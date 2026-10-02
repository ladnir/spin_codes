from fractions import Fraction as Q
from math import comb,factorial
from copy import deepcopy
from types import SimpleNamespace
from itertools import product
import unittest
from unittest.mock import patch
import numpy as np
from flint import arb,arb_mat,ctx,fmpq_mat,fmpq
import regional_count as regional
import scalar_cover as sc
import birth_classes
from occupancy_model import placement as reference_placement
from test_scalar_cover import endpoint


class OutwardRegionalCountTests(unittest.TestCase):
    def test_proposal_local_reuse_never_replaces_fresh_verification(self):
        ctx.prec=192;model=SimpleNamespace(data={'windows':32})
        local=[arb_mat([[1,0],[0,1]])];caps=[arb(1)];refined=[local[0]*2]
        witness={'parameters':['1/8','0','0'],'regional_lazy_density_through':6}
        with patch('occupancy_birth_classes.outward',return_value=local) as build, \
                patch('lazy_density.actual_caps',return_value=caps) as census, \
                patch('lazy_density.candidate',return_value=refined) as candidate, \
                patch('lazy_density.refine_actual',return_value=refined) as fresh:
            for _ in range(2):
                result=regional.local_operators(model,witness,_proposal_scratch=regional.proposal_cache(model))
                self.assertIs(result,refined)
            self.assertEqual(build.call_count,1);self.assertEqual(census.call_count,1)
            self.assertEqual(candidate.call_count,2);fresh.assert_not_called()
            # Verification does not accept the proposal's numeric scratch.
            regional.local_operators(model,witness)
            self.assertEqual(build.call_count,2);self.assertEqual(fresh.call_count,1)
            changed=dict(witness,parameters=['1/9','0','0'])
            regional.local_operators(model,changed,_proposal_scratch=regional.proposal_cache(model))
            self.assertEqual(build.call_count,3);self.assertEqual(census.call_count,2)
            ctx.prec=256
            regional.local_operators(model,changed,_proposal_scratch=regional.proposal_cache(model))
            self.assertEqual(build.call_count,4);self.assertEqual(census.call_count,3)
            model.data=dict(model.data)
            regional.local_operators(model,changed,_proposal_scratch=regional.proposal_cache(model))
            self.assertEqual(build.call_count,5);self.assertEqual(census.call_count,4)
            self.assertLessEqual(len(regional.proposal_cache(model)),2)

    def test_proposal_count_weights_reuse_only_identical_domains_and_duals(self):
        ctx.prec=192
        model=SimpleNamespace(data={},features=[Q(1,4)],active=[1],q_min=1,tilt=Q(3,16),
            weights=lambda cell,tilt:([Q(3,4),Q(1,4)],None))
        cell=(Q(1,5),Q(3,10));parts=[((Q(0),Q(1,4)),(Q(0),Q(0),Q(0)))]
        witness={'parameters':['1/8','0','0'],'variance_dual':['0','0','0'],'regional_direct_counts':True}
        checked=dict(witness,regional_count_parts=[dict(interval=['0','1/4'],dual=['0','0','0'],mgf_witnesses=[])])
        scratch=regional.proposal_cache(model)
        with patch.object(sc,'G',4),patch('variance_partition.factor',return_value=Q(1)), \
                patch.object(regional,'count_mass_caps',return_value=[arb(1)/5]*5) as count:
            scale,first=regional.proposal_count_weights(model,cell,witness,parts,checked,scratch)
            self.assertEqual(scale,1)
            changed=dict(witness,parameters=['1/4','0','0'],regional_exact_zero=True)
            scale,again=regional.proposal_count_weights(model,cell,changed,parts,checked,scratch)
            self.assertIs(first,again);self.assertEqual(count.call_count,1)
            other=deepcopy(checked);other['regional_count_parts'][0]['mgf_witnesses']=[{'tilt':'1','dual':['0','0','0']}]
            regional.proposal_count_weights(model,cell,witness,parts,other,scratch)
            self.assertEqual(count.call_count,2)
            regional.proposal_count_weights(model,(Q(1,5),Q(1,3)),witness,parts,other,scratch)
            self.assertEqual(count.call_count,3)
            ctx.prec=256
            regional.proposal_count_weights(model,cell,witness,parts,other,scratch)
            self.assertEqual(count.call_count,4)
            model.q_min=2
            regional.proposal_count_weights(model,cell,witness,parts,other,scratch)
            self.assertEqual(count.call_count,5)
            self.assertEqual(list(scratch),['count_weights'])

    def test_exact_zero_dispatch_is_opt_in_and_recomputed(self):
        ctx.prec=192;data={'windows':32};model=SimpleNamespace(data=data)
        local=[arb_mat([[1,0],[0,1]])];refined=[local[0]*arb(1)]
        witness={'parameters':['1/8','0','0']}
        with patch('occupancy_birth_classes.outward',return_value=local), \
                patch('occupancy_birth_classes.refine_zero',return_value=refined) as refine:
            self.assertIs(regional.local_operators(model,witness),local)
            refine.assert_not_called()
            for _ in range(2):
                self.assertIs(regional.local_operators(model,dict(witness,regional_exact_zero=True)),refined)
            self.assertEqual(refine.call_count,2)
            self.assertIs(refine.call_args.args[0],data)
            self.assertIs(refine.call_args.args[1],local)
            for invalid in (1,'true',None):
                with self.assertRaises(ValueError):regional.local_operators(model,dict(witness,regional_exact_zero=invalid))

    def test_feedback_class_selector_requires_a_complete_positive_integer_interval(self):
        self.assertIsNone(regional.feedback_class_interval({},32))
        self.assertEqual(regional.feedback_class_interval(dict(
            regional_feedback_classes_from=3,regional_feedback_classes_through=32),32),(3,32))
        for start,end in ((None,32),(3,None),(0,32),(3,33),(4,3),(True,32),(3,'32')):
            selected={key:value for key,value in zip(
                ('regional_feedback_classes_from','regional_feedback_classes_through'),(start,end)) if value is not None}
            with self.assertRaises(ValueError):regional.feedback_class_interval(selected,32)
        for selected in ({'regional_feedback_uniform_classes':True},
                {'regional_feedback_uniform_replace':True},
                {'regional_feedback_uniform_classes':1},
                {'regional_feedback_uniform_replace':'true'}):
            with self.assertRaises(ValueError):regional.feedback_class_interval(selected,32)

    def test_feedback_class_allocation_is_opt_in_fresh_and_after_other_refinements(self):
        ctx.prec=192;data={'windows':32};model=SimpleNamespace(data=data)
        local=[arb_mat([[1,0],[0,1]])];refined=[local[0]*2];attached={'windows':32,'attached':True}
        witness={'parameters':['1/8','0','0'],'regional_lazy_density_through':6,
            'regional_feedback_classes_from':3,'regional_feedback_classes_through':32}
        with patch('occupancy_birth_classes.outward',return_value=local), \
                patch('lazy_density.refine_actual',return_value=refined), \
                patch('fiber_density.attach',return_value=attached) as attach, \
                patch('fiber_density.candidate',return_value=refined) as candidate:
            regional.local_operators(model,witness);regional.local_operators(model,witness)
            attach.assert_called_once_with(data)
            self.assertEqual(candidate.call_count,2)
            self.assertIs(candidate.call_args.args[0],attached)
            self.assertIs(candidate.call_args.args[1],refined)
            self.assertEqual(candidate.call_args.args[3:],(3,32))
            self.assertEqual(candidate.call_args.kwargs,{'allocation':'classes'})
            regional.local_operators(SimpleNamespace(data=data),witness)
            self.assertEqual(attach.call_count,2)
            regional.local_operators(model,{'parameters':['1/8','0','0']})
            self.assertEqual(candidate.call_count,3)
            with patch('fiber_density.replace_uniform_classes',return_value=refined) as replacement:
                regional.local_operators(model,dict(witness,regional_feedback_uniform_classes=True,
                    regional_feedback_uniform_replace=True))
                self.assertIs(replacement.call_args.args[1],refined)
                self.assertIs(replacement.call_args.args[2],local)
                self.assertEqual(replacement.call_args.args[4:],(3,32))

    def test_joint_census_is_fresh_per_model_and_precedes_density_refinement(self):
        ctx.prec=192
        model=SimpleNamespace(data={'windows':32});local=[arb_mat([[1,0],[0,1]])]
        joint_local=[arb_mat([[2,0],[0,2]])];checked=[{'occupancy':0}];order=[]
        witness={'parameters':['1/8','0','0'],'regional_joint_return_through':3,'regional_lazy_density_through':6}
        def joint(data,operators,census,z):
            self.assertIs(operators,local);self.assertIs(census,checked);order.append('joint');return joint_local
        def density(data,operators,z,cutoff):
            self.assertIs(operators,joint_local);order.append('density');return joint_local
        with patch('occupancy_birth_classes.outward',return_value=local), \
             patch('return_moment.actual_census',return_value=checked) as build, \
             patch('return_moment.refine_class_returns',side_effect=joint), \
             patch('lazy_density.refine_actual',side_effect=density):
            regional.local_operators(model,witness);regional.local_operators(model,witness)
            self.assertEqual(build.call_count,1)
            regional.local_operators(SimpleNamespace(data=model.data),witness)
            self.assertEqual(build.call_count,2)
            self.assertEqual(order,['joint','density']*3)

    def test_local_density_witness_dispatches_fresh_refinement(self):
        ctx.prec=192
        model=SimpleNamespace(data={'windows':32})
        local=[arb_mat([[1,0],[0,1]])]
        witness={'parameters':['1/8','0','0']}
        with patch('occupancy_birth_classes.outward',return_value=local), \
             patch('lazy_density.refine_actual',return_value=local) as refine:
            self.assertIs(regional.local_operators(model,witness),local)
            refine.assert_not_called()
            selected=dict(witness,regional_lazy_density_through=6)
            self.assertIs(regional.local_operators(model,selected),local)
            refine.assert_called_once()
            self.assertEqual(refine.call_args.args[0],model.data)
            self.assertIs(refine.call_args.args[1],local)
            self.assertEqual(refine.call_args.args[3],6)

    def test_scalar_proposal_records_direct_count_mode(self):
        components=[('a',Q(1),sc.probabilities(Q(1,2)),1)]
        model=sc.Model(components,{'windows':32},0,inner=birth_classes,
            variance_shuffle=True,variance_bins=2,regional_count=True)
        cell=(Q(1,4),Q(1,4)+Q(1,4096))
        witness={'parameters':['1','0','0']}
        with patch.object(model,'propose_with',return_value=(0,witness)), \
             patch('variance_partition.propose',return_value=(-10,witness)), \
             patch.object(regional,'propose',side_effect=lambda m,c,w:(-100,w)) as proposed:
            score,saved=model.proposal(cell)
            self.assertEqual(score,-100)
            proposed.assert_called_once()
            self.assertIs(saved['regional_direct_counts'],True)
            self.assertIs(saved['regional_exact_zero'],True)
            self.assertIs(saved['regional_tilted_atom'],True)
            self.assertNotIn('regional_direct_counts',witness)

    def test_scalar_proposal_can_retain_lazy_density_without_discarding_baseline(self):
        components=[('a',Q(1),sc.probabilities(Q(1,2)),1)]
        model=sc.Model(components,{'windows':32},0,inner=birth_classes,
            variance_shuffle=True,variance_bins=2,regional_count=True)
        cell=(Q(1,4),Q(1,4)+Q(1,4096));witness={'parameters':['1','0','0']}
        for refinement,combined in ((-100,-100),(10,10),(-50,-100)):
            def proposal(m,c,w):
                return (combined if 'regional_joint_return_through' in w else
                        refinement if 'regional_lazy_density_through' in w else -30),w
            with patch.object(model,'propose_with',return_value=(0,witness)), \
                 patch('variance_partition.propose',return_value=(-10,witness)), \
                 patch.object(regional,'propose',side_effect=proposal) as proposed:
                score,saved=model.proposal(cell)
                self.assertEqual(score,min(-30,refinement,combined))
                self.assertEqual(proposed.call_count,2 if refinement==-100 else 3 if combined==-100 else 8)
                self.assertEqual('regional_lazy_density_through' in saved,score<-30)
                self.assertEqual('regional_joint_return_through' in saved,refinement==-50)
                if score<-30:self.assertEqual(saved['regional_lazy_density_through'],6)

    def test_scalar_proposal_selects_classified_feedback_only_when_better(self):
        components=[('a',Q(1),sc.probabilities(Q(1,2)),1)]
        model=sc.Model(components,{'windows':32},0,inner=birth_classes,
            variance_shuffle=True,variance_bins=2,regional_count=True)
        cell=(Q(1,4),Q(1,4)+Q(1,4096));witness={'parameters':['1','0','0']}
        for classified_score in (-125,42):
            def proposal(m,c,w):
                return (classified_score if 'regional_feedback_classes_from' in w else -30),w
            with patch.object(model,'propose_with',return_value=(0,witness)), \
                    patch('variance_partition.propose',return_value=(-10,witness)), \
                    patch.object(regional,'propose',side_effect=proposal) as proposed:
                score,saved=model.proposal(cell)
                self.assertEqual(score,min(-30,classified_score))
                self.assertEqual(proposed.call_count,4 if classified_score<0 else 8)
                if classified_score<0:
                    self.assertEqual(saved['regional_feedback_classes_from'],3)
                    self.assertEqual(saved['regional_feedback_classes_through'],32)
                else:self.assertNotIn('regional_feedback_classes_from',saved)

    def test_scalar_proposal_can_select_complete_uniform_row_replacement(self):
        components=[('a',Q(1),sc.probabilities(Q(1,2)),1)]
        model=sc.Model(components,{'windows':32},0,inner=birth_classes,
            variance_shuffle=True,variance_bins=2,regional_count=True)
        cell=(Q(1,4),Q(1,4)+Q(1,4096));witness={'parameters':['1','0','0']}
        with patch.object(model,'propose_with',return_value=(0,witness)), \
                patch('variance_partition.propose',return_value=(-10,witness)), \
                patch.object(regional,'propose',side_effect=lambda m,c,w:(
                    -120 if w.get('regional_feedback_uniform_replace') else -30,w)) as proposed:
            score,saved=model.proposal(cell)
            self.assertEqual(score,-120);self.assertEqual(proposed.call_count,5)
            self.assertIs(saved['regional_feedback_uniform_classes'],True)
            self.assertIs(saved['regional_feedback_uniform_replace'],True)

    def test_scalar_regional_tilt_search_retains_baseline_and_exact_parameters(self):
        components=[('a',Q(1),sc.probabilities(Q(1,2)),1)]
        model=sc.Model(components,{'windows':32},0,inner=birth_classes,
            variance_shuffle=True,variance_bins=2,regional_count=True)
        cell=(Q(1,4),Q(1,4)+Q(1,4096));witness={'parameters':['1/8','2','3']}
        for improved in (True,False):
            tried=[]
            def proposal(m,c,w):
                lam=Q(w['parameters'][0]);tried.append(lam)
                self.assertEqual(w['parameters'][1:],['2','3'])
                return (-120 if improved and lam==Q(1,10) else -30 if lam==Q(1,8) else 10),w
            with patch.object(model,'propose_with',return_value=(0,witness)), \
                    patch('variance_partition.propose',return_value=(-10,witness)), \
                    patch.object(regional,'propose',side_effect=proposal):
                score,saved=model.proposal(cell)
            self.assertEqual(score,-120 if improved else -30)
            self.assertEqual(Q(saved['parameters'][0]),Q(1,10) if improved else Q(1,8))
            self.assertEqual(tried[-2:] if improved else tried[-3:],
                [Q(9,80),Q(1,10)] if improved else [Q(9,80),Q(1,10),Q(39,320)])
            self.assertEqual(witness,{'parameters':['1/8','2','3']})

    def test_binomial_atom_maximum_over_mean_interval(self):
        ctx.prec=192
        for n in (4,7,16):
            for cell in ((Q(1,10),Q(4,5)),(Q(1,4),Q(1,4)),(Q(49,100),Q(51,100))):
                bounds=regional.binomial_interval_upper(n,cell)
                for j,value in enumerate(bounds):
                    p=min(cell[1],max(cell[0],Q(j,n)))
                    exact=comb(n,j)*p**j*(1-p)**(n-j)
                    self.assertGreaterEqual(endpoint(value),exact)
                    self.assertLess(endpoint(value)-exact,Q(1,10**40))
        for n,cell in ((0,(Q(1,4),Q(1,2))),(4,(Q(0),Q(1,2))),(4,(Q(3,4),Q(1,2)))):
            with self.assertRaises(ValueError):regional.binomial_interval_upper(n,cell)

    def test_polynomial_matrix_power_matches_ordered_placement(self):
        ctx.prec=192
        exact=[fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]]),fmpq_mat([[2,1],[0,1]])]
        local=[arb_mat([[int(m[i,j]) for j in range(2)] for i in range(2)]) for m in exact]
        for epochs in (1,2,3,5,8):
            actual=regional.placement(local,epochs)
            expected=reference_placement(exact,epochs,2,fmpq_mat,lambda x:x,maximum_groups=2*epochs)
            for matrix,truth in zip(actual,expected):
                for i in range(2):
                    for j in range(2):
                        value=Q(str(truth[i,j]))
                        self.assertGreaterEqual(endpoint(matrix[i,j]),value)
                        self.assertLess(endpoint(matrix[i,j])-value,Q(1,10**35))

    def test_binomial_recurrence_and_checked_mgf_caps(self):
        ctx.prec=192;features=list(map(Q,('0','1/5','1/2','4/5','1')));active=[0,1,1,1,1]
        for probabilities in ([features[1]]*4,[features[i] for i in (0,1,3,4)],
                              [features[i] for i in (1,2,2,3)],[features[i] for i in (0,0,4,4)]):
            mean=sum(probabilities)/4;var=sum(p*(1-p) for p in probabilities)/4
            cell=(mean-Q(1,100),mean+Q(1,100));interval=(max(Q(0),var-Q(1,100)),min(Q(1,4),var+Q(1,100)))
            q=sum(p>0 for p in probabilities)
            witnesses=regional.propose_mgf(features,active,cell,interval,q,4)
            self.assertTrue(witnesses)
            caps=regional.count_ratios(features,active,cell,interval,q,4,arb(10),witnesses)
            refined=regional.count_ratios(features,active,cell,interval,q,4,arb(10),
                [dict(w,tilted_atom=True) for w in witnesses])
            variance_witnesses=regional.propose_mgf(features,active,cell,interval,q,4,
                tilted_atom=True,tilted_variance=True)
            finite=regional.count_ratios(features,active,cell,interval,q,4,arb(10),variance_witnesses)
            direct=regional.count_mass_caps(features,active,cell,interval,q,4,arb(10),variance_witnesses)
            fallback=regional.count_mass_caps(features,active,cell,interval,q,4,arb(10),[])
            masses=[Q(1)]
            for p in probabilities:
                new=[Q(0)]*(len(masses)+1)
                for j,mass in enumerate(masses):new[j]+=mass*(1-p);new[j+1]+=mass*p
                masses=new
            for j,(value,mass) in enumerate(zip(regional.binomial_masses(4,mean),masses)):
                binomial=comb(4,j)*mean**j*(1-mean)**(4-j)
                self.assertTrue(value.contains(regional.aq(binomial)))
                self.assertGreaterEqual(endpoint(caps[j])*binomial,mass)
                self.assertGreaterEqual(endpoint(refined[j])*binomial,mass)
                self.assertLessEqual(endpoint(refined[j]),endpoint(caps[j]))
                self.assertGreaterEqual(endpoint(finite[j])*binomial,mass)
                self.assertLessEqual(endpoint(finite[j]),endpoint(refined[j]))
                self.assertGreaterEqual(endpoint(direct[j]),mass)
                self.assertLessEqual(endpoint(direct[j]),endpoint(fallback[j]))
                self.assertLessEqual(endpoint(direct[j]),1)
            for witness in witnesses:
                # A higher-precision truth avoids comparing two equally
                # loose upper endpoints when the actual bound is equality.
                ctx.prec=384
                t=regional.aq(Q(witness['tilt'])).exp();truth=arb(1)
                for p in probabilities:truth*=regional.aq(1-p)+regional.aq(p)*t
                ctx.prec=192
                upper=regional.mgf_upper(features,active,cell,interval,q,4,witness).exp()
                self.assertGreaterEqual(endpoint(upper),-endpoint(-truth))
            for witness in variance_witnesses:
                lower=regional.tilted_variance_lower(features,active,cell,interval,q,4,witness)
                ctx.prec=384;e=regional.aq(Q(witness['tilt'])).exp()
                truth=sum((regional.aq(p*(1-p))*e/(regional.aq(1-p)+regional.aq(p)*e)**2
                    for p in probabilities),arb(0))
                self.assertLessEqual(endpoint(lower),endpoint(truth))
                ctx.prec=192

    def test_invalid_dual_sign(self):
        with self.assertRaises(ValueError):
            regional.mgf_upper([Q(1,2)],[1],(Q(1,2),Q(1,2)),(Q(1,4),Q(1,4)),1,1,
                {'tilt':'1','dual':['0','0','1']})
        with self.assertRaises(ValueError):
            regional.tilted_variance_lower([Q(1,2)],[1],(Q(1,2),Q(1,2)),(Q(1,4),Q(1,4)),1,1,
                {'tilt':'1','dual':['0','0','0'],'tilted_variance_dual':['0','0','-1']})

    def test_tilted_atom_against_exact_rational_count_laws(self):
        ctx.prec=192
        features=list(map(Q,('0','1/5','1/2','4/5','1')))
        for probabilities in (features,[Q(1,2)]*8,[Q(0)]*4+[Q(1)]*4):
            v=sum(p*(1-p) for p in probabilities)/len(probabilities)
            for exponential in (Q(1,4),Q(1,2),Q(1),Q(2),Q(4)):
                t=regional.aq(exponential).log();masses=[Q(1)]
                for p in probabilities:
                    tilted=p*exponential/(1-p+p*exponential)
                    new=[Q(0)]*(len(masses)+1)
                    for j,m in enumerate(masses):new[j]+=m*(1-tilted);new[j+1]+=m*tilted
                    masses=new
                # Use an exactly representable log tilt proposal and compute
                # the truth at its actual (nearby) tilt with higher precision.
                exact_t=Q(str(float(t)))
                bound=regional.tilted_atom_upper(features,v,len(probabilities),exact_t)
                ctx.prec=384;e=regional.aq(exact_t).exp();actual=[arb(1)]
                for p in probabilities:
                    tilted=regional.aq(p)*e/(regional.aq(1-p)+regional.aq(p)*e)
                    new=[arb(0)]*(len(actual)+1)
                    for j,m in enumerate(actual):new[j]+=m*(1-tilted);new[j+1]+=m*tilted
                    actual=new
                self.assertGreaterEqual(endpoint(bound),max(-endpoint(-m) for m in actual))
                self.assertAlmostEqual(float(max(actual)),float(max(masses)),places=12)
                ctx.prec=192

    def test_regional_sum_covers_exact_compositions_and_replays_saved_duals(self):
        ctx.prec=192;n=4;regions=2
        ps=list(map(Q,('0','1/4','3/4','1')))
        components=[(str(i),Q(i+1),sc.probabilities(p),int(i!=0)) for i,p in enumerate(ps)]
        local_exact=[fmpq_mat([[1,0],[0,1]]),fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]])]
        local=[arb_mat([[int(m[i,j]) for j in range(2)] for i in range(2)]) for m in local_exact]
        region_exact=reference_placement(local_exact,2,2,fmpq_mat,lambda x:x,maximum_groups=n)
        ordered=regional.placement(local,2)
        with patch.object(sc,'G',n),patch.object(sc,'REGIONS',regions),patch.object(sc,'PACKETS',n*regions):
            model=sc.Model(components,{'windows':32},0,q_min=2,tilt=Q(2,3),inner=birth_classes,
                variance_shuffle=True,variance_bins=2,regional_count=True)
            cell=(Q(1,4),Q(3,4));cs,_,_=model.family(model.tilt)
            witness=dict(tilt=str(model.tilt),parameters=['1','0','0'],
                variance_dual=list(map(str,model.shuffle_dual(cell))),variance_partition=[
                    dict(interval=['0','1/8'],dual=['0','0','0']),
                    dict(interval=['1/8','1/4'],dual=['0','0','0'])])
            with patch('occupancy_birth_classes.outward',return_value=local),patch.object(regional,'placement',return_value=ordered):
                upper,checked=regional.outward(model,cell,witness)
                with patch.object(regional,'propose_mgf',side_effect=AssertionError('replay must not solve LPs')):
                    self.assertEqual(upper,model.outward(cell,checked))
                expected=Q(0)
                for counts in product(range(n+1),repeat=4):
                    if sum(counts)!=n or sum(a*c for a,c in zip(model.active,counts))<2:continue
                    mean=sum(f*c for f,c in zip(model.features,counts))/n
                    if not cell[0]<=mean<=cell[1]:continue
                    masses=[Q(1)];coefficient=Q(factorial(n))
                    for c,weight,p in zip(counts,cs,model.features):
                        coefficient*=weight**c/factorial(c)
                        for _ in range(c):
                            new=[Q(0)]*(len(masses)+1)
                            for j,mass in enumerate(masses):new[j]+=mass*(1-p);new[j+1]+=mass*p
                            masses=new
                    matrix=fmpq_mat(2,2)
                    for j,(mass,value) in enumerate(zip(masses,region_exact)):
                        matrix+=value*fmpq(str(mass/model.tilt**j))
                    power=matrix**regions
                    expected+=coefficient*sum(Q(str(power[0,j])) for j in range(2))
                self.assertGreaterEqual(endpoint(upper),expected)
                direct,checked_direct=regional.outward(model,cell,dict(witness,regional_direct_counts=True))
                self.assertGreaterEqual(endpoint(direct),expected)
                self.assertLessEqual(endpoint(direct),endpoint(upper))
                with patch.object(regional,'propose_mgf',side_effect=AssertionError('replay must not solve LPs')):
                    self.assertEqual(direct,model.outward(cell,checked_direct))
                numeric=np.array([[[float(m[i,j]) for j in range(2)] for i in range(2)] for m in ordered])
                with patch('regional_count_probe.scaled_placement',return_value=(numeric,np.zeros(len(numeric)))):
                    proposed,_=regional.propose(model,cell,checked_direct)
                    self.assertAlmostEqual(proposed,float(direct.log()/arb(2).log()),places=10)
                    with patch.object(regional,'count_mass_caps',side_effect=AssertionError('proposal should reuse identical weights')):
                        repeated,_=regional.propose(model,cell,checked_direct)
                    self.assertEqual(repeated,proposed)
                    # Numeric proposal caches cannot feed the outward bound.
                    with patch.object(regional,'proposal_cache',side_effect=AssertionError('outward must regenerate')):
                        self.assertEqual(model.outward(cell,checked_direct),direct)
                with self.assertRaises(ValueError):
                    model.outward(cell,dict(checked_direct,regional_direct_counts='yes'))
                with self.assertRaises(ValueError):
                    model.outward(cell,dict(checked_direct,regional_exact_zero='yes'))
                for cutoff in (-1,33,True,'6'):
                    with self.assertRaises(ValueError):
                        model.outward(cell,dict(checked_direct,regional_lazy_density_through=cutoff))
                for cutoff in (-1,5,True,'3'):
                    with self.assertRaises(ValueError):
                        model.outward(cell,dict(checked_direct,regional_joint_return_through=cutoff))
                for broken in ([],checked['regional_count_parts'][::-1]):
                    with self.assertRaises(ValueError):model.outward(cell,dict(checked,regional_count_parts=broken))
                changed=deepcopy(checked);changed['regional_count_parts'][0]['dual'][0]='1'
                with self.assertRaises(ValueError):model.outward(cell,changed)
                changed=deepcopy(checked);changed['regional_count_parts'][0]['mgf_witnesses'][0]['dual'][2]='1'
                with self.assertRaises(ValueError):model.outward(cell,changed)
                model.regional_count=False
                with self.assertRaises(ValueError):model.outward(cell,checked)

    def test_invalid_regional_mode(self):
        components=[('a',Q(1),sc.probabilities(Q(1,2)),1)]
        for kwargs in ({'regional_count':'true'},{'regional_count':True},
                       {'regional_count':True,'variance_shuffle':True,'variance_bins':2}):
            with self.assertRaises(ValueError):sc.Model(components,{},0,**kwargs)


if __name__=='__main__':unittest.main()
