from fractions import Fraction as Q
from itertools import product
from math import comb
from unittest.mock import patch
from types import SimpleNamespace
import unittest
import numpy as np
from flint import arb,ctx
import scalar_cover as sc
from mixture import group_components,G
from heterogeneous import tilted_reference
from packet_census import census,rank


def endpoint(value):
    mantissa,exponent=value.upper().man_exp()
    return Q(int(mantissa))*Q(2)**int(exponent)


class GF16Tests(unittest.TestCase):
    def test_proposal_stop_only_controls_search_effort(self):
        model=SimpleNamespace(tilt=Q(1,4),proposal_stop_bits=66)
        base=(-75,{})
        model.propose_with=lambda cell,tilt:base
        self.assertEqual(sc.Model.proposal(model,(Q(1,32),Q(1,32))),base)
        self.assertEqual(sc.proposal_cutoff(None),-90)
        for value in (True,0,-1,1.5,'66'):
            model.proposal_stop_bits=value
            with self.assertRaises(ValueError):sc.Model.proposal(model,(Q(0),Q(1)))

    def test_repropose_all_uses_only_partition_and_fresh_candidates(self):
        saved=dict(leaves={'0':dict(witness={'old':True},proposal=-1000)},
                   unresolved={'1':{}},visited=50)
        class Trial:
            root=(Q(0),Q(1))
            def __init__(self):self.proposed=[];self.checked=[]
            def proposal(self,cell):
                self.proposed.append(cell);return -100,{'fresh':list(map(str,cell))}
            def outward(self,cell,witness):
                if 'old' in witness:raise AssertionError('old witnesses must not be evaluated or imported')
                self.checked.append(cell);return arb(2)**-99
        model=Trial();result=sc.run(model,2,5,60,resume=saved,coalesce=False,repropose=True)
        self.assertEqual(set(result['leaves']),{'0','1'})
        self.assertEqual(result['unresolved'],{})
        self.assertEqual(model.proposed,model.checked)
        self.assertEqual(len(model.checked),2)
        self.assertEqual(result['visited'],52)
        self.assertTrue(saved['leaves']['0']['witness']['old'])
        for option in (1,None,'yes'):
            with self.assertRaises(ValueError):sc.run(model,2,5,60,resume=saved,repropose=option)
        with self.assertRaises(ValueError):sc.run(model,2,5,60,repropose=True)

    def test_repropose_all_cli_requires_saved_partition_without_coalescing(self):
        for argv in (['--repropose-all'],['--resume','unused.json','--repropose-all']):
            with patch('sys.argv',['scalar_cover.py',*argv]),patch('sys.stderr'),self.assertRaises(SystemExit) as caught:
                sc.main()
            self.assertEqual(caught.exception.code,2)

    def test_repropose_all_does_not_keep_old_acceptance_when_new_attempt_fails(self):
        saved=dict(leaves={'0':dict(witness={'old':True},proposal=-1000)},
                   unresolved={'1':{}},visited=50)
        class Trial:
            root=(Q(0),Q(1))
            def proposal(self,cell):return 1,{}
            def outward(self,cell,witness):raise AssertionError('no successful fresh proposal exists')
        model=Trial();result=sc.run(model,1,5,60,resume=saved,coalesce=False,repropose=True)
        self.assertEqual(result['leaves'],{})
        self.assertEqual(set(result['unresolved']),{'00','01','1'})
        cells=sc.partition(model,result['leaves'],result['unresolved'])
        self.assertEqual(sum(hi-lo for lo,hi in cells.values()),1)

    def test_retarget_claim_is_explicit_and_does_not_mutate_record(self):
        record=dict(threshold=193986,updates=3,leaves={'':{'proposal':-1000}})
        self.assertEqual(sc.retarget_claim(record),(193986,3))
        self.assertEqual(sc.retarget_claim(record,updates=4),(193986,4))
        self.assertEqual(sc.retarget_claim(record,distance='99/1000'),(207618,3))
        self.assertEqual(sc.retarget_claim(record,updates=2,distance='1/10'),(209715,2))
        self.assertEqual(record,dict(threshold=193986,updates=3,leaves={'':{'proposal':-1000}}))
        for value in (True,1,5,'4'):
            with self.assertRaises(ValueError):sc.retarget_claim(record,updates=value)
        for value in ('0','1','-1/10'):
            with self.assertRaises(ValueError):sc.retarget_claim(record,distance=value)

    def test_retarget_cli_requires_resume_and_separate_output(self):
        for argv in (['--retarget-updates','4'],
                     ['--replay','unused.json','--retarget-distance','1/10','--output','new.json'],
                     ['--resume','same.json','--retarget-updates','4','--output','same.json']):
            with patch('sys.argv',['scalar_cover.py',*argv]),patch('sys.stderr'),self.assertRaises(SystemExit) as caught:
                sc.main()
            self.assertEqual(caught.exception.code,2)

    def test_retarget_resume_rechecks_changed_claim_not_saved_scores(self):
        saved=dict(threshold=188743,updates=3,visited=2,
            leaves={'0':{'witness':{'limit':3},'proposal':-1000},
                    '1':{'witness':{'limit':4},'proposal':-1000}},unresolved={})
        class Trial:
            root=(Q(0),Q(1))
            def __init__(self):
                self.threshold,self.updates=sc.retarget_claim(saved,updates=4,distance='99/1000')
                self.checked=[]
            def outward(self,cell,witness):
                self.checked.append((cell,self.threshold,self.updates))
                return arb(2)**(-100 if self.updates<=witness['limit'] else 0)
            def proposal(self,cell):return 1,{}
        model=Trial();result=sc.run(model,1,4,70,resume=saved,coalesce=False)
        self.assertEqual(set(result['leaves']),{'1'})
        self.assertEqual(len(model.checked),2)
        self.assertTrue(all(row[1:]==(207618,4) for row in model.checked))
        self.assertEqual(sum(hi-lo for lo,hi in sc.partition(model,result['leaves'],result['unresolved']).values()),1)
        self.assertTrue(result['unresolved'])

    def test_rational_row_comparison(self):
        caps=[1,0,0,0,70,0,0,0,1]
        scale=Q(1001,1000);bias=Q(21,50)
        rows=sc.envelope(caps,scale*2**6,bias)
        for w,count in enumerate(caps):
            bound=sum(c*p**w*(1-p)**(8-w) for c,p in rows)
            self.assertGreaterEqual(bound,Q(count,comb(8,w)))
        with patch.object(sc,'authenticated_caps',return_value=caps):
            self.assertEqual(sc.actual_components(6,bias,scale),group_components(rows))
            self.assertEqual(sc.actual_components(6,bias),group_components(sc.envelope(caps,2**6,bias)))
        for scale in (0,-1):
            with self.assertRaises(ValueError):sc.actual_components(6,bias,scale)

    def test_packet_census(self):
        self.assertEqual(rank([1,2,3,4]),3)
        data=census([1<<i for i in range(8)])
        self.assertEqual(data['pair_zero_preimages_out_of_225'],{0:1})
        self.assertEqual(data['pair_largest_atom_out_of_225'],{1:1})
        self.assertEqual(data['worst_case_atom_caps'],{1:'1/15',2:'1/225'})
        data=census([1,2,4,8]*2)
        self.assertEqual(data['pair_zero_preimages_out_of_225'],{15:1})
        self.assertEqual(data['pair_largest_atom_out_of_225'],{15:1})

    def test_exact_field_action(self):
        for a in range(1,16):
            self.assertEqual(sorted(sc.multiply(a,x) for x in range(16)),list(range(16)))
            for x,y in product(range(16),repeat=2):self.assertEqual(sc.multiply(a,x^y),sc.multiply(a,x)^sc.multiply(a,y))
        for x in range(1,16):self.assertEqual(sorted(sc.multiply(a,x) for a in range(1,16)),list(range(1,16)))

    def test_activity_distribution(self):
        for p in (Q(0),Q(1,5),Q(1)):
            probs=sc.probabilities(p)
            self.assertEqual(sum(probs),1)
            self.assertEqual(probs[0],1-p)
            for w in range(1,5):self.assertEqual(probs[w],p*sum(x.bit_count()==w for x in range(1,16))/15)

    def test_binary_shuffle_comparison_after_field_action(self):
        # Exhaust all 16^2 packet outputs. The field conditional law is the
        # same on both sides, so only binary activity pays a shuffle loss.
        for p,r,t in product((Q(0),Q(2,5),Q(1)),(Q(0),Q(3,5),Q(1)),(Q(1,8),Q(1),Q(2))):
            reference,factor=tilted_reference([(1-p,p),(1-r,r)],[1,1],[1,t])
            def mass(a,x):return 1-a if x==0 else a/15
            for x,y in product(range(16),repeat=2):
                actual=(mass(p,x)*mass(r,y)+mass(r,x)*mass(p,y))/2
                self.assertLessEqual(actual,factor*mass(reference[1],x)*mass(reference[1],y))

    def test_gf_inner_against_exhaustive_reference(self):
        ctx.prec=192
        images=[sum(((s>>i)&1)*v for i,v in enumerate((1,6,120))) for s in range(8)]
        columns=[1,2,4,3,5,7,6,1]
        data=sc.kernel.prepare(images,columns,3,2)
        z=Q(3,4)
        for activity in (Q(0),Q(1,5),Q(1)):
            matrix=sc.kernel.outward_at_z(data,sc.probabilities(activity),arb(3)/4)
            for state in range(8):
                actual=[Q(0),Q(0)]
                for x in range(256):
                    probability=Q(1)
                    for offset in (0,4):probability*=activity/15 if ((x>>offset)&15) else 1-activity
                    feedback=0
                    for b,c in enumerate(columns):
                        if x>>b&1:feedback^=c
                    weighted=probability*z**(images[state]^x).bit_count()
                    if state==0:actual[int(feedback!=0)]+=weighted
                    else:
                        for refreshed in range(1,8):
                            chance=Q(int(refreshed==state),4)+Q(3,28)
                            actual[int(refreshed!=feedback)]+=weighted*chance
                for j in range(2):self.assertGreaterEqual(endpoint(matrix[int(state!=0),j]),actual[j])

    def test_scalar_bound_and_alternate_weights(self):
        ctx.prec=192
        components=group_components([(Q(1),p) for p in (Q(0),Q(1),Q(1,2),Q(2,5),Q(3,5))])
        data=sc.kernel.prepare(list(range(8)),[1,2,4,3]*32,3)
        model=sc.Model(components,data,10485)
        self.assertEqual(sum(c for c,_,_ in model.components),sum(row[1] for row in components))
        for tilt in (model.tilt,Q(1,8),Q(1)):
            for i,(c,p,a) in enumerate(model.components):
                if not a:continue
                point=Q(65,G)*model.features[i];cell=(max(model.minimum,point-Q(1,65536)),min(Q(1),point+Q(1,65536)))
                weights,duals=model.weights(cell,tilt)
                z=1-p+p*tilt
                expected=(1-Q(65,G)+Q(65,G)*(1-p)/z,Q(65,G)*p/z)
                for x,y in zip(expected,weights):self.assertLessEqual(x,y)
        cell=(Q(1,20),Q(51,1000))
        for tilt in (model.tilt,Q(1,8)):
            score,witness=model.propose_with(cell,tilt)
            bound=model.outward(cell,witness)
            self.assertAlmostEqual(score,float(bound.log()/arb(2).log()),places=5)
        del witness['weights_dual']
        with self.assertRaises(ValueError):model.outward(cell,witness)

    def test_work_limit_keeps_full_partition(self):
        class Failing:
            root=(Q(0),Q(1))
            def proposal(self,cell):return 1,{}
        model=Failing();result=sc.run(model,9,3,70)
        self.assertTrue(result['unresolved'])
        cells=sc.partition(model,result['leaves'],result['unresolved'])
        self.assertEqual(sum(hi-lo for lo,hi in cells.values()),1)

    def test_variance_shuffle_witness_and_fallback(self):
        ctx.prec=192
        components=group_components([(Q(1),p) for p in (Q(0),Q(1),Q(1,2),Q(2,5),Q(3,5))])
        data=sc.kernel.prepare(list(range(8)),[1,2,4,3]*32,3)
        model=sc.Model(components,data,10485,q_min=129,tilt=Q(1,8),variance_shuffle=True)
        cell=(Q(1,25),Q(401,10000))
        loss=model.comparison_loss(cell,model.tilt)
        self.assertLess(loss,model.loss(cell))
        self.assertEqual(model.comparison_loss(cell,Q(1)),model.loss(cell))
        for tilt in (model.tilt,Q(1)):
            score,witness=model.propose_with(cell,tilt)
            bound=model.outward(cell,witness)
            self.assertAlmostEqual(score,float(bound.log()/arb(2).log()),places=5)
            if tilt==model.tilt:
                self.assertIn('variance_dual',witness)
                wrong=dict(witness,variance_dual=['0','-1'])
                with self.assertRaises(ValueError):model.outward(cell,wrong)
                del witness['variance_dual']
                with self.assertRaises(ValueError):model.outward(cell,witness)

    def test_resume_rechecks_witnesses_and_keeps_unresolved_cells(self):
        class Trial:
            root=(Q(0),Q(1))
            def __init__(self):self.checked=[]
            def outward(self,cell,witness):
                self.checked.append(cell);return arb(2)**(-100 if witness['valid'] else 0)
            def proposal(self,cell):return 1,{}
        model=Trial()
        saved=dict(leaves={'00':dict(witness={'valid':True}),'01':dict(witness={'valid':False})},
                   unresolved={'1':dict(cell=['bogus coordinates are ignored'])},visited=5)
        result=sc.run(model,1,4,70,resume=saved)
        self.assertEqual(set(result['leaves']),{'00'})
        self.assertEqual(len(model.checked),2)
        cells=sc.partition(model,result['leaves'],result['unresolved'])
        self.assertEqual(sum(hi-lo for lo,hi in cells.values()),1)
        self.assertEqual(result['visited'],6)
        with self.assertRaises(ValueError):
            sc.run(model,1,4,70,resume=dict(leaves=saved['leaves'],unresolved={},visited=5))

    def test_resume_coalesces_only_unresolved_siblings(self):
        class Trial:
            root=(Q(0),Q(1))
            def __init__(self):self.probed=[]
            def outward(self,cell,witness):return arb(2)**-100
            def proposal(self,cell):self.probed.append(cell);return -100,{}
        model=Trial()
        saved=dict(leaves={'1':dict(witness={})},unresolved={p:{} for p in ('000','001','01')},visited=7)
        result=sc.run(model,1,4,70,resume=saved)
        self.assertEqual(model.probed,[(Q(0),Q(1,2))])
        self.assertEqual(set(result['leaves']),{'0','1'});self.assertFalse(result['unresolved'])

    def test_outer_dual_stationarity_and_boundary(self):
        from scipy.special import logsumexp
        features=np.array([0.,.02,.1,.4]);logs=np.array([0.,-10.,-12.,-20.]);active=np.array([0,1,1,1])
        for target,occupancy,positive in [(.04,.047,True),(.04,.7,True),(.0001,.047,False),(.1,1.,True)]:
            eta,mu=sc.outer_dual(logs,features,active,target,occupancy,positive)
            values=logs+eta*features+mu*active;prob=np.exp(values-logsumexp(values))
            if 0<abs(eta)<19999:self.assertAlmostEqual(float(prob@features),target,places=10)
            if 0<mu<19999:self.assertAlmostEqual(float(prob@active),occupancy,places=10)
            self.assertGreaterEqual(mu,0)
            self.assertGreaterEqual(eta if positive else -eta,0)

    def test_resume_can_keep_subdivisions_without_trusting_saved_bounds(self):
        class Trial:
            root=(Q(0),Q(1))
            def __init__(self):self.probed=[];self.checked=[]
            def outward(self,cell,witness):
                self.checked.append(cell);return arb(2)**(-100 if witness.get('valid',True) else 0)
            def proposal(self,cell):self.probed.append(cell);return -100,{}
        model=Trial()
        saved=dict(leaves={'1':dict(witness={'valid':False})},
            unresolved={p:{} for p in ('000','001','01')},visited=7)
        result=sc.run(model,2,4,70,resume=saved,coalesce=False)
        self.assertEqual(model.probed,[(Q(1,2),Q(1)),(Q(1,4),Q(1,2))])
        self.assertEqual(model.checked[0],(Q(1,2),Q(1)))
        self.assertEqual(set(result['leaves']),{'1','01'})
        self.assertEqual(set(result['unresolved']),{'000','001'})
        cells=sc.partition(model,result['leaves'],result['unresolved'])
        self.assertEqual(sum(hi-lo for lo,hi in cells.values()),1)
        self.assertEqual(result['visited'],9)
        with self.assertRaises(ValueError):sc.run(model,1,4,70,resume=saved,coalesce='no')


if __name__=='__main__':unittest.main()
