from fractions import Fraction as Q
from math import factorial
from unittest.mock import patch
import unittest

from flint import arb,arb_mat,ctx

from mixture_cover import outer_log_bound,children,cover,replay,convex_hull,separated,CoverModel,G,L,linear_upper,resolve_threshold
from probe import aq


class MixtureCoverTests(unittest.TestCase):
    def test_parent_witness_is_recomputed_for_each_child(self):
        class Toy:
            def empty(self,cell):return False
            def proposal(self,cell):return -20,['parent']
            def outward(self,cell,witness):
                return arb(2)**(-20 if cell[1]-cell[0]>Q(1,2) else -100)
        instance=Toy()
        with patch('builtins.print'),patch.object(instance,'proposal',wraps=instance.proposal) as propose:
            result=cover(instance,max_cells=3,reuse_witnesses=True)
        self.assertEqual(propose.call_count,1)
        self.assertFalse(result['unresolved'])
        self.assertEqual(replay(instance,result),arb(2)**-99)
        with patch('builtins.print'):
            partial=cover(instance,max_cells=1,reuse_witnesses=True)
            resumed=cover(instance,max_cells=2,resume=partial,reuse_witnesses=True)
        self.assertEqual(set(partial['unresolved']),{'0','1'})
        self.assertTrue(all(leaf['witness']==['parent'] for leaf in partial['unresolved'].values()))
        self.assertEqual(replay(instance,resumed),arb(2)**-99)
        with self.assertRaises(ValueError):cover(instance,reuse_witnesses=True,screen_only=True)

    def test_affine_weight_majorants_for_all_small_compositions(self):
        features=[(Q(1),Q(0),Q(0)),(Q(1,4),Q(1,2),Q(1,4)),(Q(0),Q(0),Q(1))]
        active=[0,1,1];values=[Q(2),Q(5,3),Q(1,4)]
        cell=(Q(1,8),Q(3,8),Q(1,8),Q(3,4));rho=Q(1,3)
        for dual in ((0,0,0),(2,-3,1),(-5,7,3)):
            upper=linear_upper(values,features,active,rho,cell,dual)
            for a in range(7):
                for b in range(7-a):
                    f=[Q(a,6),Q(b,6),Q(6-a-b,6)]
                    x=sum(w*p[1] for w,p in zip(f,features));y=sum(w*p[2] for w,p in zip(f,features))
                    if f[1]+f[2]<rho or not cell[0]<=x<=cell[1] or not cell[2]<=y<=cell[3]:continue
                    self.assertLessEqual(sum(w*v for w,v in zip(f,values)),upper)

    def test_reference_normalization_is_per_packet_not_per_bit(self):
        ctx.prec=192
        components=[('zero',Q(1),Q(1),Q(0),Q(0),0),
                    ('single',Q(1),Q(0),Q(1),Q(0),1),
                    ('double',Q(1),Q(0),Q(0),Q(1),1)]
        model_=CoverModel(components,None,1,[1,Q(1,2),Q(1,3)],threshold=0)
        cell=(Q(1,4),Q(1,4),Q(1,4),Q(1,4));scale=sum(model_.upper_weights(cell))
        lam=Q(1,100)
        with patch('mixture_cover.iid_kernel.outward',return_value=arb_mat([[1,0],[0,1]])):
            bound=model_.outward(cell,[lam,0,0,0])
        expected=(G*L*aq(scale).log()+L*aq(model_.cell_loss(cell)).log()
                  +G*sum(map(aq,model_.coefficients),arb(0)).log()).exp()
        self.assertTrue(bound>=expected)
        self.assertTrue(bound<expected*(1+arb(2)**-150))

    def test_saved_threshold_is_explicit_and_cannot_be_silently_raised(self):
        self.assertEqual(resolve_threshold(None),209715)
        self.assertEqual(resolve_threshold(None,{}),209715)
        self.assertEqual(resolve_threshold(None,{'threshold':167772}),167772)
        self.assertEqual(resolve_threshold(167772,{},resume=True),167772)
        self.assertEqual(resolve_threshold(167772,{'threshold':167772}),167772)
        for value in (-1,2*G*L,1.5,True):
            with self.assertRaises(ValueError):resolve_threshold(value)
        with self.assertRaises(ValueError):resolve_threshold(167772,{})
        with self.assertRaises(ValueError):resolve_threshold(209715,{'threshold':167772},resume=True)
        self.assertEqual(resolve_threshold(209715,{'threshold':167772},retarget=True),209715)
        with self.assertRaises(ValueError):resolve_threshold(None,{'threshold':167772},retarget=True)
        with self.assertRaises(ValueError):resolve_threshold(209715,retarget=True)

    def test_lower_threshold_improves_a_fixed_outward_witness(self):
        ctx.prec=192
        components=[('zero',Q(1),Q(1),Q(0),Q(0),0),
                    ('single',Q(1),Q(0),Q(1),Q(0),1),
                    ('double',Q(1),Q(0),Q(0),Q(1),1)]
        high=CoverModel(components,None,1,[1,Q(1,2),Q(1,3)],threshold=209715)
        low=CoverModel(components,None,1,[1,Q(1,2),Q(1,3)],threshold=167772)
        cell=(Q(1,4),)*4;lam=Q(1,100)
        with patch('mixture_cover.iid_kernel.outward',return_value=arb_mat([[1,0],[0,1]])):
            old=high.outward(cell,[lam,0,0,0]);new=low.outward(cell,[lam,0,0,0])
        self.assertTrue(new<old)
        expected=aq(lam*(high.threshold-low.threshold))
        self.assertLess(abs(float(old.log()-new.log()-expected)),1e-40)

    def test_exact_projection_pruning_retains_boundaries(self):
        polygon=convex_hull([(Q(1,9),Q(0)),(Q(1),Q(0)),(Q(0),Q(1)),(Q(1,2),Q(1,2))])
        self.assertEqual(len(polygon),3)
        self.assertTrue(separated((Q(0),Q(1,100),Q(0),Q(1,100)),polygon))
        self.assertFalse(separated((Q(1,9),Q(1,9),Q(0),Q(0)),polygon))
        self.assertFalse(separated((Q(1,4),Q(1,2),Q(1,4),Q(1,2)),polygon))

    def test_multinomial_dual_bounds_whole_cell(self):
        ctx.prec=192
        coeff=[Q(1),Q(3,2),Q(5,3),Q(2)]
        features=[(Q(1),Q(0),Q(0)),(Q(1,4),Q(1,2),Q(1,4)),
                  (Q(0),Q(1,2),Q(1,2)),(Q(0),Q(0),Q(1))]
        active=[0,1,1,1];n=6;q_min=2
        for cell in ((Q(0),Q(1),Q(0),Q(1)),(Q(1,8),Q(3,8),Q(1,8),Q(3,4))):
            actual=Q(0)
            for a in range(n+1):
                for b in range(n-a+1):
                    for c in range(n-a-b+1):
                        counts=[a,b,c,n-a-b-c]
                        if n-a<q_min:continue
                        x=sum(k*f[1] for k,f in zip(counts,features))/n
                        y=sum(k*f[2] for k,f in zip(counts,features))/n
                        if not cell[0]<=x<=cell[1] or not cell[2]<=y<=cell[3]:continue
                        term=Q(factorial(n))
                        for k,w in zip(counts,coeff):term*=w**k/factorial(k)
                        actual+=term
            for dual in ((0,0,0),(2,-3,1),(-5,7,3)):
                upper=outer_log_bound(coeff,features,active,n,q_min,cell,dual).exp()
                self.assertTrue(upper>=arb(actual.numerator)/actual.denominator)

    def test_partition_replay_and_rejection(self):
        class Toy:
            def empty(self,cell):return cell[0]>=Q(1,2)
            def proposal(self,cell):return (-20 if cell[1]-cell[0]>Q(1,2) else -100),[1,0,0,0]
            def outward(self,cell,witness):return arb(2)**-100
        model=Toy()
        with patch('builtins.print'):record=cover(model,max_cells=10)
        self.assertFalse(record['unresolved'])
        self.assertEqual(replay(model,record),arb(2)**-100)
        original=record['leaves'];record['leaves']={'0':original['0']}
        with self.assertRaises(ValueError):replay(model,record)
        record['leaves']={'0':original['0'],'1':original['1'],'':original['0']}
        with self.assertRaises(ValueError):replay(model,record)
        self.assertEqual(children((Q(0),Q(1),Q(0),Q(1)))[0],(Q(0),Q(1,2),Q(0),Q(1)))

    def test_retarget_requeues_failed_old_leaves_without_losing_coverage(self):
        class Toy:
            def empty(self,cell):return False
            def proposal(self,cell):return -100,['new']
            def outward(self,cell,witness):return arb(2)**(-100 if witness==['new'] or cell[0]>=Q(1,2) else -10)
        instance=Toy()
        old={'leaves':{'0':{'witness':['old']},'1':{'witness':['old']}},'unresolved':{},'upper_dyadic':[0,0]}
        with patch('builtins.print'),patch.object(instance,'proposal',wraps=instance.proposal) as propose:
            result=cover(instance,max_cells=1,resume=old,retarget=True)
        self.assertEqual(propose.call_count,1)
        self.assertFalse(result['unresolved'])
        self.assertEqual(result['leaves']['0']['witness'],['new'])
        self.assertEqual(result['leaves']['1']['witness'],['old'])
        self.assertEqual(replay(instance,result),arb(2)**-99)
        self.assertEqual(old['leaves']['0']['witness'],['old'])
        with patch('builtins.print'):partial=cover(instance,max_cells=0,resume=old,retarget=True)
        self.assertEqual(set(partial['unresolved']),{'0'})
        with self.assertRaises(ValueError):replay(instance,partial)
        with self.assertRaises(ValueError):cover(instance,retarget=True)

    def test_resume_reconstructs_partition_and_recomputes_bounds(self):
        class Toy:
            def empty(self,cell):return False
            def proposal(self,cell):return (-20 if cell[1]-cell[0]>Q(1,2) else -100),[1,0,0,0]
            def outward(self,cell,witness):return arb(2)**-100
        instance=Toy()
        with patch('builtins.print'):partial=cover(instance,max_cells=2)
        self.assertEqual(set(partial['leaves']),{'0'})
        self.assertEqual(set(partial['unresolved']),{'1'})
        # Neither the cached sum nor stored unresolved coordinates are evidence.
        partial['upper_dyadic']=[0,0];partial['unresolved']['1']['cell']=['0']*4
        with patch('builtins.print'),patch.object(instance,'outward',wraps=instance.outward) as compute:
            result=cover(instance,max_cells=2,resume=partial)
        self.assertEqual(compute.call_count,2)
        self.assertFalse(result['unresolved']);self.assertEqual(replay(instance,result),arb(2)**-99)
        bad={**partial,'unresolved':{}}
        with self.assertRaises(ValueError):cover(instance,resume=bad)
        bad={**partial,'unresolved':{'0':{},'1':{}}}
        with self.assertRaises(ValueError):cover(instance,resume=bad)
        bad={**partial,'leaves':{'0':{'empty':True}}}
        with self.assertRaises(ValueError):cover(instance,resume=bad)
        bad={**partial,'screen_only':True}
        with self.assertRaises(ValueError):cover(instance,resume=bad)


if __name__=='__main__':unittest.main()
