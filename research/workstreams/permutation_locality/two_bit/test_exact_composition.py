from fractions import Fraction as Q
from itertools import product
from unittest.mock import patch
import unittest

from flint import arb,ctx

from exact_composition import boundary_compositions
from affine_mixture import AffineModel,sparse_anchor_candidates
from row_mixture import envelope,pair_components


class ExactCompositionTests(unittest.TestCase):
    def test_sparse_density_anchor_proposal_includes_disappearance(self):
        self.assertEqual(sparse_anchor_candidates((3891,204,1)),((3891,204,0),(3891,204,1)))
        self.assertEqual(sparse_anchor_candidates((4096,0,0)),((4096,0,0),))

    def test_complete_lists_match_every_small_composition(self):
        features=[(Q(0),)*3,(Q(1,7),Q(0),Q(0)),(Q(3,13),Q(1,13),Q(0)),
                  (Q(1,5),Q(0),Q(1,2)),(Q(0),Q(1),Q(0))]
        active=[0,1,1,1,1];slots=4;q_min=2;reference=[]
        for counts in product(range(slots+1),repeat=len(features)):
            if sum(counts)!=slots or sum(n*b for n,b in zip(counts,active))<q_min:continue
            mean=tuple(sum(n*f[j] for n,f in zip(counts,features))/slots for j in range(3))
            reference.append((counts,mean))
        cells=[(Q(0),Q(1))*3]
        cells.extend(tuple(v for x in p for v in (x,x)) for _,p in reference)
        cells.extend((Q(i,28),Q(i+1,28),Q(j,52),Q(j+1,52),Q(0),Q(1,10))
                     for i in range(7) for j in range(6))
        for cell in cells:
            expected={n for n,p in reference if all(lo<=x<=hi for x,lo,hi in zip(p,cell[::2],cell[1::2]))}
            actual=boundary_compositions(features,active,cell,slots,q_min,list_limit=1000,node_limit=100000)
            self.assertIsNotNone(actual)
            self.assertEqual(set(actual),expected)
            self.assertEqual(len(actual),len(expected))
        self.assertIsNone(boundary_compositions(features,active,cells[0],slots,q_min,list_limit=1))
        self.assertIsNone(boundary_compositions(features,active,cells[0],slots,q_min,node_limit=1))

    def test_production_obstruction_and_replay_gaps(self):
        ctx.prec=192
        components=pair_components(envelope([1,2,3,2,1],4,Q(2,5)))
        instance=AffineModel(components,None,401,[1,Q(1,4),Q(1,4)])
        raw=(Q(3677,262144),Q(1839,131072),Q(15,262144),Q(1,16384),Q(0),Q(1,262144))
        cell=instance.integer_cell(raw)
        expected=tuple({'0,0':3695,'0,3':400,'4,4':1}.get(c[0],0) for c in components)
        counts=boundary_compositions(instance.features,instance.active,cell,4096,401)
        self.assertEqual(counts,(expected,))
        row={'counts':list(expected),'parameters':['1/10','1/4','1/2'],'capped':True}
        witness={'compositions':[row]}
        with patch('affine_mixture.mixture_probe.outward',return_value=arb(2)**-100):
            self.assertEqual(instance.outward(raw,witness),arb(2)**-100)
        with self.assertRaises(ValueError):instance.outward(raw,{'compositions':[]})
        with self.assertRaises(ValueError):instance.outward(raw,{'compositions':[row,row]})
        extra={**row,'counts':[4096]+[0]*(len(components)-1)}
        self.assertEqual(instance.composition_rows(cell,{'compositions':[row,extra]}),[row])

    def test_single_central_row_does_not_exhaust_the_enumeration_limit(self):
        components=pair_components(envelope([1,2,3,2,1],4,Q(2,5)))
        instance=AffineModel(components,None,401,[1,Q(1,4),Q(1,4)])
        raw=(Q(1837,131072),Q(3675,262144),Q(31,262144),Q(1,8192),Q(31,262144),Q(1,8192))
        counts=boundary_compositions(instance.features,instance.active,instance.integer_cell(raw),4096,401)
        self.assertEqual(counts,())
        self.assertTrue(instance.empty(raw))

    def test_two_central_rows_at_shifted_boundary_are_infeasible(self):
        components=pair_components(envelope([1,2,3,2,1],4,Q(2,5)))
        instance=AffineModel(components,None,409,[1,Q(1,4),Q(1,4)])
        raw=tuple(map(Q,['1887/131072','3775/262144','39/262144','5/32768','63/262144','1/4096']))
        cell=instance.integer_cell(raw)
        self.assertIsNone(boundary_compositions(instance.features,instance.active,cell,4096,409,term_limit=10))
        self.assertEqual(boundary_compositions(instance.features,instance.active,cell,4096,409),())
        self.assertTrue(instance.empty(raw))


if __name__=='__main__':unittest.main()
