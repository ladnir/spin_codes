from fractions import Fraction as Q
from itertools import product
import unittest

from integer_projection import possibly_reachable,integer_projection_empty,boundary_composition_empty


class IntegerProjectionTests(unittest.TestCase):
    def test_boundary_enumeration_agrees_with_small_complete_compositions(self):
        features=[(Q(0),)*3,(Q(1,7),Q(0),Q(0)),(Q(3,13),Q(1,13),Q(0)),
                  (Q(1,5),Q(0),Q(1,2)),(Q(0),Q(1),Q(0))]
        active=[0,1,1,1,1];slots=4;q_min=2;points=set()
        for counts in product(range(slots+1),repeat=len(features)):
            if sum(counts)!=slots or sum(n*b for n,b in zip(counts,active))<q_min:continue
            points.add(tuple(sum(n*f[j] for n,f in zip(counts,features))/slots for j in range(3)))
        cells=[tuple(v for x in p for v in (x,x)) for p in points]
        cells.extend((Q(i,28),Q(i+1,28),Q(j,52),Q(j+1,52),Q(0),Q(1,10))
                     for i in range(7) for j in range(6))
        for cell in cells:
            expected=not any(all(lo<=x<=hi for x,lo,hi in zip(p,cell[::2],cell[1::2])) for p in points)
            self.assertEqual(boundary_composition_empty(features,active,cell,slots,q_min),expected)
        cell=(Q(1,7),Q(3,14),Q(1,20),Q(3,40),Q(0),Q(0))
        self.assertFalse(boundary_composition_empty(features,active,cell,slots,q_min,node_limit=1))

    def test_reachability_matches_small_exhaustive_reference(self):
        values=(Q(1,7),Q(2,5),Q(2,3));slots=4
        sums={sum(k*v for k,v in zip(counts,values))
              for counts in product(range(slots+1),repeat=len(values)) if sum(counts)<=slots}
        for i in range(49):
            for width in (0,1,3):
                lo,hi=Q(i,21),Q(i+width,21)
                self.assertEqual(possibly_reachable(values,lo,hi,slots),any(lo<=v<=hi for v in sums))
        self.assertTrue(possibly_reachable((Q(2),Q(5)),Q(3),Q(3),1,node_limit=1))
        self.assertTrue(possibly_reachable((Q(1,100),),Q(199,100),Q(2),1000))

    def test_projection_never_discards_small_integer_compositions(self):
        features=[(Q(0),Q(0),Q(0)),(Q(1,7),Q(0),Q(0)),
                  (Q(3,13),Q(1,13),Q(0)),(Q(1,2),Q(1,2),Q(1,2))]
        slots=5
        for counts in product(range(slots+1),repeat=len(features)):
            if sum(counts)!=slots:continue
            mean=[sum(n*f[j] for n,f in zip(counts,features))/slots for j in range(3)]
            for epsilon in (Q(0),Q(1,10000)):
                cell=tuple(v for x in mean for v in (max(Q(0),x-epsilon),x+epsilon))
                self.assertFalse(integer_projection_empty(features,cell,slots))

    def test_failed_production_cell_has_no_possible_double_count(self):
        values=(Q(1),Q(1,2),Q(2,5),Q(3,5),Q(1,7),Q(2,19),Q(3,16),Q(1,13),Q(6,43),Q(9,37))
        # This is already infeasible without excluding the central components.
        self.assertFalse(possibly_reachable(values,Q(5,32),Q(41,256),4096))
        features=[(Q(0),Q(0),Q(0)),(Q(1,7),Q(0),Q(0)),
                  (Q(3,13),Q(1,13),Q(0)),(Q(5,19),Q(2,19),Q(1,2)),
                  (Q(13,43),Q(6,43),Q(0)),(Q(12,37),Q(9,37),Q(0))]
        cell=(Q(14701,1048576),Q(7351,524288),Q(5,131072),Q(41,1048576),Q(0),Q(1,1048576))
        self.assertTrue(integer_projection_empty(features,cell,4096))

    def test_joint_boundary_rejects_a_cell_that_passes_each_projection(self):
        from lifted_mixture import LiftedModel
        from row_mixture import envelope,pair_components
        components=pair_components(envelope([1,2,3,2,1],4,Q(2,5)))
        instance=LiftedModel(components,None,401,[1,Q(1,4),Q(1,4)])
        raw=(Q(459,32768),Q(3673,262144),Q(7,131072),Q(15,262144),Q(0),Q(1,262144))
        cell=instance.integer_cell(raw)
        self.assertFalse(integer_projection_empty(instance.features,cell,4096))
        self.assertTrue(boundary_composition_empty(instance.features,instance.active,cell,4096,401))
        self.assertTrue(instance.empty(raw))

    def test_next_obstruction_really_contains_an_integer_comparison(self):
        from lifted_mixture import LiftedModel
        from row_mixture import envelope,pair_components
        components=pair_components(envelope([1,2,3,2,1],4,Q(2,5)))
        instance=LiftedModel(components,None,401,[1,Q(1,4),Q(1,4)])
        cell=(Q(3677,262144),Q(1839,131072),Q(15,262144),Q(1,16384),Q(0),Q(1,262144))
        point=((Q(400,7)+Q(12,37))/4096,Q(9,37)/4096,Q(0))
        self.assertTrue(all(lo<=x<=hi for x,lo,hi in zip(point,cell[::2],cell[1::2])))
        self.assertFalse(instance.empty(cell))


if __name__=='__main__':unittest.main()
