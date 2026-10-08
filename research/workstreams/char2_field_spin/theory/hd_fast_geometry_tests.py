"""Small exact checks and recorded high-distance obstruction regressions."""
from fractions import Fraction as F
from itertools import combinations, product
import unittest

from hd_fast_screens import window_gate
from high_distance_window_obstruction import lower_bound, mean_row_dimension


class GeometryTests(unittest.TestCase):
    def test_deterministic_dimension_inequality_exhaustive(self):
        # Every possible family of fixed-size regional subsets, not samples.
        b,k,rows,window = 4,3,4,2
        sets = tuple(combinations(range(rows),window))
        for regions in product(sets,repeat=b):
            hits = [sum(i in region for region in regions) for i in range(rows)]
            dimension = sum(max(0,h-(b-k)) for h in hits)
            self.assertGreaterEqual(dimension,b*window-(b-k)*rows)
            for excluded in range(rows+1):
                retained = sum(max(0,h-(b-k)) for h in hits[:rows-excluded])
                self.assertGreaterEqual(retained,dimension-k*excluded)

    def test_excluded_row_moment(self):
        b,k,rows,window,excluded = 3,2,3,1,1
        choices = tuple(combinations(range(rows),window))
        total = 0
        for regions in product(choices,repeat=b):
            total += sum(max(0,sum(i in r for r in regions)-(b-k))
                         for i in range(rows-excluded))
        self.assertEqual(F(total,len(choices)**b),
                         (rows-excluded)*mean_row_dimension(b,k,F(window,rows)))

    def test_expanding_high_rate_failure_is_deterministic(self):
        for args in ((192,128,4,4,(5,5,6)),(256,192,4,4,(6,))):
            gate = window_gate(*args)
            self.assertEqual(gate['failure_probability_lower_bound'],'1')
            self.assertGreater(gate['deterministic_kernel_dimension_lower_bound'],90000)
            self.assertLess(F(gate['window'],gate['rows']),F(2,5))

    def test_rectangular_nine_output_failure_survives_shortening(self):
        gate = window_gate(256,144,8,8,(9,),18)
        self.assertEqual((gate['rows'],gate['window']),(7296,2912))
        self.assertGreater(gate['failure_probability_lower_bound_float'],.002)
        self.assertLess(gate['failure_probability_lower_bound_float'],.0021)
        self.assertGreater(gate['expected_shortened_dimension'],gate['exit_capacity'])

    def test_next_rectangular_gate_is_only_inconclusive(self):
        gate = window_gate(256,136,16,8,(17,),0)
        self.assertEqual(gate['failure_probability_lower_bound'],'0')
        self.assertEqual(gate['deterministic_kernel_dimension_lower_bound'],0)


if __name__ == '__main__':
    unittest.main()
