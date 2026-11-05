"""Independent exact trace, cyclic-run, box-count, and rounding checks."""
from fractions import Fraction as F
from itertools import product
import unittest

from flint import arb, ctx

from hd_shifted_general_count import transfer_value as four_marker_value
from hd_shifted_pairs_count import circular_run_count, exact_profiles
from hd_shifted_run_count import (transfer_value, box_count_upper,
                                  box_log2_upper, stationary_markers)


class RunCountTests(unittest.TestCase):
    def test_trace_matches_exhaustive_indexed_binary_cycles(self):
        for pairs in range(1, 6):
            for x, y in ((F(2,3),F(3,5)),(F(3,2),F(5,4)),(F(1),F(1))):
                expected = F(0)
                for bits in product((0,1), repeat=2*pairs):
                    weight = sum(bits)
                    runs = sum(not bits[i-1] and bits[i] for i in range(2*pairs))
                    expected += x**weight*y**runs
                self.assertEqual(transfer_value(pairs,x,y), expected)
                self.assertEqual(transfer_value(pairs,x,y),
                                 four_marker_value(pairs,(x*y,x*y,1/y,1/y)))

    def test_run_formula_matches_profile_projection_including_constants(self):
        for pairs in range(1, 9):
            counts = {}
            for (r,a,b), value in exact_profiles(pairs).items():
                key = r,r-a-b
                counts[key] = counts.get(key,0)+value
            for r in range(2*pairs+1):
                for runs in range(pairs+1):
                    self.assertEqual(counts.get((r,runs),0),
                                     circular_run_count(2*pairs,r,runs))

    def test_box_bound_covers_every_small_box_once(self):
        for pairs in range(1, 4):
            exact = exact_profiles(pairs)
            intervals = lambda n:[(lo,hi) for lo in range(n+1) for hi in range(lo,n+1)]
            for ranges in product(intervals(2*pairs),intervals(pairs),intervals(pairs)):
                count = sum(value for point,value in exact.items()
                            if all(lo<=n<=hi for n,(lo,hi) in zip(point,ranges)))
                for x,y in ((F(2,3),F(3,5)),(F(3,2),F(5,4))):
                    upper = box_count_upper(pairs,ranges,x,y)
                    self.assertGreaterEqual(upper,count,(pairs,ranges,x,y))
                    if not count:
                        self.assertEqual(upper,0)

    def test_stationary_markers_are_exact_and_log_is_directed(self):
        self.assertEqual(stationary_markers(4096,256,56),(F(155,197),F(49,24625)))
        ranges = ((2,5),(0,2),(0,2))
        count = box_count_upper(4,ranges,F(2,3),F(3,5))
        upper = box_log2_upper(4,ranges,F(2,3),F(3,5),precision=256)
        with ctx.workprec(384):
            exact = arb(count).log()/arb(2).log()
            self.assertTrue(upper>=exact)
            self.assertTrue(upper.is_exact())

    def test_invalid_inputs_and_empty_box(self):
        with self.assertRaises(TypeError):
            transfer_value(4,.5,F(1))
        with self.assertRaises(ValueError):
            transfer_value(4,0,1)
        with self.assertRaises(ValueError):
            stationary_markers(4,2,2)
        ranges = ((2,2),(1,1),(1,1))
        self.assertEqual(box_count_upper(4,ranges,F(1),F(1)),0)
        with self.assertRaises(ValueError):
            box_log2_upper(4,ranges,F(1),F(1))


if __name__ == '__main__':
    unittest.main()
