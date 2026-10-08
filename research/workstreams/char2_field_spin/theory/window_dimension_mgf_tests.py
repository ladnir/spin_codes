"""Exact small checks of NA moment domination and convex smoothing."""

from fractions import Fraction as F
from itertools import product
import unittest

from window_dimension_mgf import moment, profile_moment, exact_fixed_profile_moment


class WindowMGFTests(unittest.TestCase):
    def test_exhaustive_subset_moment_domination(self):
        for lengths in product((1,2), repeat=3):
            for parity in (1,2):
                for mark in (F(3,2), F(2), F(4)):
                    actual = exact_fixed_profile_moment(3, lengths, parity, mark)
                    marginal = profile_moment(parity, [F(w,3) for w in lengths], mark)
                    self.assertLessEqual(actual, marginal**3)

    def test_convex_smoothing(self):
        for numerator in product(range(5), repeat=3):
            probabilities = [F(n,4) for n in numerator]
            for parity in (1,2):
                for mark in (F(3,2), F(3)):
                    actual = profile_moment(parity, probabilities, mark)
                    smooth = moment(3, parity, sum(probabilities)/3, mark)
                    self.assertLessEqual(actual, smooth)

    def test_endpoints_and_normalization(self):
        for parity in range(4):
            self.assertEqual(moment(4,parity,F(0),F(3)),1)
            self.assertEqual(moment(4,parity,F(1),F(3)),3**(4-parity))
            self.assertEqual(moment(4,parity,F(2,3),F(1)),1)

    def test_window_length_count(self):
        import math
        for used in range(1,5):
            for budget in range(used,8):
                actual = sum(sum(lengths)<=budget for lengths in product(range(1,budget+1), repeat=used))
                self.assertEqual(actual, math.comb(budget,used))


if __name__ == '__main__':
    unittest.main()
