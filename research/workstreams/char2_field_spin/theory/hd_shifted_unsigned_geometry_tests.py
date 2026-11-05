"""Exact DP and lane-witness checks for the unsigned support domain."""
from itertools import product
import unittest

from hd_shifted_general_count import feasible_witness
from hd_shifted_pairs_count import exact_profiles, profile_count_upper
from hd_shifted_unsigned_geometry import (unsigned_witness,
    unsigned_is_feasible, unsigned_box_witness, unsigned_box_is_empty,
    unsigned_box_all_feasible)


class UnsignedGeometryTests(unittest.TestCase):
    def test_every_small_profile_matches_exact_cycle_dp(self):
        for pairs in range(1, 13):
            exact = exact_profiles(pairs)
            for r, a, b in product(range(2*pairs+1), range(pairs+1), range(pairs+1)):
                actual = unsigned_witness(pairs, r, a, b)
                self.assertEqual(actual is not None, (r, a, b) in exact,
                                 (pairs, r, a, b))
                self.assertEqual(profile_count_upper(pairs, r, a, b) > 0,
                                 actual is not None, (pairs, r, a, b))
                if actual is not None:
                    even, odd = actual
                    self.assertEqual(even+odd, r)
                    self.assertIsNotNone(feasible_witness(pairs, even, odd, a, b),
                                         (pairs, r, a, b, actual))

    def test_every_small_box_matches_exact_cycle_dp(self):
        for pairs in range(1, 5):
            exact = tuple(exact_profiles(pairs))
            weights = [(lo, hi) for lo in range(2*pairs+1) for hi in range(lo, 2*pairs+1)]
            full = [(lo, hi) for lo in range(pairs+1) for hi in range(lo, pairs+1)]
            for ranges in product(weights, full, full):
                expected = [point for point in exact
                            if all(lo <= n <= hi for n, (lo, hi) in zip(point, ranges))]
                witness = unsigned_box_witness(pairs, ranges)
                self.assertEqual(witness is None, not expected, (pairs, ranges))
                if witness is not None:
                    self.assertIn(witness, expected)
                volume = 1
                for lo, hi in ranges:
                    volume *= hi-lo+1
                self.assertEqual(unsigned_box_all_feasible(pairs, ranges),
                                 len(expected) == volume, (pairs, ranges))

    def test_constants_and_alternating_cycle(self):
        for pairs in (1, 2, 4096):
            self.assertEqual(unsigned_witness(pairs, 0, 0, 0), (0, 0))
            self.assertEqual(unsigned_witness(pairs, 2*pairs, pairs, pairs), (pairs, pairs))
            self.assertEqual(unsigned_witness(pairs, pairs, 0, 0), (pairs, 0))
            for full in (1, pairs-1):
                if 0 < full < pairs:
                    self.assertFalse(unsigned_is_feasible(pairs, 2*full, full, full))

    def test_forbidden_first_weight_advances_once(self):
        self.assertIsNone(unsigned_box_witness(8, ((6, 6), (3, 3), (3, 3))))
        self.assertEqual(unsigned_box_witness(8, ((6, 8), (3, 3), (3, 3))), (7, 3, 3))

    def test_zero_exclusion_and_outside_domain(self):
        self.assertEqual(unsigned_box_witness(8, ((-10, 20), (-4, 9), (-4, 9))), (0, 0, 0))
        self.assertEqual(unsigned_box_witness(8, ((-10, 20), (-4, 9), (-4, 9)),
                                             exclude_zero=True), (1, 0, 0))
        zero = ((0, 0), (0, 0), (0, 0))
        self.assertTrue(unsigned_box_is_empty(8, zero, exclude_zero=True))
        self.assertFalse(unsigned_box_all_feasible(8, zero, exclude_zero=True))
        self.assertTrue(unsigned_box_is_empty(8, ((1, 2), (2, 4), (0, 8))))
        self.assertFalse(unsigned_box_all_feasible(8, ((-1, 0), (0, 0), (0, 0))))

    def test_bad_inputs(self):
        for pairs in (0, -1, 1.5):
            with self.assertRaises(ValueError):
                unsigned_is_feasible(pairs, 0, 0, 0)
        self.assertFalse(unsigned_is_feasible(4, 1.5, 0, 0))
        self.assertFalse(unsigned_is_feasible(4, -1, 0, 0))
        self.assertFalse(unsigned_is_feasible(4, 9, 4, 4))
        for ranges in (((1, 0), (0, 0), (0, 0)), ((0, 1), (0, 0)),
                       ((0, .5), (0, 0), (0, 0))):
            with self.assertRaises(ValueError):
                unsigned_box_witness(4, ranges)


if __name__ == '__main__':
    unittest.main()
