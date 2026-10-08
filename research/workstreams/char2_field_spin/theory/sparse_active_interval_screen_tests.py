"""Exact monotonicity and independent replay tests for floating interval bins."""

from fractions import Fraction as F
from math import exp, log, log2
import unittest

import numpy as np

import superregular_return_screen as old
from superregular_return_screen_tests import exact_region, matrix_product
from sparse_active_interval_screen import (ANCHORS, FROZEN_INTERVALS, FixedMarks, anchor_zones,
                                            interpolated_anchors, log_coefficients, marks,
                                            regional_from_coefficients, subdivide)


class SparseActiveIntervals(unittest.TestCase):
    def test_regional_and_global_monotonicity_exact(self):
        y, tau, rows = F(9, 10), F(7), 32
        regions = [exact_region(active, y, tau, rows, 8) for active in range(5)]
        for lower, upper in zip(regions, regions[1:]):
            self.assertTrue(all(lower[i][j] <= upper[i][j] for i in range(2) for j in range(2)))
            lower_cube = matrix_product(matrix_product(lower, lower), lower)
            upper_cube = matrix_product(matrix_product(upper, upper), upper)
            self.assertLessEqual(sum(lower_cube[0]), sum(upper_cube[0]))

    def test_coefficient_reuse_matches_existing_regional_power(self):
        y, tau, rows = .97, 7., 64
        coefficients = log_coefficients(y, tau, rows, 5)
        for active in range(6):
            actual = regional_from_coefficients(coefficients, active, rows, log(tau))
            expected = old.region(active, y, tau, rows, 8)
            np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=3e-13)

    def test_anchor_conversion_preserves_y(self):
        anchor, endpoint, point, n = 16, 21, (7.18, log(7)), 1540096
        original = marks(anchor, point, n)
        moved = marks(endpoint, (point[0]+log(anchor/endpoint), point[1]), n)
        self.assertAlmostEqual(original[0], moved[0], places=15)
        self.assertEqual(original[1], moved[1])
        self.assertNotEqual(original[0], marks(endpoint, point, n)[0])

    def test_interval_bounds_every_integer_and_accounts_union(self):
        fixed = FixedMarks(2, (1.5, log(7)), 5, k_log=4, length=3, dimension=2)
        interval = fixed.interval(1, 5)
        self.assertAlmostEqual(interval.per_class_margin_bits-interval.union_margin_bits, log2(5))
        for active in range(1, 6):
            point = (1.5+log(2/active), log(7))
            individual_log2 = old.bound(active, point, 4, 3, 2, 8)
            self.assertLessEqual(individual_log2, -interval.per_class_margin_bits+1e-11)
        expected_argmax = max(range(1, 6), key=fixed.source.__getitem__)
        self.assertEqual(interval.source_maximizer, expected_argmax)

    def test_partition_is_complete_and_positive(self):
        zones = anchor_zones(1, 20, {1: None, 4: None, 16: None})
        self.assertEqual([a for lo, hi, _ in zones for a in range(lo, hi+1)], list(range(1, 21)))
        class Fake:
            def interval(self, lo, hi):
                from sparse_active_interval_screen import Interval
                return Interval(lo, hi, 1, 200/(hi-lo+1), 200/(hi-lo+1), lo)
        parts = subdivide(Fake(), 1, 20, 100)
        self.assertEqual([a for part in parts for a in range(part.lower, part.upper+1)], list(range(1, 21)))
        self.assertTrue(all(part.union_margin_bits >= 100 for part in parts))

    def test_frozen_fixture_and_interpolation(self):
        expanded = interpolated_anchors(ANCHORS)
        self.assertEqual([a for lo, hi, _ in FROZEN_INTERVALS for a in range(lo, hi+1)], list(range(1, 601)))
        self.assertEqual(len(FROZEN_INTERVALS), 80)
        self.assertTrue(all(anchor in expanded for _, _, anchor in FROZEN_INTERVALS))
        self.assertTrue(all(expanded[anchor] == point for anchor, point in ANCHORS.items()))
        for anchor, point in expanded.items():
            logy, logtau = marks(anchor, point, 1540096)
            self.assertLess(logy, 0)
            self.assertTrue(0 <= logtau <= log(255))


if __name__ == '__main__':
    unittest.main()
