"""Exact labelled-subset controls only; no production operator/gate."""
from fractions import Fraction as F
from itertools import product
from math import log2
import unittest

import numpy as np

from hd_two_sum_gf16_cross_coverage import (covered_groups, coverage_geometry,
    coverage_transition, subset_weights, unrestricted_groups)


class CrossCoverageTests(unittest.TestCase):
    def test_coverage_against_all_labelled_group_subsets(self):
        for cross in range(4):
            weights = [F(j+2, 3) for j in range(cross+1)]
            for groups in range(4):
                exact = F(0)
                for masks in product(range(1 << cross), repeat=groups):
                    union, mass = 0, F(1)
                    for mask in masks:
                        union |= mask
                        mass *= weights[mask.bit_count()]
                    if union == (1 << cross)-1:
                        exact += mass
                actual = covered_groups([log2(x) for x in weights], groups)
                if exact:
                    self.assertAlmostEqual(actual, log2(exact), places=13)
                else:
                    self.assertEqual(actual, -np.inf)

    def test_fixed_subset_optional_coordinates(self):
        k, cross, live, arity = 2, 2, 1, 3
        weights = [F(w+1) for w in range(6)]
        b0, bx, marker = F(2, 3), F(3, 4), F(5, 2)
        actual = subset_weights([log2(x) for x in weights], k, cross, live,
            log2(b0), log2(bx), log_marker=log2(marker), arity=arity)
        for chosen in range(cross+1):
            exact = F(0)
            for mask in range(1 << (k+live)):
                zero_count = (mask & ((1 << k)-1)).bit_count()
                live_count = (mask >> k).bit_count()
                weight = zero_count+live_count+chosen
                exact += weights[weight]*b0**zero_count*bx**chosen*marker**min(weight, arity)
            self.assertAlmostEqual(actual[chosen], log2(exact), places=13)

    def test_zero_cross_matches_unrestricted_power(self):
        for groups in (0, 1, 3, 18):
            self.assertEqual(covered_groups([2.5], groups), groups*2.5)
            self.assertEqual(unrestricted_groups([2.5], groups), groups*2.5)

    def test_coverage_is_positive_subsum_not_subtracted_upper_bounds(self):
        values = [-1., -3., -2., -5.]
        for groups in (1, 2, 4):
            self.assertLessEqual(covered_groups(values, groups), unrestricted_groups(values, groups)+1e-12)
        self.assertEqual(covered_groups([0., -np.inf, -np.inf], 4), -np.inf)

    def test_no_extra_label_multiplicity(self):
        # Uniform subset weights: each of c coordinates appears at least
        # once among r groups in (2^r−1) ways, independently.
        for cross, groups in ((1, 2), (3, 2), (4, 3)):
            self.assertAlmostEqual(covered_groups([0.]*(cross+1), groups),
                                   cross*log2(2**groups-1), places=13)

    def test_guards_and_cached_geometry(self):
        self.assertIs(coverage_geometry(3), coverage_geometry(3))
        self.assertTrue(np.all(np.isneginf(np.tril(coverage_transition([0., 1., 2.]), -1)
            [np.tril_indices(3, -1)])))
        for invalid in (-1, 17, True):
            with self.assertRaises(ValueError):
                coverage_geometry(invalid)
        with self.assertRaises(ValueError):
            subset_weights([0.]*5, 2, 2, 1, -1., -1., arity=2)


if __name__ == '__main__':
    unittest.main()
