"""Tiny positive-cut checks; no production component diagnostic."""
from itertools import product
from math import log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_pre_toeplitz_zero_component as target


class ZeroComponentTests(unittest.TestCase):
    def test_zeta_against_labelled_subset_enumeration(self):
        for group in (4, 8):
            weights = np.array([-.3*mask.bit_count()+.01*mask for mask in range(1 << group)])
            weights[0] = -np.inf
            actual = target.subset_zeta_log(weights)
            for allowed in range(1 << group):
                expected = np.logaddexp2.reduce([value for mask, value in enumerate(weights) if mask & ~allowed == 0])
                self.assertAlmostEqual(actual[allowed], expected, places=12)

    def test_refined_component_against_direct_group_mask_sum(self):
        group, groups, log_b, c, ell = 4, 2, -.4, .2, .3
        weights = np.array([-.5*mask.bit_count()+.02*mask for mask in range(1 << group)])
        prepared = target.prepare_zero_component(weights, log_b)
        result = target.zero_component(prepared, log2(c), groups, ell, columns=4, minimum=3)
        terms = []
        for masks in product(range(1 << group), repeat=groups):
            term = sum(weights[mask]+ell*bool(mask)-3*ell/4 for mask in masks)
            for bit in range(group):
                occupancy = sum(mask >> bit & 1 for mask in masks)
                term += log_b*occupancy+(log2(c) if occupancy else 0.)
            terms.append(term)
        self.assertAlmostEqual(result['refined_column_log2'], np.logaddexp2.reduce(terms), places=12)
        self.assertLessEqual(result['refined_column_log2'], result['ordinary_column_log2'])

    def test_cut_endpoints_and_empty_mass(self):
        weights = np.array([-.4*mask.bit_count() for mask in range(16)])
        prepared = target.prepare_zero_component(weights, -.3)
        full = target.zero_component(prepared, 0., 2, columns=4, minimum=3)
        self.assertAlmostEqual(full['refined_column_log2'], full['ordinary_column_log2'])
        empty = target.zero_component(prepared, -np.inf, 2, .5, columns=4, minimum=3)
        self.assertAlmostEqual(empty['refined_column_log2'], 2*(weights[0]-3*.5/4))
        self.assertEqual(empty['dominant_allowed_mask'], 0)

    def test_all00_matching_invariance(self):
        group = 8
        weights = np.array([[-.2*w-.05*f for f in range(5)] for w in range(9)])
        results = []
        for matching in (tuple((2*i, 2*i+1) for i in range(4)), tuple((i, i+4) for i in range(4))):
            table = target.schedule.mask_weight_table(weights, matching)
            prepared = target.prepare_zero_component(table, -.3)
            results.append(target.zero_component(prepared, -1.5, 3, columns=8, minimum=5)['refined_column_log2'])
        self.assertAlmostEqual(*results, places=12)

    def test_max_plus_retains_column_boundary_state(self):
        matrices = [np.array([[0., 2.], [3., -.5]]), np.array([[.7, -.2], [1.3, .1]])]
        indices = [np.array([[0, 1], [2, 3]]), np.array([[4, 5], [6, 7]])]
        value, trace, final = target.largest_state_word(list(zip(matrices, indices)))
        candidates = [(matrices[0][0, middle]+matrices[1][middle, last], middle, last)
                      for middle, last in product(range(2), repeat=2)]
        expected, middle, last = max(candidates)
        self.assertEqual(value, expected)
        self.assertEqual(final, last)
        self.assertEqual(trace, [indices[0][0, middle], indices[1][middle, last]])

    def test_guards_and_partial_result(self):
        with self.assertRaises(ValueError):
            target.log_complement(.1)
        with self.assertRaises(ValueError):
            target.subset_zeta_log(np.zeros(15))
        for budget in (True, 0, 31):
            with self.assertRaises(ValueError):
                target.run(budget)
        with patch.object(target, 'load_saved', side_effect=TimeoutError('fixture')):
            result = target.run()
        self.assertFalse(result['probability_lower_bound'])
        self.assertFalse(result['whole_probability_upper_bound'])
        self.assertFalse(result['complete_requested_controls'])


if __name__ == '__main__':
    unittest.main()
