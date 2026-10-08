"""Exact combinatorial checks for the grouped-cut diagnostic, not a certificate."""
import itertools
import unittest

import numpy as np

from hd_affine_cut import grouped_safe_counts
from hd_grouped_cut_search import (change_delta, finite_score, flip_delta,
                                   local_interior_options, parity_endpoint_seed,
                                   weighted_degrees)
from hd_two_state_regional_strides import ROWS


class GroupedCutTests(unittest.TestCase):
    def test_all_small_flip_and_swap_deltas(self):
        strides = (1, 3, 5)
        for word in range(256):
            selected = np.array([(word >> i) & 1 for i in range(8)], dtype=bool)
            original = selected.copy()
            before = np.array(grouped_safe_counts(selected, 4, strides))
            for row in range(8):
                delta = flip_delta(selected, row, strides)
                selected[row] = ~selected[row]
                self.assertTrue(np.array_equal(before+delta,
                    grouped_safe_counts(selected, 4, strides)))
                selected[row] = ~selected[row]
            for remove in np.flatnonzero(selected):
                for add in np.flatnonzero(~selected):
                    delta = change_delta(selected, (remove, add), strides)
                    self.assertTrue(np.array_equal(selected, original))
                    altered = selected.copy()
                    altered[[remove, add]] = ~altered[[remove, add]]
                    self.assertTrue(np.array_equal(before+delta,
                        grouped_safe_counts(altered, 4, strides)))

    def test_weighted_marginals(self):
        strides = (1, 3, 5)
        weights = np.array((2, 3, 7))
        for word in range(256):
            selected = np.array([(word >> i) & 1 for i in range(8)], dtype=bool)
            degrees = weighted_degrees(selected, weights, strides)
            for row in range(8):
                expected = np.dot(weights, flip_delta(selected, row, strides))
                self.assertEqual(expected, (1-2*int(selected[row]))*degrees[row])

    def test_interior_reduction(self):
        expected = {(0, 0): (0, 0, 1), (0, 1): (0, 1, 2),
                    (1, 0): (0, 1, 2), (1, 1): (0, 1, 3)}
        for endpoints, values in expected.items():
            self.assertEqual(local_interior_options(*endpoints), values)

    def test_interval_dp_against_all_subsets(self):
        n = 6
        for bits in itertools.product((0, 1), repeat=n):
            for numerator, denominator in ((1, 2), (1, 12), (3, 1)):
                best = 0
                for mask in range(1 << n):
                    chosen = [(mask >> i) & 1 for i in range(n)]
                    intervals = sum(chosen[i] and (i == 0 or not chosen[i-1])
                                    for i in range(n))
                    safe = sum(a*b for a, b in zip(bits, chosen))
                    unsafe = sum(chosen)-safe
                    gain, cost = safe-2*intervals, unsafe+2*intervals
                    best = max(best, denominator*gain-numerator*cost)
                self.assertEqual(finite_score(bits, numerator, denominator), best)

    def test_full_group_single_swap_barrier(self):
        strides = (1, 3, 5)
        weights = np.array((2, 3, 7))
        for word in range(1, 15):
            selected = np.repeat([(word >> i) & 1 for i in range(4)], 4).astype(bool)
            def score(bits):
                return np.dot(weights, grouped_safe_counts(bits, 4, strides))
            before = score(selected)
            for remove in np.flatnonzero(selected):
                for add in np.flatnonzero(~selected):
                    altered = selected.copy()
                    altered[[remove, add]] = ~altered[[remove, add]]
                    self.assertLessEqual(score(altered), before)

    def test_opposite_parity_endpoint_identity(self):
        for rows in (ROWS//4, 9830, 13107, ROWS//2, 19660, 22937, 3*ROWS//4):
            selected = parity_endpoint_seed(rows)
            self.assertEqual(int(selected.sum()), rows)
            self.assertEqual(grouped_safe_counts(selected, 4), (rows-ROWS//8,)*64)


if __name__ == "__main__":
    unittest.main()
