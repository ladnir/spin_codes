"""Regression checks for the sound formulas used by the FLOAT XOR2 screen."""
from itertools import product
from math import comb, log2
import unittest

import numpy as np

from hd_grouped_subfield_screen import coefficients
from hd_xor_column_screen import (column_density_bits, matrix_polynomial_power,
                                  marked_comparison_masks, one_group_region_words,
                                  profile_arrays, thinned_coefficients, thinned_matrices)
from hd_xor_column_probe import uniform_column_support_counts


class XorColumnScreenTests(unittest.TestCase):
    def test_full_space_density(self):
        self.assertAlmostEqual(column_density_bits(8, 1.), 0., places=12)
        self.assertAlmostEqual(column_density_bits(1, 1.), 56., places=12)
        for active in range(1, 9):
            self.assertGreaterEqual(column_density_bits(active, .75), -1e-12)

    def test_marked_kernel_at_one(self):
        for active in (1, 2, 7, 8):
            actual = [(weight, sum(2**value for size, value in terms))
                      for weight, terms in marked_comparison_masks(active)]
            counts = uniform_column_support_counts(256, 8, (1 << active)-1)
            expected = [(mask.bit_count(), count/256**active)
                        for mask, count in enumerate(counts) if count]
            self.assertEqual([item[0] for item in actual], [item[0] for item in expected])
            np.testing.assert_allclose([item[1] for item in actual], [item[1] for item in expected], rtol=1e-14)

    def test_forced_nonzero_prevents_quiet_single_input(self):
        words = one_group_region_words(.999)
        self.assertEqual(words[0][0][0], 0.)
        for bit in range(8):
            self.assertTrue(np.isneginf(words[1 << bit][0][0]))
            self.assertTrue(np.isfinite(words[1 << bit][0][1]))

    def test_no_thinning_matches_control(self):
        actual = thinned_coefficients(.5, .125, 1.)
        control = coefficients(.5, .125, h=4, quiet_inverse=True)
        finite = np.isfinite(actual) & np.isfinite(control)
        self.assertTrue(np.array_equal(np.isfinite(actual), np.isfinite(control)))
        # The control rounds each profile entry upward to 13 mantissa bits.
        self.assertLessEqual(float(np.max(actual[finite]-control[finite])), 1e-10)
        self.assertLess(float(np.max(control[finite]-actual[finite])), .0001)

    def test_thinning_against_direct_profile_average(self):
        occupancy, _ = profile_arrays()
        lookup = {tuple(row): index for index, row in enumerate(occupancy)}
        _, _, control = thinned_matrices(.75, 1.)
        probability = .625
        _, _, thinned = thinned_matrices(.75, probability)
        for upper in ((0, 0, 0, 2), (0, 1, 1, 2), (1, 2, 3, 4)):
            expected = np.zeros((2, 2))
            for lower in product(*(range(value+1) for value in upper)):
                weight = 1.
                for high, low in zip(upper, lower):
                    weight *= comb(high, low)*probability**low*(1-probability)**(high-low)
                expected += weight*control[lookup[tuple(sorted(lower))]]
            actual = thinned[lookup[upper]]
            np.testing.assert_allclose(actual[0], expected[0], rtol=2e-13)
            np.testing.assert_allclose(actual[1, 1], expected[1, 1], rtol=2e-13)
            # min(E[return cap], E[live cap]) can exceed E[min(...)].
            self.assertGreaterEqual(actual[1, 0]*(1+1e-12), expected[1, 0])

    def test_bivariate_matrix_polynomial(self):
        terms = (((0, 0), 1), ((1, 0), 2), ((0, 1), 3), ((1, 1), 5))
        local = [(powers, ((log2(value), -np.inf), (-np.inf, log2(value))))
                 for powers, value in terms]
        actual = matrix_polynomial_power(local, 2, 3)
        expected = np.zeros((4, 4), dtype=int)
        for sequence in product(terms, repeat=3):
            left = sum(term[0][0] for term in sequence)
            right = sum(term[0][1] for term in sequence)
            value = 1
            for term in sequence:
                value *= term[1]
            expected[left, right] += value
        np.testing.assert_allclose(np.exp2(actual[..., 0, 0]), expected, rtol=1e-14)
        self.assertTrue(np.all(np.isneginf(actual[..., 0, 1])))

    def test_noncommuting_matrix_polynomial(self):
        terms = (((0, 0), np.array([[1, 2], [0, 1]], dtype=int)),
                 ((1, 0), np.array([[2, 0], [3, 1]], dtype=int)),
                 ((0, 1), np.array([[1, 1], [1, 0]], dtype=int)),
                 ((1, 1), np.array([[0, 1], [2, 3]], dtype=int)))
        with np.errstate(divide='ignore'):
            local = [(powers, np.log2(matrix)) for powers, matrix in terms]
        actual = matrix_polynomial_power(local, 2, 3)
        expected = np.zeros((4, 4, 2, 2), dtype=int)
        for sequence in product(terms, repeat=3):
            left = sum(term[0][0] for term in sequence)
            right = sum(term[0][1] for term in sequence)
            value = np.eye(2, dtype=int)
            for term in sequence:
                value = value @ term[1]
            expected[left, right] += value
        np.testing.assert_allclose(np.exp2(actual), expected, rtol=1e-14)


if __name__ == '__main__':
    unittest.main()
