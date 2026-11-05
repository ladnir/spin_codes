import unittest
from math import comb

import numpy as np

from hd_shifted_good_pair_screen import (
    good_geometry, good_pair_coefficients, good_step_coefficients,
    low_impurity_frequency, matrix_product, paired_frequency, square_good_coefficients,
)
from hd_shifted_general_screen import general_coefficients


class GoodPairScreenTests(unittest.TestCase):
    def test_exact_good_geometry_marginals(self):
        impurity, full, counts, _ = good_geometry()
        self.assertEqual(len(counts), 425)
        self.assertEqual(int(sum(impurity == 0)), 153)
        self.assertEqual(int(sum(impurity == 1)), 272)
        for degree in range(33):
            self.assertEqual(int(sum(counts[(impurity == 0) & (full == degree)])), comb(32, degree))
            self.assertEqual(int(sum(counts[(impurity == 1) & (full == degree)])),
                             32*comb(31, degree) if degree <= 31 else 0)

    def test_alpha_one_equals_ordered_coefficient_product(self):
        paired = good_pair_coefficients(.5, 1.)
        reference = square_good_coefficients(good_step_coefficients(.5, 1.))
        np.testing.assert_allclose(paired, reference, rtol=2e-12, atol=1e-280)

    def test_good_coefficients_match_full_oriented_polynomial(self):
        odds = (-2., -3., -1.)
        full = general_coefficients(.5, .1, odds)
        good = good_step_coefficients(.5, .1).reshape((2, 33, 4))
        full_marker = np.exp2(np.arange(33)*odds[2])[:, None]
        np.testing.assert_allclose(full[:, 0, 0, :].T, good[0]*full_marker,
                                   rtol=2e-12, atol=1e-280)
        np.testing.assert_allclose(full[:, 1, 0, :].T, good[1]*full_marker*np.exp2(odds[0]),
                                   rtol=2e-12, atol=1e-280)
        np.testing.assert_allclose(full[:, 0, 1, :].T, good[1]*full_marker*np.exp2(odds[1]),
                                   rtol=2e-12, atol=1e-280)

    def test_before_alpha_is_positive_and_improves_bound(self):
        paired = good_pair_coefficients(.5, .1)
        reference = square_good_coefficients(good_step_coefficients(.5, .1))
        self.assertTrue(np.all(paired >= 0))
        self.assertTrue(np.all(paired <= reference*(1+2e-12)+1e-280))
        self.assertTrue(np.any(paired < .9*reference))
        self.assertLess(float(np.sum(paired)), float(np.sum(reference)))

    def test_frequency_decomposition_retains_noncommuting_order(self):
        dirty = np.array((2, 0, 1, 3), dtype=complex)[:, None]
        good = np.array((1, 2, 3, 1), dtype=complex)[:, None]
        self.assertFalse(np.array_equal(matrix_product(good, dirty), matrix_product(dirty, good)))
        result = paired_frequency(dirty, good, matrix_product(good, good))
        np.testing.assert_array_equal(result, matrix_product(dirty+good, dirty+good))

    def test_orientation_markers_and_degree_64_fold(self):
        coefficients = np.zeros((3, 65, 2, 2))
        coefficients[2, 64] = ((1, 2), (3, 4))
        frequencies = low_impurity_frequency(coefficients, (-1., -2., 0.), 64, (2, 2, 64))
        roots = np.exp(-2j*np.pi*np.arange(64)/64)
        singleton = .5*roots[:, None]+.25*roots[None, :]
        expected = np.array((1, 2, 3, 4))[:, None, None]*singleton[None, :, :]**2
        np.testing.assert_allclose(frequencies[:, :, :, 0], expected, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(frequencies[:, :, :, 17], expected, rtol=1e-14, atol=1e-14)
        truncated = low_impurity_frequency(coefficients, (-1., -2., 0.), 64, (0, 2, 64))
        expected = np.array((1, 2, 3, 4))[:, None]*.25**2*roots[None, :]**2
        np.testing.assert_allclose(truncated[:, 7, :, 0], expected, rtol=1e-14, atol=1e-14)


if __name__ == '__main__':
    unittest.main()
