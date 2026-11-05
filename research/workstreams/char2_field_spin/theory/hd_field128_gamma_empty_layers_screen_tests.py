"""Tiny exact polynomial comparisons; no K20 gate or alias bank."""

from fractions import Fraction as F
from math import comb, log2
import unittest

import numpy as np

import hd_field128_gamma_empty_layers_screen as mode
from hd_field128_gamma_empty_boundary_screen_tests import polynomial_product


class EmptyLayersTests(unittest.TestCase):
    def example(self):
        poly = {(s, f): (s+1, f+2, s+f+3, 2*s+f+1)
                for s in range(3) for f in range(3-s)}
        logs = np.full((4, 3, 3), -np.inf)
        for (s, f), row in poly.items():
            logs[:, s, f] = np.log2(row)
        return poly, logs

    def test_transformed_polynomial_and_ordered_power(self):
        poly, logs = self.example()
        whole = polynomial_product(poly, poly)
        for marker in (F(1, 2), F(1), F(2)):
            layers = mode.empty_layer_logs(logs, log2(marker), slots=1)
            for empty in range(3):
                expected = [sum(F(poly[2-empty-f, f][entry])*marker**f
                                for f in range(3-empty)) for entry in range(4)]
                np.testing.assert_allclose(np.exp2(layers[:, empty]), list(map(float, expected)), rtol=2e-14)
            for empty in range(5):
                actual, diag = mode.boundary.coefficient_logs(layers, 2, empty, log_marker=1.)
                expected = [sum(F(row[entry])*marker**f for (s, f), row in whole.items()
                                if s+f == 4-empty) for entry in range(4)]
                np.testing.assert_allclose(np.exp2(actual), list(map(float, expected)), rtol=2e-14)
                self.assertFalse(diag['nonfinite'])

    def test_normalized_mixed_bound_dominates_exact_coefficients(self):
        poly, logs = self.example()
        whole = polynomial_product(poly, poly)
        for empty in range(5):
            for full in range(5-empty):
                single = 4-empty-full
                r = single+2*full
                norm = comb(4, empty)*comb(4-empty, full)
                expected = np.array(whole[single, full], dtype=float)/norm
                for eta in (-2., 0., 2.):
                    estimate, point, diag = mode.regional_estimate(logs, r, full, eta, .5,
                                                                  slots=1, exponent=2)
                    self.assertTrue(np.all(np.exp2(estimate) >= expected*(1-1e-12)))
                    self.assertTrue(np.all(point >= estimate-1e-12))
                    self.assertFalse(diag['directed'])

    def test_empty_zero_layer_equals_boundary_polynomial_evaluation(self):
        unused, logs = self.example()
        top = mode.boundary.top_logs(logs, slots=1)
        for eta in (-3., 0., 3.):
            expected = np.logaddexp2.reduce(top+eta*np.arange(3), axis=1)
            np.testing.assert_allclose(mode.empty_layer_logs(logs, eta, slots=1)[:, 0], expected)

    def test_original_prefix_not_complemented(self):
        prefix, forbidden = mode.original_weight_prefix(F(1, 16), F(1, 4), 5120)
        self.assertEqual(forbidden, 838860)
        self.assertAlmostEqual(prefix, .25*(16384*5120-128+4*838860))
        other, unused = mode.original_weight_prefix(F(1, 16), F(1, 4), 3072)
        self.assertEqual(prefix-other, .25*16384*2048)

    def test_rejects_nontriangular_or_invalid_histograms(self):
        unused, logs = self.example()
        logs[:, 2, 2] = 0.
        with self.assertRaises(ValueError):
            mode.empty_layer_logs(logs, 0., slots=1)
        unused, logs = self.example()
        with self.assertRaises(ValueError):
            mode.regional_estimate(logs, 8, 0, 0., slots=1, exponent=2)
        with self.assertRaises(ValueError):
            mode.empty_layer_logs(logs, float('inf'), slots=1)


if __name__ == '__main__':
    unittest.main()
