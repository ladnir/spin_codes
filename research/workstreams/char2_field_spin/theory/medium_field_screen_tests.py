"""Exact local regressions and independent finite-composition checks."""
from fractions import Fraction as F
from math import comb, exp, log
import unittest

import numpy as np

import medium_field_screen as screen
from superregular_return_screen import local_matrix
from symmetric_density_screen import return_bound, bound as old_density_bound
from state_capacity_reference import matrices, mul


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


class MediumFieldTests(unittest.TestCase):
    def test_exact_q256_packet_regression(self):
        for state in (2, 4, 6, 8):
            for y, tau in ((F(1, 3), F(3)), (F(7, 8), F(19)), (F(1), F(255))):
                for w in range(17):
                    expected = local_matrix(w, y, tau, state)
                    if w:
                        expected[1][0] = return_bound(w, y, tau, state)
                    self.assertEqual(screen.packet(w, y, tau, 256, state), expected)

    def test_exact_large_field_factors_and_domain(self):
        q = 2**32
        for state in (2, 4, 8):
            factors = screen.root_factors(F(q-1), q, state)
            self.assertEqual(factors[0], 1)
            self.assertTrue(all(isinstance(v, F) for v in factors))
            self.assertTrue(all(0 < v < 1 for v in factors[1:]))
            for z in range(1, state):
                self.assertLessEqual(factors[z], F(q, q-1)**state/F(q-1)**z)
        for q_bad in (16, 2**16, 2**128):
            with self.assertRaises(ValueError):
                screen.packet(1, F(1, 2), F(2), q_bad)
        with self.assertRaises(ValueError):
            screen.packet(1, F(1, 2), F(256), 256)
        screen.packet(1, F(1, 2), F(256), q)

    def test_same_tilt_larger_field_only_changes_root_factor(self):
        for state in (2, 4, 8):
            for w in (0, 1, 4, 8, 16):
                small = screen.packet(w, F(3, 4), F(19), 256, state)
                large = screen.packet(w, F(3, 4), F(19), 2**32, state)
                self.assertEqual(small[0], large[0])
                self.assertEqual(small[1][0], large[1][0])
                self.assertLessEqual(large[1][1], small[1][1])

    def test_fraction_region_against_logarithmic_convolution(self):
        rows, active, state, q = 32, 3, 4, 2**32
        y, tau = F(9, 10), F(1000)
        coefficients = [[[F(1), F(0)], [F(0), F(1)]]]
        for _ in range(rows//16):
            following = [[[F(0), F(0)], [F(0), F(0)]] for _ in range(active+1)]
            for h, previous in enumerate(coefficients):
                for w in range(active+1-h):
                    contribution = multiply(previous, screen.packet(w, y, tau, q, state))
                    for i in range(2):
                        for j in range(2):
                            following[h+w][i][j] += comb(16, w)*contribution[i][j]
            coefficients = following
        actual = np.exp(screen.log_coefficients(active, float(y), float(tau), rows, q, state))
        np.testing.assert_allclose(actual, np.asarray(coefficients, dtype=float), rtol=2e-13)
        regional = [[sum(F(comb(active, h), comb(rows, h))*tau**h*coefficients[h][i][j]
                         for h in range(active+1)) for j in range(2)] for i in range(2)]
        logs = screen.regional_from_coefficients(np.log(np.maximum(actual, 1e-300)),
                                                active, rows, log(float(tau)))
        np.testing.assert_allclose(np.exp(logs), np.asarray(regional, dtype=float), rtol=2e-13)

    def test_global_density_q256_regression(self):
        point = (-.7, -.3, log(11), .2)
        for state in (4, 6, 8):
            expected = old_density_bound(point, 10000, 16, 64, 48, state, joint=True)
            actual = screen.density_bound(point, 10000, 256, state, 64, 48, 16)
            self.assertAlmostEqual(actual, expected, places=8)

    def test_active_formula_against_direct_matrix_product(self):
        point, active, q, state = (2., log(1000)), 2, 2**32, 4
        g = screen.geometry(6, 8, 6)
        logy = -exp(point[0])*active/g['n']
        coefficients = screen.log_coefficients(active, exp(logy), exp(point[1]), g['rows'], q, state)
        regional = np.exp(screen.regional_from_coefficients(coefficients, active, g['rows'], point[1]))
        vector = np.array((1., 0.))
        for _ in range(8):
            vector = vector @ regional
        expected = (log(comb(g['rows'], active))-(2*active+1)*point[1]
                    -g['threshold']*logy+log(float(vector.sum())))/log(2)
        self.assertAlmostEqual(screen.active_bound(point, active, q, state, 8, 6, 6), expected, places=10)

    def test_restriction_retains_ca_zero(self):
        for phase in range(4):
            a, c = matrices(phase, True)
            for state in (2, 4, 6, 8):
                for i in range(state):
                    for j in range(state):
                        value = 0
                        for h in range(16):
                            value ^= mul(c[i][h], a[h][j])
                        self.assertEqual(value, 0)


if __name__ == '__main__':
    unittest.main()
