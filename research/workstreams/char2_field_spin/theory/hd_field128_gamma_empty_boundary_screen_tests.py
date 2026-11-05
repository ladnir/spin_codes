"""Tiny independent polynomial identities; no q128 bank or selected-weight gate."""

from math import comb, log2
import unittest

import numpy as np

import hd_field128_gamma_empty_boundary_screen as mode


def multiply(a, b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2)) for i in range(2) for j in range(2))


def polynomial_product(a, b):
    out = {}
    for (s, f), x in a.items():
        for (t, g), y in b.items():
            key = s+t, f+g
            z = multiply(x, y)
            out[key] = tuple(v+w for v, w in zip(out.get(key, (0,)*4), z))
    return out


class BoundaryTests(unittest.TestCase):
    def test_top_degree_matches_independent_ordered_composition(self):
        slots = 2
        matrices = [np.array([[j+1, 2*j+1], [3*j+2, j+3]], dtype=float) for j in range(5)]
        for alpha in (1., .5):
            composed = mode.provider.all_pair_coefficients(np.log2(matrices), alpha, slots=slots)
            top = np.exp2(mode.top_logs(composed, slots=slots))
            expected = np.zeros((4, 5))
            wrong = np.zeros((4, 5))
            for i in range(3):
                for j in range(3):
                    count = comb(2, i)*comb(2, j)
                    expected[:, i+j] += count*(matrices[2+i]@matrices[2+j]).ravel()**alpha
                    wrong[:, i+j] += count*((matrices[2+i]**alpha)@(matrices[2+j]**alpha)).ravel()
            np.testing.assert_allclose(top, expected, rtol=2e-14)
            if alpha != 1:
                self.assertGreater(float(np.max(abs(expected-wrong))), 1)

    def test_boundary_coefficient_equals_full_bivariate_power(self):
        slots = 1
        matrices = [(1, 2, 3, 4), (2, 1, 4, 3), (3, 4, 1, 2)]
        geometry = [(0, 0), (1, 0), (0, 1)]
        paired = {}
        for s, f in geometry:
            for t, g in geometry:
                key = s+t, f+g
                row = multiply(matrices[s+2*f], matrices[t+2*g])
                paired[key] = tuple(a+b for a, b in zip(paired.get(key, (0,)*4), row))
        whole = polynomial_product(paired, paired)
        logs = np.array([[log2(paired[2-f, f][entry]) for f in range(3)] for entry in range(4)])
        for f in range(5):
            expected = whole[4-f, f]
            for tilt in (-3., 0., 2.):
                answer, diagnostic = mode.coefficient_logs(logs, 2, f, log_marker=tilt)
                np.testing.assert_allclose(np.exp2(answer), expected, rtol=2e-14)
                self.assertFalse(diagnostic['nonfinite'])
                marker = mode.positive_marker_logs(logs, 2, f, tilt)
                self.assertTrue(np.all(marker >= answer-1e-12))

    def test_reversal_and_scaling_preserve_noncommuting_power(self):
        logs = np.log2(np.array([[1, 2, 4], [3, 5, 1], [2, 1, 7], [4, 2, 3]], dtype=float))
        for degree in range(7):
            a, unused = mode.coefficient_logs(logs, 3, degree, log_marker=1.25, reverse=False)
            b, diag = mode.coefficient_logs(logs, 3, degree, log_marker=1.25, reverse=True)
            np.testing.assert_allclose(a, b, atol=2e-13)
            self.assertEqual(diag['polynomial_reversed'], degree > 6-degree)

    def test_underflow_and_zero_are_explicit_not_certificates(self):
        logs = np.array([[0., -2000.]]*4)
        answer, diag = mode.coefficient_logs(logs, 1, 1, reverse=False)
        self.assertTrue(diag['input_underflow_terms'])
        self.assertTrue(diag['zero_coefficient_ambiguity'])
        self.assertTrue(diag['estimated_zeros_are_not_proved'])
        self.assertFalse(diag['directed'])
        self.assertTrue(np.all(np.isneginf(answer)))

    def test_guards_and_alpha_one_factorization(self):
        logs = np.zeros((4, 3))
        for exponent, degree, eta in ((0, 1, 0.), (2, 5, 0.), (2, 1, float('nan'))):
            with self.assertRaises(ValueError):
                mode.coefficient_logs(logs, exponent, degree, log_marker=eta)
        with self.assertRaises(ValueError):
            mode.top_logs(np.zeros((4, 3, 3)), slots=2)
        matrices = np.array([[[j+1, 1], [2, j+2]] for j in range(5)], dtype=float)
        top = np.exp2(mode.top_logs(mode.provider.all_pair_coefficients(np.log2(matrices), 1., slots=2), slots=2))
        t = .7
        left = sum(comb(2, i)*t**i*matrices[2+i] for i in range(3))
        np.testing.assert_allclose((top*np.array([t**k for k in range(5)])).sum(axis=1).reshape(2, 2),
                                   left@left, rtol=2e-14)


if __name__ == '__main__':
    unittest.main()
