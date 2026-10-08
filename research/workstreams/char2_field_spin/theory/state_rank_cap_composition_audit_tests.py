"""Independent exact checks of rank-cap cone composition; no certificate."""

from fractions import Fraction as F
from math import comb, exp
import unittest

import numpy as np

from mask_regional_screen import normalized_power
from rank_cap_regional_screen import region
from state_rank_cap_transfer import local_matrix


def cone_generators():
    return [(F(1),)+(F(0),)*8] + [
        tuple(F(int(i >= j)) for i in range(9)) for j in range(1, 9)]


def assert_cone(test, matrix):
    """The nine rays generate the entire required continuation cone."""
    for generator in cone_generators():
        values = matrix @ np.asarray(generator, dtype=object)
        test.assertTrue(all(value >= 0 for value in values))
        test.assertTrue(all(values[i] <= values[i+1] for i in range(1, 8)))


def sequential_coefficients(packets, degree, y, tau):
    """Raw positive matrix polynomial, with independent sequential products."""
    local = [np.asarray(local_matrix(w, y, tau), dtype=object)*comb(16, w)
             for w in range(min(degree, 16)+1)]
    coefficients = [np.eye(9, dtype=object)]
    for _ in range(packets):
        following = [np.zeros((9, 9), dtype=object)
                     for _ in range(min(degree+1, len(coefficients)+len(local)-1))]
        for left_degree, left in enumerate(coefficients):
            for right_degree, right in enumerate(local[:len(following)-left_degree]):
                following[left_degree+right_degree] += left @ right
        coefficients = following
    return coefficients


class RankCapCompositionAudit(unittest.TestCase):
    def test_normalized_coefficients_equal_exact_sequential_polynomial(self):
        y, tau, degree, packets = F(9, 10), F(7, 2), 2, 3
        local = np.asarray([local_matrix(w, y, tau) for w in range(degree+1)], dtype=object)
        actual = normalized_power(local, 16, packets, degree, exact=True)
        raw = sequential_coefficients(packets, degree, y, tau)
        for weight in range(degree+1):
            expected = raw[weight]*F(1, comb(16*packets, weight))
            self.assertTrue(np.array_equal(actual[weight], expected))
            assert_cone(self, actual[weight])

    def test_regional_scaling_matches_exact_positive_coefficients(self):
        y, tau, active, rows = F(9, 10), F(7, 2), 2, 48
        raw = sequential_coefficients(rows//16, active, y, tau)
        expected = sum((matrix*F(comb(active, weight), comb(rows, weight))*tau**weight
                        for weight, matrix in enumerate(raw)), np.zeros((9, 9), dtype=object))
        scaled, log_scale = region(active, float(y), float(tau), rows)
        np.testing.assert_allclose(scaled*exp(log_scale), np.asarray(expected, dtype=float),
                                   rtol=5e-13, atol=1e-300)
        assert_cone(self, expected)
        # A positive regional mixture and repeated regions preserve the cone.
        assert_cone(self, expected @ expected)

    def test_uniform_abstract_rank_caps_need_not_be_exact_distributions(self):
        # The proof consumes tail domination, not entrywise rank domination.
        # These unequal-output-order products test the full continuation cone.
        y, tau = F(3, 4), F(19, 2)
        matrices = [np.asarray(local_matrix(w, y, tau), dtype=object) for w in (0, 1, 2, 9)]
        for left, right in zip(matrices, reversed(matrices)):
            assert_cone(self, left @ right)


if __name__ == '__main__':
    unittest.main()
