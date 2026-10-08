"""Independent small polynomial checks for the weight-cut diagnostic."""

import math
import unittest

import numpy as np
from scipy.special import logsumexp

from cauchy_regional_rank_screen import log_region
from cauchy_weight_cut_screen import log_power, region_coefficients


class WeightCutTests(unittest.TestCase):
    def test_noncommuting_positive_polynomial(self):
        local = np.array([[[2., 1.], [3., 2.]], [[4., 7.], [1., 5.]]])
        exact = np.eye(2)[None, :, :]
        for _ in range(5):
            result = np.zeros((len(exact)+1, 2, 2))
            for i, left in enumerate(exact):
                for j, right in enumerate(local):
                    result[i+j] += left @ right
            exact = result
        np.testing.assert_allclose(np.exp(log_power(np.log(local), 5, 5)), exact)

    def test_truncation_does_not_change_low_coefficients(self):
        local = np.log(np.array([[[1., 2.], [3., 4.]], [[5., 1.], [2., 1.]]]))
        np.testing.assert_allclose(log_power(local, 7, 3), log_power(local, 7, 7)[:4])

    def test_region_sum_matches_existing_operator(self):
        for active in (1, 2, 3):
            for rows in (16, 32):
                actual = logsumexp(region_coefficients(active, .91, 17., rows), axis=0)
                expected = log_region(active, .91, 17., rows, joint_returns=True)
                np.testing.assert_allclose(actual, expected, atol=1e-12)

    def test_scalar_binomial_weight_cut(self):
        # Diagonal entries reproduce (1+2x)^6 exactly, including the tail.
        local = np.full((2, 2, 2), -math.inf)
        local[0, 0, 0] = local[0, 1, 1] = 0.
        local[1, 0, 0] = local[1, 1, 1] = math.log(2.)
        result = log_power(local, 6, 6)
        expected = sum(math.comb(6, w)*2**w for w in range(4, 7))
        self.assertAlmostEqual(math.exp(logsumexp(result[4:, 0, 0])), expected)


if __name__ == '__main__':
    unittest.main()
