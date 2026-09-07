"""Cross-check every polynomial coefficient, including extreme scaled tails."""
import unittest

import numpy as np

import occupation_positive_poly_v1 as positive
import occupation_refresh_v1 as log_reference


class PositiveTests(unittest.TestCase):
    def test_random_extreme_coefficients(self):
        rng = np.random.default_rng(192)
        for spread in (1., 100., 1000.):
            left = rng.uniform(-spread, spread, size=(13, 4, 4))
            right = rng.uniform(-spread, spread, size=(9, 4, 4))
            left[:, 1, 2] = -np.inf
            right[0, 0, 3] = -np.inf
            expected = log_reference.polynomial_product(left, right, 19)
            actual = positive.polynomial_product(left, right, 19)
            np.testing.assert_allclose(actual, expected, rtol=0, atol=3e-12)

    def test_all_finite_field_region_coefficients(self):
        prepared = log_reference.Epochs(4, 2, {2: 2, 4: 1}, [1, 0, 2, 0, 1])
        for lam in (.001, .5, 5., 100.):
            epoch = prepared.at(lam)
            for length, maximum in ((16, 16), (256, 64), (1024, 128)):
                expected = log_reference.region_logs(epoch, 4, length, maximum)
                actual = positive.region_logs(epoch, 4, length, maximum)
                np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-8)


if __name__ == '__main__':
    unittest.main()
