"""Independent implementation checks for native Q1 coefficients."""
import unittest

import numpy as np

import activation_q1 as reference
from activation_q1_native import Q1Kernel


class NativeQ1Test(unittest.TestCase):
    def test_dense_matrices_and_endpoints(self):
        random = np.random.default_rng(17234)
        zero, one = [np.log(random.uniform(.001, .8, size=(3,3,3))) for _ in range(2)]
        kernel = Q1Kernel()
        for length in (1,2,7,32):
            actual = kernel.coefficients(zero, one, length)
            expected = reference.coefficient_logs(zero, one, length)
            np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-12)

    def test_activation_matrices_and_extreme_log_values(self):
        spectrum = {32:126,64:1}
        kernel = Q1Kernel()
        for epochs in (1,8,65536):
            regions = reference.region_logs(*reference.epoch_logs(64,7,spectrum,np.exp([-16.,-7.,-1.,0.])),epochs)
            actual = kernel.coefficients(*regions, 128)
            expected = reference.coefficient_logs(*regions, 128)
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-9)

    def test_invalid_geometry_is_rejected(self):
        kernel = Q1Kernel()
        with self.assertRaises(ValueError):
            kernel.coefficients(np.zeros((1,3,3)), np.zeros((1,3,3)), 0)


if __name__ == '__main__':
    unittest.main()
