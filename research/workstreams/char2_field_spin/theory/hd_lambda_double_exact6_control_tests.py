import unittest
from math import log2

import numpy as np

import hd_lambda_double_exact6_control as control


class ExactSixControlTests(unittest.TestCase):
    def test_dependencies(self):
        control.dependencies()

    def test_raw_multiplicity_outside_root(self):
        table = np.full((65536, 2, 2), -36.)
        result = control.coefficients({(3, 9): 255, (3, 10): 255, (0, 0): 1}, table)
        np.testing.assert_allclose(result[3], log2(510)-2, rtol=0, atol=1e-14)
        np.testing.assert_allclose(result[0], -2, rtol=0, atol=1e-14)
        self.assertTrue(np.all(np.isneginf(result[1])))

    def test_exact_zero_and_degree_six(self):
        table = np.full((65536, 2, 2), -np.inf)
        table[0, 0, 0] = 0
        table[7] = -18
        result = control.coefficients({(0, 0): 1, (6, 7): 255**6}, table)
        self.assertEqual(result[0, 0, 0], 0)
        self.assertTrue(np.isneginf(result[0, 0, 1]))
        np.testing.assert_allclose(result[6], 6*log2(255)-1, rtol=0, atol=1e-14)

    def test_invalid_counts(self):
        table = np.zeros((65536, 2, 2))
        for counts in ({(7, 0): 1}, {(1, 65536): 1}, {(1, 0): 0}, {(1, 0): 1.0}):
            with self.assertRaises(ValueError): control.coefficients(counts, table)


if __name__ == '__main__':
    unittest.main()
