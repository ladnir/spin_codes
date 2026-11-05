from fractions import Fraction as F
from itertools import product
import unittest

import numpy as np

from hd_bit_syndrome_activity_control import log_difference, marked_columns


class ActivityTests(unittest.TestCase):
    def test_scalar_majorant(self):
        for r, y, k in product(range(1, 10), (F(1), F(3, 2), F(5, 2)),
                               (F(1), F(101, 100), F(3), F(256))):
            self.assertLessEqual((1+y*(k-1))**r, 1+y**r*(k**r-1))

    def test_log_difference(self):
        result = log_difference(np.log2([1, 3, 7, 2**50]),
                                np.array([-np.inf, 0, np.log2(7), 50-1e-10]))
        self.assertEqual(result[0], 0)
        self.assertAlmostEqual(result[1], 1)
        self.assertEqual(result[2], -np.inf)
        self.assertTrue(np.isfinite(result[3]))
        with self.assertRaises(ArithmeticError):
            log_difference([0], [1e-12])

    def test_marker_matrix(self):
        base = np.array([[1., 2.], [3., 4.]])
        column = np.stack((base+1, base+3))
        for y in (F(1), F(2), F(5, 2)):
            actual = np.exp2(marked_columns(np.log2(column), np.log2(base), y, 3))
            expected = base+float(y)**3*(column-base)
            np.testing.assert_allclose(actual, expected, rtol=1e-14)

    def test_minimum_column_floor(self):
        # Exact toyword: onlywords with>=2nonzerocolumns survive.
        for y in (F(1), F(3, 2), F(5, 2)):
            kernels = (F(1, 5), F(2, 3), F(4, 7))
            exact = sum(np.prod([kernels[i] for i in range(3) if mask>>i&1])
                        for mask in range(8) if mask.bit_count() >= 2)
            marked = np.prod([1+y*k for k in kernels])/y**2
            self.assertLessEqual(exact, marked)


if __name__ == '__main__':
    unittest.main()
