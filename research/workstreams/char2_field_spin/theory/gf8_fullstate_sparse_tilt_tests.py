"""Exact mark inequality and compatibility tests, not a certificate."""
from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np

from gf8_fullstate_screen import sparse_bound
from gf8_fullstate_sparse_tilt import fixed_coefficients, evaluate
from sparse_active_interval_screen import regional_from_coefficients


class SparseTiltTests(unittest.TestCase):
    def test_exact_nonnegative_support_restriction(self):
        # Arbitrary nonnegative coefficients, including deliberately invalid
        # low-weight terms that the marked relaxation must still pay for.
        for coefficients in itertools.product(range(3), repeat=5):
            for minimum in range(5):
                actual = sum(coefficients[minimum:])
                for mark in (F(1), F(3,2), F(7)):
                    bound = sum(c*mark**(w-minimum) for w,c in enumerate(coefficients))
                    self.assertLessEqual(actual, bound)

    def test_zero_mark_matches_old_bound(self):
        for state in (8,16):
            point = (6.2, math.log(6))
            coeff,g,tau,base = fixed_coefficients(3,point,state,k_log=12)
            self.assertAlmostEqual(evaluate(coeff,g,3,tau,base,0),
                                   sparse_bound(3,point,state,k_log=12),places=10)

    def test_region_mark_is_coefficient_mark(self):
        coeff,g,tau,_ = fixed_coefficients(3,(6.2,math.log(6)),8,k_log=12)
        mark = math.log(1.8)
        directly = coeff+np.arange(len(coeff))[:,None,None]*mark
        np.testing.assert_allclose(regional_from_coefficients(directly,3,g['rows'],tau),
                                   regional_from_coefficients(coeff,3,g['rows'],tau+mark),rtol=1e-13)


if __name__ == '__main__':
    unittest.main()
