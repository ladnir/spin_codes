"""Check the native four-state kernel against independent implementations."""
import math
import unittest
from fractions import Fraction

import numpy as np

import activation_q1_refresh as reference
from activation_refresh_native import RefreshKernel
from test_bch_growth import exact_refresh
from test_activation_q1 import positive_regions, exact_coefficients


class NativeRefreshTest(unittest.TestCase):
    def test_matches_positive_rationals(self):
        kernel=RefreshKernel()
        for z in (Fraction(1,3),Fraction(3,4)):
            matrices=reference.epoch_logs(4,2,{2:2,4:1},np.array([-math.log(float(z))]))
            for epochs in (1,3,9):
                expected=exact_coefficients(*positive_regions(*exact_refresh(z),epochs),4)
                actual=kernel.coefficients(*matrices,epochs,4)[0]
                np.testing.assert_allclose(actual,[math.log(float(v)) for v in expected],atol=3e-13,rtol=2e-13)

    def test_full_length_coefficients_match_python(self):
        kernel=RefreshKernel()
        matrices=reference.epoch_logs(4,2,{2:2,4:1},np.array([1e-6,.003,.7]))
        for block,epochs in ((8,1),(128,64),(512,256)):
            np.testing.assert_allclose(kernel.coefficients(*matrices,epochs,block),
                reference.coefficient_logs(*matrices,epochs,block),atol=2e-9,rtol=2e-12)

    def test_invalid_geometry_rejected(self):
        kernel=RefreshKernel()
        with self.assertRaises(ValueError):
            kernel.coefficients(np.zeros((1,3,3)),np.zeros((1,3,3)),1,8)


if __name__=='__main__': unittest.main()
