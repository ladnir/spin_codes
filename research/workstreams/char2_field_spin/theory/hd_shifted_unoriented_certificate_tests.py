import unittest
from fractions import Fraction as F
from itertools import product
from unittest.mock import patch

import numpy as np
from flint import arb,ctx

from hd_shifted_alias_certificate import interval_fft
from hd_shifted_general_certificate import general_coefficients
from hd_shifted_good_pair_certificate import interval_polynomial
from hd_shifted_coverage import RegionalAliases
from hd_shifted_unoriented_certificate import collapse_coefficients,unsigned_components,global_box_receipt


class UnorientedCertificateTests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.workprec(256)
        self.precision.__enter__()

    def tearDown(self):
        self.precision.__exit__(None,None,None)

    def test_orientation_count_is_divided_out_before_maximum(self):
        single = [[arb(0) for _ in range(33**3)] for _ in range(4)]
        for entry in range(4):
            single[entry][2*33] = arb(2)
            single[entry][33**2+33] = arb(10)
            single[entry][2*33**2] = arb(3)
        actual = collapse_coefficients(single)
        self.assertEqual(actual[2,0],(arb(5),)*4)
        self.assertNotIn((3,0),actual)

    def test_zero_and_one_singleton_need_no_extra_orientation_loss(self):
        single = general_coefficients(F(1,2),F(1),(F(1),)*3)
        actual = collapse_coefficients(single)
        for s in range(2):
            for f in range(33-s):
                for entry in range(4):
                    expected = single[entry][s*33**2+f]
                    self.assertTrue(actual[s,f][entry]>=expected.lower())
                    self.assertTrue(abs(actual[s,f][entry]-expected)<arb('1e-40'))

    def test_unsigned_markers_do_not_introduce_two_to_the_singleton_count(self):
        single = {(1,2):(arb(3),)*4,(2,1):(arb(5),)*4}
        paired = {(2,64):(arb(7),)*4}
        good,dirty,composed = unsigned_components(single,paired,(F(2),F(3)))
        self.assertEqual(good[1,2],(arb(54),)*4)
        self.assertEqual(dirty[2,1],(arb(60),)*4)
        self.assertEqual(composed[2,64],(arb(28)*arb(3)**64,)*4)

    def test_two_dimensional_polynomial_folds_degree_64(self):
        coeff = {(0,0):(arb(2),)*4,(2,64):(arb(3),)*4}
        array = interval_polynomial(coeff,(arb(1),)*4,64,dimensions=2)
        self.assertEqual(array[0][0].shape,(4,64,64))
        for entry in range(4):
            self.assertLessEqual(array[0][0][entry,2,0],3)
            self.assertGreaterEqual(array[0][1][entry,2,0],3)

    def test_two_dimensional_alias_normalization_and_wrapping(self):
        raw = np.arange(4*4**2,dtype=float).reshape(4,4,4)/64
        zeros = np.zeros_like(raw)
        powered = interval_fft(((raw.copy(),raw.copy()),(zeros.copy(),zeros.copy())),(1,2))
        atlas = RegionalAliases.from_power(powered,(arb(1),arb(1)),arb(1))
        actual = atlas.numerator(((3,5),(1,1)))
        for entry in range(4):
            expected = max(raw[entry,i,1] for i in (3,0,1))
            self.assertTrue(actual[entry]>=arb(float(expected)))

    def test_directed_global_box_dominates_every_small_interior_profile(self):
        single = {(0,0):tuple(map(arb,(1,2,3,4)))}
        paired = {(0,0):tuple(map(arb,(5,6,7,8)))}
        ranges = ((1023,1025),(127,129),(127,129))
        markers = ((F(2,3),F(3,5)),(F(4,5),F(2,5)))
        with patch('hd_shifted_unoriented_certificate.general_coefficients',return_value=None), \
             patch('hd_shifted_unoriented_certificate.collapse_coefficients',return_value=single), \
             patch('hd_shifted_unoriented_certificate.good_coefficients',return_value=paired):
            result = global_box_receipt(ranges,F(1,2),F(1,8),markers)
            upper = F(result['log2_upper_rational'])
            self.assertEqual(result['integer_profile_triples'],27)
            for point in product(*(range(lo,hi+1) for lo,hi in ranges)):
                actual = global_box_receipt(tuple((n,n) for n in point),F(1,2),F(1,8),markers)
                self.assertLessEqual(F(actual['log2_upper_rational']),upper)
        with self.assertRaises(ValueError):
            global_box_receipt(((1,2),(2,3),(0,0)),F(1,2),F(1,8),markers)


if __name__=='__main__':
    unittest.main()
