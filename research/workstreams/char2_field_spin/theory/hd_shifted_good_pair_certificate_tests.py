import unittest
from fractions import Fraction as F

import numpy as np
from flint import arb,ctx

from hd_shifted_alias_certificate import _cadd,interval_fft,selected_inverse
from hd_shifted_general_certificate import general_coefficients
from hd_shifted_good_pair_certificate import (components,good_coefficients,
    interval_polynomial,matrix_product,inverse_alias_box)


class DirectedGoodPairTests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.workprec(256)
        self.precision.__enter__()

    def tearDown(self):
        self.precision.__exit__(None,None,None)

    def test_degree_64_is_folded_with_outward_sum(self):
        coeff = {(0,0,0):tuple(map(arb,(1,2,3,4))),
                 (0,0,64):tuple(map(arb,(5,6,7,8)))}
        values = interval_polynomial(coeff,(arb(1),)*4,64)
        for entry,total in enumerate((6,8,10,12)):
            self.assertLessEqual(values[0][0][entry,0,0,0],total)
            self.assertGreaterEqual(values[0][1][entry,0,0,0],total)

    def test_interval_product_retains_order(self):
        left = interval_polynomial({(0,0,0):tuple(map(arb,(1,2,0,1)))},(arb(1),)*4,2)
        right = interval_polynomial({(0,0,0):tuple(map(arb,(1,0,3,1)))},(arb(1),)*4,2)
        actual = matrix_product(left,right)
        for entry,target in enumerate((7,2,3,1)):
            self.assertLessEqual(actual[0][0][entry,0,0,0],target)
            self.assertGreaterEqual(actual[0][1][entry,0,0,0],target)

    def test_three_components_match_constant_polynomial(self):
        g = {(0,0,0):tuple(map(arb,(1,2,0,1)))}
        d = {(0,0,0):tuple(map(arb,(1,0,3,1)))}
        # C=G^2. Then C+GD+DG+DD=(G+D)^2, whose coefficient is exact.
        c = {(0,0,0):tuple(map(arb,(1,4,0,1)))}
        transformed = [interval_fft(interval_polynomial(v,(arb(1),)*4,2),(1,2,3))
                       for v in (g,d,c)]
        g,d,c = transformed
        result = _cadd(c,_cadd(matrix_product(g,d),_cadd(matrix_product(d,g),matrix_product(d,d))))
        selected = selected_inverse(result,(0,0,0))
        for entry,target in enumerate((10,8,12,10)):
            self.assertLessEqual(selected[0][0][entry],target)
            self.assertGreaterEqual(selected[0][1][entry],target)

    def test_markers_and_target_truncation(self):
        single = [[arb(0) for _ in range(33**3)] for _ in range(4)]
        for entry in range(4):
            single[entry][33**2] = arb(entry+1)
            single[entry][33] = arb(entry+2)
        paired = {(2,64):list(map(arb,(1,2,3,4)))}
        g,d,c = components(single,paired,(F(2),F(3),F(1)),(1,2,64))
        self.assertEqual(g[1,0,0][0],arb(2))
        self.assertEqual(g[0,1,0][0],arb(6))
        self.assertEqual(c[1,1,64][0],arb(12))
        self.assertEqual(c[0,2,64][0],arb(9))
        self.assertNotIn((2,0,64),c)

    def test_inverse_neighborhood_caps_exact_coefficients_and_wraps(self):
        values = (np.arange(4*4**3,dtype=float).reshape(4,4,4,4)%17)/16
        zeros = np.zeros_like(values)
        for ranges in (((1,2),(2,3),(0,2)),((3,5),(0,8),(4,4))):
            transformed = interval_fft(((values.copy(),values.copy()),(zeros.copy(),zeros.copy())),(1,2,3))
            caps = inverse_alias_box(transformed,ranges)
            indices = [np.unique(np.arange(lo,hi+1)%4) for lo,hi in ranges]
            exact = np.max(values[np.ix_(np.arange(4),*indices)],axis=(1,2,3))
            self.assertTrue(np.all(caps[0][0]<=exact))
            self.assertTrue(np.all(exact<=caps[0][1]))

    def test_alpha_one_agrees_with_independent_single_step_polynomial(self):
        single = general_coefficients(F(1,2),F(1),(F(1),)*3)
        actual = good_coefficients(F(1,2),F(1))
        reference = {}
        for i in range(2):
            for j in range(33-i):
                left = tuple(row[i*33**2+j] for row in single)
                for k in range(2):
                    for ell in range(33-k):
                        right = tuple(row[k*33**2+ell] for row in single)
                        key = i+k,j+ell
                        if key not in reference:
                            reference[key] = [arb(0) for _ in range(4)]
                        for r in range(2):
                            for s in range(2):
                                reference[key][2*r+s] += left[2*r]*right[s]+left[2*r+1]*right[2+s]
        for key,values in reference.items():
            for value,expected in zip(actual[key],values):
                self.assertTrue((value-expected).contains(0),(key,value,expected))


if __name__=='__main__':
    unittest.main()
