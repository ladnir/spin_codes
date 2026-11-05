"""Tiny native-ACB tests; no construction-size certificate gate."""

from fractions import Fraction as F
from math import comb
import unittest

from flint import acb,arb,ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_acb_alias_certificate as native
from hd_field128_lift_alias_certificate_tests import multiply,power


class NativeAliasTests(unittest.TestCase):
    def test_native_inverse_normalization(self):
        entries = tuple([acb(F((i+3*j)%11,7).numerator)/7 for i in range(16)] for j in range(4))
        with ctx.workprec(192):
            restored = native.dft2(native.dft2(entries,4),4,inverse=True)
            for got,want in zip(restored,entries):
                for actual,expected in zip(got,want):
                    self.assertTrue(actual.contains(expected))

    def test_noncommuting_pairing_and_degree64_alias(self):
        good,dirty,composed = (1,2,3,4),(2,5,7,1),(11,13,17,19)
        expected = power(tuple(composed[i]+multiply(good,dirty)[i]+multiply(dirty,good)[i]
                               +multiply(dirty,dirty)[i] for i in range(4)),64)
        with ctx.workprec(192):
            atlas,fallback,diagnostic = native.component_alias(
                {(0,0):tuple(map(arb,good))},{(0,0):tuple(map(arb,dirty))},
                {(0,64):tuple(map(arb,composed))},64,base.Budget(20))
            actual = atlas.numerator(((0,0),(0,0)))
            for got,positive,want in zip(actual,fallback,expected):
                self.assertTrue(got.upper()>=arb(want).lower())
                self.assertTrue(positive>=arb(want))
                self.assertTrue(got/arb(want)<arb(1)+arb(2)**-100)
            self.assertEqual(diagnostic['transform_precision'],192)

    def test_positive_residues_and_tiny_coefficients(self):
        matrix = (F(1,4),F(1,8),F(1,16),F(1,4))
        final = power(matrix,128)
        with ctx.workprec(256):
            dirty = {(0,0):tuple(map(base.rational,matrix)),
                     (1,0):tuple(base.rational(x/4) for x in matrix)}
            atlas,_,_ = native.component_alias({},dirty,{},64,base.Budget(20))
            for residue in (0,1,2,31,63):
                scalar = sum((F(comb(128,n),4**n) for n in range(residue,129,64)),F(0))
                for got,entry in zip(atlas.numerator(((residue,residue),(0,0))),final):
                    self.assertTrue(got.upper()>=base.rational(scalar*entry).lower())
            # Arithmetic range is independent of binary64's smallest subnormal.
            tiny = tuple(base.rational(x)*arb(2)**-2000 for x in matrix)
            atlas,_,_ = native.component_alias({}, {},{(0,0):tiny},64,base.Budget(20))
            for got,want in zip(atlas.numerator(((0,0),(0,0))),power(matrix,64)):
                self.assertTrue(got>0)
                self.assertTrue(got.upper()>=(base.rational(want)*arb(2)**-128000).lower())

    def test_source_closure_has_both_backends(self):
        pins = native.source_pins()
        for name in ('hd_field128_lift_acb_alias_certificate.py','hd_field128_lift_alias_certificate.py',
                     'hd_field128_lift_certificate.py','hd_field128_lift_tail_certificate.py'):
            self.assertIn(name,pins)
        self.assertNotIn('hd_field128_lift_alias_screen.py',pins)


if __name__=='__main__':
    unittest.main()
