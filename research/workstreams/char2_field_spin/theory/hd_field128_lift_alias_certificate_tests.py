"""Tiny exact-reference tests; no full-field numerical certificate gate."""

from fractions import Fraction as F
from math import comb
import unittest

from flint import arb,ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_alias_certificate as alias


def multiply(a,b):
    return tuple(a[2*r]*b[c]+a[2*r+1]*b[2+c] for r in range(2) for c in range(2))


def power(a,n):
    result = (1,0,0,1)
    while n:
        if n&1:
            result = multiply(result,a)
        n >>= 1
        if n:
            a = multiply(a,a)
    return result


class DirectedAliasTests(unittest.TestCase):
    def test_similarity_exceeds_binary64_range(self):
        with ctx.workprec(256):
            for exponents in ((1400,1300,1250,1399),(1399,1250,1300,1400),
                              (-1500,-1700,-1600,-1501)):
                entries = tuple(arb(2)**n for n in exponents)
                vector,rho = alias.similarity(entries)
                self.assertTrue(rho>0 and rho.is_finite())
                self.assertTrue(vector[1]>0 and vector[1].is_finite())
                for r in range(2):
                    row = (entries[2*r]*vector[0]+entries[2*r+1]*vector[1])/vector[r]
                    self.assertTrue(rho*rho>=row)

    def test_component_split_and_exact_marker_tilts(self):
        with ctx.workprec(256):
            single = {(0,0):tuple(map(arb,(1,2,3,4))),
                      (1,1):tuple(map(arb,(2,3,4,5))),
                      (2,0):tuple(map(arb,(3,4,5,6)))}
            composed = {(2,64):tuple(map(arb,(5,6,7,8)))}
            markers = (F(1,4),F(1,2))
            good,dirty,paired = alias.components(single,composed,markers)
            self.assertEqual(set(good),{(0,0),(1,1)})
            self.assertEqual(set(dirty),{(2,0)})
            self.assertEqual(good[1,1],tuple(arb(x)/8 for x in (2,3,4,5)))
            self.assertEqual(paired[2,64],tuple(arb(x)*arb(2)**-68 for x in (5,6,7,8)))
            good,dirty,paired = alias.components(single,None,markers)
            self.assertFalse(good or paired)
            self.assertEqual(set(dirty),set(single))

    def test_constant_noncommuting_pair_enclosure(self):
        good,dirty,composed = (1,2,3,4),(2,5,7,1),(11,13,17,19)
        pair = tuple(composed[i]+multiply(good,dirty)[i]+multiply(dirty,good)[i]
                     +multiply(dirty,dirty)[i] for i in range(4))
        expected = power(pair,64)
        with ctx.workprec(256):
            atlas,fallback,_ = alias.component_alias(
                {(0,0):tuple(map(arb,good))},{(0,0):tuple(map(arb,dirty))},
                {(0,0):tuple(map(arb,composed))},64,base.Budget(10))
            actual = atlas.numerator(((0,0),(0,0)))
            for got,positive,want in zip(actual,fallback,expected):
                self.assertTrue(got.upper()>=arb(want).lower())
                self.assertTrue(positive>=arb(want))
                self.assertLess(float(got/arb(want)),1.000001)

    def test_complete_positive_residue_sums_not_selected_coefficients(self):
        matrix = (F(1,4),F(1,8),F(1,16),F(1,4))
        final = power(matrix,128)
        with ctx.workprec(256):
            dirty = {(0,0):tuple(map(base.rational,matrix)),
                     (1,0):tuple(base.rational(x/4) for x in matrix)}
            atlas,_,_ = alias.component_alias({},dirty,{},64,base.Budget(10))
            for residue in (0,1,2,31,63):
                scalar = sum((F(comb(128,n),4**n) for n in range(residue,129,64)),F(0))
                actual = atlas.numerator(((residue,residue),(0,0)))
                for got,entry in zip(actual,final):
                    self.assertTrue(got.upper()>=base.rational(scalar*entry).lower())
            # Degree64 is folded, not discarded at the transform boundary.
            composed = {(0,64):tuple(map(base.rational,matrix))}
            atlas,_,_ = alias.component_alias({}, {},composed,64,base.Budget(10))
            for got,want in zip(atlas.numerator(((0,0),(0,0))),power(matrix,64)):
                self.assertTrue(got.upper()>=base.rational(want).lower())

    def test_normalized_bank_min_and_independent_positive_fallback(self):
        class Atlas:
            def __init__(self,values):
                self.values = tuple(map(arb,values))
            def numerator(self,ranges):
                return self.values
        bank = object.__new__(alias.AliasBank)
        bank.markers = ((F(1,4),F(1,8)),(F(1,8),F(1,2)))
        bank.atlases = [Atlas((1,8,2,6)),Atlas((3,2,6,1))]
        bank.numerators = [tuple(map(arb,(4,4,4,4))),tuple(map(arb,(4,4,4,4)))]
        with ctx.workprec(256):
            result = bank.regional(3)
            for full,row in result.items():
                single = 3-2*full
                count = comb(base.PAIRS,full)*comb(base.PAIRS-full,single)
                candidates = [tuple(min(a,b)/base.rational(count*m[0]**single*m[1]**full)
                                    for a,b in zip(atlas.values,fallback))
                              for m,atlas,fallback in zip(bank.markers,bank.atlases,bank.numerators)]
                self.assertEqual(row,tuple(min(x[i].upper() for x in candidates) for i in range(4)))

    def test_pin_closure_and_invalid_grid(self):
        pins = alias.source_pins()
        for name in ('hd_field128_lift_alias_certificate.py','hd_field128_lift_certificate.py',
                     'hd_field128_lift_tail_certificate.py','hd_shifted_alias_certificate.py',
                     'hd_shifted_good_pair_certificate.py','hd_shifted_coverage.py'):
            self.assertIn(name,pins)
        self.assertNotIn('hd_field128_lift_alias_screen.py',pins)
        with self.assertRaises(ValueError):
            alias.component_alias({}, {},{},32,base.Budget(10))


if __name__=='__main__':
    unittest.main()
