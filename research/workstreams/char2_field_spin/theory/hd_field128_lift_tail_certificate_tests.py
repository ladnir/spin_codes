"""Tiny structural/directed tests; no q128 full-range numerical gate."""

from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

from flint import arb,ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_tail_certificate as tail
from hd_shifted_cover_backend import histogram_vertices,profile_count
from hd_randomized_outer_proof import Parameters,ball_log2


def multiply(a,b):
    return tuple(a[2*r]*b[c]+a[2*r+1]*b[2+c] for r in range(2) for c in range(2))


class DirectedTailTests(unittest.TestCase):
    def test_support_bounds_dominate_exact_indexed_cycles(self):
        with ctx.workprec(256):
            for pairs in (2,3,4):
                for weight in range(1,2*pairs+1):
                    matrices = {f:(F(3,4)+F(f,100),F(1,8),F(1,32),F(1,2))
                                for f in tail.full_counts(pairs,weight)}
                    exact = F(0)
                    for bits in product((0,1),repeat=2*pairs):
                        if sum(bits)!=weight:
                            continue
                        a = sum(bits[2*i]*bits[2*i+1] for i in range(pairs))
                        b = sum(bits[2*i+1]*bits[(2*i+2)%(2*pairs)] for i in range(pairs))
                        pair = multiply(matrices[a],matrices[b])
                        run = multiply(pair,pair)
                        exact += run[0]+run[1]
                    regions = {f:tuple(map(base.rational,row)) for f,row in matrices.items()}
                    for ratio in (F(1,8),F(1),F(8)):
                        for method in ('marginal','runs'):
                            bound = tail.support_sum(pairs,weight,regions,ratio,method=method,repetitions=4)
                            self.assertTrue(bound>=base.rational(exact))

    def test_histogram_box_geometry_replays_old(self):
        for ranges in (((1,5),(0,2),(0,3)),((4090,4100),(0,4096),(0,4096)),
                       ((8191,8192),(4095,4096),(4095,4096)),((2,2),(1,1),(1,1))):
            self.assertEqual(tail.histogram_vertices(ranges),histogram_vertices(ranges))
            self.assertEqual(tail.profile_count(ranges),profile_count(ranges))

    def test_regional_bank_min_after_normalization(self):
        bank = object.__new__(tail.PositiveBank)
        bank.markers = ((F(1,4),F(1,8)),(F(1,8),F(1,2)))
        bank.numerators = [tuple(map(arb,(1,8,2,6))),tuple(map(arb,(3,2,6,1)))]
        with ctx.workprec(256):
            result = bank.regional(3)
            for f in tail.full_counts(base.PAIRS,3):
                s = 3-2*f
                count = comb(base.PAIRS,f)*comb(base.PAIRS-f,s)
                expected = []
                for entry in range(4):
                    candidates = [row[entry]/base.rational(count*m[0]**s*m[1]**f)
                                  for row,m in zip(bank.numerators,bank.markers)]
                    expected.append(min(x.upper() for x in candidates))
                self.assertEqual(result[f],tuple(expected))

    def test_positive_pair_preserves_order(self):
        with ctx.workprec(256):
            # Noncommuting matrices exercise both GD and DG independently.
            good,dirty,composed = (1,2,3,4),(2,5,7,1),(11,13,17,19)
            single = {(0,0):tuple(map(arb,good)),(2,0):tuple(map(arb,dirty))}
            pair = {(0,0):tuple(map(arb,composed))}
            actual = tail.positive_pair(single,pair,(F(1),F(1)))
            expected = tuple(composed[i]+multiply(good,dirty)[i]+multiply(dirty,good)[i]
                             +multiply(dirty,dirty)[i] for i in range(4))
            self.assertEqual(actual,tuple(map(arb,expected)))

    def test_terminal_q256_replays_exact_ball(self):
        result = tail.terminal_receipt(q=256,first=7649,precision=256)
        with ctx.workprec(256):
            parameters = Parameters(field=256,state=16,limit=734003)
            ball = ball_log2(parameters,mark=F(result['witness']['mark']))
            expected = ball-7649*(96*8+32*tail.log2arb(arb(255)))-tail.log2arb(arb(255))
            self.assertGreaterEqual(F(result['log2_upper_rational']),base.rational_ceiling(expected))
        self.assertEqual(result['weight_range'],[7649,8192])
        self.assertTrue(result['directed'])
        self.assertTrue(result['no_support_subset_union'])
        self.assertFalse(result['whole_code_certificate'])

    def test_source_closure_and_scope_validation(self):
        pins = tail.source_pins()
        for name in ('hd_field128_lift_tail_certificate.py','hd_field128_lift_certificate.py',
                     'hd_common_refresh_principal.py','hd_grouped_quiet.py'):
            self.assertIn(name,pins)
        self.assertNotIn('hd_field128_lift_tail_screen.py',pins)
        for first in (0,8193):
            with self.assertRaises(ValueError):
                tail.terminal_receipt(q=256,first=first)
        with self.assertRaises(ValueError):
            tail.contraction(tuple(map(arb,(1,2,3,4))),F(0))


if __name__=='__main__':
    unittest.main()
