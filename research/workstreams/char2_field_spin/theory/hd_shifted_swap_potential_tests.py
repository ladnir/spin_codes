"""Exact toy checks for the common-potential matching-marginal bound."""
from fractions import Fraction as F
from math import comb
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx

from hd_shifted_pairs_count import exact_profiles, one_pairing_count
from hd_shifted_swap_potential import (contraction_arb, full_counts, marginal_sum_arb,
    marginal_log_screen, optimize_potential_screen, positive_vector,
    weighted_transfer_arb,weighted_count_upper_arb,weighted_marginal_sum_arb,
    directed_bank_regions,checkpoint_receipt,SAVED_RUN_WITNESSES)
from hd_shifted_alias_certificate import rational


def multiply(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


class PotentialTests(unittest.TestCase):
    def test_both_matching_marginals(self):
        for pairs in range(1, 7):
            profiles = exact_profiles(pairs)
            for weight in range(2*pairs+1):
                for full in full_counts(pairs, weight):
                    expected = one_pairing_count(pairs, weight, full)
                    self.assertEqual(sum(n for (r,a,b), n in profiles.items()
                                         if r == weight and a == full), expected)
                    self.assertEqual(sum(n for (r,a,b), n in profiles.items()
                                         if r == weight and b == full), expected)

    def test_cauchy_without_independent_matchings(self):
        for pairs in range(2, 7):
            profiles = exact_profiles(pairs)
            for weight in range(1, 2*pairs):
                scalars = {f: F((f+1)**2+2, f+3) for f in full_counts(pairs, weight)}
                exact = sum(n*scalars[a]**3*scalars[b]**3
                            for (r,a,b), n in profiles.items() if r == weight)
                upper = sum(one_pairing_count(pairs, weight, f)*value**6
                            for f, value in scalars.items())
                self.assertLessEqual(exact, upper)

    def test_ordered_noncommuting_products(self):
        with ctx.workprec(192):
            pairs, repetitions, vector = 5, 6, (F(3), F(2))
            profiles = exact_profiles(pairs)
            for weight in range(1, 2*pairs):
                regions = {f: tuple(map(F, (1+f, 2, 3+f, 2+2*f)))
                           for f in full_counts(pairs, weight)}
                exact = F(0)
                for (r,a,b), count in profiles.items():
                    if r != weight:
                        continue
                    step = multiply(regions[a], regions[b])
                    run = (F(1), F(0), F(0), F(1))
                    for _ in range(repetitions//2):
                        run = multiply(run, step)
                    exact += count*(run[0]+run[1])
                upper = marginal_sum_arb(pairs, weight, regions, vector, repetitions=repetitions)
                self.assertGreaterEqual(upper, arb(exact.numerator)/exact.denominator)

    def test_one_common_vector_required(self):
        # Individually perfect Perron contractions cannot be multiplied when
        # their vectors differ: each nilpotent matrix has spectral radius0,
        # but their ordered product has a positive eigenvalue.
        first, second = (0, 2, 0, 0), (0, 0, 2, 0)
        self.assertEqual(multiply(first, second), (4, 0, 0, 0))
        self.assertEqual(contraction_arb(first, (1,1)), arb(2))
        self.assertEqual(contraction_arb(second, (1,1)), arb(2))

    def test_weighted_cycle_transfer_exactly(self):
        with ctx.workprec(256):
            for pairs in range(1,7):
                profiles = exact_profiles(pairs)
                for t,s,f in ((F(2),F(3),F(5)),(F(2,3),F(4,5),F(7,6))):
                    exact = sum(n*s**(r-2*a)*f**a*t**(a+b)
                                for (r,a,b),n in profiles.items())
                    actual = weighted_transfer_arb(pairs,t,s,f)
                    self.assertTrue(actual.contains(arb(exact.numerator)/exact.denominator))
                    for weight in range(2*pairs+1):
                        for full in full_counts(pairs,weight):
                            exact_count = sum(n*t**(a+b) for (r,a,b),n in profiles.items()
                                              if r==weight and a==full)
                            upper = weighted_count_upper_arb(pairs,weight,full,t,(s,f))
                            self.assertGreaterEqual(upper,arb(exact_count.numerator)/exact_count.denominator)

    def test_weighted_cauchy_with_noncommuting_matrices(self):
        with ctx.workprec(256):
            pairs,weight,steps,vector = 5,6,6,(F(2),F(3))
            regions = {f: tuple(map(F,(1+f,2,3+f,2+2*f))) for f in full_counts(pairs,weight)}
            exact = F(0)
            for (r,a,b),n in exact_profiles(pairs).items():
                if r!=weight:
                    continue
                step = multiply(regions[a],regions[b])
                run = (F(1),F(0),F(0),F(1))
                for _ in range(steps//2):
                    run = multiply(run,step)
                exact += n*(run[0]+run[1])
            markers = {f:(F(1+f,3),F(2+f,5)) for f in regions}
            for tilt in (F(1,3),F(1),F(7,2)):
                upper = weighted_marginal_sum_arb(pairs,weight,regions,vector,tilt,markers,repetitions=steps)
                self.assertGreaterEqual(upper,arb(exact.numerator)/exact.denominator)

    def test_float_matches_exact_scalar_fixture(self):
        pairs, weight = 6, 5
        regions = {f: (F(f+1), F(0), F(0), F(f+1)) for f in full_counts(pairs,weight)}
        logs = {f: np.array([[np.log2(f+1), -np.inf], [-np.inf, np.log2(f+1)]])
                for f in regions}
        with ctx.workprec(192):
            exact = marginal_sum_arb(pairs, weight, regions, (1,1), repetitions=4)
            expected = float(exact.log()/arb(2).log())
        self.assertAlmostEqual(marginal_log_screen(pairs, weight, logs, 0, repetitions=4), expected)
        result = optimize_potential_screen(pairs, weight, logs, repetitions=4)
        self.assertAlmostEqual(result['log2_marginal_screen'], expected)

    def test_invalid_inputs(self):
        for vector in ((0,1), (1,-1), (1.,1)):
            with self.assertRaises(ValueError):
                positive_vector(vector)
        with self.assertRaises(ValueError):
            full_counts(3,7)
        with self.assertRaises(ValueError):
            contraction_arb((-1,0,0,1), (1,1))

    def test_bank_minimum_follows_histogram_normalization(self):
        rows = ((F(2),F(3)),(F(5),F(7)))
        raw = {rows[0]:(4,1,5,2),rows[1]:(2,6,1,3)}
        class Atlas:
            def __init__(self,row):
                self.row = row
            def numerator(self,target):
                singleton,full = target[0][0],target[1][0]
                normalizer = comb(4096,singleton)*comb(4096-singleton,full)
                tilt = rational(self.row[0])**singleton*rational(self.row[1])**full
                return tuple(normalizer*tilt*x for x in raw[self.row])
        with patch('hd_shifted_swap_certificate.swap_unmarked',return_value=(None,None)) as coefficients, \
             patch('hd_shifted_unoriented_certificate.unsigned_components',side_effect=lambda a,b,row:(row,)), \
             patch('hd_shifted_good_pair_certificate.component_power',side_effect=lambda row,*a,**k:(row,None,None)) as powering, \
             patch('hd_shifted_coverage.RegionalAliases.from_power',side_effect=lambda row,*a:Atlas(row)):
            regions = directed_bank_regions(2,F(1,2),F(1,3),(*rows,rows[0]))
        self.assertEqual(coefficients.call_count,1)
        self.assertEqual(powering.call_count,2)
        with ctx.workprec(256):
            for region in regions.values():
                for value,expected in zip(region,(2,1,1,2)):
                    self.assertGreaterEqual(value,arb(expected))
                    self.assertLess(value,arb(expected)+arb(2)**-100)

    def test_saved_replay_uses_exact_literals_without_screening(self):
        self.assertEqual({r:len(w['markers']) for r,w in SAVED_RUN_WITNESSES.items()},
                         {128:9,512:19,1024:22})
        with patch('hd_shifted_swap_potential.run_bank_receipt',return_value='receipt') as replay:
            for weight in SAVED_RUN_WITNESSES:
                self.assertEqual(checkpoint_receipt(weight,precision=384),'receipt')
                args,kwargs = replay.call_args
                self.assertEqual(args[0],weight)
                self.assertTrue(all(isinstance(x,F) for x in args[1:3]))
                self.assertTrue(all(isinstance(x,F) and x>0 for row in args[3] for x in row))
                self.assertTrue(all(isinstance(x,F) and x>0 for x in args[4]))
                self.assertEqual(kwargs['precision'],384)


if __name__ == '__main__':
    unittest.main()
