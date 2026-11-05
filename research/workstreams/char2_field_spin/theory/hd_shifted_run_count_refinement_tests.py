"""Exact enumeration checks; no large construction replay."""
from fractions import Fraction as F
from itertools import product
from math import log2
import unittest

import numpy as np
from flint import arb,ctx

from hd_shifted_pairs_count import exact_profiles,circular_run_count
from hd_shifted_run_potential import run_sum_arb
from hd_shifted_swap_potential import full_counts
from hd_shifted_run_count_refinement import (conditional_profile_count,
    exact_profile_count,run_partition_sum_arb,run_partition_log_screen,
    exact_joint_sum_arb,exact_joint_log_screen)


class RunCountRefinementTests(unittest.TestCase):
    def test_last_exit_count_against_indexed_bit_enumeration(self):
        for pairs in range(1,7):
            counts = {}
            for bits in product((0,1),repeat=2*pairs):
                r = sum(bits)
                a = sum(bits[2*i]*bits[2*i+1] for i in range(pairs))
                b = sum(bits[2*i+1]*bits[(2*i+2)%(2*pairs)] for i in range(pairs))
                key = r,a,b,sum(bits[::2])
                counts[key] = counts.get(key,0)+1
            for r in range(2*pairs+1):
                for a in full_counts(pairs,r):
                    for b in full_counts(pairs,r):
                        for even in range(max(0,r-pairs),min(pairs,r)+1):
                            self.assertEqual(conditional_profile_count(pairs,r,a,b,even),
                                             counts.get((r,a,b,even),0),(pairs,r,a,b,even))
                        self.assertEqual(exact_profile_count(pairs,r,a,b),
                                         sum(n for (rr,aa,bb,_),n in counts.items() if (rr,aa,bb)==(r,a,b)))

    def test_conditioning_on_runs_partitions_without_rotations(self):
        for pairs in range(1,7):
            profiles = exact_profiles(pairs)
            for r in range(2*pairs+1):
                for k in range(pairs+1):
                    self.assertEqual(sum(n for (rr,a,b),n in profiles.items() if rr==r and r-a-b==k),
                                     circular_run_count(2*pairs,r,k))

    def test_uniform_allocation_of_run_mass_is_false(self):
        profiles = exact_profiles(4)
        total = circular_run_count(8,4,2)
        self.assertEqual(total,36)
        self.assertEqual([profiles[4,a,2-a] for a in range(3)],[2,32,2])
        self.assertGreater(profiles[4,1,1],total//3)

    def test_q_norm_bounds_exact_correlated_sum_and_tightens_q1(self):
        with ctx.workprec(256):
            for pairs in range(1,6):
                profiles = exact_profiles(pairs)
                for r in range(2*pairs+1):
                    scalar = {f:F(f+3,f+1) for f in full_counts(pairs,r)}
                    regions = {f:(v,0,0,v) for f,v in scalar.items()}
                    exact = sum(n*scalar[a]**3*scalar[b]**3 for (rr,a,b),n in profiles.items() if rr==r)
                    original = run_sum_arb(pairs,r,regions,(1,1),repetitions=6)
                    previous = None
                    for q in (1,2,4,8):
                        value = run_partition_sum_arb(pairs,r,regions,(1,1),q=q,repetitions=6)
                        self.assertGreaterEqual(value,arb(exact.numerator)/exact.denominator)
                        if q==1:
                            self.assertLess(abs(value-original),abs(original)*arb(2)**-180)
                        if previous is not None:
                            self.assertLessEqual(value,previous+abs(previous)*arb(2)**-180)
                        previous = value

    def test_screen_matches_directed_and_maximum_is_a_bound(self):
        with ctx.workprec(256):
            pairs,r = 6,6
            regions = {f:(f+1,0,0,f+1) for f in full_counts(pairs,r)}
            logs = {f:(log2(f+1),-np.inf,-np.inf,log2(f+1)) for f in regions}
            for q in (1,2,4,8):
                exact = run_partition_sum_arb(pairs,r,regions,(1,1),q=q,repetitions=6)
                screen = run_partition_log_screen(pairs,r,logs,0,q=q,repetitions=6)
                self.assertAlmostEqual(screen,float(exact.log()/arb(2).log()),places=10)
            actual = sum(n*(a+1)**3*(b+1)**3 for (rr,a,b),n in exact_profiles(pairs).items() if rr==r)
            self.assertGreaterEqual(run_partition_log_screen(pairs,r,logs,0,q=None,repetitions=6),log2(actual))

    def test_factorized_moments_equal_exact_correlated_scalar_sum(self):
        with ctx.workprec(256):
            for pairs in range(1,8):
                profiles = exact_profiles(pairs)
                for r in range(2*pairs+1):
                    regions = {f:(F(f+2,f+1),0,0,F(f+2,f+1)) for f in full_counts(pairs,r)}
                    exact = sum(n*F(a+2,a+1)**3*F(b+2,b+1)**3
                                for (rr,a,b),n in profiles.items() if rr==r)
                    value = exact_joint_sum_arb(pairs,r,regions,(1,1),repetitions=6)
                    expected = arb(exact.numerator)/exact.denominator
                    self.assertGreaterEqual(value,expected)
                    self.assertLess(abs(value-expected),abs(expected)*arb(2)**-180)
                    logs = {f:(log2(F(f+2,f+1)),-np.inf,-np.inf,log2(F(f+2,f+1))) for f in regions}
                    screen = exact_joint_log_screen(pairs,r,logs,0,repetitions=6)
                    self.assertAlmostEqual(screen,log2(exact),places=10)

    def test_factorized_count_uniform_caps_and_noncommuting_products(self):
        from math import comb
        def mul(x,y):
            return tuple(x[2*i]*y[j]+x[2*i+1]*y[2+j] for i in range(2) for j in range(2))
        with ctx.workprec(256):
            profiles = exact_profiles(5)
            for r in range(11):
                identity = {f:(1,0,0,1) for f in full_counts(5,r)}
                value = exact_joint_sum_arb(5,r,identity,(1,1),repetitions=6)
                self.assertLess(abs(value-comb(10,r)),arb(2)**-180)
                regions = {f:tuple(map(F,(1+f,2,3,4+f))) for f in full_counts(5,r)}
                exact = F(0)
                for (rr,a,b),count in profiles.items():
                    if rr!=r:
                        continue
                    pair = mul(regions[a],regions[b])
                    run = mul(mul(pair,pair),pair)
                    exact += count*(run[0]+run[1])
                for vector in ((2,3),(3,2)):
                    value = exact_joint_sum_arb(5,r,regions,vector,repetitions=6)
                    self.assertGreaterEqual(value,arb(exact.numerator)/exact.denominator)


if __name__=='__main__':
    unittest.main()
