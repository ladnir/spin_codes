"""Exact support enumeration independently checks the run convolution bound."""
from fractions import Fraction as F
import unittest

import numpy as np
from flint import arb,ctx

from hd_shifted_pairs_count import exact_profiles,circular_run_count,profile_count_upper
from hd_shifted_swap_potential import full_counts
from hd_shifted_run_potential import run_sum_arb,run_log_screen,log_joint_counts


class RunPotentialTests(unittest.TestCase):
    def test_every_joint_class_is_bounded_without_independence(self):
        for pairs in range(1,7):
            for (weight,a,b),count in exact_profiles(pairs).items():
                self.assertLessEqual(count,circular_run_count(2*pairs,weight,weight-a-b))

    def test_convolution_majorizes_exact_correlated_support_sum(self):
        with ctx.workprec(192):
            for pairs in range(1,7):
                profiles = exact_profiles(pairs)
                for weight in range(2*pairs+1):
                    scalars = {f:F(f+3,f+1) for f in full_counts(pairs,weight)}
                    regions = {f:(q,0,0,q) for f,q in scalars.items()}
                    exact = sum(n*scalars[a]**3*scalars[b]**3
                                for (r,a,b),n in profiles.items() if r==weight)
                    upper = run_sum_arb(pairs,weight,regions,(1,1),repetitions=6)
                    self.assertGreaterEqual(upper,arb(exact.numerator)/exact.denominator)

    def test_float_matches_directed_scalar_fixture(self):
        pairs,weight = 6,6
        regions = {f:(f+1,0,0,f+1) for f in full_counts(pairs,weight)}
        logs = {f:(np.log2(f+1),-np.inf,-np.inf,np.log2(f+1)) for f in regions}
        with ctx.workprec(192):
            exact = run_sum_arb(pairs,weight,regions,(1,1),repetitions=6)
            expected = float(exact.log()/arb(2).log())
        self.assertAlmostEqual(run_log_screen(pairs,weight,logs,0,repetitions=6),expected)

    def test_noncommuting_products_with_both_potential_orderings(self):
        def multiply(left,right):
            a,b,c,d = left
            e,f,g,h = right
            return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h
        with ctx.workprec(192):
            profiles = exact_profiles(5)
            for weight in range(1,10):
                regions = {f:tuple(map(F,(1+f,2,3+f,2+2*f))) for f in full_counts(5,weight)}
                exact = F(0)
                for (r,a,b),count in profiles.items():
                    if r!=weight:
                        continue
                    pair = multiply(regions[a],regions[b])
                    run = multiply(multiply(pair,pair),pair)
                    exact += count*(run[0]+run[1])
                for vector in ((3,2),(2,3)):
                    bound = run_sum_arb(5,weight,regions,vector,repetitions=6)
                    self.assertGreaterEqual(bound,arb(exact.numerator)/exact.denominator)

    def test_full_full_nonconstant_profile_contributes_zero(self):
        # Only the endpoint contraction is nonzero. Both matchings full
        # would force a constant support, impossible at this weight.
        regions = {f:(0,0,0,0) for f in full_counts(6,6)}
        regions[3] = (1,0,0,1)
        self.assertEqual(run_sum_arb(6,6,regions,(1,1)),arb(0))

    def test_bad_repetition_count_is_rejected(self):
        for repetitions in (0,3,4.5):
            with self.assertRaises(ValueError):
                run_sum_arb(4,4,{},(1,1),repetitions=repetitions)

    def test_joint_count_screen_matches_integer_reference(self):
        from math import log2
        for pairs in range(1,7):
            for weight in range(2*pairs+1):
                counts = log_joint_counts(pairs,weight)
                for i,a in enumerate(full_counts(pairs,weight)):
                    for j,b in enumerate(full_counts(pairs,weight)):
                        value = profile_count_upper(pairs,weight,a,b)
                        self.assertEqual(counts[i,j],log2(value) if value else -np.inf)


if __name__=='__main__':
    unittest.main()
