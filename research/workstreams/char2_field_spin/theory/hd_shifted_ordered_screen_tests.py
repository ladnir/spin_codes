"""Ordered matrix screens against exact small positive matrix products."""
from fractions import Fraction as F
from math import log2
import unittest

import numpy as np

from hd_shifted_ordered_screen import ordered_log_screen
from hd_shifted_pairs_count import circular_run_count,profile_count_upper
from hd_shifted_run_potential import run_log_screen
from hd_shifted_swap_potential import full_counts


def multiply(left,right):
    a,b,c,d = left
    e,f,g,h = right
    return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h


class OrderedScreenTests(unittest.TestCase):
    def test_all_pairs_match_exact_reference_with_both_counts(self):
        for pairs in range(2,6):
            for weight in range(1,2*pairs):
                regions = {f:tuple(map(F,(f+1,2,3+f,2*f+1))) for f in full_counts(pairs,weight)}
                logs = {f:np.array(list(map(log2,values))) for f,values in regions.items()}
                for joint in (False,True):
                    value = F(0)
                    for a in regions:
                        for b in regions:
                            pair = multiply(regions[a],regions[b])
                            run = multiply(multiply(pair,pair),pair)
                            count = (profile_count_upper(pairs,weight,a,b) if joint else
                                     circular_run_count(2*pairs,weight,weight-a-b))
                            value += count*(run[0]+run[1])
                    result = ordered_log_screen(pairs,weight,logs,repetitions=6,joint_counts=joint)
                    self.assertAlmostEqual(result['log2_support_sum_screen'],log2(value),places=11)
                    self.assertFalse(result['directed'])

    def test_common_potential_majorizes_ordered_products(self):
        regions = {f:np.log2(np.array((f+1,2,3+f,2*f+1))) for f in full_counts(6,5)}
        direct = ordered_log_screen(6,5,regions,repetitions=8)['log2_support_sum_screen']
        for ratio in (-2,0,2):
            self.assertGreaterEqual(run_log_screen(6,5,regions,ratio,repetitions=8),direct)

    def test_zero_entries_and_constant_support(self):
        for weight in (0,8):
            regions = {weight//2:np.array((0,-np.inf,-np.inf,0))}
            self.assertEqual(ordered_log_screen(4,weight,regions)['log2_support_sum_screen'],0)

    def test_odd_or_zero_repetition_count_rejected(self):
        for repetitions in (0,3,4.5):
            with self.assertRaises(ValueError):
                ordered_log_screen(4,4,{},repetitions=repetitions)


if __name__=='__main__':
    unittest.main()
