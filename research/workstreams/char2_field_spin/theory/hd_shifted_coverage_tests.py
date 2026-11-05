import unittest
from fractions import Fraction as F
from itertools import product

import numpy as np
from flint import arb,ctx

from hd_shifted_alias_certificate import interval_fft
from hd_shifted_coverage import (RegionalAliases,SupportCounts,certainly_empty,cover)
from hd_shifted_general_count import ProfileBox,box_count_upper,exact_profiles


def points(ranges):
    return set(product(*(range(lo,hi+1) for lo,hi in ranges)))


class CoverageTests(unittest.TestCase):
    def test_cached_aliases_include_wraparound_and_restore_similarity(self):
        with ctx.workprec(256):
            raw = np.arange(4*4**3,dtype=float).reshape(4,4,4,4)/256
            zeros = np.zeros_like(raw)
            powered = interval_fft(((raw.copy(),raw.copy()),(zeros.copy(),zeros.copy())),(1,2,3))
            atlas = RegionalAliases.from_power(powered,(arb(1),arb(2)),arb(1))
            for ranges in (((1,2),(2,3),(0,2)),((3,5),(0,8),(4,4))):
                actual = atlas.numerator(ranges)
                residues = [np.unique(np.arange(lo,hi+1)%4) for lo,hi in ranges]
                exact = np.max(raw[np.ix_(np.arange(4),*residues)],axis=(1,2,3))
                for i in range(4):
                    expected = arb(float(exact[i]))*arb((1,2)[i//2])/((1,2)[i%2])
                    self.assertTrue(actual[i]>=expected)
            self.assertFalse(atlas.upper.flags.writeable)

    def test_cached_counts_match_uncached_reference(self):
        markers = (F(1,2),F(3,2),F(2,3),F(3))
        cache = SupportCounts(5,markers)
        for ranges in (((0,5),)*4,((1,3),(2,4),(0,2),(1,3)),((2,2),(2,2),(1,1),(1,1))):
            box = ProfileBox(*ranges)
            self.assertEqual(cache.upper(box),box_count_upper(5,box,markers))

    def test_empty_filter_never_removes_small_feasible_profiles(self):
        for pairs in range(1,6):
            for profile in exact_profiles(pairs):
                box = ProfileBox(*((n,n) for n in profile))
                self.assertFalse(certainly_empty(pairs,box),profile)
                larger = ProfileBox(*((max(0,n-1),min(pairs,n+1)) for n in profile))
                self.assertFalse(certainly_empty(pairs,larger),profile)

    @staticmethod
    def evaluate(box):
        count = len(points(box.ranges))
        # A synthetic checker that only closes sufficiently small leaves.
        return dict(box=box.ranges,log2_upper_rational=str(-100 if count<=4 else 0),directed=True)

    def check_partition(self,result):
        seen = set()
        leaves = [item['box'] for item in result['accepted']]+result['empty']+result['unresolved']
        for ranges in leaves:
            current = points(ranges)
            self.assertFalse(seen&current)
            seen |= current
        self.assertEqual(seen,points(result['root']))

    def test_disjoint_complete_cover_and_directed_union(self):
        root = ProfileBox((1,3),(0,3),(0,3),(0,3))
        result = cover(root,self.evaluate,pairs=3,max_evaluations=1000)
        self.check_partition(result)
        self.assertTrue(result['complete_local_cover'])
        self.assertFalse(result['whole_code_certificate'])
        with ctx.workprec(512):
            ceiling = F(result['accepted_union_log2_upper'])
            self.assertTrue(arb(2)**(arb(ceiling.numerator)/ceiling.denominator)
                            >= len(result['accepted'])*arb(2)**-100)

    def test_budget_preserves_every_unresolved_box(self):
        root = ProfileBox((1,3),(0,3),(0,3),(0,3))
        result = cover(root,self.evaluate,pairs=3,max_evaluations=5)
        self.check_partition(result)
        self.assertFalse(result['complete_local_cover'])
        self.assertEqual(result['evaluations'],5)

    def test_nonclosing_singletons_remain_explicit(self):
        root = ProfileBox((1,1),(0,0),(0,0),(0,0))
        result = cover(root,lambda box:dict(box=box.ranges,log2_upper_rational='0',directed=True),pairs=3)
        self.assertEqual(result['unresolved'],[root.ranges])
        self.assertFalse(result['complete_local_cover'])

    def test_rejects_float_receipt_mismatched_box_and_zero_message(self):
        root = ProfileBox((1,1),(0,0),(0,0),(0,0))
        for bad in (dict(box=root.ranges,log2_upper_rational='-100',directed=False),
                    dict(box=((2,2),)*4,log2_upper_rational='-100',directed=True)):
            with self.assertRaises(ValueError):
                cover(root,lambda box:bad,pairs=3)
        with self.assertRaises(ValueError):
            cover(ProfileBox((0,1),(0,1),(0,1),(0,1)),self.evaluate,pairs=3)


if __name__=='__main__':
    unittest.main()
