"""Exact counting and positivity checks, not a distance certificate."""
from fractions import Fraction as F
from math import comb
import unittest
import hd_field128_dense_gamma_exclusions_v43 as candidate


class Tests(unittest.TestCase):
    def test_partition_count(self):
        # Independent closed sum over v=t+unused; not merely the j32 shape.
        row = candidate.exclusion_count(62)
        count = sum(comb(64, v)*(2**v-1)*sum(comb(64-v, e) for e in range(16-v, 32))
                    for v in (1, 2))
        self.assertEqual(row['ordered_partitions'], count)
        self.assertIn([32, 31, 31, comb(64, 32)*32**2], row['shapes'])
        self.assertEqual(F(row['bad_gamma_upper']),
                         F(16*count+8*comb(64, 32), (1 << 128)-(1 << 64)))

    def test_targeted_rank_and_nested_counts(self):
        last = set()
        for cutoff in (62, 61, 60, 59, 58):
            row = candidate.exclusion_count(cutoff)
            shapes = {tuple(x[:3]) for x in row['shapes']}
            self.assertTrue(last <= shapes)
            for j, c, e in shapes:
                t = j-c
                self.assertTrue(1 <= t <= 64-cutoff)
                self.assertTrue(e < 32 and c+e >= cutoff)
                self.assertTrue(min(c, 64-c, 32-t) < 32)
                self.assertEqual(candidate.prior.rank_floor(j, c, e), t-1)
            last = shapes


if __name__ == '__main__':
    unittest.main()
