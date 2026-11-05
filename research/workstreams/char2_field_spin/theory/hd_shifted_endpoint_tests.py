"""Exact small swap laws, alpha chronology, quiet enumerator, and directed replay."""
from fractions import Fraction as F
from itertools import product
import unittest

from flint import arb, ctx

from hd_common_refresh_packets import local_moments
from hd_shifted_alias_certificate import rational
from hd_shifted_endpoint import (swap_occupancy_counts,fractional_average,
                                fixed_lane_quiet_coefficient,fixed_lane_obstruction,
                                swap_endpoint_receipt)


class EndpointTests(unittest.TestCase):
    def test_fair_swap_counts_exhaustive(self):
        for half in range(1,5):
            expected = {}
            for bits in product((0,1),repeat=2*half):
                a,b = sum(bits[:half]),sum(bits[half:])
                occupancy = tuple(sorted((a,half-a,b,half-b)))
                expected[occupancy] = expected.get(occupancy,0)+1
            self.assertEqual(dict(swap_occupancy_counts(half)),expected)
        self.assertEqual(sum(n for _,n in swap_occupancy_counts()),2**32)

    def test_only_two_alternating_nonempty_independent_sets(self):
        for pairs in range(1,7):
            supports = []
            for bits in product((0,1),repeat=2*pairs):
                if sum(bits)==pairs and all(not(bits[i] and bits[(i+1)%(2*pairs)])
                                           for i in range(2*pairs)):
                    supports.append(bits)
            self.assertEqual(set(supports),{(0,1)*pairs,(1,0)*pairs})

    def test_fractional_power_precedes_average(self):
        with ctx.workprec(256):
            result = fractional_average(((F(1,2),(F(1,4),)*4),
                                         (F(1,2),(F(9,4),)*4)),F(1,2))
            for entry in result:
                self.assertTrue(entry.contains(1))
                self.assertTrue(entry < rational(F(5,4)).sqrt())
            with self.assertRaises(ValueError):
                fractional_average(((F(1,3),(1,1,1,1)),),F(1,2))

    def test_quiet_enumerator_matches_exact_packet_moment(self):
        for z in (F(1,32),F(2,5),F(3,4)):
            polynomial = sum(fixed_lane_quiet_coefficient(w)*z**w for w in range(33))
            self.assertEqual(polynomial,local_moments(z,q=256,h=4)[0][16]**2)

    def test_uniform_comparator_obstruction_is_directed_and_positive(self):
        first = fixed_lane_obstruction(precision=256)
        second = fixed_lane_obstruction(precision=384)
        self.assertEqual(first['log2_lower_rational'],second['log2_lower_rational'])
        self.assertGreater(F(first['log2_lower_rational']),324696)
        self.assertFalse(first['actual_distance_refutation'])

    def test_256_and_384_bit_swap_receipts(self):
        first = swap_endpoint_receipt(precision=256)
        second = swap_endpoint_receipt(precision=384)
        self.assertEqual(first['log2_upper_rational'],second['log2_upper_rational'])
        self.assertEqual(F(first['log2_upper_rational']),F(-95969876597769549,10**12))
        self.assertEqual(first['exact_support_count'],2)
        self.assertFalse(first['whole_code_certificate'])
        with self.assertRaises(TypeError):
            swap_endpoint_receipt(.03125,F(1,4))


if __name__ == '__main__':
    unittest.main()
