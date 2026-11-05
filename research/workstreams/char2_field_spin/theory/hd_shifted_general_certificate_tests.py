"""Exact input and composition checks for general directed alias receipts."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import factorial
import unittest
from unittest.mock import patch

from flint import arb, ctx

from hd_shifted_alias_certificate import odd_coefficients, rational
from hd_shifted_general_certificate import LIMIT, general_coefficients, receipt, truncate_coefficients
from hd_shifted_general_count import joint_count_upper, matching_histograms


class ConstantPackets(dict):
    def items(self):
        return ((occupancy, (F(1),)*4)
                for occupancy in combinations_with_replacement(range(17), 4))


class GeneralCertificateTests(unittest.TestCase):
    def test_three_distinct_rational_markers(self):
        markers = (F(1, 8), F(1, 16), F(1, 32))
        with ctx.workprec(192), patch('hd_shifted_alias_certificate.packet_table', return_value=ConstantPackets()):
            coefficients = general_coefficients(F(1, 2), F(1), markers)
            for a, b, c in ((0, 0, 0), (1, 2, 3), (7, 9, 16), (32, 0, 0), (0, 32, 0)):
                count = factorial(32)//(factorial(32-a-b-c)*factorial(a)*factorial(b)*factorial(c))
                exact = rational(count*markers[0]**a*markers[1]**b*markers[2]**c)
                for row in coefficients:
                    self.assertTrue(row[a*33**2+b*33+c].contains(exact))

    def test_equal_marker_api_preserves_previous_coefficients(self):
        with ctx.workprec(128), patch('hd_shifted_alias_certificate.packet_table', return_value=ConstantPackets()):
            old = odd_coefficients(F(1, 2), F(1, 2), F(1, 8), F(1, 16))
            new = general_coefficients(F(1, 2), F(1, 2), (F(1, 8), F(1, 8), F(1, 16)))
        for first, second in zip(old, new):
            for left, right in zip(first, second):
                self.assertEqual(str(left), str(right))

    def test_truncation_keeps_axes_separate_and_input_unchanged(self):
        coeff = [[arb(i+1) for i in range(27)] for _ in range(4)]
        reduced = truncate_coefficients(coeff, (0, 1, 2))
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    index = a*9+b*3+c
                    self.assertTrue(coeff[0][index] == arb(index+1))
                    self.assertTrue(reduced[0][index] == (arb(index+1) if a == 0 and b <= 1 else arb(0)))

    def test_receipt_uses_both_histograms_and_keeps_region_order(self):
        seen = []
        def region(z, alpha, histogram, markers, **options):
            seen.append((histogram, markers))
            values = (1, 2, 3, 1) if len(seen) == 1 else (2, 0, 1, 3)
            return tuple(arb(value)/2**10000 for value in values), {'directed': True}
        profile = (12, 15, 8, 9)
        markers = ((F(1, 8), F(1, 16), F(1, 32)), (F(1, 7), F(1, 15), F(1, 31)))
        with patch('hd_shifted_general_certificate.regional', side_effect=region):
            actual = receipt(profile, F(1, 2), F(1, 10), markers, support_count='upper')
        self.assertEqual([item[0] for item in seen], list(matching_histograms(4096, *profile)))
        self.assertEqual([item[1] for item in seen], list(markers))
        self.assertEqual(actual['distinct_regional_replays'], 2)
        self.assertTrue(actual['accepted'])
        self.assertFalse(actual['whole_code_certificate'])
        # Compute the ordered power with exact integers, independently of
        # the Arb matrix-power helper used by the production composition.
        base, power = (4, 6, 7, 3), (1, 0, 0, 1)
        for _ in range(128):
            power = tuple(power[2*r]*base[s]+power[2*r+1]*base[2+s]
                          for r in range(2) for s in range(2))
        with ctx.workprec(256):
            log2 = lambda value: value.log()/arb(2).log()
            expected = (log2(arb(joint_count_upper(4096, *profile)))
                        + rational(F(1, 10))*(27*(2048-128*log2(arb(255)))
                                               -log2(arb(255))+LIMIT)
                        + log2(arb(power[0]+power[1]))-2560000)
            self.assertTrue(arb(actual['log2_upper_enclosure']).overlaps(expected))

    def test_equal_regions_reused_without_independence_claim(self):
        with patch('hd_shifted_general_certificate.regional', return_value=((arb(2)**-10000,)*4, {'directed': True})) as region:
            actual = receipt(support_count='upper')
        self.assertEqual(region.call_count, 1)
        self.assertEqual(actual['distinct_regional_replays'], 1)
        self.assertFalse(actual['full_under_either_matching'])


if __name__ == '__main__':
    unittest.main()
