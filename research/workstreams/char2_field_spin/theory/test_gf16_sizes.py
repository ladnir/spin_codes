"""Geometry, source-tilt, and new-state transfer checks for GF16 size work."""

from collections import Counter
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction as F
import math
import unittest

from batched_transfer import multiply
from gf16_certificate import interval_local_envelope as baseline_local
from gf16_certificate import upper_regions as baseline_regions
from gf16_sizes import (geometry, interval_local_envelope, singleton_weights,
                        source_eta, upper_regions, upper_tilted_source, validate_witnesses,
                        middle_upper, configuration, verify_middle_entry)
from gf16_size_witnesses import K18_S6
from q1_batched_screen import fixed_support_multiplicity
from subfield_certificate import Interval


class GF16SizeTests(unittest.TestCase):
    def test_geometry_and_thresholds(self):
        self.assertEqual(geometry(16), (8192, 131072, 13107))
        self.assertEqual(geometry(18), (32768, 524288, 52428))
        self.assertEqual(geometry(20), (131072, 2097152, 209715))
        with self.assertRaises(ValueError):
            geometry(17)

    def test_preserved_k16_provider(self):
        z = Interval.rational(17, 19)
        previous, new = baseline_local(z), interval_local_envelope(z, 4)
        for w in range(17):
            for row in range(2):
                for col in range(2):
                    a, b = previous[w][row][col], new[w][row][col]
                    self.assertLessEqual(max(a.lo, b.lo), min(a.hi, b.hi))
                    self.assertLess(abs(F(a.hi)-F(b.hi)), F(1, 10**70))
        a, b = baseline_regions(z, 2), upper_regions(z, 2, 8192, 4)
        for left, right in zip(a, b):
            for x, y in zip(left, right):
                self.assertLess(abs(F(x)-F(y)), max(F(x), F(y), F(1, 10**100))*F(1, 10**70))

    def test_all_state_sizes_have_ca_zero_and_direct_singleton_counts(self):
        mul = lambda a, b: multiply(a, b, 0x13)
        for state in range(1, 9):
            rows = []
            for j in range(16):
                powers = [1]
                for _ in range(state-1):
                    powers.append(mul(powers[-1], j))
                rows.append(powers)
            for left in range(state):
                for right in range(state):
                    total = 0
                    for row in rows:
                        total ^= mul(row[left], row[right])
                    self.assertEqual(total, 0)
            direct = []
            for active in range(16):
                count = 0
                for out in range(16):
                    value = int(active == out)
                    for p in range(state):
                        value ^= mul(rows[active][p], rows[out][p])
                    count += value != 0
                direct.append(count)
            self.assertEqual(tuple(direct), singleton_weights(state))
        self.assertEqual(Counter(singleton_weights(4)), Counter({15: 13, 14: 2, 16: 1}))
        with self.assertRaises(ValueError):
            singleton_weights(9)

    def test_affine_source_envelope_exactly(self):
        for rho in [F(871, 1000), F(1), F(3), F(567, 100)]:
            eta = source_eta(rho)
            for w in range(9, 17):
                actual = fixed_support_multiplicity(16, 8, w, 16)
                self.assertLessEqual(actual, eta*rho**w*16**(w-8))

    def test_tilted_source_rounds_up_against_rationals(self):
        regional = [[F(1, 7), F(2, 9), F(1, 11), F(2, 5)],
                    [F(1, 17), F(1, 13), F(2, 23), F(1, 7)],
                    [F(1, 19), F(1, 21), F(1, 29), F(1, 13)]]
        upper = [[Interval.rational(x.numerator, x.denominator).hi for x in row] for row in regional]
        rho = F(871, 1000)
        matrix = [sum(math.comb(2, h)*16**h*rho**h*regional[h][j] for h in range(3)) for j in range(4)]
        vector = [F(1), F(0)]
        for _ in range(16):
            vector = [vector[0]*matrix[0]+vector[1]*matrix[2],
                      vector[0]*matrix[1]+vector[1]*matrix[3]]
        expected = sum(vector)*source_eta(rho)**2/16**16
        actual = F(upper_tilted_source(upper, 2, rho))
        self.assertGreaterEqual(actual, expected)
        self.assertLess(actual-expected, expected*F(1, 10**70))

    def test_fixed_witness_coverage_and_mutations(self):
        self.assertEqual(validate_witnesses(K18_S6), 32768)
        self.assertEqual(configuration(18, 6)['bad_weight_threshold'], 52428)
        self.assertEqual(configuration(18, 6)['refresh'], 'uniform-gl6-pre-output')
        variants = []
        bad = deepcopy(K18_S6)
        bad['sparse'] = bad['sparse'][1:]
        variants.append(bad)
        bad = deepcopy(K18_S6)
        first = bad['middle'][0]
        bad['middle'][0] = (first[0]+1,)+first[1:]
        variants.append(bad)
        bad = deepcopy(K18_S6)
        bad['dense_first'] += 1
        variants.append(bad)
        bad = deepcopy(K18_S6)
        first = bad['sparse'][0]
        bad['sparse'][0] = (first[0], first[1], '1')
        variants.append(bad)
        for bad in variants:
            with self.assertRaises(ValueError):
                validate_witnesses(bad)

    def test_middle_directed_witness_is_finite(self):
        first, last, x, z, rho = K18_S6['middle'][0]
        bound = middle_upper(18, 6, first, F(x), F(z), rho)
        self.assertTrue(bound.hi.is_finite())
        self.assertGreater(bound.lo, 0)
        self.assertLess(F(bound.hi), F(1, 1 << 80))

    def test_middle_displayed_endpoint_tampering(self):
        replay = dict(first=224, last=239, endpoint_upper=['0.125', '0.25'])
        verify_middle_entry(replay, replay, Decimal(4))
        for values in [['0', '0.25'], ['NaN', '0.25'], ['Infinity', '0.25'], ['0.125']]:
            supplied = dict(replay, endpoint_upper=values)
            with self.assertRaises(ValueError):
                verify_middle_entry(supplied, replay, Decimal(4))
        supplied = dict(replay, endpoint_upper=['0.125', '0.5'])
        with self.assertRaises(ValueError):
            verify_middle_entry(supplied, replay, Decimal(4))
        verify_middle_entry(supplied, replay, Decimal(8))


if __name__ == '__main__':
    unittest.main()
