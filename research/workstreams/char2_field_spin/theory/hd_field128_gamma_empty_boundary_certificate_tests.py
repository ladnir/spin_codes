"""Tiny independent directed-polynomial tests; no K20 bank or gate."""

from fractions import Fraction as F
from itertools import product
from math import comb
from unittest.mock import patch
import unittest

from flint import arb, ctx

import hd_field128_gamma_empty_boundary_certificate as mode


def multiply(a, b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2)) for i in range(2) for j in range(2))


def poly_product(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            z = multiply(x, y)
            out[i+j] = tuple(v+w for v, w in zip(out.get(i+j, (F(0),)*4), z))
    return out


class BoundaryTests(unittest.TestCase):
    def test_exact_rational_oracle_tilt_truncation_and_reversal(self):
        with ctx.workprec(256):
            slots, power = 2, 3
            raw = {j: tuple(F(x, 8) for x in (j+1, 2*j+1, j+2, 3*j+1)) for j in range(5)}
            table = {j: tuple(mode.base.rational(x) for x in row) for j, row in raw.items()}
            paired = {}
            for i in range(3):
                for j in range(3):
                    row = multiply(raw[2+i], raw[2+j])
                    count = comb(2, i)*comb(2, j)
                    paired[i+j] = tuple(v+count*x for v, x in zip(paired.get(i+j, (F(0),)*4), row))
            exact = {0: (F(1), F(0), F(0), F(1))}
            for unused in range(power):
                exact = poly_product(exact, paired)
            for full in range(13):
                expected = tuple(x/comb(12, full) for x in exact[full])
                for marker, reverse in ((F(1), False), (F(1, 3), False), (F(1, 3), True)):
                    actual, meta = mode.regional_boundary(table, F(1), full, slots=slots,
                        paired_power=power, marker=marker, reverse=reverse)
                    for got, want in zip(actual, expected):
                        target = mode.base.rational(want)
                        self.assertTrue(got >= target.upper())
                        self.assertTrue(abs(got-target) < arb(2)**-210)
                    self.assertEqual(meta['extracted_degree'], min(full, 12-full) if reverse else full)

    def test_fractional_raw_product_precedes_alpha(self):
        with ctx.workprec(256):
            table = {j: tuple(arb(v) for v in (j+1, j+2, 2*j+1, j+3)) for j in range(5)}
            actual = mode.top_coefficients(table, F(1, 2), slots=2)
            for degree in range(5):
                expected = [arb(0) for unused in range(4)]
                wrong = [arb(0) for unused in range(4)]
                for i in range(3):
                    j = degree-i
                    if not 0 <= j <= 2:
                        continue
                    count = comb(2, i)*comb(2, j)
                    raw = multiply(table[2+i], table[2+j])
                    bad = multiply(tuple(x.sqrt() for x in table[2+i]), tuple(x.sqrt() for x in table[2+j]))
                    for entry in range(4):
                        expected[entry] += count*raw[entry].sqrt()
                        wrong[entry] += count*bad[entry]
                for entry in range(4):
                    got = actual[entry][degree]
                    self.assertTrue(got >= expected[entry].lower())
                    self.assertTrue(abs(got-expected[entry]) < arb(2)**-200)
                    self.assertTrue(abs(got-wrong[entry]) > arb('0.1'))

    def test_profile_transfer_and_boundary_count_tiny_words(self):
        with ctx.workprec(256):
            matrix = tuple(map(F, ('1/2', '1/4', '1/8', '1/2')))
            powered = (F(1), F(0), F(0), F(1))
            for unused in range(6):
                powered = multiply(powered, matrix)
            support, transfer = mode.profile_value(tuple(mode.base.rational(x) for x in matrix), 7, repetitions=6)
            self.assertEqual(transfer, mode.base.rational(powered[0]+powered[1]))
            self.assertEqual(support, mode.base.rational(7*(powered[0]+powered[1])))
        for pairs in range(2, 5):
            for weight in range(pairs, 2*pairs+1):
                f, actual = weight-pairs, 0
                for word in product((0, 1), repeat=2*pairs):
                    if sum(word) != weight:
                        continue
                    a = sum(word[2*i]*word[2*i+1] for i in range(pairs))
                    b = sum(word[2*i+1]*word[(2*i+2)%(2*pairs)] for i in range(pairs))
                    actual += a == b == f
                count = mode.midpoint.MidpointCounts(pairs, weight, max_singletons=64).count(f, f)
                self.assertGreaterEqual(count, actual)
                if weight == pairs:
                    self.assertEqual(actual, 2)

    def test_public_receipt_scope_prefix_and_incomplete_failure(self):
        class Counts:
            def __init__(self, pairs, weight, **unused):
                self.pairs, self.weight = pairs, weight
            def count(self, a, b):
                self_args = (self.pairs, self.weight, a, b)
                if self_args != (4096, 5120, 1024, 1024):
                    raise AssertionError(self_args)
                return 7
            def diagnostics(self):
                return dict(synthetic_count_fixture=True)
        matrix = tuple(arb(x) for x in ('0.5', '0.25', '0.125', '0.5'))
        with patch.object(mode.provider.frozen.sparse, 'packet_table', return_value={}) as packet, \
                patch.object(mode, 'regional_boundary', return_value=(matrix, {'synthetic': True})), \
                patch.object(mode.midpoint, 'MidpointCounts', Counts):
            result = mode.boundary_receipt('1/2', '1/8', 5120, max_seconds=30)
        packet.assert_called_once()
        self.assertEqual(result['count_upper'], 7)
        self.assertEqual(result['profile_pairs_evaluated'], 1)
        self.assertTrue(result['complete_boundary_profile'])
        self.assertFalse(result['complete_selected_weight'])
        self.assertFalse(result['covers_all_supports_at_weight'])
        self.assertFalse(result['whole_code_certificate'])
        with ctx.workprec(256):
            alpha, z = F(1, 8), F(1, 2)
            prefix = mode.base.rational(alpha)*(5120*mode.tail.log2arb(mode.base.rational(mode.base.outer_density()))
                -mode.tail.log2arb(arb((1 << 128)-1))-int(F(2, 5)*(1 << 21))*mode.tail.log2arb(mode.base.rational(z)))
            self.assertTrue(arb(result['original_weight_outer_prefix_enclosure']).overlaps(prefix))
        with patch.object(mode.provider.frozen.sparse, 'packet_table', side_effect=TimeoutError('tiny fixture')):
            incomplete = mode.boundary_receipt('1/2', '1/8', 5120, max_seconds=30)
        self.assertFalse(incomplete['complete_boundary_profile'])
        self.assertFalse(incomplete['directed'])
        self.assertIsNone(incomplete['log2_upper_rational'])
        self.assertEqual(incomplete['profile_pairs_evaluated'], 0)

    def test_guards_and_exit_source_check(self):
        for z, alpha, weight in ((.5, '1/8', 5120), ('1/2', '0', 5120), ('1/2', '1/8', 4095)):
            with self.assertRaises(ValueError):
                mode.boundary_receipt(z, alpha, weight)
        with patch.object(mode, 'check_source_pins', side_effect=[None, RuntimeError('changed')]), \
                patch.object(mode.provider.frozen.sparse, 'packet_table', side_effect=TimeoutError('tiny fixture')):
            with self.assertRaisesRegex(RuntimeError, 'changed'):
                mode.boundary_receipt('1/2', '1/8', 5120)


if __name__ == '__main__':
    unittest.main()
