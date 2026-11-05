"""Exact rational checks of the directed fixed-A matrix polynomial."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import unittest

from flint import arb, ctx, fmpq

import hd_xor_fixed_mask_floor_certificate as certificate
from hd_xor_fixed_mask_floor import comparison_counts
from hd_xor_paired_gate import CONSTRUCTION


def multiply(left, right):
    return tuple(sum(left[2*i+k]*right[2*k+j] for k in range(2))
                 for i in range(2) for j in range(2))


def directed(matrix):
    return tuple(arb(value.numerator)/value.denominator for value in matrix)


IDENTITY = (F(1), F(0), F(0), F(1))


class FixedMaskFloorCertificateTests(unittest.TestCase):
    def assertRationalInInterval(self, interval, value):
        rational = fmpq(value.numerator, value.denominator)
        # arb.contains(fmpq) first rounds the rational to a ball and can
        # reject that larger ball although the exact rational lies inside.
        self.assertLessEqual(interval.lower().fmpq(), rational)
        self.assertGreaterEqual(interval.upper().fmpq(), rational)

    def assertEncloses(self, actual, expected):
        for interval, value in zip(actual, expected):
            self.assertRationalInInterval(interval, value)

    def test_labelled_mask_products_and_exact_split(self):
        # Noncommuting, non-dyadic entries expose both order and rounding.
        regions = tuple(tuple(tuple(F(1+7*i+3*b+j, 29+j)
                                    for j in range(4)) for b in range(2))
                        for i in range(4))
        exact = []
        for mask in range(16):
            value = IDENTITY
            for index, table in enumerate(regions):
                value = multiply(value, table[(mask >> index) & 1])
            exact.append(value)
        with ctx.workprec(192):
            operators = certificate.mask_operators_arb(
                tuple(tuple(directed(matrix) for matrix in table) for table in regions))
            for actual, expected in zip(operators, exact):
                self.assertEncloses(actual, expected)
            for active_mask in range(1, 16):
                counts = comparison_counts(active_mask, group=4, q=4)
                normalizer = 4**active_mask.bit_count()
                empty, nonempty = certificate.fixed_column_split_arb(
                    operators, active_mask, group=4, q=4)
                self.assertEncloses(empty, tuple(x/normalizer for x in exact[0]))
                self.assertEncloses(nonempty, tuple(
                    sum(counts[mask]*exact[mask][entry] for mask in range(1, 16))/normalizer
                    for entry in range(4)))

    def test_positive_convolution_matches_every_exact_word(self):
        empty = (F(1, 3), F(1, 7), F(2, 5), F(1, 11))
        nonempty = (F(1, 13), F(2, 3), F(1, 17), F(1, 5))
        self.assertNotEqual(multiply(empty, nonempty), multiply(nonempty, empty))
        for columns in range(6):
            expected = [[F(0) for _ in range(4)] for _ in range(columns+1)]
            for word in product((0, 1), repeat=columns):
                value = IDENTITY
                for flag in word:
                    value = multiply(value, (empty, nonempty)[flag])
                for entry in range(4):
                    expected[sum(word)][entry] += value[entry]
            with ctx.workprec(192):
                actual = certificate.column_polynomial_arb(
                    directed(empty), directed(nonempty), columns)
                for interval, exact in zip(actual, expected):
                    self.assertEncloses(interval, exact)
                for floor in range(columns+1):
                    tail = certificate.tail_mass_arb(actual, floor)
                    value = sum(row[0]+row[1] for row in expected[floor:])
                    self.assertRationalInInterval(tail, value)

    def test_zero_column_elimination_and_exact_probability_mass(self):
        with ctx.workprec(192):
            operators = [directed(IDENTITY)]*16
            empty, nonempty = certificate.fixed_column_split_arb(operators, 7, 4, 4)
            values = certificate.column_polynomial_arb(empty, nonempty, 4)
            self.assertEncloses(values[0], tuple(x/F(64**4) for x in IDENTITY))
            self.assertTrue(certificate.tail_mass_arb(values, 0).contains(1))
            expected = 1-F(1, 64**4)
            self.assertRationalInInterval(certificate.tail_mass_arb(values, 1), expected)

    def test_fresh_dual_receipts_agree_on_all219_ceilings_and_scope(self):
        directory = Path(__file__).resolve().parent
        receipts = [json.loads((directory/f'hd_xor_fixed_mask_floor_r1_receipt{bits}.json').read_text())
                    for bits in (256, 384)]
        source_hash = sha256(Path(certificate.__file__).read_bytes()).hexdigest()
        for bits, row in zip((256, 384), receipts):
            self.assertEqual(row['precision_bits'], bits)
            self.assertEqual(row['source_sha256'], source_hash)
            self.assertEqual(row['construction'], CONSTRUCTION)
            self.assertEqual(row['mark'], certificate.SAVED_MARK)
            self.assertEqual(row['active_groups'], 1)
            self.assertEqual(row['source_supports'], 8192)
            self.assertEqual(row['source_matching_counts'], [0, 0])
            self.assertEqual(row['active_constituent_counts'], list(range(3, 9)))
            self.assertEqual(row['nonempty_column_floor'], 17)
            self.assertEqual(row['projective_normalizations'], 1)
            self.assertTrue(row['directed'])
            self.assertFalse(row['whole_code_certificate'])
            self.assertFalse(row['comparison_is_actual_rs_law'])
            self.assertFalse(row['extra_row_set_binomial_factor'])
            self.assertEqual([item['active_row_mask'] for item in row['entries']],
                             [mask for mask in range(256) if mask.bit_count() >= 3])
            self.assertEqual(len(row['entries']), 219)
            self.assertLess(F(row['failure_log2_upper']), -40)
        for key in ('entries', 'arity_summaries', 'failure_log2_upper'):
            self.assertEqual(receipts[0][key], receipts[1][key])
        self.assertEqual(receipts[0]['failure_log2_upper'], '-82398007651003/500000000000')


if __name__ == '__main__':
    unittest.main()
