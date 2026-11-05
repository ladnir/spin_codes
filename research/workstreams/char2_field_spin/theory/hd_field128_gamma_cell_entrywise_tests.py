"""Synthetic matrix/marker tests only; no distance gates or bank construction."""

from fractions import Fraction as F
from itertools import product
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_cell_entrywise as c


def synthetic_bank():
    markers = ((F(1, 2), F(1, 3)), (F(2, 3), F(3, 5)))
    return SimpleNamespace(markers=markers,
                           log_markers=tuple(tuple(c.tail.log2arb(c.base.rational(value)) for value in row)
                                             for row in markers),
                           constant=arb(-100), slope=arb(3),
                           numerator=lambda index, rect: tuple(arb(x) for x in ((1, 9, 4, 3), (5, 2, 8, 1))[index]))


def normalizer(weight, full, marker):
    single = weight-2*full
    count = c.base.comb(c.base.PAIRS, full)*c.base.comb(c.base.PAIRS-full, single)
    return count*c.base.rational(marker[0])**single*c.base.rational(marker[1])**full


class EntrywiseCellTests(unittest.TestCase):
    def test_singleton_mixed_markers_match_noncommuting_point_oracle(self):
        with ctx.workprec(256):
            bank, cell = synthetic_bank(), c.Cell('interior', ((9, 9), (2, 2), (3, 3)))
            matrices = (((1, 9, 4, 3), (5, 2, 8, 1)), ((3, 8, 1, 5), (2, 4, 6, 2)))
            def numerator(index, rectangle):
                full = rectangle[1][0]
                return tuple(arb(value)*normalizer(9, full, bank.markers[index])
                             for value in matrices[full-2][index])
            bank.numerator = numerator
            result = c._entrywise_bound(bank, cell, 37, c.base.Budget())
            balanced = c._entrywise_bound(bank, cell, 37, c.base.Budget(), balanced_density=True)
            left, right = tuple(map(arb, (1, 2, 4, 1))), tuple(map(arb, (2, 4, 1, 2)))
            run = c.base.matrix_power(c.base.matrix_product(left, right), 128)
            reverse = c.base.matrix_power(c.base.matrix_product(right, left), 128)
            self.assertFalse((run[0]+run[1]).overlaps(reverse[0]+reverse[1]))
            expected = c.tail.log2arb(37*(run[0]+run[1]))+bank.constant+9*bank.slope
            self.assertGreaterEqual(F(result['log2_upper_rational']), c.base.rational_ceiling(expected))
            self.assertLess(F(result['log2_upper_rational'])-c.base.rational_ceiling(expected), F(1, 10**9))
            self.assertEqual(result['marker_indices_by_parity_entry'], [[0, 1, 0, 1], [1, 1, 0, 1]])
            self.assertEqual(result['marker_entries_evaluated'], 16)
            self.assertLess(abs(F(balanced['log2_upper_rational'])-F(result['log2_upper_rational'])), F(1, 10**9))
            self.assertEqual(balanced['marker_indices_by_parity_entry'], result['marker_indices_by_parity_entry'])
            self.assertTrue(balanced['balanced_density'])
            self.assertFalse(balanced['outer_prefix_maximized_separately'])
            self.assertEqual(F(balanced['prefix_log2_upper_rational']), -100)
            fixed = []
            for even, odd in product(range(2), repeat=2):
                product_matrix = c.base.matrix_product(tuple(map(arb, matrices[0][even])),
                                                      tuple(map(arb, matrices[1][odd])))
                powered = c.base.matrix_power(product_matrix, 128)
                fixed.append(c.tail.log2arb(37*(powered[0]+powered[1]))+bank.constant+9*bank.slope)
            self.assertLess(c.base.rational(F(result['log2_upper_rational'])), min(value.lower() for value in fixed))

    def test_all_vertex_normalizers_bound_interior_and_both_faces(self):
        with ctx.workprec(256):
            bank = synthetic_bank()
            for cell in (c.Cell('interior', ((7, 10), (0, 4), (1, 4))),
                         c.Cell('even_face', ((4, 7), (1, 4))),
                         c.Cell('odd_face', ((4, 7), (1, 4)))):
                for balanced in (False, True):
                    result = c._entrywise_bound(bank, cell, 100, c.base.Budget(), balanced_density=balanced)
                    for point in product(*(range(lo, hi+1) for lo, hi in cell.ranges)):
                        if not cell.contains(point):
                            continue
                        weight, even, odd = cell.profile(point)
                        for parity, full in enumerate((even, odd)):
                            for entry in range(4):
                                index = result['marker_indices_by_parity_entry'][parity][entry]
                                value = bank.numerator(index, None)[entry]/normalizer(weight, full, bank.markers[index])
                                if balanced:
                                    value *= arb(2)**(bank.slope*weight/256)
                                ceiling = c.base.rational(F(result['regional_log2_upper_rational'][parity][entry]))
                                self.assertLessEqual(c.tail.log2arb(value).upper(), ceiling)
                        prefix = bank.constant if balanced else bank.constant+bank.slope*weight
                        self.assertLessEqual(prefix.upper(), c.base.rational(F(result['prefix_log2_upper_rational'])))
                    self.assertEqual(len(result['vertices']), len(cell.vertices()))
                # Each broad-cell mode is checked as its own complete bound;
                # selecting one mode never substitutes a pointwise comparison.

    def test_empty_zero_and_timeout_return_no_partial_claim(self):
        bank = synthetic_bank()
        empty = c.Cell('interior', ((1, 1), (1, 1), (0, 0)))
        result = c._entrywise_bound(bank, empty, 0, c.base.Budget())
        self.assertTrue(result['zero'])
        self.assertTrue(result['integer_empty'])
        self.assertIsNone(result['log2_upper_rational'])
        cell = c.Cell('interior', ((9, 9), (2, 2), (3, 3)))
        with self.assertRaises(ValueError):
            c._entrywise_bound(bank, cell, 0, c.base.Budget())
        with self.assertRaises(TimeoutError):
            c._entrywise_bound(bank, cell, 1, c.base.Budget(1e-12))
        bank.numerator = lambda index, rect: (arb(0),)*4
        self.assertTrue(c._entrywise_bound(bank, cell, 1, c.base.Budget())['zero'])

    def test_strict_provider_count_pins_and_public_receipt_metadata(self):
        fake, cell = synthetic_bank(), c.Cell('interior', ((9, 9), (2, 2), (3, 3)))
        with self.assertRaises(ValueError):
            c.entrywise_bound(fake, cell, 1, c.base.Budget())
        bank = object.__new__(c.provider.AliasBank)
        bank.__dict__.update(fake.__dict__, precision=256, pins={'mock-only': 'not-a-gate'})
        bank.check_pins = lambda: None
        bank.metadata = lambda: dict(construction='mock only', directed=True, whole_code_certificate=False)
        with patch.object(c.counts, 'support_count', return_value=10):
            with self.assertRaises(ValueError):
                c.entrywise_bound(bank, cell, 9, c.base.Budget())
            result = c.entrywise_bound(bank, cell, 10, c.base.Budget())
        self.assertTrue(result['support_count_verified'])
        self.assertTrue(result['complete_cell'])
        self.assertEqual(result['checked_projection_support_count'], 10)
        self.assertEqual(result['bank_source_sha256'], {'mock-only': 'not-a-gate'})
        self.assertFalse(result['whole_code_certificate'])
        pins = c.source_pins()
        for name in ('hd_field128_gamma_cell_entrywise.py', 'hd_field128_gamma_cell_counts.py',
                     'hd_field128_gamma_cell_geometry.py', 'hd_field128_gamma_all_pair_certificate.py'):
            self.assertIn(name, pins)
        with patch.object(c.counts, 'support_count', return_value=10), \
             patch.object(c.profiles, 'check_source_pins', side_effect=RuntimeError('changed source')):
            with self.assertRaises(RuntimeError):
                c.entrywise_bound(bank, cell, 10, c.base.Budget())

    def test_scope_and_nonfinite_normalizers_fail_closed(self):
        bank, budget = synthetic_bank(), c.base.Budget()
        tiny = c.Cell('interior', ((3, 3), (0, 0), (1, 1)), pairs=4)
        with self.assertRaises(ValueError):
            c._entrywise_bound(bank, tiny, 1, budget)
        cell = c.Cell('interior', ((9, 9), (2, 2), (3, 3)))
        with self.assertRaises(ValueError):
            c._entrywise_bound(bank, cell, 1, budget, balanced_density=1)
        bank.log_markers = ((arb('nan'), arb(1)),)*2
        with self.assertRaises(ArithmeticError):
            c._entrywise_bound(bank, cell, 1, budget)


if __name__ == '__main__':
    unittest.main()
