"""Tiny table/oracle and mocked-provider tests; no actual bank or distance gate."""

from fractions import Fraction as F
from itertools import product
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_cell_regional_table as c
import hd_field128_gamma_cell_cover as cover
from hd_shifted_pairs_count import exact_profiles, profile_count_upper


def regions(pairs, weight):
    return {f: tuple(c.base.rational(F((f+2)*(entry+1), 101+7*f+entry)) for entry in range(4))
            for f in c.tail.full_counts(pairs, weight)}


def mock_bank():
    bank = object.__new__(c.provider.AliasBank)
    bank.precision, bank.constant, bank.slope = 256, arb(-1000), arb(2)
    bank.construction, bank.state, bank.q = 'synthetic-test-only', 30, 1 << 128
    bank.distance, bank.forbidden = F(2, 5), 838860
    bank.z, bank.alpha, bank.grid, bank.bits = F(1, 2), F(1, 8), 128, 80
    bank.markers = ((F(1, 2), F(1, 3)), (F(2, 3), F(3, 5)))
    bank.pins, bank.check_pins = {'mock': 'not-a-gate'}, lambda: None
    return bank


class RegionalTableTests(unittest.TestCase):
    def test_every_sparse_interval_matches_direct_componentwise_max(self):
        with ctx.workprec(256):
            for pairs in range(1, 8):
                for weight in range(1, 2*pairs):
                    rows = regions(pairs, weight)
                    table = c.RegionalMaxTable(pairs, weight, rows, c.base.Budget())
                    for low in rows:
                        for high in range(low, table.last+1):
                            expected = tuple(max(rows[f][entry].upper() for f in range(low, high+1))
                                             for entry in range(4))
                            self.assertEqual(table.query(low, high), expected)
                    self.assertEqual(len(table.levels), len(rows).bit_length())

    def test_exact_locked_projections_interior_faces_empty_and_rejection(self):
        for pairs in range(2, 7):
            for weight in range(1, 2*pairs):
                for root in cover.selected_weight_roots(weight, pairs):
                    intervals = c.full_intervals(root, weight)
                    points = [root.profile(point) for point in product(
                        *(range(low, high+1) for low, high in root.ranges)) if root.contains(point)]
                    if not points:
                        self.assertIsNone(intervals)
                    else:
                        expected = tuple((min(point[parity+1] for point in points),
                                          max(point[parity+1] for point in points)) for parity in (0, 1))
                        self.assertEqual(intervals, expected)
        for cell in (c.Cell('interior', ((5, 6), (0, 2), (0, 2))),
                     c.Cell('even_face', ((2, 3), (1, 2))),
                     c.Cell('odd_face', ((3, 3), (1, 2)))):
            with self.assertRaises(ValueError):
                c.full_intervals(cell, 5)
        self.assertIsNone(c.full_intervals(c.Cell('interior', ((3, 3), (2, 2), (0, 0))), 3))

    def test_singletons_match_frozen_profile_arithmetic_and_broad_cells_dominate(self):
        with ctx.workprec(256):
            pairs, weight, prefix = 5, 4, arb(-100)
            rows = regions(pairs, weight)
            # These rows are not commuting; order is checked by the direct oracle.
            rows[0], rows[2] = tuple(map(arb, (1, 2, 4, 1))), tuple(map(arb, (2, 4, 1, 2)))
            table = c.RegionalMaxTable(pairs, weight, rows, c.base.Budget())
            actual_counts = exact_profiles(pairs)
            for root in cover.selected_weight_roots(weight, pairs):
                count = c.counts.support_count(root, c.base.Budget())
                bound = c._cell_bound(table, root, prefix, count, c.base.Budget())
                actual_sum = arb(0)
                for point in product(*(range(lo, hi+1) for lo, hi in root.ranges)):
                    if not root.contains(point):
                        continue
                    r, a, b = root.profile(point)
                    support = profile_count_upper(pairs, r, a, b)
                    cell = c.Cell(root.kind, tuple((x, x) for x in point), pairs)
                    row = c._cell_bound(table, cell, prefix, support, c.base.Budget())
                    direct = c.base.matrix_power(c.base.matrix_product(rows[a], rows[b]), 128)
                    actual_sum += actual_counts.get((r, a, b), 0)*(direct[0]+direct[1])
                    value = c.tail.log2arb(support*(direct[0]+direct[1]))+prefix
                    self.assertGreaterEqual(F(row['log2_upper_rational']), c.base.rational_ceiling(value))
                    self.assertLess(F(row['log2_upper_rational'])-c.base.rational_ceiling(value), F(1, 10**9))
                    self.assertGreaterEqual(F(bound['log2_upper_rational']), F(row['log2_upper_rational']))
                actual_log = c.tail.log2arb(actual_sum)+prefix
                self.assertGreaterEqual(F(bound['log2_upper_rational']), c.base.rational_ceiling(actual_log))

    def test_constructor_calls_frozen_normalization_once_and_driver_recognizes_subclass(self):
        bank = mock_bank()
        with ctx.workprec(256):
            rows = {f: (c.base.rational(F(1, 4)),)*4 for f in range(3)}
        with patch.object(c.restored, 'RegionalCache') as cache:
            cache.return_value.regional.return_value = rows
            evaluator = c.Evaluator(bank, 4)
            for root in cover.selected_weight_roots(4):
                result = evaluator.receipt(root, failure_bits=64)
                self.assertTrue(result['complete_cell'])
                self.assertTrue(cover._receipt_pass(result, root, 64))
                self.assertEqual(result['regional_table_builds'], 1)
                self.assertEqual(result['regional_table_source_sha256'], evaluator.pins)
                self.assertEqual(result['witness']['markers'], [['1/2', '1/3'], ['2/3', '3/5']])
            cache.assert_called_once_with(bank)
            self.assertEqual(cache.return_value.regional.call_count, 1)
            self.assertEqual(cache.return_value.regional.call_args.args[0], 4)
            result = cover.cover(evaluator, 4, max_nodes=10, max_seconds=30., cell_seconds=5.)
            self.assertTrue(result['complete_selected_weight'])
            self.assertTrue(result['directed'])
        self.assertIsInstance(evaluator, c.numeric.Evaluator)
        for cell in (c.Cell('interior', ((3, 3), (0, 1), (0, 1))),
                     c.Cell('even_face', ((2, 3), (1, 1)))):
            with self.assertRaises(ValueError):
                evaluator.receipt(cell)
        with self.assertRaises(ValueError):
            evaluator.receipt(cover.selected_weight_roots(4)[0], marker_pairs=((0, 0),))

    def test_failed_build_and_cell_are_closed_and_complete_keys_are_required(self):
        bank = mock_bank()
        with patch.object(c.restored, 'RegionalCache') as cache:
            cache.return_value.regional.side_effect = TimeoutError('mock build deadline')
            with self.assertRaises(TimeoutError):
                c.Evaluator(bank, 4)
        with ctx.workprec(256):
            rows = regions(4096, 4)
        with patch.object(c.restored, 'RegionalCache') as cache:
            cache.return_value.regional.return_value = rows
            evaluator = c.Evaluator(bank, 4)
        result = evaluator.receipt(cover.selected_weight_roots(4)[0], max_seconds=1e-12)
        self.assertFalse(result['complete_cell'])
        self.assertFalse(result['directed'])
        self.assertIsNone(result['log2_upper_rational'])
        self.assertFalse(result['meets_cell_budget'])
        with self.assertRaises(ValueError):
            c.RegionalMaxTable(4096, 4, {0: rows[0]}, c.base.Budget())
        with self.assertRaises(ValueError):
            evaluator.table.query(-1, 0)
        with self.assertRaises(ValueError):
            c.Evaluator(SimpleNamespace(), 4)
        for weight in (0, 8192, True):
            with self.assertRaises(ValueError):
                c.Evaluator(bank, weight)
        pins = c.source_pins()
        for name in ('hd_field128_gamma_cell_regional_table.py', 'hd_field128_lift_alias_cover.py',
                     'hd_field128_gamma_cell_counts.py', 'hd_field128_gamma_cell_certificate.py'):
            self.assertIn(name, pins)
        with patch.object(c.profiles, 'check_source_pins', side_effect=RuntimeError('changed helper')):
            with self.assertRaises(RuntimeError):
                evaluator.check_pins()


if __name__ == '__main__':
    unittest.main()
