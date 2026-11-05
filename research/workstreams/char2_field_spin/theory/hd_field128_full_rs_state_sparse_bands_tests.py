"""Tiny q256 replay and mocked timeout checks; no full-field distance gate."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_full_rs_state_sparse_bands as bands


class StateSparseBandTests(unittest.TestCase):
    def test_shared_table_replays_each_tiny_selected_weight_and_union(self):
        z, alpha = F(7, 8), F(1, 4)
        rows = []
        result = bands.batch_receipt(1, 3, z, alpha, maximum=4, q=256, state=24,
                                     max_seconds=30, checkpoint=rows.append)
        self.assertEqual(result['regional_tables_built'], 1)
        self.assertEqual(result['coefficient_maximum'], 4)
        self.assertEqual(result['coefficient_marker'], '1/1024')
        self.assertEqual([row['active_groups'] for row in rows], [1, 2, 3])
        self.assertTrue(result['complete_requested_range'])
        self.assertEqual(result['pending_weight_ranges'], [])
        for row in rows:
            frozen = bands.certificate.weight_receipt(row['active_groups'], z, alpha, q=256,
                                                      state=24, distance=F(2, 5), max_seconds=30)
            self.assertLess(abs(F(row['log2_upper_rational'])-F(frozen['log2_upper_rational'])), F(1, 10**9))
            self.assertEqual(row['profiles'], frozen['profiles'])
            self.assertEqual(row['construction'], bands.local.construction(24))
            self.assertEqual(row['witness']['coefficient_maximum'], 4)
            self.assertTrue(row['directed'])
            self.assertFalse(row['whole_code_certificate'])
        with ctx.workprec(256):
            mass = sum((arb(2)**bands.base.rational(F(row['log2_upper_rational'])) for row in rows), arb(0))
            exact = bands.base.rational_ceiling(mass.log()/arb(2).log())
        self.assertEqual(F(result['log2_upper_rational']), exact)
        self.assertTrue(result['no_weight_interpolation'])

    def test_timeout_preserves_completed_row_and_explicit_pending_range(self):
        def row(weight, *args, **kwargs):
            if weight == 3:
                raise TimeoutError('mock row deadline')
            return dict(active_groups=weight, profiles=weight, log2_upper_rational='-80',
                        log2_upper_enclosure='-80', covers_all_supports_at_weight=True,
                        covers_all_even_odd_splits=True)
        saved = []
        with patch.object(bands.certificate, 'regional_directed', return_value={}) as builder, \
             patch.object(bands, '_weight_from_regions', side_effect=row):
            result = bands.batch_receipt(2, 4, F(7, 8), F(1, 4), maximum=64, checkpoint=saved.append)
        builder.assert_called_once()
        self.assertEqual([item['active_groups'] for item in saved], [2])
        self.assertEqual(result['completed_weight_range'], [2, 2])
        self.assertEqual(result['pending_weight_ranges'], [[3, 4]])
        self.assertFalse(result['every_weight_evaluated'])
        self.assertFalse(result['covers_all_supports_in_range'])
        self.assertIsNone(result['log2_upper_rational'])
        self.assertEqual(result['failure']['active_groups'], 3)
        self.assertEqual(result['failure']['phase'], 'weight')
        self.assertLessEqual(F(result['completed_union_log2_upper_rational']), -79)

    def test_timeout_before_table_reports_every_weight_pending(self):
        with patch.object(bands.certificate, 'regional_directed') as builder:
            result = bands.batch_receipt(1, 4, F(7, 8), F(1, 4), max_seconds=1e-12)
        builder.assert_not_called()
        self.assertEqual(result['regional_tables_built'], 0)
        self.assertEqual(result['weight_receipts'], [])
        self.assertEqual(result['pending_weight_ranges'], [[1, 4]])
        self.assertEqual(result['failure']['phase'], 'regional-table')

    def test_validation_and_exclusive_hash_pinned_row_files(self):
        for first, last, maximum in ((0, 2, 2), (3, 2, 4), (1, 5, 4), (1, 65, 65)):
            with self.assertRaises(ValueError):
                bands.batch_receipt(first, last, F(7, 8), F(1, 4), maximum=maximum)
        with self.assertRaises(ValueError):
            bands.batch_receipt(1, 1, .875, F(1, 4))
        with TemporaryDirectory() as directory:
            path = Path(directory)/'weight_0002.json'
            receipt = bands.write_new(path, {'active_groups': 2, 'scope': 'mock'})
            self.assertEqual(receipt['sha256'], sha256(path.read_bytes()).hexdigest())
            self.assertEqual(json.loads(path.read_bytes())['active_groups'], 2)
            with self.assertRaises(FileExistsError):
                bands.write_new(path, {})
        for name in ('hd_field128_full_rs_state_sparse_bands.py', 'hd_field128_full_rs_state_certificate.py',
                     'hd_field128_full_rs_state_screen.py', 'hd_field128_lift_certificate.py'):
            self.assertIn(name, bands.source_pins())


if __name__ == '__main__':
    unittest.main()
