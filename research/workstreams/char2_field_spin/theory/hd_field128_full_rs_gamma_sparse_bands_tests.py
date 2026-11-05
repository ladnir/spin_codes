"""Synthetic-region arithmetic and mocked fail-closed checks; no distance gates."""

from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_full_rs_gamma_sparse_bands as bands


def synthetic_regions(maximum):
    return {(r, full): tuple(bands.base.rational(F((entry+1)*(full+1), 257+r+entry)) for entry in range(4))
            for r in range(maximum+1) for full in range(r//2+1)}


def fake_row(weight, *args, **kwargs):
    return dict(active_groups=weight, profiles=weight, log2_upper_rational='-80',
                log2_upper_enclosure='-80', covers_all_supports_at_weight=True,
                covers_all_even_odd_splits=True)


class GammaSparseBandTests(unittest.TestCase):
    def test_shared_table_replays_frozen_per_weight_arithmetic_and_union(self):
        z, alpha, saved = F(7, 8), F(1, 4), []
        with ctx.workprec(256):
            regions = synthetic_regions(4)
        with patch.object(bands.certificate, 'regional_directed', return_value=regions) as builder:
            result = bands.batch_receipt(1, 3, z, alpha, maximum=4, state=30, checkpoint=saved.append)
        builder.assert_called_once()
        self.assertEqual(builder.call_args.args, (z, alpha, 4))
        self.assertEqual(builder.call_args.kwargs['state'], 30)
        self.assertEqual(result['regional_tables_built'], 1)
        self.assertEqual([row['active_groups'] for row in saved], [1, 2, 3])
        self.assertEqual(result['coefficient_marker'], '1/1024')
        self.assertTrue(result['every_integer_weight_evaluated'])
        self.assertTrue(result['complete_requested_range'])
        self.assertTrue(result['no_weight_interpolation'])
        for row in saved:
            with patch.object(bands.certificate, 'regional_directed', return_value=regions):
                control = bands.certificate.weight_receipt(row['active_groups'], z, alpha,
                                                           state=30, distance=F(2, 5), max_seconds=10)
            self.assertEqual(row['log2_upper_rational'], control['log2_upper_rational'])
            self.assertEqual(row['profiles'], control['profiles'])
            self.assertEqual(row['construction'], bands.local.construction(30))
            self.assertEqual(row['feedback'], 'gamma A transpose')
            self.assertEqual(row['witness']['coefficient_maximum'], 4)
            self.assertFalse(row['all_pair_before_alpha'])
            self.assertFalse(row['whole_code_certificate'])
        with ctx.workprec(256):
            mass = sum((arb(2)**bands.base.rational(F(row['log2_upper_rational'])) for row in saved), arb(0))
            expected = bands.base.rational_ceiling(mass.log()/arb(2).log())
        self.assertEqual(F(result['log2_upper_rational']), expected)

    def test_every_public_state_metadata_and_exact_provider_dispatch(self):
        for state in bands.certificate.STATES:
            with patch.object(bands.certificate, 'regional_directed', return_value={}), \
                 patch.object(bands, '_weight_from_regions', side_effect=fake_row) as evaluate:
                result = bands.batch_receipt(1, 1, F(1, 2), F(1, 4), state=state)
            self.assertEqual(result['construction'], bands.local.construction(state))
            self.assertEqual(result['state_dimension'], state)
            self.assertEqual(result['gamma_extension_degree'], 16)
            self.assertEqual(result['outer_measure'], 'invertible-toeplitz')
            self.assertEqual(evaluate.call_args.kwargs['q'], 1 << 128)
            self.assertFalse(result['old_numerical_bounds_reused'])
        self.assertIs(bands._weight_from_regions, bands.frozen_driver._weight_from_regions)

    def test_timeout_preserves_immediate_rows_and_pending_range_without_band_bound(self):
        saved = []
        def stop(weight, *args, **kwargs):
            if weight == 3:
                self.assertEqual([row['active_groups'] for row in saved], [2])
                raise TimeoutError('mock row deadline')
            return fake_row(weight)
        with patch.object(bands.certificate, 'regional_directed', return_value={}) as builder, \
             patch.object(bands, '_weight_from_regions', side_effect=stop):
            result = bands.batch_receipt(2, 4, F(7, 8), F(1, 4), maximum=64, checkpoint=saved.append)
        builder.assert_called_once()
        self.assertEqual(result['completed_weight_range'], [2, 2])
        self.assertEqual(result['pending_weight_ranges'], [[3, 4]])
        self.assertEqual(result['status'], 'incomplete')
        self.assertFalse(result['every_integer_weight_evaluated'])
        self.assertFalse(result['covers_all_supports_in_range'])
        self.assertFalse(result['covers_all_even_odd_splits'])
        self.assertIsNone(result['log2_upper_rational'])
        self.assertEqual(result['failure']['phase'], 'weight')
        self.assertEqual(result['failure']['active_groups'], 3)
        self.assertLessEqual(F(result['completed_union_log2_upper_rational']), -79)

    def test_timeout_before_build_and_arithmetic_failure_are_closed(self):
        with patch.object(bands.certificate, 'regional_directed') as builder:
            result = bands.batch_receipt(1, 4, F(7, 8), F(1, 4), max_seconds=1e-12)
        builder.assert_not_called()
        self.assertEqual(result['regional_tables_built'], 0)
        self.assertEqual(result['pending_weight_ranges'], [[1, 4]])
        self.assertEqual(result['weight_receipts'], [])
        self.assertIsNone(result['log2_upper_rational'])
        with patch.object(bands.certificate, 'regional_directed', side_effect=ArithmeticError('mock invalid bound')):
            result = bands.batch_receipt(1, 2, F(1, 2), F(1, 4))
        self.assertEqual(result['failure']['error_type'], 'ArithmeticError')
        self.assertIsNone(result['completed_union_log2_upper_rational'])

    def test_scope_exactness_exclusive_files_and_transitive_pin_guard(self):
        for first, last, maximum in ((0, 2, 2), (3, 2, 4), (1, 5, 4), (1, 65, 65)):
            with self.assertRaises(ValueError):
                bands.batch_receipt(first, last, F(7, 8), F(1, 4), maximum=maximum)
        for kwargs in (dict(state=32), dict(bits=8), dict(precision=64), dict(paired=1), dict(failure_bits=-1)):
            with self.assertRaises(ValueError):
                bands.batch_receipt(1, 1, F(7, 8), F(1, 4), **kwargs)
        with self.assertRaises(ValueError):
            bands.batch_receipt(1, 1, .875, F(1, 4))
        with TemporaryDirectory() as directory:
            path = Path(directory)/'weight_0002.json'
            receipt = bands.write_new(path, {'active_groups': 2, 'scope': 'mock only'})
            self.assertEqual(receipt['sha256'], sha256(path.read_bytes()).hexdigest())
            with self.assertRaises(FileExistsError):
                bands.write_new(path, {})
        pins = bands.source_pins()
        for name in ('hd_field128_full_rs_gamma_sparse_bands.py', 'hd_field128_full_rs_gamma_certificate.py',
                     'hd_field128_full_rs_gamma_screen.py', 'hd_field128_full_rs_state_sparse_bands.py',
                     'hd_shifted_pairs_count.py'):
            self.assertIn(name, pins)
        bands.check_pins(pins)
        changed = dict(pins)
        changed['hd_field128_full_rs_gamma_sparse_bands.py'] = 'changed'
        with self.assertRaises(RuntimeError):
            bands.check_pins(changed)

    def test_incomplete_cli_checkpoints_summary_and_exits_nonzero(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)/'new-output'
            argv = ['gamma-sparse-bands', '--first', '1', '--last', '2', '--z', '1/2',
                    '--alpha', '1/4', '--output-dir', str(output)]
            def stop(weight, *args, **kwargs):
                if weight == 2:
                    raise TimeoutError('mock stop')
                return fake_row(weight)
            with patch('sys.argv', argv), patch.object(bands.certificate, 'regional_directed', return_value={}), \
                 patch.object(bands, '_weight_from_regions', side_effect=stop), redirect_stdout(StringIO()):
                with self.assertRaises(SystemExit) as error:
                    bands.main()
            self.assertEqual(error.exception.code, 2)
            row = json.loads((output/'weight_0001.json').read_bytes())
            result = json.loads((output/'summary.json').read_bytes())
            self.assertEqual(row['active_groups'], 1)
            self.assertEqual(row['source_sha256'], result['source_sha256'])
            self.assertIsNone(result['log2_upper_rational'])
            self.assertEqual(result['pending_weight_ranges'], [[2, 2]])


if __name__ == '__main__':
    unittest.main()
