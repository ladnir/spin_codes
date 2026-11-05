"""Synthetic compact artifacts only; no coefficient or FFT construction."""

from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import ctx

import hd_field128_gamma_near_face_index_update as update
from hd_field128_gamma_near_face_band_audit_tests import synthetic_bank
from hd_field128_gamma_near_face_band_tests import FakeCache


def audited_fixture(root):
    bank = synthetic_bank()
    metadata = bank.metadata()
    metadata['witness']['alpha'] = '1/'+str(1 << 128)
    bank.metadata = lambda: deepcopy(metadata)
    with ctx.workprec(256):
        bank.constant = update.audited.outer_prefix(metadata, 0)
        bank.slope = update.audited.outer_prefix(metadata, 1)-bank.constant
    FakeCache.fail_weight = None
    with patch.object(update.band, 'RegionalCache', FakeCache):
        result = update.band.band_receipt(bank, 1, 3, output_dir=root/'band', max_weights=3, max_singletons=2)
    report = update.audited.audit_band(root/'band'/'summary.json')
    pin = update.band.write_new(root/'audit.json', report)
    return result, report, pin


def prior_index(root, weight=9):
    receipt = update.band.write_new(root/'prior-receipt.json', {'test_only': True})
    event = dict(receipt, weight_range=[weight, weight], log2_upper_rational='-100')
    value = dict(schema='field128-gamma-partial-cover-index-v9',
                 **{key: update.audited.EXPECTED[key] for key in ('construction', 'field_order', 'source_dimension',
                     'output_length', 'target_distance', 'whole_code_certificate')},
                 entries=[event], receipt_count=1, covered_weight_ranges=[[weight, weight]], covered_weight_count=1,
                 pending_weight_ranges=update.index.ranges(set(range(1, 8193))-{weight}), pending_weight_count=8191,
                 covered_events_union_log2_upper_rational='-99', aggregation_source_sha256={}, producer_receipt_counts={})
    return update.band.write_new(root/'prior.json', value)


class UpdateTests(unittest.TestCase):
    def test_slice_preserves_context_and_rows_with_fresh_subset_union(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            result, report, pin = audited_fixture(root)
            value = update.slice_receipt(pin['path'], 2, 3)
            self.assertEqual(value['context'], result['context'])
            self.assertEqual(value['weight_receipts'], result['weight_receipts'][1:])
            self.assertEqual(value['excluded_parent_weights'], [[1, 1]])
            self.assertFalse(value['parent_event_included'])
            self.assertFalse(value['whole_code_certificate'])
            self.assertEqual(value['requested_failure_bits'], 64)
            expected = update.index.union_ceiling([update.fraction(x['log2_upper_rational'])
                                                   for x in result['weight_receipts'][1:]], 384)
            self.assertEqual(value['log2_upper_rational'], str(expected))
            for lo, hi in ((0, 2), (2, 4), (3, 2)):
                with self.assertRaises(ValueError):
                    update.slice_receipt(pin['path'], lo, hi)

    def test_v10_accepts_loader_byte_count_and_rejects_overlap(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            result, report, pin = audited_fixture(root)
            self.assertIn('bytes', report['context'])
            self.assertNotIn('bytes', result['context'])
            prior = prior_index(root)
            events = [('compact-band', root/'band'/'summary.json', pin['path'])]
            value = update.extend_index(prior['path'], events, version=10)
            self.assertEqual(value['covered_weight_ranges'], [[1, 3], [9, 9]])
            self.assertEqual(value['added_weight_count'], 3)
            self.assertEqual(value['receipt_count'], 2)
            self.assertFalse(value['whole_code_certificate'])
            newer = update.band.write_new(root/'v10.json', value)
            with self.assertRaisesRegex(ValueError, 'overlaps'):
                update.extend_index(newer['path'], events, version=11)

    def test_derived_slice_is_recomputed_not_trusted(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            result, report, pin = audited_fixture(root)
            value = update.slice_receipt(pin['path'], 2, 3)
            saved = update.band.write_new(root/'slice.json', value)
            loader = update.index.Validator(update.F(2, 5))
            event, closure = update._event('saved-slice', saved['path'], None, loader)
            self.assertEqual(event['weight_range'], [2, 3])
            value['log2_upper_rational'] = '-999999'
            bad = update.band.write_new(root/'bad-slice.json', value)
            with self.assertRaisesRegex(ValueError, 'Derived slice differs'):
                update._event('saved-slice', bad['path'], None, loader)

    def test_changed_row_bytes_and_audit_scope_are_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            result, report, pin = audited_fixture(root)
            bad = deepcopy(report)
            bad['target_distance'] = '9/20'
            bad_pin = update.band.write_new(root/'bad-audit.json', bad)
            with self.assertRaises(ValueError):
                update.slice_receipt(bad_pin['path'], 2, 3)
            Path(result['weight_receipts'][0]['path']).write_bytes(b'test mutation')
            with self.assertRaisesRegex(RuntimeError, 'artifact changed'):
                update.slice_receipt(pin['path'], 2, 3)


if __name__ == '__main__':
    unittest.main()
