"""Tiny synthetic contracts; no production index or coefficient bank is evaluated."""

from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v16_tests as fixtures
import hd_field128_gamma_selected_index_v17 as mode


def parent_fixture(folder, *, source_conflict=False):
    event = fixtures.write(folder/'prior.json', dict(synthetic_fixture=True))
    audit = fixtures.write(folder/'prior_audit.json', dict(synthetic_prior_audit=True))
    entries = [dict(event, weight_range=[r, r], log2_upper_rational='-100') for r in range(1, 114)]
    entries[-1].update(zero=True, log2_upper_rational=None, audit=audit)
    parent = dict(schema='field128-gamma-partial-cover-index-v16',
        construction=mode.midpoint.previous.EXPECTED['construction'], field_order=1 << 128,
        source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        whole_code_certificate=False, directed_aggregation=True, entries=entries,
        covered_weight_ranges=[[1, 113]], covered_weight_count=113,
        pending_weight_ranges=[[114, 8192]], pending_weight_count=8079, receipt_count=113,
        covered_events_union_log2_upper_rational='-90',
        producer_receipt_counts={'hd_field128_gamma_midpoint_weight_bands': 8},
        latest_audited_multiweight_events=[entries[-1]],
        aggregation_source_sha256=({'hd_field128_gamma_midpoint_weight_bands.py': 'wrong'}
                                   if source_conflict else {}))
    pin = fixtures.write(folder/'parent.json', parent)
    return Path(pin['path']), pin['sha256']


class IndexTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.folder = Path(self.temporary.name)
        self.parent, self.parent_hash = parent_fixture(self.folder/'old')

    def tearDown(self):
        self.temporary.cleanup()

    def extend(self, paths):
        with patch.object(mode, 'PARENT_SHA256', self.parent_hash):
            return mode.extend_v16(self.parent, paths)

    def test_retains_all_113_entries_and_reaggregates_two_new_bands(self):
        first, unused = fixtures.band_fixture(self.folder/'first', 114, 115)
        second, unused = fixtures.band_fixture(self.folder/'second', 117, 118)
        raw = self.parent.read_bytes()
        prior = json.loads(raw)
        result = self.extend([first, second])
        self.assertEqual(result['schema'], 'field128-gamma-partial-cover-index-v17')
        self.assertEqual(result['entries'][:113], prior['entries'])
        self.assertEqual(self.parent.read_bytes(), raw)
        self.assertEqual(result['receipt_count'], 115)
        self.assertEqual(result['added_event_count'], 2)
        self.assertEqual(result['added_weight_count'], 4)
        self.assertEqual(result['covered_weight_ranges'], [[1, 115], [117, 118]])
        self.assertEqual(result['pending_weight_ranges'], [[116, 116], [119, 8192]])
        expected = mode.index.union_ceiling(mode.frozen._bounds(result['entries']), 384)
        self.assertEqual(result['covered_events_union_log2_upper_rational'], str(expected))
        self.assertEqual(result['producer_receipt_counts']['hd_field128_gamma_midpoint_weight_bands'], 10)
        self.assertEqual(result['parent_index']['sha256'], self.parent_hash)
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['coefficient_replay'])
        self.assertIn('hd_field128_gamma_selected_index_v17.py', result['aggregation_source_sha256'])
        self.assertIn('hd_field128_gamma_selected_index_v16.py', result['aggregation_source_sha256'])

    def test_zero_added_event_preserves_null_semantics(self):
        path, unused = fixtures.band_fixture(self.folder/'zero', 114, 115, zero=True)
        result = self.extend([path])
        self.assertTrue(result['added_events_zero'])
        self.assertIsNone(result['added_events_union_log2_upper_rational'])
        self.assertEqual(result['covered_weight_count'], 115)
        self.assertTrue(result['entries'][-1]['zero'])
        self.assertIsNone(result['entries'][-1]['log2_upper_rational'])
        self.assertEqual(result['covered_events_union_log2_upper_rational'],
                         str(mode.index.union_ceiling([F(-100)]*112, 384)))

    def test_explicit_limit_accepts_64_and_rejects_invalid_lists(self):
        # Isolate the adapter's 64-entry limit; real band validation is covered
        # above and inherited unchanged from the reviewed v16 validator.
        paths = [self.folder/f'audit{i}.json' for i in range(64)]
        for i, path in enumerate(paths):
            fixtures.write(path, dict(synthetic_fixture=True, i=i))
        def event(path, loader, pins):
            unused, pin = loader.load(path, role='synthetic-limit-test')
            r = 114+paths.index(path)
            return dict(pin, weight_range=[r, r], weights_evaluated=1,
                        zero=False, log2_upper_rational='-200')
        with patch.object(mode.frozen, '_band_event', side_effect=event):
            result = self.extend(paths)
        self.assertEqual(result['added_event_count'], 64)
        self.assertEqual(result['receipt_count'], 177)
        self.assertEqual(result['covered_weight_ranges'], [[1, 177]])
        for invalid in ([], paths+[paths[0]], str(paths[0]), iter(paths)):
            with self.assertRaisesRegex(ValueError, 'explicit audited bands'):
                self.extend(invalid)

    def test_overlap_duplicate_and_hard_parent_pin_rejected(self):
        old, unused = fixtures.band_fixture(self.folder/'overlap', 113, 114)
        new, unused = fixtures.band_fixture(self.folder/'new', 114, 115)
        for paths in ([old], [new, new]):
            with self.assertRaisesRegex(ValueError, 'overlaps'):
                self.extend(paths)
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            mode.extend_v16(self.parent, [new])

    def test_inherited_source_and_new_receipt_corruption_rejected(self):
        valid, pin = fixtures.band_fixture(self.folder/'valid', 114, 115)
        self.parent, self.parent_hash = parent_fixture(self.folder/'conflict', source_conflict=True)
        with self.assertRaisesRegex(ValueError, 'Conflicting source pin'):
            self.extend([valid])
        self.parent, self.parent_hash = parent_fixture(self.folder/'clean')
        bad, unused = fixtures.band_fixture(self.folder/'partial', 114, 115,
                                            mutate=lambda row: row.update(complete_requested_band=False))
        with self.assertRaises(ValueError):
            self.extend([bad])
        read = Path.read_bytes
        target = Path(pin['path'])
        def changed(path):
            raw = read(path)
            return raw+b' ' if path == target else raw
        with patch.object(Path, 'read_bytes', changed):
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                self.extend([valid])


if __name__ == '__main__':
    unittest.main()
