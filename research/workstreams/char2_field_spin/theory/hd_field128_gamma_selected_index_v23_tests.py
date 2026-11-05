"""Synthetic v23 extension tests; no production index, bank, or remote job."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v16_tests as fixtures
import hd_field128_gamma_selected_index_v22_typed_tests as typed
import hd_field128_gamma_selected_index_v23 as mode


def parent_fixture(folder, mutate=None):
    event = fixtures.write(folder/'prior.json', dict(synthetic_fixture=True))
    audit = fixtures.write(folder/'prior_audit.json', dict(synthetic_prior_audit=True))
    entries = [dict(event, weight_range=[r, r], log2_upper_rational='-100') for r in range(1, 241)]
    entries.extend([dict(event, weight_range=[241, 7598], log2_upper_rational='-100'),
                    dict(event, weight_range=[8192, 8192], zero=True,
                         log2_upper_rational=None, audit=audit)])
    parent = dict(schema='field128-gamma-partial-cover-index-v22-typed',
        construction=mode.midpoint.previous.EXPECTED['construction'], field_order=1 << 128,
        source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        whole_code_certificate=False, directed_aggregation=True, entries=entries,
        covered_weight_ranges=[[1, 7598], [8192, 8192]], covered_weight_count=7599,
        pending_weight_ranges=[[7599, 8191]], pending_weight_count=593, receipt_count=242,
        covered_events_union_log2_upper_rational='-90', aggregation_source_sha256={},
        producer_receipt_counts={'hd_field128_gamma_midpoint_weight_bands': 37})
    if mutate:
        mutate(parent)
    pin = fixtures.write(folder/'parent.json', parent)
    return Path(pin['path']), pin['sha256']


class IndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pins = mode.source_pins()

    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.folder = Path(self.temporary.name)
        self.parent, self.parent_hash = parent_fixture(self.folder/'old')
        self.pins_patch = patch.object(mode, 'source_pins', side_effect=lambda: dict(self.pins))
        self.pins_patch.start()

    def tearDown(self):
        self.pins_patch.stop()
        self.temporary.cleanup()

    def extend(self, paths):
        with patch.object(mode, 'PARENT_SHA256', self.parent_hash):
            return mode.extend_v22(self.parent, paths)

    def fake_band(self, name, first=7599, last=None, *, bound='-200', zero=False):
        pin = fixtures.write(self.folder/name, dict(weight_range=[first, first if last is None else last],
            zero=zero, log2_upper_rational=None if zero else bound))
        return Path(pin['path'])

    @staticmethod
    def fake_event(path, loader, pins, source_inputs):
        row, pin = loader.load(path, role='synthetic-adapter-mechanics-only')
        lo, hi = row['weight_range']
        return dict(pin, **row, weights_evaluated=hi-lo+1)

    def test_real_typed_v1_v2_validation_and_unchanged_parent_prefix(self):
        first, unused = typed.typed_fixture(self.folder/'v1', 7599, 7600)
        second, unused = typed.typed_fixture(self.folder/'v2', 7602, 7603, zero=True, version=2)
        raw = self.parent.read_bytes()
        prior = json.loads(raw)
        result = self.extend([first, second])
        self.assertEqual(result['schema'], 'field128-gamma-partial-cover-index-v23')
        self.assertEqual(result['entries'][:242], prior['entries'])
        self.assertEqual(json.dumps(result['entries'][:242]), json.dumps(prior['entries']))
        self.assertEqual(self.parent.read_bytes(), raw)
        self.assertEqual(result['receipt_count'], 244)
        self.assertEqual(result['covered_weight_count'], 7603)
        self.assertEqual(result['pending_weight_count'], 589)
        self.assertEqual(result['covered_weight_ranges'], [[1, 7600], [7602, 7603], [8192, 8192]])
        self.assertEqual(result['added_weight_count'], 4)
        self.assertEqual(result['producer_receipt_counts']['hd_field128_gamma_midpoint_weight_bands'], 39)
        self.assertEqual(result['covered_events_union_log2_upper_rational'],
            str(mode.index.union_ceiling(mode.frozen._bounds(result['entries']), 384)))
        self.assertEqual(result['covered_events_failure_bits'], 54)
        self.assertTrue(result['prior_semantic_audits_inherited'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['coefficient_replay'])
        self.assertEqual(result['parent_index']['sha256'], self.parent_hash)
        self.assertTrue(result['entries'][-1]['zero'])
        self.assertIn('hd_field128_gamma_selected_index_v23.py', result['aggregation_source_sha256'])
        self.assertIn('hd_field128_gamma_selected_index_v22_typed.py', result['aggregation_source_sha256'])

    def test_input_count_limit_and_full_coverage_stays_partial(self):
        paths = [self.fake_band(f'event{i}.json', 7599+i) for i in range(64)]
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            result = self.extend(paths)
        self.assertEqual(result['added_event_count'], 64)
        for invalid in ([], paths+[paths[0]], str(paths[0]), iter(paths)):
            with self.assertRaisesRegex(ValueError, 'explicit audited bands'):
                self.extend(invalid)
        all_missing = [self.fake_band(f'fill{lo}.json', lo, min(8191, lo+63))
                       for lo in range(7599, 8192, 64)]
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            result = self.extend(all_missing)
        self.assertEqual(result['covered_weight_count'], 8192)
        self.assertEqual(result['pending_weight_ranges'], [])
        self.assertFalse(result['whole_code_certificate'])
        self.assertIn('no whole-code certificate', result['scope'])

    def test_parent_hash_and_parent_full_claim_rejected(self):
        new = self.fake_band('new.json')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            mode.extend_v22(self.parent, [new])
        for i, mutate in enumerate((lambda p: p.update(whole_code_certificate=True),
                lambda p: p.update(receipt_count=241), lambda p: p['entries'].pop(),
                lambda p: p.update(covered_weight_count=7600),
                lambda p: p.update(covered_events_union_log2_upper_rational='-100'),
                lambda p: p['aggregation_source_sha256'].update(
                    {'hd_field128_gamma_midpoint_weight_bands.py': 'wrong'}))):
            self.parent, self.parent_hash = parent_fixture(self.folder/f'bad{i}', mutate)
            with self.assertRaises(ValueError):
                self.extend([new])

    def test_overlap_duplicate_and_aggregate_budget_rejected(self):
        old = self.fake_band('old.json', 7598, 7599)
        new = self.fake_band('new.json', 7599, 7600)
        threshold = self.fake_band('threshold.json', bound='-54')
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            for paths in ([old], [new, new]):
                with self.assertRaisesRegex(ValueError, 'overlaps'):
                    self.extend(paths)
            with self.assertRaisesRegex(ValueError, 'union fails 54-bit'):
                self.extend([threshold])

    def test_only_explicit_typed_audits_and_no_new_full_claim(self):
        legacy, unused = fixtures.band_fixture(self.folder/'legacy', 7599, 7600)
        with self.assertRaisesRegex(ValueError, 'typed v1/v2'):
            self.extend([legacy])
        for i, mutate in enumerate((lambda r: r.update(whole_code_certificate=True),
                lambda r: r.update(complete_requested_band=False),
                lambda r: r.update(typed_provenance_checked=False),
                lambda r: r.update(nested_source_maps_checked=False))):
            path, unused = typed.typed_fixture(self.folder/f'badtyped{i}', 7599, 7600,
                                               version=2, mutate=mutate)
            with self.assertRaises(ValueError):
                self.extend([path])

    def test_zero_only_and_invalid_zero_semantics(self):
        path = self.fake_band('zero.json', zero=True)
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            result = self.extend([path])
        self.assertTrue(result['added_events_zero'])
        self.assertIsNone(result['added_events_union_log2_upper_rational'])
        def invalid(*args):
            row = self.fake_event(*args)
            row['log2_upper_rational'] = '-200'
            return row
        with patch.object(mode, '_band_event', side_effect=invalid):
            with self.assertRaisesRegex(ValueError, 'Zero event'):
                self.extend([path])

    def test_closing_input_and_source_changes_rejected(self):
        path = self.fake_band('new.json')
        read = Path.read_bytes
        calls = 0
        def changed_at_close(p):
            nonlocal calls
            raw = read(p)
            if p == self.parent:
                calls += 1
                if calls > 1:
                    return raw+b' '
            return raw
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            with patch.object(Path, 'read_bytes', changed_at_close):
                with self.assertRaisesRegex(RuntimeError, 'immutable band artifact changed'):
                    self.extend([path])
            bad = dict(self.pins)
            bad[Path(mode.__file__).name] = '0'*64
            with patch.object(mode, 'source_pins', return_value=bad):
                with self.assertRaisesRegex(RuntimeError, 'pinned weight-band source changed'):
                    self.extend([path])


if __name__ == '__main__':
    unittest.main()
