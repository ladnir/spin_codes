"""Tiny synthetic contracts; no production index or coefficient bank is evaluated."""

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v16_tests as fixtures
import hd_field128_gamma_selected_index_v22_typed as mode


def parent_fixture(folder, *, mutate=None):
    event = fixtures.write(folder/'prior.json', dict(synthetic_fixture=True))
    audit = fixtures.write(folder/'prior_audit.json', dict(synthetic_prior_audit=True))
    entries = [dict(event, weight_range=[r, r], log2_upper_rational='-100') for r in range(1, 228)]
    entries[-1].update(zero=True, log2_upper_rational=None, audit=audit)
    parent = dict(schema='field128-gamma-partial-cover-index-v21',
        construction=mode.midpoint.previous.EXPECTED['construction'], field_order=1 << 128,
        source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
        whole_code_certificate=False, directed_aggregation=True, entries=entries,
        covered_weight_ranges=[[1, 227]], covered_weight_count=227,
        pending_weight_ranges=[[228, 8192]], pending_weight_count=7965, receipt_count=227,
        covered_events_union_log2_upper_rational='-90',
        producer_receipt_counts={'hd_field128_gamma_midpoint_weight_bands': 22},
        latest_audited_multiweight_events=[entries[-1]], aggregation_source_sha256={})
    if mutate:
        mutate(parent)
    pin = fixtures.write(folder/'parent.json', parent)
    return Path(pin['path']), pin['sha256']


def typed_fixture(folder, first=228, last=229, *, zero=False, mutate=None, version=1):
    legacy, summary = fixtures.band_fixture(folder, first, last, zero=zero)
    report = json.loads(legacy.read_bytes())
    helper_path = Path(mode.previous.__file__).resolve()
    helper = dict(path=str(helper_path), sha256=sha256(helper_path.read_bytes()).hexdigest(),
                  kind='python-source')
    provider = mode.audited if version == 1 else mode.audited_v2
    report.update(schema=f'gamma-midpoint-weight-band-typed-independent-audit-v{version}',
        source_sha256=provider.source_pins(),
        checked_source_sha256={**provider.source_pins(), helper_path.name: helper['sha256']},
        checked_non_json_file_pins=[helper], typed_provenance_checked=True,
        source_and_input_bytes_checked_at_exit=True, cell_transfer_recomputed=False,
        execution_driver_receipt=None,
        frozen_numerical_auditor_sha256=mode.audited.PREVIOUS_SHA256)
    if version == 2:
        report.update(nested_source_maps_checked=True,
                      frozen_typed_v1_auditor_sha256=mode.audited_v2.TYPED_V1_SHA256)
    if mutate:
        mutate(report)
    pin = fixtures.write(folder/'typed_audit.json', report)
    return Path(pin['path']), summary


class IndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Compute the real immutable closure once for the synthetic tests.
        cls.pins = mode.source_pins()

    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.folder = Path(self.temporary.name)
        self.parent, self.parent_hash = parent_fixture(self.folder/'old')
        self.pin_patch = patch.object(mode, 'source_pins', side_effect=lambda: dict(self.pins))
        self.pin_patch.start()

    def tearDown(self):
        self.pin_patch.stop()
        self.temporary.cleanup()

    def extend(self, paths):
        with patch.object(mode, 'PARENT_SHA256', self.parent_hash):
            return mode.extend_v21(self.parent, paths)

    def fake_band(self, name, first, last=None, *, bound='-200', zero=False):
        return Path(fixtures.write(self.folder/name, dict(
            weight_range=[first, first if last is None else last],
            zero=zero, log2_upper_rational=None if zero else bound))['path'])

    @staticmethod
    def fake_event(path, loader, pins, source_inputs):
        # Only adapter mechanics use this stub; the first test exercises the
        # unchanged real v16 validator with synthetic complete cell receipts.
        row, pin = loader.load(path, role='synthetic-band-adapter-test')
        lo, hi = row['weight_range']
        return dict(pin, **row, weights_evaluated=hi-lo+1)

    def test_real_frozen_validator_retains_227_events_and_both_zero_semantics(self):
        first, unused = typed_fixture(self.folder/'first', 228, 229)
        zero, unused = fixtures.band_fixture(self.folder/'zero', 231, 232, zero=True)
        raw = self.parent.read_bytes()
        prior = json.loads(raw)
        result = self.extend([first, zero])
        self.assertEqual(result['schema'], 'field128-gamma-partial-cover-index-v22-typed')
        self.assertEqual(result['entries'][:227], prior['entries'])
        self.assertEqual(self.parent.read_bytes(), raw)
        self.assertEqual(result['receipt_count'], 229)
        self.assertEqual(result['added_event_count'], 2)
        self.assertEqual(result['added_weight_count'], 4)
        self.assertEqual(result['covered_weight_ranges'], [[1, 229], [231, 232]])
        self.assertEqual(result['pending_weight_ranges'], [[230, 230], [233, 8192]])
        self.assertTrue(result['entries'][-1]['zero'])
        self.assertIsNone(result['entries'][-1]['log2_upper_rational'])
        expected = mode.index.union_ceiling(mode.frozen._bounds(result['entries']), 384)
        self.assertEqual(result['covered_events_union_log2_upper_rational'], str(expected))
        self.assertEqual(result['covered_events_failure_bits'], 54)
        self.assertTrue(result['meets_covered_events_target'])
        self.assertEqual(result['producer_receipt_counts']['hd_field128_gamma_midpoint_weight_bands'], 24)
        self.assertEqual(result['parent_index']['sha256'], self.parent_hash)
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['coefficient_replay'])
        self.assertFalse(result['normalized_regional_replay'])
        for version in (16, 17, 18, 19, 20, 21, 22):
            self.assertIn(f'hd_field128_gamma_selected_index_v{version}.py',
                          result['aggregation_source_sha256'])

    def test_zero_only_addition_and_invalid_zero_log(self):
        good = self.fake_band('zero.json', 228, 229, zero=True)
        bad = self.fake_band('badzero.json', 230)
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            result = self.extend([good])
            self.assertTrue(result['added_events_zero'])
            self.assertIsNone(result['added_events_union_log2_upper_rational'])
            self.assertEqual(result['covered_events_union_log2_upper_rational'],
                             str(mode.index.union_ceiling([F(-100)]*226, 384)))
        def invalid(path, loader, pins, source_inputs):
            row = self.fake_event(path, loader, pins, source_inputs)
            row['zero'] = True
            return row
        with patch.object(mode, '_band_event', side_effect=invalid):
            with self.assertRaisesRegex(ValueError, 'Zero event'):
                self.extend([bad])

    def test_explicit_64_band_limit(self):
        paths = [self.fake_band(f'audit{i}.json', 228+i) for i in range(64)]
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            result = self.extend(paths)
        self.assertEqual(result['added_event_count'], 64)
        self.assertEqual(result['receipt_count'], 291)
        self.assertEqual(result['covered_weight_ranges'], [[1, 291]])
        for invalid in ([], paths+[paths[0]], str(paths[0]), iter(paths)):
            with self.assertRaisesRegex(ValueError, 'explicit audited bands'):
                self.extend(invalid)

    def test_complete_integer_coverage_does_not_promote_whole_code_claim(self):
        def almost_complete(parent):
            parent['entries'][-1]['weight_range'] = [228, 8192]
            parent.update(covered_weight_ranges=[[1, 226], [228, 8192]],
                          covered_weight_count=8191, pending_weight_ranges=[[227, 227]],
                          pending_weight_count=1)
        self.parent, self.parent_hash = parent_fixture(
            self.folder/'almost_complete', mutate=almost_complete)
        last = self.fake_band('last.json', 227)
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            result = self.extend([last])
        self.assertEqual(result['covered_weight_ranges'], [[1, 8192]])
        self.assertEqual(result['covered_weight_count'], 8192)
        self.assertEqual(result['pending_weight_count'], 0)
        self.assertEqual(result['pending_weight_ranges'], [])
        self.assertFalse(result['whole_code_certificate'])
        self.assertIn('no whole-code certificate', result['scope'])

    def test_overlap_duplicate_pin_and_aggregate_budget_rejected(self):
        old = self.fake_band('overlap.json', 227, 228)
        new = self.fake_band('new.json', 228, 229)
        threshold = self.fake_band('threshold.json', 230, bound='-54')
        with patch.object(mode, '_band_event', side_effect=self.fake_event):
            for paths in ([old], [new, new]):
                with self.assertRaisesRegex(ValueError, 'overlaps'):
                    self.extend(paths)
            # A -54 admitted event plus nonzero prior events exceeds the total
            # -54 budget; the adapter must reject rather than round it down.
            with self.assertRaisesRegex(ValueError, 'union fails 54-bit'):
                self.extend([threshold])
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            mode.extend_v21(self.parent, [new])

    def test_parent_count_coverage_union_and_source_tampering(self):
        new = self.fake_band('new.json', 228)
        mutations = [
            lambda p: p.update(receipt_count=226),
            lambda p: p['entries'].pop(),
            lambda p: p.update(covered_weight_count=228),
            lambda p: p.update(covered_events_union_log2_upper_rational='-100'),
            lambda p: p['aggregation_source_sha256'].update(
                {'hd_field128_gamma_midpoint_weight_bands.py': 'wrong'}),
        ]
        for i, mutate in enumerate(mutations):
            self.parent, self.parent_hash = parent_fixture(self.folder/f'bad{i}', mutate=mutate)
            with self.assertRaises(ValueError):
                self.extend([new])

    def test_real_partial_receipt_and_changed_source_or_input_fail_closed(self):
        partial, unused = fixtures.band_fixture(self.folder/'partial', 228, 229,
            mutate=lambda row: row.update(complete_requested_band=False))
        with self.assertRaises(ValueError):
            self.extend([partial])
        new = self.fake_band('new.json', 228)
        original = Path.read_bytes
        prior_event = self.parent.parent/'prior.json'
        def changed(path):
            raw = original(path)
            return raw+b' ' if path == prior_event else raw
        with patch.object(Path, 'read_bytes', changed):
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                self.extend([new])
        # Exercise the unconditional final source-byte guard, not only merge
        # conflicts: one closure pin is deliberately incorrect in this mock.
        bad_pins = dict(self.pins)
        bad_pins['hd_field128_gamma_selected_index_v22_typed.py'] = '0'*64
        with patch.object(mode, 'source_pins', return_value=bad_pins):
            with patch.object(mode, '_band_event', side_effect=self.fake_event):
                with self.assertRaisesRegex(RuntimeError, 'pinned weight-band source changed'):
                    self.extend([new])


    def test_typed_identity_helper_kind_and_pin_are_strict(self):
        changes = [
            lambda r: r.update(schema='unreviewed-audit'),
            lambda r: r.update(typed_provenance_checked=False),
            lambda r: r.update(frozen_numerical_auditor_sha256='bad'),
            lambda r: r.update(source_sha256={}),
            lambda r: r['checked_non_json_file_pins'][0].update(kind='markdown-source'),
            lambda r: r['checked_non_json_file_pins'][0].update(sha256='0'*64),
            lambda r: r['checked_source_sha256'].pop('hd_field128_gamma_selected_index_v22.py'),
            lambda r: r['checked_non_json_file_pins'][0].update(path='relative.py'),
        ]
        for i, mutate in enumerate(changes):
            path, unused = typed_fixture(self.folder/f'typed_bad{i}', mutate=mutate)
            with self.assertRaises(ValueError, msg=str(i)):
                self.extend([path])

    def test_typed_zero_and_closing_byte_checks(self):
        path, unused = typed_fixture(self.folder/'typed_zero', zero=True)
        raw_parent = self.parent.read_bytes()
        result = self.extend([path])
        self.assertTrue(result['added_events_zero'])
        self.assertIsNone(result['added_events_union_log2_upper_rational'])
        self.assertEqual(result['non_json_input_file_pins_verified'], 1)
        self.assertEqual(result['entries'][-1]['receipt_kind'],
                         'typed-audited-midpoint-multiweight-cell-band')
        self.assertEqual(self.parent.read_bytes(), raw_parent)
        target = Path(mode.previous.__file__).resolve()
        original = Path.read_bytes
        calls = 0
        def changed_at_close(p):
            nonlocal calls
            raw = original(p)
            if p.resolve() == target:
                calls += 1
                if calls >= 2:
                    return raw+b' '
            return raw
        with patch.object(Path, 'read_bytes', changed_at_close):
            with self.assertRaises((ValueError, RuntimeError)):
                self.extend([path])

    def test_typed_json_and_byte_lanes_cannot_be_confused(self):
        # A JSON input in the source-byte lane must fail before union use.
        path, summary = typed_fixture(self.folder/'json_source',
            mutate=lambda r: r['checked_non_json_file_pins'].__setitem__(0,
                dict(path=str((self.folder/'json_source'/'summary.json').resolve()),
                     sha256='0'*64, kind='python-source')))
        with self.assertRaisesRegex(ValueError, 'Unsupported non-JSON'):
            self.extend([path])
        # Conversely, a Python helper in the JSON lane is not silently skipped.
        path, unused = typed_fixture(self.folder/'source_json',
            mutate=lambda r: r['checked_file_pins'].append(r['checked_non_json_file_pins'][0]))
        with self.assertRaises((ValueError, json.JSONDecodeError)):
            self.extend([path])

    def test_v1_v2_and_legacy_audits_have_explicit_distinct_identities(self):
        first, unused = typed_fixture(self.folder/'v1', 228, 229)
        second, unused = typed_fixture(self.folder/'v2', 230, 231, version=2)
        third, unused = fixtures.band_fixture(self.folder/'legacy', 232, 233, zero=True)
        result = self.extend([first, second, third])
        self.assertEqual(result['covered_weight_ranges'], [[1, 233]])
        self.assertEqual(result['added_event_count'], 3)
        self.assertEqual(result['non_json_input_file_pins_verified'], 1)
        self.assertIn(Path(mode.audited_v2.__file__).name, result['aggregation_source_sha256'])
        for i, change in enumerate((
                lambda r: r.update(nested_source_maps_checked=False),
                lambda r: r.update(frozen_typed_v1_auditor_sha256='bad'),
                lambda r: r.update(source_sha256=mode.audited.source_pins()))):
            path, unused = typed_fixture(self.folder/f'bad_v2_{i}', version=2, mutate=change)
            with self.assertRaises(ValueError):
                self.extend([path])


if __name__ == '__main__':
    unittest.main()
