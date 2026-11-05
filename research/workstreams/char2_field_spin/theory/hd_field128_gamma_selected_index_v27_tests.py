"""Memory-backed bookkeeping fixtures: no files, coefficients, or proof runs."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v27 as mode


class IndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pins = mode.source_pins()

    def setUp(self):
        self.root = (mode.DIRECTORY/'__v27_memory_only__').resolve()
        self.files = {}
        self.read_counts = {}
        self.change_on_close = None
        self.original_read = Path.read_bytes
        self.original_stat = Path.stat
        self.original_is_file = Path.is_file
        for name, method in [('read_bytes', self.read), ('stat', self.stat), ('is_file', self.is_file)]:
            def dispatch(path, *args, method=method, **kwargs):
                return method(path, *args, **kwargs)
            mocked = patch.object(Path, name, dispatch)
            mocked.start(); self.addCleanup(mocked.stop)
        patcher = patch.object(mode, 'source_pins', side_effect=lambda: dict(self.pins))
        patcher.start(); self.addCleanup(patcher.stop)
        self.event = self.put('prior.json', dict(synthetic=True))
        self.prior_audit = self.put('prior_audit.json', dict(synthetic_audit=True))
        source_path = self.root/'prior.py'
        raw = b'# Synthetic inherited source.\n'
        self.files[str(source_path)] = raw
        self.source = dict(path=str(source_path), sha256=sha256(raw).hexdigest(), kind='python-source')
        entries = [dict(self.event, weight_range=[r, r], log2_upper_rational='-100') for r in range(1, 252)]
        entries += [dict(self.event, weight_range=[252, 7927], log2_upper_rational='-100', audit=self.prior_audit),
                    dict(self.event, weight_range=[8192, 8192], zero=True, log2_upper_rational=None)]
        self.parent = dict(schema='field128-gamma-partial-cover-index-v26',
            construction=mode.midpoint.previous.EXPECTED['construction'], field_order=1 << 128,
            source_dimension=1 << 20, output_length=1 << 21, target_distance='2/5',
            whole_code_certificate=False, directed_aggregation=True, entries=entries,
            receipt_count=253, covered_weight_count=7928, pending_weight_count=264,
            covered_weight_ranges=[[1, 7927], [8192, 8192]], pending_weight_ranges=[[7928, 8191]],
            covered_events_union_log2_upper_rational='-90', aggregation_source_sha256={},
            checked_non_json_file_pins=[self.source], producer_receipt_counts={},
            latest_explicit_selected_events=[dict(old_selected=True)],
            new_selected_rows_and_prefix_reaggregated=True)
        self.parent_pin = self.put('parent.json', self.parent)

    def read(self, path):
        key = str(path)
        if key not in self.files:
            return self.original_read(path)
        self.read_counts[key] = self.read_counts.get(key, 0)+1
        raw = self.files[key]
        return raw+b' ' if key == self.change_on_close and self.read_counts[key] > 1 else raw

    def stat(self, path, *args, **kwargs):
        if str(path) in self.files:
            return SimpleNamespace(st_size=len(self.files[str(path)]))
        return self.original_stat(path, *args, **kwargs)

    def is_file(self, path):
        return str(path) in self.files or self.original_is_file(path)

    def put(self, name, value):
        path = self.root/name
        raw = (json.dumps(value, sort_keys=True)+'\n').encode()
        self.files[str(path)] = raw
        return dict(path=str(path), sha256=sha256(raw).hexdigest())

    def event_at(self, first=7928, last=None, *, bound='-100', zero=False):
        last = first if last is None else last
        pin = self.put(f'band{first}_{last}.json', dict(weight_range=[first, last],
            zero=zero, log2_upper_rational=None if zero else bound))
        return pin['path']

    @staticmethod
    def fake_event(path, loader, pins, source_inputs):
        row, pin = loader.load(path, role='memory-only-band-contract')
        lo, hi = row['weight_range']
        return dict(pin, **row, weights_evaluated=hi-lo+1)

    def extend(self, paths):
        with patch.object(mode, 'PARENT_SHA256', self.parent_pin['sha256']):
            return mode.extend_v26(self.parent_pin['path'], paths)

    def test_parent_prefix_pins_fresh_union_and_typed_delegation(self):
        raw = self.files[self.parent_pin['path']]
        path = self.event_at(7928, 7930)
        with patch.object(mode.bands, '_band_event', side_effect=self.fake_event) as validator:
            result = self.extend([path])
        validator.assert_called_once()
        self.assertEqual(result['entries'][:253], self.parent['entries'])
        self.assertEqual(self.files[self.parent_pin['path']], raw)
        self.assertEqual((result['receipt_count'], result['covered_weight_count'], result['pending_weight_count']),
                         (254, 7931, 261))
        self.assertEqual(result['parent_index']['sha256'], self.parent_pin['sha256'])
        self.assertIn(self.source, result['checked_non_json_file_pins'])
        self.assertEqual(result['aggregation_source_sha256'][Path(mode.__file__).name],
                         self.pins[Path(mode.__file__).name])
        self.assertEqual(result['covered_events_union_log2_upper_rational'],
                         str(mode.index.union_ceiling(mode.frozen._bounds(result['entries']), 384)))
        self.assertEqual(result['latest_explicit_selected_events'], [])
        self.assertEqual(result['new_audit_schemas_allowed'], list(mode.bands.AUDIT_SCHEMAS))
        self.assertFalse(result['new_selected_rows_and_prefix_reaggregated'])
        for field in ('whole_code_certificate', 'coefficient_replay', 'normalized_regional_replay', 'transfer_recomputed'):
            self.assertFalse(result[field])

    def test_limits_and_full_coverage_never_claim_whole_code(self):
        for invalid in ([], [None]*65, 'audit.json', iter([])):
            with self.assertRaisesRegex(ValueError, 'one to 64'):
                self.extend(invalid)
        with patch.object(mode.bands, '_band_event', side_effect=self.fake_event):
            many = self.extend([self.event_at(r) for r in range(7928, 7992)])
            self.assertEqual(many['added_event_count'], 64)
            result = self.extend([self.event_at(r, min(r+63, 8191)) for r in range(7928, 8192, 64)])
        self.assertEqual(result['covered_weight_count'], 8192)
        self.assertEqual(result['pending_weight_ranges'], [])
        self.assertFalse(result['whole_code_certificate'])
        self.assertIn('no whole-code certificate', result['scope'])

    def test_wrong_parent_metadata_coverage_union_and_source_rejected(self):
        path = self.event_at()
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            mode.extend_v26(self.parent_pin['path'], [path])
        good = deepcopy(self.parent)
        mutations = [lambda p: p.update(whole_code_certificate=True),
            lambda p: p.update(receipt_count=252), lambda p: p['entries'].pop(),
            lambda p: p.update(pending_weight_ranges=[]),
            lambda p: p.update(covered_events_union_log2_upper_rational='-100'),
            lambda p: p['aggregation_source_sha256'].update({Path(mode.__file__).name: '0'*64})]
        for mutate in mutations:
            self.parent = deepcopy(good); mutate(self.parent)
            self.parent_pin = self.put('parent.json', self.parent)
            with self.assertRaises(ValueError):
                self.extend([path])

    def test_overlap_and_global_budget_rejected(self):
        path = self.event_at()
        with patch.object(mode.bands, '_band_event', side_effect=self.fake_event):
            for paths in ([path, path], [self.event_at(7927, 7928)]):
                with self.assertRaisesRegex(ValueError, 'overlaps'):
                    self.extend(paths)
            with self.assertRaisesRegex(ValueError, '54-bit budget'):
                self.extend([self.event_at(bound='-54')])

    def test_zero_events_have_null_logs_and_positive_union_only(self):
        path = self.event_at(zero=True)
        with patch.object(mode.bands, '_band_event', side_effect=self.fake_event):
            result = self.extend([path])
        self.assertTrue(result['added_events_zero'])
        self.assertIsNone(result['added_events_union_log2_upper_rational'])
        def invalid(*args):
            result = self.fake_event(*args); result['log2_upper_rational'] = '-100'
            return result
        with patch.object(mode.bands, '_band_event', side_effect=invalid):
            with self.assertRaisesRegex(ValueError, 'Zero event'):
                self.extend([path])

    def test_untyped_audit_rejected_before_numerical_validation(self):
        path = self.put('not_typed.json', dict(schema='gamma-midpoint-selected-receipts-audit-v1'))['path']
        with self.assertRaisesRegex(ValueError, 'typed v1/v2'):
            self.extend([path])

    def test_input_source_mutations_and_closing_checks_fail_closed(self):
        path = self.event_at()
        with patch.object(mode.bands, '_band_event', side_effect=self.fake_event):
            self.files[self.source['path']] += b'# changed\n'
            with self.assertRaisesRegex(ValueError, 'source byte pin mismatch'):
                self.extend([path])
            self.files[self.source['path']] = b'# Synthetic inherited source.\n'
            self.read_counts.clear(); self.change_on_close = self.parent_pin['path']
            with self.assertRaisesRegex(RuntimeError, 'immutable band artifact changed'):
                self.extend([path])
            self.read_counts.clear(); self.change_on_close = self.source['path']
            with self.assertRaisesRegex(ValueError, 'Pinned non-JSON input changed'):
                self.extend([path])
            self.change_on_close = None
            bad = dict(self.pins); bad[Path(mode.__file__).name] = '0'*64
            with patch.object(mode, 'source_pins', return_value=bad):
                with self.assertRaisesRegex(RuntimeError, 'pinned weight-band source changed'):
                    self.extend([path])


if __name__ == '__main__':
    unittest.main()
