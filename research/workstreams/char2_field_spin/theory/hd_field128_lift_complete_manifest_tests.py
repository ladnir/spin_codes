"""Synthetic scalar/coverage checks; no regional coefficient or distance gate."""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_lift_complete_manifest as manifest
from hd_field128_lift_partial_index_tests import alias_receipt


index = manifest.index


def full_fixture(closures, pin):
    """Mock coefficient evidence; only deterministic cover logic is exercised."""
    provenance = dict(pin, kind='synthetic-mock')
    events = [dict(active_groups=r, log2_upper_rational='-100', provenance=dict(provenance, row_index=r-1))
              for r in range(1, index.TERMINAL_RANGE[0])]
    selected = [dict(row, used_in_union=True, passes_individual_target=True,
                     compatible_candidate_count=1) for row in events]
    terminal = dict(active_groups=None, log2_upper_rational='-200', provenance=provenance,
                    weight_range=index.TERMINAL_RANGE, used_in_union=True,
                    passes_individual_target=True, compatible_candidate_count=1)
    events.append(dict(manifest.event_core(terminal), weight_range=index.TERMINAL_RANGE))
    return dict(index.Validator(manifest.DISTANCE).full(), schema='field128-lift-partial-index-v1',
                indexed_weight_domain=[1, index.GROUPS], all_weights_covered=True,
                passing_weight_count=index.GROUPS, passing_weight_ranges=[[1, index.GROUPS]],
                unresolved_weight_ranges=[], unresolved_weights=[], missing_weight_ranges=[], missing_weights=[],
                known_above_target_weight_ranges=[], known_above_target_weights=[],
                whole_code_threshold_conditions_met=True, selected_union_meets_aggregate_target=True,
                selected_per_weight=selected, terminal_event=terminal, union_events=events,
                union_event_count=len(events), terminal_union_event_count=1,
                individual_event_failure_bits=80, aggregate_failure_bits=40, aggregation_precision=256,
                selected_union_log2_upper_rational='-87', explicit_inputs=[pin], receipt_file_pins=[pin],
                source_sha256=closures[manifest.INDEX], checked_producer_closures={},
                checked_source_sha256=closures[manifest.INDEX])


class CompleteManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reviewed = manifest.reviewed_closures()

    def setUp(self):
        self.tmp = TemporaryDirectory(prefix='literal35-manifest-test-')
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.source = self.save('mock-receipt.json', dict(scope='mock scalar evidence only'))
        raw = self.source.read_bytes()
        self.pin = dict(path=str(self.source.resolve()), sha256=sha256(raw).hexdigest(), bytes=len(raw))

    def save(self, name, value):
        path = self.directory/name
        path.write_text(json.dumps(value), encoding='utf-8')
        return path

    def fixture(self):
        return full_fixture(self.reviewed[0], self.pin)

    def test_complete_mock_reconstruction_union_and_honest_scope(self):
        row = self.fixture()
        path = self.save('index.json', row)
        original = path.read_bytes()
        with patch.object(manifest, 'reviewed_closures', return_value=self.reviewed), \
             patch.object(index, 'build_index', return_value=row) as reconstruct:
            result = manifest.build_manifest(path)
        reconstruct.assert_called_once_with([str(self.source.resolve())], distance=F(7, 20),
                                             target_bits=80, aggregate_bits=40, precision=256)
        self.assertEqual(path.read_bytes(), original)
        self.assertTrue(result['whole_code_certificate'])
        self.assertFalse(result['full_coefficient_replay'])
        self.assertFalse(result['verification']['receipt_authentication'])
        self.assertFalse(result['original_index_whole_code_certificate'])
        self.assertEqual(result['complete_weight_count'], 8192)
        self.assertEqual(result['terminal_union_event_count'], 1)
        self.assertEqual(result['union_event_count'], 5854)
        self.assertEqual([entry['precision'] for entry in result['directed_union_replays']], [256, 384])
        self.assertLess(F(result['selected_union_log2_upper_rational']), -87)

    def test_real_partial_index_reconstruction_is_rejected(self):
        receipt = self.save('one-row.json', alias_receipt(1, [-100]))
        partial = index.build_index([receipt])
        path = self.save('partial.json', partial)
        with patch.object(manifest, 'reviewed_closures', return_value=self.reviewed):
            with self.assertRaisesRegex(ValueError, 'complete index coverage'):
                manifest.build_manifest(path)

    def test_duplicate_terminal_and_tail_overlap_rejected(self):
        row = self.fixture()
        row['union_events'].append(deepcopy(row['union_events'][-1]))
        row['union_event_count'] += 1
        row['terminal_union_event_count'] += 1
        with self.assertRaisesRegex(ValueError, 'counted twice'):
            manifest.validate_cover(row)
        row = self.fixture()
        tail_row = dict(active_groups=5854, log2_upper_rational='-300', provenance=self.pin,
                        used_in_union=True, passes_individual_target=True)
        row['selected_per_weight'].append(tail_row)
        row['union_events'].append(manifest.event_core(tail_row))
        with self.assertRaisesRegex(ValueError, 'overlapping'):
            manifest.validate_cover(row)

    def test_fake_complete_flag_does_not_fill_gap(self):
        row = self.fixture()
        row['union_events'].pop(100)
        with self.assertRaisesRegex(ValueError, 'does not cover all'):
            manifest.validate_cover(row)

    def test_unreviewed_and_conflicting_source_pins_rejected(self):
        row = self.fixture()
        row['checked_producer_closures']['unreviewed.py'] = {}
        with self.assertRaisesRegex(ValueError, 'Unreviewed'):
            manifest.validate_reviewed_sources(row, *self.reviewed)
        row = self.fixture()
        row['checked_source_sha256'] = dict(row['checked_source_sha256'], extra='0'*64)
        with self.assertRaisesRegex(ValueError, 'source union'):
            manifest.validate_reviewed_sources(row, *self.reviewed)
        original = index.source_closure
        def altered(name):
            return dict(original(name), **{name: '0'*64})
        with patch.object(index, 'source_closure', side_effect=altered):
            with self.assertRaisesRegex(ValueError, 'Reviewed producer'):
                manifest.reviewed_closures()

    def test_stale_index_wrong_construction_and_input_hash_rejected(self):
        row = self.fixture()
        stale = deepcopy(row)
        stale['union_events'][0]['log2_upper_rational'] = '-999'
        with patch.object(manifest, 'reviewed_closures', return_value=self.reviewed), \
             patch.object(index, 'build_index', return_value=row):
            with self.assertRaisesRegex(ValueError, 'differs from reconstruction'):
                manifest.build_manifest(self.save('stale.json', stale))
        wrong = self.fixture()
        wrong['construction'] = 'full-RS128-48-not-the-literal-construction'
        with self.assertRaisesRegex(ValueError, 'construction'):
            manifest.validate_cover(wrong)
        row = self.fixture()
        row['explicit_inputs'][0] = dict(row['explicit_inputs'][0], sha256='0'*64)
        with patch.object(manifest, 'reviewed_closures', return_value=self.reviewed):
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                manifest.build_manifest(self.save('wrong-pin.json', row))

    def test_union_replay_exact_multiplicity_and_ceiling(self):
        for precision in manifest.PRECISIONS:
            row = manifest.replay_union([F(-100), F(-100)], precision)
            self.assertLessEqual(F(-99), F(row['log2_upper_rational']))
            self.assertLessEqual(F(row['log2_upper_rational']), F(-99)+F(1, 10**12))
            self.assertEqual(row['distinct_ceilings'], 1)
            with ctx.workprec(precision):
                self.assertTrue(index.rational_arb(F(row['log2_upper_rational'])) >= arb(row['log2_upper_enclosure']).upper())
        row = self.fixture()
        row['individual_event_failure_bits'] = 40
        for item in row['selected_per_weight']:
            item['log2_upper_rational'] = '-45'
        for item in row['union_events'][:-1]:
            item['log2_upper_rational'] = '-45'
        with patch.object(manifest, 'reviewed_closures', return_value=self.reviewed), \
             patch.object(index, 'build_index', return_value=row):
            with self.assertRaisesRegex(ValueError, 'fails the40bit'):
                manifest.build_manifest(self.save('weak.json', row))


if __name__ == '__main__':
    unittest.main()
