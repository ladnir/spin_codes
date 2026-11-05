"""Tiny saved-artifact tests; no coefficient or FFT construction."""

from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import ctx

import hd_field128_gamma_compact_index_extend as extender
from hd_field128_gamma_near_face_index_update_tests import audited_fixture, prior_index
from hd_field128_gamma_near_face_band_tests import FakeCache


class ExtendTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.summary, cls.audit, cls.audit_pin = audited_fixture(cls.root)
        prior = prior_index(cls.root)
        cls.parent = json.loads(Path(prior['path']).read_bytes())
        cls.parent.update(schema='field128-gamma-partial-cover-index-v12', directed_aggregation=True)
        cls.parent_pin = extender.band.write_new(cls.root/'v12.json', cls.parent)
        cls.events = [(cls.root/'band'/'summary.json', cls.audit_pin['path'])]

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_version_independent_append_reuses_frozen_resolver(self):
        original = extender.frozen._event
        with patch.object(extender.frozen, '_event', wraps=original) as event:
            result = extender.extend_compact_index(self.parent_pin['path'], self.events)
        event.assert_called_once()
        self.assertEqual(event.call_args.args[0], 'compact-band')
        self.assertEqual(result['schema'], 'field128-gamma-partial-cover-index-v13')
        self.assertEqual(result['covered_weight_ranges'], [[1, 3], [9, 9]])
        self.assertEqual(result['covered_weight_count'], 4)
        self.assertEqual(result['pending_weight_count'], 8188)
        self.assertEqual(result['added_event_count'], 1)
        self.assertEqual(result['added_weight_count'], 3)
        self.assertFalse(result['whole_code_certificate'])
        self.assertIn(Path(extender.__file__).name, result['aggregation_source_sha256'])

    def test_duplicate_new_events_rejected(self):
        with self.assertRaisesRegex(ValueError, 'overlaps'):
            extender.extend_compact_index(self.parent_pin['path'], self.events*2)

    def test_parent_scope_counts_and_source_conflicts_rejected(self):
        for number, changes in enumerate(({'schema': 'field128-gamma-partial-cover-index-v11'},
                {'target_distance': '9/20'}, {'pending_weight_count': 0},
                {'aggregation_source_sha256': {Path(extender.__file__).name: '0'*64}})):
            value = deepcopy(self.parent)
            value.update(changes)
            pin = extender.band.write_new(self.root/f'bad-parent-{number}.json', value)
            with self.assertRaises(ValueError):
                extender.extend_compact_index(pin['path'], self.events)
        for events in ([], self.events*33):
            with self.assertRaises(ValueError):
                extender.extend_compact_index(self.parent_pin['path'], events)

    def test_exit_source_check_cannot_be_skipped(self):
        with patch.object(extender.band.selected, 'check_source_pins', side_effect=RuntimeError('source mutation')):
            with self.assertRaisesRegex(RuntimeError, 'source mutation'):
                extender.extend_compact_index(self.parent_pin['path'], self.events)

    def test_inline_proposal_resolves_original_bytes_and_frozen_row_auditor(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            context = json.loads(Path(self.summary['context']['path']).read_bytes())
            source = extender.band.write_new(root/'literal-source.json', {'synthetic_literal_source': True})
            proposal = dict(schema='field128-gamma-exact-proposal-v1',
                construction=extender.frozen.audited.EXPECTED['construction'], state_dimension=30,
                target_distance='2/5', calibration_active_groups=1, marker_count=1, source=source,
                witness={key: context['witness'][key] for key in ('z', 'alpha', 'markers', 'grid', 'precision', 'packet_mantissa_bits')},
                directed=False, whole_code_certificate=False)
            proposal_pin = extender.band.write_new(root/'proposal.json', proposal)
            bank = object.__new__(extender.band.provider.AliasBank)
            bank.precision, bank.construction = 256, context['construction']
            bank.metadata, bank.check_pins = lambda: deepcopy(context), lambda: None
            with ctx.workprec(256):
                bank.constant = extender.frozen.audited.outer_prefix(context, 0)
                bank.slope = extender.frozen.audited.outer_prefix(context, 1)-bank.constant
            FakeCache.fail_weight = None
            with patch.object(extender.band, 'RegionalCache', FakeCache):
                summary = extender.band.band_receipt(bank, 1, 2, output_dir=root/'band', witness_source=proposal)
            original = Path(summary['context']['path']).read_bytes()
            with patch.object(extender.frozen.audited, 'audit_rows', wraps=extender.frozen.audited.audit_rows) as rows:
                report = extender.audit_inline_band(root/'band'/'summary.json', proposal_pin['path'])
            self.assertEqual(rows.call_count, 2)
            self.assertEqual(report['audit_result'], 'PASS')
            self.assertEqual(report['resolved_inline_proposal']['sha256'], proposal_pin['sha256'])
            self.assertEqual(Path(summary['context']['path']).read_bytes(), original)
            saved = extender.band.write_new(root/'audit.json', report)
            event, pins = extender._event(root/'band'/'summary.json', saved['path'], extender.index.Validator(extender.F(2, 5)))
            self.assertEqual(event['weight_range'], [1, 2])
            self.assertEqual(pins, extender.source_pins())
            proposal['calibration_active_groups'] = 2
            different = extender.band.write_new(root/'different-proposal.json', proposal)
            with self.assertRaisesRegex(ValueError, 'differs from immutable proposal'):
                extender.audit_inline_band(root/'band'/'summary.json', different['path'])


if __name__ == '__main__':
    unittest.main()
