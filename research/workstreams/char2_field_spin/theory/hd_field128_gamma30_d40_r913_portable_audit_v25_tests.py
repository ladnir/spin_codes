"""Saved-fixture and mocked arithmetic tests; no profile-count or coefficient run."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import hd_field128_gamma30_d40_r913_portable_audit_v25 as mode

RUN = mode.HERE/'hd_field128_gamma30_d40_weights913_913_midpoint_v25_20261010'


class AdapterTests(unittest.TestCase):
    def test_explicit_theory_relative_resolution(self):
        link = dict(path='some_run/proposal.json', sha256='a'*64)
        self.assertEqual(mode.portable_path(link), mode.HERE/'some_run/proposal.json')
        for name in ('../p.json', '/tmp/p.json', 'C:/p.json', 'a\\p.json', 'a//p.json', './p.json'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                mode.portable_path(dict(link, path=name))

    def test_complete_metadata_and_source_pin_guards(self):
        record = json.loads((RUN/'summary.json').read_bytes())
        mode.public_metadata(record)
        for key, value in (('active_groups', 914), ('profile_pairs_evaluated', 208848),
                           ('complete_selected_weight', False), ('log2_upper_rational', '-53')):
            changed = dict(record, **{key: value})
            with self.subTest(key=key), self.assertRaises(ValueError):
                mode.public_metadata(changed)
        with self.assertRaises(ValueError):
            mode.source_maps({'source_sha256': {mode.NAME: '0'*64}}, {})

    def test_context_reconstructs_without_receipt_rewrite(self):
        record = json.loads((RUN/'summary.json').read_bytes())
        bank = json.loads((RUN/'bank.json').read_bytes())
        saved = json.loads((RUN/'context.json').read_bytes())
        context = mode.expected_context(record, bank)
        self.assertEqual(context, saved['context'])
        digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        self.assertEqual(digest, record['context_sha256'])
        changed = deepcopy(bank); changed['witness']['markers'][0][0] = '1'
        self.assertNotEqual(mode.expected_context(record, changed), context)

    def test_full_adapter_uses_frozen_arithmetic_once_and_preserves_bytes(self):
        before = (RUN/'summary.json').read_bytes()
        stub = dict(active_groups=913, even_rows=457, profile_pairs=208849,
                    saved_log2_upper_rational='-58945660269253/1000000000000')
        with patch.object(mode.frozen, 'audit_record', return_value=stub) as arithmetic:
            result = mode.audit_summary(RUN/'summary.json')
        arithmetic.assert_called_once()
        self.assertEqual((RUN/'summary.json').read_bytes(), before)
        self.assertTrue(result['external_checkpoint_context_reconstructed'])
        self.assertEqual(result['external_rows_checked'], 457)
        self.assertEqual(result['selected_weights'], [913])
        self.assertFalse(result['coefficient_replay'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertIn(str((RUN/'row_0000.json').relative_to(mode.HERE)).replace('\\','/'),
                      result['external_input_sha256'])

    def test_changed_external_row_and_budget_fail_before_arithmetic(self):
        original = mode.frozen.index.Validator.load
        def changed(loader, path, **options):
            value, pin = original(loader, path, **options)
            if Path(path).name == 'row_0000.json':
                value = dict(value, active_groups=914)
            return value, pin
        with patch.object(mode.frozen.index.Validator, 'load', changed), \
                patch.object(mode.frozen, 'audit_record', side_effect=AssertionError('No arithmetic')):
            with self.assertRaises(ValueError):
                mode.audit_summary(RUN/'summary.json')
        for seconds in (0, 121, True):
            with self.assertRaises(ValueError):
                mode.audit_summary(RUN/'summary.json', max_seconds=seconds)


if __name__ == '__main__':
    unittest.main()
