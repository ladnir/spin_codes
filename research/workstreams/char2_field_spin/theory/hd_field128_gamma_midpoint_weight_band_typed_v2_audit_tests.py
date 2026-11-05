"""Tiny v2 nested-source compatibility contracts; no production arithmetic."""
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import hd_field128_gamma_midpoint_weight_band_typed_v2_audit as mode
import hd_field128_gamma_midpoint_weight_band_typed_audit_tests as prior_tests


class TypedV2Tests(unittest.TestCase):
    def test_recursive_maps_preserve_leaf_paths(self):
        with tempfile.TemporaryDirectory() as root:
            p = prior_tests.source_file(Path(root), 'helper.py', 'VALUE = 1\n')
            nested = dict(replay_source_sha256=dict(provider={p['path']: p['sha256']},
                                                   band=dict(deeper={p['path']: p['sha256']})))
            sources = {}
            mode._provenance(nested, mode.Inputs(), sources)
            self.assertEqual(sources, {p['path']: p['sha256']})
            with self.assertRaises(ValueError):
                mode._provenance(dict(source_sha256={p['path']: '0'*64}), mode.Inputs(), sources)

    def test_malformed_nested_leaves_and_maps_rejected(self):
        bad = [[], None, 1, 'abc', dict(provider=[]), dict(provider=dict(file=1)),
               dict(provider=dict(file=None)), dict(provider=dict(file='not-a-digest')),
               dict(provider=dict(file='A'*64)), {'': {}}, {3: {}}]
        for value in bad:
            with self.assertRaises(ValueError, msg=repr(value)):
                mode._provenance(dict(replay_source_sha256=value), mode.Inputs(), {})

    def test_nested_json_and_python_link_types_and_tamper(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            helper = prior_tests.source_file(directory, 'helper.py', 'VALUE = 1\n')
            child = mode.producer.profiles.write_new(directory/'child.json',
                dict(driver=helper, replay_source_sha256=dict(provider={helper['path']: helper['sha256']})))
            inputs, sources = mode.Inputs(), {}
            mode._provenance(child, inputs, sources)
            self.assertEqual(len(inputs.pins), 1)
            self.assertEqual(len(inputs.non_json_pins), 1)
            target = Path(helper['path'])
            original = Path.read_bytes
            with patch.object(Path, 'read_bytes', lambda p: original(p)+b' ' if p == target else original(p)):
                with self.assertRaises(ValueError):
                    inputs.check()

    def test_explicit_source_identity_and_unchanged_frozen_numerics(self):
        pins = mode.source_pins()
        mode.check_source_pins(pins)
        self.assertEqual(pins[Path(mode.typed_v1.__file__).name], mode.TYPED_V1_SHA256)
        self.assertEqual(pins[Path(mode.previous.__file__).name], mode.PREVIOUS_SHA256)
        self.assertIs(mode.previous.audit_numeric, mode.typed_v1.previous.audit_numeric)
        self.assertIs(mode.Inputs, mode.typed_v1.Inputs)
        record = prior_tests.fixtures.synthetic()
        check = mode.previous.audit_numeric(record, 4, mode.base.Budget())
        self.assertEqual(check['union'], record['log2_upper_rational'])

    def test_public_resolver_v2_identity_with_nested_map(self):
        result = prior_tests.fixtures.synthetic()
        metadata = prior_tests.fixtures.production_metadata()
        context = dict(metadata.pop('context'), **result['context'])
        result.update(metadata, context=context)
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            code = prior_tests.source_file(directory, 'helper.py', 'VALUE = 1\n')
            proposal_pin = mode.producer.profiles.write_new(directory/'proposal.json',
                dict(witness={key: context['witness'][key] for key in ('z', 'alpha', 'markers')},
                     source_code=code, replay_source_sha256=dict(provider={code['path']: code['sha256']})))
            context['witness_source'] = proposal_pin
            result['context_sha256'] = digest = sha256(mode.canonical(context)).hexdigest()
            for row in result['nodes']:
                row['receipt']['context_sha256'] = digest
            result['checkpoints'] = len(result['nodes'])+1
            mode.producer.profiles.write_new(directory/'context.json', dict(context=context,
                context_sha256=digest, directed=False, complete_requested_band=False, whole_code_certificate=False))
            for sequence in range(result['checkpoints']):
                changed = result['nodes'] if not sequence else [result['nodes'][sequence-1]]
                mode.producer.profiles.write_new(directory/f'node_{sequence:06d}.json',
                    dict(sequence=sequence, context_sha256=digest, changed_nodes=changed,
                         complete_requested_band=False, covers_all_supports_in_band=False, directed=False,
                         batch_source_verified_at_exit=False, whole_code_certificate=False, log2_upper_rational=None))
            driver_pin = mode.producer.profiles.write_new(directory/'driver.json',
                dict(driver=code, witness_source=proposal_pin, requested_weight_range=result['requested_weight_range'],
                     cell_failure_bits=result['cell_failure_bits'], aggregate_failure_bits=result['requested_failure_bits'],
                     numerical_provider=Path(mode.producer.__file__).name))
            mode.producer.profiles.write_new(directory/'summary.json', result)
            numeric = mode.previous.audit_numeric
            # Test-only small-domain injection; production calls the frozen P4096 numerical validator.
            with patch.object(mode.previous, 'audit_numeric', lambda record, pairs, budget: numeric(record, 4, budget)):
                report = mode.audit_summary(directory/'summary.json')
            self.assertEqual(report['schema'], 'gamma-midpoint-weight-band-typed-independent-audit-v2')
            self.assertEqual(report['execution_driver_receipt'], driver_pin)
            self.assertEqual(report['checked_non_json_file_pins'], [dict(code, kind='python-source')])
            self.assertTrue(all(Path(x['path']).suffix == '.json' for x in report['checked_file_pins']))
            self.assertEqual(report['saved_band_log2_upper_rational'], report['reaggregated_cell_log2_upper_rational'])
            self.assertTrue(report['nested_source_maps_checked'])
            self.assertEqual(report['frozen_typed_v1_auditor_sha256'], mode.TYPED_V1_SHA256)
            self.assertFalse(report['coefficient_replay_in_this_audit'])
            self.assertFalse(report['cell_transfer_recomputed'])
            self.assertFalse(report['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
