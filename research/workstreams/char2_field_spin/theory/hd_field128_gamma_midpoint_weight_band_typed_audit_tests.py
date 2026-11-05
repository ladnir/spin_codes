"""Tiny typed-provenance and unchanged numerical-contract tests."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import hd_field128_gamma_midpoint_weight_band_typed_audit as mode
import hd_field128_gamma_midpoint_weight_band_audit_tests as fixtures


def source_file(directory, name, text):
    path = directory/name
    path.write_text(text, encoding='utf-8')
    return dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest())


class TypedAuditTests(unittest.TestCase):
    def test_typed_recursive_links_and_no_source_execution(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            code = source_file(directory, 'helper.py', "raise RuntimeError('must not execute')\n")
            note = source_file(directory, 'proof.md', '# A mathematical note\n')
            linked = mode.producer.profiles.write_new(directory/'linked.json', dict(code=code, note=note))
            inputs, sources = mode.Inputs(), {}
            mode._provenance(dict(first=linked, again=linked), inputs, sources)
            self.assertEqual(len(inputs.pins), 1)
            self.assertEqual(len(inputs.non_json_pins), 2)
            self.assertEqual(set(sources.values()), {code['sha256'], note['sha256']})
            inputs.check()
            code_path = Path(code['path'])
            original = Path.read_bytes
            with patch.object(Path, 'read_bytes', lambda p: original(p)+b' ' if p == code_path else original(p)):
                with self.assertRaises(ValueError):
                    inputs.check()

    def test_extension_content_mismatch_and_tamper_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            good = source_file(directory, 'source.py', 'VALUE = 1\n')
            inputs = mode.Inputs()
            with self.assertRaises(ValueError):
                inputs.load_source(good['path'], '0'*64)
            for name, text in (('disguised.py', '{"value": 1}'), ('disguised.md', '[1,2]'),
                               ('other.txt', 'text'), ('invalid.py', 'not python ???')):
                pin = source_file(directory, name, text)
                with self.assertRaises((ValueError, SyntaxError)):
                    mode._provenance(pin, mode.Inputs(), {})
            bad = source_file(directory, 'source.json', 'VALUE = 1\n')
            with self.assertRaises(ValueError):
                mode.Inputs().load(bad['path'])
            with self.assertRaises(ValueError):
                mode.Inputs().load(good['path'])
            scalar = source_file(directory, 'scalar.json', '42')
            with self.assertRaises(ValueError):
                mode.Inputs().load(scalar['path'])

    def test_conflicting_links_and_orchestration_source_maps(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            pin = source_file(directory, 'helper.py', 'VALUE = 1\n')
            inputs, sources = mode.Inputs(), {}
            mode._provenance(pin, inputs, sources)
            with self.assertRaises(ValueError):
                mode._provenance(dict(pin, sha256='0'*64), inputs, sources)
            mode._provenance(dict(orchestration_source_sha256={pin['path']: pin['sha256']}), inputs, sources)
            with self.assertRaises(ValueError):
                mode._provenance(dict(orchestration_source_sha256={pin['path']: '0'*64}), inputs, sources)

    def test_frozen_closure_and_numeric_oracle_unchanged(self):
        pins = mode.source_pins()
        mode.check_source_pins(pins)
        self.assertEqual(pins[Path(mode.previous.__file__).name], mode.PREVIOUS_SHA256)
        self.assertIn(Path(mode.__file__).name, pins)
        result = fixtures.synthetic()
        report = mode.previous.audit_numeric(result, 4, mode.base.Budget())
        self.assertEqual(report['union'], result['log2_upper_rational'])
        bad = deepcopy(result)
        bad['nodes'][0]['receipt']['count_upper'] += 1
        with self.assertRaises(ValueError):
            mode.previous.audit_numeric(bad, 4, mode.base.Budget())
        with self.assertRaises(TimeoutError):
            mode.previous.audit_numeric(result, 4, mode.base.Budget(1e-12))

    def test_public_resolver_source_link_and_execution_receipt(self):
        result = fixtures.synthetic()
        metadata = fixtures.production_metadata()
        context = dict(metadata.pop('context'), **result['context'])
        result.update(metadata, context=context)
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            code = source_file(directory, 'helper.py', 'VALUE = 1\n')
            proposal_pin = mode.producer.profiles.write_new(directory/'proposal.json',
                dict(witness={key: context['witness'][key] for key in ('z', 'alpha', 'markers')},
                     source_code=code))
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
            self.assertEqual(report['schema'], 'gamma-midpoint-weight-band-typed-independent-audit-v1')
            self.assertEqual(report['execution_driver_receipt'], driver_pin)
            self.assertEqual(report['checked_non_json_file_pins'], [dict(code, kind='python-source')])
            self.assertTrue(all(Path(x['path']).suffix == '.json' for x in report['checked_file_pins']))
            self.assertEqual(report['saved_band_log2_upper_rational'], report['reaggregated_cell_log2_upper_rational'])
            self.assertFalse(report['coefficient_replay_in_this_audit'])
            self.assertFalse(report['cell_transfer_recomputed'])
            self.assertFalse(report['whole_code_certificate'])

    def test_no_frozen_source_mutation(self):
        self.assertEqual(sha256(Path(mode.previous.__file__).read_bytes()).hexdigest(), mode.PREVIOUS_SHA256)


if __name__ == '__main__':
    unittest.main()

