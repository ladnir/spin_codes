"""Tiny independent count, tree/union, metadata, and byte-link audit tests."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb

import hd_field128_gamma_midpoint_weight_band_audit as mode
import hd_field128_gamma_midpoint_weight_bands_tests as fixtures
from hd_field128_gamma_exact_profile_counts import joint_count


def synthetic():
    table = fixtures.tables(4, 2, 6)
    evaluator = fixtures.fake_evaluator(table, split_to_points=True)
    result = mode.producer.cover(evaluator, max_nodes=1000, cell_failure_bits=10)
    context = result['context']
    context.update(root_cells=result['root_cells'],
        parity_weights={str(eps): [r for r in table.weights if r % 2 == eps] for eps in (0, 1)},
        regional_histograms_prepared=table.regional_histograms)
    result['context_sha256'] = digest = sha256(mode.canonical(context)).hexdigest()
    # Forced split caps are synthetic, so use a valid enclosure of their zero log.
    for row in result['nodes']:
        row['receipt']['context_sha256'] = digest
    return result


def production_metadata():
    context = dict(mode.EXPECTED)
    context.update(source_sha256=mode.producer.source_pins(),
        bank_source_sha256=mode.producer.provider.source_pins(),
        witness=dict(z='1/2', alpha='1/4', markers=[['1/2', '1/3']], precision=256,
                     grid=128, packet_mantissa_bits=80, state_dimension=30, inner_width=64,
                     gamma_extension_degree=16, gamma_base_field_order=256,
                     all_pair_before_alpha=True, joint_rank_return_cap=True))
    return dict(context=context, schema='gamma-midpoint-weight-band-cover-v1', directed=True,
                numerical_backend='actual-all-pair-bank', whole_code_certificate=False,
                batch_source_verified_at_entry=True, batch_source_verified_at_exit=True)


class AuditTests(unittest.TestCase):
    def test_independent_prefix_counts_every_tiny_rectangle(self):
        for pairs in range(1, 7):
            for r in range(2*pairs+1):
                own = mode.RectangleCounts(pairs, r, mode.base.Budget())
                scalar = mode.producer.boxes.midpoint.MidpointCounts(pairs, r, max_singletons=0)
                last = r//2
                for a in range(last+1):
                    for b in range(last+1):
                        self.assertEqual(own.rectangle(a, a, b, b), scalar.midpoint_upper(a, b))
                        self.assertGreaterEqual(own.rectangle(a, a, b, b), joint_count(pairs, r, a, b))
                for al in range(-1, last+1):
                    for ah in range(al, last+2):
                        for bl, bh in ((-1, last+1), (0, 0), (last//2, last)):
                            expected = sum(scalar.midpoint_upper(a, b) for a in range(al, ah+1)
                                           for b in range(bl, bh+1))
                            self.assertEqual(own.rectangle(al, ah, bl, bh), expected)

    def test_complete_tree_counts_and_fresh_union(self):
        result = synthetic()
        audited = mode.audit_numeric(result, 4, mode.base.Budget())
        self.assertEqual(audited['weights'], [2, 3, 4, 5, 6])
        self.assertEqual(audited['cells_checked'], len(result['nodes']))
        self.assertEqual(audited['union'], result['log2_upper_rational'])
        self.assertEqual(audited['positive_count_zero_cells'], 0)

    def test_corrupt_counts_partition_weights_and_enclosures_rejected(self):
        result = synthetic()
        for mutation in ('count', 'split', 'weight', 'enclosure', 'union', 'partial', 'duplicate'):
            bad = deepcopy(result)
            if mutation == 'count':
                bad['nodes'][0]['receipt']['count_upper'] += 1
            elif mutation == 'split':
                next(row for row in bad['nodes'] if row['status'] == 'split')['split']['integer_cut'] += 1
            elif mutation == 'weight':
                bad['requested_weights'].pop()
            elif mutation == 'enclosure':
                next(row for row in bad['nodes'] if not row['receipt']['zero'])['receipt']['log2_upper_enclosure'] = '1000000'
            elif mutation == 'union':
                bad['log2_upper_rational'] = '-1000000'
            elif mutation == 'partial':
                bad['nodes'][-1]['status'] = 'pending'
            else:
                bad['nodes'].append(deepcopy(bad['nodes'][0]))
            with self.assertRaises(ValueError, msg=mutation):
                mode.audit_numeric(bad, 4, mode.base.Budget())

    def test_production_scope_and_source_closure(self):
        record = production_metadata()
        mode._metadata(record)
        for key, value in (('state_dimension', 24), ('target_distance', '9/20'),
                           ('normalized_every_original_weight', False), ('bank_source_sha256', {})):
            bad = deepcopy(record)
            bad['context'][key] = value
            with self.assertRaises(ValueError):
                mode._metadata(bad)
        pins = mode.source_pins()
        mode.check_source_pins(pins)
        pins[Path(mode.__file__).name] = '0'*64
        with self.assertRaises(ValueError):
            mode.check_source_pins(pins)

    def test_byte_link_and_context_fail_closed_without_arithmetic(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'input.json'
            pin = mode.producer.profiles.write_new(path, dict(value=1))
            inputs = mode.Inputs()
            self.assertEqual(inputs.load(path, pin['sha256'])[0], dict(value=1))
            with self.assertRaises(ValueError):
                inputs.load(path, '0'*64)
            original = Path.read_bytes
            with patch.object(Path, 'read_bytes', lambda p: original(p)+b' ' if p == path else original(p)):
                with self.assertRaises(ValueError):
                    inputs.check()
        result = synthetic()
        result['context']['requested_weights'] = [2, 6]
        with self.assertRaisesRegex(ValueError, 'Context hash'):
            mode.audit_numeric(result, 4, mode.base.Budget())

    def test_timeout_and_zero_count_edges(self):
        with self.assertRaises(TimeoutError):
            mode.audit_numeric(synthetic(), 4, mode.base.Budget(1e-12))
        for pairs in (1, 4):
            self.assertEqual(mode.RectangleCounts(pairs, pairs, mode.base.Budget()).rectangle(0, 0, 0, 0), 2)
            self.assertEqual(mode.RectangleCounts(pairs, 0, mode.base.Budget()).rectangle(0, 0, 0, 0), 1)
        result = synthetic()
        zero = next(row for row in result['nodes'] if row['receipt']['zero'])
        self.assertEqual(zero['receipt']['count_upper'], 0)
        zero['receipt']['log2_upper_rational'] = '0'
        with self.assertRaises(ValueError):
            mode.audit_numeric(result, 4, mode.base.Budget())

    def test_public_file_resolver_with_explicit_tiny_numeric_injection(self):
        result = synthetic()
        metadata = production_metadata()
        context = dict(metadata.pop('context'), **result['context'])
        result.update(metadata, context=context)
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            proposal_pin = mode.producer.profiles.write_new(directory/'proposal.json',
                dict(witness={key: context['witness'][key] for key in ('z', 'alpha', 'markers')}))
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
            mode.producer.profiles.write_new(directory/'summary.json', result)
            numeric = mode.audit_numeric
            # This injection is test-only: public production audit always uses P4096.
            with patch.object(mode, 'audit_numeric', lambda record, pairs, budget: numeric(record, 4, budget)):
                report = mode.audit_summary(directory/'summary.json')
            self.assertEqual(report['audit_result'], 'PASS')
            self.assertEqual(report['saved_band_log2_upper_rational'], report['reaggregated_cell_log2_upper_rational'])
            self.assertFalse(report['coefficient_replay_in_this_audit'])
            self.assertFalse(report['cell_transfer_recomputed'])
            self.assertEqual(len(report['checked_file_pins']), result['checkpoints']+3)


if __name__ == '__main__':
    unittest.main()
