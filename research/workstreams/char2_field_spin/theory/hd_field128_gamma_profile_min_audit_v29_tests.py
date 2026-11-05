"""Memory-only audit contracts and tiny Arb rows; no bank construction or files."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
import unittest
from unittest.mock import patch

from flint import arb, ctx
import hd_field128_gamma_profile_min_audit_v29 as audit
import hd_field128_gamma_midpoint_certificate as cert


def encoded(obj):
    return (json.dumps(obj, sort_keys=True, allow_nan=False)+'\n').encode()


class MemoryInputs:
    def __init__(self, data):
        self.data, self.pins = data, {}
    def load(self, name, expected=None):
        raw = encoded(self.data[name])
        digest = sha256(raw).hexdigest()
        audit.require(expected is None or digest == expected, 'memory pin mismatch')
        self.pins[name] = digest
        return deepcopy(self.data[name])
    pinned = audit.Inputs.pinned
    def check(self):
        for name, digest in self.pins.items():
            audit.require(sha256(encoded(self.data[name])).hexdigest() == digest, 'memory input changed')


def fixture():
    p, r, data = audit.producer, 2, {}
    def put(name, obj):
        data[name] = deepcopy(obj)
        return dict(path=name, sha256=sha256(encoded(obj)).hexdigest())
    pins = p.source_pins(cert)
    witness = dict(z='1/2', alpha='1/2', markers=[['1', '1']], grid=128, state_dimension=30)
    proposals, banks = [], []
    with ctx.workprec(256):
        # __new__ + metadata only: never call a numerical bank constructor.
        bank = cert.provider.AliasBank.__new__(cert.provider.AliasBank)
        bank.z = bank.alpha = F(1, 2)
        bank.markers, bank.state, bank.q, bank.field_bits = ((F(1), F(1)),), 30, 1 << 128, 128
        bank.distance, bank.forbidden, bank.outer_measure = F(2, 5), 838860, 'invertible-toeplitz'
        bank.precision, bank.bits, bank.grid, bank.diagnostics = 256, 80, 128, []
        bank.construction, bank.pins = p.literals.CONSTRUCTION, cert.provider.source_pins()
        density = cert.tail.log2arb(cert.base.rational(cert.base.outer_density(bank.q, bank.outer_measure)))
        exponent = -cert.base.rational(bank.alpha)*(cert.tail.log2arb(arb(bank.q-1))
            + bank.forbidden*cert.tail.log2arb(cert.base.rational(bank.z))) + cert.base.rational(bank.alpha)*density*r
        prefactor = cert.base.upper(arb(2)**exponent)
        for i in range(2):
            original = dict(construction=bank.construction, active_groups=r, state_dimension=30,
                target_distance='2/5', witness=witness,
                source={'path': 'deliberately-unread-historical.json', 'sha256': '0'*64})
            proposal = dict(source=put(f'p{i}.json', original), witness=deepcopy(witness), source_marker_key='markers')
            proposals.append(proposal)
            metadata = dict(bank.metadata(), proposal_source=proposal['source'], exact_input_witness=witness,
                source_marker_key='markers', original_weight=r, fresh_bank_rebuilt=True,
                normalized_regional_full_range=[0, 1], normalized_regional_matrix_count=2,
                outer_prefactor=cert.frozen.frozen.support_log(prefactor),
                outer_prefactor_formula='exp2(bank.constant + bank.slope * original_weight)',
                historical_origin_verified=False, old_numerical_bounds_reused=False)
            banks.append(dict(index=i, proposal=proposal, metadata=put(f'result/bank_{i:02d}.json', metadata)))
        counts = cert.get_counts(4096, r, None, 64, cert.base.Budget())
        matrix = {f: (arb(2)**-4000, arb(0), arb(0), arb(2)**-4000) for f in counts.full}
        value, partial = p.directed_rows([matrix, matrix], [prefactor, prefactor], 4096, r, counts, cert)
        raw = cert.frozen.frozen.support_log(value)
    fresh = p.aggregate_rows(partial['rows'], counts.full, cert)
    ceiling = str(max(F(raw['log2_upper_rational']), F(fresh['log2_upper_rational'])))
    common = dict(construction=bank.construction, field_order=1 << 128, state_dimension=30,
        target_distance='2/5', active_groups=r, source_sha256=pins, directed=True,
        whole_code_certificate=False, historical_origin_verified=False, old_numerical_bounds_reused=False)
    context = dict(common, schema=p.SCHEMA, banks=banks, precision=256, packet_rounding_bits=80, grid=128,
        repetitions=256, paired_region_power=128, expected_even_rows=2, expected_profile_pairs=4,
        count_provenance=counts.diagnostics(), common_count_cap_outside_alpha=True,
        own_outer_prefactor_inside_each_completed_profile_term=True, entrywise_cross_witness_minimum=False,
        original_weight_used_for_all_histograms_and_prefixes=True)
    digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    put('result/context.json', dict(context=context, context_sha256=digest,
        complete_selected_weight=False, whole_code_certificate=False))
    put('result/driver.json', dict(schema='gamma-profile-min-directed-driver-v29', source_sha256=pins,
        proposals=proposals, active_groups=r, hard_wall_seconds=600, cooperative_seconds=540.,
        memory_bytes=4 << 30, table_bytes=64 << 20, minimum_free_bytes=20 << 30,
        maximum_output_bytes=1 << 20, maximum_simultaneous_banks=1, fresh_banks_required=True,
        directed=False, whole_code_certificate=False))
    checkpoints = [put(f'result/row_{i:06d}.json', dict(row, context_sha256=digest, active_groups=r,
        directed=True, scope='One fully weighted row only; not a selected-weight bound'))
        for i, row in enumerate(partial['rows'])]
    record = dict(common, schema=p.SCHEMA, status='complete', failure=None, pending_stage=None,
        complete_selected_weight=True, covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
        outer_prefactors_already_included=True, requested_failure_bits=48, proposals=proposals,
        banks=banks, context_sha256=digest, row_checkpoints=checkpoints, **partial,
        log2_upper_rational=ceiling, result=dict(zero=False, log2_upper_rational=ceiling,
            direct_union_256=raw, saved_row_union_384=fresh, aggregation_precision=384,
            meets_requested_selected_weight_target=True))
    return data, put('result/summary.json', record)


class AuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, cls.pin = fixture()

    def test_complete_literal_only_memory_lifecycle(self):
        loader = MemoryInputs(deepcopy(self.data))
        with patch.object(audit, 'Inputs', return_value=loader):
            report = audit.audit_summary(self.pin['path'], self.pin['sha256'])
        self.assertEqual(report['audit_result'], 'PASS')
        self.assertEqual(report['profile_pairs'], 4)
        self.assertEqual(report['requested_failure_bits'], 48)
        self.assertTrue(report['meets_requested_selected_weight_target'])
        self.assertFalse(report['historical_origin_verified'])
        self.assertFalse(report['second_outer_prefix_applied'])
        self.assertFalse(report['transfer_or_winner_replay'])
        self.assertNotIn('deliberately-unread-historical.json', loader.pins)

    def test_literal_prefix_and_bank_scope_corruption(self):
        original = self.data['result/bank_00.json']
        proposal = self.data['result/summary.json']['proposals'][0]
        for key, value in (('original_weight', 3), ('fresh_bank_rebuilt', False),
                           ('field_order', 256), ('all_pair_before_alpha', False)):
            wrong = deepcopy(original); wrong[key] = value
            with self.assertRaises(ValueError):
                audit.verify_bank(wrong, proposal, 2, cert)
        wrong = deepcopy(original); wrong['witness']['alpha'] = '1/3'
        with self.assertRaises(ValueError):
            audit.verify_bank(wrong, proposal, 2, cert)
        wrong = deepcopy(original); wrong['outer_prefactor']['log2_upper_rational'] = '0'
        with self.assertRaises(ValueError):
            audit.verify_bank(wrong, proposal, 2, cert)

    def test_row_partition_counts_target_and_union_corruptions(self):
        context = self.data['result/context.json']['context']
        original = self.data['result/summary.json']
        for mutate in (lambda d: d['rows'].reverse(),
                       lambda d: d['rows'][0].update(profiles_with_positive_count=999),
                       lambda d: d['rows'][0].update(outer_prefactors_already_included=False),
                       lambda d: d['result'].update(meets_requested_selected_weight_target=False),
                       lambda d: d.update(log2_upper_rational=str(F(d['log2_upper_rational'])-1))):
            wrong = deepcopy(original); mutate(wrong)
            with self.assertRaises(ValueError):
                audit.audit_rows(wrong, context, cert, cert.base.Budget())

    def test_pin_checkpoint_context_and_partial_rejection(self):
        for name, mutate in (
            ('result/row_000000.json', lambda d: d.update(even_full=1)),
            ('result/context.json', lambda d: d['context'].update(repetitions=128)),
            ('result/summary.json', lambda d: d.update(complete_selected_weight=False))):
            data = deepcopy(self.data); mutate(data[name])
            pin = dict(self.pin)
            if name == self.pin['path']:
                pin['sha256'] = sha256(encoded(data[name])).hexdigest()
            with patch.object(audit, 'Inputs', return_value=MemoryInputs(data)):
                with self.assertRaises(ValueError):
                    audit.audit_summary(pin['path'], pin['sha256'])


if __name__ == '__main__':
    unittest.main()
