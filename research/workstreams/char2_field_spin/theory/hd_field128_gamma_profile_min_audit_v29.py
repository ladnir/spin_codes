"""Literal-only saved-data audit of the distinct profile-min producer.

Recheck fresh execution/source links and exact literals, not calibration
ancestry. Rows already include each chosen witness's outer prefactor. This
audit checks counts/partition and their 384-bit positive union, not coefficients,
regional matrices, per-profile transfers or the producer's winner selection.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx
import hd_field128_gamma_profile_min_certificate_v29 as producer

HERE, NAME = Path(__file__).resolve().parent, Path(__file__).name
PRODUCER_SHA = 'b79bd9edbeaab13f6adebb39cc2d73c16a784d453b54e96c73208c3f7d45d5e3'
SCHEMA = 'gamma-profile-min-literal-only-saved-audit-v29'


def require(test, message):
    if not test:
        raise ValueError(message)


def same(actual, expected, label):
    require(isinstance(actual, dict), label+' must be an object')
    for key, value in expected.items():
        require(actual.get(key) == value and type(actual.get(key)) is type(value), label+': '+key)


def source_pins():
    require(sha256((HERE/producer.NAME).read_bytes()).hexdigest() == PRODUCER_SHA, 'Wrong producer bytes')
    cert = producer.helpers()
    pins = producer.source_pins(cert)
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return cert, dict(sorted(pins.items()))


def check_sources(pins):
    for name, digest in pins.items():
        require(sha256((HERE/name).read_bytes()).hexdigest() == digest, 'Source changed: '+name)


class Inputs:
    def __init__(self):
        self.pins = {}

    def load(self, name, expected=None):
        name = producer.literals.proposal_name(name)
        path = HERE/name
        require(path.resolve().is_relative_to(HERE.resolve()), 'Escaping input path')
        require(path.stat().st_size <= 1 << 20, 'Input exceeds 1 MiB')
        raw = path.read_bytes()
        digest = sha256(raw).hexdigest()
        require(expected is None or expected == digest, 'Input byte pin mismatch: '+name)
        require(name not in self.pins or self.pins[name] == digest, 'Input changed during audit')
        self.pins[name] = digest
        obj = json.loads(raw)
        require(isinstance(obj, dict), 'Expected a JSON object')
        return obj

    def pinned(self, pin):
        require(isinstance(pin, dict) and set(pin) == {'path', 'sha256'}, 'Malformed file pin')
        require(isinstance(pin['sha256'], str) and len(pin['sha256']) == 64
                and all(c in '0123456789abcdef' for c in pin['sha256']), 'Malformed SHA-256')
        return self.load(pin['path'], pin['sha256'])

    def check(self):
        for name, digest in self.pins.items():
            require(sha256((HERE/name).read_bytes()).hexdigest() == digest, 'Input changed at exit')


def bound(value, cert):
    """Validate stored upward logarithmic enclosures and explicit zero semantics."""
    require(isinstance(value, dict) and type(value.get('zero')) is bool, 'Missing zero flag')
    if value['zero']:
        require(value.get('log2_upper_rational') is None and value.get('log2_upper_enclosure') is None,
                'Zero bound has a logarithm')
        return None
    text = value['log2_upper_rational']
    require(isinstance(text, str) and str(F(text)) == text, 'Noncanonical ceiling')
    enclosure = arb(value['log2_upper_enclosure'])
    require(enclosure.is_finite() and cert.base.rational(F(text)) >= enclosure.upper(),
            'Saved ceiling understates its enclosure')
    return F(text)


def verify_bank(bank, proposal, r, cert):
    witness = proposal['witness']
    expected = dict(construction=producer.literals.CONSTRUCTION, field_order=1 << 128, field_bits=128,
        state_dimension=30, inner_width=64, target_distance='2/5', output_length=1 << 21,
        source_dimension=1 << 20, outer_groups=8192, forbidden_weight=838860,
        outer_measure='invertible-toeplitz', projective_divisor=str((1 << 128)-1),
        pair_swaps=True, state_persists_across_regions=True, full_field_random_coefficients=True,
        initial_state='zero', final_state='unrestricted', gamma_extension_degree=16,
        gamma_base_field_order=256, pair_slots_per_region=4096, regions=256,
        packets_per_region=128, composed_pairs_per_region=64,
        all_pair_before_alpha=True, counts_outside_alpha=True, fresh_bank_rebuilt=True,
        no_coefficient_truncation=True, original_weight=r, directed=True, whole_code_certificate=False,
        historical_origin_verified=False, old_numerical_bounds_reused=False,
        proposal_source=proposal['source'], exact_input_witness=witness,
        source_marker_key=proposal['source_marker_key'],
        normalized_regional_full_range=[0, r//2], normalized_regional_matrix_count=r//2+1,
        outer_prefactor_formula='exp2(bank.constant + bank.slope * original_weight)',
        source_sha256=cert.provider.source_pins(), bank_source_sha256=cert.provider.source_pins())
    same(bank, expected, 'fresh bank')
    same(bank['witness'], dict(witness, precision=256, packet_mantissa_bits=80,
                              all_pair_before_alpha=True), 'bank literal witness')
    with ctx.workprec(256):
        z, alpha = map(F, (witness['z'], witness['alpha']))
        q = 1 << 128
        density = cert.tail.log2arb(cert.base.rational(cert.base.outer_density(q, 'invertible-toeplitz')))
        slope = cert.base.rational(alpha)*density
        constant = -cert.base.rational(alpha)*(cert.tail.log2arb(arb(q-1))
            + 838860*cert.tail.log2arb(cert.base.rational(z)))
        expected_prefix = cert.frozen.frozen.support_log(cert.base.upper(arb(2)**(constant+slope*r)))
    require(bank['outer_prefactor'] == expected_prefix, 'Own original-r prefix differs from exact literals')


def audit_rows(record, context, cert, budget, *, pairs=4096):
    r = record['active_groups']
    full = cert.tail.full_counts(pairs, r)
    counts = cert.get_counts(pairs, r, None, 64, budget)
    same(context, dict(expected_even_rows=len(full), expected_profile_pairs=len(full)**2,
                      count_provenance=counts.diagnostics()), 'profile domain/count policy')
    rows = record['rows']
    require(len(rows) == len(full), 'Missing complete rows')
    positive, winners = 0, [0]*len(record['proposals'])
    with ctx.workprec(384):
        for expected, row in zip(full, rows):
            budget.check()
            actual_positive = sum(counts.count(expected, other) > 0 for other in full)
            same(row, dict(even_full=expected, odd_full_range=[full.start, full.stop-1],
                profile_pairs_evaluated=len(full), profiles_with_positive_count=actual_positive,
                outer_prefactors_already_included=True, complete_even_profile_row=True,
                covers_all_odd_profiles_for_this_even_profile=True, complete_selected_weight=False,
                whole_code_certificate=False), 'complete weighted row')
            row_wins = row['witness_winning_profile_counts']
            require(isinstance(row_wins, list) and len(row_wins) == len(winners)
                and all(type(n) is int and n >= 0 for n in row_wins)
                and sum(row_wins) == actual_positive, 'Winner counts differ from the common count domain')
            bound(row['weighted_row_bound'], cert)
            winners = [a+b for a, b in zip(winners, row_wins)]
            positive += actual_positive
        same(record, dict(completed_even_rows=len(full), pending_even_full_range=None,
            profile_pairs_evaluated=len(full)**2, profiles_with_positive_count=positive,
            witness_winning_profile_counts=winners), 'complete selected profile domain')
        result = record['result']
        fresh = producer.aggregate_rows(rows, full, cert)
        same(result, dict(zero=fresh['zero'], aggregation_precision=384), 'union precision/zero')
        raw, saved, rebuilt = (bound(obj, cert) for obj in
                              (result['direct_union_256'], result['saved_row_union_384'], fresh))
        require((raw is None) == (saved is None) == (rebuilt is None), 'Union zero mismatch')
        if rebuilt is None:
            require(record['log2_upper_rational'] is None and result['log2_upper_rational'] is None,
                    'Zero selected union has a ceiling')
            admitted = True
        else:
            ceiling = F(record['log2_upper_rational'])
            require(str(ceiling) == record['log2_upper_rational'] == result['log2_upper_rational']
                    and ceiling >= max(raw, saved, rebuilt), 'Published union understates a row union')
            admitted = ceiling <= -record['requested_failure_bits']
        require(result['meets_requested_selected_weight_target'] is admitted and admitted,
                'Selected-weight requested target not met')
    return dict(even_rows=len(full), profile_pairs=len(full)**2, positive_count_profiles=positive,
                witness_winning_profile_counts=winners, fresh_saved_row_union_384=fresh)


def audit_summary(name, expected_sha256, *, max_seconds=120.):
    cert, pins = source_pins()
    inputs, budget = Inputs(), cert.base.Budget(max_seconds)
    producer_pins = producer.source_pins(cert)
    record = inputs.pinned(dict(path=name, sha256=expected_sha256))
    directory = Path(name).parent.as_posix()
    try:
        check_sources(pins)
        r = producer.literals.validate_target(record['active_groups'])
        common = dict(construction=producer.literals.CONSTRUCTION, field_order=1 << 128,
            state_dimension=30, target_distance='2/5', active_groups=r, source_sha256=producer_pins,
            directed=True, whole_code_certificate=False, historical_origin_verified=False,
            old_numerical_bounds_reused=False)
        same(record, dict(common, schema=producer.SCHEMA, status='complete', failure=None, pending_stage=None,
            complete_selected_weight=True, covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
            outer_prefactors_already_included=True), 'selected summary')
        target = record['requested_failure_bits']
        require(type(target) is int and 40 <= target <= 4096, 'Expected at least 40-bit admission')
        proposals = record['proposals']
        require(isinstance(proposals, list) and 2 <= len(proposals) <= 4, 'Wrong witness count')
        seen = set()
        for proposal in proposals:
            original = inputs.pinned(proposal['source'])
            require(proposal['source']['path'] not in seen, 'Duplicate proposal input')
            seen.add(proposal['source']['path'])
            witness, key = producer.literals.project_literals(original, r)
            same(proposal, dict(witness=witness, source_marker_key=key), 'exact literal-only projection')
        driver = inputs.load(directory+'/driver.json')
        same(driver, dict(schema='gamma-profile-min-directed-driver-v29', source_sha256=producer_pins,
            proposals=proposals, active_groups=r, hard_wall_seconds=600, cooperative_seconds=540.,
            memory_bytes=4 << 30, table_bytes=64 << 20, minimum_free_bytes=20 << 30,
            maximum_output_bytes=1 << 20, maximum_simultaneous_banks=1, fresh_banks_required=True,
            directed=False, whole_code_certificate=False), 'execution driver')
        packed = inputs.load(directory+'/context.json')
        context = packed['context']
        digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        require(digest == packed['context_sha256'] == record['context_sha256'], 'Context digest mismatch')
        same(packed, dict(complete_selected_weight=False, whole_code_certificate=False), 'context scope')
        same(context, dict(common, schema=producer.SCHEMA, banks=record['banks'], precision=256,
            packet_rounding_bits=80, grid=128, repetitions=256, paired_region_power=128,
            common_count_cap_outside_alpha=True, own_outer_prefactor_inside_each_completed_profile_term=True,
            entrywise_cross_witness_minimum=False, original_weight_used_for_all_histograms_and_prefixes=True), 'context')
        require(len(record['banks']) == len(proposals), 'Not every witness has a fresh bank')
        for i, (entry, proposal) in enumerate(zip(record['banks'], proposals)):
            budget.check()
            same(entry, dict(index=i, proposal=proposal), 'bank identity')
            require(entry['metadata']['path'] == directory+f'/bank_{i:02d}.json', 'Wrong bank file')
            verify_bank(inputs.pinned(entry['metadata']), proposal, r, cert)
        rows, checkpoints = record['rows'], record['row_checkpoints']
        require(len(checkpoints) == len(rows), 'Missing row checkpoints')
        for i, (row, pin) in enumerate(zip(rows, checkpoints)):
            budget.check()
            require(pin['path'] == directory+f'/row_{i:06d}.json', 'Wrong checkpoint path/order')
            checkpoint = inputs.pinned(pin)
            expected = dict(row, context_sha256=digest, active_groups=r, directed=True,
                scope='One fully weighted row only; not a selected-weight bound')
            require(checkpoint == expected, 'Checkpoint differs from completed row')
        numeric = audit_rows(record, context, cert, budget)
        budget.check()
        return dict(schema=SCHEMA, audit_result='PASS', input_summary=dict(path=name, sha256=expected_sha256),
            context_sha256=digest, construction=producer.literals.CONSTRUCTION, active_groups=r,
            state_dimension=30, field_order=1 << 128, target_distance='2/5', producer_sha256=PRODUCER_SHA,
            log2_upper_rational=record['log2_upper_rational'], **numeric,
            requested_failure_bits=target, meets_requested_selected_weight_target=True,
            complete_selected_weight=True, covers_all_matching_profiles=True,
            exact_literal_bank_context_matches=True, original_r_prefixes_recomputed=True,
            historical_origin_verified=False, coefficient_replay_in_this_audit=False,
            normalized_regional_replay=False, transfer_or_winner_replay=False,
            row_partition_and_count_policy_checked=True, directed_union_precision=384,
            second_outer_prefix_applied=False, whole_code_certificate=False,
            checked_file_pins=[dict(path=p, sha256=h) for p, h in sorted(inputs.pins.items())],
            source_sha256=pins,
            scope='Saved-data and producer-source audit; no coefficient, region, transfer or historical calibration replay')
    finally:
        check_sources(pins)
        inputs.check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--output-name', default='independent_profile_min_audit.json')
    args = parser.parse_args()
    with producer.portable.remote_limits():
        directory = HERE/Path(producer.literals.proposal_name(args.summary)).parent
        producer.portable.reserve_guard(directory, 1 << 20)
        report = audit_summary(args.summary, args.expected_sha256)
        outputs = producer.portable.Outputs(directory, lambda obj, unused: obj)
        pin = outputs.write(args.output_name, report)
        outputs.check()
    print(json.dumps(dict(audit_result='PASS', audit=pin, active_groups=report['active_groups'],
                         log2_upper_rational=report['log2_upper_rational'])), flush=True)


if __name__ == '__main__':
    main()
