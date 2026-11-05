"""Complete ideal-setup coverage from reviewed evidence."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import re
import shutil
from flint import arb, ctx
import hd_field128_gamma_selected_index_v31 as selected

HERE = Path(__file__).resolve().parent
CONSTRUCTION = 'toeplitz-fullfield-full-rs64-30-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'
INDEX = 'hd_field128_gamma30_d40_partial_cover_v31_20261010.json'
SELECTION = 'hd_field128_gamma30_d40_v31_selection_20261010.json'
PARENT = 'hd_field128_gamma30_d40_partial_cover_v30_20261010.json'
PARENT_CHECK = 'hd_field128_gamma30_d40_partial_cover_v30_20261010_independent_check.json'
TERMINAL = 'hd_field128_gamma30_d40_terminal384_20261010.json'
PINS = {
    INDEX: 'b98f340ba669ae9712be7ff9ac401af5bcccb7e83d93bbccf5025f37173664f0',
    SELECTION: '065dc8d8432cf70d1e4e02eb7cbadf3efc326155db78207f5192b598808bc6c9',
    PARENT: '875c572f4ae84b5bc139a4c21be9b0b4e18aed51ab1f288c45c125a260aab88f',
    PARENT_CHECK: '22357f681f345f94f13d72394df478984f8756d5b68fbaf0db917d0b6e9b604a',
    TERMINAL: '014aef17603f952a4cc5d4e8a93cc0281a7c80d0bff0bc889ec64757adde8c91',
    'hd_field128_gamma_selected_index_v31.py': 'a98322ce10ae59a98bd58a8ba6a634058402dd40111fab5d7ccf5bf901dc7328',
    'hd_field128_gamma_selected_index_v31_tests.py': '7efa1e843a513a7d916755341b9029989b36fd939a2ecdc7bf192f4766824604',
    'hd_field128_gamma30_d40_v31_record_20261010.py': '84ba1877bb100c210b0c90afd75e20c64dd17477d488e946d22eb54ff58a3091',
    'hd_field128_gamma_terminal_certificate.py': 'cf8f589bcc00a675c2a5754d89149623c801aee3fad0bb026df6d39c66c0bd80',
}
NOTES = {
    '../HD_FIELD128_LIFT_PROOF.md': '48b24b48a7ca7ea2d8f2bf90ba0bf4f798ff91f9708f446b54d9d1de072347d1',
    '../HD_FIELD128_GAMMA_RETURN.md': '0a09a2fce55aa3a6272656e37dce75ff584b798d535f7da560ac5cf8a5c2068c',
    '../HD_FIELD128_ALL_PAIR_COMPOSITION.md': '2bd44dc6f83603177cef082e62a3d60d0c9bd12f0f2bc5efd9552551c058a11a',
    '../HD_FIELD128_PROFILE_MIN.md': '82791dd0a61c0a2b3fd6aaac3d405423612bfb6223d5a7a2b446c5c6d3b48eff',
}
ADDITIONS = ((521, 547), (549, 575), (917, 917), (918, 918), (919, 919),
             (1146, 1146), (1147, 1147), (1148, 1148), (1149, 1149))
LIMIT, RESERVE = 32 << 10, 20 << 30


def require(test, message):
    if not test:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def relative(name):
    require(isinstance(name, str) and '\\' not in name, 'Expected relative POSIX path')
    p = PurePosixPath(name)
    require(p.parts and not p.is_absolute() and all(x not in ('.', '..') and ':' not in x for x in p.parts),
            'Unsafe relative path')
    path = (HERE / p).resolve()
    require(path != HERE and path.is_relative_to(HERE), 'Path escapes theory')
    return path


def source_pins():
    pins = selected.source_pins()
    for name, digest in {**PINS, **NOTES}.items():
        if name.endswith('.json'):
            continue
        require(sha256((HERE / name).read_bytes()).hexdigest() == digest, 'Frozen source changed: '+name)
        if name in pins:
            require(pins[name] == digest, 'Source conflict')
        pins[name] = digest
    for name in (Path(__file__).name, 'hd_field128_gamma_completion_v31_tests.py'):
        pins[name] = sha256((HERE / name).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def positive_union(entries):
    """One term per event, not per weight."""
    with ctx.workprec(384):
        total, ln2 = arb(0), arb(2).log()
        positive = 0
        for entry in entries:
            zero = entry.get('zero', False)
            require(type(zero) is bool, 'Invalid zero flag')
            if zero:
                require(entry['log2_upper_rational'] is None, 'Zero event has a logarithm')
                continue
            value = F(entry['log2_upper_rational'])
            total += (ln2 * arb(value.numerator) / value.denominator).exp()
            positive += 1
        require(positive > 0, 'Expected at least the positive terminal event')
        value = total.log() / ln2
        ceiling = F(int((value.upper() * 10**12).upper().ceil().fmpz()), 10**12)
        return str(ceiling), str(value)


def verify_objects(saved, parent, terminal, rebuilt):
    require(canonical(saved) == canonical(rebuilt), 'Helper reconstruction differs')
    common = dict(construction=CONSTRUCTION, field_order=1 << 128, source_dimension=1 << 20,
                  output_length=1 << 21, target_distance='2/5', whole_code_certificate=False)
    for record in (saved, parent, terminal):
        require(all(record.get(k) == v for k, v in common.items()), 'Ensemble/dimensions/distance mismatch')
    require(saved['schema'] == 'field128-gamma-partial-cover-index-v31-composite-selected', 'Wrong index schema')
    require(saved['whole_code_certificate'] is False and saved['completion_audit_performed'] is False,
            'Partial index was relabeled')
    require(len(parent['entries']) == 265 and len(saved['entries']) == 274, 'Wrong event count')
    require(canonical(saved['entries'][:265]) == canonical(parent['entries']), 'Parent prefix/order changed')
    covered, identities = set(), set()
    for entry in saved['entries']:
        bounds = entry['weight_range']
        require(isinstance(bounds, list) and len(bounds) == 2 and all(type(x) is int for x in bounds), 'Bad range')
        lo, hi = bounds
        require(1 <= lo <= hi <= 8192, 'Range outside nonzero outer support domain')
        weights = set(range(lo, hi+1))
        require(not covered & weights, 'Overlapping weight ranges')
        covered |= weights
        identity = (entry['path'], entry['sha256'])
        require(identity not in identities, 'Duplicate event identity')
        identities.add(identity)
    require(covered == set(range(1, 8193)), 'Incomplete active-group coverage')
    require({tuple(x['weight_range']) for x in saved['entries'][265:]} == set(ADDITIONS), 'Unexpected new events')
    for key, value in dict(receipt_count=274, covered_weight_count=8192, pending_weight_count=0,
        covered_weight_ranges=[[1, 8192]], pending_weight_ranges=[], added_event_count=9,
        added_weight_count=61, aggregation_precision_bits=384, covered_events_failure_bits=40,
        directed_aggregation=True, meets_covered_events_target=True, terminal_counted_once=True).items():
        require(saved.get(key) == value, 'Index counter/scope differs: '+key)
    tails = [entry for entry in saved['entries'] if entry['weight_range'][1] >= 6678]
    require(len(tails) == 1 and tails[0]['weight_range'] == [6678, 8192], 'Terminal suffix must occur once')
    tail = tails[0]
    require(tail['path'] == TERMINAL and tail['sha256'] == PINS[TERMINAL], 'Terminal receipt binding differs')
    require(tail['log2_upper_rational'] == terminal['log2_upper_rational'], 'Terminal ceiling changed')
    terminal_expected = dict(state_dimension=30, inner_width=64, forbidden_weight=838860,
        weight_range=[6678, 8192], gamma_extension_degree=16, gamma_base_field_order=256,
        feedback='gamma A transpose', CA_zero=True, fixed_nodes='GF256 additive subspace0..63',
        state_nodes='128..157', refresh='D2 H_30 D1 with independent full-field nonzero diagonals',
        initial_state='zero', final_state='unrestricted', state_persists_across_regions=True,
        outer_independent_of_route_and_inner=True, full_field_random_coefficients=True,
        outer_groups=8192, outer_randomizers_per_group=64, pair_swaps=True,
        outer_randomizer='independent uniform invertible4x4 Toeplitz per group and column',
        method='inverse-image-terminal-range', directed=True,
        covers_all_weights_in_range_simultaneously=True, no_union_over_tail_weights=True)
    require(all(terminal.get(k) == v for k, v in terminal_expected.items()), 'Wrong terminal ensemble or forbidden threshold')
    require(838860 == (2 * (1 << 21)) // 5 and F(838861, 1 << 21) > F(2, 5), 'Threshold arithmetic')
    bound, enclosure = positive_union(saved['entries'])
    require(F(bound) <= -40, 'Actual event union misses the 40-bit target')
    require(F(bound) == F(saved['covered_events_union_log2_upper_rational']), 'Saved event union differs')
    return dict(event_count=274, positive_event_count=sum(not x.get('zero', False) for x in saved['entries']),
        covered_active_group_counts=[1, 8192], covered_weight_count=8192, pending_weight_count=0,
        nonoverlapping_weight_ranges=True, probability_events_assumed_disjoint=False,
        parent_entries_and_order_preserved=265, terminal_event_count=1, terminal_weight_range=[6678, 8192],
        independent_positive_union_precision_bits=384, actual_union_log2_upper_rational=bound,
        actual_union_log2_upper_enclosure=enclosure, meets_global_40_bit_target=True)


def theorem():
    return dict(field='GF(2^128) with a fixed GF256 embedding', source_dimension=1 << 20,
        output_length=1 << 21, rate='1/2', state_dimension=30, inner_width=64, forbidden_weight=838860,
        minimum_distance_lower_bound=838861, relative_symbol_distance_lower_bound='838861/2097152',
        setup_distribution=[
            '8192 groups of four RS[64,32]; uniform invertible4x4 Toeplitz per group/column.',
            '256 uniform permutations of4096 pair slots, fair swaps, alternating cyclic matchings.',
            'Nonzero uniform full-field D1,D2 entries per packet. All setup choices independent and fixed across messages.'],
        fixed_inner='A=RS[64,30] from GF256 nodes128..157 to0..63; C=A^T, CA=0; H30 is the totally nonsingular Lagrange map.',
        feedback='Fix gamma of degree16 over GF256. Y=X+AQ; Q_next=D2 H30 D1(Q+gamma C X).',
        state_boundary='Initial zero, continuous across regions, final unrestricted.',
        exhaustion='r counts active outer groups, not message symbols. Injectivity excludes nonzero messages at r=0.',
        linearity='Injective outer and causal identity feedthrough give F-linear dimension K.',
        bad_event='There exists a nonzero F-message whose encoded output has symbol Hamming weight at most838860.',
        probability_space='One ideal setup; exists-message quantification inside the failure event.',
        conclusion='Under reviewed inequalities and audit evidence, failure probability <= actual event union <=2^-40.',
        counting_conventions=['All original-r ordered profiles, including odd/axes/dense cases; only counts use complementation.',
            'Counts outside alpha; projective q-1 divisor included once.',
            'Witness minima follow each full prefactor and chronological transfer; no witness-union factor.',
            'One terminal event; no random-gamma failure term; Toeplitz conditioning already charged.'])


def check(index_sha256, selection_sha256):
    require(index_sha256 == PINS[INDEX] and selection_sha256 == PINS[SELECTION], 'Unexpected final index/selection')
    inputs, pins = {}, source_pins()
    def read(name, digest):
        require(isinstance(digest, str) and re.fullmatch('[0-9a-f]{64}', digest), 'Invalid explicit SHA256')
        path = relative(name)
        raw = path.read_bytes()
        require(sha256(raw).hexdigest() == digest, 'Changed compact input: '+name)
        inputs[path] = raw
        return json.loads(raw) if name.endswith('.json') else None
    saved, choice = read(INDEX, index_sha256), read(SELECTION, selection_sha256)
    parent, terminal = read(PARENT, PINS[PARENT]), read(TERMINAL, PINS[TERMINAL])
    read(PARENT_CHECK, PINS[PARENT_CHECK])
    require(choice['schema'] == 'gamma-v31-explicit-remote-evidence-selection-v1', 'Wrong selection schema')
    expected_pins = dict(parent=PARENT, parent_independent_check=PARENT_CHECK,
        validator='hd_field128_gamma_selected_index_v31.py', validator_tests='hd_field128_gamma_selected_index_v31_tests.py',
        recorder='hd_field128_gamma30_d40_v31_record_20261010.py')
    for key, name in expected_pins.items():
        require(choice[key] == dict(path=name, sha256=PINS[name]), 'Wrong assembly pin: '+key)
        read(name, PINS[name])
    original = read(**dict(name=choice['input_selection']['path'], digest=choice['input_selection']['sha256']))
    specs = choice['event_specs']
    require(original['event_specs'] == specs and len(specs) == 9, 'Selection proposal differs')
    require(choice['whole_code_certificate'] is False and choice['completion_audit_performed'] is False,
            'Selection was relabeled')
    for spec in specs:
        require(set(spec) == {'audit', 'archive_receipt'}, 'Unexpected specification')
        for pin in spec.values():
            require(set(pin) == {'path', 'sha256'}, 'Unexpected input pin')
            read(pin['path'], pin['sha256'])
    rebuilt = selected.extend_v30(HERE / PARENT, HERE / PARENT_CHECK, specs)
    checked = verify_objects(saved, parent, terminal, rebuilt)
    for closure in (saved['aggregation_source_sha256'], terminal['source_sha256']):
        selected.frozen._merge(pins, closure)
    selected.frozen.producer.check_source_pins(pins)
    require(source_pins().items() <= pins.items(), 'Checker source changed')
    require(all(path.read_bytes() == raw for path, raw in inputs.items()), 'Input changed at exit')
    return dict(schema='gamma-ideal-setup-distance-completion-v31', audit_result='PASS',
        construction=CONSTRUCTION, **checked, theorem=theorem(),
        completion_audit_performed=True, whole_code_certificate=True, partial_index_whole_code_certificate=False,
        exact_frozen_helper_object_reconstructed=True, same_ensemble_and_threshold_evidence='Parent audit, new-audit dispatch, terminal metadata; producer validity inherited.',
        source_sha256=dict(sorted(pins.items())),
        compact_input_pins=[dict(path=p.relative_to(HERE).as_posix(), sha256=sha256(raw).hexdigest(), bytes=len(raw))
                            for p, raw in sorted(inputs.items())],
        source_and_input_bytes_checked_at_exit=True, reviewed_mathematical_inequalities_inherited=True,
        saved_data_audits_inherited=True, remote_archive_checks_inherited=True,
        coefficient_replay=False, normalized_regional_replay=False, transfer_recomputed=False,
        remote_archive_revalidation=False, remote_network_access=False, historical_execution_authenticated=False,
        concrete_seed_certified=False, binary_hamming_distance_claim=False, encoder_timing_claim=False,
        scope='Ideal GF(2^128) setup and symbol distance; conditional on pinned reviewed producer and audit evidence.')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index-sha256', required=True)
    parser.add_argument('--selection-sha256', required=True)
    parser.add_argument('--output', required=True, help='New theory-relative manifest filename')
    args = parser.parse_args(argv)
    target = relative(args.output)
    require(not target.exists(), 'Refusing to overwrite an existing artifact')
    require(shutil.disk_usage(HERE).free - LIMIT >= RESERVE, 'Preserve20GiB after planned output')
    report = check(args.index_sha256, args.selection_sha256)
    raw = (json.dumps(report, indent=2, allow_nan=False)+'\n').encode()
    require(len(raw) <= LIMIT and shutil.disk_usage(HERE).free-len(raw) >= RESERVE, 'Output budget')
    selected.frozen.producer.check_source_pins(report['source_sha256'])
    for pin in report['compact_input_pins']:
        require(sha256(relative(pin['path']).read_bytes()).hexdigest() == pin['sha256'], 'Input changed before publication')
    with target.open('xb') as stream:
        stream.write(raw)
    require(target.read_bytes() == raw, 'Published manifest bytes differ')
    print(json.dumps(dict(path=target.name, bytes=len(raw), sha256=sha256(raw).hexdigest(),
        actual_union_log2_upper_rational=report['actual_union_log2_upper_rational'],
        covered_weight_count=8192, pending_weight_count=0, whole_code_certificate=True)), flush=True)


if __name__ == '__main__':
    main()
