"""Final M22 audit. No solver, statistical test, or frozen-file overwrite.

Rebuilds the exact joint LP witness, independently recomputes all 92 original
one-row transfers, checks deterministic caps and higher-occupation aggregates,
and recursively checks retained evidence. The cell-DP replay is a separate
documented command. --record creates the final audit receipt once.
"""
import argparse
import itertools
import json
import math
import sys
from fractions import Fraction as Q
from pathlib import Path

sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
import numpy as np
from flint import arb, ctx
from audit_bch_q1_full_arb import decode, encode, rational
from audit_bch_closure_envelope import audit as envelope_audit
from bch_joint_shell_objective import audit as joint_audit
from certify_bch_m25_closure import deterministic_caps
from export_constant_weight_delsarte import attach_qsopt_certificate
from run_higher_endpoint_preflight import sha, write_new

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / 'generated'
OUTPUT = GEN / 'bch256_m22_unconditional_closure_audit.json'
JOINT = GEN / 'shift_rank_oa29_joint'


def independent_m22(s):
    """Normalized support recurrence; region polynomial binary powering.

    This uses neither regions() nor support_moments() from the primary
    unnormalized transfer implementation. All arithmetic is 256-bit Arb.
    """
    t = arb(s.numerator) / s.denominator
    bit = (1 + (-t).exp()) / 2
    d = bit / (1 << 22)
    f = bit - d
    zero = [[arb(1), arb(0)], [d, f]]
    one = [[d, f], [d, f]]

    def mul(a, b):
        return [[sum((a[i][k] * b[k][j] for k in range(2)), arb(0))
                 for j in range(2)] for i in range(2)]

    def add(a, b):
        return [[a[i][j] + b[i][j] for j in range(2)] for i in range(2)]

    def polymul(a, b):
        return [mul(a[0], b[0]), add(mul(a[0], b[1]), mul(a[1], b[0]))]

    ident = [[arb(1), arb(0)], [arb(0), arb(1)]]
    empty = [[arb(0), arb(0)], [arb(0), arb(0)]]
    power, result, exponent = [zero, one], [ident, empty], 8192
    while exponent:
        if exponent & 1:
            result = polymul(result, power)
        exponent >>= 1
        if exponent:
            power = polymul(power, power)
    rz = result[0]
    ra = [[v / 8192 for v in row] for row in result[1]]
    rows = [[arb(1), arb(0)]]
    for length in range(1, 257):
        updated = []
        for w in range(length + 1):
            item = [arb(0), arb(0)]
            if w < length:
                item = [sum((rows[w][k] * rz[k][j] for k in range(2)), arb(0))
                        * (length - w) / length for j in range(2)]
            if w:
                item = [item[j] + sum((rows[w-1][k] * ra[k][j] for k in range(2)), arb(0))
                        * w / length for j in range(2)]
            updated.append(item)
        rows = updated
    return [sum(row) * (209716 * t).exp() * 8192 for row in rows]


def build():
    cache = {}
    hashed = {}

    def checked(path):
        path = path.resolve()
        if path in cache:
            return cache[path]
        value = json.loads(path.read_text())
        cache[path] = value
        aliases = {'diagnostic': GEN / 'random_inner_oa15_majorant_diagnostic.json',
                   'johnson_w52': GEN / 'johnson_n256_w52_d38.json',
                   'johnson_w54': GEN / 'johnson_n256_w54_d38.json'} if path.name == 'random_inner_threshold_outward.json' else {}

        def walk(obj):
            if isinstance(obj, dict):
                for name, expected in obj.get('source_sha256', {}).items():
                    dependency = aliases.get(name, ROOT / name)
                    if path.name == 'bch256_q2_envelope_diagnostic.json':
                        assert Path(name).name == name
                        dependency = GEN / 'bch256_q2_envelope_diagnostic_sources' / name
                    dependency = dependency.resolve()
                    actual = sha(dependency)
                    assert actual == expected, (str(path), name)
                    hashed[dependency] = actual
                    if dependency.suffix == '.json':
                        checked(dependency)
                for key, child in obj.items():
                    if key != 'source_sha256':
                        walk(child)
            elif isinstance(obj, list):
                for child in obj:
                    walk(child)
        walk(value)
        return value

    full = checked(JOINT / 'audit.json')
    assert joint_audit('prepare_bch_shift_rank_oa29_probe', JOINT) == full
    print('Exact 1163-row joint LP and reconstructed algebra passed', flush=True)
    envelope = checked(GEN / 'bch256_closure_deterministic_envelope.json')
    assert envelope_audit() == envelope
    assert not any(envelope[k] for k in ('orbit_lower_bounds_used', 'statistical_caps_used', 'lattice_rounding_used'))
    rank = checked(GEN / 'bch256_shift_rank_q30_refined.json')
    assert rank['Qdual_extended_distance_lower'] >= 30
    # This independently proved distance also justifies all inherited OA21 rows;
    # no published dual-distance table is needed as a mathematical assumption.
    for w in (52, 54):
        path = GEN / f'johnson_n256_w{w}_d38.json'
        stored = checked(path)
        assert (stored['n'], stored['weight'], stored['minimum_distance']) == (256, w, 38)
        assert attach_qsopt_certificate(stored.copy(), 256, w, 19, path.with_suffix('.sol')) == stored
    caps = deterministic_caps()

    transfer = checked(GEN / 'bch256_q1_full_arb_transfer.json')
    assert (transfer['parameters']['outer_rows'], transfer['parameters']['memory_bits'],
            transfer['parameters']['output_bits'], transfer['parameters']['distance']) == (8192, 22, 1 << 21, 209716)
    assert set(map(int, transfer['coefficient_rows'])) == {w for w, c in caps.items() if c}
    groups = {}
    for name, row in transfer['coefficient_rows'].items():
        assert int(name) == row['weight']
        s = decode(row['surprisal'])
        assert s > 0
        groups.setdefault(s, []).append(row)
    ctx.prec = 256
    shell_count = 0
    for s, rows in groups.items():
        values = independent_m22(s)
        for row in rows:
            upper = min(Q(8192), rational(values[row['weight']].upper()))
            assert upper <= decode(row['coefficient_upper']), row['weight']
            shell_count += 1
    assert shell_count == 92
    print('All 92 M22 one-row Chernoff coefficients independently checked', flush=True)

    r = Q((1 << 20) + 1, 1 << 20)
    higher = []
    for q, field, cap_field in ((2, 'conditional_pair_upper', 'spectrum_envelope'),
                              (3, 'conditional_upper', 'weight_class_caps')):
        stem = 'bch256_q2_positive_outward' if q == 2 else 'bch256_q3_capped_outward'
        receipt = checked(GEN / (stem + '.json'))
        params = receipt['parameters']
        assert (params['memory_bits'], params['distance_cutoff'], params['occupation'], params['outer_rows'],
                params['outer_length'], params['message_bits'], params['output_bits']) == (22, 209716, q, 8192, 256, 1 << 20, 1 << 21)
        old_caps = {int(w): c for w, c in receipt[cap_field].items() if int(w) and c}
        if q == 2:
            assert set(old_caps) == {w for w, c in caps.items() if c}
            assert all(caps[w] <= r * old_caps[w] for w in old_caps)
        else:
            assert set(old_caps) == set(range(38, 72, 2))
            assert old_caps[70] == (1 << 128) - 1
            assert all(caps[w] <= r * old_caps[w] for w in range(38, 70, 2))
        with np.load(GEN / (stem + '.npz')) as data:
            array = data[field]
            value = math.comb(8192, q) * sum((math.prod(old_caps[w] for w in ws) * Q.from_float(float(array[ws]))
                for ws in itertools.product(sorted(old_caps), repeat=q)), Q(0))
        assert value == decode(receipt[f'Q{q}_upper'])
        higher.append(value)
    low_check = checked(GEN / 'bch256_low_occupation_crosschecks.json')
    assert low_check['Q2_zero_weight_axis_checks'] == 514
    assert low_check['Q3_two_zero_weight_axis_checks'] == 213
    assert low_check['both_symmetrized_aggregates_below_2_to_minus_60']

    tail = checked(GEN / 'bch256_occupation_tail_outward.json')
    params = dict(outer_rows=8192, outer_length=256, message_bits=1 << 20,
                  output_bits=1 << 21, distance_cutoff=209716, memory_bits=22)
    assert tail['parameters'] == params == full['parameters']
    assert [row['q'] for row in tail['rows']] == list(range(4, 8193))
    for row in tail['rows']:
        value = decode(row['log2_upper'])
        assert row['dyadic_upper_exponent'] == -(-value.numerator // value.denominator) <= -58
    exponent = min(row['dyadic_upper_exponent'] for row in tail['rows'])
    tail_sum = Q(sum(1 << (row['dyadic_upper_exponent'] - exponent) for row in tail['rows']), 1 << -exponent)
    assert tail_sum == decode(tail['tail_upper']) < Q(1, 1 << 83)
    crosscheck = checked(GEN / 'bch256_tail_transfer_crosschecks.json')
    assert crosscheck['all_8189_occupations_checked']
    assert crosscheck['unnormalized_region_count_checks'] == 83
    assert crosscheck['direct_full_length_power_checks'] == 8106
    assert checked(GEN / 'bch256_tail_envelope_checks.json')['exact_prefix_comparisons_checked'] == 5397
    inflation = 1 / (1 - 8192 * (r - 1))
    assert inflation == Q(128, 127)
    objective = checked(JOINT / 'objective.json')
    higher_total = inflation * (sum(higher, Q(0)) + tail_sum)
    assert higher_total == decode(objective['higher_occupation_upper'])

    cells = checked(GEN / 'bch256_q1_activity_cells_outward.json')
    coefficients = {int(w): decode(row['coefficient_upper']) for w, row in transfer['coefficient_rows'].items()}
    coefficients.update({row['weight']: decode(row['coefficient_upper']) for row in cells['shells']})
    assert {row['weight'] for row in cells['shells']} == {38, 40, 42}
    excluded = {38, 40, 42, 214, 216, 218}
    rest = higher_total + sum((caps[w] * c for w, c in coefficients.items() if w not in excluded), Q(0))
    assert rest == decode(full['constant_rest_upper'])
    total = decode(full['paired_shells_upper']) + rest
    assert total == decode(full['full_M22_first_moment_upper']) < Q(61, 100 * (1 << 40)) < Q(1, 1 << 40)
    assert full['original_M22_target_closed'] and not full['statistical_shell_assumptions_used']
    print('All occupations reaggregated; exact bound below 0.61 times 2^-40', flush=True)

    # Direct inputs plus recursive retained evidence; arrays absent from older
    # receipt manifests are explicitly included in this final receipt.
    direct = [Path(__file__), JOINT / 'audit.json', JOINT / 'objective.json',
              GEN / 'bch256_closure_deterministic_envelope.json',
              GEN / 'bch256_shift_rank_q30_refined.json',
              GEN / 'bch256_q1_activity_cells_outward.json', GEN / 'bch256_q1_activity_cells_outward.npz',
              GEN / 'bch256_q1_full_arb_transfer.json', GEN / 'bch256_low_occupation_crosschecks.json',
              GEN / 'bch256_tail_transfer_crosschecks.json', GEN / 'bch256_tail_envelope_checks.json']
    for path in direct:
        hashed[path.resolve()] = sha(path)
    return dict(classification='Final exact deterministic-spectrum ideal RandomStepConv-M22 closure audit',
                parameters=params, original_M22_target_closed=True,
                full_first_moment_upper=encode(total), margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
                exact_bound_below_0_61_times_2_to_minus_40=True,
                joint_LP_rows_rebuilt_and_verified=1163, independent_M22_transfer_shells_checked=92,
                all_8192_occupations_covered=True, JSON_evidence_files_checked=len(cache),
                statistical_shell_assumptions_used=False, empirical_orbit_lower_bounds_assumed=False,
                external_dual_distance_table_required=False, published_exact_endpoint_spectra_used=True,
                cell_DP_replay_command='python -B code/certify_bch_q1_activity_cells.py --verify',
                higher_transfer_replay_scope='Retained independent arithmetic certificates hash-checked; Q2/Q3 and all tail bounds exactly reaggregated, not recomputed here',
                paper_equation19_transfer_established=False, RM2Sub_transfer_established=False,
                PRG_replacement_established=False, eleven_percent_distance_claimed=False,
                source_sha256={str(p.relative_to(ROOT)): h for p, h in sorted(hashed.items())})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    if args.record:
        assert not OUTPUT.exists()
    value = build()
    if args.record:
        write_new(OUTPUT, value)
    else:
        assert value == json.loads(OUTPUT.read_text())
    print(json.dumps({k: v for k, v in value.items() if k not in ('source_sha256', 'full_first_moment_upper')}, indent=2))
