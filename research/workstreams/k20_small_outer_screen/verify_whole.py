"""Independent authentication and exact union audit of K20 small-outer receipts.

Uses only the standard library, without importing any numerical producer.
Rechecks fixed maps, local receipt structure, every component's provenance,
coverage, deterministic minimum selection, and all 8192 dyadic summands.
Sparse shell folds are recomputed with exact rationals. State censuses and
matrix endpoints are authenticated, not independently regenerated here.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import comb, isclose, log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MAP = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
DEFAULT_RECEIPT = ROOT / 'tmp/k20-small-outer-whole-v1.json'
BETA = F(1 << 256, 65535**8)
CONTINUITY = 'retained_between_every_physical_step_and_region'
COMMON = dict(K=1048576, N=2097152, threshold=209715, distance='1/10',
    groups=8192, group_dimension=128, group_output_bits=256, regions=64,
    state_bits=22, physical_t=64, macros_per_region=256,
    physical_steps_per_macro=2, zero_initial_state=True, final_flush=False)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fields(record, expected, context):
    require(isinstance(record, dict), f'{context}: object required')
    require(all(type(record.get(key)) is type(value) and record[key] == value
                for key, value in expected.items()), f'{context}: geometry or scope mismatch')


def dyadic(pair):
    require(isinstance(pair, list) and len(pair) == 2, 'dyadic pair required')
    n, e = pair
    require(type(n) is int and n >= 0 and type(e) is int,
            'nonnegative integer dyadic pair required')
    return F(n << e) if e >= 0 else F(n, 1 << -e)


def exact_union(encoded):
    require(isinstance(encoded, dict) and set(encoded) == {'numerator_hex', 'exponent'},
            'hexadecimal exact union required')
    n, e = encoded['numerator_hex'], encoded['exponent']
    require(type(n) is str and n.startswith('0x') and type(e) is int,
            'hexadecimal numerator and integer exponent required')
    numerator = int(n, 16)
    require(numerator > 0 and hex(numerator) == n, 'positive canonical hexadecimal numerator required')
    return dyadic([numerator, e])


def sum_dyadics(pairs):
    """Align to one denominator before summing; retain even the smallest term."""
    pairs = tuple(pairs)
    require(pairs, 'nonempty dyadic union required')
    for pair in pairs:
        require(isinstance(pair, list) and len(pair) == 2
                and type(pair[0]) is int and pair[0] > 0 and type(pair[1]) is int,
                'positive integer dyadic pair required')
    exponent = min(pair[1] for pair in pairs)
    numerator = sum(n << (e-exponent) for n, e in pairs)
    return dyadic([numerator, exponent])


def margin(value):
    require(value > 0, 'positive union required')
    return log2(value.denominator) - log2(value.numerator)


def display(record, key, value):
    require(type(record.get(key)) in (int, float)
            and isclose(record[key], margin(value), rel_tol=0, abs_tol=1e-8),
            f'{key}: displayed margin mismatch')


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


class SourceAudit:
    def __init__(self):
        self.observed = {}
        self.locals = {}

    def digest(self, path):
        path = str(Path(path).resolve())
        if path not in self.observed:
            self.observed[path] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        return self.observed[path]

    def pins(self, pins):
        require(isinstance(pins, dict) and pins, 'nonempty source pins required')
        for path, digest in pins.items():
            require(type(path) is str and Path(path).is_absolute()
                    and type(digest) is str and len(digest) == 64,
                    'absolute source path and SHA256 required')
            require(self.digest(path) == digest, f'source pin mismatch: {path}')

    def required(self, pins, path):
        path = str(Path(path).resolve())
        require(pins.get(path) == self.digest(path), f'required source pin missing: {path}')

    def inherited(self, pins, origin):
        require(all(pins.get(path) == digest for path, digest in origin.items()),
                'inherited source pin missing or conflicting')

    def stable(self):
        for path, digest in self.observed.items():
            require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest,
                    f'source changed during audit: {path}')


def audit_map(record, sources, pins):
    fields(record, dict(t=64, s=22, map_sha256=MAP, distribution='uniform_gl',
        return_denominator=4194303, zero_initial_state=True, final_flush=False,
        state_continuity=CONTINUITY, full_state_census=True, state_count=4194304), 'map')
    rows, columns = record.get('expansion_rows_hex'), record.get('feedback_columns')
    require(isinstance(rows, list) and len(rows) == 22 and isinstance(columns, list)
            and len(columns) == 64, 'complete map arrays required')
    rows = [int(value, 16) for value in rows]
    require(all(0 <= value < 1 << 64 for value in rows)
            and all(type(value) is int and 0 <= value < 1 << 22 for value in columns),
            'map coordinates out of range')
    payload = dict(bits=22, width=64, expansion_rows=list(map(hex, rows)), feedback_columns=columns)
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    require(digest == MAP, 'explicit map arrays mismatch')
    sources.pins(record['source_sha256'])
    sources.inherited(pins, record['source_sha256'])


def audit_local(path, sources, inherited_pins, expected_map):
    path = str(Path(path).resolve())
    sources.required(inherited_pins, path)
    if path in sources.locals:
        local = sources.locals[path]
        require(local['map_record'] == expected_map, 'local and component map records differ')
        sources.inherited(inherited_pins, local['source_sha256'])
        return local
    local = load_json(path)
    fields(local, dict(schema='k18-certified22-outward-physical-export-1',
                       whole_code_certificate=False), 'local receipt')
    require(type(local.get('precision')) is int and local['precision'] >= 192,
            'local precision must be at least 192')
    require(local['map_record'] == expected_map, 'local and component map records differ')
    require(F(local['tilt']) > 0 and 0 <= F(local['local_selection_activity']) <= 1,
            'invalid local tilt or row-selection activity')
    sources.pins(local['source_sha256'])
    sources.inherited(inherited_pins, local['source_sha256'])
    audit_map(local['map_record'], sources, local['source_sha256'])
    fields(local['physical_metadata'], dict(bits=22, windows=16, distance=16,
        map_sha256=MAP, distribution='uniform_gl', birth_density='capped', birth_degree=8,
        birth_class_levels=[16, 24, 28, 32, 36, 40, 48, 64], map_images={'length':4194304},
        columns={'length':64}, records={'shape':[25,5], 'dtype':'int64'},
        birth_class_census={'shape':[8,25], 'dtype':'int64'}), 'local physical metadata')
    raw = local['matrices']
    require(isinstance(raw, list) and len(raw) == 17 and all(isinstance(matrix, list)
        and len(matrix) == 11 and all(isinstance(row, list) and len(row) == 11 for row in matrix)
        for matrix in raw), 'physical family must have 17 matrices of shape 11 by 11')
    for matrix in raw:
        for row in matrix:
            for pair in row:
                dyadic(pair)
    require(dyadic(raw[0][0][0]) == 1 and all(dyadic(pair) == 0 for pair in raw[0][0][1:]),
            'zero-input zero-state transition mismatch')
    choices = local['birth_row_choices']
    require(isinstance(choices, list) and len(choices) == 17
            and all(type(index) is int and index >= 0 for index in choices),
            'complete physical birth-row choices required')
    sources.locals[path] = local
    return local


def coverage(record, minimum=3):
    qs = record.get('occupancy_covered')
    values = record.get('occupancy_uppers')
    require(isinstance(qs, list) and qs and all(type(q) is int and minimum <= q <= 8192 for q in qs)
            and qs == sorted(set(qs)), 'explicit sorted unique occupancy coverage required')
    require(isinstance(values, dict) and set(values) == set(map(str, qs)), 'component coverage mismatch')
    require(all(dyadic(pair) > 0 for pair in values.values()), 'positive component endpoints required')
    return values


@lru_cache(maxsize=1)
def exact_shells():
    """Independent inclusion-exclusion for the GF(2^16) MDS[16,8] shells."""
    counts = [F(0) for _ in range(65)]
    power = [1]
    symbol = [0, 60, 1350, 13500, 50625]
    for h in range(1, 17):
        next_power = [0] * (len(power) + 4)
        for i, a in enumerate(power):
            for j in range(1, 5):
                next_power[i+j] += a * symbol[j]
        power = next_power
        if h < 9:
            continue
        words = comb(16, h) * sum((-1)**j * comb(h, j) * (65536**(h-8-j)-1)
                                  for j in range(h-8))
        for v, coefficient in enumerate(power):
            counts[v] += F(words * coefficient, 65535**h)
    require(counts[0] == 0 and sum(counts) == (1 << 128)-1, 'exact shell normalization failed')
    return tuple(counts)


def sparse_fold(record, occupancy):
    counts = exact_shells()
    digest = hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest()
    require(record.get('count_sha256') == digest, 'exact shell hash mismatch')
    if occupancy == 1:
        raw = record['support_probability_uppers']
        require(isinstance(raw, list) and len(raw) == 65, 'all q1 supports required')
        probabilities = list(map(dyadic, raw))
        require(all(0 < p <= 1 for p in probabilities), 'q1 support bound outside (0,1]')
        total = 8192 * sum((count * p for count, p in zip(counts, probabilities)), F(0))
    else:
        raw = record['support_pair_probability_uppers']
        require(isinstance(raw, list) and len(raw) == 65 and all(isinstance(row, list)
                and len(row) == v+1 for v, row in enumerate(raw)), 'all triangular q2 support pairs required')
        total = F(0)
        for v, row in enumerate(raw):
            for u, pair in enumerate(row):
                p = dyadic(pair)
                require(0 < p <= 1, 'q2 support bound outside (0,1]')
                total += (1 if u == v else 2) * counts[u] * counts[v] * p
        total *= comb(8192, 2)
    endpoint = dyadic(record[f'q{occupancy}_upper'])
    # Small slack permits the producer's outward rounding of rational shell
    # coefficients and the sum, while detecting any material fold alteration.
    require(total <= endpoint <= total * (1 + F(1, 1 << (record['precision']-16))),
            f'q{occupancy} endpoint differs from exact shell fold')


def audit_component(record, sources):
    fields(record, dict(COMMON, whole_code_certificate=False), 'component')
    require(type(record.get('precision')) is int and record['precision'] >= 192,
            'component precision must be at least 192')
    pins = record['source_sha256']
    sources.pins(pins)
    audit_map(record['map_record'], sources, pins)
    schema = record.get('schema')
    if schema in ('finite-packet-rs-state-q1-1', 'finite-packet-rs-state-q2-1', 'rs16-state-ladder-tail-1'):
        fields(record, dict(physical_steps=32768, macro_t=128, macro_steps=16384,
            inner_distribution='uniform_gl', inner_group='GL(22,2)', physical_updates_independent=True,
            outer_symbol_group='GL(16,2)', state_continuity=CONTINUITY, fresh_computation=True,
            count_kind='exact_expected_shells'), 'ordinary component')
        require(record.get('count_sha256') == hashlib.sha256(
            json.dumps(list(map(str, exact_shells()))).encode()).hexdigest(), 'outer shell identity mismatch')
        require(isinstance(record.get('tilts'), list) and record['tilts']
                and all(F(t) > 0 for t in record['tilts']), 'positive component tilts required')
        if schema != 'rs16-state-ladder-tail-1':
            q = 1 if schema == 'finite-packet-rs-state-q1-1' else 2
            fields(record, dict(occupancy_covered=[q], all_supports_covered=True,
                               all_two_group_support_pairs_covered=q == 2), 'sparse coverage')
            fields(record['count_premises'], dict(outer='four parallel GF16 RS[16,8] words',
                label_mixing='independent uniform GL16 per aligned four-nibble symbol',
                independent_setups_between_groups=True, exact_expected_shell_counts=True,
                independent_uniform_group_packet_shuffles=True,
                independent_uniform_regional_shuffles=True), 'sparse counting premises')
            sources.required(pins, HERE.parent / 'k16_design/packet_rs_state_sparse.py')
            sparse_fold(record, q)
            return {str(q):record[f'q{q}_upper']}
        fields(record, dict(return_denominator=4194303, birth_density='capped',
            every_requested_occupancy_checked=True, every_shell_checked=True, beta=str(BETA)), 'ordinary tail')
        require(record.get('method') in ('exact', 'fugacity'), 'ordinary tail method required')
        if record['method'] == 'exact':
            require(record.get('placement_backend') in ('linear', 'binary'), 'exact placement backend required')
        else:
            require(record.get('marker_probabilities') and all(0 < F(p) <= 1
                    for p in record['marker_probabilities']), 'valid ordinary marker probabilities required')
        sources.required(pins, HERE.parent / 'k16_design/packet_rs_state_ladder.py')
        return coverage(record)
    fields(record, dict(slots_per_region=8192, continuous_state=True), 'fractional component')
    sources.required(pins, HERE.parent / 'k18_margin60/replay_fractional.py')
    sources.required(pins, HERE.parent / 'k16_design/packet_regional_power.py')
    values = coverage(record)
    if schema == 'k20-small-outer-fractional-range-1':
        fields(record, dict(macro_slots=32, precision=256), 'fractional range')
        require(record['occupancy_covered'] == list(range(min(map(int, values)), max(map(int, values))+1)),
                'fractional range must cover every integer in its interval')
        require(0 < F(record['alpha']) <= 1 and F(record['tilt']) > 0, 'invalid fractional witness')
        path = record['local_receipt']
        require(sources.digest(path) == record['local_sha256'], 'local receipt hash mismatch')
        local = audit_local(path, sources, pins, record['map_record'])
        require(F(local['tilt']) == F(record['tilt']), 'component and local tilts differ')
        sources.required(pins, HERE / 'range_replay.py')
        sources.required(pins, HERE / 'outward_screen.py')
        return values
    require(schema == 'k20-small-outer-fractional-dense-complete-1', 'unsupported component schema')
    fields(record, dict(packet_bits=4, beta=str(BETA), every_requested_occupancy_checked=True,
        all_occupancies_covered=False, local_census_recomputed_by_this_replay=False), 'fractional dense')
    first, last = record['occupancy_interval']
    require(record['occupancy_covered'] == list(range(first, last+1)), 'dense interval coverage mismatch')
    alphas, markers = list(map(F, record['alphas'])), list(map(F, record['marker_probabilities']))
    require(alphas and all(0 < a <= 1 for a in alphas) and markers
            and all(0 < p <= 1 for p in markers), 'invalid dense witnesses')
    locals_by_path = {}
    for path in record['local_receipts']:
        local = audit_local(path, sources, pins, record['map_record'])
        locals_by_path[str(Path(path).resolve())] = local
    require(locals_by_path, 'dense local receipts required')
    choices = record['occupancy_choices']
    require(set(choices) == set(values), 'dense choices must cover every occupancy')
    for q, choice in choices.items():
        path = str(Path(choice['local_receipt']).resolve())
        p = F(choice['marker_probability'])
        require(path in locals_by_path and F(choice['tilt']) == F(locals_by_path[path]['tilt'])
                and F(choice['alpha']) in alphas and p in markers and (p != 1 or q == '8192'),
                f'invalid dense witness at q{q}')
    total = sum_dyadics(values.values())
    require(exact_union(record['exact_union']) == total, 'dense exact union mismatch')
    display(record, 'margin_bits_display', total)
    sources.required(pins, HERE / 'dense_complete.py')
    sources.required(pins, HERE.parent / 'k18_margin60/SLACK_REVIEW.md')
    return values


def audit_record(record):
    sources = SourceAudit()
    fields(record, dict(COMMON, schema='k20-small-outer-whole-1', minimum_distance=209716,
        slots_per_region=8192, macro_slots=32, continuous_state=True, map_sha256=MAP,
        all_occupancies_covered=True), 'whole receipt')
    target = record.get('target_bits')
    require(type(target) is int and target in (40, 50), 'target must be 40 or 50 bits')
    pins = record['source_sha256']
    sources.pins(pins)
    sources.required(pins, HERE / 'assemble_whole.py')
    sources.required(pins, HERE / 'METHODS.md')
    components = record.get('components')
    require(isinstance(components, list) and components, 'component receipts required')
    best, selected, component_maps = {}, {}, []
    for index, entry in enumerate(components):
        path = entry['path']
        require(Path(path).is_absolute(), 'absolute component path required')
        require(sources.digest(path) == entry['sha256'], 'component file hash mismatch')
        sources.required(pins, path)
        component = load_json(path)
        require(entry['schema'] == component.get('schema'), 'component schema mismatch')
        values = audit_component(component, sources)
        sources.inherited(pins, component['source_sha256'])
        component_maps.append(component['map_record'])
        for q, pair in values.items():
            if q not in best or dyadic(pair) < dyadic(best[q]):
                best[q], selected[q] = pair, index
    require(all(value == component_maps[0] for value in component_maps), 'component map records differ')
    require(record.get('map_record') == component_maps[0], 'whole and component map records differ')
    expected_keys = set(map(str, range(1, 8193)))
    require(set(best) == expected_keys and set(record.get('occupancy_uppers', {})) == expected_keys
            and set(record.get('selected_components', {})) == expected_keys, 'whole coverage must be q1..8192')
    for q in expected_keys:
        require(record['occupancy_uppers'][q] == best[q], f'wrong minimum endpoint selection at q{q}')
        require(type(record['selected_components'][q]) is int
                and record['selected_components'][q] == selected[q], f'wrong first-minimum component selection at q{q}')
    total = sum_dyadics(best.values())
    tail = sum_dyadics(pair for q, pair in best.items() if int(q) >= 3)
    require(exact_union(record['exact_union']) == total, 'whole exact union mismatch')
    require(exact_union(record['exact_tail_union']) == tail, 'tail exact union mismatch')
    passes = total < F(1, 1 << target)
    require(type(record.get('target_met')) is bool and record['target_met'] == passes
            and type(record.get('whole_code_certificate')) is bool
            and record['whole_code_certificate'] == passes, 'target or certificate claim differs from exact comparison')
    display(record, 'margin_bits_display', total)
    display(record, 'tail_margin_bits_display', tail)
    sources.stable()
    return dict(schema='k20-small-outer-independent-whole-audit-1', verified=True,
        authenticated_file_count=len(sources.observed), component_count=len(components),
        authenticated_occupancies=8192, authenticated_local_receipts=len(sources.locals),
        exact_selections_match=True, exact_union_matches=True, exact_tail_matches=True,
        target_bits=target, target_met=passes, whole_code_certificate=passes,
        margin_bits=margin(total), tail_margin_bits=margin(tail),
        sparse_shell_folds_independently_recomputed=True, local_census_recomputed=False,
        matrix_endpoints_independently_recomputed=False,
        scope='Independent source and component authentication, explicit-map identity, geometry, '
              'complete occupancy coverage, sparse shell folds, exact selection and union audit. '
              'Local censuses and matrix endpoints are authenticated, not independently regenerated.')


def verify(path):
    raw = Path(path).read_bytes()
    result = audit_record(json.loads(raw))
    result['receipt_sha256'] = hashlib.sha256(raw).hexdigest()
    require(Path(path).read_bytes() == raw, 'whole receipt changed during audit')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', nargs='?', type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    require(args.output is None or not args.output.exists(), 'fresh audit output required')
    result = verify(args.receipt)
    print(json.dumps(result, indent=2))
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, indent=2)
            handle.write('\n')


if __name__ == '__main__':
    main()
