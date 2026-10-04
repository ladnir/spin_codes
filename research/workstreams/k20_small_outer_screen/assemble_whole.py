"""Authenticate components, minimize each occupancy bound, and sum exactly.

Only complete coverage of all 8192 occupancies can produce a whole-code
certificate. This is an assembly of retained and fresh outward endpoints,
not an independent numerical recomputation of those endpoints.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAP = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
COMMON = dict(K=1048576, N=2097152, threshold=209715, state_bits=22,
    groups=8192, group_dimension=128, group_output_bits=256, regions=64,
    physical_t=64, zero_initial_state=True, final_flush=False)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def dyadic(pair):
    require(isinstance(pair, (list, tuple)) and len(pair) == 2
        and all(type(x) is int for x in pair) and pair[0] > 0
        and abs(pair[1]) <= 10_000_000, 'positive finite dyadic pair required')
    m, e = pair
    zeros = (m & -m).bit_length()-1
    return m >> zeros, e+zeros


def less(a, b):
    a, b = dyadic(a), dyadic(b)
    tops = a[0].bit_length()+a[1], b[0].bit_length()+b[1]
    if tops[0] != tops[1]:
        return tops[0] < tops[1]
    e = min(a[1], b[1])
    return (a[0] << (a[1]-e)) < (b[0] << (b[1]-e))


def exact_sum(values):
    buckets = {}
    for value in values:
        m, e = dyadic(value)
        buckets[e] = buckets.get(e, 0)+m
    require(bool(buckets), 'nonempty union required')
    exponents = sorted(buckets, reverse=True)
    e, m = exponents[0], buckets[exponents[0]]
    for ne in exponents[1:]:
        m = (m << (e-ne))+buckets[ne]
        e = ne
    return dyadic((m, e))


def encode(value):
    m, e = dyadic(value)
    return dict(numerator_hex=hex(m), exponent=e)


def margin(value):
    m, e = dyadic(value)
    return -log2(m)-e


def merge_pins(target, pins):
    require(isinstance(pins, dict) and bool(pins), 'source pins required')
    for path, expected in pins.items():
        require(path not in target or target[path] == expected, f'conflicting source pin: {path}')
        require(digest(path) == expected, f'source pin mismatch: {path}')
        target[path] = expected


def component(record):
    require(all(record.get(k) == v for k, v in COMMON.items()), 'component construction or geometry mismatch')
    require(record.get('whole_code_certificate') is False, 'occupancy component required')
    require(record['map_record']['map_sha256'] == MAP, 'component map mismatch')
    require(record.get('precision', 0) >= 192, 'outward precision >=192 required')
    schema = record.get('schema')
    if schema in ('finite-packet-rs-state-q1-1', 'finite-packet-rs-state-q2-1'):
        q = 1 if schema.endswith('q1-1') else 2
        require(record.get('occupancy_covered') == [q]
            and record.get('all_supports_covered') is True
            and record.get('all_two_group_support_pairs_covered') is (q == 2)
            and record.get('inner_distribution') == 'uniform_gl'
            and record.get('outer_symbol_group') == 'GL(16,2)'
            and record.get('count_kind') == 'exact_expected_shells', 'sparse scope mismatch')
        values = {str(q): record[f'q{q}_upper']}
    elif schema == 'k20-small-outer-fractional-range-1':
        require(record.get('continuous_state') is True
            and record.get('slots_per_region') == 8192
            and record.get('macros_per_region') == 256
            and record.get('macro_slots') == 32
            and record.get('physical_steps_per_macro') == 2
            and 0 < Fraction(record['alpha']) <= 1
            and Fraction(record['tilt']) > 0, 'fractional range premises mismatch')
        require(digest(record['local_receipt']) == record['local_sha256'], 'local receipt mismatch')
        values = record['occupancy_uppers']
    elif schema == 'rs16-state-ladder-tail-1':
        require(record.get('every_requested_occupancy_checked') is True
            and record.get('every_shell_checked') is True
            and record.get('inner_distribution') == 'uniform_gl'
            and record.get('outer_symbol_group') == 'GL(16,2)'
            and record.get('state_continuity') == 'retained_between_every_physical_step_and_region'
            and record.get('method') in ('exact', 'fugacity'), 'ordinary tail premises mismatch')
        values = record['occupancy_uppers']
    elif schema == 'k20-small-outer-fractional-dense-complete-1':
        require(record.get('continuous_state') is True
            and record.get('slots_per_region') == 8192
            and record.get('macros_per_region') == 256
            and record.get('physical_steps_per_macro') == 2
            and record.get('every_requested_occupancy_checked') is True
            and Fraction(record['beta']) == Fraction(1 << 256, 65535**8),
            'dense fractional premises mismatch')
        require(record.get('alphas') and all(0 < Fraction(a) <= 1 for a in record['alphas'])
            and record.get('marker_probabilities')
            and all(0 < Fraction(p) <= 1 for p in record['marker_probabilities']),
            'valid dense witnesses required')
        for path in record['local_receipts']:
            require(record['source_sha256'].get(path) == digest(path), 'dense local receipt pin missing')
        values = record['occupancy_uppers']
        require(set(record['occupancy_choices']) == set(values), 'dense witness coverage mismatch')
        for q, choice in record['occupancy_choices'].items():
            require(choice['local_receipt'] in record['local_receipts']
                and choice['alpha'] in record['alphas']
                and choice['marker_probability'] in record['marker_probabilities']
                and (Fraction(choice['marker_probability']) < 1 or q == '8192'),
                'invalid selected dense witness')
    else:
        raise ValueError(f'unrecognized component schema: {schema}')
    qs = record.get('occupancy_covered')
    require(isinstance(qs, list) and qs and qs == sorted(set(qs))
        and all(type(q) is int and 1 <= q <= 8192 for q in qs)
        and set(values) == {str(q) for q in qs}, 'explicit component coverage mismatch')
    for pair in values.values():
        dyadic(pair)
    return values


def assemble(paths, output, target_bits=50):
    output = Path(output)
    require(not output.exists() and target_bits in (40, 50), 'fresh output and supported target required')
    values, selected, pins, refs = {}, {}, {}, []
    map_record = None
    for index, path in enumerate(paths):
        path = Path(path).resolve()
        record = json.loads(path.read_text(encoding='utf-8'))
        current = component(record)
        if map_record is None:
            map_record = record['map_record']
        require(map_record == record['map_record'], 'full map records differ')
        merge_pins(pins, record['source_sha256'])
        merge_pins(pins, {str(path): digest(path)})
        refs.append(dict(path=str(path), sha256=digest(path), schema=record['schema']))
        for q, pair in current.items():
            if q not in values or less(pair, values[q]):
                values[q], selected[q] = pair, index
    missing = sorted(set(range(1,8193)) - {int(q) for q in values})
    require(not missing, f'incomplete coverage: {len(missing)} missing occupancies; first {missing[:20]}')
    union = exact_sum(values.values())
    tail = exact_sum(v for q, v in values.items() if int(q) >= 3)
    met = less(union, (1, -target_bits))
    merge_pins(pins, {str(path.resolve()): digest(path)
                      for path in (Path(__file__), HERE/'METHODS.md')})
    result = dict(COMMON, schema='k20-small-outer-whole-1', distance='1/10',
        minimum_distance=209716, slots_per_region=8192, macro_slots=32,
        macros_per_region=256, physical_steps_per_macro=2, continuous_state=True,
        map_sha256=MAP, map_record=map_record, components=refs, source_sha256=pins,
        occupancy_uppers={str(q): values[str(q)] for q in range(1,8193)},
        selected_components={str(q): selected[str(q)] for q in range(1,8193)},
        exact_union=encode(union), exact_tail_union=encode(tail),
        margin_bits_display=margin(union), tail_margin_bits_display=margin(tail),
        target_bits=target_bits, target_met=met, whole_code_certificate=met,
        all_occupancies_covered=True, component_endpoints_recomputed_by_assembler=False,
        scope='Whole-code setup-failure bound for the fixed ideal independent ensemble. '
              'Authenticated component endpoints, exact per-occupancy minima and exact union; '
              'no claim about the minimum distance of every seed.')
    merge_pins({}, pins)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k:result[k] for k in ('target_met','margin_bits_display','tail_margin_bits_display')}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--target-bits', type=int, choices=[40,50], default=50)
    parser.add_argument('components', nargs='+', type=Path)
    args = parser.parse_args()
    assemble(args.components, args.output, args.target_bits)
