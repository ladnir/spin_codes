"""Independent whole-receipt check with sequential regional row propagation.

Reuses the audited positive primitives and placement recurrence. The complete
state calculation differs from the producer's binary matrix powers: apply
the regional matrix to the initial row64 times, then sum terminal states.
No producer source is modified and no fresh local census is performed.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from math import comb, log2
from pathlib import Path
from time import monotonic

import numpy as np
import certify_wider24 as producer
import global_outward as go
import outward_positive as op
import scaled_positive as sp
import q1_outward as q1

HERE = Path(__file__).resolve().parent
GEOMETRY = dict(K=65536,N=131072,rate='1/2',outer_groups=256,group_dimension=256,
    group_output_bits=512,symbol_bits=32,regions=64,slots_per_region=256,
    packet_bits=8,physical_t=64,state_bits=24,steps_per_region=32,
    cutoff=13107,zero_initial_state=True,final_flush=False,continuous_state=True)


def decode(record):
    """Decode integer hexadecimal dyadics without any floating conversion."""
    if type(record.get('exponent')) is not int or not isinstance(record.get('numerator_hex'), str):
        raise ValueError('integer exponent and integer hexadecimal numerator required')
    numerator = int(record['numerator_hex'], 16)
    if numerator < 0:
        raise ValueError('nonnegative dyadic required')
    return F(numerator)*F(2)**record['exponent']


def encode(value):
    value = F(value)
    if value < 0 or value.denominator & (value.denominator-1):
        raise ValueError('nonnegative dyadic required')
    return dict(numerator_hex=hex(value.numerator), exponent=1-value.denominator.bit_length())


def margin(value):
    value = F(value)
    return log2(value.denominator)-log2(value.numerator) if value else float('inf')


def cap(value):
    return max(F(1, 1 << 200), min(F(1), F(value)))


def check_endpoint_union(records, recorded_union, *, groups=256, target=40):
    if set(records) != {str(q) for q in range(1, groups+1)}:
        raise ValueError('all occupancy endpoints required exactly once')
    endpoints = {int(q): decode(item) for q, item in records.items()}
    if any(not 0 <= value <= 1 for value in endpoints.values()):
        raise ValueError('probability endpoints required')
    total = sum(endpoints.values(), F(0))
    if total != decode(recorded_union):
        raise ArithmeticError('stored union is not the exact endpoint sum')
    if total > F(1, 1 << target):
        raise ArithmeticError('exact stored union exceeds target')
    return endpoints, total


def authenticate(pins):
    if not pins:
        raise ValueError('nonempty source pins required')
    for name, digest in pins.items():
        path = Path(name)
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ArithmeticError(f'pin mismatch: {name}')


def sequential_moment(regional, regions=64):
    if type(regions) is not int or regions < 1:
        raise ValueError('positive region count required')
    size = regional.value.shape[0]
    initial = np.zeros((1, size))
    initial[0, 0] = 1.
    row = sp.Scaled(initial)
    for _ in range(regions):
        row = sp.matmul(row, regional)
    return sp.matmul(row, sp.Scaled(np.ones((size, 1))))


def sequential_bound(regional, q, z, alpha, *, groups=256, regions=64,
                     cutoff=13107, beta=None, a=None):
    """Full row-state bound; scalar prefactor follows the same proved formula."""
    z, alpha = F(z), F(alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1 or not 1 <= q <= groups:
        raise ValueError('valid witness and occupancy required')
    beta = F(1 << 512, ((1 << 32)-1)**8) if beta is None else F(beta)
    a = ((1+z)/2)**8 if a is None else F(a)
    scalar = sp.integer_power(beta, q)
    scalar = sp.multiply(scalar, sp.integer_power(1/z, cutoff))
    scalar = sp.multiply(scalar, sp.integer_power(a, regions*q))
    scalar = sp.fractional_power(scalar, alpha)
    bound = sp.multiply(sequential_moment(regional, regions), scalar)
    return sp.multiply(bound, comb(groups, q))


def check_selection(saved, endpoints):
    trials = {str(Path(t['macro_receipt']).resolve()):t for t in saved['trials']}
    choices = saved['selected_witnesses']
    if set(choices) != {str(q) for q in range(1, 257)}:
        raise ValueError('complete witness selection required')
    for q in range(2, 257):
        all_bounds = []
        for trial in trials.values():
            if set(trial['endpoints']) != {str(k) for k in range(2, 257)}:
                raise ValueError('complete trial endpoint collection required')
            value = decode(trial['endpoints'][str(q)])
            if not 0 <= value <= 1:
                raise ValueError('trial probability endpoint out of range')
            all_bounds.append(value)
        best = min([F(1)]+all_bounds)
        if endpoints[q] != best:
            raise ArithmeticError('selected endpoint is not the recorded minimum')
        choice = choices[str(q)]
        if choice == 'trivial probability cap':
            if best != 1:
                raise ArithmeticError('nontrivial endpoint lacks a selected witness')
        elif (str(Path(choice).resolve()) not in trials
              or decode(trials[str(Path(choice).resolve())]['endpoints'][str(q)]) != best):
            raise ArithmeticError('selected macro does not supply the selected endpoint')
    support_trials = [tuple(map(decode, t['q1_support_uppers'])) for t in saved['trials']
                      if 'q1_support_uppers' in t]
    if not support_trials or any(len(v) != 65 for v in support_trials):
        raise ValueError('complete q1 support witnesses required')
    minimum = tuple(min(v) for v in zip((F(1),)*65, *support_trials))
    if tuple(map(decode, saved['q1_selected_support_uppers'])) != minimum:
        raise ArithmeticError('q1 support minimum differs from recorded trials')
    if endpoints[1] != cap(q1.combine_supports(minimum)):
        raise ArithmeticError('selected q1 endpoint differs from exact shell combination')
    return trials


def own_pins():
    pins = producer.source_pins()
    for path in (Path(__file__).resolve(), HERE/'test_verify_whole.py'):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def run(args):
    if args.output.exists():
        raise ValueError('fresh verification output required')
    op.check_runtime()
    source, output = args.receipt.resolve(), args.output.resolve()
    saved = json.loads(source.read_text(encoding='utf-8'))
    if (saved.get('schema') != 'packet8-wider24-whole-code-outward-1'
            or saved.get('geometry') != GEOMETRY or saved.get('target_margin_bits') != 40
            or not saved.get('source_pins_verified_at_finish')
            or not saved.get('whole_code_certificate') or not saved.get('all_occupancies_checked')
            or not saved.get('outward_under_stated_arithmetic_contract')):
        raise ValueError('completed same-construction outward receipt required')
    authenticate(saved['source_sha256'])
    pins = own_pins()
    for name, digest in producer.source_pins().items():
        if saved['source_sha256'].get(name) != digest:
            raise ArithmeticError('whole receipt omits a current producer dependency')
    producer.merge_pins(pins, saved['source_sha256'])
    pins[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
    authenticate(pins)
    endpoints, stored_total = check_endpoint_union(saved['selected_endpoints'], saved['exact_union'])
    trials = check_selection(saved, endpoints)
    selected = defaultdict(list)
    for q in range(2, 257):
        choice = saved['selected_witnesses'][str(q)]
        if choice == 'trivial probability cap':
            raise ArithmeticError('accepted receipt unexpectedly relies on probability1')
        selected[str(Path(choice).resolve())].append(q)

    started = monotonic()
    independent = {}
    records = []
    # Choose the certified local corresponding to the known near-one witness.
    q1_local, q1_z = None, None
    q1_distance = None
    for name, qs in selected.items():
        path = Path(name)
        if pins.get(str(path)) != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ArithmeticError('selected macro is not authenticated by the whole receipt')
        macro, active, z, alpha, macro_saved = producer.load_macro(path)
        distance = abs(float(z)-.9900498337491681)
        if q1_distance is None or distance < q1_distance:
            q1_local, q1_z, q1_distance = active, z, distance
        regional = go.placement(macro, max(qs), epochs=8, windows=32)
        largest_difference = 0.
        for q in qs:
            value = cap(sp.scalar_fraction(sequential_bound(regional[q], q, z, alpha)))
            independent[q] = value
            largest_difference = max(largest_difference, abs(margin(value)-margin(endpoints[q])))
        records.append(dict(macro_receipt=name, occupancies=qs,
            largest_display_margin_difference=largest_difference))
        print(f'sequential replay: {len(qs)} selected q values, alpha={alpha}, '
              f'z={float(z):.7g}; elapsed={monotonic()-started:.2f}s', flush=True)

    # The .01 local may have no selected q>=2; inspect other authenticated
    # trials to select the nearest witness without rerunning any local census.
    for name in trials:
        if name in selected:
            continue
        path = Path(name)
        if pins.get(str(path)) != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ArithmeticError('q1 macro is not authenticated by the whole receipt')
        _, active, z, _, _ = producer.load_macro(path)
        distance = abs(float(z)-.9900498337491681)
        if q1_distance is None or distance < q1_distance:
            q1_local, q1_z, q1_distance = active, z, distance
    if q1_local is None:
        raise ArithmeticError('no authenticated local witness for q1')
    support = q1.q1_support_upper(q1_local, q1_z)
    independent[1] = cap(q1.combine_supports(support))
    if set(independent) != set(range(1, 257)):
        raise ArithmeticError('independent replay missed an occupancy')
    total = sum(independent.values(), F(0))
    if total > F(1, 1 << 40):
        raise ArithmeticError('independent outward union fails target')
    authenticate(pins)
    result = dict(schema='packet8-wider24-independent-sequential-replay-1',
        whole_receipt=str(source), source_sha256=pins, source_pins_verified_at_finish=True,
        recorded_exact_union_verified=True, all_occupancies_recomputed=True,
        regional_composition='64 sequential row-vector updates; no state resets',
        arithmetic_scope='same authenticated local/macros and positive primitives; independent global row composition',
        q1_witness_z=str(q1_z), q1_support_uppers=[encode(x) for x in support],
        q1_margin_bits_display=margin(independent[1]),
        endpoints={str(q):encode(independent[q]) for q in range(1, 257)},
        exact_union=encode(total), exact_union_at_most_2_minus40=True,
        combined_margin_bits_display=margin(total), producer_margin_bits_display=margin(stored_total),
        difference_bits_display=margin(total)-margin(stored_total),
        trials=records, elapsed_seconds=monotonic()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(f'INDEPENDENT EXACT UNION PASS: {margin(total):.12f}bits, '
          f'q1={margin(independent[1]):.12f}', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, default=HERE/'whole_wider24_v1.json')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())

