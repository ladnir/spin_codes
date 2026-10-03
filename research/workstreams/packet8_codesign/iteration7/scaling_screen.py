"""Bounded outward K18/K20 screens; never a whole-code certificate.

Only global geometry changes. Authenticated actual24 active local operators
and G4 macros are imported read-only from iteration5. The exact dyadic bounds
cover only the occupancies explicitly listed in each result.
"""
from __future__ import annotations

import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/'iteration5'
sys.path.insert(0, str(OLD))
import certify_wider24 as cert
import global_outward as go
import q1_outward as q1
import scaled_positive as sp


def geometry(power):
    if type(power) is not int or power not in (16, 18, 20):
        raise ValueError('screen geometry is K16, K18, or K20')
    k = 1 << power
    groups = k//256
    return dict(K=k, N=2*k, groups=groups, regions=64, packet_bits=8,
        physical_t=64, state_bits=24, physical_epochs=groups//8,
        macro_epochs=groups//32, macro_windows=32, cutoff=(2*k)//10,
        zero_initial_state=True, continuous_state=True, final_flush=False)


def regional_family(macro, limit, geo):
    return go.placement(macro, limit, epochs=geo['macro_epochs'], windows=32)


def occupancy(regional, q, z, alpha, geo):
    return go.occupancy_upper(regional, q, z, alpha, groups=geo['groups'],
        regions=geo['regions'], cutoff=geo['cutoff'])


def one_group(active, z, geo):
    supports = q1.q1_support_upper(active, z, regions=geo['regions'],
        epochs=geo['physical_epochs'], cutoff=geo['cutoff'])
    return supports, q1.combine_supports(supports, groups=geo['groups'])


def probability_endpoint(value):
    """Keep diagnostic raw margins separate from compact probability caps."""
    exact = sp.scalar_fraction(value)
    return q1.dyadic_upper(max(Fraction(1,1 << 200),min(Fraction(1),exact)))


def own_pins():
    return {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), HERE/'test_scaling_screen.py')}


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    if not args.q or min(args.q) < 2 or args.seconds <= 0:
        raise ValueError('nonempty q>=2 and positive time budget required')
    source = json.loads(args.whole.read_text(encoding='utf-8'))
    if not (source.get('whole_code_certificate') and cert.checked_pins(source['source_sha256'])):
        raise ValueError('authenticated original K16 certificate required')
    pins = own_pins()
    cert.merge_pins(pins, source['source_sha256'])
    pins[str(args.whole.resolve())] = hashlib.sha256(args.whole.read_bytes()).hexdigest()
    paths = [] if args.only_extra else [Path(trial['macro_receipt']) for trial in source['trials']]
    paths += args.extra_macros
    if not paths:
        raise ValueError('at least one macro required')
    paths.sort(key=lambda p: -float(Fraction(json.loads(p.read_text())['z'])))
    started = monotonic()
    receipt = dict(schema='packet8-wider24-scaling-SCREEN-1',
        whole_code_certificate=False, all_occupancies_checked=False,
        outward_partial_bounds_under_stated_arithmetic_contract=True,
        probability_failure_not_established=True, source_sha256=pins,
        source_pins_verified_at_finish=False, results={}, trials=[],
        budget_seconds=args.seconds, stopped_for_budget=False)
    qset = sorted(set(args.q))
    best = {power: {q: Fraction(1) for q in qset} for power in args.powers}
    raw_best = {power: {q: float('-inf') for q in qset} for power in args.powers}
    supports = {power: tuple(Fraction(1) for _ in range(65)) for power in args.powers}
    choices = {power: {} for power in args.powers}
    seen = set()
    cert.op.check_runtime()
    for path in paths:
        saved_header = json.loads(path.read_text(encoding='utf-8'))
        z = Fraction(saved_header['z'])
        # Selecting the witness subset is not a numerical proof operation.
        import math
        theta = -math.log(float(z))
        if theta > args.max_theta + 1e-10:
            continue
        macro, active, z, alpha, saved = cert.load_macro(path)
        cert.merge_pins(pins, saved['source_sha256'])
        pins[str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
        for power in args.powers:
            if monotonic()-started >= args.seconds:
                receipt['stopped_for_budget'] = True
                break
            geo = geometry(power)
            if max(qset) > geo['groups']:
                raise ValueError('q exceeds group count')
            beginning = monotonic()
            regional = regional_family(macro, max(qset), geo)
            trial = dict(K_power=power, macro_receipt=str(path.resolve()),
                z=str(z), alpha=str(alpha), theta_display=theta, endpoints={}, margins={})
            for q in qset:
                upper = occupancy(regional[q], q, z, alpha, geo)
                bound = probability_endpoint(upper)
                margin = go.display_margin(upper)
                trial['endpoints'][str(q)] = cert.dyadic_record(bound)
                trial['margins'][str(q)] = margin
                if margin > raw_best[power][q]:
                    raw_best[power][q] = margin
                    choices[power][q] = str(path.resolve())
                best[power][q] = min(best[power][q], bound)
            if not args.no_q1 and theta <= .0600000001 and (power, z) not in seen:
                single_support, single_bound = one_group(active, z, geo)
                supports[power] = q1.min_supports(supports[power], single_support)
                trial['q1_support_uppers'] = [cert.dyadic_record(x) for x in single_support]
                trial['q1_margin_display'] = cert.bits(single_bound)
                seen.add((power, z))
            trial['elapsed_seconds'] = monotonic()-beginning
            receipt['trials'].append(trial)
            if not cert.checked_pins(pins):
                raise ArithmeticError('source pins changed during screen')
            print(f'K{power} theta={theta:.5g} alpha={alpha}: '
                  f'{trial["margins"]}; q1={trial.get("q1_margin_display")}; '
                  f'point={trial["elapsed_seconds"]:.2f}s total={monotonic()-started:.2f}s', flush=True)
        if receipt['stopped_for_budget']:
            break
    for power in args.powers:
        geo = geometry(power)
        item = dict(geometry=geo, tested_q=qset,
            selected_endpoints={str(q): cert.dyadic_record(best[power][q]) for q in qset},
            selected_raw_margins_display=raw_best[power], selected_witnesses=choices[power])
        if any(p == power for p, z in seen):
            value = q1.combine_supports(supports[power], groups=geo['groups'])
            item.update(q1_endpoint=cert.dyadic_record(value), q1_margin_display=cert.bits(value),
                q1_selected_support_uppers=[cert.dyadic_record(x) for x in supports[power]])
        receipt['results'][str(power)] = item
    if not cert.checked_pins(pins):
        raise ArithmeticError('final source authentication failed')
    receipt.update(source_pins_verified_at_finish=True, elapsed_seconds=monotonic()-started)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--whole', type=Path, default=OLD/'whole_wider24_v1.json')
    parser.add_argument('--powers', type=int, nargs='+', default=[18,20])
    parser.add_argument('--q', type=int, nargs='+', default=[2,3,4,8,16])
    parser.add_argument('--max-theta', type=float, default=.2)
    parser.add_argument('--seconds', type=float, default=240)
    parser.add_argument('--no-q1', action='store_true')
    parser.add_argument('--extra-macros', type=Path, nargs='*', default=[])
    parser.add_argument('--only-extra', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
