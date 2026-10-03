"""Bounded baseline16 g4 occupancy curve; no map or encoder changes.

Reuse each powered macro family across requested occupancies. Check every
selected winner using a fresh all-logarithmic regional recurrence. Both
evaluations remain floating proposals, not outward-rounded certificates.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, isfinite, log
from pathlib import Path
from time import monotonic

import numpy as np

import fine_grouped_gate as fine

HERE = Path(__file__).resolve().parent
prior, maps, marked = fine.prior, fine.maps, fine.marked

SCHEDULE = {
    16: ('.04', '.06', '.08'),
    64: ('.2', '.25', '.3'),
    2: ('.005', '.01', '.02', '.04'),
    4: ('.005', '.01', '.02', '.04'),
    8: ('.005', '.01', '.02', '.04'),
    32: ('.08', '.12', '.16'),
    96: ('.35', '.4', '.45'),
    160: ('.6', '.7', '.8'),
    256: ('1', '1.2', '1.4'),
    384: ('1.6', '1.8', '2'),
    512: ('2.2',),
}


def source_pins():
    result = fine.source_pins()
    path = Path(__file__).resolve()
    result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def margin_from_region(region, q, theta, alpha, beta_log, cutoff):
    moment = prior.logarithmic.log_power_matrix(region, 32)
    margin = -(log(comb(512, q))+float(alpha)*(q*beta_log+float(theta)*cutoff)+moment)/log(2)
    if not isfinite(margin):
        raise FloatingPointError('nonfinite curve margin')
    return margin, moment


def run(args):
    qs = tuple(args.qs)
    alphas = tuple(Fraction(value) for value in args.alphas)
    if (args.output.exists() or not qs or len(set(qs)) != len(qs)
            or any(q not in SCHEDULE for q in qs) or not alphas
            or len(set(alphas)) != len(alphas) or min(alphas) <= 0 or max(alphas) > 1):
        raise ValueError('fresh output, distinct scheduled q and valid alpha values required')
    started, pins = monotonic(), source_pins()
    data, record = maps.prepare('byte_native')
    shape = prior.geometry(data)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    degrees = {}
    for q in qs:
        for theta in map(Fraction, SCHEDULE[q]):
            degrees[theta] = max(q, degrees.get(theta, 0))
    cache = {}
    result = dict(schema='packet8-g4-curve-proposal-1', map_record=record, geometry=shape,
        construction_changed=False, group_steps=4, proposal_only=True,
        whole_code_certificate=False, has_outward_endpoints=False,
        schedule={str(q): list(SCHEDULE[q]) for q in qs}, alphas=list(map(str, alphas)),
        source_sha256=pins, source_pins_verified_at_finish=False,
        envelope=envelope.metadata(), trials=[], winners={}, failed_tilts={})
    def save():
        if source_pins() != pins:
            raise ArithmeticError('source changed during curve screen')
        result['elapsed_seconds'] = monotonic()-started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    for q in qs:
        candidates = []
        for theta in map(Fraction, SCHEDULE[q]):
            if theta not in cache:
                try:
                    weighted, emission, diagnostic = prior.moments(data, np.exp(-float(theta)))
                    potential = marked.potential_operators(prior.operators(weighted, emission), 0.)
                    tuples = fine.tuple_products(potential, 4)
                    cache[theta] = {}
                    for alpha in alphas:
                        local = fine.fine_operators(tuples, alpha)
                        regional, backend = prior.placement(local, degrees[theta], epochs=16, windows=32)
                        cache[theta][alpha] = (local, regional, backend, diagnostic)
                    del weighted, emission, potential, tuples
                except (ArithmeticError, FloatingPointError) as error:
                    result['failed_tilts'][str(theta)] = str(error)
                    cache[theta] = None
                    save()
                    print(f'theta={theta}: numerical rejection: {error}', flush=True)
            if cache[theta] is None:
                continue
            for alpha in alphas:
                local, regional, backend, diagnostic = cache[theta][alpha]
                margin, moment = margin_from_region(regional[q], q, theta, alpha, beta_log, shape['cutoff'])
                trial = dict(q=q, tilt=str(theta), alpha=str(alpha), margin_bits=margin,
                    log_moment=moment, backend=backend, local_diagnostics=diagnostic)
                candidates.append(trial)
                result['trials'].append(trial)
                print(f'q={q} theta={theta} alpha={alpha}: {margin:.6f} bits', flush=True)
        if candidates:
            winner = dict(max(candidates, key=lambda trial: trial['margin_bits']))
            theta, alpha = Fraction(winner['tilt']), Fraction(winner['alpha'])
            local = cache[theta][alpha][0]
            logs, backend = prior.placement(local, q, epochs=16, windows=32, force_log=True)
            checked, checked_moment = margin_from_region(logs[q], q, theta, alpha, beta_log, shape['cutoff'])
            difference = checked-winner['margin_bits']
            if abs(difference) > 5e-7:
                raise ArithmeticError(f'all-log winner disagrees at q={q}: {difference} bits')
            winner.update(all_log_margin_bits=checked, all_log_difference_bits=difference,
                all_log_backend=backend, powered_macro_operators=local.tolist())
            result['winners'][str(q)] = winner
            print(f'WINNER q={q}: {checked:.6f} bits, theta={theta}, alpha={alpha}; '
                f'all-log delta={difference:.3g}; elapsed={monotonic()-started:.2f}s', flush=True)
        save()
    result['source_pins_verified_at_finish'] = True
    save()
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qs', nargs='+', type=int, default=list(SCHEDULE))
    parser.add_argument('--alphas', nargs='+', default=['.2', '.4', '.6', '1'])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
