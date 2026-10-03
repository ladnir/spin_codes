"""One bounded H2 route-conditioning gate; no new encoder construction."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import numpy as np

import gate
import cap_counts


def marked(local, nu, cap=2, *, packet_bits=8):
    if type(cap) is not int or not 1 <= cap < len(local) or nu < 0:
        raise ValueError('positive packet cap and nonnegative tilt required')
    result = gate.potential_operators(local, 0., packet_bits=packet_bits)
    result *= np.exp(float(nu)*np.minimum(np.arange(len(local)), cap))[:, None, None]
    return result


def sources():
    result = gate.source_pins()
    for path in (Path(__file__).resolve(), Path(cap_counts.__file__).resolve(), gate.HERE/'test_cap_gate.py'):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    if args.output.exists():
        raise ValueError('a fresh output is required')
    tilts, nus = tuple(map(Fraction, args.tilts)), tuple(map(Fraction, args.nus))
    if min(tilts) <= 0 or min(nus) < 0:
        raise ValueError('positive weight tilts and nonnegative route tilts required')
    q, h, x = args.q, args.h, Fraction(args.chernoff_x)
    witness = cap_counts.rational_witness(q, h, x, cap=2, target_bits=60)
    if not witness['exact_integer_check_passed']:
        raise ValueError('the requested route threshold must have an exact60-bit witness')
    bad_margin = witness['certified_margin_bits_floor']
    started, pins = monotonic(), sources()
    name = 'byte_native_A_scaled' if args.scaled else 'byte_native'
    data, record = gate.maps.prepare(name)
    beta = gate.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
    beta_log = log(beta.numerator)-log(beta.denominator)
    result = dict(schema='packet8-cap2-conditioned-proposal-1', map_record=record,
        unchanged_baseline_construction=not args.scaled, q=q, h=h, cap=2,
        geometry=gate.prior.geometry(data), bad_route_witness=witness,
        bad_route_bound_exact=True, good_message_bound_outward=False,
        proposal_only=True, whole_code_certificate=False, source_sha256=pins,
        source_pins_verified_at_finish=False, input_tilts=list(map(str, tilts)),
        route_tilts=list(map(str, nus)), trials=[], best_margin_bits=-float('inf'))
    for tilt in tilts:
        weighted, emission, _ = gate.prior.moments(data, np.exp(-float(tilt)))
        local = gate.prior.operators(weighted, emission)
        for nu in nus:
            local_marked = marked(local, float(nu))
            regional, backend = gate.prior.placement(local_marked, q, epochs=64, windows=8)
            moment = gate.prior.logarithmic.log_power_matrix(regional[q], 32)
            exponent = log(comb(512, q))+q*beta_log+float(tilt)*13107-float(nu)*h+moment
            good_margin = -exponent/log(2)
            margin = -float(np.logaddexp(exponent, -bad_margin*log(2)))/log(2)
            trial = dict(tilt=str(tilt), route_tilt=str(nu), backend=backend,
                marked_log_moment=moment, good_message_margin_bits=good_margin,
                combined_margin_bits=margin)
            result['trials'].append(trial)
            if margin > result['best_margin_bits']:
                result['best_margin_bits'], result['best_choice'] = margin, trial
            if pins != sources():
                raise ArithmeticError('source changed during proposal')
            result['elapsed_seconds'] = monotonic()-started
            args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
            print(f'H2 {name} theta={tilt} nu={nu}: good={good_margin:.6f}; '
                  f'best={result["best_margin_bits"]:.6f}; elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q', type=int, default=119)
    parser.add_argument('--h', type=int, default=2559)
    parser.add_argument('--chernoff-x', default='9/35')
    parser.add_argument('--scaled', action='store_true')
    parser.add_argument('--tilts', nargs='+', default=['.3', '.4', '.5'])
    parser.add_argument('--nus', nargs='+', default=['.5', '1', '1.5', '2', '3'])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
