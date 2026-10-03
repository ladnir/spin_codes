"""Floating route-conditioned gate for the unchanged byte-native encoder.

For q selected outer groups, H counts physical steps containing at least one
POTENTIAL packet from these groups, regardless of its eventual byte label.
The setup event requires H>=h for every q-group subset. Its failure union
has C(G,q) terms, not C(G,q)*beta^q. On good routes, exp(nu*(H-h)) bounds
the indicator for nu>=0. The remaining message comparison uses beta^q.
No route is reset or changed, and no claim is an outward certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, isfinite, log
from pathlib import Path
import sys
from time import monotonic

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import maps
import screen as prior
import route_counts


def potential_operators(local, nu, *, packet_bits=8):
    """Thin j potential slots, then mark whether their physical step is hit."""
    local = np.asarray(local, dtype=float)
    if local.ndim != 3 or local.shape[1] != local.shape[2] or nu < 0:
        raise ValueError('positive square local family and nonnegative route tilt required')
    W = len(local)-1
    p = 1-2.**(-packet_bits)
    mixed = np.empty_like(local)
    for j in range(W+1):
        mixed[j] = sum(comb(j, k)*p**k*(1-p)**(j-k)*local[k] for k in range(j+1))
        if j:
            mixed[j] *= np.exp(float(nu))
    if not np.isfinite(mixed).all() or np.any(mixed < 0):
        raise FloatingPointError('invalid marked potential-slot family')
    return mixed


def threshold(value):
    q, h, eta = value.split(':')
    q, h, eta = int(q), int(h), Fraction(eta)
    if not 2 <= q <= 512 or not 32*((q+7)//8) <= h <= 32*min(64, q) or eta <= 0:
        raise ValueError('feasible q:h:positive_eta required')
    return q, h, eta


def source_pins():
    result = prior.sources('byte_native')
    for path in (Path(__file__).resolve(), HERE/'test_gate.py', Path(route_counts.__file__).resolve()):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    tilts = tuple(Fraction(x) for x in args.tilts)
    nus = tuple(Fraction(x) for x in args.nus)
    thresholds = tuple(threshold(x) for x in args.thresholds)
    if (args.output.exists() or not tilts or min(tilts) <= 0 or not nus or min(nus) < 0
            or len(set(tilts)) != len(tilts) or len(set(nus)) != len(nus)
            or len(set(q for q, _, _ in thresholds)) != len(thresholds)):
        raise ValueError('fresh output, distinct tilts and distinct occupancy thresholds required')
    started, pins = monotonic(), source_pins()
    data, map_record = maps.prepare('byte_native')
    shape = prior.geometry(data)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    bad = {q: route_counts.bad_route_margin(q, h, float(eta), regions=32) for q, h, eta in thresholds}
    result = dict(schema='packet8-route-conditioned-proposal-1', geometry=shape,
        map_record=map_record, construction_changed=False, proposal_only=True,
        whole_code_certificate=False, has_outward_endpoints=False,
        source_sha256=pins, source_pins_verified_at_finish=False,
        envelope=envelope.metadata(), input_tilts=list(map(str, tilts)),
        route_tilts=list(map(str, nus)), thresholds=[dict(q=q, h=h, eta=str(eta),
            bad_route_margin_bits=bad[q]) for q, h, eta in thresholds],
        bad_route_union_factor='C(512,q), with no beta^q',
        good_message_union_factor='C(512,q)*beta^q',
        q1_checked=False, all_occupancies_checked=False, trials=[], best={}, choices={})
    maximum = max(q for q, _, _ in thresholds)
    for tilt in tilts:
        weighted, emission, local_diagnostic = prior.moments(data, np.exp(-float(tilt)))
        local = prior.operators(weighted, emission)
        old_regional, _ = prior.placement(local, maximum, epochs=64, windows=8)
        reference = {q: prior.logarithmic.log_power_matrix(prior.uniform_mixture(old_regional, q), 32)
                     for q, _, _ in thresholds}
        for nu in nus:
            potential = potential_operators(local, float(nu))
            regional, backend = prior.placement(potential, maximum, epochs=64, windows=8)
            trial = dict(tilt=str(tilt), route_tilt=str(nu), backend=backend,
                local_diagnostics=local_diagnostic, witnesses={})
            for q, h, eta in thresholds:
                moment = prior.logarithmic.log_power_matrix(regional[q], 32)
                if nu == 0 and abs(moment-reference[q]) > 2e-8:
                    raise ArithmeticError('potential-slot thinning disagrees with prior active-slot mixture')
                exponent = log(comb(512, q))+q*beta_log+float(tilt)*13107-float(nu)*h+moment
                good_margin = -exponent/log(2)
                total_margin = -float(np.logaddexp(-bad[q]*log(2), exponent))/log(2)
                if not isfinite(good_margin) or not isfinite(total_margin):
                    raise FloatingPointError('invalid route-conditioned margin')
                trial['witnesses'][str(q)] = dict(good_message_margin_bits=good_margin,
                    bad_route_margin_bits=bad[q], combined_margin_bits=total_margin,
                    marked_log_moment=moment, old_log_moment=reference[q],
                    improvement_over_unconditioned_bits=(reference[q]-moment+float(nu)*h)/log(2))
                if str(q) not in result['best'] or total_margin > result['best'][str(q)]:
                    result['best'][str(q)] = total_margin
                    result['choices'][str(q)] = dict(tilt=str(tilt), route_tilt=str(nu),
                        good_message_margin_bits=good_margin, bad_route_margin_bits=bad[q])
            result['trials'].append(trial)
            if source_pins() != pins:
                raise ArithmeticError('source changed during proposal')
            result['elapsed_seconds'] = monotonic()-started
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
            print(f'tilt={tilt} nu={nu}: best={result["best"]}; '
                  f'elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--thresholds', nargs='+', default=[
        '64:1057:1.5993122190268654', '119:1483:1.8952931340202324', '128:1531:1.9687464380734174'])
    parser.add_argument('--tilts', nargs='+', default=['.3', '.4', '.5', '.6', '.8', '1.0'])
    parser.add_argument('--nus', nargs='+', default=['0', '.5', '1', '1.5', '2', '2.5', '3', '4'])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
