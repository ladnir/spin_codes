"""Fractional first moments over routing geometry; floating proposals only.

Condition on the chronological potential-occupancy sequence J, not on exact
positions. Average within-step positions and uniform byte labels first.
For 0<alpha<=1, min(1,B F(J)) <= B^alpha F(J)^alpha. Subadditivity bounds
F(J)^alpha by the product of entrywise-powered, already-thinned operators.
The group-subset union C(512,q) remains outside the fractional power.
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
BASE = HERE.parent
sys.path.insert(0, str(BASE))
import maps
import screen as prior
sys.path.insert(0, str(BASE/'route_conditioning'))
import gate as marked
import cap_counts


def fractional_operators(local, alpha, *, nu=0., cap=2, packet_bits=8):
    """Thin potential slots, then apply the pathwise fractional majorant."""
    alpha, nu = float(alpha), float(nu)
    if (not isfinite(alpha) or not 0 < alpha <= 1 or not isfinite(nu) or nu < 0
            or type(cap) is not int or not 1 <= cap < len(local)):
        raise ValueError('0<alpha<=1, nonnegative finite nu and a valid integer cap required')
    potential = marked.potential_operators(local, 0., packet_bits=packet_bits)
    powered = potential**alpha
    costs = np.minimum(np.arange(len(local)), cap)
    powered *= np.exp(nu*costs)[:, None, None]
    if not np.isfinite(powered).all() or np.any(powered < 0):
        raise FloatingPointError('nonnegative finite fractional operators required')
    return powered


def source_pins():
    result = prior.sources('byte_native')
    for path in (Path(__file__).resolve(), HERE/'test_fractional_gate.py',
            Path(marked.__file__).resolve(), Path(cap_counts.__file__).resolve()):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    tilts = tuple(Fraction(x) for x in args.tilts)
    alphas = tuple(Fraction(x) for x in args.alphas)
    nus = tuple(Fraction(x) for x in args.nus)
    if (args.output.exists() or not tilts or min(tilts) <= 0 or not alphas
            or min(alphas) <= 0 or max(alphas) > 1 or not nus or min(nus) < 0
            or any(len(set(values)) != len(values) for values in (tilts, alphas, nus))
            or not 2 <= args.q <= 512):
        raise ValueError('fresh output, distinct valid tilt grids and q in2..512 required')
    started, pins = monotonic(), source_pins()
    data, record = maps.prepare('byte_native')
    shape = prior.geometry(data)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    route_witness = None
    if max(nus) > 0:
        route_witness = cap_counts.rational_witness(args.q, args.h, args.route_x,
            cap=args.cap, target_bits=60)
        if not route_witness['exact_integer_check_passed']:
            raise ValueError('the requested optional route threshold fails its exact60-bit gate')
    result = dict(schema='packet8-fractional-geometry-proposal-1', geometry=shape,
        map_record=record, construction_changed=False, proposal_only=True,
        whole_code_certificate=False, has_outward_endpoints=False,
        q=args.q, input_tilts=list(map(str, tilts)), fractional_powers=list(map(str, alphas)),
        route_tilts=list(map(str, nus)), route_witness=route_witness,
        power_order='potential-label mixture first; entrywise alpha power second',
        clipping_condition='chronological potential occupancies, with positions averaged conditionally',
        message_factor='C(512,q)*beta^(q*alpha)', subset_union_is_not_powered=True,
        source_sha256=pins, source_pins_verified_at_finish=False,
        envelope=envelope.metadata(), trials=[], best=None, best_choice=None)
    for tilt in tilts:
        weighted, emission, diagnostics = prior.moments(data, np.exp(-float(tilt)))
        local = prior.operators(weighted, emission)
        potential = marked.potential_operators(local, 0.)
        reference_regional, _ = prior.placement(potential, args.q, epochs=64, windows=8)
        reference = prior.logarithmic.log_power_matrix(reference_regional[args.q], 32)
        for alpha in alphas:
            for nu in nus:
                powered = fractional_operators(local, alpha, nu=nu, cap=args.cap)
                regional, backend = prior.placement(powered, args.q, epochs=64, windows=8)
                moment = prior.logarithmic.log_power_matrix(regional[args.q], 32)
                if alpha == 1 and nu == 0 and abs(moment-reference) > 2e-8:
                    raise ArithmeticError('alpha1 does not reproduce exact-q potential placement')
                exponent = (log(comb(512, args.q))
                    +float(alpha)*(args.q*beta_log+float(tilt)*shape['cutoff'])
                    -float(nu)*args.h+moment)
                combined = exponent
                if nu > 0:
                    route_bits = route_witness['certified_margin_bits_floor']
                    combined = float(np.logaddexp(-route_bits*log(2), exponent))
                margin = -combined/log(2)
                if not isfinite(margin):
                    raise FloatingPointError('nonfinite fractional gate result')
                trial = dict(tilt=str(tilt), alpha=str(alpha), route_tilt=str(nu),
                    good_message_margin_bits=-exponent/log(2), combined_margin_bits=margin,
                    fractional_log_moment=moment, alpha1_unmarked_log_moment=reference,
                    backend=backend, local_diagnostics=diagnostics)
                result['trials'].append(trial)
                if result['best'] is None or margin > result['best']:
                    result['best'], result['best_choice'] = margin, trial
                if source_pins() != pins:
                    raise ArithmeticError('source changed during proposal')
                result['elapsed_seconds'] = monotonic()-started
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
                print(f'tilt={tilt} alpha={alpha} nu={nu}: margin={margin:.6f}; '
                    f'best={result["best"]:.6f}; elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q', type=int, default=119)
    parser.add_argument('--tilts', nargs='+', default=['.3', '.4', '.5'])
    parser.add_argument('--alphas', nargs='+', default=['.25', '.4', '.55', '.7', '.85', '1'])
    parser.add_argument('--nus', nargs='+', default=['0'])
    parser.add_argument('--cap', type=int, default=2)
    parser.add_argument('--h', type=int, default=2559)
    parser.add_argument('--route-x', default='9/35')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
