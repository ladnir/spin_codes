"""Retain physical occupancies while grouping hidden comparison-state paths.

For each tuple (j1,...,jg), multiply its already-thinned physical operators,
then take an entrywise fractional power. Average tuples conditional on their
total potential occupancy only AFTER that power. Hypergeometric weights and
the optional H2 exponential marker are not raised to the fractional power.
Floating research proposals only; no encoder or map change.
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

import grouped_fractional_gate as coarse

HERE = Path(__file__).resolve().parent
physical, prior, maps, marked = coarse.physical, coarse.prior, coarse.maps, coarse.marked


def tuple_products(potential, group_steps, *, cap=2):
    """Chronological products and integer multiplicities for fine counts."""
    potential = np.asarray(potential, dtype=float)
    if (type(group_steps) is not int or group_steps not in (1, 2, 4)
            or potential.ndim != 3 or potential.shape[1] != potential.shape[2]
            or len(potential) < 2 or not np.isfinite(potential).all()
            or np.any(potential < 0) or type(cap) is not int or not 1 <= cap < len(potential)):
        raise ValueError('positive family, group size1/2/4 and valid integer cap required')
    W, size = len(potential)-1, potential.shape[1]
    counts = np.arange(W+1)
    costs = np.minimum(counts, cap)
    choices = np.array([comb(W, j) for j in counts], dtype=float)
    matrices = np.eye(size)[None]
    totals, marks, multiplicities = np.zeros(1, dtype=int), np.zeros(1, dtype=int), np.ones(1)
    for _ in range(group_steps):
        matrices = (matrices[:, None]@potential[None]).reshape(-1, size, size)
        totals = (totals[:, None]+counts[None]).reshape(-1)
        marks = (marks[:, None]+costs[None]).reshape(-1)
        multiplicities = (multiplicities[:, None]*choices[None]).reshape(-1)
    denominator = np.array([comb(W*group_steps, int(j)) for j in totals], dtype=float)
    weights = multiplicities/denominator
    for j in range(W*group_steps+1):
        if abs(weights[totals == j].sum()-1) > 3e-14:
            raise ArithmeticError('conditional fine-count weights do not sum to one')
    return dict(products=matrices, totals=totals, costs=marks, weights=weights,
        windows=W*group_steps, group_steps=group_steps, cap=cap)


def fine_operators(tuples, alpha, *, nu=0.):
    alpha, nu = float(alpha), float(nu)
    if not isfinite(alpha) or not 0 < alpha <= 1 or not isfinite(nu) or nu < 0:
        raise ValueError('0<alpha<=1 and finite nonnegative nu required')
    products = tuples['products']
    powered = products**alpha
    # Neither the route multiplicities nor the H_c marker belong inside alpha.
    weighted = (tuples['weights']*np.exp(nu*tuples['costs']))[:, None, None]*powered
    local = np.zeros((tuples['windows']+1, products.shape[1], products.shape[2]))
    np.add.at(local, tuples['totals'], weighted)
    if not np.isfinite(local).all() or np.any(local < 0):
        raise FloatingPointError('invalid fine-grouped operator')
    return local


def source_pins():
    result = coarse.source_pins()
    for path in (Path(__file__).resolve(), HERE/'test_fine_grouped_gate.py'):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    tilts = tuple(Fraction(x) for x in args.tilts)
    alphas = tuple(Fraction(x) for x in args.alphas)
    nus = tuple(Fraction(x) for x in args.nus)
    groups = tuple(args.groups)
    if (args.output.exists() or not tilts or min(tilts) <= 0 or not alphas
            or min(alphas) <= 0 or max(alphas) > 1 or not nus or min(nus) < 0
            or not groups or any(g not in (1, 2, 4) for g in groups)
            or any(len(set(v)) != len(v) for v in (tilts, alphas, nus, groups))
            or not 2 <= args.q <= 512):
        raise ValueError('fresh output, distinct valid grids and q in2..512 required')
    started, pins = monotonic(), source_pins()
    data, record = maps.prepare('byte_native')
    shape = prior.geometry(data)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    witness = None
    if max(nus) > 0:
        witness = physical.cap_counts.rational_witness(args.q, args.h,
            args.route_x, cap=2, target_bits=60)
        if not witness['exact_integer_check_passed']:
            raise ArithmeticError('route threshold fails its exact60-bit witness')
    result = dict(schema='packet8-fine-grouped-fractional-proposal-1', geometry=shape,
        map_record=record, construction_changed=False, proposal_only=True,
        whole_code_certificate=False, has_outward_endpoints=False,
        q=args.q, input_tilts=list(map(str, tilts)), fractional_powers=list(map(str, alphas)),
        group_steps=list(groups), route_tilts=list(map(str, nus)), route_witness=witness,
        alpha1_regression_at_first_tilt=True,
        power_order='potential labels, ordered tuple product, entrywise alpha, tuple average',
        route_marker_order='exp(nu*sum_i min(j_i,2)) after power; weights never powered',
        clipping_condition='physical potential counts; hidden state paths grouped before power',
        message_factor='C(512,q)*beta^(q*alpha)', subset_union_is_not_powered=True,
        source_sha256=pins, source_pins_verified_at_finish=False,
        envelope=envelope.metadata(), trials=[], best=None, best_choice=None)
    for tilt_index, tilt in enumerate(tilts):
        weighted, emission, diagnostics = prior.moments(data, np.exp(-float(tilt)))
        physical_local = prior.operators(weighted, emission)
        potential = marked.potential_operators(physical_local, 0.)
        for group in groups:
            tuples = tuple_products(potential, group)
            powers = alphas+(Fraction(1),) if tilt_index == 0 and 1 not in alphas else alphas
            for nu in nus:
                physical_marked = potential*np.exp(float(nu)*np.minimum(np.arange(9), 2))[:, None, None]
                reference_regional, _ = prior.placement(physical_marked, args.q, epochs=64, windows=8)
                reference = prior.logarithmic.log_power_matrix(reference_regional[args.q], 32)
                for alpha in powers:
                    local = fine_operators(tuples, alpha, nu=nu)
                    regional, backend = prior.placement(local, args.q,
                        epochs=64//group, windows=8*group)
                    moment = prior.logarithmic.log_power_matrix(regional[args.q], 32)
                    regression_error = moment-reference if alpha == 1 else None
                    if alpha == 1 and abs(regression_error) > 5e-8:
                        raise ArithmeticError('fine-grouped alpha1 does not regress marked physical placement')
                    exponent = (log(comb(512, args.q))+float(alpha)*(
                        args.q*beta_log+float(tilt)*shape['cutoff'])-float(nu)*args.h+moment)
                    combined = exponent
                    if nu > 0:
                        combined = float(np.logaddexp(exponent,
                            -witness['certified_margin_bits_floor']*log(2)))
                    margin = -combined/log(2)
                    if not isfinite(margin):
                        raise FloatingPointError('nonfinite fine-grouped result')
                    trial = dict(tilt=str(tilt), alpha=str(alpha), group_steps=group,
                        route_tilt=str(nu), good_message_margin_bits=-exponent/log(2),
                        margin_bits=margin, fractional_log_moment=moment,
                        alpha1_physical_log_moment=reference, alpha1_regression_error=regression_error,
                        backend=backend, local_diagnostics=diagnostics)
                    result['trials'].append(trial)
                    if result['best'] is None or margin > result['best']:
                        result['best'], result['best_choice'] = margin, trial
                    if source_pins() != pins:
                        raise ArithmeticError('source changed during fine-grouped proposal')
                    result['elapsed_seconds'] = monotonic()-started
                    args.output.parent.mkdir(parents=True, exist_ok=True)
                    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
                    print(f'g={group} tilt={tilt} alpha={alpha} nu={nu}: margin={margin:.6f}; '
                        f'best={result["best"]:.6f}; elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q', type=int, default=119)
    parser.add_argument('--tilts', nargs='+', default=['.4', '.5', '.6'])
    parser.add_argument('--alphas', nargs='+', default=['.3', '.4', '.5'])
    parser.add_argument('--groups', nargs='+', type=int, default=[2])
    parser.add_argument('--nus', nargs='+', default=['0'])
    parser.add_argument('--h', type=int, default=2559)
    parser.add_argument('--route-x', default='9/35')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
