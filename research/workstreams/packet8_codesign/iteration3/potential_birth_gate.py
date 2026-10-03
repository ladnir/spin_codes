"""Rebase comparison birth families after potential-label thinning.

The positive comparison moment is unchanged for every physical occupancy
sequence. A common row-stochastic change-of-basis matrix Q satisfies
Tnew_j Q = Q Told_j, e0 Q=e0, and Q 1=1. Fractional path bounds can improve
because a potential-j birth occupies one coordinate instead of a mixture.
This does not identify the filled comparison closure with the physical law.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import numpy as np

import fine_grouped_gate as fine

HERE = Path(__file__).resolve().parent
prior, maps, marked = fine.prior, fine.maps, fine.marked

POINTS = {16: ('.06', '.4'), 64: ('.25', '.4'), 119: ('.475', '.35')}


def rebase_potential(potential):
    """Replace each potential-j birth mixture by its own normalized family."""
    old = np.asarray(potential, dtype=float)
    if (old.ndim != 3 or old.shape[1] != old.shape[2] or old.shape[1] != len(old)+1
            or not np.isfinite(old).all() or np.any(old < 0)
            or np.any(old[:, 1:, 2:] != 0) or np.any(old[:, 0, 1] != 0)
            or np.any(old[0, 0, 2:] != 0)):
        raise ValueError('finite positive delta0/U/physical-birth closure required')
    W, size = len(old)-1, old.shape[1]
    mass = old[1:, 0, 2:].sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError('each nonempty potential occupancy needs positive birth mass')
    change = np.eye(size)
    change[2:, :] = 0
    change[2:, 2:] = old[1:, 0, 2:]/mass[:, None]
    if not np.allclose(change.sum(axis=1), 1., rtol=0, atol=3e-15):
        raise ArithmeticError('birth change of basis is not row stochastic')
    transformed = change@old
    new = np.zeros_like(old)
    new[:, :, :2] = transformed[:, :, :2]
    for j in range(1, W+1):
        new[j, 0, j+1] = mass[j-1]
    left, right = new@change, change@old
    positive = right > 0
    if (np.any(left[~positive] != 0)
            or not np.allclose(left[positive], right[positive], rtol=2e-13, atol=0)):
        raise ArithmeticError('comparison-family intertwining failed')
    relative_error = float(np.max(np.abs(left[positive]/right[positive]-1)))
    return new, change, dict(intertwining_max_relative_error=relative_error,
        stochastic_max_absolute_error=float(np.max(np.abs(change.sum(axis=1)-1))),
        comparison_moment_preserved=True, physical_kernel_not_claimed_exact=True)


def source_pins():
    result = fine.source_pins()
    for path in (Path(__file__).resolve(), HERE/'test_potential_birth_gate.py'):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def evaluate(local, q, theta, alpha, beta_log, cutoff, *, group=4, force_log=False):
    regional, backend = prior.placement(local, q, epochs=64//group,
        windows=8*group, force_log=force_log)
    moment = prior.logarithmic.log_power_matrix(regional[q], 32)
    margin = -(log(comb(512, q))+float(alpha)*(q*beta_log+float(theta)*cutoff)+moment)/log(2)
    return margin, moment, backend


def run(args):
    if args.output.exists() or not args.qs or len(set(args.qs)) != len(args.qs) or any(q not in POINTS for q in args.qs):
        raise ValueError('fresh output and distinct supported points required')
    started, pins = monotonic(), source_pins()
    data, record = maps.prepare('byte_native')
    shape = prior.geometry(data)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    result = dict(schema='packet8-potential-birth-basis-proposal-1', map_record=record,
        geometry=shape, construction_changed=False, group_steps=4,
        proposal_only=True, whole_code_certificate=False, has_outward_endpoints=False,
        source_sha256=pins, source_pins_verified_at_finish=False,
        moment_identity='Tnew_j Q=Q Told_j; e0 Q=e0; Q1=1',
        comparison_not_physical_kernel=True, envelope=envelope.metadata(), trials=[])
    for q in args.qs:
        theta, alpha = map(Fraction, POINTS[q])
        weighted, emission, diagnostic = prior.moments(data, np.exp(-float(theta)))
        old = marked.potential_operators(prior.operators(weighted, emission), 0.)
        new, change, algebra = rebase_potential(old)
        old_group = fine.fine_operators(fine.tuple_products(old, 4), alpha)
        new_group = fine.fine_operators(fine.tuple_products(new, 4), alpha)
        old_margin, old_moment, old_backend = evaluate(old_group, q, theta, alpha, beta_log, shape['cutoff'])
        new_margin, new_moment, new_backend = evaluate(new_group, q, theta, alpha, beta_log, shape['cutoff'])
        checked, _, check_backend = evaluate(new_group, q, theta, alpha, beta_log,
            shape['cutoff'], force_log=True)
        if abs(checked-new_margin) > 5e-7:
            raise ArithmeticError('fresh all-log check failed')
        old_alpha1 = evaluate(old, q, theta, 1, beta_log, shape['cutoff'], group=1)
        new_alpha1 = evaluate(new, q, theta, 1, beta_log, shape['cutoff'], group=1)
        if abs(old_alpha1[1]-new_alpha1[1]) > 5e-8:
            raise ArithmeticError('alpha1 global comparison moment changed')
        trial = dict(q=q, tilt=str(theta), alpha=str(alpha), old_margin_bits=old_margin,
            new_margin_bits=new_margin, gain_bits=new_margin-old_margin,
            old_backend=old_backend, new_backend=new_backend, check_backend=check_backend,
            all_log_margin_bits=checked, all_log_difference_bits=checked-new_margin,
            alpha1_log_moment_difference=new_alpha1[1]-old_alpha1[1],
            change_of_basis=change.tolist(), algebra_checks=algebra,
            powered_macro_operators=new_group.tolist(), local_diagnostics=diagnostic)
        result['trials'].append(trial)
        if source_pins() != pins:
            raise ArithmeticError('source changed during potential-birth gate')
        result['elapsed_seconds'] = monotonic()-started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(f'q={q} theta={theta} alpha={alpha}: old={old_margin:.6f}, '
            f'new={new_margin:.6f}, gain={new_margin-old_margin:.6f}; '
            f'elapsed={monotonic()-started:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qs', nargs='+', type=int, default=[16, 64, 119])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
