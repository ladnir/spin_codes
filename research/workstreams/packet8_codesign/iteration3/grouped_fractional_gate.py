"""Average short runs conditional on total occupancy, then take path powers.

The construction is unchanged. A macro contains g consecutive physical
steps. Given j potential slots in that macro, Q_j is the ordered conditional
average of the already-thinned physical operators. Apply entrywise alpha
to Q_j, not to its summands. All evaluations are floating proposals.
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

import fractional_gate as physical

HERE = Path(__file__).resolve().parent
prior, marked, maps = physical.prior, physical.marked, physical.maps


def grouped_operators(potential, group_steps):
    """E[ordered physical operator product | macro potential count j]."""
    potential = np.asarray(potential, dtype=float)
    if (type(group_steps) is not int or group_steps < 1 or potential.ndim != 3
            or potential.shape[1] != potential.shape[2] or len(potential) < 2
            or not np.isfinite(potential).all() or np.any(potential < 0)):
        raise ValueError('positive square family and positive integer group size required')
    windows = len(potential)-1
    logs, backend = prior.placement(potential, group_steps*windows,
        epochs=group_steps, windows=windows)
    finite = np.isfinite(logs)
    if np.any(logs[finite] < log(np.finfo(float).tiny)):
        raise FloatingPointError('macro coefficients require a log-domain power implementation')
    grouped = np.exp(logs)
    if not np.isfinite(grouped).all() or np.any(grouped < 0):
        raise FloatingPointError('invalid macro operator family')
    return grouped, backend


def powered_operators(grouped, alpha):
    alpha = float(alpha)
    if not isfinite(alpha) or not 0 < alpha <= 1:
        raise ValueError('finite alpha in (0,1] required')
    return np.asarray(grouped, dtype=float)**alpha


def source_pins():
    result = physical.source_pins()
    for path in (Path(__file__).resolve(), HERE/'test_grouped_fractional_gate.py'):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    tilts = tuple(Fraction(x) for x in args.tilts)
    alphas = tuple(Fraction(x) for x in args.alphas)
    groups = tuple(args.groups)
    if (args.output.exists() or not tilts or min(tilts) <= 0 or not alphas
            or min(alphas) <= 0 or max(alphas) > 1 or not groups
            or any(g not in (1, 2, 4) for g in groups)
            or any(len(set(v)) != len(v) for v in (tilts, alphas, groups))
            or not 2 <= args.q <= 512):
        raise ValueError('fresh output, distinct valid grids and q in 2..512 required')
    started, pins = monotonic(), source_pins()
    data, record = maps.prepare('byte_native')
    shape = prior.geometry(data)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    result = dict(schema='packet8-grouped-fractional-proposal-1', geometry=shape,
        map_record=record, construction_changed=False, proposal_only=True,
        whole_code_certificate=False, has_outward_endpoints=False,
        q=args.q, input_tilts=list(map(str, tilts)), fractional_powers=list(map(str, alphas)),
        group_steps=list(groups), alpha1_regression_at_first_tilt=True,
        power_order='potential labels, macro conditional average, then entrywise alpha',
        clipping_condition='chronological macro potential counts; within-macro positions averaged',
        message_factor='C(512,q)*beta^(q*alpha)', subset_union_is_not_powered=True,
        source_sha256=pins, source_pins_verified_at_finish=False,
        envelope=envelope.metadata(), trials=[], best=None, best_choice=None)
    for tilt_index, tilt in enumerate(tilts):
        weighted, emission, diagnostics = prior.moments(data, np.exp(-float(tilt)))
        local = prior.operators(weighted, emission)
        potential = marked.potential_operators(local, 0.)
        reference_regional, _ = prior.placement(potential, args.q, epochs=64, windows=8)
        reference = prior.logarithmic.log_power_matrix(reference_regional[args.q], 32)
        for group in groups:
            grouped, macro_backend = grouped_operators(potential, group)
            powers = alphas+(Fraction(1),) if tilt_index == 0 and 1 not in alphas else alphas
            for alpha in powers:
                powered = powered_operators(grouped, alpha)
                regional, backend = prior.placement(powered, args.q,
                    epochs=64//group, windows=8*group)
                moment = prior.logarithmic.log_power_matrix(regional[args.q], 32)
                regression_error = moment-reference if alpha == 1 else None
                if alpha == 1 and abs(regression_error) > 5e-8:
                    raise ArithmeticError('grouped alpha1 does not regress physical placement')
                exponent = log(comb(512, args.q))+float(alpha)*(
                    args.q*beta_log+float(tilt)*shape['cutoff'])+moment
                margin = -exponent/log(2)
                if not isfinite(margin):
                    raise FloatingPointError('nonfinite grouped fractional result')
                trial = dict(tilt=str(tilt), alpha=str(alpha), group_steps=group,
                    margin_bits=margin, fractional_log_moment=moment,
                    alpha1_physical_log_moment=reference, alpha1_regression_error=regression_error,
                    macro_backend=macro_backend, backend=backend, local_diagnostics=diagnostics)
                result['trials'].append(trial)
                if result['best'] is None or margin > result['best']:
                    result['best'], result['best_choice'] = margin, trial
                if source_pins() != pins:
                    raise ArithmeticError('source changed during grouped proposal')
                result['elapsed_seconds'] = monotonic()-started
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
                print(f'g={group} tilt={tilt} alpha={alpha}: margin={margin:.6f}; '
                    f'best={result["best"]:.6f}; elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q', type=int, default=119)
    parser.add_argument('--tilts', nargs='+', default=['.4', '.5', '.6'])
    parser.add_argument('--alphas', nargs='+', default=['.3', '.4', '.5', '.6'])
    parser.add_argument('--groups', nargs='+', type=int, default=[2, 4])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
