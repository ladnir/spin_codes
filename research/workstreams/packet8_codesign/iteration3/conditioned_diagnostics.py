"""H2-tilted comparison diagnostics: potential slots, active bytes, and state.

The probability measure below is the normalized positive comparison path
sum, including its H2 exponential marker. It is not the actual encoder's
failure distribution, nor the exact law conditioned on H2>=2559.
"""
from __future__ import annotations

import os
for _variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_variable] = '1'

import argparse
import hashlib
import json
from math import comb, log
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from scipy.special import gammaln, logsumexp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'route_conditioning'))
import cap_gate


def batch_moments(families, q, *, epochs=64, windows=8, regions=32):
    """Batch independent floating perturbations of the ordered recurrence.

    Each occupancy and each perturbation has its own scale. Subnormal
    contributions may underflow; no probability law is changed. Every real
    diagnostic checks baseline and selected perturbations against the frozen
    all-logarithmic evaluator before reporting finite-difference counts.
    """
    families = np.asarray(families, dtype=float)
    B, terms, size, other = families.shape
    if (terms != windows+1 or size != other or q > epochs*windows
            or np.any(families < 0) or not np.isfinite(families).all()):
        raise ValueError('positive square operator families and feasible geometry required')
    scales = families.max(axis=(2, 3))
    if np.any(scales <= 0):
        raise ValueError('positive operator for every occupancy required')
    with np.errstate(under='ignore', over='raise', divide='raise', invalid='raise'):
        normalized = families/scales[:, :, None, None]
        logs = np.log(scales)
        current = np.broadcast_to(np.eye(size), (B, 1, size, size)).copy()
        current_logs = np.zeros((B, 1))
        for epoch in range(1, epochs+1):
            limit, total, previous = min(q, epoch*windows), epoch*windows, (epoch-1)*windows
            contributions, base = [], np.full((B, limit+1), -np.inf)
            for k in range(min(windows, limit)+1):
                js = np.arange(k, min(limit, k+current.shape[1]-1)+1)
                old = js-k
                weight = (log(comb(windows, k))+gammaln(previous+1)-gammaln(old+1)
                    -gammaln(previous-old+1)-gammaln(total+1)+gammaln(js+1)+gammaln(total-js+1))
                exponent = weight[None, :]+current_logs[:, old]+logs[:, k, None]
                base[:, js] = np.maximum(base[:, js], exponent)
                contributions.append((k, js, old, exponent))
            following = np.zeros((B, limit+1, size, size))
            for k, js, old, exponent in contributions:
                following[:, js] += np.exp(exponent-base[:, js])[:, :, None, None]*(
                    current[:, old] @ normalized[:, k, None])
            scale = following.max(axis=(2, 3))
            if np.any(scale <= 0) or not np.isfinite(scale).all():
                raise FloatingPointError('nonpositive batched placement scale')
            current, current_logs = following/scale[:, :, None, None], base+np.log(scale)
    matrix = current[:, q]
    logged = np.full_like(matrix, -np.inf)
    np.log(matrix, where=matrix > 0, out=logged)
    logged += current_logs[:, q, None, None]
    row = np.full((B, size), -np.inf)
    row[:, 0] = 0.
    for _ in range(regions):
        row = logsumexp(row[:, :, None]+logged, axis=1)
    totals = logsumexp(row, axis=1)
    return totals, np.exp(row-totals[:, None])


def make_perturbations(local, nu, *, epsilon=1e-4, packet_bits=8):
    """Scale active operators before thinning; potential operators afterward."""
    base = cap_gate.marked(local, nu, cap=2, packet_bits=packet_bits)
    W = len(local)-1
    descriptions = []
    for j in range(W+1):
        descriptions.append((f'potential_{j}', 'potential', (j, slice(None), slice(None))))
    for k in range(W+1):
        descriptions.append((f'active_{k}', 'active', (k, slice(None), slice(None))))
    descriptions.extend([
        ('zero_source', 'active', (slice(None), 0, slice(None))),
        ('nonzero_source', 'active', (slice(None), slice(1, None), slice(None))),
        ('uniform_source', 'active', (slice(None), 1, slice(None))),
        ('birth_source', 'active', (slice(None), slice(2, None), slice(None))),
        ('empty_zero_to_zero', 'active', (0, 0, 0)),
        ('active_zero_feedback_zero_to_zero', 'active', (slice(1, None), 0, 0)),
        ('birth_transition', 'active', (slice(1, None), 0, slice(2, None))),
        ('return_transition', 'active', (slice(1, None), slice(1, None), 0)),
    ])
    families = [base]
    for _, domain, selector in descriptions:
        for sign in (-1, 1):
            changed = (base if domain == 'potential' else local).copy()
            changed[selector] *= np.exp(sign*epsilon)
            if domain == 'active':
                changed = cap_gate.marked(changed, nu, cap=2, packet_bits=packet_bits)
            families.append(changed)
    return np.asarray(families), [name for name, _, _ in descriptions]


def logarithmic_moment(local, q, *, epochs=64, windows=8, regions=32):
    regional, _ = cap_gate.gate.prior.placement(local, q, epochs=epochs, windows=windows, force_log=True)
    return cap_gate.gate.prior.logarithmic.log_power_matrix(regional[q], regions)


def checked_pins(pins):
    return all(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest for path, digest in pins.items())


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    pins = cap_gate.sources()
    for path in (Path(__file__).resolve(), HERE/'test_conditioned_diagnostics.py'):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    results = dict(schema='H2-conditioned-comparison-diagnostics-1',
        proposal_only=True, whole_code_certificate=False, actual_failure_distribution=False,
        exact_good_event_conditional_distribution=False, source_sha256=pins,
        source_pins_verified_at_finish=False, epsilon=args.epsilon,
        geometry=dict(q=119, regions=32, slots_per_region=512, physical_steps_per_region=64,
            windows=8, total_physical_steps=2048, zero_initial_state=True, final_flush=False,
            continuous_state=True), route_marker=dict(cap=2, nu=3.5, h=2559), cases={})
    started = monotonic()
    for name in args.cases:
        case_start = monotonic()
        if name == 'baseline16':
            data, record = cap_gate.gate.maps.prepare('byte_native')
            theta = .4
            weighted, emission, _ = cap_gate.gate.prior.moments(data, np.exp(-theta))
            local = cap_gate.gate.prior.operators(weighted, emission)
        elif name == 'actual24':
            source = args.receipt.resolve()
            receipt = json.loads(source.read_text(encoding='utf-8'))
            if not receipt['source_pins_verified_at_finish'] or not checked_pins(receipt['source_sha256']):
                raise ValueError('complete authenticated24-bit source required')
            pins.update(receipt['source_sha256'])
            pins[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
            local, record, theta = np.asarray(receipt['local_operator_matrices']), receipt['map_record'], receipt['tilt']
            if record['state_bits'] != 24 or theta != .5 or local.shape != (9, 10, 10):
                raise ValueError('the specified actual24-bit point is required')
        else:
            raise ValueError('unknown case')
        families, labels = make_perturbations(local, 3.5, epsilon=args.epsilon)
        moments, terminals = batch_moments(families, 119)
        estimates = (moments[2::2]-moments[1::2])/(2*args.epsilon)
        counts = dict(zip(labels, map(float, estimates)))
        potential = [counts[f'potential_{j}'] for j in range(9)]
        active = [counts[f'active_{k}'] for k in range(9)]
        # Check selected positive and negative perturbations independently in
        # log space, including active thinning and the rare zero-feedback cell.
        checks = [0, 1, 2, 1+2*labels.index('active_0'), 2+2*labels.index('active_0'),
                  1+2*labels.index('active_zero_feedback_zero_to_zero'),
                  2+2*labels.index('active_zero_feedback_zero_to_zero')]
        errors = {}
        for index in checks:
            expected = logarithmic_moment(families[index], 119)
            errors[str(index)] = float(moments[index]-expected)
            if abs(errors[str(index)]) > 3e-8:
                raise ArithmeticError('batched scaled diagnostic differs from log evaluator')
        residuals = dict(potential_steps=sum(potential)-2048,
            active_steps=sum(active)-2048, potential_packets=sum(j*x for j, x in enumerate(potential))-32*119,
            state_partition=counts['zero_source']+counts['nonzero_source']-2048,
            nonzero_partition=counts['uniform_source']+counts['birth_source']-counts['nonzero_source'],
            zero_partition=counts['empty_zero_to_zero']+counts['active_zero_feedback_zero_to_zero']
                +counts['birth_transition']-counts['zero_source'],
            state_flow=counts['birth_transition']-counts['return_transition']-(1-float(terminals[0, 0])))
        if max(map(abs, residuals.values())) > 5e-4:
            raise ArithmeticError('diagnostic conservation checks failed')
        beta = cap_gate.gate.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
        fixed = log(comb(512, 119))+119*(log(beta.numerator)-log(beta.denominator))+theta*13107-3.5*2559
        results['cases'][name] = dict(theta=theta, map_record=record, no_large_state_census=True,
            marked_log_moment=float(moments[0]), good_message_margin_bits=-(fixed+float(moments[0]))/log(2),
            potential_step_histogram=potential, active_step_histogram=active,
            expected_H1=2048-potential[0], expected_H2=sum(min(j, 2)*x for j, x in enumerate(potential)),
            expected_active_packets=sum(k*x for k, x in enumerate(active)),
            expected_potential_packets=sum(j*x for j, x in enumerate(potential)),
            state_and_transition_counts={key: value for key, value in counts.items()
                if not key.startswith(('potential_', 'active_')) or key == 'active_zero_feedback_zero_to_zero'},
            terminal_probabilities=terminals[0].tolist(), conservation_residuals=residuals,
            scaled_vs_log_checks=errors, elapsed_seconds=monotonic()-case_start)
        if not checked_pins(pins):
            raise ArithmeticError('source changed during diagnostic')
        results['elapsed_seconds'] = monotonic()-started
        args.output.write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
        print(f'{name}: margin={results["cases"][name]["good_message_margin_bits"]:.6f}; '
            f'H1={results["cases"][name]["expected_H1"]:.3f}; '
            f'H2={results["cases"][name]["expected_H2"]:.3f}; '
            f'counts={results["cases"][name]["state_and_transition_counts"]}; '
            f'elapsed={monotonic()-case_start:.2f}s', flush=True)
    results['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', nargs='+', choices=['baseline16', 'actual24'], default=['baseline16', 'actual24'])
    parser.add_argument('--epsilon', type=float, default=1e-4)
    parser.add_argument('--receipt', type=Path, default=HERE.parent/'larger_state/trajectory_v1.json')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
