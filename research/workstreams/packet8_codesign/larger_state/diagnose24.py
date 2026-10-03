"""One-point path diagnostics for the24-bit comparison operator, not the code.

Every restricted path below is a contribution to the positive comparison
expression. None is a lower bound on the actual code's failure probability.
"""
import argparse
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import screen24
import numpy as np
from scipy.special import logsumexp


def path_statistics(matrix, regions):
    size = len(matrix)
    forward = np.full((regions+1, size), -np.inf)
    forward[0, 0] = 0.
    backward = np.zeros((regions+1, size))
    for k in range(regions):
        forward[k+1] = logsumexp(forward[k, :, None]+matrix, axis=0)
        backward[k+1] = logsumexp(matrix+backward[k, None, :], axis=1)
    total = float(logsumexp(forward[-1]))
    transitions = np.zeros_like(matrix)
    for k in range(regions):
        transitions += np.exp(forward[k, :, None]+matrix+backward[regions-1-k, None, :]-total)
    value = forward[0].copy()
    parents = []
    for k in range(regions):
        choices = value[:, None]+matrix
        parent = np.argmax(choices, axis=0)
        value = choices[parent, np.arange(size)]
        parents.append(parent)
    final = int(np.argmax(value))
    maximum = float(value[final])
    path = [final]
    for parent in reversed(parents):
        path.append(int(parent[path[-1]]))
    path.reverse()
    return dict(log_moment=total, terminal_probabilities=np.exp(forward[-1]-total).tolist(),
        expected_boundary_transitions=transitions.tolist(), maximum_path=path,
        maximum_path_log_moment=maximum, full_over_maximum_path_log2=(total-maximum)/log(2))


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    started = monotonic()
    pins = screen24.sources()
    pins[str(Path(__file__).resolve())] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    data = screen24.maps24.make_maps()
    counts = screen24.maps24.census(data)
    local, diagnostics = screen24.local_operators(data, counts, np.exp(-args.tilt), progress=True)
    regional, _ = screen24.prior.placement(local, args.q, epochs=64, windows=8)
    matrix = screen24.prior.uniform_mixture(regional, args.q)
    statistics = path_statistics(matrix, 32)
    beta = screen24.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
    fixed = log(comb(512, args.q))+args.q*(log(beta.numerator)-log(beta.denominator))+args.tilt*13107
    margin = lambda value: -(fixed+float(value))/log(2)
    statistics['full_margin_bits'] = margin(statistics['log_moment'])
    statistics['maximum_path_margin_bits'] = margin(statistics['maximum_path_log_moment'])
    statistics['always_zero_regional_boundary_margin_bits'] = margin(32*matrix[0, 0])
    statistics['uniform_after_first_region_margin_bits'] = margin(matrix[0, 1]+31*matrix[1, 1])
    statistics['zero_or_uniform_boundaries_margin_bits'] = margin(
        screen24.prior.logarithmic.log_power_matrix(matrix[:2, :2], 32))
    perturbed = []
    delta = 1e-4
    for direction in (-1, 1):
        changed = local.copy()
        changed[:, 1:, 0] *= np.exp(direction*delta)
        transfer, _ = screen24.prior.placement(changed, args.q, epochs=64, windows=8)
        moment = screen24.prior.logarithmic.log_power_matrix(screen24.prior.uniform_mixture(transfer, args.q), 32)
        perturbed.append(float(moment))
    statistics['expected_physical_return_count_finite_difference'] = (perturbed[1]-perturbed[0])/(2*delta)
    labels = ['zero', 'uniform_nonzero']+[f'birth_{j}' for j in range(1, 9)]
    statistics['maximum_path_labels'] = [labels[i] for i in statistics['maximum_path']]
    result = dict(schema='actual24-bit-comparison-path-diagnostic-1',
        proposal_only=True, whole_code_certificate=False, actual_failure_lower_bound=False,
        q=args.q, tilt=args.tilt, map_record=data['record'], coordinate_labels=labels,
        source_sha256=pins, source_pins_verified_at_finish=False,
        local_operator_matrices=local.tolist(), regional_mixture_log_matrix=matrix.tolist(),
        local_diagnostics=diagnostics, statistics=statistics, elapsed_seconds=monotonic()-started)
    current = screen24.sources()
    current[str(Path(__file__).resolve())] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if current != pins:
        raise ArithmeticError('source changed during diagnostic')
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(statistics, indent=2), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q', type=int, default=119)
    parser.add_argument('--tilt', type=float, default=.5)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
