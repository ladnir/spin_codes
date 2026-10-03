"""Tilted comparison-path occupancy diagnostics from an authenticated receipt."""
import argparse
import hashlib
import json
from math import log
from pathlib import Path
from time import monotonic

import screen24
import numpy as np
from scipy.special import gammaln, logsumexp


def moment(local, q):
    regional, _ = screen24.prior.placement(local, q, epochs=64, windows=8)
    mixture = screen24.prior.uniform_mixture(regional, q)
    return screen24.prior.logarithmic.log_power_matrix(mixture, 32), regional, mixture


def marked_matrix_expectation(matrix, marked, *, regions=32):
    """Derivative of e0 Q^regions 1 for a nonnegative marked matrix."""
    size = len(matrix)
    forward = np.full((regions+1, size), -np.inf)
    forward[0, 0] = 0
    backward = np.zeros((regions+1, size))
    for i in range(regions):
        forward[i+1] = logsumexp(forward[i, :, None]+matrix, axis=0)
        backward[i+1] = logsumexp(matrix+backward[i, None, :], axis=1)
    total = float(logsumexp(forward[-1]))
    result = 0.
    for i in range(regions):
        result += np.exp(logsumexp(forward[i, :, None]+marked+backward[regions-1-i, None, :])-total)
    return float(result)


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    started = monotonic()
    receipt = json.loads(args.input.read_text(encoding='utf-8'))
    if not receipt['source_pins_verified_at_finish']:
        raise ValueError('finished source-authenticated input required')
    for name, expected in receipt['source_sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'stale source: {name}')
    local, q = np.asarray(receipt['local_operator_matrices']), receipt['q']
    baseline, regional, matrix = moment(local, q)
    if abs(baseline-receipt['statistics']['log_moment']) > 1e-8:
        raise ArithmeticError('input local moments do not reproduce the receipt')
    expected = []
    delta = 1e-4
    for j in range(9):
        values = []
        for sign in (-1, 1):
            changed = local.copy()
            changed[j] *= np.exp(sign*delta)
            values.append(moment(changed, q)[0])
        expected.append(float((values[1]-values[0])/(2*delta)))
    occupancies = np.arange(q+1)
    weights = (gammaln(q+1)-gammaln(occupancies+1)-gammaln(q-occupancies+1)
               +occupancies*log(255)-q*log(256))
    marked = logsumexp(regional[1:q+1]+weights[1:, None, None]
                      +np.log(occupancies[1:])[:, None, None], axis=0)
    packets = marked_matrix_expectation(matrix, marked)
    if abs(sum(expected)-2048) > 2e-4 or abs(sum(j*x for j, x in enumerate(expected))-packets) > 3e-4:
        raise ArithmeticError('occupancy derivatives failed total/packet checks')
    result = dict(schema='actual24-bit-comparison-occupancy-1', q=q, tilt=receipt['tilt'],
        proposal_only=True, whole_code_certificate=False, actual_failure_distribution=False,
        input_receipt=str(args.input.resolve()), input_receipt_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        input_source_sha256=receipt['source_sha256'],
        source_sha256={str(Path(__file__).resolve()): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        finite_difference_delta=delta, expected_physical_steps_by_occupancy=expected,
        total_expected_steps=sum(expected), expected_occupied_packets=packets,
        expected_occupied_packets_from_step_derivatives=sum(j*x for j, x in enumerate(expected)),
        unconditional_uniform_expected_packets=32*q*255/256,
        mean_packets_per_occupied_step=packets/(2048-expected[0]),
        elapsed_seconds=monotonic()-started)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
