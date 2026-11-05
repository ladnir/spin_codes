"""Bounded FLOAT mark search for explicitly named, frozen research providers.

This only proposes a selected-weight witness; it does not certify a distance.
"""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib import import_module
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

from scipy.optimize import minimize


PROVIDERS = {
    'four-packet': 'hd_field128_lift_four_packet_screen',
    'full-rs': 'hd_field128_full_rs_screen',
}


def proposal(provider, r, z, alpha, distance, evaluations, seconds):
    if not (type(r) is int and 1 <= r <= 64 and 0 < z < 1 and 0 < alpha <= 1
            and 1 <= evaluations <= 1000 and isfinite(seconds) and seconds > 0):
        raise ValueError('Invalid selected-weight proposal parameters')
    started, history = monotonic(), []
    best, exhausted = None, False

    def objective(point):
        nonlocal best
        remaining = seconds - (monotonic() - started)
        if remaining <= 0:
            raise TimeoutError('Proposal wall-time budget exhausted')
        ell, exponent = point
        if not 1e-8 <= ell <= 64 or not 1e-6 <= exponent <= 1:
            return 1e100
        row = provider.weight_screen(r, 2.**-ell, float(exponent),
                                     distance=F(distance), max_seconds=remaining)
        history.append({key: row[key] for key in
                        ('z', 'alpha', 'log2_upper_screen', 'worst_profile')})
        if best is None or row['log2_upper_screen'] < best['log2_upper_screen']:
            best = row
        return row['log2_upper_screen']

    try:
        minimize(objective, [-log2(z), alpha], method='Nelder-Mead',
                 options=dict(maxfev=evaluations, xatol=1e-7, fatol=1e-5))
    except TimeoutError:
        exhausted = True
    if best is None:
        raise TimeoutError('No complete selected-weight result')
    best.update(proposal_evaluations=len(history), proposal_history=history,
                resource_limit_reached=exhausted, max_seconds=seconds,
                max_evaluations=evaluations, elapsed_seconds=monotonic()-started)
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--provider', choices=PROVIDERS, required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--z', type=float, required=True)
    parser.add_argument('--alpha', type=float, required=True)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--max-evaluations', type=int, default=48)
    parser.add_argument('--max-seconds', type=float, default=240.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new file')
    provider = import_module(PROVIDERS[args.provider])
    result = proposal(provider, args.r, args.z, args.alpha, F(args.distance),
                      args.max_evaluations, args.max_seconds)
    result['source_sha256'] = provider.source_pins()
    path = Path(__file__)
    result['source_sha256'][path.name] = sha256(path.read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('source_sha256', 'proposal_history')}, allow_nan=False))


if __name__ == '__main__':
    main()
