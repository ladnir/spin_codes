"""Retune one common state potential without refitting a frozen alias bank.

All calculations here propose witnesses in FLOAT. The saved exact z,
alpha, and coefficient markers remain unchanged. Directed replay must
reconstruct the complete bound before any new interval is certified.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import isfinite, log2

from scipy.optimize import minimize_scalar

from hd_shifted_potential_band_neighbor import fixed_screen_from_artifact, checked_window
from hd_shifted_run_count_screen import run_log_screen_float


def retune_potential(artifact, weight, *, bounds=(-32., 0.), outer_measure=None):
    """Return a fixed bank with one exactified common potential, plus diagnostics."""
    if (type(weight) is not int or not 1 <= weight <= 8192
            or len(bounds) != 2 or not all(map(isfinite, bounds))
            or not -512 <= bounds[0] < bounds[1] <= 512):
        raise ValueError('An admissible central weight and bounded log-potential interval are required')
    screen = fixed_screen_from_artifact(artifact, outer_measure=outer_measure)
    screen.source_witness = deepcopy(screen.source_witness)
    regions = screen.regional(weight)
    objective = lambda ratio: run_log_screen_float(4096, weight, regions, ratio)
    initial = screen.log_ratio
    optimized = minimize_scalar(objective, method='bounded', bounds=bounds,
                                options={'xatol': 1.e-9, 'maxiter': 80})
    candidates = (initial, bounds[0], bounds[1], float(optimized.x))
    scored = [(objective(ratio), ratio) for ratio in candidates]
    if not all(isfinite(value) and isfinite(ratio) for value, ratio in scored):
        raise ArithmeticError('Nonfinite common-potential proposal')
    _, ratio = min(scored)
    literal = str(F(2.**ratio))
    # Use the ratio of the exactified binary64 potential in subsequent screens.
    screen.log_ratio = log2(F(literal))
    screen.source_witness['potential'] = literal
    screen.source_witness.pop('alpha_and_potential', None)
    screen.source_witness['alpha_preserved_from_source'] = True
    diagnostics = dict(central_weight=weight, initial_log2_ratio=initial,
        selected_log2_ratio=screen.log_ratio, evaluations=int(optimized.nfev),
        initial_support_sum=objective(initial), selected_support_sum=objective(screen.log_ratio),
        log2_ratio_bounds=list(bounds), new_marker_fits=0, directed=False,
        requires_fresh_directed_replay=True)
    screen.source_witness['potential_retuning'] = diagnostics
    return screen, diagnostics


if __name__ == '__main__':
    import argparse
    import json
    from hashlib import sha256
    from pathlib import Path
    from unittest.mock import patch
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--central-weight', type=int, required=True)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--max-weights', type=int, default=33)
    parser.add_argument('--threshold', type=float, default=-60.)
    parser.add_argument('--outer-measure')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checked_window(args.first, args.last, args.max_weights)
    if not isfinite(args.threshold):
        parser.error('A finite threshold is required')
    if args.output.exists():
        parser.error('Output exists; frozen artifacts are never overwritten')
    raw = args.source.read_bytes()
    screen, diagnostics = retune_potential(json.loads(raw), args.central_weight,
                                           outer_measure=args.outer_measure)
    with patch('hd_shifted_potential_band_screen.run_log_screen', run_log_screen_float):
        result = screen.band(args.first, args.last, max_weights=args.max_weights,
                             threshold=args.threshold)
    result.update(source_artifact_sha256=sha256(raw).hexdigest(), potential_retuning=diagnostics,
                  count_method='FLOAT-log-factorial-cyclic-run-counts',
                  requires_fresh_directed_replay=True)
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, allow_nan=False, separators=(',', ':'))
        output.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()), potential_retuning=diagnostics,
                         failed_weights=result['failed_weights'],
                         union_log2_upper_screen=result['union_log2_upper_screen'], directed=False)), flush=True)
