"""Bounded fixed-marker FLOAT z/alpha search using cached midpoint counts."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_gamma_midpoint_enrich as enrich


mode, provider, HERE = enrich.mode, enrich.provider, enrich.HERE


def source_pins():
    pins = enrich.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run_search(source_path, weight, *, max_evaluations=40, max_seconds=120., checkpoint=None):
    if (type(weight) is not int or not 2 <= weight < 4096 or type(max_evaluations) is not int
            or not 1 <= max_evaluations <= 40 or not isfinite(max_seconds) or not 0 < max_seconds <= 120):
        raise ValueError('Expected weight2..4095, at most40 evaluations and120 seconds')
    path = Path(source_path).resolve()
    raw, pins = path.read_bytes(), source_pins()
    artifact = json.loads(raw)
    if weight not in artifact.get('requested_weights', [artifact.get('active_groups')]):
        raise ValueError('Weight must be explicitly calibrated in the source')
    state, distance = artifact['state_dimension'], F(artifact['target_distance'])
    if artifact['construction'] != provider.gamma.construction(state) or not 0 < distance < F(1, 2):
        raise ValueError('Wrong source construction or distance')
    original = artifact['witness']
    markers, marker_key = enrich.positive_markers(original)
    odds = [tuple(log2(float(value)) for value in pair) for pair in markers]
    literals = [list(map(str, pair)) for pair in markers]
    z_exact, alpha_exact = F(original['z']), F(original['alpha'])
    z0, alpha0 = float(z_exact), float(alpha_exact)
    if F(z0) != z_exact or F(alpha0) != alpha_exact:
        raise ValueError('This FLOAT search requires exactly representable initial binary64 marks')
    if not 2.**-1022 <= z0 <= 2.**-1e-8 or not 2.**-32 <= alpha0 <= 1:
        raise ValueError('Initial marks lie outside the bounded logarithmic chart')
    grid = original.get('grid', 128)
    if grid not in (64, 128, 256):
        raise ValueError('Unsupported grid')
    for name, digest in artifact.get('source_sha256', {}).items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting source pin: '+name)
        pins[name] = digest
    provenance = dict(path=str(path), sha256=sha256(raw).hexdigest())
    def check_pins():
        if path.read_bytes() != raw:
            raise RuntimeError('Source proposal changed')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Search source changed: '+name)
    check_pins()
    started, budget = monotonic(), provider.old.Budget(max_seconds)
    history, failures, cache, best = [], [], {}, None
    attempted, counts = 0, None
    initial = np.array((-log2(z0), log2(alpha0)))
    simplex = np.array((initial, (min(1022., initial[0]+.5), initial[1]),
                        (initial[0], initial[1]+.25 if initial[1] <= -.25 else initial[1]-.25)))
    if simplex[1, 0] == initial[0]:
        simplex[1, 0] = initial[0]-.5

    def objective(point):
        nonlocal attempted, best
        budget.check()
        ell, exponent = map(float, point)
        if not 1e-8 <= ell <= 1022 or not -32 <= exponent <= 0:
            return 1e300
        z = z0 if ell == initial[0] else 2.**-ell
        alpha = alpha0 if exponent == initial[1] else 2.**exponent
        key = z, alpha
        if key in cache:
            return cache[key]
        if attempted >= max_evaluations:
            raise TimeoutError('Completed/attempted bank limit reached')
        number = attempted
        attempted += 1
        try:
            coeff = provider.AliasCoefficients(z, alpha, state=state, _budget=budget)
            bank = provider.AliasBankScreen(z, alpha, odds, 0., state=state, distance=distance, grid=grid,
                        floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
            regions = bank.regional(weight)
            row = mode.support_sum(4096, weight, regions, counts=counts, budget=budget)
            prefix = bank.constant+bank.alpha*bank.density*weight
            score = prefix+row['log2_support_sum_screen']
            if not isfinite(score):
                raise ArithmeticError('Nonfinite complete midpoint objective')
            row.update(schema='gamma-midpoint-mark-trial-FLOAT-v1', evaluation_index=number,
                       active_groups=weight, state_dimension=state, construction=bank.metadata()['construction'],
                       target_distance=str(distance), source=provenance, source_sha256=dict(sorted(pins.items())),
                       log2_outer_prefix=prefix, log2_upper_screen=score,
                       witness=dict(z=str(F(z)), alpha=str(F(alpha)), markers=literals, grid=grid,
                                    state_dimension=state, marker_log_odds_exact=[[str(F(v)) for v in pair] for pair in odds]),
                       regional_alias_diagnostics=dict(bank.regional_diagnostics),
                       regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                       coefficient_table_sha256=sha256(coeff.composed.tobytes()).hexdigest(),
                       cached_count_provenance=counts.diagnostics(), complete_selected_weight_screen=True,
                       common_potential_used=False, neutral_unused_log_potential=0,
                       potential_role='new neutral interface placeholder; not inherited and not used by objective',
                       old_numerical_bounds_reused=False, elapsed_seconds=monotonic()-started,
                       directed=False, selected_weight_certificate=False, whole_code_certificate=False)
            check_pins()
            link = checkpoint(row) if checkpoint else None
            history.append(dict(evaluation_index=number, z=str(F(z)), alpha=str(F(alpha)), log2_upper_screen=score,
                                largest_profile=row['largest_profile'], checkpoint=link,
                                regional_alias_diagnostics=row['regional_alias_diagnostics']))
            if best is None or score < best['log2_upper_screen']:
                best = row
            cache[key] = score
            return score
        except (TimeoutError, ArithmeticError) as error:
            failures.append(dict(evaluation_index=number, z=str(F(z)), alpha=str(F(alpha)),
                                 error_type=type(error).__name__, message=str(error), complete_selected_weight_screen=False))
            if isinstance(error, TimeoutError):
                raise
            cache[key] = 1e300
            return cache[key]
    try:
        counts = mode.MidpointCounts(4096, weight, max_singletons=64, budget=budget)
        fitted = minimize(objective, initial, method='Nelder-Mead',
                          options=dict(maxfev=max_evaluations, initial_simplex=simplex, xatol=1e-7, fatol=1e-5))
        outcome = dict(success=bool(fitted.success), message=str(fitted.message), function_calls=int(fitted.nfev))
    except TimeoutError as error:
        outcome = dict(success=False, message=str(error), timeout=True)
    check_pins()
    return dict(schema='gamma-midpoint-fixed-marker-mark-search-FLOAT-v1', active_groups=weight,
                state_dimension=state, construction=provider.gamma.construction(state), target_distance=str(distance),
                initial_witness=original, source=provenance, source_sha256=dict(sorted(pins.items())),
                source_marker_key=marker_key, fixed_positive_marker_literals=literals,
                initial_z_preserved_exactly=str(z_exact), initial_alpha_preserved_exactly=str(alpha_exact),
                initial_simplex=[[str(F(v)) for v in row] for row in simplex],
                marker_log_conversion='log2(float(saved exact positive rational)); unchanged across every trial',
                best=best, history=history, failures=failures, attempted_banks=attempted, completed_trials=len(history),
                counts_built=int(counts is not None), count_cutoff=64, optimizer=outcome, max_evaluations=max_evaluations,
                max_seconds=max_seconds, elapsed_seconds=monotonic()-started, common_potential_used=False,
                support_method='complete midpoint-profile sum with all-pair chronological transfer',
                old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
                whole_code_certificate=False, scope='One selected-weight FLOAT mark search; no contiguous-band certificate')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--max-evaluations', type=int, default=40)
    parser.add_argument('--max-seconds', type=float, default=120.)
    args = parser.parse_args()
    args.output_dir.mkdir(exist_ok=False)
    def checkpoint(row):
        pin = provider.search.write_new(args.output_dir/f"evaluation_{row['evaluation_index']:03d}.json", row)
        print(json.dumps(dict(evaluation=row['evaluation_index'], z=row['witness']['z'], alpha=row['witness']['alpha'],
                              log2_upper_screen=row['log2_upper_screen'], largest_profile=row['largest_profile'],
                              elapsed_seconds=row['elapsed_seconds'])), flush=True)
        return pin
    result = run_search(args.source, args.r, max_evaluations=args.max_evaluations, max_seconds=args.max_seconds, checkpoint=checkpoint)
    pin = provider.search.write_new(args.output_dir/'summary.json', result)
    print(json.dumps(dict(artifact=pin, completed_trials=result['completed_trials'], optimizer=result['optimizer'],
                          best=None if result['best'] is None else result['best']['log2_upper_screen'])), flush=True)


if __name__ == '__main__':
    main()
