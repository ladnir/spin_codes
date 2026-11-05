"""Bounded FLOAT mark search at original upper-half weights; counts alone complement.

This is a separate adapter. The frozen lower-half search and all numerical
providers remain unchanged. Every trial rebuilds gamma coefficients and sums
the complete feasible midpoint-profile domain at the original weight.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_gamma_midpoint_mark_search as frozen


mode, provider, HERE = frozen.mode, frozen.provider, frozen.HERE


def source_pins():
    pins = frozen.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def checked_domain(weight):
    if type(weight) is not int or not 4096 <= weight < 8192:
        raise ValueError('Expected original weight4096..8191')
    full = provider.alias.full_counts(4096, weight)
    if list(map(int, full)) != list(range(weight-4096, weight//2+1)):
        raise ArithmeticError('Unexpected original-weight feasible domain')
    return full


def run_search(source_path, weight, *, max_evaluations=40, max_seconds=120., checkpoint=None):
    full = checked_domain(weight)
    if (type(max_evaluations) is not int or not 1 <= max_evaluations <= 40
            or not isfinite(max_seconds) or not 0 < max_seconds <= 120):
        raise ValueError('Expected at most40 evaluations and120 seconds')
    path = Path(source_path).resolve()
    raw, pins = path.read_bytes(), source_pins()
    source = json.loads(raw)
    if weight not in source.get('requested_weights', [source.get('active_groups')]):
        raise ValueError('Original weight must be explicitly calibrated in source')
    state, distance = source['state_dimension'], F(source['target_distance'])
    if source['construction'] != provider.gamma.construction(state) or not 0 < distance < F(1, 2):
        raise ValueError('Wrong source construction or distance')
    original = source['witness']
    markers, marker_key = frozen.enrich.positive_markers(original)
    literals = [list(map(str, row)) for row in markers]
    odds = [tuple(log2(float(v)) for v in row) for row in markers]
    z_exact, alpha_exact = F(original['z']), F(original['alpha'])
    z0, alpha0 = float(z_exact), float(alpha_exact)
    if F(z0) != z_exact or F(alpha0) != alpha_exact:
        raise ValueError('Initial marks must be exactly binary64-representable')
    if not 2.**-1022 <= z0 <= 2.**-1e-8 or not 2.**-32 <= alpha0 <= 1:
        raise ValueError('Initial marks outside logarithmic chart')
    grid = original.get('grid', 128)
    if grid not in (64, 128, 256):
        raise ValueError('Unsupported alias grid')
    for name, digest in source.get('source_sha256', {}).items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting source pin: '+name)
        pins[name] = digest
    provenance = dict(path=str(path), sha256=sha256(raw).hexdigest())

    def check():
        if path.read_bytes() != raw:
            raise RuntimeError('Source changed')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Helper source changed: '+name)

    check()
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
            raise TimeoutError('Attempted bank limit reached')
        number = attempted
        attempted += 1
        try:
            coeff = provider.AliasCoefficients(z, alpha, state=state, _budget=budget)
            bank = provider.AliasBankScreen(z, alpha, odds, 0., state=state, distance=distance,
                        grid=grid, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
            regions = bank.regional(weight)
            if regions.shape != (len(full), 2, 2):
                raise ArithmeticError('Original-weight regional domain mismatch')
            row = mode.support_sum(4096, weight, regions, counts=counts, budget=budget)
            if row['profile_pairs_evaluated'] != len(full)**2:
                raise ArithmeticError('Incomplete original-weight profile sum')
            prefix = bank.constant+bank.alpha*bank.density*weight
            score = prefix+row['log2_support_sum_screen']
            if not isfinite(score):
                raise ArithmeticError('Nonfinite complete objective')
            row.update(schema='gamma-dense-midpoint-mark-trial-FLOAT-v1', evaluation_index=number,
                active_groups=weight, state_dimension=state, construction=bank.metadata()['construction'],
                target_distance=str(distance), source=provenance, source_sha256=dict(sorted(pins.items())),
                feasible_full_range=[int(full[0]), int(full[-1])], canonical_count_weight=8192-weight,
                count_full_shift=weight-4096, count_complement_only=True,
                original_weight_used_for_histograms_markers_matrices_and_prefix=True,
                log2_outer_prefix=prefix, log2_upper_screen=score,
                witness=dict(z=str(F(z)), alpha=str(F(alpha)), markers=literals, grid=grid,
                    state_dimension=state, marker_log_odds_exact=[[str(F(v)) for v in pair] for pair in odds]),
                regional_alias_diagnostics=dict(bank.regional_diagnostics),
                regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                coefficient_table_sha256=sha256(coeff.composed.tobytes()).hexdigest(),
                cached_count_provenance=counts.diagnostics(), complete_selected_weight_screen=True,
                common_potential_used=False, old_numerical_bounds_reused=False,
                elapsed_seconds=monotonic()-started, directed=False,
                selected_weight_certificate=False, whole_code_certificate=False)
            check()
            link = checkpoint(row) if checkpoint else None
            history.append(dict(evaluation_index=number, z=str(F(z)), alpha=str(F(alpha)),
                log2_upper_screen=score, largest_profile=row['largest_profile'], checkpoint=link,
                regional_alias_diagnostics=row['regional_alias_diagnostics']))
            if best is None or score < best['log2_upper_screen']:
                best = row
            cache[key] = score
            return score
        except (TimeoutError, ArithmeticError, MemoryError) as error:
            failures.append(dict(evaluation_index=number, z=str(F(z)), alpha=str(F(alpha)),
                error_type=type(error).__name__, message=str(error), complete_selected_weight_screen=False))
            if isinstance(error, (TimeoutError, MemoryError)):
                raise
            cache[key] = 1e300
            return cache[key]

    try:
        counts = mode.MidpointCounts(4096, weight, max_singletons=64, budget=budget)
        fitted = minimize(objective, initial, method='Nelder-Mead',
            options=dict(maxfev=max_evaluations, initial_simplex=simplex, xatol=1e-7, fatol=1e-5))
        outcome = dict(success=bool(fitted.success), message=str(fitted.message), function_calls=int(fitted.nfev))
    except (TimeoutError, MemoryError) as error:
        outcome = dict(success=False, message=str(error), resource_stop=type(error).__name__)
    check()
    return dict(schema='gamma-dense-midpoint-fixed-marker-search-FLOAT-v1', active_groups=weight,
        state_dimension=state, construction=provider.gamma.construction(state), target_distance=str(distance),
        source=provenance, source_sha256=dict(sorted(pins.items())), source_marker_key=marker_key,
        initial_witness=original, fixed_positive_marker_literals=literals,
        initial_z_preserved_exactly=str(z_exact), initial_alpha_preserved_exactly=str(alpha_exact),
        initial_simplex=[[str(F(v)) for v in row] for row in simplex],
        best=best, history=history, failures=failures, attempted_banks=attempted, completed_trials=len(history),
        counts_built=int(counts is not None), optimizer=outcome, max_evaluations=max_evaluations,
        max_seconds=max_seconds, elapsed_seconds=monotonic()-started, common_potential_used=False,
        count_complement_only=True, original_weight_used_for_histograms_markers_matrices_and_prefix=True,
        old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, scope='One complete-profile FLOAT objective per completed trial; no contiguous certificate')
