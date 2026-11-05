"""Bounded FLOAT dense calibration from fresh all-pair gamma coefficients."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_gamma_midpoint_enrich as literals


provider, mode, HERE = literals.provider, literals.mode, literals.HERE
PAIRS = 4096


def source_pins():
    pins = literals.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def histogram(pairs, weight, full):
    if any(type(v) is not int for v in (pairs, weight, full)):
        raise ValueError('Integer histogram parameters required')
    result = pairs-weight+full, weight-2*full, full
    if pairs < 1 or not 0 <= weight <= 2*pairs or min(result) < 0:
        raise ValueError('Infeasible original-weight histogram')
    return result


def target_full_counts(pairs, weight):
    lo, hi = max(0, weight-pairs), weight//2
    histogram(pairs, weight, lo)
    span = hi-lo
    return tuple(sorted({lo+min(span, offset) for offset in
                         (0, 1, 4, 16, span//4, span//2, 3*span//4, span)}))


def fit_entry(numerators, pairs, weight, full, entry, old_odds, budget, *, max_calls=128):
    """Finite Powell proposal; every objective call, including seeds, is bounded."""
    if type(entry) is not int or not 0 <= entry < 4 or not 8 <= max_calls <= 128:
        raise ValueError('Expected one of four entries and8..128 objective calls')
    empty, single, full = histogram(pairs, weight, full)
    if not full:
        raise ValueError('Dense targets require positive full count')
    if empty == 0:
        # The full regional degree is at most pairs. A fixed finite radial
        # coordinate avoids division by empty and makes no limiting claim.
        chart = 'empty-zero-finite-radial-t64'
        center = (log2(single/full),)
        encode = lambda point: (float(point[0])+64., 64.)
        projected = [(u-v,) for u, v in old_odds]
        bounds = [(-64., 64.)]
    elif single == 0:
        chart = 'single-zero-positive-x2^-64'
        center = (log2(full/empty),)
        encode = lambda point: (-64., float(point[0]))
        projected = [(v,) for u, v in old_odds]
        bounds = [(-64., 128.)]
    else:
        chart = 'positive-interior-log-markers'
        center = (log2(single/empty), log2(full/empty))
        encode = lambda point: tuple(map(float, point))
        projected = list(old_odds)
        bounds = [(-64., 128.)]*2
    def clip(point):
        return tuple(max(lo, min(hi, value)) for value, (lo, hi) in zip(point, bounds))
    center = clip(center)
    nearest = sorted(map(clip, projected), key=lambda p: sum((a-b)**2 for a, b in zip(p, center)))[:4]
    seeds = tuple(dict.fromkeys((center, *nearest)))
    best, best_point, best_value, calls = None, None, float('inf'), 0
    def objective(point):
        nonlocal best, best_point, best_value, calls
        budget.check()
        if calls >= max_calls:
            raise ArithmeticError('Optimizer exceeded its per-fit call bound')
        calls += 1
        point = tuple(map(float, point))
        odds = encode(point)
        value = float(numerators(odds)[entry]-single*odds[0]-full*odds[1])
        if not isfinite(value):
            raise ArithmeticError('Nonfinite positive fitting objective')
        if value < best_value:
            best, best_point, best_value = odds, point, value
        return value
    for point in seeds:
        objective(point)
    fitted = minimize(objective, best_point, method='Powell', bounds=bounds,
                      options=dict(maxfev=max_calls-calls, xtol=1e-8, ftol=1e-10))
    positive = tuple(F(2.**value) for value in best)
    if not all(v > 0 and isfinite(float(v)) for v in positive):
        raise ArithmeticError('Nonfinite or zero exact positive proposal')
    actual_odds = tuple(log2(float(v)) for v in positive)
    return positive, dict(active_groups=weight, histogram=dict(empty=empty, single=single, full=full),
                          entry=entry, chart=chart, chart_bounds=bounds, objective_calls=calls,
                          optimizer_success=bool(fitted.success), optimizer_message=str(fitted.message),
                          positive_objective=best_value, fitted_log_literals=[str(F(v)) for v in best],
                          positive_marker_literals=list(map(str, positive)),
                          evaluated_marker_log_literals=[str(F(v)) for v in actual_odds])


def regional_with_selection(bank, weight, budget, *, pairs=PAIRS):
    """Preserve the frozen normalization; index every entry by actual full count."""
    full = provider.alias.full_counts(pairs, weight)
    single = weight-2*full
    regions = bank.regional(weight)
    if regions.shape != (len(full), 2, 2):
        raise ArithmeticError('Original-weight feasible regional shape mismatch')
    best = np.full((4, len(full)), np.inf)
    reliable_best = np.zeros_like(best, dtype=bool)
    for odds, table in zip(bank.odds, bank.tables):
        budget.check()
        reliable = table.reliable[:, single % bank.grid, full % bank.grid]
        values = np.where(reliable, table.logs[:, single % bank.grid, full % bank.grid],
                          table.positive_numerator[:, None])-single*odds[0]-full*odds[1]
        improved = values < best
        reliable_best[improved] = reliable[improved]
        np.minimum(best, values, out=best)
    return regions, reliable_best.T.reshape(regions.shape)


def run(source_path, weight, *, max_seconds=90., max_calls=128):
    if (type(weight) is not int or not PAIRS < weight < 2*PAIRS
            or not isfinite(max_seconds) or not 0 < max_seconds <= 90
            or type(max_calls) is not int or not 8 <= max_calls <= 128):
        raise ValueError('Expected original dense weight4097..8191, at most90 seconds/128 calls per fit')
    path = Path(source_path).resolve()
    raw, pins = path.read_bytes(), source_pins()
    saved = json.loads(raw)
    construction = provider.gamma.construction(30)
    if (saved['state_dimension'] != 30 or F(saved['target_distance']) != F(2, 5)
            or saved['construction'] != construction):
        raise ValueError('Expected the gamma30/d40 construction')
    witness = saved['witness']
    if (witness.get('state_dimension', 30) != 30 or witness.get('inner_width', 64) != 64
            or int(witness.get('field_order', saved.get('field_order', 1 << 128))) != 1 << 128):
        raise ValueError('Source field or dimensions disagree')
    positives, marker_key = literals.positive_markers(witness)
    z, alpha = F(witness['z']), F(witness['alpha'])
    if not 0 < float(z) < 1 or not 0 < float(alpha) <= 1:
        raise ValueError('Invalid fixed z/alpha')
    grid = witness.get('grid', 128)
    if grid not in (64, 128, 256):
        raise ValueError('Unsupported grid')
    calibration, parent_pin = saved.get('active_groups'), None
    if calibration is None:
        proposal = saved['witness_proposal']
        parent_positive, unused = literals.positive_markers(proposal['witness'])
        if (parent_positive != positives or any(F(proposal['witness'][k]) != F(witness[k]) for k in ('z', 'alpha'))
                or F(proposal['new_target_distance']) != F(2, 5)):
            raise ValueError('Embedded proposal disagrees with exact bank literals')
        calibration, parent_pin = proposal['source_calibration_weight'], proposal['source']
    if type(calibration) is not int or not 1 <= calibration <= 8192:
        raise ValueError('Explicit source calibration required')
    targets = target_full_counts(PAIRS, weight)
    initial_count = len(positives)
    if initial_count+4*len(targets) > 128:
        raise ValueError('Retained source plus all requested entry fits exceeds128 markers')
    for name, digest in saved.get('source_sha256', {}).items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting source pin: '+name)
        pins[name] = digest
    def check():
        if path.read_bytes() != raw:
            raise RuntimeError('Exact source artifact changed')
        if parent_pin and sha256(Path(parent_pin['path']).read_bytes()).hexdigest() != parent_pin['sha256']:
            raise RuntimeError('Embedded literal-proposal source changed')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source closure changed: '+name)
    check()
    odds = [tuple(log2(float(v)) for v in pair) for pair in positives]
    started, budget = monotonic(), provider.old.Budget(max_seconds)
    cache, rounds, fits = {}, [], []
    pending, status, coefficient_hash, initial_metadata = dict(stage='initial-bank'), 'complete', None, None
    coeff, counts = None, None
    full = provider.alias.full_counts(PAIRS, weight)
    def numerators(point):
        budget.check()
        point = tuple(map(float, point))
        if point not in cache:
            cache[point] = np.asarray(provider.alias.power(coeff.pair_logs(point, _budget=budget), 64)).ravel()
        return cache[point]
    def complete_round(number, target):
        nonlocal pending
        pending = dict(stage='all-profile-sum', round=number, target_full=target, marker_count=len(positives))
        row = mode.support_sum(PAIRS, weight, regions, counts=counts, budget=budget)
        if row['profile_pairs_evaluated'] != len(full)**2:
            raise ArithmeticError('Incomplete original-weight profile sum')
        row.update(round=number, target_full=target, active_groups=weight, marker_count=len(positives),
                   log2_outer_prefix=prefix, log2_upper_screen=prefix+row['log2_support_sum_screen'],
                   regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                   selected_reliable_entries=int(reliability.sum()),
                   selected_positive_fallback_entries=int((~reliability).sum()),
                   complete_selected_weight_screen=True, elapsed_seconds=monotonic()-started)
        rounds.append(row)
        pending = None
        print(json.dumps(dict(round=number, target_full=target, marker_count=len(positives),
                              score=row['log2_upper_screen'], largest=row['largest_profile'],
                              elapsed_seconds=row['elapsed_seconds'])), flush=True)
    try:
        coeff = provider.AliasCoefficients(float(z), float(alpha), state=30, _budget=budget)
        coefficient_hash = sha256(coeff.composed.tobytes()).hexdigest()
        bank = provider.AliasBankScreen(float(z), float(alpha), odds, 0., state=30, distance=F(2, 5),
                    grid=grid, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
        initial_metadata = bank.metadata()
        regions, reliability = regional_with_selection(bank, weight, budget)
        prefix = bank.constant+bank.alpha*bank.density*weight
        del bank
        counts = mode.MidpointCounts(PAIRS, weight, max_singletons=64, budget=budget)
        complete_round(0, None)
        for number, target in enumerate(targets, 1):
            pending = dict(stage='fit', round=number, target_full=target, marker_count=len(positives))
            proposed = []
            for entry in range(4):
                positive, record = fit_entry(numerators, PAIRS, weight, target, entry, odds, budget, max_calls=max_calls)
                proposed.append(positive)
                fits.append(record)
            fresh = list(dict.fromkeys(pair for pair in proposed if pair not in positives))
            if fresh:
                pending['stage'] = 'fresh-atlases'
                fresh_odds = [tuple(log2(float(v)) for v in pair) for pair in fresh]
                bank = provider.AliasBankScreen(float(z), float(alpha), fresh_odds, 0., state=30, distance=F(2, 5),
                            grid=grid, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
                added, reliable = regional_with_selection(bank, weight, budget)
                del bank
                improved = added < regions
                reliability[improved] = reliable[improved]
                np.minimum(regions, added, out=regions)
                positives.extend(fresh)
                odds.extend(fresh_odds)
            complete_round(number, target)
    except (TimeoutError, MemoryError) as error:
        status = 'timeout' if isinstance(error, TimeoutError) else 'memory-error'
    if coeff is not None and sha256(coeff.composed.tobytes()).hexdigest() != coefficient_hash:
        raise RuntimeError('Shared fresh coefficient table changed')
    check()
    return dict(schema='gamma-original-weight-dense-calibration-FLOAT-v1', status=status,
                construction=construction, state_dimension=30, target_distance='2/5', active_groups=weight,
                feasible_full_range=[int(full[0]), int(full[-1])], requested_full_targets=list(targets),
                source=dict(path=str(path), sha256=sha256(raw).hexdigest()), embedded_literal_source=parent_pin,
                source_calibration_active_groups=calibration, source_marker_key=marker_key,
                source_sha256=dict(sorted(pins.items())), source_role='exact literals only; no old numerical bounds',
                initial_marker_count=initial_count, marker_count=len(positives),
                initial_bank_metadata=initial_metadata, coefficient_table_sha256=coefficient_hash,
                fresh_coefficient_table_built=coeff is not None, shared_fresh_coefficients=True,
                witness=dict(z=str(z), alpha=str(alpha), markers=[list(map(str, p)) for p in positives], grid=grid,
                             marker_log_odds_exact=[[str(F(v)) for v in p] for p in odds], state_dimension=30),
                fits=fits, rounds=rounds, last_complete_round=rounds[-1] if rounds else None,
                pending_round=pending, witness_may_include_incomplete_round=pending is not None,
                completed_target_count=max(0, len(rounds)-1), counts_built=int(counts is not None),
                max_seconds=max_seconds, max_calls_per_fit=max_calls, max_fits=4*len(targets),
                distinct_positive_points=len(cache), elapsed_seconds=monotonic()-started,
                conversion='Original positive rationals retained; fitted exp2 binary64 values frozen as exact rational proposals',
                count_complement_only=True, original_weight_used_for_histograms_markers_matrices_and_prefix=True,
                common_potential_used=False, old_numerical_bounds_reused=False, no_weight_interpolation=True,
                directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                scope='Completed original-weight FLOAT sums only; no finite-boundary-limit or certificate claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--weight', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=90.)
    parser.add_argument('--max-calls', type=int, default=128)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = run(args.source, args.weight, max_seconds=args.max_seconds, max_calls=args.max_calls)
    print(json.dumps(dict(status=result['status'], artifact=provider.search.write_new(args.output, result))), flush=True)


if __name__ == '__main__':
    main()
