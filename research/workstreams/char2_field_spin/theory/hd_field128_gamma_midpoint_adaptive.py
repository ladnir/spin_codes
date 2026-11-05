"""Bounded FLOAT adaptive markers at lower-half original weights.

Rebuild one initial bank, then only one atlas per new marker. Every completed
addition is scored with the full midpoint-profile sum; old bounds are unused.
"""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_gamma_midpoint_dense_adaptive as dense

frozen = dense.frozen


provider, mode, HERE, PAIRS = frozen.provider, frozen.mode, frozen.HERE, frozen.PAIRS


def source_pins():
    pins = dense.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


choose_target, regional_winners = dense.choose_target, dense.regional_winners


def fit_entry(numerators, pairs, weight, full, entry, old_odds, budget, *, max_calls=128):
    """Reuse the dense fitter except for the missing full-count-zero axis."""
    if full != 0:
        return frozen.fit_entry(numerators, pairs, weight, full, entry, old_odds,
                                budget, max_calls=max_calls)
    empty, single, full = frozen.histogram(pairs, weight, full)
    if (type(entry) is not int or not 0 <= entry < 4 or type(max_calls) is not int
            or not 8 <= max_calls <= 128 or empty == 0 or single == 0):
        raise ValueError('Expected a nontrivial full-zero axis and bounded fit')
    bounds = [(-64., 128.)]
    center = max(-64., min(128., log2(single/empty)))
    projected = [max(-64., min(128., float(u))) for u, v in old_odds]
    nearest = sorted(projected, key=lambda u: abs(u-center))[:4]
    seeds = tuple(dict.fromkeys((center, *nearest)))
    best, best_value, calls = None, float('inf'), 0

    def objective(point):
        nonlocal best, best_value, calls
        budget.check()
        if calls >= max_calls:
            raise ArithmeticError('Optimizer exceeded its per-fit call bound')
        calls += 1
        odds = float(point[0]), -64.
        value = float(numerators(odds)[entry]-single*odds[0])
        if not isfinite(value):
            raise ArithmeticError('Nonfinite positive fitting objective')
        if value < best_value:
            best, best_value = odds, value
        return value

    for seed in seeds:
        objective((seed,))
    fitted = minimize(objective, [best[0]], method='Powell', bounds=bounds,
                      options=dict(maxfev=max_calls-calls, xtol=1e-8, ftol=1e-10))
    positive = tuple(F(2.**v) for v in best)
    if not all(v > 0 and isfinite(float(v)) for v in positive):
        raise ArithmeticError('Nonfinite or zero exact positive proposal')
    actual_odds = tuple(log2(float(v)) for v in positive)
    return positive, dict(active_groups=weight, histogram=dict(empty=empty, single=single, full=0),
                          entry=entry, chart='full-zero-positive-y2^-64', chart_bounds=bounds,
                          objective_calls=calls, optimizer_success=bool(fitted.success),
                          optimizer_message=str(fitted.message), positive_objective=best_value,
                          fitted_log_literals=[str(F(v)) for v in best],
                          positive_marker_literals=list(map(str, positive)),
                          evaluated_marker_log_literals=[str(F(v)) for v in actual_odds])

def run(source_path, weight, *, max_seconds=90., max_fits=17, max_calls=128, target_bits=80):
    if (type(weight) is not int or not 2 <= weight < PAIRS
            or not isfinite(max_seconds) or not 0 < max_seconds <= 90
            or type(max_fits) is not int or not 0 <= max_fits <= 32
            or type(max_calls) is not int or not 8 <= max_calls <= 128
            or type(target_bits) is not int or not 0 <= target_bits <= 4096):
        raise ValueError('Expected original weight 2..4095 and bounded time/fits/calls/target')
    path = Path(source_path).resolve()
    raw, pins = path.read_bytes(), source_pins()
    source = json.loads(raw)
    construction = provider.gamma.construction(30)
    if (source.get('active_groups') != weight or source['state_dimension'] != 30
            or F(source['target_distance']) != F(2, 5) or source['construction'] != construction):
        raise ValueError('Source must explicitly calibrate this original gamma30/d40 weight')
    original = source['witness']
    if (original.get('state_dimension', 30) != 30 or original.get('inner_width', 64) != 64
            or int(original.get('field_order', source.get('field_order', 1 << 128))) != 1 << 128):
        raise ValueError('Source field or dimensions disagree')
    positives, marker_key = frozen.literals.positive_markers(original)
    z, alpha = F(original['z']), F(original['alpha'])
    if (not 0 < float(z) < 1 or not 0 < float(alpha) <= 1
            or F(float(z)) != z or F(float(alpha)) != alpha):
        raise ValueError('Fixed initial z/alpha must be valid exact binary64 values')
    grid = original.get('grid', 128)
    if grid not in (64, 128):
        raise ValueError('This memory-bounded adapter supports grids 64 and 128')
    for name, digest in source.get('source_sha256', {}).items():
        if name in pins and pins[name] != digest:
            raise ValueError('Source conflict: '+name)
        pins[name] = digest
    def check():
        if path.read_bytes() != raw:
            raise RuntimeError('Source witness changed')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Helper source changed: '+name)
    check()
    initial_count = len(positives)
    odds = [tuple(log2(float(v)) for v in row) for row in positives]
    full = provider.alias.full_counts(PAIRS, weight)
    started, budget = monotonic(), provider.old.Budget(max_seconds)
    cache, history, fits, attempts = {}, [], [], {}
    coeff, counts, coefficient_hash, initial_metadata, completed_winners = None, None, None, None, None
    pending, status = dict(stage='initial-bank'), 'fit-limit'
    def numerators(point):
        budget.check()
        point = tuple(map(float, point))
        if point not in cache:
            cache[point] = np.asarray(provider.alias.power(coeff.pair_logs(point, _budget=budget), 64)).ravel()
        return cache[point]
    def complete(stage):
        nonlocal pending, completed_winners
        pending = dict(stage='all-profile-sum', fit_index=stage, marker_count=len(positives))
        row = mode.support_sum(PAIRS, weight, regions, counts=counts, budget=budget)
        if row['profile_pairs_evaluated'] != len(full)**2:
            raise ArithmeticError('Incomplete original-weight profile sum')
        row.update(step=len(history), fit_index=stage, active_groups=weight, marker_count=len(positives),
                   log2_outer_prefix=prefix, log2_upper_screen=prefix+row['log2_support_sum_screen'],
                   regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                   selected_reliable_entries=int(reliability.sum()),
                   selected_positive_fallback_entries=int((~reliability).sum()),
                   complete_selected_weight_screen=True, elapsed_seconds=monotonic()-started)
        completed_winners = winners.copy()
        history.append(row)
        pending = None
        print(json.dumps(dict(step=row['step'], fit_index=stage, markers=len(positives),
                              score=row['log2_upper_screen'], largest=row['largest_profile'],
                              elapsed_seconds=row['elapsed_seconds'])), flush=True)
    try:
        coeff = provider.AliasCoefficients(float(z), float(alpha), state=30, _budget=budget)
        coefficient_hash = sha256(coeff.composed.tobytes()).hexdigest()
        bank = provider.AliasBankScreen(float(z), float(alpha), odds, 0., state=30, distance=F(2, 5),
                    grid=grid, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
        initial_metadata = bank.metadata()
        regions, reliability, winners = regional_winners(bank, weight, budget)
        prefix = bank.constant+bank.alpha*bank.density*weight
        del bank
        counts = mode.MidpointCounts(PAIRS, weight, max_singletons=64, budget=budget)
        complete(None)
        for number in range(max_fits):
            if history[-1]['log2_upper_screen'] <= -target_bits:
                status = 'FLOAT-target-met'
                break
            if len(positives) >= 128:
                status = 'marker-limit'
                break
            target, entry = choose_target(history[-1]['largest_profile'], attempts, pairs=PAIRS, weight=weight)
            pending = dict(stage='fit', fit_index=number, target_full=target, entry=entry,
                           marker_count=len(positives))
            positive, record = fit_entry(numerators, PAIRS, weight, target, entry, odds, budget,
                                                max_calls=max_calls)
            attempts[target, entry] = attempts.get((target, entry), 0)+1
            record.update(fit_index=number, selected_from_complete_step=history[-1]['step'],
                          selected_from_profile=dict(history[-1]['largest_profile']))
            fits.append(record)
            if positive in positives:
                record['addition_status'] = 'duplicate-not-added'
                pending = None
                continue
            record['addition_status'] = 'pending-atlas'
            pending['stage'] = 'fresh-single-marker-atlas'
            fresh_odds = tuple(log2(float(v)) for v in positive)
            bank = provider.AliasBankScreen(float(z), float(alpha), [fresh_odds], 0., state=30, distance=F(2, 5),
                        grid=grid, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
            added, reliable = frozen.regional_with_selection(bank, weight, budget)
            del bank
            improved = added < regions
            reliability[improved], winners[improved] = reliable[improved], len(positives)
            np.minimum(regions, added, out=regions)
            positives.append(positive)
            odds.append(fresh_odds)
            record['addition_status'] = 'applied-awaiting-complete-sum'
            complete(number)
            record['addition_status'] = 'complete-scored-addition'
        if history[-1]['log2_upper_screen'] <= -target_bits:
            status = 'FLOAT-target-met'
        elif len(positives) == 128:
            status = 'marker-limit'
    except (TimeoutError, MemoryError) as error:
        status = 'timeout' if isinstance(error, TimeoutError) else 'memory-error'
    if coeff is not None and sha256(coeff.composed.tobytes()).hexdigest() != coefficient_hash:
        raise RuntimeError('Shared fresh coefficient table changed')
    check()
    completed_count = history[-1]['marker_count'] if history else initial_count
    retained = sorted(set(map(int, completed_winners.ravel()))) if completed_winners is not None else None
    return dict(schema='gamma-original-weight-adaptive-FLOAT-v1', status=status,
                construction=construction, state_dimension=30, target_distance='2/5', active_groups=weight,
                source=dict(path=str(path), sha256=sha256(raw).hexdigest()), source_sha256=dict(sorted(pins.items())),
                source_marker_key=marker_key, source_role='exact literals only; all initial numerical caps rebuilt',
                initial_marker_count=initial_count, marker_count=completed_count, applied_marker_count=len(positives),
                witness=dict(z=str(z), alpha=str(alpha), markers=[list(map(str, p)) for p in positives[:completed_count]],
                             marker_log_odds_exact=[[str(F(v)) for v in p] for p in odds[:completed_count]],
                             grid=grid, state_dimension=30),
                unscored_applied_markers=[list(map(str, p)) for p in positives[completed_count:]],
                witness_scope='Exact marker prefix of last complete sum; no incomplete addition exported in witness',
                feasible_full_range=[int(full[0]), int(full[-1])], count_complement_only=True,
                original_weight_used_for_histograms_markers_matrices_and_prefix=True,
                initial_bank_metadata=initial_metadata, coefficient_table_sha256=coefficient_hash,
                fresh_coefficients_once=coeff is not None, counts_built=int(counts is not None),
                history=history, fits=fits, last_complete_step=history[-1] if history else None, pending_step=pending,
                winning_marker_indices_at_last_complete=retained,
                possible_winner_only_count=len(retained) if retained is not None else None,
                pruning_performed=False, pruning_requires_fresh_original_weight_replay=True,
                max_seconds=max_seconds, max_fits=max_fits, max_calls_per_fit=max_calls, max_markers=128,
                target_bits=target_bits, distinct_positive_points=len(cache), elapsed_seconds=monotonic()-started,
                common_potential_used=False, old_numerical_bounds_reused=False, no_weight_interpolation=True,
                directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                retained_numerical_storage='One coefficient table, one initial atlas bank then one fresh atlas at a time, O(r) regional minima and midpoint counts',
                scope='Complete original-weight FLOAT sums only; no convergence or certificate claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--weight', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=90.)
    parser.add_argument('--max-fits', type=int, default=17)
    parser.add_argument('--max-calls', type=int, default=128)
    parser.add_argument('--target-bits', type=int, default=80)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = run(args.source, args.weight, max_seconds=args.max_seconds, max_fits=args.max_fits,
                 max_calls=args.max_calls, target_bits=args.target_bits)
    print(json.dumps(dict(status=result['status'], artifact=provider.search.write_new(args.output, result))), flush=True)


if __name__ == '__main__':
    main()
