"""Bounded FLOAT marker enrichment, rebuilding every initial cap from literals."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_gamma_midpoint_counts as mode


provider = mode.near.provider
HERE = Path(__file__).resolve().parent


def source_pins():
    pins = mode.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def positive_markers(witness):
    key = 'exact_positive_marker_proposals' if 'exact_positive_marker_proposals' in witness else 'markers'
    rows = witness.get(key)
    if not isinstance(rows, list) or not 1 <= len(rows) <= 128 or any(len(row) != 2 for row in rows):
        raise ValueError('Expected1..128 explicit positive marker pairs; log-only sources are not accepted')
    markers = [tuple(map(F, row)) for row in rows]
    if not all(value > 0 and isfinite(float(value)) and float(value) > 0 for row in markers for value in row):
        raise ValueError('Markers must be positive and finite in binary64')
    return markers, key


def regional_with_selection(bank, weight, budget):
    regions = bank.regional(weight)
    full = np.arange(weight//2+1)
    single = weight-2*full
    best, reliable_best = np.full((4, len(full)), np.inf), np.zeros((4, len(full)), dtype=bool)
    for odds, table in zip(bank.odds, bank.tables):
        budget.check()
        reliable = table.reliable[:, single % bank.grid, full % bank.grid]
        values = np.where(reliable, table.logs[:, single % bank.grid, full % bank.grid],
                          table.positive_numerator[:, None])-single*odds[0]-full*odds[1]
        improved = values < best
        reliable_best[improved] = reliable[improved]
        np.minimum(best, values, out=best)
    return regions, reliable_best.T.reshape(regions.shape)


def run(source_path, weights, *, target_bits=54, max_seconds=90., max_rounds=12, max_markers=128):
    """Reuse one coefficient table; each completed round evaluates every weight."""
    weights = tuple(weights)
    if (not 1 <= len(weights) <= 8 or len(set(weights)) != len(weights)
            or any(type(r) is not int or not 2 <= r < 4096 for r in weights)):
        raise ValueError('Expected1..8 distinct integer weights in2..4095')
    if (not isfinite(max_seconds) or not 0 < max_seconds <= 90 or type(max_rounds) is not int or not 0 <= max_rounds <= 12
            or type(max_markers) is not int or not 1 <= max_markers <= 128
            or type(target_bits) is not int or not 0 <= target_bits <= 4096):
        raise ValueError('Invalid bounded run limits')
    path = Path(source_path).resolve()
    raw, pins = path.read_bytes(), source_pins()
    saved = json.loads(raw)
    witness = saved['witness']
    state, distance = saved['state_dimension'], F(saved['target_distance'])
    construction = provider.gamma.construction(state)
    if saved.get('construction', construction) != construction or not 0 < distance < F(1, 2):
        raise ValueError('Wrong construction or distance')
    if (witness.get('state_dimension', state) != state or witness.get('inner_width', 64) != 64
            or int(witness.get('field_order', saved.get('field_order', 1 << 128))) != 1 << 128):
        raise ValueError('Witness field or dimensions disagree')
    z, alpha = F(witness['z']), F(witness['alpha'])
    positives, marker_key = positive_markers(witness)
    if not (0 < float(z) < 1 and 0 < float(alpha) <= 1 and len(positives) <= max_markers):
        raise ValueError('Invalid fixed marks or initial marker count')
    grid = witness.get('grid', 128)
    if grid not in (64, 128, 256):
        raise ValueError('Unsupported alias grid')
    for name, digest in saved.get('source_sha256', {}).items():
        if name in pins and pins[name] != digest:
            raise ValueError('Source closure conflict: '+name)
        pins[name] = digest
    def check_pins():
        if path.read_bytes() != raw:
            raise RuntimeError('Source artifact changed')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source closure changed: '+name)
    check_pins()
    odds = [tuple(log2(float(value)) for value in pair) for pair in positives]
    initial_count, initial_odds = len(positives), list(odds)
    started, budget = monotonic(), provider.old.Budget(max_seconds)
    history, fit_history, targets, cache = [], [], [], {}
    pending, status, metadata, coefficient_hash = dict(round=0, stage='initial-bank', completed=[], pending_weights=list(weights)), 'round-limit', None, None
    counts, regions, reliability, prefixes = {}, {}, {}, {}

    def numerators(point):
        budget.check()
        point = tuple(map(float, point))
        if point not in cache:
            cache[point] = np.asarray(provider.alias.power(coeff.pair_logs(point, _budget=budget), 64)).ravel()
        return cache[point]

    def complete_round(number, target):
        nonlocal pending
        pending = dict(round=number, stage='complete-profile-sums', marker_count=len(positives),
                       completed=[], pending_weights=list(weights), target=target)
        for weight in weights:
            budget.check()
            if weight not in counts:
                counts[weight] = mode.MidpointCounts(4096, weight, max_singletons=64, budget=budget)
            result = mode.support_sum(4096, weight, regions[weight], counts=counts[weight], budget=budget)
            result.update(active_groups=weight, log2_outer_prefix=prefixes[weight],
                          log2_upper_screen=prefixes[weight]+result['log2_support_sum_screen'],
                          regional_logs_sha256=sha256(regions[weight].tobytes()).hexdigest(),
                          selected_reliable_entries=int(reliability[weight].sum()),
                          selected_positive_fallback_entries=int((~reliability[weight]).sum()),
                          complete_selected_weight_screen=True)
            pending['completed'].append(result)
            pending['pending_weights'].remove(weight)
        row = dict(pending, elapsed_seconds=monotonic()-started, complete_all_requested_weights=True)
        history.append(row)
        pending = None
        print(json.dumps(dict(round=number, markers=len(positives), target=target,
                              scores={r['active_groups']: r['log2_upper_screen'] for r in row['completed']},
                              elapsed_seconds=row['elapsed_seconds'])), flush=True)
        return max(row['completed'], key=lambda item: item['log2_upper_screen'])

    def fit_target(weight, full):
        single, empty = weight-2*full, 4096-weight+full
        rows, proposed = [], []
        occupancy = (log2(single/empty) if single else -64., log2(full/empty) if full else -64.)
        seeds = [tuple(np.clip(point, -64., 8.)) for point in [occupancy, *odds]]
        seeds = [(x if single else -64., y if full else -64.) for x, y in seeds]
        for entry in range(4):
            best = min(seeds, key=lambda point: numerators(point)[entry]-single*point[0]-full*point[1])
            best_value = float(numerators(best)[entry]-single*best[0]-full*best[1])
            calls = 0
            def objective(raw_point):
                nonlocal best, best_value, calls
                budget.check()
                calls += 1
                point = (float(raw_point[0]), -64.) if not full else (-64., float(raw_point[0])) if not single else tuple(map(float, raw_point))
                value = float(numerators(point)[entry]-single*point[0]-full*point[1])
                if value < best_value:
                    best, best_value = point, value
                return value
            initial = [best[0]] if not full else [best[1]] if not single else best
            fitted = minimize(objective, initial, method='Powell', bounds=[(-64., 8.)]*len(initial),
                              options=dict(maxfev=256, xtol=1e-8, ftol=1e-10))
            positive = tuple(F(2.**value) for value in best)
            rows.append(dict(active_groups=weight, single=single, full=full, entry=entry, calls=calls,
                             fitted_log_literals=list(map(lambda value: str(F(value)), best)),
                             positive_marker_literals=list(map(str, positive)), positive_objective=best_value,
                             optimizer_success=bool(fitted.success)))
            proposed.append(positive)
        return proposed, rows

    try:
        coeff = provider.AliasCoefficients(float(z), float(alpha), state=state, _budget=budget)
        coefficient_hash = sha256(coeff.composed.tobytes()).hexdigest()
        bank = provider.AliasBankScreen(float(z), float(alpha), odds, 0., state=state, distance=distance,
                    grid=grid, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
        metadata = bank.metadata()
        for weight in weights:
            regions[weight], reliability[weight] = regional_with_selection(bank, weight, budget)
            prefixes[weight] = bank.constant+bank.alpha*bank.density*weight
        worst = complete_round(0, None)
        for number in range(1, max_rounds+1):
            if worst['log2_upper_screen'] <= -target_bits:
                status = 'FLOAT-target-met'
                break
            weight, profile = worst['active_groups'], worst['largest_profile']
            candidates = [profile['even_full'], profile['odd_full']]
            full = next((f for f in candidates if (weight, f) not in targets), candidates[0])
            target = dict(active_groups=weight, full=full)
            pending = dict(round=number, stage='fit', target=target, completed=[], pending_weights=list(weights))
            proposed, fits = fit_target(weight, full)
            fit_history.extend(fits)
            fresh = list(dict.fromkeys(pair for pair in proposed if pair not in positives))
            if not fresh or len(positives)+len(fresh) > max_markers:
                status = 'no-new-markers' if not fresh else 'marker-limit'
                pending = None
                break
            pending['stage'] = 'fresh-atlases'
            fresh_odds = [tuple(log2(float(value)) for value in pair) for pair in fresh]
            bank = provider.AliasBankScreen(float(z), float(alpha), fresh_odds, 0., state=state, distance=distance,
                        grid=grid, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
            for weight in weights:
                added, reliable = regional_with_selection(bank, weight, budget)
                improved = added < regions[weight]
                reliability[weight][improved] = reliable[improved]
                regions[weight] = np.minimum(regions[weight], added)
            positives.extend(fresh)
            odds.extend(fresh_odds)
            targets.append((target['active_groups'], full))
            worst = complete_round(number, target)
        if history and max(row['log2_upper_screen'] for row in history[-1]['completed']) <= -target_bits:
            status = 'FLOAT-target-met'
    except TimeoutError:
        status = 'timeout'
    if coefficient_hash is not None and sha256(coeff.composed.tobytes()).hexdigest() != coefficient_hash:
        raise RuntimeError('Fixed coefficient table changed')
    check_pins()
    return dict(schema='gamma-fixed-marks-midpoint-enrichment-FLOAT-v1', construction=construction,
                state_dimension=state, target_distance=str(distance), requested_weights=list(weights), status=status,
                source=dict(path=str(path), sha256=sha256(raw).hexdigest()), source_marker_key=marker_key,
                source_active_groups=saved.get('active_groups'), source_role='exact positive literals only; no saved numerical caps',
                source_sha256=dict(sorted(pins.items())), coefficient_table_sha256=coefficient_hash,
                original_coefficients_reused=coefficient_hash is not None, initial_bank_rebuilt=metadata is not None,
                old_numerical_bounds_reused=False,
                bank_metadata=metadata, initial_marker_count=initial_count, marker_count=len(positives),
                witness=dict(z=str(z), alpha=str(alpha), markers=[list(map(str, pair)) for pair in positives], grid=grid,
                             marker_log_odds_exact=[[str(F(v)) for v in pair] for pair in odds],
                             floor_policy='positive-fallback', state_dimension=state),
                original_marker_log_odds_exact=[[str(F(v)) for v in pair] for pair in initial_odds],
                source_log_witness_literals=witness.get('marker_log_odds_exact'),
                conversion='Saved positive rationals are retained exactly; FLOAT uses log2(float(rational)). New exp2 binary64 outputs are exact rational proposals.',
                rounds=history, last_complete_round=history[-1] if history else None, pending_round=pending,
                pending_weights=list(weights) if pending is None and not history else pending['pending_weights'] if pending else [],
                witness_may_include_incomplete_round=pending is not None, fits=fit_history, distinct_positive_points=len(cache),
                target_bits=target_bits, max_seconds=max_seconds, max_rounds=max_rounds, max_markers=max_markers,
                optimizer='Powell on positive numerator minus histogram log tilt', optimizer_bounds=[-64., 8.],
                max_objective_calls_per_fit=256, elapsed_seconds=monotonic()-started,
                common_potential_used=False, all_profiles_retained_for_completed_weights=True,
                directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                scope='FLOAT calibration of listed weights only; no interpolated or contiguous coverage claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--weights', type=int, nargs='+', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--max-seconds', type=float, default=90.)
    parser.add_argument('--max-rounds', type=int, default=12)
    parser.add_argument('--max-markers', type=int, default=128)
    parser.add_argument('--target-bits', type=int, default=54)
    args = parser.parse_args()
    if Path(args.output).exists():
        raise FileExistsError(args.output)
    result = run(args.source, args.weights, target_bits=args.target_bits, max_seconds=args.max_seconds,
                 max_rounds=args.max_rounds, max_markers=args.max_markers)
    pin = provider.search.write_new(args.output, result)
    print(json.dumps(dict(status=result['status'], artifact=pin, marker_count=result['marker_count'])), flush=True)


if __name__ == '__main__':
    main()
