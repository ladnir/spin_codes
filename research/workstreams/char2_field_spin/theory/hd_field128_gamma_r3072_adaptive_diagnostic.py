"""Bounded, fixed-mark FLOAT marker enrichment at one selected weight.

Saved regional caps are retained and improved by fresh four-marker atlases.
Every completed round evaluates all profiles. No coefficient certificate or
probability bound is inherited from FLOAT arithmetic.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_gamma_midpoint_counts as mode


screen = mode.near.provider
HERE = Path(__file__).resolve().parent
SOURCE = HERE/'hd_field128_gamma30_d40_r3072_fresh_midpoint_float_20261010.json'
FIRST = HERE/'hd_field128_gamma_r3072_f1463_four_marker_diagnostic_20261010.json'
OUTPUT = HERE/'hd_field128_gamma30_d40_r3072_adaptive_midpoint_float_20261010.json'


def run():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    sources = []
    for path in (SOURCE, FIRST):
        raw = path.read_bytes()
        sources.append((path, raw, json.loads(raw)))
    saved, first = sources[0][2], sources[1][2]
    pins = mode.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    for name, value in saved['source_sha256'].items():
        if sha256((HERE/name).read_bytes()).hexdigest() != value:
            raise RuntimeError('Saved source closure changed: '+name)
        if name in pins and pins[name] != value:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = value
    weight, state, distance = 3072, 30, F(2, 5)
    z, alpha = F(saved['witness']['z']), F(saved['witness']['alpha'])
    if (saved['active_groups'], saved['state_dimension'], saved['target_distance']) != (weight, state, str(distance)):
        raise ValueError('Wrong saved construction or selected weight')
    if (first['z'], first['alpha'], first['regional_full']) != (str(z), str(alpha), 1463):
        raise ValueError('Initial diagnostic differs from retained marks')
    regions = np.array(saved['regions'], dtype=float)
    if regions.shape != (1537, 2, 2) or sha256(regions.tobytes()).hexdigest() != saved['regional_logs_sha256']:
        raise ValueError('Saved regional table does not match its content pin')
    odds = [tuple(float(F(x)) for x in row) for row in saved['witness']['marker_log_odds_exact']]
    positives = [tuple(F(2.**x) for x in row) for row in odds]
    if len(odds) != 76:
        raise ValueError('Expected original 76-marker bank')
    started, budget = monotonic(), screen.old.Budget(180.)
    coeff = screen.AliasCoefficients(float(z), float(alpha), state=state, _budget=budget)
    counts = mode.MidpointCounts(4096, weight, max_singletons=64, budget=budget)
    prefix = float(saved['log2_outer_prefix'])
    fits = list(saved['fits'])
    cache, history, targets = {}, [], []
    pending, status = None, 'round-limit'
    # 0: inherited saved cap (reliability identity was not stored); 1: a fresh
    # reliable alias won; 2: a fresh positive fallback won.
    origin = np.zeros(regions.shape, dtype=np.int8)

    def numerators(point):
        budget.check()
        key = tuple(map(float, point))
        if key not in cache:
            cache[key] = np.asarray(screen.alias.power(coeff.pair_logs(key, _budget=budget), 64)).ravel()
        return cache[key]

    def complete_round(index, target, new_fits):
        result = mode.support_sum(4096, weight, regions, counts=counts, budget=budget)
        result['log2_upper_screen'] = prefix+result['log2_support_sum_screen']
        result.update(round=index, target_full=target, marker_count=len(odds),
                      new_fits=new_fits, elapsed_seconds=monotonic()-started,
                      regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                      retained_saved_cap_entries=int((origin == 0).sum()),
                      selected_fresh_reliable_entries=int((origin == 1).sum()),
                      selected_fresh_positive_entries=int((origin == 2).sum()))
        history.append(result)
        print(json.dumps(dict(round=index, markers=len(odds), target=target,
                              log2_upper_screen=result['log2_upper_screen'],
                              worst=result['largest_profile'], elapsed=result['elapsed_seconds'])), flush=True)
        return result

    def fit_target(full):
        single, empty = weight-2*full, 4096-weight+full
        output, proposed = [], []
        for entry in range(4):
            budget.check()
            entry_fits = [row for row in fits if row['entry'] == entry]
            nearest = sorted(entry_fits, key=lambda row: abs(row['full']-full))[:4]
            seeds = []
            for row in nearest:
                x, y = [float(F(v)) for v in row['marker_logs_exact']]
                sf, ff, ef = row['single'], row['full'], 4096-weight+row['full']
                if single and sf:
                    x += log2(single/empty/(sf/ef))
                if full and ff:
                    y += log2(full/empty/(ff/ef))
                seeds.append((x if single else -64., y if full else -64.))
            seeds.append((log2(single/empty) if single else -64., log2(full/empty) if full else -64.))
            best = min(seeds, key=lambda point: numerators(point)[entry]-single*point[0]-full*point[1])
            best_value = float(numerators(best)[entry]-single*best[0]-full*best[1])
            calls = 0
            def objective(raw):
                nonlocal best, best_value, calls
                budget.check()
                calls += 1
                point = (float(raw[0]), -64.) if not full else (-64., float(raw[0])) if not single else tuple(map(float, raw))
                value = float(numerators(point)[entry]-single*point[0]-full*point[1])
                if value < best_value:
                    best, best_value = point, value
                return value
            initial = [best[0]] if not full else [best[1]] if not single else best
            result = minimize(objective, initial, method='Powell', bounds=[(-64., 8.)]*len(initial),
                              options=dict(maxfev=256, xtol=1e-8, ftol=1e-10))
            positive = tuple(F(2.**x) for x in best)
            row = dict(entry=entry, full=full, single=single,
                       marker_logs_exact=[str(F(x)) for x in best],
                       marker_positive_rationals=list(map(str, positive)), objective=best_value,
                       calls=calls, optimizer_success=bool(result.success))
            output.append(row)
            proposed.append(positive)
        return proposed, output

    try:
        result = complete_round(0, None, [])
        for index in range(1, 13):
            if result['log2_upper_screen'] <= -40.:
                status = 'FLOAT-target-met'
                break
            budget.check()
            pending = dict(round=index, stage='fit')
            if index == 1:
                target = 1463
                proposed = [tuple(map(F, row)) for row in first['exact_new_positive_markers']]
                new_fits = [dict(row, full=target, single=weight-2*target) for row in first['fits']]
            else:
                worst = result['largest_profile']
                candidates = [worst['even_full'], worst['odd_full']]
                target = next((f for f in candidates if f not in targets), candidates[0])
                proposed, new_fits = fit_target(target)
            fresh = []
            for marker in proposed:
                if marker not in positives and marker not in fresh:
                    fresh.append(marker)
            if not fresh:
                status = 'no-new-markers'
                break
            if len(odds)+len(fresh) > 128:
                status = 'marker-limit'
                break
            pending.update(stage='atlas', target_full=target, marker_count=len(fresh))
            fresh_odds = [tuple(log2(float(x)) for x in marker) for marker in fresh]
            bank = screen.AliasBankScreen(float(z), float(alpha), fresh_odds, 0., state=state,
                        distance=distance, grid=128, floor_policy='positive-fallback',
                        _coefficients=coeff, _budget=budget)
            new_regions = bank.regional(weight)
            if abs(bank.constant+bank.alpha*bank.density*weight-prefix) > 1e-8:
                raise ArithmeticError('Fixed outer prefix changed')
            # Independently retain the identity of each freshly selected cap.
            full = np.arange(1537)
            single = weight-2*full
            best_logs = np.full((4, len(full)), np.inf)
            reliability = np.zeros(best_logs.shape, dtype=bool)
            # bank.regional already normalized all entries identically; its
            # winner can be recovered before that common normalizer.
            for mark, table in zip(fresh_odds, bank.tables):
                rel = table.reliable[:, single % bank.grid, full % bank.grid]
                logs = np.where(rel, table.logs[:, single % bank.grid, full % bank.grid],
                                table.positive_numerator[:, None])-single*mark[0]-full*mark[1]
                win = logs < best_logs
                reliability[win] = rel[win]
                best_logs = np.minimum(best_logs, logs)
            reliability = reliability.T.reshape(regions.shape)
            improved = new_regions < regions
            origin[improved] = np.where(reliability[improved], 1, 2)
            regions = np.minimum(regions, new_regions)
            odds.extend(fresh_odds)
            positives.extend(fresh)
            fits.extend(new_fits)
            targets.append(target)
            pending['stage'] = 'complete-profile-sum'
            result = complete_round(index, target, new_fits)
            pending = None
        else:
            if result['log2_upper_screen'] <= -40.:
                status = 'FLOAT-target-met'
    except TimeoutError:
        status = 'timeout'
    # A timed-out incomplete round does not inherit the preceding full-weight
    # result. Export its improved caps as proposals only, plus the completed
    # round's independently recorded numerical result and marker prefix.
    for name, value in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != value:
            raise RuntimeError('Source changed: '+name)
    for path, raw, _ in sources:
        if path.read_bytes() != raw:
            raise RuntimeError('Input artifact changed')
    witness = dict(saved['witness'])
    witness.update(marker_log_odds=[list(row) for row in odds],
                   marker_log_odds_exact=[[str(F(v)) for v in row] for row in odds],
                   exact_positive_marker_proposals=[list(map(str, row)) for row in positives])
    output = dict(schema='gamma30-r3072-adaptive-midpoint-float-v1', active_groups=weight,
                  state_dimension=state, target_distance=str(distance), status=status,
                  source=[dict(path=str(path), sha256=sha256(raw).hexdigest()) for path, raw, _ in sources],
                  source_sha256=pins, witness=witness, construction=saved['construction'],
                  original_marker_count=76, marker_count=len(odds), target_history=targets,
                  fixed_original_z_alpha=True, log2_outer_prefix=prefix, calibration_history=history,
                  last_completed=history[-1] if history else None, pending=pending,
                  regions=regions.tolist(), regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                  proposal_regions_may_include_pending_round=pending is not None,
                  max_seconds=180, max_rounds=12, max_markers=128, elapsed_seconds=monotonic()-started,
                  complete_rounds=len(history), distinct_positive_fit_points=len(cache),
                  source_role='saved FLOAT regional caps retained, fresh atlases combined entrywise; no directed result reused',
                  marker_proposal_rule='binary64 exp2 outputs are exact positive rational proposals, not exact exponential identities',
                  all_profiles_retained=True, common_potential_used=False, directed=False,
                  selected_weight_certificate=False, whole_code_certificate=False)
    pin = screen.search.write_new(OUTPUT, output)
    print(json.dumps(dict(stage='complete', status=status, artifact=pin,
                          last_completed=history[-1]['log2_upper_screen'] if history else None)), flush=True)


if __name__ == '__main__':
    run()
