"""Bounded same-profile FLOAT retune after changing the return cap.

Every trial retains all97 original/centered markers. At most24 trial values;
the final result is not an all-profile or directed distance certificate.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite, log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_reversed_return_v44 as prior

HERE = Path(__file__).resolve().parent
PRIOR_SHA = '485685eec3829eb01d4821b0c87cfc98cc3fc27b9ecd8328f4a827c05444341f'
LIMIT, RESERVE = 32 << 10, 20 << 30


def source_pins():
    pins = prior.source_pins()
    if pins[Path(prior.__file__).name] != PRIOR_SHA:
        raise RuntimeError('V44 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(source_name, output_name):
    with prior.prior.control.remote_limits():
        return _run(source_name, output_name)


def _run(source_name, output_name):
    import numpy as np
    from scipy.optimize import minimize
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    control, old = prior.prior.control, prior.prior.control.old
    source, output = old.prior.previous.inside(source_name), old.prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20 or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Storage/input guard')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != old.prior.SOURCE_SHA:
        raise ValueError('Pinned V34 required')
    row = old.prior.validate_source(json.loads(raw))
    w, profile = row['witness'], row['largest_profile']
    z0, a0 = float(F(w['z'])), float(F(w['alpha']))
    markers, _, _ = control.marker_extension(w['markers'], a0)
    odds = [tuple(log2(float(F(v))) for v in pair) for pair in markers]
    density = gamma.old.outer_log_density(1 << 128, 'invertible-toeplitz')
    started, budget, history, cache = monotonic(), gamma.old.Budget(80.), [], {}
    def objective(point):
        budget.check()
        key = tuple(map(float, point))
        if key in cache:
            return cache[key]
        if len(history) >= 24:
            raise RuntimeError('Exceeded24 trial cap')
        # Local chart: multiplicative changes in -log2(z) and alpha.
        if not all(isfinite(x) and abs(x) <= 1.5 for x in key):
            return 1.e9
        lz, alpha = log2(z0)*2.**key[0], a0*2.**key[1]
        z = 2.**lz
        coeff = old.Coefficients(prior.local_logs(z), z, alpha, budget)
        best = control.Best()
        for i, point_odds in enumerate(odds):
            atlas = alias.alias_table(z, alpha, point_odds, _coefficients=coeff, _budget=budget)
            value, reliable = old.selected(atlas, point_odds)
            best.add(i, value, reliable)
        prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096, density)
        score = prefix+profile['log2_count_upper']+old.prior.profile_transfer(best.logs, best.logs)
        if not isfinite(score):
            raise ArithmeticError('Nonfinite objective')
        history.append(dict(z=str(F(z)), alpha=str(F(alpha)), score=score, chart=list(key),
            regional_logs=best.logs.tolist(), reliable=best.reliable.tolist(), winners=best.winners.tolist()))
        cache[key] = score
        return score
    baseline = objective((0., 0.))
    if abs(baseline-4009.084049851517) > 1e-4:
        raise ArithmeticError('V44 baseline differs')
    answer = minimize(objective, np.zeros(2), method='Nelder-Mead', options=dict(
        maxfev=24, initial_simplex=np.array([[0., 0.], [.1, 0.], [0., .1]]),
        xatol=1.e-5, fatol=1.e-3))
    best = min(history, key=lambda item: item['score'])
    if source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Source changed')
    result = dict(schema='reversed-return-retune-v45-FLOAT', input_sha256=sha256(raw).hexdigest(),
        source_sha256=pins, markers=markers, grid=128, active_groups=4096,
        selected_profile=[2026, 2026], target_distance='9/20', history=history, best=best,
        trials=len(history), optimizer_success=bool(answer.success), optimizer_message=str(answer.message),
        log2_profile_count=profile['log2_count_upper'],
        principal_bad_gamma_upper=str(F(8*comb(64, 32), (1 << 128)-(1 << 64))),
        construction='A128 emission,A64-transpose feedback,state32,global random degree16gamma',
        original_counts_once=True, actual_returns_retained=True, markers_fixed_across_trials=True,
        conditional_scores_exclude_setup_risk=True, directed=False, all_profiles_replayed=False,
        whole_code_certificate=False, elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Storage cap')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, sha256=sha256(data).hexdigest(), bytes=len(data),
                trials=len(history), best=best, baseline=baseline, elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
