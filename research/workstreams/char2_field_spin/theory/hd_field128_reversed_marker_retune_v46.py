"""Same-profile FLOAT joint retune of z,alpha and one coefficient marker.

At most100 trials; final97-marker bank plus the optimized marker. No changes
to construction/caps, and no all-profile or directed certificate claim.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_reversed_retune_v45 as previous

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = '55bb7bfa70deea29c3cddc07043393e130c662665feaae636dfbd3492a24e260'
INPUT_SHA = 'a51b6d2fcf061f59b071a5d3f0f8d2e0f1cb79b4337a0f0aa3f02b63b5aa9f43'
LIMIT, RESERVE = 32 << 10, 20 << 30


def source_pins():
    pins = previous.source_pins()
    if pins[Path(previous.__file__).name] != PREVIOUS_SHA:
        raise RuntimeError('V45 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(source_name, output_name):
    with previous.prior.prior.control.remote_limits():
        return _run(source_name, output_name)


def _run(source_name, output_name):
    import numpy as np
    from scipy.optimize import minimize
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    control = previous.prior.prior.control
    old = control.old
    source, output = old.prior.previous.inside(source_name), old.prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > LIMIT or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Storage/input guard')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Pinned V45 required')
    artifact = json.loads(raw)
    if artifact['selected_profile'] != [2026, 2026] or artifact['active_groups'] != 4096:
        raise ValueError('Original profile required')
    seed = artifact['best']
    z0, a0 = float(F(seed['z'])), float(F(seed['alpha']))
    index = seed['winners'][0][0]
    if seed['winners'] != [[index, index], [index, index]]:
        raise ValueError('Single seed marker required')
    original_odds = tuple(log2(float(F(v))) for v in artifact['markers'][index])
    density = gamma.old.outer_log_density(1 << 128, 'invertible-toeplitz')
    count = artifact['log2_profile_count']
    started, budget, history, coefficients = monotonic(), gamma.old.Budget(80.), [], {}
    def evaluate(point):
        budget.check()
        if len(history) >= 100:
            raise RuntimeError('Trial cap')
        if not all(isfinite(x) and abs(x) <= 4 for x in point):
            return 1.e9
        z, alpha = 2.**(log2(z0)*2.**point[0]), a0*2.**point[1]
        odds = (original_odds[0]+point[2], original_odds[1]+point[3])
        key = float(z), float(alpha)
        if key not in coefficients:
            if len(coefficients) >= 8:
                coefficients.pop(next(iter(coefficients)))
            coefficients[key] = old.Coefficients(previous.prior.local_logs(z), z, alpha, budget)
        atlas = alias.alias_table(z, alpha, odds, _coefficients=coefficients[key], _budget=budget)
        logs, reliable = old.selected(atlas, odds)
        prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096, density)
        score = prefix+count+old.prior.profile_transfer(logs, logs)
        history.append(dict(chart=[float(x) for x in point], score=float(score)))
        return float(score)
    baseline = evaluate(np.zeros(4))
    if abs(baseline-seed['score']) > 1e-4:
        raise ArithmeticError('V45 replay mismatch')
    simplex = np.zeros((5, 4))
    for k, step in enumerate((.03, .03, .1, .1)):
        simplex[k+1, k] = step
    answer = minimize(evaluate, np.zeros(4), method='Nelder-Mead', options=dict(
        maxfev=98, initial_simplex=simplex, xatol=1.e-5, fatol=1.e-3))
    selected = min(history, key=lambda row: row['score'])
    point = selected['chart']
    z, alpha = 2.**(log2(z0)*2.**point[0]), a0*2.**point[1]
    new_marker = [str(F(2.**(original_odds[k]+point[k+2]))) for k in range(2)]
    markers = artifact['markers']+[new_marker]
    coeff = old.Coefficients(previous.prior.local_logs(z), z, alpha, budget)
    best = control.Best()
    for i, pair in enumerate(markers):
        odds = tuple(log2(float(F(v))) for v in pair)
        atlas = alias.alias_table(z, alpha, odds, _coefficients=coeff, _budget=budget)
        logs, reliable = old.selected(atlas, odds)
        best.add(i, logs, reliable)
    prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096, density)
    final = prefix+count+old.prior.profile_transfer(best.logs, best.logs)
    if source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Source changed')
    result = dict(schema='reversed-marker-retune-v46-FLOAT', source_sha256=pins,
        source=dict(path=source_name, sha256=INPUT_SHA), baseline=baseline,
        best=dict(z=str(F(z)), alpha=str(F(alpha)), score=final, new_marker=new_marker,
                  regional_logs=best.logs.tolist(), winners=best.winners.tolist(), reliable=best.reliable.tolist()),
        markers=markers, history=history, trials=len(history), optimizer_success=bool(answer.success),
        active_groups=4096, selected_profile=[2026, 2026], target_distance='9/20',
        same_valid_V44_return_caps=True, original_count_and_prefix_once=True,
        directed=False, all_profiles_replayed=False, whole_code_certificate=False,
        elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Storage cap')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, sha256=sha256(data).hexdigest(), bytes=len(data),
                trials=len(history), baseline=baseline, best=result['best'], elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
