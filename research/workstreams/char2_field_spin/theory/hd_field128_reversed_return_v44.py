"""Use the involution to bound return ranks in both directions (FLOAT).

L=I+gamma AC, CA=0 imply L^2=I and CL=C. Swapping input and output
interchanges active-zero set D with unused set V, preserving the syndrome.
Only rank is transferred; all probabilities use the original input size.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_dense_feedback_return_v42 as prior

HERE = Path(__file__).resolve().parent
PRIOR_SHA = '4e054473298279dd2c9db148b850144cb9b76c51930f461c96f6b53217857578'
LIMIT, Q = 32 << 10, 1 << 128


def rank_floor(j, c, e, *, state=32, degree=16):
    n = 2*state
    if not (0 <= c <= j <= n and 0 <= e <= n-j and degree >= 1):
        raise ValueError('Invalid zero/support counts')
    t, u, v = j-c, min(state, j-c), n-j-e
    if e >= state:
        return u
    h = max(0, e+u-state, min((e+1)//2, t), min(e, (t+1)//2))
    loss = min(c, v, state-h)//degree if h else 0
    return max(0, e+u-state, h-loss)


def probability(j, c, e, q, *, state=32, degree=16):
    u = min(state, j-c)
    delta = u-rank_floor(j, c, e, state=state, degree=degree)
    return F((q**delta-1)*q**(j-c-u), q**j)


def local_logs(z):
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    logs = prior.local_logs(z)
    lz, lw, log_d = log2(z), log2(1-z), 4096+gamma.old.log1mexp2(-4096)
    for j in range(1, 65):
        terms = []
        for c in range(j):
            u = min(32, j-c)
            for e in range(min(32, 65-j)):
                delta = u-rank_floor(j, c, e)
                if delta:
                    lp = (delta-c-u)*128+gamma.old.log1mexp2(-128*delta)
                    terms.append(log2(comb(j, c)*comb(64-j, e))+(64-c-e)*lz+(c+e)*lw+lp)
        cap = float(np.logaddexp2.reduce(terms)) if terms else -np.inf
        logs[j, 1, 0] = min(logs[j, 1, 0], cap-log_d)
    return logs


def source_pins():
    pins = prior.source_pins()
    if pins[Path(prior.__file__).name] != PRIOR_SHA:
        raise RuntimeError('V42 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(source_name, output_name):
    with prior.control.remote_limits():
        return _run(source_name, output_name)


def _run(source_name, output_name):
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    control, old = prior.control, prior.control.old
    source, output = old.prior.previous.inside(source_name), old.prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20 or shutil.disk_usage(HERE).free-LIMIT < prior.RESERVE:
        raise ValueError('Storage/input guard')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != old.prior.SOURCE_SHA:
        raise ValueError('Pinned V34 required')
    row = old.prior.validate_source(json.loads(raw))
    witness, profile = row['witness'], row['largest_profile']
    z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
    markers, _, _ = control.marker_extension(witness['markers'], alpha)
    prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096,
        gamma.old.outer_log_density(Q, 'invertible-toeplitz'))
    if abs(prefix-row['log2_outer_prefix']) > 1e-7 or (profile['even_full'], profile['odd_full']) != (2026, 2026):
        raise ArithmeticError('Original-r/profile guard')
    offset = prefix+profile['log2_count_upper']
    started, budget, rows = monotonic(), gamma.old.Budget(80.), []
    for name, logs in [('V42', prior.local_logs(z)), ('reversed', local_logs(z))]:
        coeff, best = old.Coefficients(logs, z, alpha, budget), control.Best()
        for i, literal in enumerate(markers):
            budget.check()
            odds = tuple(log2(float(F(x))) for x in literal)
            atlas = alias.alias_table(z, alpha, odds, _coefficients=coeff, _budget=budget)
            value, reliable = old.selected(atlas, odds)
            best.add(i, value, reliable)
        rows.append(dict(name=name, score=offset+old.prior.profile_transfer(best.logs, best.logs),
            regional_logs=best.logs.tolist(), reliable=best.reliable.tolist(), winners=best.winners.tolist(),
            local_return=logs[1:, 1, 0].tolist()))
    if abs(rows[0]['score']-4057.940908270539) > 1e-4:
        raise ArithmeticError('V42 replay mismatch')
    if source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Source changed')
    result = dict(schema='reversed-return-v44-FLOAT', source_sha256=pins,
        input_sha256=sha256(raw).hexdigest(), witness=witness, markers=markers,
        active_groups=4096, selected_profile=[2026, 2026], target_distance='9/20', rows=rows,
        principal_bad_gamma_upper=str(F(8*comb(64, 32), Q-(1 << 64))),
        additional_gamma_exclusions=False, original_input_normalization=True,
        rank_reversal_uses_involution=True, conditional_scores_exclude_setup_risk=True,
        construction='A128 emission,A64-transpose feedback,state32; other stages unchanged',
        directed=False, all_profiles_replayed=False, whole_code_certificate=False,
        elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, separators=(',', ':'), sort_keys=True, allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < prior.RESERVE:
        raise RuntimeError('Output/reserve guard')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, sha256=sha256(data).hexdigest(), bytes=len(data),
                scores=[row['score'] for row in rows], regional_logs=[row['regional_logs'] for row in rows],
                elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
