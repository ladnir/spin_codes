"""Conditional FLOAT return-cap screen for A128 / A64-transpose feedback.

Uses a proved half-block Cauchy rank floor, not generic transversality.
The setup must satisfy every active principal determinant. No certificate
or timing claim is made. Frozen local providers remain unchanged.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_marker_saddles_v40 as control

HERE = Path(__file__).resolve().parent
CONTROL_SHA = '415007705114836ae00da11fd46c89e91f489b3b5ff97a0176775b7c48a9b5dc'
LIMIT, RESERVE, Q = 32 << 10, 20 << 30, 1 << 128


def rank_floor(j, c, e, *, state=32, degree=16):
    """Uniform over all half splits, conditional on principal nonsingularity."""
    n = 2*state
    if not (0 <= c <= j <= n and 0 <= e <= n-j and degree >= 1):
        raise ValueError('Invalid zero/support counts')
    t, u = j-c, min(state, j-c)
    if e >= state:
        return u
    h = max(0, e+u-state, min((e+1)//2, t), min(e, (t+1)//2))
    loss = min(c, n-c, state-h)//degree if h else 0
    return max(0, e+u-state, h-loss)


def probability(j, c, e, q, *, state=32, degree=16):
    u = min(state, j-c)
    delta = u-rank_floor(j, c, e, state=state, degree=degree)
    return F((q**delta-1)*q**(j-c-u), q**j)


def local_logs(z):
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    logs = gamma.local_logs(z, Q, state=32).copy()
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
    pins = control.source_pins()
    if pins[Path(control.__file__).name] != CONTROL_SHA:
        raise RuntimeError('Frozen V40 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(source_name, output_name):
    with control.remote_limits():
        return _run(source_name, output_name)


def _run(source_name, output_name):
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    old = control.old
    source, output = old.prior.previous.inside(source_name), old.prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20:
        raise ValueError('Exclusive output/bounded input required')
    if shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise RuntimeError('20GiB reserve required')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != old.prior.SOURCE_SHA:
        raise ValueError('Pinned v34 result required')
    row = old.prior.validate_source(json.loads(raw))
    witness, profile = row['witness'], row['largest_profile']
    if (profile['even_full'], profile['odd_full']) != (2026, 2026):
        raise ValueError('Selected profile mismatch')
    started, budget = monotonic(), gamma.old.Budget(80.)
    z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
    before = gamma.local_logs(z, Q, state=32)
    after = local_logs(z)
    if not np.array_equal(before[:, 0, :], after[:, 0, :]) or not np.array_equal(before[:, 1, 1], after[:, 1, 1]):
        raise ArithmeticError('Only return may change')
    coefficients = [old.Coefficients(logs, z, alpha, budget) for logs in (before, after)]
    best = [control.Best(), control.Best()]
    markers, labels, centers = control.marker_extension(witness['markers'], alpha)
    for i, literal in enumerate(markers):
        odds = tuple(log2(float(F(v))) for v in literal)
        for k, producer in enumerate(coefficients):
            budget.check()
            atlas = alias.alias_table(z, alpha, odds, _coefficients=producer, _budget=budget)
            logs, reliable = old.selected(atlas, odds)
            best[k].add(i, logs, reliable)
    prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096,
        gamma.old.outer_log_density(Q, 'invertible-toeplitz'))
    if abs(prefix-row['log2_outer_prefix']) > 1e-7:
        raise ArithmeticError('Original-r prefix mismatch')
    offset = prefix+profile['log2_count_upper']
    scores = [offset+old.prior.profile_transfer(item.logs, item.logs) for item in best]
    if abs(scores[0]-old.EXPECTED_BASELINE) > 1e-4:
        raise ArithmeticError('Baseline mismatch')
    result = dict(schema='dense-feedback-return-v42-FLOAT', source=dict(path=source_name, sha256=sha256(raw).hexdigest()),
        source_sha256=pins, construction='A128 emission, A64-transpose feedback, random degree16 global gamma',
        principal_nonsingularity_required=True, state=32, target_distance='9/20', active_groups=4096,
        selected_profile=[2026, 2026], witness=witness, marker_literals=markers,
        scores=scores, regional_logs=[x.logs.tolist() for x in best], winners=[x.winners.tolist() for x in best],
        reliable=[x.reliable.tolist() for x in best], local_return_before=before[1:, 1, 0].tolist(),
        local_return_after=after[1:, 1, 0].tolist(),
        rank_method='one-half scaled Cauchy; worst split; degree16 bordered-minor multiplicity',
        log2_outer_prefix=prefix, log2_profile_count=profile['log2_count_upper'],
        counts_original_and_once=True, actual_return_retained=True, conditional_formula=True,
        setup_risk_added=False, all_profiles_replayed=False, directed=False,
        whole_code_certificate=False, encoder_implemented=False, elapsed_seconds=monotonic()-started)
    if source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Source changed')
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Storage guard')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, sha256=sha256(data).hexdigest(), bytes=len(data), scores=scores,
                elapsed_seconds=result['elapsed_seconds'], regional_logs=result['regional_logs'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
