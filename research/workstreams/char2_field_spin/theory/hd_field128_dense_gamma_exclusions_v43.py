"""Conditional FLOAT screen with explicit once-global gamma exclusions.

For selected (D,T,Z,U) partitions, exclude roots of one canonically chosen
nonzero minor numerator. Its degree <=31 implies at most16 degree16 roots.
This changes the setup hypothesis, not the sampling distribution. Risk is
reported separately; no selected-profile score is a whole-code certificate.
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
LIMIT, Q = 48 << 10, 1 << 128


def targeted(j, c, e, cutoff):
    if not 58 <= cutoff <= 62:
        raise ValueError('Cutoff must be58..62')
    t, u = j-c, min(32, j-c)
    if not 0 <= c < j <= 64 or not 0 <= e <= 64-j or e >= 32 or c+e < cutoff:
        return False
    h = max(0, e+u-32, min((e+1)//2, t), min(e, (t+1)//2))
    return h == t and prior.rank_floor(j, c, e) < t


def exclusion_count(cutoff):
    shapes, partitions = [], 0
    for j in range(1, 65):
        for c in range(j):
            for e in range(65-j):
                if targeted(j, c, e, cutoff):
                    count = comb(64, j)*comb(j, c)*comb(64-j, e)
                    shapes.append([j, c, e, count])
                    partitions += count
    numerator = 16*partitions+8*comb(64, 32)
    return dict(cutoff=cutoff, shapes=shapes, ordered_partitions=partitions,
        bad_gamma_upper=str(F(numerator, Q-(1 << 64))),
        bad_gamma_margin_bits_FLOAT=128+log2(1-2.**-64)-log2(numerator),
        principal_risk_included_once=True, one_minor_per_partition=True)


def local_logs(z, cutoff):
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    logs = prior.local_logs(z)
    lz, lw, log_d = log2(z), log2(1-z), 4096+gamma.old.log1mexp2(-4096)
    for j in range(1, 65):
        terms = []
        for c in range(j):
            u = min(32, j-c)
            for e in range(min(32, 65-j)):
                if targeted(j, c, e, cutoff):
                    continue
                delta = u-prior.rank_floor(j, c, e)
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
    # Inherit90s/1GiB/32KiB hard limits; actual output must stay below32KiB.
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
        raise ValueError('Pinned V34 input required')
    row = old.prior.validate_source(json.loads(raw))
    w, profile = row['witness'], row['largest_profile']
    z, alpha = float(F(w['z'])), float(F(w['alpha']))
    markers, _, _ = control.marker_extension(w['markers'], alpha)
    prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096,
        gamma.old.outer_log_density(Q, 'invertible-toeplitz'))
    if abs(prefix-row['log2_outer_prefix']) > 1e-7 or (profile['even_full'], profile['odd_full']) != (2026, 2026):
        raise ArithmeticError('Original-r/profile guard')
    offset = prefix+profile['log2_count_upper']
    started, budget, rows = monotonic(), gamma.old.Budget(80.), []
    for cutoff in (62, 61, 60, 59, 58):
        budget.check()
        logs = local_logs(z, cutoff)
        coeff, best = old.Coefficients(logs, z, alpha, budget), control.Best()
        for i, literal in enumerate(markers):
            budget.check()
            odds = tuple(log2(float(F(x))) for x in literal)
            atlas = alias.alias_table(z, alpha, odds, _coefficients=coeff, _budget=budget)
            value, reliable = old.selected(atlas, odds)
            best.add(i, value, reliable)
        counts = exclusion_count(cutoff)
        # Preserve the exact sum, not a bulky list of every shape.
        counts['targeted_shapes'] = len(counts.pop('shapes'))
        rows.append(dict(**counts, conditional_selected_profile_score=offset+old.prior.profile_transfer(best.logs, best.logs),
            regional_logs=best.logs.tolist(), reliable=best.reliable.tolist(), winners=best.winners.tolist(),
            return_at_occupancy32=float(logs[32, 1, 0])))
    if source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Source changed')
    result = dict(schema='dense-feedback-targeted-gamma-exclusions-v43-FLOAT', source_sha256=pins,
        input_sha256=sha256(raw).hexdigest(), source=source_name, witness=w, markers=markers,
        active_groups=4096, selected_profile=[2026, 2026], target_distance='9/20', rows=rows,
        setup='one independent global uniform gamma from GF2^128 minus GF2^64',
        construction='A128 emission,A64-transpose feedback,state32; outer,route,refresh unchanged',
        output_zero_terms_removed_only_under_explicit_good_event=True,
        gamma_risk_is_once_global=True, gamma_resampling=False, conditional_scores_exclude_setup_risk=True,
        rank_union_uses_one_canonical_minor_not_all_minors=True,
        directed=False, all_profiles_replayed=False, whole_code_certificate=False,
        elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, separators=(',', ':'), sort_keys=True, allow_nan=False)+'\n').encode()
    if len(data) > 32 << 10 or shutil.disk_usage(HERE).free-len(data) < prior.RESERVE:
        raise RuntimeError('32KiB output/reserve guard')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, sha256=sha256(data).hexdigest(), bytes=len(data), rows=rows, elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
