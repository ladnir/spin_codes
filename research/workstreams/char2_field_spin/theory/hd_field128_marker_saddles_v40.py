"""Same-witness FLOAT marker diagnosis, not a certificate.
For J0 only: region00<=2^(128alpha*c-alpha*r*128), c=max(logK_j+128j),
since every route has sum(j)=r. Region11<=max(P_j11)^(128alpha).
Also P_j11<=B*t^64, t=z+(1-z)/q, by MDS zero-set probabilities<=q^-e.
Full-return diagonals include leaving/reentering paths: these caps do NOT apply.
"""
import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
import shutil
import sys
from time import monotonic

import hd_field128_return_grid_paths_v38 as control

HERE = Path(__file__).resolve().parent
CONTROL_SHA = '1e0fc4ee6f4858bf3c6877c987831dc3bbb8d38c2aa18e1ed2117e6ba8a3252d'
LIMIT, RESERVE, Q, R = 32 << 10, 20 << 30, 1 << 128, 4096
old, Best = control.old, control.Best


def source_pins():
    pins = control.source_pins()
    if pins.get(Path(control.__file__).name) != CONTROL_SHA:
        raise RuntimeError('Frozen v38 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def marker_extension(original, alpha):
    if not isfinite(alpha) or not 0 < alpha <= 1:
        raise ValueError('Invalid alpha')
    center = (F(44, 2026), F(1))
    # Exact rational literals of the proposed binary64 quiet saddle.
    quiet = (F(float(center[0])*2.**(128*alpha)), F(2.**(256*alpha)))
    result, labels = [list(pair) for pair in original], []
    seen = {tuple(F(v) for v in pair) for pair in original}
    for name, point in (('histogram', center), ('quiet-weighted', quiet)):
        for x in (F(1, 2), F(1), F(2)):
            for y in (F(1, 2), F(1), F(2)):
                pair = point[0]*x, point[1]*y
                if pair not in seen:
                    seen.add(pair)
                    labels.append(dict(index=len(result), family=name, factors=[str(x), str(y)]))
                    result.append(list(map(str, pair)))
    return result, labels, dict(histogram=list(map(str, center)), quiet_weighted=list(map(str, quiet)))


def diagonal_caps(local, z, alpha):
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    local = np.asarray(local)
    if (local.shape != (65, 2, 2) or np.any(np.isnan(local)) or np.any(np.isposinf(local))
            or not 0 < z < 1 or not 0 < alpha <= 1):
        raise ValueError('Invalid local table/marks')
    c = float(np.max(local[:, 0, 0]+128*np.arange(65)))
    d = float(np.max(local[:, 1, 1]))
    t = float(np.logaddexp2(log2(z), log2(1-z)-128))
    log_b = gamma.old.log1mexp2(-4096)-32*gamma.old.log1mexp2(-128)
    symbolic = log_b+64*t
    return dict(quiet_residual_max=c, live_packet_max=d, live_symbolic_uniform=symbolic,
        quiet_regional=128*alpha*c-alpha*R*128,
        live_regional_from_local=128*alpha*d,
        live_regional_symbolic=128*alpha*symbolic,
        live_regional=128*alpha*min(d, symbolic))


def cap_j0(matrix, caps):
    import numpy as np
    value = np.array(matrix, dtype=float, copy=True)
    if value.shape != (2, 2) or not np.isneginf(value[1, 0]):
        raise ValueError('Independent diagonal caps apply only to exact triangular J0')
    value[0, 0] = min(value[0, 0], caps['quiet_regional'])
    value[1, 1] = min(value[1, 1], caps['live_regional'])
    return value


def snapshot(full, zero, caps, offset):
    import numpy as np
    def matrix(value):
        return [[None if np.isneginf(x) else float(x) for x in row] for row in value]
    capped = cap_j0(zero.logs, caps)
    return dict(full_return_profile=offset+old.prior.profile_transfer(full.logs, full.logs),
        optimistic_J0_profile=offset+old.prior.profile_transfer(zero.logs, zero.logs),
        optimistic_J0_capped_profile=offset+old.prior.profile_transfer(capped, capped),
        full_return_regional=matrix(full.logs), J0_regional=matrix(zero.logs), J0_capped_regional=matrix(capped),
        full_return_winners=full.winners.tolist(), J0_winners=zero.winners.tolist(),
        full_return_reliable=full.reliable.tolist(), J0_reliable=zero.reliable.tolist(),
        capped_all_quiet=offset+256*capped[0, 0],
        capped_one_birth=offset+old.prior.power(capped, 256)[0, 1])


@contextmanager
def remote_limits():
    if sys.platform != 'linux':
        raise RuntimeError('Production requires Linux')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (85, 90))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    signal.alarm(90)
    try:
        yield
    finally:
        signal.alarm(0)


def run(source_name, output_name):
    with remote_limits():
        return _run(source_name, output_name)


def _run(source_name, output_name):
    source, output = old.prior.previous.inside(source_name), old.prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20:
        raise ValueError('Exclusive output / bounded input required')
    if shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise RuntimeError('20GiB planned reserve required')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != old.prior.SOURCE_SHA:
        raise ValueError('Pinned v34 result required')
    row = old.prior.validate_source(json.loads(raw))
    witness, profile = row['witness'], row['largest_profile']
    if (len(witness['markers']) != 79 or witness['grid'] != 128
            or (profile['even_full'], profile['odd_full']) != (2026, 2026)):
        raise ValueError('Expected79markers/grid128/profile2026,2026')
    def check():
        if source.read_bytes() != raw or source_pins() != pins:
            raise RuntimeError('Input or source closure changed')
    try:
        import hd_field128_full_rs_gamma_screen as gamma
        import hd_field128_lift_alias_screen as alias
        started, budget = monotonic(), gamma.old.Budget(80.)
        z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
        local = gamma.local_logs(z, Q, state=32)
        full_coeff = old.Coefficients(local, z, alpha, budget)
        zero_coeff = old.Coefficients(old.ablate(local), z, alpha, budget)
        caps = diagonal_caps(local, z, alpha)
        markers, labels, centers = marker_extension(witness['markers'], alpha)
        prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), R,
            gamma.old.outer_log_density(Q, 'invertible-toeplitz'))
        if abs(prefix-row['log2_outer_prefix']) > 1e-7:
            raise ArithmeticError('Original-r prefix mismatch')
        offset = prefix+profile['log2_count_upper']
        full, zero, before = Best(), Best(), None
        for index, literal in enumerate(markers):
            budget.check()
            odds = tuple(log2(float(F(v))) for v in literal)
            atlas = alias.alias_table(z, alpha, odds, _coefficients=full_coeff, _budget=budget)
            value, reliable = old.selected(atlas, odds)
            full.add(index, value, reliable)
            tri = old.triangular_atlas(zero_coeff.composed, odds, full_coeff.pair_logs(odds, _budget=budget),
                                       check=budget.check)
            value, reliable = old.selected(tri, odds)
            zero.add(index, value, reliable)
            if index == 78:
                before = snapshot(full, zero, caps, offset)
                if (abs(before['full_return_profile']-old.EXPECTED_BASELINE) > 1e-4
                        or abs(before['optimistic_J0_profile']-3561.912788278656) > 1e-4):
                    raise ArithmeticError('V37 controls differ')
        budget.check()
        result = dict(schema='gamma32-marker-saddles-selected-profile-FLOAT-v40',
            source=dict(path=source_name, sha256=old.prior.SOURCE_SHA), source_sha256=pins,
            conditional_construction=old.prior.previous.state_metadata(32), active_groups=R,
            selected_profile=[2026, 2026], target_distance='9/20', witness=witness,
            marker_literals=markers, new_marker_labels=labels, saddle_centers=centers,
            original79_retained_unchanged=markers[:79] == witness['markers'],
            completed_markers=len(markers), grid=128, same_z_alpha=True,
            diagonal_caps_J0_only=caps, before=before, after=snapshot(full, zero, caps, offset),
            log2_outer_prefix=prefix, log2_profile_count=profile['log2_count_upper'],
            count_source='pinned FLOAT profile count cap; not recomputed',
            normalization='original H4096/count/prefix, each once', all_profiles_replayed=False,
            grid_stability_implies_coefficient_equality=False, symbolic_live_cap_uses_t_not_z=True,
            full_return_formula_valid_conditional=True, J0_is_valid_upper_bound=False,
            setup_risk_added_to_score=False, optimizer_run=False, directed=False,
            selected_weight_certificate=False, whole_code_certificate=False,
            arithmetic='FLOAT with positive fallback; no outward rounding',
            scope='one-profile alias-marker diagnostic; no impossibility or distance claim',
            elapsed_seconds=monotonic()-started,
            limits=dict(wall_seconds=90, address_space_bytes=1 << 30, output_bytes=LIMIT,
                        minimum_free_after_write=RESERVE))
        check()
        data = (json.dumps(old.prior.previous.json_ready(result), sort_keys=True,
                           separators=(',', ':'), allow_nan=False)+'\n').encode()
        if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
            raise RuntimeError('32KiB output /20GiB planned reserve exceeded')
        with output.open('xb') as stream:
            stream.write(data)
        if output.read_bytes() != data:
            raise RuntimeError('Exclusive output changed')
        return dict(path=output_name, sha256=sha256(data).hexdigest(), before=before,
                    after=result['after'], markers=len(markers))
    finally:
        check()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    answer = run(args.source, args.output)
    print(json.dumps({k: answer[k] for k in ('path', 'sha256', 'markers')}, allow_nan=False), flush=True)
