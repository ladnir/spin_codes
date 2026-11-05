"""All-profile min of three complete conditional FLOAT witnesses, r4096.

Minimize only complete per-profile bounds after chronological transfer and
the outer prefix. Never mix coefficients or matrix entries across z/alpha.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_reversed_allprofiles_v47 as previous

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = '1d2fb58cbe80ee4b4dc2e0b7296caa414f00b0128d198d4be0d28a8452fd2341'
LIMIT, RESERVE = 32 << 10, 20 << 30
INPUTS = [
    ('hd_field128_d45_state32_search_r4096_v34_20261010/result.json',
     '938de22bbf7e1b25fc343813c404bb50fb73564017b010c4c8efc223250d3aa9'),
    ('hd_field128_reversed_retune_v45_20261010.json',
     'a51b6d2fcf061f59b071a5d3f0f8d2e0f1cb79b4337a0f0aa3f02b63b5aa9f43'),
    ('hd_field128_reversed_marker_retune_v46_20261010.json',
     'f7314a8fcb5c627a1818db943dc94d94f4f03648fb89f13c6f5348bc654f6dd1')]


def source_pins():
    pins = previous.source_pins()
    if pins[Path(previous.__file__).name] != PREVIOUS_SHA:
        raise RuntimeError('V47 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(output_name):
    with previous.previous.previous.prior.prior.control.remote_limits():
        return _run(output_name)


def _run(output_name):
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    import hd_field128_reversed_return_v44 as caps
    control = caps.prior.control
    old, counts_module = control.old, previous.counts_module
    output = old.prior.previous.inside(output_name)
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Storage guard')
    pins, inputs = source_pins(), []
    for name, digest in INPUTS:
        raw = (HERE/name).read_bytes()
        if len(raw) > 64 << 10 or sha256(raw).hexdigest() != digest:
            raise ValueError('Pinned input mismatch')
        inputs.append(json.loads(raw))
    initial = old.prior.validate_source(inputs[0])['witness']
    a0 = float(F(initial['alpha']))
    initial_markers, _, _ = control.marker_extension(initial['markers'], a0)
    specs = [dict(z=initial['z'], alpha=initial['alpha'], markers=initial_markers)]
    specs += [dict(z=a['best']['z'], alpha=a['best']['alpha'], markers=a['markers']) for a in inputs[1:]]
    started, budget, regions, prefixes = monotonic(), gamma.old.Budget(80.), [], []
    for spec in specs:
        z, alpha = float(F(spec['z'])), float(F(spec['alpha']))
        odds = [tuple(log2(float(F(v))) for v in pair) for pair in spec['markers']]
        coefficient = old.Coefficients(caps.local_logs(z), z, alpha, budget)
        bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128, distance=F(9, 20), grid=128,
            floor_policy='positive-fallback', max_table_bytes=384 << 20, _coefficients=coefficient, _budget=budget)
        value = bank.regional(4096)
        if value.shape != (2049, 2, 2):
            raise ArithmeticError('Incomplete regional domain')
        regions.append(value)
        prefixes.append(bank.constant+bank.alpha*bank.density*4096)
        if abs(prefixes[-1]-old.prior.previous.prefactor(z, alpha, F(9, 20), 4096, bank.density)) > 1.e-7:
            raise ArithmeticError('Original-r prefix')
        del bank, coefficient
    counts = counts_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
    transfer = counts_module.near.provider.profiles.transfer
    full = np.arange(2049)
    singles = np.full(3, -np.inf)
    total, largest, seen = -np.inf, None, 0
    winners = np.zeros(3, dtype=np.int64)
    for start in range(0, 2049, 32):
        budget.check()
        left = full[start:start+32]
        count_logs = counts.log_counts(left[:, None], full[None, :], budget=budget)
        terms = []
        for i, table in enumerate(regions):
            pair = transfer.multiply(table[start:start+32, None], table[None])
            repeated = transfer.power(pair, 128)
            value = prefixes[i]+count_logs+np.logaddexp2(repeated[..., 0, 0], repeated[..., 0, 1])
            terms.append(value)
            singles[i] = np.logaddexp2(singles[i], np.logaddexp2.reduce(value.ravel()))
        terms = np.array(terms)
        selected = np.argmin(terms, axis=0)
        best = np.min(terms, axis=0)
        winners += np.bincount(selected.ravel(), minlength=3)
        total = float(np.logaddexp2(total, np.logaddexp2.reduce(best.ravel())))
        pos = np.unravel_index(np.argmax(best), best.shape)
        if largest is None or best[pos] > largest['log2_term']:
            largest = dict(even_full=int(left[pos[0]]), odd_full=int(full[pos[1]]),
                log2_term=float(best[pos]), witness=int(selected[pos]), individual_scores=terms[:, pos[0], pos[1]].tolist())
        seen += best.size
    if seen != 2049**2 or abs(singles[2]-32091.6695407331) > 1.e-4:
        raise ArithmeticError('Coverage or V47 control')
    if source_pins() != pins or any(sha256((HERE/n).read_bytes()).hexdigest() != d for n, d in INPUTS):
        raise RuntimeError('Source/input changed')
    result = dict(schema='reversed-profile-min-v48-FLOAT', inputs=INPUTS, source_sha256=pins,
        conditional_log2_upper_screen=total, individual_allprofile_scores=singles.tolist(),
        largest_profile=largest, winning_profile_counts=winners.tolist(), profiles_evaluated=seen,
        active_groups=4096, target_distance='9/20', all_profiles_retained=True,
        witness_specs=specs, min_after_complete_transfer_and_prefix=True,
        same_translated_feedback_construction=True, directed=False, whole_code_certificate=False,
        setup_risk_added=False, elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Storage cap')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, sha256=sha256(data).hexdigest(), bytes=len(data), score=total,
        individual_scores=singles.tolist(), largest_profile=largest, winning_profiles=winners.tolist(),
        profiles_evaluated=seen, elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), allow_nan=False), flush=True)
