"""Monotone complete-profile marker exchange, conditional FLOAT at r4096.

Each witness keeps its own z, alpha, and all previous positive markers.
Entrywise minima occur only within one witness. Cross-witness minima follow
the complete ordered transfer and its own outer prefix. No certificate.
"""
import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_reversed_profile_min_v48 as previous

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = 'a509e13a43bde40490bf07e961d7efb775998218a38f1aeeba2d6f8f8dc83714'
INPUT = 'hd_field128_reversed_profile_min_v48_20261010.json'
INPUT_SHA = '85dee604671a28264156c8619358e54d906182ce07cae3b0b237188f9cf204de'
LIMIT, RESERVE = 96 << 10, 20 << 30


def source_pins():
    pins = previous.source_pins()
    if pins[Path(previous.__file__).name] != PREVIOUS_SHA:
        raise RuntimeError('V48 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def literal_marker(odds):
    if len(odds) != 2 or not all(isfinite(x) and abs(x) <= 32 for x in odds):
        raise ValueError('Bounded finite marker required')
    literals = [str(F(2.**float(x))) for x in odds]
    return literals, tuple(log2(float(F(x))) for x in literals)


def marker_regional(table, odds):
    """Same normalization and positive fallback as frozen AliasBankScreen."""
    import numpy as np
    import hd_field128_lift_alias_screen as alias
    full = np.arange(2049)
    single = 4096-2*full
    norm = (alias.LOG_FACTORIAL[4096]-alias.LOG_FACTORIAL[full]
            -alias.LOG_FACTORIAL[single]-alias.LOG_FACTORIAL[full])
    grid = table.logs.shape[1]
    good = table.reliable[:, single % grid, full % grid]
    numerator = np.where(good, table.logs[:, single % grid, full % grid],
                         table.positive_numerator[:, None])
    logs = numerator-norm-single*odds[0]-full*odds[1]
    return logs.T.reshape((-1, 2, 2))


def envelope(regions, prefixes, counts, budget, *, repetitions=128):
    """Enumerate every ordered profile; never mix different witnesses."""
    import numpy as np
    transfer = previous.previous.counts_module.near.provider.profiles.transfer
    size = regions[0].shape[0]
    if (len(regions) != len(prefixes) or not regions
            or any(x.shape != (size, 2, 2) for x in regions)):
        raise ValueError('Regional domain mismatch')
    full = np.arange(size)
    total, largest, seen = -np.inf, None, 0
    winners = np.zeros(len(regions), dtype=np.int64)
    for start in range(0, size, 32):
        budget.check()
        left = full[start:start+32]
        count_logs = counts.log_counts(left[:, None], full[None, :], budget=budget)
        terms = []
        for table, prefix in zip(regions, prefixes):
            pair = transfer.multiply(table[start:start+32, None], table[None])
            repeated = transfer.power(pair, repetitions)
            terms.append(prefix+count_logs+np.logaddexp2(repeated[..., 0, 0], repeated[..., 0, 1]))
        terms = np.array(terms)
        selected, best = np.argmin(terms, axis=0), np.min(terms, axis=0)
        if np.any(np.isnan(best)) or np.any(np.isposinf(best)):
            raise ArithmeticError('Invalid profile score')
        winners += np.bincount(selected.ravel(), minlength=len(regions))
        total = float(np.logaddexp2(total, np.logaddexp2.reduce(best.ravel())))
        pos = np.unravel_index(np.argmax(best), best.shape)
        if largest is None or best[pos] > largest['log2_term']:
            largest = dict(even_full=int(left[pos[0]]), odd_full=int(full[pos[1]]),
                log2_term=float(best[pos]), witness=int(selected[pos]),
                individual_scores=terms[:, pos[0], pos[1]].tolist())
        seen += best.size
    if seen != size**2 or not isfinite(total):
        raise ArithmeticError('Incomplete profile enumeration')
    return dict(score=total, largest_profile=largest, profiles_evaluated=seen,
                winning_profile_counts=winners.tolist())


@contextmanager
def remote_limits():
    import sys
    if sys.platform != 'linux':
        raise RuntimeError('Numerical production is remote-only')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (175, 180))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    signal.alarm(180)
    try:
        yield
    finally:
        signal.alarm(0)


def run(output_name, *, rounds=6):
    with remote_limits():
        return _run(output_name, rounds)


def _run(output_name, rounds):
    import numpy as np
    from scipy.optimize import minimize
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    import hd_field128_reversed_return_v44 as caps
    old, count_module = caps.prior.control.old, previous.previous.counts_module
    output = old.prior.previous.inside(output_name)
    if not 1 <= rounds <= 6 or output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Round/output/storage guard')
    raw, pins = (HERE/INPUT).read_bytes(), source_pins()
    if len(raw) > 32 << 10 or sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Pinned V48 input required')
    artifact = json.loads(raw)
    specs = artifact['witness_specs']
    started, budget = monotonic(), gamma.old.Budget(165.)
    regions, prefixes, coefficients = [], [], []
    for spec in specs:
        z, alpha = float(F(spec['z'])), float(F(spec['alpha']))
        coefficient = old.Coefficients(caps.local_logs(z), z, alpha, budget)
        odds = [tuple(log2(float(F(v))) for v in pair) for pair in spec['markers']]
        bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128,
            distance=F(9, 20), grid=128, floor_policy='positive-fallback',
            max_table_bytes=384 << 20, _coefficients=coefficient, _budget=budget)
        region = bank.regional(4096)
        if region.shape != (2049, 2, 2):
            raise ArithmeticError('Original-r regional domain')
        regions.append(region)
        prefixes.append(bank.constant+bank.alpha*bank.density*4096)
        coefficients.append(coefficient)
        del bank
    counts = count_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
    completed = envelope(regions, prefixes, counts, budget)
    if abs(completed['score']-artifact['conditional_log2_upper_screen']) > 1.e-4:
        raise ArithmeticError('V48 full-envelope replay mismatch')
    initial, history, rejected, stop = completed, [], [], 'round_limit'
    for round_id in range(rounds):
        try:
            budget.check()
            worst = completed['largest_profile']
            index = worst['witness']
            left, right = worst['even_full'], worst['odd_full']
            spec = specs[index]
            z, alpha = float(F(spec['z'])), float(F(spec['alpha']))
            count = float(counts.log_counts(left, right, budget=budget))
            local_prefix = count+prefixes[index]
            probes = []
            def evaluate(odds):
                budget.check()
                if not all(isfinite(x) and abs(x) <= 32 for x in odds):
                    return 1.e12
                atlas = alias.alias_table(z, alpha, tuple(odds),
                    _coefficients=coefficients[index], _budget=budget)
                candidate = marker_regional(atlas, odds)
                standalone = local_prefix+old.prior.profile_transfer(candidate[left], candidate[right])
                mixed_left = np.minimum(regions[index][left], candidate[left])
                mixed_right = np.minimum(regions[index][right], candidate[right])
                mixed = local_prefix+old.prior.profile_transfer(mixed_left, mixed_right)
                probes.append(dict(odds=[float(x) for x in odds], standalone=float(standalone), mixed=float(mixed)))
                return float(standalone)
            # Seed from old markers and a finite multinomial center. Optimize
            # standalone to avoid flat regions caused by the retained minimum.
            for marker in spec['markers']:
                evaluate(tuple(log2(float(F(v))) for v in marker))
            center = (left+right)/2
            evaluate((log2(max(1., 4096-2*center)/max(1., center)), 0.))
            seed = np.array(min(probes, key=lambda p: p['standalone'])['odds'])
            simplex = np.array([seed, seed+[.1, 0.], seed+[0., .1]])
            minimize(evaluate, seed, method='Nelder-Mead', options=dict(
                initial_simplex=simplex, maxfev=70, xatol=1.e-5, fatol=1.e-4))
            selected = min(probes, key=lambda p: p['mixed'])
            marker, exact_odds = literal_marker(selected['odds'])
            atlas = alias.alias_table(z, alpha, exact_odds,
                _coefficients=coefficients[index], _budget=budget)
            candidate = marker_regional(atlas, exact_odds)
            proposed_region = np.minimum(regions[index], candidate)
            local_after = local_prefix+old.prior.profile_transfer(proposed_region[left], proposed_region[right])
            if marker in spec['markers'] or local_after >= worst['individual_scores'][index]-1.e-6:
                rejected.append(dict(round=round_id, reason='no_new_local_improvement', marker=marker,
                                     exact_literal_score=local_after, probes=len(probes)))
                stop = 'no_local_improvement'
                break
            tentative = list(regions)
            tentative[index] = proposed_region
            replay = envelope(tentative, prefixes, counts, budget)
            if replay['score'] > completed['score']+1.e-7 or replay['profiles_evaluated'] != 2049**2:
                raise ArithmeticError('Nonmonotone or incomplete round')
            # Commit both numerical regions and exported literals only after
            # a complete replay. Any TimeoutError leaves the last state intact.
            regions = tentative
            spec['markers'] = spec['markers']+[marker]
            history.append(dict(round=round_id, witness=index, target=[left, right],
                marker=marker, local_before=worst['individual_scores'][index], local_after=local_after,
                probes=len(probes), before=completed, after=replay))
            completed = replay
            print(json.dumps(dict(round=round_id, score=completed['score'], worst=completed['largest_profile'])), flush=True)
            if completed['score'] <= -80.:
                stop = 'selected_weight_FLOAT_target'
                break
        except TimeoutError:
            stop = 'budget_last_complete_round_preserved'
            break
    if source_pins() != pins or (HERE/INPUT).read_bytes() != raw:
        raise RuntimeError('Source/input changed')
    result = dict(schema='reversed-marker-exchange-v49-FLOAT', source_sha256=pins,
        input=dict(path=INPUT, sha256=INPUT_SHA), initial=initial, final=completed,
        completed_rounds=history, rejected=rejected, stop_reason=stop, witness_specs=specs,
        active_groups=4096, target_distance='9/20', conditional_on_principal_setup=True,
        original_r_for_coefficients_and_prefix=True, min_after_complete_transfer_and_prefix=True,
        old_markers_retained=True, all_profiles_retained=True,
        regional_sha256=[sha256(x.tobytes()).hexdigest() for x in regions],
        directed=False, whole_code_certificate=False, setup_risk_added=False,
        elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Output storage cap')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, bytes=len(data), sha256=sha256(data).hexdigest(),
        score=completed['score'], rounds=len(history), stop=stop, elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--rounds', type=int, default=6)
    args = parser.parse_args()
    print(json.dumps(run(args.output, rounds=args.rounds), allow_nan=False), flush=True)
