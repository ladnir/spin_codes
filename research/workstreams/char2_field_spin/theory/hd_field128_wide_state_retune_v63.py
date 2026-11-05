"""Retune V60 state40 proof witnesses; no change to its construction.

Fast V62 logs propose marks only. Every accepted witness is rebuilt from
exact V60 rational caps and checked across all 4,198,401 midpoint profiles.
These are conditional FLOAT upper-bound screens, not certificates.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
import resource
import shutil
from time import monotonic

import numpy as np
from scipy.optimize import minimize
import hd_field128_wide_state_profile_v61 as previous
import hd_field128_wide_state_logs_v62 as fast

HERE = Path(__file__).resolve().parent
INPUT = 'hd_field128_wide_state_profile_v61_20261010.json'
INPUT_SHA = '1d157fb50ac8fcd31afe8f3a503d7546684b1b727ecb7769be103b6e2cbb13c5'
CHECK = 'hd_field128_wide_state_logs_check_v62_20261010.json'
CHECK_SHA = '8de6e830c6e7d5f25c26ea5787d7d63c8ab2aa447c7ef1e03531933701d254c5'
STATE, LIMIT, RESERVE = 40, 192 << 10, 20 << 30


def source_pins():
    pins = previous.source_pins()
    for name, digest in fast.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting source identity')
        pins[name] = digest
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def decode_point(point, z0, alpha0, left0, right0):
    if len(point) != 6 or not all(isfinite(v) and abs(v) <= 4 for v in point):
        raise ValueError('Bounded six-dimensional proposal required')
    z = 2.**(log2(z0)*2.**float(point[0]))
    alpha = alpha0*2.**float(point[1])
    left = (left0[0]+point[2], left0[1]+point[3])
    right = (right0[0]+point[4], right0[1]+point[5])
    if not (0 < z < 1 and 0 < alpha <= 1
            and all(isfinite(v) and abs(v) <= 32 for v in left+right)):
        raise ValueError('Invalid mark or positive-marker logarithm')
    return float(z), float(alpha), tuple(map(float, left)), tuple(map(float, right))


def run(name, *, rounds=2):
    if HERE != Path('/tmp/spin-field128-proof-52GBmlZD/bundle/theory'):
        raise RuntimeError('Numerical search is remote-only')
    if Path(name).name != name or not name.endswith('.json') or not 1 <= rounds <= 2:
        raise ValueError('Expected result basename and at most two rounds')
    output = HERE/name
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Exclusive output and planned reserve required')
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (550, 555))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    framework, provider = previous.framework, previous.provider
    exchange = framework.control.previous.previous.exchange
    old, alias, gamma = framework.caps.prior.control.old, framework.alias, framework.gamma
    counts_module = exchange.previous.previous.counts_module
    paths = [HERE/previous.SOURCE, HERE/INPUT, HERE/CHECK]
    if any(path.stat().st_size > LIMIT for path in paths):
        raise ValueError('Oversized pinned input')
    raws = [path.read_bytes() for path in paths]
    if (sha256(raws[0]).hexdigest() != previous.SOURCE_SHA
            or sha256(raws[1]).hexdigest() != INPUT_SHA
            or sha256(raws[2]).hexdigest() != CHECK_SHA):
        raise ValueError('Frozen witness/baseline input changed')
    check = json.loads(raws[2])
    fast_pins = fast.source_pins()
    if (check['status'] != 'PASS' or check['scalar_entries_checked'] != 5460
            or any(check['source_sha256'].get(k) != v for k, v in fast_pins.items())):
        raise RuntimeError('Matching complete V62 exact-log check required')
    specs = json.loads(raws[0])['witness_specs']
    baseline = next(row for row in json.loads(raws[1])['rows'] if row['state'] == STATE)
    pins, started, budget = source_pins(), monotonic(), gamma.old.Budget(530.)
    regions, prefixes, local_hashes = [], [], []

    def build(spec):
        budget.check()
        z_exact, alpha_exact = F(spec['z']), F(spec['alpha'])
        if F(float(z_exact)) != z_exact or F(float(alpha_exact)) != alpha_exact:
            raise ValueError('Exact binary64 witness literals required')
        z, alpha = float(z_exact), float(alpha_exact)
        logs = provider.local_logs(z_exact, STATE)
        # Source exact rationals, NOT the faster proposal approximation.
        coefficient = old.Coefficients(logs, z, alpha, budget)
        odds = [tuple(log2(float(F(v))) for v in marker) for marker in spec['markers']]
        bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128,
            distance=F(9, 20), grid=128, floor_policy='positive-fallback',
            max_table_bytes=384 << 20, _coefficients=coefficient, _budget=budget)
        region = bank.regional(4096)
        if region.shape != (2049, 2, 2) or bank.odds != tuple(odds):
            raise ArithmeticError('Incomplete regional domain or changed markers')
        return region, bank.constant+bank.alpha*bank.density*4096, sha256(logs.tobytes()).hexdigest()

    for spec in specs:
        region, prefix, digest = build(spec)
        regions.append(region)
        prefixes.append(prefix)
        local_hashes.append(digest)
    if (local_hashes != baseline['local_log_sha256']
            or [sha256(v.tobytes()).hexdigest() for v in regions] != baseline['regional_sha256']):
        raise ArithmeticError('Exact V61 state40 replay mismatch')
    counts = counts_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
    completed = exchange.envelope(regions, prefixes, counts, budget)
    if abs(completed['score']-baseline['final']['score']) > 1e-7:
        raise ArithmeticError('Baseline full-envelope replay mismatch')
    initial, history, rejected, stop = completed, [], [], 'round_limit'
    density = gamma.old.outer_log_density(1 << 128, 'invertible-toeplitz')
    fast.local_logs(.5, STATE)  # Prepare exact coefficient tables once.
    print(json.dumps(dict(baseline_score=completed['score'],
                          elapsed_seconds=monotonic()-started)), flush=True)
    for round_id in range(rounds):
        try:
            budget.check()
            worst = completed['largest_profile']
            left, right, index = worst['even_full'], worst['odd_full'], worst['witness']
            seed_spec = specs[index]
            z0, a0 = float(F(seed_spec['z'])), float(F(seed_spec['alpha']))
            count = float(counts.log_counts(left, right, budget=budget))
            seed_coefficient = old.Coefficients(fast.local_logs(z0, STATE), z0, a0, budget)
            starts = []
            for marker in seed_spec['markers']:
                odds = tuple(log2(float(F(v))) for v in marker)
                atlas = alias.alias_table(z0, a0, odds, _coefficients=seed_coefficient, _budget=budget)
                regional = exchange.marker_regional(atlas, odds)
                starts.append((regional[left], regional[right], odds))
            transfer = counts_module.near.provider.profiles.transfer
            lefts, rights = np.array([row[0] for row in starts]), np.array([row[1] for row in starts])
            products = transfer.power(transfer.multiply(lefts[:, None], rights[None]), 128)
            scores = np.logaddexp2(products[..., 0, 0], products[..., 0, 1])
            pair = np.unravel_index(np.argmin(scores), scores.shape)
            left0, right0 = starts[pair[0]][2], starts[pair[1]][2]
            coefficients, trials = {(z0, a0): seed_coefficient}, []

            def evaluate(point):
                budget.check()
                try:
                    z, alpha, lo, ro = decode_point(point, z0, a0, left0, right0)
                except ValueError:
                    return 1e12
                key = z, alpha
                if key not in coefficients:
                    if len(coefficients) >= 8:
                        coefficients.pop(next(iter(coefficients)))
                    coefficients[key] = old.Coefficients(fast.local_logs(z, STATE), z, alpha, budget)
                lc = alias.alias_table(z, alpha, lo, _coefficients=coefficients[key], _budget=budget)
                rc = alias.alias_table(z, alpha, ro, _coefficients=coefficients[key], _budget=budget)
                l = exchange.marker_regional(lc, lo)[left]
                r = exchange.marker_regional(rc, ro)[right]
                prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096, density)
                score = float(prefix+count+old.prior.profile_transfer(l, r))
                trials.append(dict(point=list(map(float, point)), z=z, alpha=alpha,
                    left_odds=list(lo), right_odds=list(ro), score=score))
                return score

            # Wider coarse seed offsets test a shifted optimum before small
            # local steps. Distinct region histograms retain separate markers.
            for a in (-.5, 0., .5):
                for b in (-.5, 0., .5):
                    evaluate(np.array([a, b, 0., 0., 0., 0.]))
            seeds = sorted(trials, key=lambda row: row['score'])[:2]
            coarse_scores = [row['score'] for row in trials]
            for seed in seeds:
                point = np.asarray(seed['point'])
                simplex = np.tile(point, (7, 1))
                for k, step in enumerate((.08, .08, .15, .15, .15, .15)):
                    simplex[k+1, k] += step
                minimize(evaluate, point, method='Nelder-Mead', options=dict(
                    initial_simplex=simplex, maxfev=80, xatol=1e-5, fatol=1e-3))
            selected = min(trials, key=lambda row: row['score'])
            markers = list(seed_spec['markers'])
            for key in ('left_odds', 'right_odds'):
                marker, unused = exchange.literal_marker(selected[key])
                if marker not in markers:
                    markers.append(marker)
            candidate = dict(z=str(F(selected['z'])), alpha=str(F(selected['alpha'])), markers=markers)
            region, prefix, digest = build(candidate)
            local = float(prefix+count+old.prior.profile_transfer(region[left], region[right]))
            if local >= worst['log2_term']-1e-6:
                rejected.append(dict(round=round_id, reason='no_exact_local_improvement',
                    target=[left, right], exact_local_score=local, selected=selected,
                    trials=len(trials), coarse_scores=coarse_scores))
                stop = 'no_local_improvement'
                break
            replay = exchange.envelope(regions+[region], prefixes+[prefix], counts, budget)
            if replay['profiles_evaluated'] != 2049**2 or replay['score'] > completed['score']+1e-7:
                raise ArithmeticError('Incomplete or nonmonotone full replay')
            # Commit only a complete replay. No partial envelope or bound from
            # another construction may enter the retained result.
            regions, prefixes = regions+[region], prefixes+[prefix]
            specs, local_hashes = specs+[candidate], local_hashes+[digest]
            history.append(dict(round=round_id, target=[left, right], parent_witness=index,
                new_witness=len(specs)-1, trials=len(trials), coarse_scores=coarse_scores,
                selected=selected, exact_local_score=local, before=completed, after=replay))
            improvement = completed['score']-replay['score']
            completed = replay
            print(json.dumps(dict(round=round_id, score=completed['score'], gain=improvement,
                worst=completed['largest_profile'], elapsed_seconds=monotonic()-started)), flush=True)
            if completed['score'] < -80 or improvement < .01:
                stop = 'selected_weight_FLOAT_target' if completed['score'] < -80 else 'small_union_gain'
                break
        except TimeoutError:
            stop = 'budget_last_complete_round_preserved'
            break
    if source_pins() != pins or any(path.read_bytes() != raw for path, raw in zip(paths, raws)):
        raise RuntimeError('Source or input changed')
    result = dict(schema='wide-state40-retuned-midpoint-v63-FLOAT', source_sha256=pins,
        inputs=[dict(path=p.name, sha256=sha256(raw).hexdigest()) for p, raw in zip(paths, raws)],
        construction=provider.construction(STATE), state=STATE, target_distance='9/20',
        active_groups=4096, initial=initial, final=completed, completed_rounds=history,
        rejected=rejected, stop_reason=stop, witness_specs=specs, local_log_sha256=local_hashes,
        regional_sha256=[sha256(v.tobytes()).hexdigest() for v in regions],
        exact_rational_local_caps_for_all_retained_witnesses=True,
        fast_logs_used_only_for_proposals=True, all_profiles_retained=True,
        original_r_for_coefficients_and_prefix=True,
        witness_min_after_own_complete_transfer=True, cross_construction_minimum=False,
        conditional_on_principal_setup=True, setup_risk_added=False,
        directed=False, whole_code_certificate=False, encoder_implemented=False,
        elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Output storage cap')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=name, bytes=len(data), sha256=sha256(data).hexdigest(),
        score=completed['score'], rounds=len(history), stop=stop,
        elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--rounds', type=int, default=2)
    args = parser.parse_args()
    print(json.dumps(run(args.output, rounds=args.rounds), allow_nan=False), flush=True)
