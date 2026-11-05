"""Add two-marker FLOAT witnesses, retaining the V51 envelope at r4096.

Tune z, alpha and distinct even/odd markers for the current worst profile. Rebuild every
marker at the new parameters, then replay every profile before committing.
This is conditional one-weight numerical evidence, not a certificate.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_reversed_witness_exchange_v51 as previous

exchange = previous.previous

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = '1c637469f50e87f09b3f253c9ae89e51c774224467fff522101e690eaaab6d0e'
INPUT = 'hd_field128_reversed_witness_exchange_v51_20261010.json'
INPUT_SHA = 'e4562fa1feb10010380612c0fcc7ffdba8a0a667ba9718d2455d6c9317b54a2e'
LIMIT, RESERVE = 96 << 10, 20 << 30


def source_pins():
    pins = previous.source_pins()
    if pins[Path(previous.__file__).name] != PREVIOUS_SHA:
        raise RuntimeError('V51 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(output_name, *, rounds=3):
    with exchange.remote_limits():
        return _run(output_name, rounds)


def _run(output_name, rounds):
    import numpy as np
    from scipy.optimize import minimize
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    import hd_field128_reversed_return_v44 as caps
    old, counts_module = caps.prior.control.old, exchange.previous.previous.counts_module
    output = old.prior.previous.inside(output_name)
    if not 1 <= rounds <= 3 or output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Round/output/storage guard')
    raw, pins = (HERE/INPUT).read_bytes(), source_pins()
    if len(raw) > LIMIT or sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Pinned V51 required')
    artifact = json.loads(raw)
    specs = artifact['witness_specs']
    started, budget = monotonic(), gamma.old.Budget(165.)
    regions, prefixes = [], []
    def build(spec):
        z, alpha = float(F(spec['z'])), float(F(spec['alpha']))
        coefficient = old.Coefficients(caps.local_logs(z), z, alpha, budget)
        odds = [tuple(log2(float(F(v))) for v in pair) for pair in spec['markers']]
        bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128,
            distance=F(9, 20), grid=128, floor_policy='positive-fallback',
            max_table_bytes=384 << 20, _coefficients=coefficient, _budget=budget)
        result = bank.regional(4096)
        if result.shape != (2049, 2, 2):
            raise ArithmeticError('Incomplete regional domain')
        prefix = bank.constant+bank.alpha*bank.density*4096
        return result, prefix
    for spec in specs:
        region, prefix = build(spec)
        regions.append(region)
        prefixes.append(prefix)
    if [sha256(x.tobytes()).hexdigest() for x in regions] != artifact['regional_sha256']:
        raise ArithmeticError('V51 regional replay mismatch')
    counts = counts_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
    completed = exchange.envelope(regions, prefixes, counts, budget)
    if abs(completed['score']-artifact['final']['score']) > 1.e-6:
        raise ArithmeticError('V51 envelope replay mismatch')
    initial, history, rejected, stop = completed, [], [], 'round_limit'
    density = gamma.old.outer_log_density(1 << 128, 'invertible-toeplitz')
    for round_id in range(rounds):
        try:
            budget.check()
            worst = completed['largest_profile']
            index = worst['witness']
            left, right = worst['even_full'], worst['odd_full']
            seed_spec = specs[index]
            z0, a0 = float(F(seed_spec['z'])), float(F(seed_spec['alpha']))
            count = float(counts.log_counts(left, right, budget=budget))
            seed_coefficient = old.Coefficients(caps.local_logs(z0), z0, a0, budget)
            starts = []
            for marker in seed_spec['markers']:
                odds = tuple(log2(float(F(v))) for v in marker)
                atlas = alias.alias_table(z0, a0, odds, _coefficients=seed_coefficient, _budget=budget)
                regional = exchange.marker_regional(atlas, odds)
                starts.append((regional[left], regional[right], odds))
            transfer = counts_module.near.provider.profiles.transfer
            lefts, rights = np.array([v[0] for v in starts]), np.array([v[1] for v in starts])
            products = transfer.power(transfer.multiply(lefts[:, None], rights[None]), 128)
            scores = np.logaddexp2(products[..., 0, 0], products[..., 0, 1])
            seed_pair = np.unravel_index(np.argmin(scores), scores.shape)
            odds_left, odds_right = starts[seed_pair[0]][2], starts[seed_pair[1]][2]
            coefficients = {(z0, a0): seed_coefficient}
            trials = []
            def evaluate(point):
                budget.check()
                if not all(isfinite(v) and abs(v) <= 4 for v in point):
                    return 1.e12
                z = 2.**(log2(z0)*2.**point[0])
                alpha = a0*2.**point[1]
                left_odds = (odds_left[0]+point[2], odds_left[1]+point[3])
                right_odds = (odds_right[0]+point[4], odds_right[1]+point[5])
                if not (0 < z < 1 and 0 < alpha <= 1 and all(abs(v) <= 32 for v in left_odds+right_odds)):
                    return 1.e12
                key = float(z), float(alpha)
                if key not in coefficients:
                    if len(coefficients) >= 8:
                        coefficients.pop(next(iter(coefficients)))
                    coefficients[key] = old.Coefficients(caps.local_logs(z), z, alpha, budget)
                atlas_left = alias.alias_table(z, alpha, left_odds, _coefficients=coefficients[key], _budget=budget)
                atlas_right = alias.alias_table(z, alpha, right_odds, _coefficients=coefficients[key], _budget=budget)
                regional_left = exchange.marker_regional(atlas_left, left_odds)[left]
                regional_right = exchange.marker_regional(atlas_right, right_odds)[right]
                prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096, density)
                score = prefix+count+old.prior.profile_transfer(regional_left, regional_right)
                trials.append(dict(z=float(z), alpha=float(alpha),
                    left_odds=[float(v) for v in left_odds], right_odds=[float(v) for v in right_odds],
                    score=float(score)))
                return float(score)
            evaluate(np.zeros(6))
            simplex = np.zeros((7, 6))
            for k, step in enumerate((.04, .04, .1, .1, .1, .1)):
                simplex[k+1, k] = step
            minimize(evaluate, np.zeros(6), method='Nelder-Mead', options=dict(
                initial_simplex=simplex, maxfev=160, xatol=1.e-5, fatol=1.e-3))
            selected = min(trials, key=lambda p: p['score'])
            additions = [exchange.literal_marker(selected[key])[0] for key in ('left_odds', 'right_odds')]
            markers = list(seed_spec['markers'])
            for marker in additions:
                if marker not in markers:
                    markers.append(marker)
            candidate_spec = dict(z=str(F(selected['z'])), alpha=str(F(selected['alpha'])), markers=markers)
            # Rebuild all caps from the exact serialized literals. No old
            # regional entries or prefix are reused for this new witness.
            candidate_region, candidate_prefix = build(candidate_spec)
            exact_score = candidate_prefix+count+old.prior.profile_transfer(candidate_region[left], candidate_region[right])
            if exact_score >= worst['log2_term']-1.e-6:
                rejected.append(dict(round=round_id, target=[left, right], reason='no_local_improvement',
                                     exact_score=exact_score, trials=len(trials)))
                stop = 'no_local_improvement'
                break
            proposed_regions, proposed_prefixes = regions+[candidate_region], prefixes+[candidate_prefix]
            replay = exchange.envelope(proposed_regions, proposed_prefixes, counts, budget)
            if replay['score'] > completed['score']+1.e-7 or replay['profiles_evaluated'] != 2049**2:
                raise ArithmeticError('Nonmonotone or incomplete envelope')
            # No mutation of old arrays/specs, even on an incomplete round.
            regions, prefixes = proposed_regions, proposed_prefixes
            specs = specs+[candidate_spec]
            history.append(dict(round=round_id, parent_witness=index, new_witness=len(specs)-1,
                target=[left, right], trials=len(trials), exact_local_score=exact_score,
                selected_optimizer_score=selected['score'], before=completed, after=replay))
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
    result = dict(schema='reversed-witness-exchange-v52-FLOAT', source_sha256=pins,
        input=dict(path=INPUT, sha256=INPUT_SHA), initial=initial, final=completed,
        completed_rounds=history, rejected=rejected, stop_reason=stop, witness_specs=specs,
        active_groups=4096, target_distance='9/20', conditional_on_principal_setup=True,
        original_r_for_coefficients_and_prefix=True, min_after_complete_transfer_and_prefix=True,
        separate_positive_markers_for_even_odd_profiles=True,
        old_complete_witnesses_retained=True, all_profiles_retained=True,
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
    parser.add_argument('--rounds', type=int, default=3)
    args = parser.parse_args()
    print(json.dumps(run(args.output, rounds=args.rounds), allow_nan=False), flush=True)


