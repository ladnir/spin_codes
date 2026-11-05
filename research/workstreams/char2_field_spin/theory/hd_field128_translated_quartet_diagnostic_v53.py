"""Separate alias slack from four-packet composition at one V52 profile.

Retain all omitted occupancies through V35's positive remainder. The
optimistic boundary ceiling is a diagnostic LOWER bound on the selected
quartet bound, not an upper bound on code failure. No certificate.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_reversed_dual_marker_v52 as previous

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = '7e26a1cc5f0f4da2ea509c7e411f0ecea1613ba921aa2332b0fcee21c88d7c3f'
INPUT = 'hd_field128_reversed_dual_marker_v52_20261010.json'
INPUT_SHA = '1ec14fce6aede5f7f45d301fdeae5cb7057670a09f8e244bad52e9bd1b69b3e5'
LIMIT, RESERVE = 32 << 10, 20 << 30


def source_pins():
    pins = previous.source_pins()
    if pins[Path(previous.__file__).name] != PREVIOUS_SHA:
        raise RuntimeError('V52 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(output_name):
    # The frozen guard enforces Linux-only, 1GiB and a60s hard deadline.
    import hd_field128_return_alias_ablation_v37 as old
    with old.remote_limits():
        return _run(output_name)


def _run(output_name):
    import numpy as np
    import hd_field128_reversed_return_v44 as caps
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    import hd_field128_gamma_four_packet_positive_v35 as four
    exchange = previous.exchange
    old = caps.prior.control.old
    counts_module = exchange.previous.previous.counts_module
    output = old.prior.previous.inside(output_name)
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Output/storage guard')
    raw, pins = (HERE/INPUT).read_bytes(), source_pins()
    if len(raw) > 96 << 10 or sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Pinned V52 required')
    artifact = json.loads(raw)
    worst = artifact['final']['largest_profile']
    left, right, index = worst['even_full'], worst['odd_full'], worst['witness']
    spec = artifact['witness_specs'][index]
    z, alpha = float(F(spec['z'])), float(F(spec['alpha']))
    started, budget = monotonic(), gamma.old.Budget(50.)
    logs = caps.local_logs(z)
    coefficient = old.Coefficients(logs, z, alpha, budget)
    odds = [tuple(log2(float(F(v))) for v in pair) for pair in spec['markers']]
    bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128,
        distance=F(9, 20), grid=128, floor_policy='positive-fallback',
        max_table_bytes=384 << 20, _coefficients=coefficient, _budget=budget)
    if bank.odds != tuple(odds):
        raise ArithmeticError('Marker deduplication changed the table index')
    regions = bank.regional(4096)
    if sha256(regions.tobytes()).hexdigest() != artifact['regional_sha256'][index]:
        raise ArithmeticError('V52 regional digest changed')
    counts = counts_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
    prefix = bank.constant+bank.alpha*bank.density*4096
    count = float(counts.log_counts(left, right, budget=budget))
    def score(values):
        return prefix+count+four.profile_transfer(values[0], values[1])
    baseline = score(regions[[left, right]])
    if abs(baseline-worst['individual_scores'][index]) > 1.e-6:
        raise ArithmeticError('V52 profile control')
    individual = np.array([exchange.marker_regional(table, mark)[[left, right]]
                           for table, mark in zip(bank.tables, odds)])
    winners = np.argmin(individual, axis=0)
    frequency = np.bincount(winners.ravel(), minlength=len(odds))
    selected = list(map(int, np.argsort(-frequency, kind='stable')[:4]))
    raw_pairs = four.compose(logs[:, None], logs[None])
    best = {name: np.full((2, 2, 2), np.inf) for name in
            ('pair_positive', 'quartet_positive', 'optimistic_selected_boundary_floor')}
    details = []
    for marker_index in selected:
        mark = odds[marker_index]
        weights = four.occupancy_weights(mark)
        selected_j = four.choose_occupancies(four.pair_terms(raw_pairs, weights, alpha), 16)
        row = four.positive_partition(raw_pairs, weights, alpha, selected_j, budget.check)
        c2, rest = row['selected_two'], row['remainder_two']
        c2square = four.compose(c2, c2)
        optimistic = four.matrix_sum([c2square+alpha-1,
            four.compose(c2, rest), four.compose(rest, c2), four.compose(rest, rest)])
        if np.any(optimistic > row['four']+1.e-7):
            raise ArithmeticError('Optimistic boundary lower diagnostic failed')
        for name, matrix, exponent in (('pair_positive', row['two'], 64),
                ('quartet_positive', row['four'], 32),
                ('optimistic_selected_boundary_floor', optimistic, 32)):
            candidate = np.array([four.regional(matrix, exponent, 4096, f, mark) for f in (left, right)])
            np.minimum(best[name], candidate, out=best[name])
        details.append(dict(marker_index=marker_index, marker=spec['markers'][marker_index],
            occupancies=list(selected_j), quartets_evaluated=row['quartets_evaluated'],
            selected_quartet_fraction=np.exp2(c2square-row['old_four']).tolist(),
            log2_gain_per_quartet=(row['old_four']-row['four']).tolist()))
    scores = {name: score(value) for name, value in best.items()}
    scores.update(alias_all_markers=baseline, alias_selected_markers=score(np.min(individual[selected], axis=0)))
    if scores['quartet_positive'] > scores['pair_positive']+1.e-6:
        raise ArithmeticError('Quartet refinement worsened the positive bound')
    if source_pins() != pins or (HERE/INPUT).read_bytes() != raw:
        raise RuntimeError('Source/input changed')
    result = dict(schema='translated-quartet-diagnostic-v53-FLOAT', source_sha256=pins,
        input=dict(path=INPUT, sha256=INPUT_SHA), selected_profile=[left, right],
        witness_index=index, z=spec['z'], alpha=spec['alpha'], scores=scores, markers=details,
        active_groups=4096, target_distance='9/20', conditional_on_principal_setup=True,
        chronological_quartet_before_alpha=True, omitted_occupancies_retained_positively=True,
        optimistic_floor_is_not_a_probability_upper_bound=True,
        all_profiles_replayed=False, directed=False, whole_code_certificate=False,
        setup_risk_added=False, elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Output cap')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, bytes=len(data), sha256=sha256(data).hexdigest(),
                scores=scores, elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), allow_nan=False), flush=True)
