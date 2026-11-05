"""Full-occupancy quartet FLOAT diagnostic at one pinned V52 profile.

Each chronological four-packet product is formed before entrywise alpha.
Histogram multiplicities stay outside alpha. All 65**4 ordered occupancy
quartets are retained; no remainder or coefficient truncation is used.
The boundary floor is a lower diagnostic, not a failure-probability bound.
This is not a directed, selected-weight, or whole-code certificate.
"""
import argparse
import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite, log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_translated_quartet_diagnostic_v53 as previous
import hd_field128_gamma_four_packet_positive_v35 as four
import hd_field128_reversed_marker_exchange_v49 as guard

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = '21c1e20038101526bc1e41daf0d869deeda3565c620ac6acaf1e6456a9e37c10'
GUARD_SHA = 'bb9d01af9135af040a1b86ea8efae5e02acb6fc932070485d2ea1c03df8c4866'
INPUT53 = 'hd_field128_translated_quartet_diagnostic_v53_20261010.json'
INPUT53_SHA = '993bd2d5bdac3aa1f62c2e900810e3c479be2fe65d596980ae68785a614695d8'
INPUT52 = 'hd_field128_reversed_dual_marker_v52_20261010.json'
INPUT52_SHA = '1ec14fce6aede5f7f45d301fdeae5cb7057670a09f8e244bad52e9bd1b69b3e5'
LIMIT, RESERVE = 32 << 10, 20 << 30
PROFILE, WEIGHT, WINNER = (1613, 1727), 4096, 0


def source_pins():
    """Pin every local Python import, plus the frozen helpers' note closure."""
    pins = previous.source_pins()
    if (pins.get(Path(previous.__file__).name) != PREVIOUS_SHA
            or pins.get(Path(guard.__file__).name) != GUARD_SHA):
        raise RuntimeError('Frozen V53 or V49 changed')
    pending, seen = [Path(__file__).name], set()
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        raw = (HERE/name).read_bytes()
        digest = sha256(raw).hexdigest()
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting dependency: '+name)
        pins[name] = digest
        for node in ast.walk(ast.parse(raw)):
            modules = ([x.name for x in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                child = module.split('.')[0]+'.py'
                if (HERE/child).is_file() and child not in seen:
                    pending.append(child)
    return dict(sorted(pins.items()))


def quartet_count(size):
    if type(size) is not int or not 1 <= size <= 65:
        raise ValueError('Occupancy table size must be 1..65')
    return size**4


def full_quartets(raw_pairs, weights, alpha, *, chunk=32, check=lambda: None):
    """Return log matrices for 1..4 markers, streaming all ordered quartets.

    raw_pairs[j,k] is log2(P_j P_k), before alpha. Each row of weights is
    log2(g_j(x,y)) for one marker. Outputs have shape (markers,2,2).
    Zero entries are represented by -inf. The largest raw quartet block
    has shape (chunk,size**2,2,2); no size**4 tensor is allocated.
    """
    import numpy as np
    raw_pairs, weights = np.asarray(raw_pairs, dtype=float), np.asarray(weights, dtype=float)
    if raw_pairs.ndim != 4 or raw_pairs.shape[2:] != (2, 2):
        raise ValueError('Expected an occupancy-pair 2x2 log table')
    size = raw_pairs.shape[0]
    expected = quartet_count(size)
    if weights.ndim == 1:
        weights = weights[None, :]
    if (raw_pairs.shape[1] != size or weights.ndim != 2 or weights.shape[1] != size
            or not 1 <= weights.shape[0] <= 4 or type(chunk) is not int or not 1 <= chunk <= 32
            or not isfinite(alpha) or not 0 < alpha <= 1
            or np.any(np.isnan(raw_pairs)) or np.any(np.isposinf(raw_pairs))
            or np.any(np.isnan(weights)) or np.any(np.isposinf(weights))):
        raise ValueError('Invalid markers, chunk, alpha, or positive-log input')
    check()
    pairs = raw_pairs.reshape((-1, 2, 2))
    pair_weights = (weights[:, :, None]+weights[:, None, :]).reshape((len(weights), -1))
    two = np.array([four.matrix_sum(pair_weights[i, :, None, None]+alpha*pairs)
                    for i in range(len(weights))])
    result = np.full((len(weights), 2, 2), -np.inf)
    seen, blocks = 0, 0
    for start in range(0, len(pairs), chunk):
        check()
        stop = min(start+chunk, len(pairs))
        # The left pair precedes the right pair. Alpha acts on their full
        # two-state path sum, never on the multiplicities or separate pairs.
        raw = four.compose(pairs[start:stop, None], pairs[None])
        raw *= alpha
        for marker in range(len(weights)):
            check()
            values = (raw+pair_weights[marker, start:stop, None, None, None]
                      +pair_weights[marker, None, :, None, None])
            result[marker] = np.logaddexp2(result[marker], four.matrix_sum(values))
            del values
        seen += (stop-start)*len(pairs)
        blocks += 1
        del raw
    check()
    if seen != expected:
        raise ArithmeticError('Incomplete ordered quartet enumeration')
    old_four = four.compose(two, two)
    floor = old_four+alpha-1
    finite = np.isfinite(old_four)
    tolerance = 2.e-10*max(1., float(np.max(np.abs(old_four[finite]), initial=0.)))
    if (np.any(np.isnan(result)) or np.any(np.isposinf(result))
            or np.any(result[~finite] != -np.inf)
            or np.any(result[finite] > old_four[finite]+tolerance)
            or np.any(result[finite] < floor[finite]-tolerance)):
        raise ArithmeticError('Full quartet boundary inequalities failed')
    return dict(two=two, four=result, old_four=old_four, optimistic_floor=floor,
                ordered_quartets_per_marker=seen, raw_quartets_composed=seen,
                pair_blocks=blocks, chunk=chunk, markers=len(weights))


def validate_inputs(v52, v53):
    if (v52.get('schema') != 'reversed-witness-exchange-v52-FLOAT'
            or v53.get('schema') != 'translated-quartet-diagnostic-v53-FLOAT'):
        raise ValueError('Wrong diagnostic input schema')
    for item in (v52, v53):
        if (item.get('active_groups') != WEIGHT or item.get('target_distance') != '9/20'
                or item.get('conditional_on_principal_setup') is not True
                or item.get('directed') is not False or item.get('whole_code_certificate') is not False
                or item.get('setup_risk_added') is not False):
            raise ValueError('Wrong original-weight or conditional construction scope')
    worst = v52['final']['largest_profile']
    if (v52.get('all_profiles_retained') is not True
            or v52['final'].get('profiles_evaluated') != 2049**2
            or (worst.get('even_full'), worst.get('odd_full')) != PROFILE
            or worst.get('witness') != WINNER or v53.get('selected_profile') != list(PROFILE)
            or v53.get('witness_index') != WINNER
            or v53.get('input') != dict(path=INPUT52, sha256=INPUT52_SHA)):
        raise ValueError('Wrong pinned profile or winning witness')
    spec = v52['witness_specs'][WINNER]
    z, alpha = F(spec['z']), F(spec['alpha'])
    if (not 0 < z < 1 or not 0 < alpha <= 1 or F(float(z)) != z or F(float(alpha)) != alpha
            or (v53.get('z'), v53.get('alpha')) != (spec['z'], spec['alpha'])):
        raise ValueError('Exact witness marks differ')
    markers = spec['markers']
    if not 4 <= len(markers) <= 128:
        raise ValueError('Expected 4..128 positive marker pairs')
    for pair in markers:
        if len(pair) != 2 or any(F(v) <= 0 or not isfinite(float(F(v))) or float(F(v)) <= 0 for v in pair):
            raise ValueError('Invalid positive marker literal')
    selected = [row['marker_index'] for row in v53['markers']]
    if (len(selected) != 4 or len(set(selected)) != 4
            or any(type(i) is not int or not 0 <= i < len(markers) for i in selected)
            or any(row['marker'] != markers[i] for row, i in zip(v53['markers'], selected))):
        raise ValueError('V53 selected-marker literals differ')
    return spec, selected


def close_control(actual, expected, label):
    if not isfinite(actual) or not isfinite(expected) or abs(actual-expected) > 2.e-6:
        raise ArithmeticError('Saved FLOAT control differs: '+label)


def write_result(output, result):
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('32KiB output or planned 20GiB reserve exceeded')
    with output.open('xb') as stream:
        stream.write(data)
    if output.read_bytes() != data:
        raise RuntimeError('Exclusive output changed')
    return dict(path=output.relative_to(HERE).as_posix(), bytes=len(data), sha256=sha256(data).hexdigest())


def run(output_name):
    # V49: Linux-only, 1GiB address space, 175/180 CPU seconds, 180s alarm.
    # Its 96KiB file limit is supplemented by our 32KiB serialized cap.
    with guard.remote_limits():
        return _run(output_name)


def _run(output_name):
    import numpy as np
    import hd_field128_reversed_return_v44 as caps
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    exchange = previous.previous.exchange
    old = caps.prior.control.old
    counts_module = exchange.previous.previous.counts_module
    output = old.prior.previous.inside(output_name)
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Exclusive output and planned reserve required')
    pins, raw_inputs = source_pins(), {}
    for name, digest in ((INPUT52, INPUT52_SHA), (INPUT53, INPUT53_SHA)):
        path = HERE/name
        if path.stat().st_size > 96 << 10:
            raise ValueError('Oversized input')
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Pinned input changed: '+name)
        raw_inputs[name] = raw
    def check_inputs():
        if source_pins() != pins or any((HERE/name).read_bytes() != raw for name, raw in raw_inputs.items()):
            raise RuntimeError('Source or input changed')
    try:
        v52, v53 = (json.loads(raw_inputs[name]) for name in (INPUT52, INPUT53))
        if (v52.get('source_sha256') != previous.previous.source_pins()
                or v53.get('source_sha256') != previous.source_pins()):
            raise ValueError('Saved source closures differ')
        spec, selected = validate_inputs(v52, v53)
        z, alpha = float(F(spec['z'])), float(F(spec['alpha']))
        started, budget = monotonic(), gamma.old.Budget(165.)
        logs = caps.local_logs(z)
        coefficient = old.Coefficients(logs, z, alpha, budget)
        odds = [tuple(log2(float(F(v))) for v in pair) for pair in spec['markers']]
        bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128,
            distance=F(9, 20), grid=128, floor_policy='positive-fallback',
            max_table_bytes=384 << 20, _coefficients=coefficient, _budget=budget)
        if bank.odds != tuple(odds):
            raise ArithmeticError('Marker deduplication changed table indices')
        regions = bank.regional(WEIGHT)
        if (regions.shape != (2049, 2, 2)
                or sha256(regions.tobytes()).hexdigest() != v52['regional_sha256'][WINNER]):
            raise ArithmeticError('Pinned original-r regional control differs')
        counts = counts_module.MidpointCounts(4096, WEIGHT, max_singletons=64, budget=budget)
        count = float(counts.log_counts(*PROFILE, budget=budget))
        prefix = bank.constant+bank.alpha*bank.density*WEIGHT
        def score(values):
            return prefix+count+four.profile_transfer(values[0], values[1])
        baseline = score(regions[list(PROFILE)])
        close_control(baseline, v52['final']['largest_profile']['individual_scores'][WINNER], 'V52 winner')
        individual = np.array([exchange.marker_regional(table, mark)[list(PROFILE)]
                               for table, mark in zip(bank.tables, odds)])
        winners = np.argmin(individual, axis=0)
        frequency = np.bincount(winners.ravel(), minlength=len(odds))
        if list(map(int, np.argsort(-frequency, kind='stable')[:4])) != selected:
            raise ArithmeticError('V53 marker selection changed')
        scores = dict(alias_all_markers=baseline,
                      alias_selected_markers=score(np.min(individual[selected], axis=0)))
        # Keep only the two profile rows and four marker literals. No atlas
        # is needed during the full quartet loop.
        del bank, coefficient, regions, individual, counts
        raw_pairs = four.compose(logs[:, None], logs[None])
        best = {name: np.full((2, 2, 2), np.inf) for name in
                ('pair_positive', 'quartet_positive', 'optimistic_selected_boundary_floor',
                 'full_quartet_positive', 'optimistic_full_boundary_floor')}
        all_weights, details = [], []
        for marker_index, saved in zip(selected, v53['markers']):
            budget.check()
            mark = odds[marker_index]
            weights = four.occupancy_weights(mark)
            chosen = four.choose_occupancies(four.pair_terms(raw_pairs, weights, alpha), 16)
            row = four.positive_partition(raw_pairs, weights, alpha, chosen, budget.check)
            c2, rest = row['selected_two'], row['remainder_two']
            c2square = four.compose(c2, c2)
            optimistic = four.matrix_sum([c2square+alpha-1, four.compose(c2, rest),
                                         four.compose(rest, c2), four.compose(rest, rest)])
            if list(chosen) != saved['occupancies'] or row['quartets_evaluated'] != saved['quartets_evaluated']:
                raise ArithmeticError('V53 selected occupancy control differs')
            for name, matrix, exponent in (('pair_positive', row['two'], 64),
                    ('quartet_positive', row['four'], 32),
                    ('optimistic_selected_boundary_floor', optimistic, 32)):
                candidate = np.array([four.regional(matrix, exponent, WEIGHT, f, mark) for f in PROFILE])
                np.minimum(best[name], candidate, out=best[name])
            all_weights.append(weights)
            details.append(dict(marker_index=marker_index, marker=spec['markers'][marker_index],
                selected_occupancies=list(chosen), selected16_log_four=row['four'].tolist()))
        for name in ('pair_positive', 'quartet_positive', 'optimistic_selected_boundary_floor'):
            scores[name] = score(best[name])
        for name, expected in v53['scores'].items():
            close_control(scores[name], expected, 'V53 '+name)
        full = full_quartets(raw_pairs, np.asarray(all_weights), alpha, chunk=32, check=budget.check)
        for i, marker_index in enumerate(selected):
            mark = odds[marker_index]
            selected_four = np.asarray(details[i]['selected16_log_four'])
            if np.any(full['four'][i] > selected_four+1.e-7):
                raise ArithmeticError('Full quartet failed to dominate selected16')
            for name, matrix in (('full_quartet_positive', full['four'][i]),
                                 ('optimistic_full_boundary_floor', full['optimistic_floor'][i])):
                candidate = np.array([four.regional(matrix, 32, WEIGHT, f, mark) for f in PROFILE])
                np.minimum(best[name], candidate, out=best[name])
            details[i].update(full_log_four=full['four'][i].tolist(),
                pair_log_two=full['two'][i].tolist(), optimistic_full_log_four=full['optimistic_floor'][i].tolist(),
                gain_bits_per_quartet=(full['old_four'][i]-full['four'][i]).tolist())
        scores.update({name: score(best[name]) for name in
                       ('full_quartet_positive', 'optimistic_full_boundary_floor')})
        if (not all(map(isfinite, scores.values()))
                or scores['full_quartet_positive'] > scores['quartet_positive']+2.e-6
                or scores['optimistic_full_boundary_floor'] > scores['full_quartet_positive']+2.e-6):
            raise ArithmeticError('Invalid full-quartet comparison')
        budget.check()
        result = dict(schema='translated-full-quartet-diagnostic-v54-FLOAT', source_sha256=pins,
            inputs=[dict(path=n, sha256=sha256(raw).hexdigest()) for n, raw in raw_inputs.items()],
            construction='A128 emission,A64-transpose feedback,state32; uniform global degree16 gamma',
            active_groups=WEIGHT, target_distance='9/20', selected_profile=list(PROFILE),
            witness_index=WINNER, z=spec['z'], alpha=spec['alpha'], scores=scores, markers=details,
            log2_outer_prefix=prefix, log2_count_upper=count, count_method='MidpointCounts(max_singletons=64)',
            counts_outside_alpha=True, original_r_for_normalization_and_prefix=True,
            chronological_quartet_before_alpha=True, all_occupancies_retained=True,
            remainder_used=False, coefficient_truncation=False, pair_slots_per_packet=32,
            packets_per_region=128, four_blocks_per_region=32, regions=256,
            ordered_quartets_per_marker=full['ordered_quartets_per_marker'],
            raw_quartets_composed_once=full['raw_quartets_composed'], pair_blocks=full['pair_blocks'],
            pair_block_limit=full['chunk'], selected_markers=len(selected),
            pair_to_full_quartet_gain_bits=scores['pair_positive']-scores['full_quartet_positive'],
            full_quartet_gain_ceiling_bits=8192*(1-alpha),
            optimistic_floor_is_not_a_probability_upper_bound=True,
            conditional_on_principal_setup=True,
            principal_bad_gamma_upper=str(F(8*comb(64, 32), (1 << 128)-(1 << 64))),
            setup_risk_added=False, all_profiles_replayed=False, quartet_fft_built=False,
            directed=False, selected_weight_certificate=False, whole_code_certificate=False,
            arithmetic='FLOAT only; log summation rounding is not a directed enclosure',
            elapsed_seconds=monotonic()-started,
            limits=dict(wall_seconds=180, cooperative_seconds=165, address_space_bytes=1 << 30,
                        output_bytes=LIMIT, minimum_free_after_write=RESERVE))
        check_inputs()
        receipt = write_result(output, result)
        return dict(receipt=receipt, scores=scores, elapsed_seconds=result['elapsed_seconds'],
                    directed=False, all_profiles_replayed=False)
    finally:
        check_inputs()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), allow_nan=False), flush=True)
