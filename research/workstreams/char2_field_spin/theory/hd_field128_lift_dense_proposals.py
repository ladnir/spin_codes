"""Dense full-field FLOAT proposals with complete re-centered marker banks.

The source artifact supplies exact proposal literals and their calibration
histograms, never a probability bound. Each old marker is moved along the
entire feasible full-count interval. Boundary proxies are positive markers,
not assertions that a histogram has a positive empty category. Every FLOAT
evaluation still sums all matching profiles at its explicitly selected weight.
Exported witness artifacts are compatible with the frozen directed cover CLI.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import floor, isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.special import gammaln

import hd_field128_lift_alias_screen as alias
import hd_field128_lift_tail_screen as tail


CONSTRUCTION = 'toeplitz-fullfield-h4-alternating-swapped-pairs-k20-rate-half-v1'
OLD_NAME = 'hd_shifted_potential_band_3936_3769_4050.json'
OLD_SHA256 = '593f4db0f25aeda39c79fd26d92b565bbfcf4bcf55eb5884ad03e2ff170a2029'
DEFAULT_CANDIDATES = ((4428, 18), (4920, 20), (5412, 22), (5904, 24))
PAIRS = tail.PAIRS
ZERO_MARKER, EMPTY_PROXY = F(1, 1 << 64), F(1, 1 << 32)


def source_pins():
    pins = alias.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def exact(value):
    if isinstance(value, float):
        raise ValueError('Proposal parameters must have exact literals')
    return F(value)


def rounded(value, digits=17):
    if type(digits) is not int or not 6 <= digits <= 17:
        raise ValueError('Use 6..17 significant decimal digits')
    value = F(str(float(value)) if digits == 17 else format(float(value), f'.{digits}g'))
    if value <= 0:
        raise ValueError('Rounding lost a positive proposal')
    return value


def exact_ratio(log_ratio):
    """A short positive dyadic-scaled proposal, even outside float value range."""
    if not isfinite(log_ratio) or abs(log_ratio) > 4096:
        raise ValueError('Invalid log potential')
    exponent = floor(log_ratio)
    mantissa = F(format(2.**(log_ratio-exponent), '.12g'))
    return mantissa*(F(2)**exponent)


def log_rational(value):
    value = F(value)
    if value <= 0:
        raise ValueError('A positive rational is required')
    return log2(value.numerator)-log2(value.denominator)


def read_source(path=None):
    path = Path(path or Path(__file__).with_name(OLD_NAME)).resolve()
    raw = path.read_bytes()
    pin = dict(path=str(path), sha256=sha256(raw).hexdigest())
    if pin['sha256'] != OLD_SHA256:
        raise ValueError('The pinned old proposal artifact changed')
    artifact = json.loads(raw)
    fixed = artifact['fixed_witness']
    row = fixed['source_exact_witness']
    if (fixed.get('pair_swaps') is not True
            or fixed.get('orientation_aggregation') != 'independent-fair-swaps'
            or row['active_groups'] != 3936):
        raise ValueError('Unexpected old proposal semantics')
    markers = tuple(tuple(map(F, pair.split())) for pair in row['markers'])
    targets = tuple(map(F, row['central_targets']))
    if (len(markers) != len(targets) or not 2 <= len(markers) <= 100
            or any(len(pair) != 2 or min(pair) <= 0 for pair in markers)
            or targets != tuple(sorted(set(targets))) or targets[0] != 0 or targets[-1] != 1968):
        raise ValueError('Source markers need their full ordered calibration grid')
    return dict(weight=3936, z=F(row['z']), alpha=F(row['alpha']),
                ratio=F(row['potential']), grid=fixed['grid'], markers=markers,
                targets=targets, provenance=pin)


def histogram(weight, full):
    full = F(full)
    if type(weight) is not int or not 1 <= weight < 2*PAIRS:
        raise ValueError('Expected a nonterminal integer weight')
    result = (PAIRS-weight+full, weight-2*full, full)
    if min(result) < 0:
        raise ValueError('Infeasible calibration histogram')
    return result


def recentered_markers(source, center, *, empty_proxy=EMPTY_PROXY, zero_marker=ZERO_MARKER,
                       factors=(F(1), F(1))):
    """Retain all source calibrations; borrow endpoint tilts only as proposals."""
    empty_proxy, zero_marker = exact(empty_proxy), exact(zero_marker)
    factors = tuple(map(exact, factors))
    if empty_proxy <= 0 or zero_marker <= 0 or len(factors) != 2 or min(factors) <= 0:
        raise ValueError('Markers, proxies, and multipliers must be positive')
    old_weight, old_targets = source['weight'], source['targets']
    old_lo, old_hi = old_targets[0], old_targets[-1]
    lo, hi = F(max(0, center-PAIRS)), F(center, 2)
    histogram(center, lo)
    calibration = []
    for index, (target, marker) in enumerate(zip(old_targets, source['markers'])):
        e, s, f = histogram(old_weight, target)
        tilt = []
        for axis, category in enumerate((s, f)):
            if category:
                tilt.append(marker[axis]*e/category)
            else:
                interior = min((j for j, t in enumerate(old_targets)
                                if histogram(old_weight, t)[axis+1] > 0),
                               key=lambda j: abs(old_targets[j]-target))
                ee, ss, ff = histogram(old_weight, old_targets[interior])
                tilt.append(source['markers'][interior][axis]*ee/(ss if axis == 0 else ff))
        calibration.append(tuple(tilt))
    requests = [(lo+(target-old_lo)*(hi-lo)/(old_hi-old_lo), index, 'warped-original')
                for index, target in enumerate(old_targets)]
    for offset in (0, 1, 2, 4, 8, 16, 32):
        for target in (lo+offset, hi-offset):
            if lo <= target <= hi:
                old_target = old_lo+(target-lo)*(old_hi-old_lo)/(hi-lo)
                index = min(range(len(old_targets)), key=lambda j: abs(old_targets[j]-old_target))
                requests.append((target, index, 'endpoint-reinforcement'))
    markers, records, seen = [], [], set()
    for target, index, kind in requests:
        e, s, f = histogram(center, target)
        denominator = e if e else empty_proxy
        pair = tuple(factor*(kappa*category/denominator if category else zero_marker)
                     for factor, kappa, category in zip(factors, calibration[index], (s, f)))
        if pair in seen:
            continue
        seen.add(pair)
        markers.append(pair)
        records.append(dict(kind=kind, source_marker_index=index, target_full=str(target),
                            target_histogram=list(map(str, (e, s, f))),
                            uses_empty_proxy=not bool(e), uses_zero_marker=not bool(s*f)))
    if len(markers) > 128:
        raise ValueError('Dense bank exceeds the frozen FLOAT alias cap of128')
    return tuple(markers), records


def make_witness(source, center, exponent, *, grid=128, digits=17, factors=(F(1), F(1)),
                 empty_proxy=EMPTY_PROXY, alpha_factor=F(1)):
    if type(exponent) is not int or not 1 <= exponent <= 64 or grid not in (64, 128, 256):
        raise ValueError('Invalid mark exponent or alias grid')
    alpha_factor = exact(alpha_factor)
    z, alpha = rounded(source['z']**exponent, digits), rounded(source['alpha']/16*alpha_factor, digits)
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError('Proposal mark or fractional power is inadmissible')
    markers, records = recentered_markers(source, center, empty_proxy=empty_proxy, factors=factors)
    return dict(schema='field128-lift-alias-witness-only-v1', construction=CONSTRUCTION,
                field_order=tail.FIELD_ORDER, field_bits=128, pair_swaps=True,
                full_field_random_coefficients=True, directed=False, whole_code_certificate=False,
                status='untested exact proposal; no probability or interval coverage claim',
                source=source['provenance'], old_numerical_ceilings_reused=False,
                center=center, source_mark_exponent=exponent, source_marker_count=len(source['markers']),
                marker_count=len(markers), marker_calibrations=records,
                transform=dict(mark='rounded old_z**exponent', alpha='rounded old_alpha/16*alpha_factor',
                               significant_decimal_digits=digits, alpha_factor=str(alpha_factor),
                               markers='per-marker source histogram tilt, affine full-count warp, endpoint reinforcement',
                               empty_proxy=str(exact(empty_proxy)), zero_marker=str(ZERO_MARKER),
                               marker_factors=list(map(str, factors)), equality_claim=False),
                state_potential_ratio=str(source['ratio']),
                witness=dict(z=str(z), alpha=str(alpha), grid=grid, markers=[list(map(str, pair)) for pair in markers]))


def _selected_markers(bank, weight, budget):
    full = tail.full_counts(PAIRS, weight)
    single, empty = weight-2*full, PAIRS-weight+full
    normalizer = tail.LOG_FACTORIAL[PAIRS]-tail.LOG_FACTORIAL[empty]-tail.LOG_FACTORIAL[single]-tail.LOG_FACTORIAL[full]
    best = np.full((4, len(full)), np.inf)
    selected = np.full(best.shape, -1, dtype=int)
    reliable = np.zeros(best.shape, dtype=bool)
    for index, (odds, table) in enumerate(zip(bank.odds, bank.tables)):
        budget.check()
        raw = table.logs[:, single % bank.grid, full % bank.grid]
        valid = table.reliable[:, single % bank.grid, full % bank.grid]
        if bank.floor_policy == 'positive-fallback':
            raw = np.where(valid, raw, table.positive_numerator[:, None])
        candidate = raw-normalizer-single*odds[0]-full*odds[1]
        improve = candidate < best
        selected[improve], reliable[improve] = index, valid[improve]
        np.minimum(best, candidate, out=best)
    return full, best.T.reshape((-1, 2, 2)), selected, reliable, normalizer


def diagnose(bank, weight, log_ratio, *, budget):
    """Decompose the largest term of the same complete cyclic-run envelope."""
    full, regions, selected, reliable, normalizer = _selected_markers(bank, weight, budget)
    contraction = tail.contraction_logs(regions, log_ratio)
    length, runs = 2*PAIRS, np.arange(weight+1)
    counts = np.full(weight+1, tail.NEG)
    valid = (runs >= 1) & (runs <= min(weight, length-weight))
    c = runs[valid]
    counts[valid] = (np.log2(length/c)+(gammaln(weight)-gammaln(c)-gammaln(weight-c+1)
                    +gammaln(length-weight)-gammaln(c)-gammaln(length-weight-c+1))/tail.LN2)
    maximum, positions, count_bits = tail.NEG, None, None
    for first in range(0, len(full), 32):
        budget.check()
        run_count = weight-full[first:first+32, None]-full[None, :]
        terms = counts[run_count]+128*contraction[first:first+32, None]+128*contraction[None, :]
        local = np.unravel_index(np.argmax(terms), terms.shape)
        if float(terms[local]) > maximum:
            maximum = float(terms[local])
            positions = (first+local[0], local[1])
            count_bits = float(counts[run_count[local]])
    endpoint_details = []
    for position in positions:
        f = int(full[position])
        e, s, _ = histogram(weight, f)
        region = regions[position].ravel()
        state_rows = (float(np.logaddexp2(region[0], region[1]+log_ratio)),
                      float(np.logaddexp2(region[2]-log_ratio, region[3])))
        entries = []
        for entry in range(4):
            index = int(selected[entry, position])
            odds = bank.odds[index]
            entries.append(dict(entry=entry, marker_index=index, log2_normalized=float(region[entry]),
                                log2_histogram_count=float(normalizer[position]),
                                log2_marker_monomial=float(s*odds[0]+f*odds[1]),
                                alias_reliable=bool(reliable[entry, position])))
        endpoint_details.append(dict(full=f, single=int(s), empty=int(e),
                                     log2_contraction=float(contraction[position]),
                                     dominant_state_row=int(state_rows[1] > state_rows[0]), entries=entries))
    outer = bank.alpha*bank.density*weight
    threshold = -bank.alpha*bank.forbidden*log2(bank.z)
    projective = bank.constant-threshold
    front = max(0., -log_ratio)
    return dict(dominant_cyclic_run_summand=endpoint_details,
                decomposition=dict(outer_density=outer, forbidden_weight_mark=threshold,
                                   projective_divisor=projective, support_run_count=count_bits,
                                   regional_even=128*float(contraction[positions[0]]),
                                   regional_odd=128*float(contraction[positions[1]]), potential_front=front),
                log2_largest_summand_screen=outer+bank.constant+maximum+front,
                selected_below_floor_entries=int((~reliable).sum()),
                complete_histogram_count=len(full),
                interpretation='largest summand of a complete upper envelope; not an actual codeword or lower bound')


def screen_witness(witness, weights, *, floor_policy='positive-fallback', max_seconds=120.,
                   potential_evaluations=40, optimize_potential=True, progress=None):
    weights = tuple(dict.fromkeys(weights))
    if (not weights or len(weights) > 16 or any(type(r) is not int or not 1 <= r < 8192 for r in weights)
            or type(potential_evaluations) is not int or not 1 <= potential_evaluations <= 128):
        raise ValueError('Use1..16 explicit weights and a bounded potential optimizer')
    if (witness.get('construction') != CONSTRUCTION or witness.get('field_order') != tail.FIELD_ORDER
            or witness.get('pair_swaps') is not True or witness.get('full_field_random_coefficients') is not True):
        raise ValueError('A literal full-field construction witness is required')
    budget, pins = tail.Budget(max_seconds), source_pins()
    rows, failure = [], None
    try:
        w = witness['witness']
        markers = tuple(tuple(log_rational(exact(x)) for x in pair) for pair in w['markers'])
        bank = alias.AliasBankScreen(float(exact(w['z'])), float(exact(w['alpha'])), markers,
                                    log_rational(exact(witness['state_potential_ratio'])), grid=w['grid'],
                                    floor_policy=floor_policy, distance=F(7, 20), _budget=budget)
        for weight in weights:
            result = bank.evaluate(weight, method='runs', optimize_potential=optimize_potential,
                                   max_evaluations=potential_evaluations, _budget=budget)
            ratio = exact_ratio(result['log2_potential_ratio'])
            exported_log_ratio = log_rational(ratio)
            result['optimizer_log2_potential_ratio'] = result['log2_potential_ratio']
            support = tail.support_sum_logs(PAIRS, weight, bank.regional(weight), exported_log_ratio,
                                            method='runs', _budget=budget)
            result.update(log2_potential_ratio=exported_log_ratio, log2_support_sum_screen=support,
                          log2_upper_screen=result['log2_outer_prefix']+support)
            diagnostics = diagnose(bank, weight, exported_log_ratio, budget=budget)
            diagnostics['log2_sum_minus_largest_summand'] = result['log2_upper_screen']-diagnostics['log2_largest_summand_screen']
            export = dict(witness, state_potential_ratio=str(ratio), selected_weight=weight,
                          status='FLOAT-screened exact proposal; fresh directed replay required')
            result.update(diagnostics=diagnostics, directed_witness=export,
                          exported_ratio_is_new_exact_proposal=True)
            rows.append(result)
            if progress:
                progress(result)
    except (TimeoutError, ArithmeticError) as error:
        failure = dict(error_type=type(error).__name__, message=str(error))
    if source_pins() != pins:
        raise RuntimeError('A pinned dense proposal helper changed during execution')
    return dict(construction=CONSTRUCTION, requested_weights=list(weights), results=rows,
                pending_weights=[r for r in weights if r not in {row['active_groups'] for row in rows}],
                source_sha256=pins, failure=failure, floor_policy=floor_policy,
                elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                directed=False, whole_code_certificate=False,
                scope='explicit selected weights only; no interpolation or band coverage claim')


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(raw)
    return dict(path=str(Path(path).resolve()), sha256=sha256(raw).hexdigest())


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--candidate', type=int, nargs=2, action='append', metavar=('CENTER', 'EXPONENT'))
    parser.add_argument('--weight', type=int, action='append')
    parser.add_argument('--mode', choices=('plan', 'float'), default='plan')
    parser.add_argument('--grid', type=int, default=128)
    parser.add_argument('--digits', type=int, default=17)
    parser.add_argument('--alpha-factor', default='1')
    parser.add_argument('--x-factor', default='1')
    parser.add_argument('--y-factor', default='1')
    parser.add_argument('--empty-proxy', default=str(EMPTY_PROXY))
    parser.add_argument('--floor-policy', choices=('positive-fallback', 'legacy-cushion'), default='positive-fallback')
    parser.add_argument('--max-seconds', type=float, default=300.)
    parser.add_argument('--max-candidate-seconds', type=float, default=120.)
    parser.add_argument('--potential-evaluations', type=int, default=40)
    parser.add_argument('--fixed-potential', action='store_true')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory must be unused')
    candidates = args.candidate or DEFAULT_CANDIDATES
    if not 1 <= len(candidates) <= 8 or not isfinite(args.max_candidate_seconds) or not 0 < args.max_candidate_seconds <= 3600:
        parser.error('At most8 candidates and bounded positive per-candidate time required')
    budget = tail.Budget(args.max_seconds)
    source = read_source(args.source)
    proposals = [make_witness(source, center, exponent, grid=args.grid, digits=args.digits,
                              factors=(F(args.x_factor), F(args.y_factor)), empty_proxy=F(args.empty_proxy),
                              alpha_factor=F(args.alpha_factor)) for center, exponent in candidates]
    args.output_dir.mkdir(exist_ok=False)
    receipts, attempted, completed, failure = [], [], [], None
    for index, proposal in enumerate(proposals):
        write_new(args.output_dir/f'witness_{index:02d}.json', proposal)
    for index, proposal in enumerate(proposals if args.mode == 'float' else ()):
        try:
            budget.check()
            seconds = min(args.max_candidate_seconds, args.max_seconds-(monotonic()-budget.started))
            if seconds <= 0:
                raise TimeoutError('Dense proposal deadline exhausted before bank construction')
            weights = args.weight or [min(5853, max(4068, proposal['center']))]
            result = screen_witness(proposal, weights, floor_policy=args.floor_policy, max_seconds=seconds,
                                    potential_evaluations=args.potential_evaluations,
                                    optimize_potential=not args.fixed_potential,
                                    progress=lambda row: print(json.dumps({key: row[key] for key in
                                               ('active_groups', 'log2_upper_screen', 'diagnostics')}), flush=True))
            receipts.append(write_new(args.output_dir/f'float_{index:02d}.json', result))
            for row in result['results']:
                write_new(args.output_dir/f'witness_{index:02d}_r{row["active_groups"]}.json', row['directed_witness'])
            attempted.append(index)
            if not result['pending_weights']:
                completed.append(index)
            if result['failure']:
                failure = dict(candidate=index, detail=result['failure'])
                break
        except (TimeoutError, ArithmeticError) as error:
            failure = dict(candidate=index, error_type=type(error).__name__, message=str(error))
            break
    summary = dict(mode=args.mode, candidate_centers=[p['center'] for p in proposals],
                   requested_weights_by_candidate=[args.weight or [min(5853, max(4068, p['center']))] for p in proposals],
                   attempted_candidates=attempted, completed_candidates=completed,
                   pending_candidate_indices=[i for i in range(len(proposals)) if i not in completed],
                   receipts=receipts, failure=failure,
                   source_sha256=source_pins(), directed=False, whole_code_certificate=False,
                   scope='proposal generation and selected-weight FLOAT diagnosis only')
    write_new(args.output_dir/'summary.json', summary)
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
