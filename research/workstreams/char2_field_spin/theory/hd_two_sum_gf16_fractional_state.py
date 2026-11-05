"""Bounded FLOAT fractional-mask all-cohort G16 state-envelope control.

Condition on routing, average outer amplitudes and inner setup first, then
apply min(1,S)<=S**alpha. Subadditivity is used over mask words and state
paths. Thus packet caps are powered BEFORE route-geometry averaging; route
probabilities themselves are not powered. Dense weights are fresh moments
sum_A N_A(T)**(16*alpha), not powers of the alpha1 Holder weights.

The fixed-type Holder/Minkowski and chronological path enumeration are
reused unchanged from state_envelope. Nothing here is a directed receipt.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_closed_kernel import (branch_counts, choose, zero_one_moment,
                                     zero_two_moments)
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_two_sum_gf16_state_envelope import (affine_envelope, elementary_logs,
    kernel_from_types, logsum, matrix_power_log, path_classes, regional_curves,
    slope_candidates)
from hd_shifted_swap_sparse import power_screen, sparse_geometry, valid_degrees
from hd_xor_multigroup import LIMIT, nonzero_matrices, profile_arrays


CONTROLS = ((F(99, 100), F(1)), (F(24, 25), F(1, 2)),
            (F(9, 10), F(1, 4)), (F(3, 4), F(1, 8)),
            (F(1, 2), F(1, 16)))
FIRST, LAST, COLUMNS, REGIONS = 18, 32, 16, 16


def checked_alpha(alpha):
    alpha = F(alpha)
    if not 0 < alpha <= 1:
        raise ValueError('A positive fractional exponent at most1 is required')
    return alpha


def floor_array(log_floors, group=16):
    values = np.asarray(log_floors, dtype=float)
    if values.ndim == 0:
        values = np.full(group, float(values))
    if values.shape != (group,) or np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError('One fixed nonnegative log floor marker for each type is required')
    return values


def floor_templates(mark, alpha):
    """Five deterministic per-type templates, never adapted to a state path.

    At each anchor, Q pays the prefix/source count plus60 bits and16 types.
    Every binary64 log is frozen as its exact rational value. These values
    are proposal witnesses, not directed arithmetic or optimized markers.
    """
    alpha = checked_alpha(alpha)
    prefix = float(alpha)*(-LIMIT*log2(float(mark))-log2(255))
    empty = np.array([log2(comb(16, a))-8*a*float(alpha)*log2(255) for a in range(1, 17)])
    floors = np.array([9, 18]+[9]*14)
    result = [dict(name='all-zero', anchor=None, scale=0, log2_floor_markers=[F(0)]*16)]
    for anchor in (FIRST, LAST):
        q = (prefix+log2(comb(8192, anchor))+60)/anchor+4
        for scale in (1, 2):
            values = scale*np.maximum(0., (q+empty)/floors)
            result.append(dict(name=f'anchor{anchor}-scale{scale}', anchor=anchor, scale=scale,
                               log2_floor_markers=[F(float(value)) for value in values]))
    return result


def all_empty_log(mark, alpha, log_floors, active_groups):
    """Admitted all-empty-mask term BEFORE regional envelope intercept loss."""
    alpha = checked_alpha(alpha)
    markers = floor_array(log_floors)
    empty = np.array([log2(comb(16, a))-8*a*float(alpha)*log2(255) for a in range(1, 17)])
    floor_counts = np.array([9, 18]+[9]*14)
    prefix = float(alpha)*(-LIMIT*log2(float(mark))-log2(255))
    return float(prefix+log2(comb(8192, active_groups))
                 +active_groups*logsum(empty-floor_counts*markers))


def fractional_packet_coefficients(table, geometry, alpha, maximum, marker):
    """Positive packet coefficients; geometry weights stay outside the power."""
    alpha = checked_alpha(alpha)
    if type(maximum) is not int or not 1 <= maximum <= 64 or not 0 < marker <= 1:
        raise ValueError('Maximum1..64 and positive coefficient marker required')
    result = np.zeros((4, maximum+1, maximum//2+1))
    powered = {}
    for single, full, count, occupancy in geometry:
        if (not 0 <= single+2*full <= maximum or min(single, full) < 0 or count <= 0):
            raise ValueError('Positive complete packet-geometry counts required')
        if occupancy not in powered:
            matrix = np.asarray(table[occupancy])
            if matrix.shape != (4,) or np.any(~np.isfinite(matrix)) or np.any(matrix < 0):
                raise ValueError('Four nonnegative finite packet caps required')
            powered[occupancy] = matrix**float(alpha)
        result[:, single, full] += ((count/2**single)*marker**(single+full)
                                    *powered[occupancy])
    return result


def normalize_regional_coefficients(powered, maximum, marker, pairs=4096):
    """The exact support-geometry count is NOT raised to alpha."""
    if (np.shape(powered) != (4, maximum+1, maximum//2+1)
            or np.any(~np.isfinite(powered)) or np.any(np.asarray(powered) < 0)
            or not 0 < marker <= 1 or type(pairs) is not int or pairs < maximum):
        raise ValueError('Finite nonnegative regional coefficients required')
    result = {}
    with np.errstate(divide='ignore'):
        for single, full in valid_degrees(maximum):
            normalizer = (log2(comb(pairs, single)*comb(pairs-single, full))
                          +(single+full)*log2(marker))
            result[single+2*full, full] = (
                np.log2(powered[:, single, full])-normalizer).reshape(2, 2)
    return result


def fractional_regions(z, alpha, maximum=LAST):
    """All NZ-amplitude R_alpha[m,f], using M_packet**alpha before averaging."""
    alpha = checked_alpha(alpha)
    if not 0 < z < 1 or type(maximum) is not int or not 1 <= maximum <= 64:
        raise ValueError('Mark in(0,1) and maximum1..64 required')
    _, _, matrices = nonzero_matrices(float(z))
    occupancies, _ = profile_arrays()
    table = {tuple(occ): matrix.ravel() for occ, matrix in zip(occupancies, matrices)}
    marker = maximum/4096
    base = fractional_packet_coefficients(table, sparse_geometry(maximum), alpha, maximum, marker)
    return normalize_regional_coefficients(power_screen(base, 128, maximum), maximum, marker)


def fractional_power_sum_log(arity, weight, alpha, *, group=16, columns=16):
    """log2 sum_{|A|=arity} N_A(T)**(columns*alpha), including overlap first.

    This evaluates the proved closed counting kernel, not an interpolation
    in alpha. Positive real powers are allowed; every branch base and
    multiplicity is an exact integer. Zero bases contribute exactly zero.
    """
    alpha = checked_alpha(alpha)
    if (group not in (8, 16) or type(columns) is not int or columns < 1
            or type(arity) is not int or type(weight) is not int
            or not 0 <= arity <= group or not 0 <= weight <= group):
        raise ValueError('Checked group, arity, target weight and column count required')
    power = columns*float(alpha)
    if arity == group:
        return power*weight*log2(255)
    if arity == group-1:
        # Quiet/hole coincide for a missing output coordinate; hole/full
        # coincide for a present coordinate. Add those counts before pow.
        terms = [(group-weight, zero_one_moment(weight))]
        if weight:
            terms.append((weight, 255**(weight-1)-zero_one_moment(weight-1)))
    else:
        terms = [(choose(group-weight, arity-weight), zero_two_moments(weight))]
        inside = arity-group+weight
        if 0 <= inside <= arity:
            terms.append((choose(weight, inside), branch_counts(arity, inside, group)[2]))
        inside += 1
        if 0 <= inside <= arity:
            terms.append(((group-weight)*choose(weight, inside),
                          branch_counts(arity, inside, group)[1]))
    logs = [log2(count)+power*log2(value) for count, value in terms if count and value]
    return float(logsum(logs)) if logs else -np.inf


def fractional_dense_weights(alpha):
    """W_a,alpha(T)^16=sum_A N_A(T)^(16alpha)/255^(8a alpha)."""
    alpha = checked_alpha(alpha)
    return {arity: np.array([(fractional_power_sum_log(arity, weight, alpha)
                             -8*arity*float(alpha)*log2(255))/16
                            for weight in range(17)]) for arity in range(3, 17)}


def fractional_type_weights(log_marks, dense_weights, log_floors, alpha, *,
                            columns=16, dimension=8, minimum=9, q=256):
    """Fixed-type column masses after fractionalization and floor marking.

    The two singleton source categories remain separate. Double full-mask
    mass is (q-1)*(q-1-g), and each hole has mass(q-1), before alpha.
    Row choices are unpowered. Floor markers are introduced after alpha
    and remain fixed across every path and every outer column.
    """
    alpha = checked_alpha(alpha)
    log_marks = np.asarray(log_marks)
    group = log_marks.shape[-1]
    log_floors = floor_array(log_floors, group)
    if (not 2 <= group < q-1 or set(dense_weights) != set(range(3, group+1))
            or type(columns) is not int or not 1 <= minimum <= columns
            or not 1 <= dimension <= columns):
        raise ValueError('Complete type weights and fixed nonnegative floor logs required')
    alpha = float(alpha)
    elementary = elementary_logs(log_marks)
    full, logq = elementary[..., group], log2(q-1)
    single_floor, pair_floor = log_floors[:2]
    single = (log2(group)-dimension*alpha*logq-minimum*single_floor)/columns
    single += np.logaddexp2(0., alpha*logq+single_floor+full)
    pair_terms = np.stack((np.zeros_like(full), 1+alpha*logq+pair_floor+full,
        2*pair_floor+np.logaddexp2(alpha*(logq+log2(q-1-group))+full,
                                  alpha*logq+elementary[..., group-1])), axis=-1)
    pair = (log2(comb(group, 2))-2*dimension*alpha*logq-2*minimum*pair_floor)/columns
    pair += logsum(pair_terms, axis=-1)
    values = [single, pair]
    for arity in range(3, group+1):
        log_floor = log_floors[arity-1]
        weights = np.asarray(dense_weights[arity])
        if weights.shape != (group+1,) or np.any(np.isnan(weights)) or np.any(np.isposinf(weights)):
            raise ValueError('Complete finite-or-zero dense cardinality weights required')
        value = np.logaddexp2(weights[0], log_floor+logsum(
            elementary[..., 1:]+weights[1:], axis=-1))
        values.append(value-minimum*log_floor/columns)
    return np.stack(values, axis=-1)


def evaluate(curves, weights, mark, alpha, log_floors, slopes):
    alpha = checked_alpha(alpha)
    intercepts, slack = affine_envelope(curves, slopes)
    classes = path_classes(REGIONS)
    counts = np.asarray([row[2] for row in classes])
    marks = np.asarray([np.repeat(slopes, row[2]) for row in classes])
    typed = fractional_type_weights(marks, weights, log_floors, alpha)
    path_a = counts @ intercepts
    prefix = float(alpha)*(-LIMIT*log2(float(mark))-log2(255))
    entries = []
    for r in range(FIRST, LAST+1):
        kernel = kernel_from_types(classes, path_a, typed, r)
        power = matrix_power_log(kernel, COLUMNS)
        entries.append(dict(active_groups=r, failure_log2_screen=float(
            prefix+log2(comb(8192, r))+logsum(power[0])),
            all_empty_mask_log2_screen=all_empty_log(mark, alpha, log_floors, r)))
    return dict(mark=str(F(mark)), alpha=str(alpha),
        log2_floor_markers=[str(F(float(value))) for value in floor_array(log_floors)],
        log2_intercepts=intercepts.tolist(), log2_slopes=np.asarray(slopes).tolist(),
        envelope_contacts=[np.flatnonzero(row < 1e-8).tolist() for row in slack],
        envelope_minimum_slack=float(np.min(slack)),
        zero_occupancy_envelope_log2=intercepts.tolist(), entries=entries,
        failure_log2_screen=float(logsum([row['failure_log2_screen'] for row in entries])))


def refine_one_round(curves, weights, baseline, check_time):
    """Only the globally best baseline gets one bounded coordinate round."""
    _, candidates = slope_candidates(curves)
    best = baseline
    memo = {tuple(baseline['log2_slopes']): baseline}
    for index in range(4):
        slopes = np.asarray(best['log2_slopes'])
        for candidate in candidates[index]:
            check_time()
            trial = slopes.copy()
            trial[index] = candidate
            key = tuple(trial)
            if key not in memo:
                memo[key] = evaluate(curves, weights, F(baseline['mark']), F(baseline['alpha']),
                                     list(map(F, baseline['log2_floor_markers'])), trial)
            row = memo[key]
            if row['failure_log2_screen'] < best['failure_log2_screen']:
                best = row
    return dict(best, coordinate_rounds=1, slope_evaluations=len(memo),
        slope_candidates=[list(row) for row in candidates],
        slope_trials=[dict(log2_slopes=row['log2_slopes'], failure_log2_screen=row['failure_log2_screen'])
                      for row in memo.values()])


def run(*, max_seconds=30, refine=True, progress=None):
    if (type(max_seconds) not in (int, float) or not np.isfinite(max_seconds)
            or not 0 < max_seconds <= 30 or type(refine) is not bool):
        raise ValueError('Explicit budget at most30s and Boolean refinement flag required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start >= max_seconds:
            raise TimeoutError('Fractional all-cohort FLOAT control exhausted its budget')
    result = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        forbidden_weight=LIMIT, active_group_band=[FIRST, LAST], all_cohort_splits=True,
        all_constituent_row_sets=True, all_source_supports=True, source_counts='C(8192,r)',
        total_columns=COLUMNS, regions_per_column=REGIONS, regional_amplitudes='uniform-nonzero',
        pair_swaps=True, fixed_route_before_fractionalization=True,
        packet_power_before_geometry_average=True, route_normalizer_powered=False,
        outer_and_inner_average_before_fractionalization=True,
        fractional_subadditivity_over_mask_words_and_state_paths=True,
        dense_weights='(sum_A N_A(T)^(16*alpha)/255^(8*a*alpha))^(1/16)',
        dense_alpha1_weights_powered=False, rowset_choices_powered=False,
        projective_normalizations=1, projective_inside_alpha=True,
        endpoint_potential_factors=0, common_potential_used=False,
        initial_state=0, terminal_vector=[1, 1], fixed_types_across_columns=True,
        type_markers_fixed_across_paths_and_columns=True, floor_markers_after_fractionalization=True,
        singleton_row_floor=9, pair_aggregate_row_floor=18, dense_nonempty_column_floor=9,
        fixed_controls=[[str(z), str(a)] for z, a in CONTROLS],
        floor_template_names=['all-zero', 'anchor18-scale1', 'anchor18-scale2', 'anchor32-scale1', 'anchor32-scale2'],
        floor_templates_per_control=5, baseline_evaluation_limit=25,
        floor_template_rule='Q=(prefix+log2(C(8192,r0))+60)/r0+4; logx_a=scale*max(0,(Q+logC(16,a)-8a*alpha*log255)/floor_a)',
        all_empty_diagnostic_before_envelope_intercept=True,
        baseline_coordinate_rounds=0, maximum_refined_baselines=int(refine), maximum_refinement_rounds=int(refine),
        mark_alpha_fits=0, floor_marker_fits=0, directed=False, old_receipts_transferred=False,
        proposal_requires_fresh_directed_replay=True, whole_code_certificate=False,
        complete_requested_band=False, compute_budget_seconds=max_seconds,
        banks=[], witnesses=[], raw_regional_banks=0)
    banks = []
    try:
        for mark, alpha in CONTROLS:
            check_time()
            curves = regional_curves(fractional_regions(mark, alpha, LAST), LAST)
            weights = fractional_dense_weights(alpha)
            check_time()
            initial, _ = slope_candidates(curves)
            banks.append((curves, weights))
            result['raw_regional_banks'] += 1
            result['banks'].append(dict(mark=str(mark), alpha=str(alpha),
                curves=[[None if not np.isfinite(v) else float(v) for v in row] for row in curves],
                zero_entries_are_none=True, initial_log2_slopes=initial.tolist(),
                dense_log_weights={a: [None if not np.isfinite(v) else float(v) for v in row]
                                   for a, row in weights.items()}))
            for template in floor_templates(mark, alpha):
                check_time()
                row = evaluate(curves, weights, mark, alpha, template['log2_floor_markers'], initial)
                row.update(bank_index=len(banks)-1, coordinate_rounds=0, slope_evaluations=1,
                           floor_template=template['name'], floor_anchor=template['anchor'], floor_scale=template['scale'])
                result['witnesses'].append(row)
                if progress:
                    progress(dict(stage='baseline', completed=len(result['witnesses']),
                        mark=str(mark), alpha=str(alpha), floor_template=template['name'],
                        failure_log2_screen=row['failure_log2_screen'], elapsed_seconds=perf_counter()-start))
        check_time()
        baseline_index = min(range(len(result['witnesses'])),
                             key=lambda i: result['witnesses'][i]['failure_log2_screen'])
        result['best_baseline_index'] = baseline_index
        if refine:
            baseline = result['witnesses'][baseline_index]
            row = refine_one_round(*banks[baseline['bank_index']], baseline, check_time)
            row.update(bank_index=baseline['bank_index'], refined_baseline_index=baseline_index)
            result['witnesses'].append(row)
            if progress:
                progress(dict(stage='one-round-refinement', baseline_index=baseline_index,
                    failure_log2_screen=row['failure_log2_screen'], evaluations=row['slope_evaluations'],
                    elapsed_seconds=perf_counter()-start))
        check_time()
        entries = []
        for index, r in enumerate(range(FIRST, LAST+1)):
            witness = min(range(len(result['witnesses'])),
                          key=lambda i: result['witnesses'][i]['entries'][index]['failure_log2_screen'])
            entries.append(dict(result['witnesses'][witness]['entries'][index], witness_index=witness))
        result.update(status='completed', complete_requested_band=True, entries=entries,
            failure_log2_screen=float(logsum([row['failure_log2_screen'] for row in entries])),
            witness_selection_per_whole_weight=True, extra_witness_factor=False)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-start
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--no-refine', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen proposal artifacts must not be overwritten')
    result = run(max_seconds=args.max_seconds, refine=not args.no_refine,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('banks', 'witnesses', 'entries')}, allow_nan=False))
