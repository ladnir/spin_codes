"""FLOAT all-cohort control retaining the two regional states.

For each fixed type vector, row-set Holder bounds give one column matrix.
Holder over the SAME type vector across16 columns gives an entrywise
ell16 norm. Minkowski then sums the16-region state paths. Type floor
markers remain fixed across every state path and column.

No common potential is used. Positive exponential regional envelopes can
overestimate exact zero entries; their zero-occupancy slack is reported.
This isolated diagnostic does not change any provider or certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_factorized import forced_nonzero_regions, dense_arity_logs
from hd_two_sum_gf16_single_dense_matrix import log_matmul
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_xor_multigroup import LIMIT


MARKS = (F(99, 100), F(49, 50), F(24, 25))
FLOOR_LOGS = (0, 1, 4, 8)
FIRST, LAST, COLUMNS, REGIONS = 18, 32, 16, 16


def logsum(values, axis=None):
    return np.logaddexp2.reduce(np.asarray(values), axis=axis)


@lru_cache(None)
def path_classes(length=16):
    """Exact counts of binary paths by endpoints and four transition counts."""
    if type(length) is not int or not 1 <= length <= 16:
        raise ValueError('Path length1..16 required')
    values = {(i, i, (0, 0, 0, 0)): 1 for i in range(2)}
    for _ in range(length):
        following = {}
        for (initial, state, counts), multiplicity in values.items():
            for final in range(2):
                target = list(counts)
                target[2*state+final] += 1
                key = (initial, final, tuple(target))
                following[key] = following.get(key, 0)+multiplicity
        values = following
    return tuple((*key, multiplicity) for key, multiplicity in sorted(values.items()))


def regional_curves(regions, maximum):
    """Four log-entry maxima over all matching counts, with exact zeros kept."""
    required = {(m, f) for m in range(maximum+1) for f in range(m//2+1)}
    if (type(maximum) is not int or maximum < 1 or set(regions) != required):
        raise ValueError('Complete regional matching histograms required')
    values = []
    for m in range(maximum+1):
        matrices = np.asarray([regions[m, f] for f in range(m//2+1)])
        if (matrices.shape != (m//2+1, 2, 2) or np.any(np.isnan(matrices))
                or np.any(np.isposinf(matrices))):
            raise ValueError('Finite-or-zero log regional matrices required')
        values.append(np.max(matrices, axis=0).ravel())
    result = np.asarray(values).T
    if np.any(~np.any(np.isfinite(result), axis=1)):
        raise ValueError('Each entry must be nonzero somewhere in the requested range')
    return result


def affine_envelope(curves, slopes):
    """Smallest intercept for each supplied slope, checked at EVERY occupancy."""
    curves, slopes = np.asarray(curves), np.asarray(slopes)
    if (curves.ndim != 2 or curves.shape[0] != 4 or slopes.shape != (4,)
            or np.any(~np.isfinite(slopes)) or np.any(np.isnan(curves))
            or np.any(np.isposinf(curves))):
        raise ValueError('Four complete log curves and four finite slopes required')
    intercepts = np.max(curves-slopes[:, None]*np.arange(curves.shape[1]), axis=1)
    if np.any(~np.isfinite(intercepts)):
        raise ValueError('Each regional entry needs a positive majorant')
    slack = intercepts[:, None]+slopes[:, None]*np.arange(curves.shape[1])-curves
    if np.min(slack) < -1e-10:
        raise ArithmeticError('Complete affine envelope check failed')
    return intercepts, slack


def slope_candidates(curves):
    """At most16 deterministic slopes/entry, not a continuous optimizer.

    Use all secants among selected finite occupancies, plus slope0. The
    first/last finite secant initializes each coordinate. Intercepts still
    use every occupancy, including those not selected for secants.
    """
    maximum = curves.shape[1]-1
    candidates, initial = [], []
    for curve in curves:
        finite = np.flatnonzero(np.isfinite(curve))
        if not len(finite):
            raise ValueError('A nonzero regional entry is required')
        selected = sorted(set([int(finite[0]), int(finite[-1])]
            +[m for m in (0, 1, 2, maximum//4, maximum//2, maximum)
              if np.isfinite(curve[m])]))
        # Bound the selected set to six points, retaining both endpoints.
        if len(selected) > 6:
            selected = [selected[0], *selected[1:5], selected[-1]]
        slopes = {0.}
        for i, left in enumerate(selected):
            for right in selected[i+1:]:
                slopes.add(float((curve[right]-curve[left])/(right-left)))
        candidates.append(tuple(sorted(slopes)))
        initial.append(0. if finite[0] == finite[-1] else
            float((curve[finite[-1]]-curve[finite[0]])/(finite[-1]-finite[0])))
    return np.asarray(initial), tuple(candidates)


def elementary_logs(log_marks):
    """Log elementary symmetric polynomials, including each labelled bit once."""
    log_marks = np.asarray(log_marks)
    if log_marks.ndim < 1 or np.any(~np.isfinite(log_marks)):
        raise ValueError('Finite positive marker logs required')
    length = log_marks.shape[-1]
    result = np.full(log_marks.shape[:-1]+(length+1,), -np.inf)
    result[..., 0] = 0.
    for index in range(length):
        result[..., 1:index+2] = np.logaddexp2(result[..., 1:index+2],
            result[..., :index+1]+log_marks[..., index, None])
    return result


def type_weights(log_marks, dense_weights, log_floor, *, columns=16,
                 dimension=8, minimum=9, q=256):
    """Per-column logs for types1,2,3..g, with fixed whole-word floor charges.

    Type1/2 row choices enter through columns-th roots. Dense weights
    already contain their fixed-row-set sum and whole-word normalization.
    All nonempty dense masks receive the SAME marker before path summation.
    """
    log_marks = np.asarray(log_marks)
    group = log_marks.shape[-1]
    if (not 2 <= group < q-1 or set(dense_weights) != set(range(3, group+1))
            or type(columns) is not int or not 1 <= minimum <= columns
            or not 1 <= dimension <= columns or not np.isfinite(log_floor) or log_floor < 0):
        raise ValueError('Complete type weights and one fixed nonnegative floor log required')
    e = elementary_logs(log_marks)
    full, logq = e[..., group], log2(q-1)
    single = (log2(group)-dimension*logq-minimum*log_floor)/columns
    single += np.logaddexp2(0., logq+log_floor+full)
    pair_terms = np.stack((np.zeros_like(full), 1+logq+log_floor+full,
        logq+2*log_floor+np.logaddexp2(log2(q-1-group)+full, e[..., group-1])), axis=-1)
    pair = (log2(comb(group, 2))-2*dimension*logq-2*minimum*log_floor)/columns
    pair += logsum(pair_terms, axis=-1)
    values = [single, pair]
    for arity in range(3, group+1):
        weights = np.asarray(dense_weights[arity])
        if (weights.shape != (group+1,) or np.any(np.isnan(weights)) or np.any(np.isposinf(weights))):
            raise ValueError('Complete finite-or-zero dense cardinality weights required')
        value = np.logaddexp2(weights[0], log_floor+logsum(e[..., 1:]+weights[1:], axis=-1))
        values.append(value-minimum*log_floor/columns)
    return np.stack(values, axis=-1)


def kernel_from_types(classes, log_path_a, type_logs, active_groups, columns=16):
    """Minkowski over paths, then factor the fixed ordered group-type vector."""
    log_path_a, type_logs = np.asarray(log_path_a), np.asarray(type_logs)
    if (log_path_a.shape != (len(classes),) or type_logs.ndim != 2
            or type_logs.shape[0] != len(classes) or type(active_groups) is not int
            or active_groups < 1 or type(columns) is not int or columns < 1
            or np.any(~np.isfinite(log_path_a)) or np.any(np.isnan(type_logs))
            or np.any(np.isposinf(type_logs))):
        raise ValueError('Complete path/type table and positive group/column counts required')
    terms = (log_path_a+active_groups/columns*logsum(columns*type_logs, axis=-1)
             +np.log2([row[3] for row in classes]))
    matrix = np.full((2, 2), -np.inf)
    for i in range(2):
        for j in range(2):
            matrix[i, j] = logsum(terms[[row[0] == i and row[1] == j for row in classes]])
    return matrix


def matrix_power_log(matrix, exponent):
    if type(exponent) is not int or exponent < 1:
        raise ValueError('Positive matrix power required')
    result = np.array([[0., -np.inf], [-np.inf, 0.]])
    base = np.asarray(matrix)
    while exponent:
        if exponent & 1:
            result = log_matmul(result, base)
        exponent >>= 1
        if exponent:
            base = log_matmul(base, base)
    return result


def evaluate(curves, dense_weights, mark, log_floor, slopes):
    intercepts, slack = affine_envelope(curves, slopes)
    classes = path_classes(REGIONS)
    counts = np.asarray([row[2] for row in classes])
    marks = np.asarray([np.repeat(slopes, row[2]) for row in classes])
    typed = type_weights(marks, dense_weights, log_floor)
    path_a = counts @ intercepts
    prefix = -LIMIT*log2(float(mark))-log2(255)
    entries = []
    for r in range(FIRST, LAST+1):
        kernel = kernel_from_types(classes, path_a, typed, r)
        power = matrix_power_log(kernel, COLUMNS)
        entries.append(dict(active_groups=r, failure_log2_screen=float(
            prefix+log2(comb(8192, r))+logsum(power[0]))))
    return dict(mark=str(mark), log2_floor_marker=int(log_floor),
        log2_intercepts=[float(x) for x in intercepts], log2_slopes=[float(x) for x in slopes],
        envelope_contacts=[np.flatnonzero(row < 1e-8).tolist() for row in slack],
        envelope_minimum_slack=float(np.min(slack)),
        zero_occupancy_envelope_log2=[float(x) for x in intercepts],
        entries=entries, failure_log2_screen=float(logsum([row['failure_log2_screen'] for row in entries])))


def selected_search(curves, weights, mark, log_floor, check_time, rounds=3):
    """Bounded deterministic coordinate selection from at most16 secants/entry."""
    if type(rounds) is not int or not 0 <= rounds <= 3:
        raise ValueError('At most three deterministic coordinate rounds')
    slopes, candidates = slope_candidates(curves)
    memo = {}
    def objective(values):
        check_time()
        key = tuple(float(x) for x in values)
        if key not in memo:
            memo[key] = evaluate(curves, weights, mark, log_floor, np.asarray(key))
        return memo[key]
    best = objective(slopes)
    for _ in range(rounds):
        changed = False
        for index in range(4):
            previous = slopes[index]
            for candidate in candidates[index]:
                trial = slopes.copy()
                trial[index] = candidate
                row = objective(trial)
                if row['failure_log2_screen'] < best['failure_log2_screen']:
                    best = row
            slopes = np.asarray(best['log2_slopes'])
            changed |= bool(slopes[index] != previous)
        if not changed:
            break
    return dict(best, evaluations=len(memo), coordinate_rounds_limit=rounds,
        slope_candidates=[list(row) for row in candidates],
        exact_float_slopes=[str(F(x)) for x in best['log2_slopes']])


def run(*, curves_only=False, rounds=3, max_seconds=30, progress=None):
    if (type(curves_only) is not bool or type(rounds) is not int or not 0 <= rounds <= 3
            or type(max_seconds) not in (int, float) or not np.isfinite(max_seconds)
            or not 0 < max_seconds <= 30):
        raise ValueError('Explicit at-most30s budget and at-most3 rounds required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start >= max_seconds:
            raise TimeoutError('Bounded state-envelope control exhausted its budget')
    result = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        forbidden_weight=LIMIT, active_group_band=[FIRST, LAST], all_cohort_splits=True,
        all_constituent_row_sets=True, all_source_supports=True, source_counts='C(8192,r)',
        total_columns=COLUMNS, regions_per_column=REGIONS, regional_alpha=1,
        regional_amplitudes='independent-uniform-nonzero', pair_swaps=True,
        initial_state=0, terminal_vector=[1, 1], common_potential_used=False,
        endpoint_potential_factors=0, projective_normalizations=1,
        fixed_types_across_columns=True, type_holder_exponent=16,
        path_minkowski_after_type_holder=True, type_markers_fixed_across_paths_and_columns=True,
        singleton_row_floor=9, pair_aggregate_row_floor=18, dense_nonempty_column_floor=9,
        rowset_multiplicity_inside_type_weights=True, extra_cohort_factor=False,
        exact_binary_path_classes=len(path_classes()), fixed_marks=list(map(str, MARKS)),
        fixed_log2_floor_markers=list(FLOOR_LOGS), mark_fits=0, floor_marker_fits=0,
        directed=False, proposal_requires_fresh_directed_replay=True, whole_code_certificate=False,
        complete_requested_band=False, compute_budget_seconds=max_seconds,
        curves_only=curves_only, raw_regional_banks=0, banks=[], witnesses=[])
    try:
        weights, pins = dense_arity_logs()
        result['exact_dense_moment_sha256'] = pins
        for mark in MARKS:
            check_time()
            curves = regional_curves(forced_nonzero_regions(float(mark), LAST), LAST)
            check_time()
            result['raw_regional_banks'] += 1
            initial, candidates = slope_candidates(curves)
            result['banks'].append(dict(mark=str(mark),
                curves=[[None if not np.isfinite(v) else float(v) for v in row] for row in curves],
                zero_entries_are_none=True, initial_log2_slopes=initial.tolist(),
                slope_candidates=[list(row) for row in candidates]))
            if not curves_only:
                for log_floor in FLOOR_LOGS:
                    row = selected_search(curves, weights, mark, log_floor, check_time, rounds)
                    result['witnesses'].append(row)
                    if progress:
                        progress(dict(completed_witnesses=len(result['witnesses']), mark=str(mark),
                            log2_floor_marker=log_floor, failure_log2_screen=row['failure_log2_screen'],
                            evaluations=row['evaluations'], envelope_contacts=row['envelope_contacts'],
                            elapsed_seconds=perf_counter()-start))
        check_time()
        result['status'] = 'curves-complete' if curves_only else 'completed'
        if not curves_only:
            selected = []
            for index, r in enumerate(range(FIRST, LAST+1)):
                witness = min(range(len(result['witnesses'])), key=lambda w:
                    result['witnesses'][w]['entries'][index]['failure_log2_screen'])
                selected.append(dict(result['witnesses'][witness]['entries'][index], witness_index=witness))
            result.update(complete_requested_band=True, entries=selected,
                failure_log2_screen=float(logsum([row['failure_log2_screen'] for row in selected])),
                witness_selection_per_weight=True, extra_witness_factor=False)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-start
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--curves-only', action='store_true')
    parser.add_argument('--rounds', type=int, default=3)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite an isolated FLOAT control')
    data = run(curves_only=args.curves_only, rounds=args.rounds, max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    data['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(data, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({key: value for key, value in data.items()
                     if key not in ('witnesses', 'banks', 'entries', 'exact_dense_moment_sha256')}, allow_nan=False))
