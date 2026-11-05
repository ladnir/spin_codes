"""Three bounded FLOAT controls for the G16 low-total-arity stratum.

At each fixed active-group weight r, only A_tot<=B=A0-1 is claimed. A
deterministic0<t<=1 contributes t**(a/16) to type-a column weights and
t**(-B) once outside alpha. Types, markers and floor witnesses stay fixed
across every chronological path and all16 columns. No high-side result or
continuous active-group band is inferred by this helper.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_gf16_fractional_state as fractional
import hd_two_sum_gf16_state_envelope as state


ANCHORS = ((18, 247, F(691, 700)), (24, 317, F(172, 175)), (32, 463, F(171, 175)))
ALPHA = F(1, 5)
LOG_T = (0, -1, -4, -16, -64)
MAXIMUM = 32
FRACTIONAL_SHA256 = 'bbeca202d77780ea325bbc70e959e73e929beec7c0341c2d80d85023c96273d8'
STATE_SHA256 = '5107b945ab18415191c48a3434a354c101d8fc5aeda66ec83a3fd2a5c3498af7'


def checked_dependencies():
    pins = {}
    for module, expected in ((fractional, FRACTIONAL_SHA256), (state, STATE_SHA256)):
        path = Path(module.__file__)
        actual = sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Reviewed dependency changed: {path.name}')
        pins[path.name] = actual
    return pins


def checked_low_scope(r, a0, mark, alpha=ALPHA):
    if ((r, a0, F(mark)) not in ANCHORS or F(alpha) != ALPHA
            or type(r) is not int or type(a0) is not int):
        raise ValueError('Only the three frozen low-arity anchors are in scope')
    return a0-1


def floor_templates(r, mark, alpha=ALPHA):
    """Zero plus scale1/2 at THIS active-group anchor; no type adaptation."""
    if r not in (18, 24, 32) or F(alpha) != ALPHA:
        raise ValueError('Frozen anchor and exponent required')
    prefix = float(alpha)*(-fractional.LIMIT*log2(float(mark))-log2(255))
    empty = np.array([log2(comb(16, a))-8*a*float(alpha)*log2(255) for a in range(1, 17)])
    counts = np.array([9, 18]+[9]*14)
    q = (prefix+log2(comb(8192, r))+60)/r+4
    result = [dict(name='all-zero', scale=0, log2_floor_markers=[F(0)]*16)]
    for scale in (1, 2):
        values = scale*np.maximum(0., (q+empty)/counts)
        result.append(dict(name=f'anchor{r}-scale{scale}', scale=scale,
                           log2_floor_markers=[F(float(value)) for value in values]))
    return result


def peak_tail_slopes(curves):
    """Initial peak-to-tail secants; candidates never use the empty endpoint.

    Zero and all secants among first nonempty, peak, midpoint and tail are
    included (at most7 candidates per entry). Intercepts still maximize
    over EVERY occupancy, including zero and every unselected position.
    """
    curves = np.asarray(curves)
    if (curves.ndim != 2 or curves.shape[0] != 4 or curves.shape[1] < 2
            or np.any(np.isnan(curves)) or np.any(np.isposinf(curves))):
        raise ValueError('Complete finite-or-zero four-entry curves required')
    initial, candidates, points = [], [], []
    for curve in curves:
        finite = np.flatnonzero(np.isfinite(curve[1:]))+1
        if not len(finite):
            raise ValueError('Every entry must have a nonempty-support value')
        peak = int(finite[np.argmax(curve[finite])])
        tail = int(finite[-1])
        middle = int(finite[len(finite)//2])
        selected = sorted({int(finite[0]), peak, middle, tail})
        choices = {0.}
        for i, left in enumerate(selected):
            for right in selected[i+1:]:
                choices.add(float((curve[right]-curve[left])/(right-left)))
        initial.append(0. if peak == tail else float((curve[tail]-curve[peak])/(tail-peak)))
        candidates.append(tuple(sorted(choices)))
        points.append(dict(nonempty_peak=peak, tail=tail, secant_positions=selected))
    return np.asarray(initial), tuple(candidates), points


def tilt_type_logs(type_logs, log_t, columns=16):
    """Whole-word factor t**a, distributed over columns; NOT t**(alpha*a)."""
    values = np.asarray(type_logs)
    if (values.ndim != 2 or values.shape[1] != 16 or type(columns) is not int or columns < 1
            or not np.isfinite(log_t) or log_t > 0
            or np.any(np.isnan(values)) or np.any(np.isposinf(values))):
        raise ValueError('Complete type logs and a nonpositive arity-marker log required')
    return values+np.arange(1, 17)*(float(log_t)/columns)


def prepared_envelope(curves, weights, floors, slopes):
    intercepts, slack = state.affine_envelope(curves, slopes)
    classes = state.path_classes(16)
    counts = np.asarray([row[2] for row in classes])
    marks = np.asarray([np.repeat(slopes, row[2]) for row in classes])
    typed = fractional.fractional_type_weights(marks, weights, floors, ALPHA)
    return dict(classes=classes, path_a=counts@intercepts, typed=typed,
                intercepts=intercepts, slack=slack, slopes=np.asarray(slopes), floors=list(floors))


def evaluate_prepared(prepared, r, a0, mark, log_t):
    bound = checked_low_scope(r, a0, mark)
    if log_t not in LOG_T:
        raise ValueError('Only the five frozen total-arity marker logs are permitted')
    tilted = tilt_type_logs(prepared['typed'], log_t)
    kernel = state.kernel_from_types(prepared['classes'], prepared['path_a'], tilted, r)
    power = state.matrix_power_log(kernel, 16)
    prefix = float(ALPHA)*(-fractional.LIMIT*log2(float(mark))-log2(255))
    probability_log = (prefix+log2(comb(8192, r))-bound*log_t
                       +float(state.logsum(power[0])))
    empty = np.array([log2(comb(16, a))-8*a*float(ALPHA)*log2(255) for a in range(1, 17)])
    floors = fractional.floor_array(prepared['floors'])
    empty -= np.array([9, 18]+[9]*14)*floors
    empty += np.arange(1, 17)*log_t
    empty_bound = prefix+log2(comb(8192, r))-bound*log_t+r*float(state.logsum(empty))
    return dict(active_groups=r, high_arity_cutoff=a0, low_total_active_constituents_max=bound,
        mark=str(mark), alpha=str(ALPHA), log2_total_arity_marker=log_t,
        total_arity_marker=str(F(1, 2**(-log_t))),
        outside_total_arity_log_factor=-bound*log_t,
        log2_floor_markers=[str(F(value)) for value in prepared['floors']],
        log2_slopes=prepared['slopes'].tolist(), log2_intercepts=prepared['intercepts'].tolist(),
        envelope_contacts=[np.flatnonzero(row < 1e-8).tolist() for row in prepared['slack']],
        envelope_minimum_slack=float(np.min(prepared['slack'])),
        zero_occupancy_envelope_log2=prepared['intercepts'].tolist(),
        all_empty_mask_log2_screen=float(empty_bound), failure_log2_screen=float(probability_log))


def one_joint_round(curves, weights, anchor, baseline, templates, candidates, check_time):
    """One coordinate pass: logt, floor template, then four entry slopes."""
    r, a0, mark = anchor
    best = baseline
    memo = {}
    def objective(slopes, template_index, log_t):
        check_time()
        key = (tuple(float(x) for x in slopes), template_index, log_t)
        if key not in memo:
            prepared = prepared_envelope(curves, weights, templates[template_index]['log2_floor_markers'], slopes)
            row = evaluate_prepared(prepared, r, a0, mark, log_t)
            row.update(floor_template_index=template_index, floor_template=templates[template_index]['name'])
            memo[key] = row
        return memo[key]
    for log_t in LOG_T:
        row = objective(best['log2_slopes'], best['floor_template_index'], log_t)
        if row['failure_log2_screen'] < best['failure_log2_screen']:
            best = row
    for template_index in range(3):
        row = objective(best['log2_slopes'], template_index, best['log2_total_arity_marker'])
        if row['failure_log2_screen'] < best['failure_log2_screen']:
            best = row
    for coordinate in range(4):
        slopes = np.asarray(best['log2_slopes'])
        for slope in candidates[coordinate]:
            trial = slopes.copy()
            trial[coordinate] = slope
            row = objective(trial, best['floor_template_index'], best['log2_total_arity_marker'])
            if row['failure_log2_screen'] < best['failure_log2_screen']:
                best = row
    return dict(best, joint_coordinate_rounds=1, joint_evaluations=len(memo),
        joint_trials=[dict(log2_slopes=row['log2_slopes'], floor_template_index=row['floor_template_index'],
            log2_total_arity_marker=row['log2_total_arity_marker'], failure_log2_screen=row['failure_log2_screen'])
            for row in memo.values()])


def run(*, max_seconds=30, refine_worst=False, progress=None):
    if (type(max_seconds) not in (int, float) or not np.isfinite(max_seconds)
            or not 0 < max_seconds <= 30 or type(refine_worst) is not bool):
        raise ValueError('Explicit at-most30s budget and Boolean optional refinement required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start >= max_seconds:
            raise TimeoutError('Low-total-arity anchor control exhausted its total budget')
    result = dict(construction=fractional.CONSTRUCTION, field_order=256, forbidden_weight=fractional.LIMIT,
        anchors=[dict(active_groups=r, high_arity_cutoff=a0, low_total_active_constituents_max=a0-1,
                      mark=str(z), alpha=str(ALPHA)) for r, a0, z in ANCHORS],
        source_counts='C(8192,r)', all_source_supports=True, all_cohorts_within_low_total_arity=True,
        all_constituent_row_sets_within_low_total_arity=True, continuous_weight_band_covered=False,
        fractional_mask_mass=True, packet_power_before_geometry_average=True, route_normalizer_powered=False,
        total_arity_marker_outside_alpha=True, total_arity_marker_per_type_column='t^(a/16)',
        fixed_types_across_columns=True, type_markers_fixed_across_paths_and_columns=True,
        independently_resampled_rowsets=False, maximum_regional_occupancy=MAXIMUM,
        envelope_scope='all feasible matching histograms and0<=m<=anchor_r',
        initial_slopes='nonempty peak-to-tail; zero if no nontrivial tail',
        fixed_log2_total_arity_markers=list(LOG_T), floor_templates_per_anchor=3,
        baseline_evaluation_limit=45, maximum_extra_refined_anchors=int(refine_worst),
        maximum_extra_joint_coordinate_rounds=int(refine_worst), fixed_z_alpha=True,
        regional_banks_rebuilt=0, directed=False, whole_code_certificate=False,
        requires_fresh_directed_replay=True, high_side_included=False,
        complete_requested_anchors=False, compute_budget_seconds=max_seconds, banks=[], witnesses=[])
    banks, template_sets, candidate_sets, best_indices = [], [], [], []
    try:
        check_time()
        result['dependency_pins'] = checked_dependencies()
        weights = fractional.fractional_dense_weights(ALPHA)
        for anchor_index, (r, a0, mark) in enumerate(ANCHORS):
            check_time()
            all_curves = state.regional_curves(fractional.fractional_regions(mark, ALPHA, MAXIMUM), MAXIMUM)
            curves = all_curves[:, :r+1]
            initial, candidates, points = peak_tail_slopes(curves)
            templates = floor_templates(r, mark)
            banks.append(curves)
            template_sets.append(templates)
            candidate_sets.append(candidates)
            result['regional_banks_rebuilt'] += 1
            result['banks'].append(dict(active_groups=r, mark=str(mark), alpha=str(ALPHA),
                curves=[[None if not np.isfinite(x) else float(x) for x in row] for row in all_curves],
                zero_entries_are_none=True, envelope_maximum=r, initial_log2_slopes=initial.tolist(),
                slope_candidates=[list(row) for row in candidates], initial_slope_points=points))
            indices = []
            for template_index, template in enumerate(templates):
                prepared = prepared_envelope(curves, weights, template['log2_floor_markers'], initial)
                for log_t in LOG_T:
                    check_time()
                    row = evaluate_prepared(prepared, r, a0, mark, log_t)
                    row.update(anchor_index=anchor_index, floor_template_index=template_index,
                               floor_template=template['name'], joint_coordinate_rounds=0)
                    indices.append(len(result['witnesses']))
                    result['witnesses'].append(row)
            best_index = min(indices, key=lambda i: result['witnesses'][i]['failure_log2_screen'])
            best_indices.append(best_index)
            if progress:
                row = result['witnesses'][best_index]
                progress(dict(stage='anchor-baselines', active_groups=r, completed_anchors=anchor_index+1,
                    failure_log2_screen=row['failure_log2_screen'], log2_total_arity_marker=row['log2_total_arity_marker'],
                    floor_template=row['floor_template'], elapsed_seconds=perf_counter()-start))
        check_time()
        result['initial_best_witness_indices'] = list(best_indices)
        if refine_worst:
            worst = max(range(3), key=lambda i: result['witnesses'][best_indices[i]]['failure_log2_screen'])
            row = one_joint_round(banks[worst], weights, ANCHORS[worst],
                result['witnesses'][best_indices[worst]], template_sets[worst], candidate_sets[worst], check_time)
            row.update(anchor_index=worst, refined_baseline_index=best_indices[worst])
            result['refined_anchor_index'] = worst
            best_indices[worst] = len(result['witnesses'])
            result['witnesses'].append(row)
        check_time()
        rows = [dict(result['witnesses'][index], witness_index=index) for index in best_indices]
        result.update(status='completed', complete_requested_anchors=True,
            selected_anchor_bounds=rows, selected_witness_indices=best_indices,
            failure_log2_screen=float(state.logsum([row['failure_log2_screen'] for row in rows])),
            deterministic_witness_selection_per_anchor=True, extra_witness_factor=False)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-start
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--refine-worst', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen artifacts must not be overwritten')
    result = run(max_seconds=args.max_seconds, refine_worst=args.refine_worst,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('banks', 'witnesses', 'selected_anchor_bounds')}, allow_nan=False))
