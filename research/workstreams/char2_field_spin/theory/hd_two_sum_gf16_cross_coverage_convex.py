"""ONE global convex Chernoff marker, then matched cross-coverage controls.

The failed adjacent-coefficient driver stays pinned and unchanged. This
wrapper minimizes the full positive word polynomial on log2(X)>=0 using
bracketing and48 bisections. The resulting marker is shared by every path,
column, and both controls. No provider, regional fit, or witness grid.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_gf16_cross_coverage_control as old


DIRECTORY = Path(__file__).resolve().parent
OLD_SOURCE_SHA = 'f3b54eddb375d533119527c11f04b83ffaefb28b14c72c2955188cdc324f4d4a'
FAILED_ARTIFACT = 'hd_two_sum_gf16_cross_coverage_gate_float.json'
FAILED_ARTIFACT_SHA = '4286c359d475aa8d4253fac094d74eccaca65a904b5d2678ca16a07d08ea39dd'


def optimal_marker(word, floor, *, bisections=48, max_log_marker=512., check_time=lambda: None):
    """Minimize log2(sum a_d*2**(ell*d))-floor*ell on ell>=0.

    The derivative is the tilted mean degree minus floor. The objective
    is convex; its second derivative is ln(2) times the tilted variance.
    FLOAT only. An absent finite stationary bracket fails explicitly.
    """
    coefficients = np.asarray(word, dtype=float)
    if (coefficients.ndim != 1 or not len(coefficients) or np.any(np.isnan(coefficients))
            or np.any(np.isposinf(coefficients)) or not np.any(np.isfinite(coefficients))
            or type(floor) is not int or floor < 0
            or type(bisections) is not int or not 1 <= bisections <= 64
            or not np.isfinite(max_log_marker) or not 1 <= max_log_marker <= 512):
        raise ValueError('A positive finite log polynomial and bounded optimizer parameters required')
    degree = np.flatnonzero(np.isfinite(coefficients)).astype(float)
    coefficients = coefficients[degree.astype(int)]
    evaluations = 0
    def evaluate(marker):
        nonlocal evaluations
        check_time()
        values = coefficients+marker*(degree-floor)
        peak = float(np.max(values))
        scaled = np.exp2(values-peak)
        total = float(np.sum(scaled))
        derivative = float(np.dot(scaled, degree-floor)/total)
        evaluations += 1
        return peak+float(np.log2(total)), derivative
    at_zero, derivative = evaluate(0.)
    if derivative >= 0:
        return 0., dict(method='global-convex-mean-degree', endpoint_minimum=True,
            objective_log2=at_zero, mean_degree=derivative+floor, derivative=derivative,
            bracket=[0., 0.], bisections=0, evaluations=evaluations,
            max_log_marker=max_log_marker, finite_marker=True)
    if np.max(degree) <= floor:
        raise ValueError('No finite stationary marker: all supported degrees are at or below the floor')
    low, high = 0., 1.
    _, high_derivative = evaluate(high)
    while high_derivative < 0:
        low = high
        if high >= max_log_marker:
            raise ValueError('Stationary marker exceeds the explicit finite bracket cap')
        high = min(2*high, max_log_marker)
        _, high_derivative = evaluate(high)
    for _ in range(bisections):
        middle = (low+high)/2
        _, derivative = evaluate(middle)
        if derivative < 0:
            low = middle
        else:
            high = middle
    marker = (low+high)/2
    objective, derivative = evaluate(marker)
    if not np.isfinite(marker) or marker < 0 or not np.isfinite(np.exp2(marker)):
        raise ArithmeticError('The selected X must be finite and at least one')
    return marker, dict(method='global-convex-mean-degree', endpoint_minimum=False,
        objective_log2=objective, mean_degree=derivative+floor, derivative=derivative,
        bracket=[low, high], bisections=bisections, evaluations=evaluations,
        max_log_marker=max_log_marker, finite_marker=True)


def read_inputs():
    if sha256((DIRECTORY/'hd_two_sum_gf16_cross_coverage_control.py').read_bytes()).hexdigest() != OLD_SOURCE_SHA:
        raise ValueError('Frozen original coverage driver changed')
    return old.read_inputs()


def run(*, quiet=False, max_seconds=30, progress=None):
    if (type(quiet) is not bool or type(max_seconds) not in (int, float)
            or not np.isfinite(max_seconds) or not 0 < max_seconds <= 30):
        raise ValueError('Fixed quiet/full scope and explicit budget at most30 seconds required')
    reference, data, measure = read_inputs()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Convex-marker coverage control exhausted its budget')
    result = dict(active_groups=18, arity=8, source_full_pair_max=1,
        mark='691/700', alpha=1, scope='quiet-only' if quiet else 'complete homogeneous arity8/F1 class',
        source_count=data['source_count'], all_source_sets_in_stratum=True,
        all_regional_state_paths_included=True, fixed_rowset_holder_exponent=16,
        raw_weights_exclude_rs_normalizer=True, extra_rowset_factor=False,
        rs_normalizer_exponent=1152, projective_normalizations=1,
        source_count_once=True, whole_word_surrogate_floor=1296,
        surrogate='min(output_weight,8)', individual_row_floors=False,
        initial_state=0, terminal_vector=[1, 1], matrix_order='chronological',
        coverage_by_positive_dp=True, signed_inclusion_exclusion=False,
        marker_source_scope='full-outer full-state ordered word polynomial',
        fixed_marker_count=1, per_path_markers=False, new_moment_candidates=0,
        raw_regional_banks_built=0, regional_parameter_fits=0,
        marker_selection='ONE global convex scalar minimization;48 bisections',
        requested_controls=2, prefix_log2=reference['prefix_log2'], state_measure=measure,
        frozen_original_source_sha256=OLD_SOURCE_SHA,
        reference_artifact=old.BASE_ARTIFACT, reference_artifact_sha256=old.BASE_ARTIFACT_SHA,
        preserved_failed_marker_artifact=FAILED_ARTIFACT,
        preserved_failed_marker_artifact_sha256=FAILED_ARTIFACT_SHA,
        compute_budget_seconds=max_seconds, directed=False,
        whole_code_certificate=False, complete_requested_controls=False, controls=[])
    try:
        full_weights = old.base.raw_holder_weights()
        matrix, _ = old.base.column_matrix(full_weights, old.base.GROUPS, measure,
            arity=old.base.ARITY, check_time=check_time)
        word = old.base.matrix_poly_power(matrix, old.base.COLUMNS, check_time)
        initial_row = np.logaddexp2.reduce(word[0], axis=0)
        polynomial = reference['prefix_log2']+float(np.logaddexp2.reduce(initial_row[old.base.FLOOR:]))
        expected = next(row['log2_screen'] for row in reference['controls'] if row['mode'] == 'polynomial')
        if abs(polynomial-expected) > 1e-7:
            raise ArithmeticError('Frozen full matrix-polynomial tail did not reproduce')
        marker, diagnostics = optimal_marker(initial_row, old.base.FLOOR, check_time=check_time)
        result.update(log2_marker=marker, marker_log_as_exact_float=str(F(marker)),
            scalar_marker_as_exact_float=str(F(float(np.exp2(marker)))),
            optimizer=diagnostics, reproduced_full_polynomial_log2=polynomial)
        if progress:
            progress(dict(stage='global-convex-marker', log2_marker=marker,
                expected_scalar_log2=reference['prefix_log2']+diagnostics['objective_log2'],
                mean_degree=diagnostics['mean_degree']))
        weights = old.base.raw_holder_weights(quiet=True) if quiet else full_weights
        covered, unrestricted, costs = old.paired_columns(weights, old.base.GROUPS, measure,
            arity=old.base.ARITY, log_marker=marker, check_time=check_time)
        original, _ = old.base.column_matrix(weights, old.base.GROUPS, measure,
            arity=old.base.ARITY, log_marker=marker, check_time=check_time)
        if not np.allclose(unrestricted, original, rtol=0., atol=1e-8):
            raise ArithmeticError('Matched unrestricted scalar matrix did not reproduce')
        for label, scalar_matrix in (('unrestricted', unrestricted), ('cross-covered', covered)):
            check_time()
            power = old.base.log_matrix_power(scalar_matrix, old.base.COLUMNS)
            scalar = float(np.logaddexp2.reduce(power[0]))-old.base.FLOOR*marker
            if label == 'unrestricted' and not quiet and abs(scalar-diagnostics['objective_log2']) > 1e-7:
                raise ArithmeticError('Scalar matrix evaluation disagrees with full polynomial evaluation')
            row = dict(control=label, offdiagonal_empty_zeros_retained=label == 'cross-covered',
                column_log2_matrix=scalar_matrix.tolist(), tail_upper_log2=scalar,
                log2_screen=reference['prefix_log2']+scalar, **costs)
            result['controls'].append(row)
            if progress:
                progress(row)
        result['coverage_gain_bits'] = result['controls'][0]['log2_screen']-result['controls'][1]['log2_screen']
        check_time()
        result['complete_requested_controls'] = True
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    answer = run(quiet=args.quiet, max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=answer['complete_requested_controls'],
        output_sha256=sha256(args.output.read_bytes()).hexdigest())))
