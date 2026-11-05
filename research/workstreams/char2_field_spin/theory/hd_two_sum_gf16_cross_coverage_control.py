"""One saved full-word marker; positive cross coverage versus matched control.

No provider, slope change, parameter grid, or per-path marker choice.
Recompute the frozen full matrix polynomial once to read floor/floor+1.
Both scalar controls then use that one fixed marker. Optional quiet scope
uses the SAME full-word marker. Each completed stage survives timeout.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, expm1, log, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_matrix_moment as base
from hd_two_sum_gf16_cross_coverage import (covered_groups, subset_weights,
    unrestricted_groups)


DIRECTORY = Path(__file__).resolve().parent
BASE_SOURCE_SHA = 'bc192ebdafc9abec59f7977e660c5d4741737eb62cf479c920fdbb83a060ce67'
COVER_SOURCE_SHA = 'e161d506e040abcf8446d931342a9649b892ec62530701db26713b85204f82a6'
BASE_ARTIFACT = 'hd_two_sum_matrix_moment_gate_float.json'
BASE_ARTIFACT_SHA = '0465d355fd56c5c0cdfabac206a377ed74f62bf638b78af354bb4a147a9b5783'


def read_inputs():
    for name, expected in (('hd_two_sum_matrix_moment.py', BASE_SOURCE_SHA),
                           ('hd_two_sum_gf16_cross_coverage.py', COVER_SOURCE_SHA)):
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: '+name)
    raw = (DIRECTORY/BASE_ARTIFACT).read_bytes()
    if sha256(raw).hexdigest() != BASE_ARTIFACT_SHA:
        raise ValueError('Frozen full ordered-matrix control changed')
    reference = json.loads(raw)
    required = dict(active_groups=18, arity=8, source_full_pair_max=1,
        mark='691/700', alpha='1', mode='both', complete_requested_controls=True,
        extra_rowset_factor=False, raw_weights_exclude_rs_normalizer=True,
        source_sha256=BASE_SOURCE_SHA)
    if any(reference.get(key) != value for key, value in required.items()):
        raise ValueError('Exact full homogeneous a8/F1 reference required')
    data, curves = base.load_saved()
    if not np.all(np.isneginf(curves[1:3, 0])):
        raise ArithmeticError('Both offdiagonal occupancy-zero entries must be exactly zero')
    if data['source_count'] != reference['source_count']:
        raise ValueError('Exact source count mismatch')
    old = data['two_atom_envelope']
    measure = base.entry_measure(curves, old['log2_probability'], old['log2_nonzero_atom'])
    if measure['exact_float_parameters'] != reference['state_measure']['exact_float_parameters']:
        raise ArithmeticError('Frozen four-atom parameters did not reproduce')
    return reference, data, measure


def paired_columns(weights, groups, measure, *, arity, log_marker,
                   check_time=lambda: None):
    """Return covered/unrestricted matrices under identical path weights."""
    weights = np.asarray(weights, dtype=float)
    size = len(weights)-1
    if (weights.ndim != 1 or not 1 <= size <= 16 or type(groups) is not int or groups < 1
            or type(arity) is not int or not 1 <= arity <= size
            or not np.isfinite(log_marker) or log_marker < 0
            or np.any(np.isnan(weights)) or np.any(np.isposinf(weights))):
        raise ValueError('Valid fixed group/region parameters and scalar X>=1 required')
    logp = measure['log2_probability']
    if not np.isfinite(logp) or logp > 0:
        raise ValueError('The00 atom probability must belong to(0,1]')
    log_other = -np.inf if logp == 0 else log2(-expm1(logp*log(2)))
    covered = np.full((2, 2), -np.inf)
    unrestricted = np.full((2, 2), -np.inf)
    cached = {}
    terms = 0
    for initial, final, counts, multiplicity in base.path_classes(size):
        check_time()
        n00, n01, n10, n11 = counts
        state_mass = (log2(multiplicity)+n01*measure['log2_scale01']
                      +n10*measure['log2_scale10']+n11*measure['log2_scale11'])
        for active00 in range(n00+1):
            if logp == 0 and active00 != n00:
                continue
            atom_mass = (log2(comb(n00, active00))+active00*logp
                         +((n00-active00)*log_other if active00 < n00 else 0.))
            key = active00, n01+n10, n11
            if key not in cached:
                check_time()
                coefficients = subset_weights(weights, *key,
                    measure['log2_atom00'], measure['log2_atom_cross'],
                    log_marker=log_marker, arity=arity)
                all_masks = unrestricted_groups(coefficients, groups)
                complete = covered_groups(coefficients, groups)
                if complete > all_masks+1e-9:
                    raise ArithmeticError('Coverage exceeded the matched positive sum')
                cached[key] = complete, all_masks
            complete, all_masks = cached[key]
            factor = state_mass+atom_mass
            covered[initial, final] = np.logaddexp2(covered[initial, final], complete+factor)
            unrestricted[initial, final] = np.logaddexp2(unrestricted[initial, final], all_masks+factor)
            terms += 1
    return covered, unrestricted, dict(path_classes=len(base.path_classes(size)),
        positive_path_atom_terms=terms, distinct_coverage_geometries=len(cached))


def run(*, quiet=False, max_seconds=30, progress=None):
    if (type(quiet) is not bool or type(max_seconds) not in (int, float)
            or not np.isfinite(max_seconds) or not 0 < max_seconds <= 30):
        raise ValueError('Fixed quiet/full scope and explicit budget at most30 seconds required')
    reference, data, measure = read_inputs()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Fixed-marker coverage control exhausted its budget')
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
        raw_regional_banks_built=0, parameter_fits=0, requested_controls=2,
        reference_artifact=BASE_ARTIFACT, reference_artifact_sha256=BASE_ARTIFACT_SHA,
        frozen_base_source_sha256=BASE_SOURCE_SHA, coverage_source_sha256=COVER_SOURCE_SHA,
        state_measure=measure, prefix_log2=reference['prefix_log2'],
        compute_budget_seconds=max_seconds, directed=False,
        whole_code_certificate=False, complete_requested_controls=False, controls=[])
    try:
        full_weights = base.raw_holder_weights()
        matrix, _ = base.column_matrix(full_weights, base.GROUPS, measure,
            arity=base.ARITY, check_time=check_time)
        word = base.matrix_poly_power(matrix, base.COLUMNS, check_time)
        initial_row = np.logaddexp2.reduce(word[0], axis=0)
        poly_bound = reference['prefix_log2']+float(np.logaddexp2.reduce(initial_row[base.FLOOR:]))
        expected = next(row['log2_screen'] for row in reference['controls'] if row['mode'] == 'polynomial')
        if abs(poly_bound-expected) > 1e-7:
            raise ArithmeticError('Frozen full matrix-polynomial tail did not reproduce')
        marker = base.adjacent_floor_marker(initial_row, base.FLOOR)
        result.update(log2_marker=marker, marker_log_as_exact_float=str(F(marker)),
            marker_rule='max(0,log(e0 coeff1296 ones)-log(e0 coeff1297 ones)) of FULL ordered word',
            adjacent_full_coefficients=initial_row[base.FLOOR:base.FLOOR+2].tolist(),
            reproduced_full_polynomial_log2=poly_bound)
        if progress:
            progress(dict(stage='full-word-marker', log2_marker=marker, reproduced_bound=poly_bound))
        check_time()
        weights = base.raw_holder_weights(quiet=True) if quiet else full_weights
        covered, unrestricted, costs = paired_columns(weights, base.GROUPS, measure,
            arity=base.ARITY, log_marker=marker, check_time=check_time)
        original, _ = base.column_matrix(weights, base.GROUPS, measure,
            arity=base.ARITY, log_marker=marker, check_time=check_time)
        if not np.allclose(unrestricted, original, rtol=0., atol=1e-8):
            raise ArithmeticError('Matched unrestricted scalar operator did not reproduce')
        for label, scalar_matrix in (('unrestricted', unrestricted), ('cross-covered', covered)):
            check_time()
            power = base.log_matrix_power(scalar_matrix, base.COLUMNS)
            scalar = float(np.logaddexp2.reduce(power[0]))-base.FLOOR*marker
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
