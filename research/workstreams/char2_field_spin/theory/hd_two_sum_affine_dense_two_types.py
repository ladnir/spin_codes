"""One frozen-witness FLOAT gate for all21 affine h2 dense arity pairs.

Only r17/source(8,0), cohort(0,15,2). The two fixed row sets range over
all219^2 ordered assignments. This is not a directed certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_affine_dense import affine_class_moments, ARITY_CLASSES, CONSTRUCTION, LABELS
from hd_two_sum_affine_dense_two import two_mask_masses, occupancy_vectors, branch_caps, arity_pairs
from hd_xor_dense_holder_classes import class_log_weights
from hd_xor_dense_two_feature import positive_bivariate_tail
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import column_log_coefficients
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


DIRECTORY = Path(__file__).resolve().parent
FIT_NAME = 'hd_two_sum_affine_dense_h2_a8_8_0_fit24.json'
FIT_SHA256 = '358a73269e326fb50fcc845d4c2511e9cdcb71b6e9f741dbc499dec25dfbd071'
Z = F(1116892707587883, 1125899906842624)
LOGP = F(-40)


def frozen_witness():
    raw = (DIRECTORY/FIT_NAME).read_bytes()
    if sha256(raw).hexdigest() != FIT_SHA256:
        raise ValueError('The approved frozen h2 fit artifact changed')
    source = json.loads(raw)
    if (source['construction'] != CONSTRUCTION or source['cohort_counts'] != [0, 15, 2]
            or source['source_matching_counts'] != [8, 0] or source['dense_arities'] != [8, 8]
            or source['field_order'] != 256 or source['labels'] != list(LABELS)
            or source['forbidden_weight'] != LIMIT or source['total_columns'] != 32
            or source['winning_literal']['mark'] != str(Z)
            or source['winning_literal']['log2_potential'] != str(LOGP)
            or source['winning_witness_hits_boundary'] is not True
            or source['directed'] is not False):
        raise ValueError('Frozen fit identity or scope mismatch')
    return source


def type_partition():
    """Source placements, not per-column multiplicities or extra Holder mass."""
    classes = dict(ARITY_CLASSES)
    rows = []
    for a, b, multiplicity in arity_pairs():
        rows.append(dict(dense_arities=[a, b], dense_type_multiplicity=multiplicity,
            fixed_row_sets_left=list(classes[a]), fixed_row_sets_right=list(classes[b]),
            ordered_fixed_dense_rowset_assignments=multiplicity*len(classes[a])*len(classes[b])))
    if (len(rows) != 21 or sum(row['ordered_fixed_dense_rowset_assignments'] for row in rows) != 219**2
            or sorted(mask for _, masks in ARITY_CLASSES for mask in masks)
                != [mask for mask in range(256) if mask.bit_count() >= 3]):
        raise ArithmeticError('Incomplete fixed-row-set type partition')
    return rows


def run(*, max_seconds=120, progress=None):
    if not 0 < max_seconds <= 120:
        raise ValueError('One frozen gate with at most120 seconds is authorized')
    source = frozen_witness()
    started = perf_counter()
    entries = []
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Bounded all21 affine h2 FLOAT budget exhausted')
    def emit(row):
        if progress:
            progress(dict(row, elapsed_seconds=perf_counter()-started))
    source_count = matching_joint_counts(8192, 17)[8, 0]
    result = dict(construction=CONSTRUCTION, field_order=256, field_polynomial='0x11b',
        labels=list(LABELS), forbidden_weight=LIMIT, active_groups=17, cohort_counts=[0, 15, 2],
        source_matching_counts=[8, 0], source_classes=1, source_supports=source_count,
        source_colorings=comb(17, 2), sparse_pair_assignments=28**15, pair_swaps=True,
        total_columns=32, sparse_support_floor=510, dense_nonempty_column_floor=34,
        dense_floor_is_aggregate_relaxation=True, sparse_normalizer_exponent=480, holder_exponent=32,
        expected_arity_pairs=21, expected_ordered_fixed_dense_rowset_assignments=219**2,
        dense_rowsets_charged_once_inside_holder=True, extra_dense_beta_factor=False,
        extra_dense_rowset_factor=False, fixed_rows_across_columns=True, fixed_source_geometry=True,
        iid_rowset_resampling_assumed=False, affine_mask_kernel_rebuilt=True, old_xor_mask_counts_used=False,
        labelled_dense_occupancies_retained=True, sparse_hole_histogram_exact=True,
        potential_quotient_after_hole_and_dense_mixtures=True, one_global_mark_and_potential=True,
        regional_matrices_are_upper_bounds=True, ordered_regional_matrices=True,
        pair_multiplicity_strata_used=False, resource_maximum_used=False,
        projective_normalizations=1, endpoint_potential_factors=1, alpha=1,
        mark=str(Z), log2_potential=str(LOGP), witness_banks=1, fitted_new_witnesses=0,
        witness_fit_filename=FIT_NAME, witness_fit_sha256=FIT_SHA256,
        witness_origin_hits_boundary=True, witness_origin_boundary_coordinates=['mark-lower', 'log2-potential-upper'],
        witness_origin_bounds=source['bounds'], directed=False, old_certificates_transferred=False,
        whole_code_certificate=False, complete_source_coverage=False,
        complete_dense_rowset_coverage=False, complete_selected_cohort_and_source=False,
        compute_budget_seconds=max_seconds, entries=entries)
    try:
        weights, moment_pins = {}, {}
        for arity, masks in ARITY_CLASSES:
            check_time()
            nums, denominator = affine_class_moments(masks)
            logs = class_log_weights(nums, denominator)
            weight = np.exp2(logs)
            if np.any(~np.isfinite(weight)) or np.any((weight == 0) & np.isfinite(logs)):
                raise ArithmeticError('FLOAT Holder weight overflow/underflow is not permitted')
            weights[arity] = weight
            moment_pins[arity] = sha256(json.dumps([list(nums), denominator], separators=(',', ':')).encode()).hexdigest()
        kernel_seconds = perf_counter()-started
        raw = pair_regions_screen(float(Z), 17)
        regions = ordered_region_tables(raw, 17, 8, 0)
        vectors = occupancy_vectors(regions, float(LOGP), check_time=check_time, progress=emit)
        operator_seconds = perf_counter()-started-kernel_seconds
        prefix = (log2(source_count)+log2(comb(17, 2))+15*log2(28)-log2(255)
                  -LIMIT*log2(float(Z))-float(LOGP)-480*log2(255))
        for row in type_partition():
            check_time()
            a, b = row['dense_arities']
            masses = np.asarray(two_mask_masses(weights[a], weights[b]), dtype=float)
            if np.any(~np.isfinite(masses)) or np.any(masses < 0):
                raise ArithmeticError('Invalid FLOAT positive mask convolution')
            with np.errstate(divide='ignore'):
                mass_logs = np.log2(masses)
            caps = branch_caps(vectors, mass_logs, float(LOGP))
            branches = np.array([column_log_coefficients(caps[:, e], 15) for e in range(3)])
            tail, peak = positive_bivariate_tail(branches, 32, 510, 34)
            check_time()
            row.update(failure_log2_screen=prefix+tail+log2(row['dense_type_multiplicity']), **peak,
                exact_moment_sha256_left=moment_pins[a], exact_moment_sha256_right=moment_pins[b])
            entries.append(row)
            emit(dict(stage='arity-pair-complete', completed=len(entries), total=21,
                dense_arities=row['dense_arities'], failure_log2_screen=row['failure_log2_screen'], **peak))
        # The a8 x a8 cell is a fresh FLOAT control of the pinned winner.
        if abs(entries[-1]['failure_log2_screen']-source['winning_literal']['failure_log2_screen']) > 1e-8:
            raise ArithmeticError('Frozen a8 x a8 result did not reproduce')
        result.update(status='completed', complete_dense_rowset_coverage=True,
            complete_selected_cohort_and_source=True, ordered_fixed_dense_rowset_assignments=219**2,
            evaluated_arity_pairs=len(entries),
            failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])),
            failed_arity_pairs=sum(row['failure_log2_screen'] >= -40 for row in entries),
            kernel_seconds=kernel_seconds, operator_seconds=operator_seconds,
            occupancy_vectors=3**8, occupied_banks=16, retained_vector_bank_doubles=int(vectors.size))
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error), evaluated_arity_pairs=len(entries))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite frozen all21 FLOAT evidence')
    result = run(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False))
