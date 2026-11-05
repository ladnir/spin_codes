"""One frozen affine h2 arity-pair FLOAT gate with labelled occupancy sums.

The initial executable control is a8 x a8, r17/source(8,0), fifteen
two-row groups. It is a triage subclass, not all219^2 dense assignments.
No fitting, old-XOR kernels, or directed certificate are used.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_affine_dense import affine_class_moments, ARITY_CLASSES, CONSTRUCTION, LABELS
from hd_two_sum_affine_dense_frozen import Z, LOGP
from hd_xor_dense_holder_classes import class_log_weights
from hd_xor_dense_two_feature import positive_bivariate_tail
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import column_log_coefficients
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


def two_mask_masses(left, right):
    """Positive convolution by labelled ternary occupancy and nonempty count.

    Integer/Fraction inputs remain exact. Ordered dense group assignments
    enter here; no symmetry multiplier is implicit in this operation.
    """
    size = len(left)
    length = size.bit_length()-1
    if (len(right) != size or size != 1 << length or length < 1
            or any(x < 0 for x in left) or any(x < 0 for x in right)):
        raise ValueError('Two complete nonnegative labelled mask laws required')
    codes = [sum(((mask >> i) & 1)*3**i for i in range(length)) for mask in range(size)]
    result = [[0 for _ in range(3**length)] for _ in range(3)]
    for first, a in enumerate(left):
        if a:
            for second, b in enumerate(right):
                if b:
                    result[int(first != 0)+int(second != 0)][codes[first]+codes[second]] += a*b
    return result


def offset_hole_suffix(regions, log_potential, occupied, dense_groups=2, check_time=None):
    """Shared labelled dense suffixes; positive EGF over sparse hole labels."""
    tables = np.asarray(regions, dtype=float)
    if (tables.ndim != 4 or tables.shape[2:] != (2, 2) or not 1 <= len(tables) <= 8
            or type(dense_groups) is not int or dense_groups not in (1, 2)
            or type(occupied) is not int or not 0 <= occupied < tables.shape[1]-dense_groups
            or not np.isfinite(log_potential) or np.any(np.isnan(tables))
            or np.any(np.isposinf(tables))):
        raise ValueError('Complete labelled regional tables and bounded dense offsets required')
    base = dense_groups+1
    log_binomials = [np.array([log2(comb(k, h)) for k in range(h, occupied+1)])
                    for h in range(occupied+1)]
    values = np.full((1, occupied+1, 2), -np.inf)
    values[0, 0] = [0., log_potential]
    for table in reversed(tables):
        if check_time:
            check_time()
        following = np.full((base*len(values), occupied+1, 2), -np.inf)
        for bit in range(base):
            target = following[bit::base]
            for holes in range(occupied+1):
                prior = values[:, :occupied+1-holes]
                local = table[occupied-holes+bit]
                action = np.logaddexp2(local[None, None, :, 0]+prior[:, :, 0, None],
                                      local[None, None, :, 1]+prior[:, :, 1, None])
                block = target[:, holes:]
                np.logaddexp2(block, action+log_binomials[holes][None, :, None], out=block)
        values = following
    return values


def occupancy_vectors(regions, log_potential, sparse_groups=15, q=256, *, check_time=None, progress=None):
    if type(sparse_groups) is not int or not 0 <= sparse_groups <= 15 or q-1 <= len(regions):
        raise ValueError('At most15 sparse groups and positive full-support probability required')
    kinds = categories(sparse_groups)
    indices = {kind: index for index, kind in enumerate(kinds)}
    result = np.full((len(kinds), 3**len(regions), 2), -np.inf)
    for occupied in range(sparse_groups+1):
        vectors = offset_hole_suffix(regions, log_potential, occupied, check_time=check_time)
        for double in range(occupied+1):
            probabilities = np.array([log2(comb(double, k))
                +(double-k)*log2(q-1-len(regions))-double*log2(q-1)
                for k in range(double+1)])
            result[indices[occupied-double, double]] = np.logaddexp2.reduce(
                vectors[:, :double+1]+probabilities[None, :, None], axis=1)
        if progress:
            progress(dict(stage='ternary-hole-bank', completed=occupied+1, total=sparse_groups+1))
    return result


def arity_pairs():
    """Unordered types with multiplicity for their two source placements."""
    return tuple((a, b, 1 if a == b else 2) for a in range(3, 9) for b in range(a, 9))


def branch_caps(vectors, mass_logs, log_potential):
    vectors, masses = np.asarray(vectors), np.asarray(mass_logs)
    if (vectors.ndim != 3 or vectors.shape[2] != 2 or masses.shape != (3, vectors.shape[1])
            or not np.isfinite(log_potential) or np.any(np.isnan(vectors)) or np.any(np.isposinf(vectors))
            or np.any(np.isnan(masses)) or np.any(np.isposinf(masses))):
        raise ValueError('Complete labelled occupancy vectors and positive mass logs required')
    mixed = np.stack([np.logaddexp2.reduce(vectors+mass[None, :, None], axis=1)
                      for mass in masses], axis=1)
    return np.max(mixed-np.array([0., log_potential]), axis=2)


def run_one(arity_left=8, arity_right=8, *, max_seconds=120, progress=None):
    if (arity_left, arity_right) != (8, 8) or not 0 < max_seconds <= 120:
        raise ValueError('This executable triage uses only a8 x a8 and at most120s')
    started = perf_counter()
    last = {}
    def emit(row):
        last.update(row)
        if progress:
            progress(dict(row, elapsed_seconds=perf_counter()-started))
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Bounded affine h2 FLOAT budget exhausted')
    classes = dict(ARITY_CLASSES)
    source_count = matching_joint_counts(8192, 17)[8, 0]
    metadata = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        forbidden_weight=LIMIT, cohort_counts=[0, 15, 2], active_groups=17,
        source_matching_counts=[8, 0], source_supports=source_count, source_colorings=comb(17, 2),
        sparse_pair_assignments=28**15, dense_arities=[8, 8], dense_type_multiplicity=1,
        ordered_fixed_dense_rowset_assignments=1, all_dense_rowset_assignments=219**2,
        complete_dense_rowset_coverage=False, total_columns=32, sparse_support_floor=510,
        dense_nonempty_column_floor=34, dense_floor_is_aggregate_relaxation=True,
        sparse_normalizer_exponent=480, holder_exponent=32, pair_swaps=True,
        dense_rowsets_charged_once_inside_holder=True, extra_dense_beta_factor=False,
        extra_dense_rowset_factor=False, fixed_rows_across_columns=True, fixed_source_geometry=True,
        iid_rowset_resampling_assumed=False, affine_mask_kernel_rebuilt=True, old_xor_mask_counts_used=False,
        labelled_dense_occupancies_retained=True, sparse_hole_histogram_exact=True,
        potential_quotient_after_hole_and_dense_mixtures=True, one_global_mark_and_potential=True,
        regional_matrices_are_upper_bounds=True, ordered_regional_matrices=True,
        pair_multiplicity_strata_used=False, resource_maximum_used=False,
        projective_normalizations=1, endpoint_potential_factors=1, alpha=1,
        mark=str(Z), log2_potential=str(LOGP), witness_banks=1, fitted_new_witnesses=0,
        directed=False, old_certificates_transferred=False, whole_code_certificate=False,
        compute_budget_seconds=max_seconds)
    try:
        nums, denominator = affine_class_moments(classes[8])
        logs = class_log_weights(nums, denominator)
        weights = np.exp2(logs)
        if np.any(~np.isfinite(weights)) or np.any((weights == 0) & np.isfinite(logs)):
            raise ArithmeticError('FLOAT mask weight overflow/underflow is not permitted')
        masses = np.asarray(two_mask_masses(weights, weights), dtype=float)
        with np.errstate(divide='ignore'):
            mass_logs = np.log2(masses)
        raw = pair_regions_screen(float(Z), 17)
        tables = ordered_region_tables(raw, 17, 8, 0)
        check_time()
        vectors = occupancy_vectors(tables, float(LOGP), check_time=check_time, progress=emit)
        caps = branch_caps(vectors, mass_logs, float(LOGP))
        operator_seconds = perf_counter()-started
        branches = np.array([column_log_coefficients(caps[:, e], 15) for e in range(3)])
        check_time()
        tail, peak = positive_bivariate_tail(branches, 32, 510, 34)
        check_time()
        prefix = (log2(source_count)+log2(comb(17, 2))+15*log2(28)
            -log2(255)-LIMIT*log2(float(Z))-float(LOGP)-480*log2(255))
        metadata.update(status='completed', failure_log2_screen=prefix+tail, **peak,
            operator_seconds=operator_seconds, occupancy_vectors=3**8, occupied_banks=16,
            retained_vector_bank_doubles=int(vectors.size),
            affine_moment_sha256=sha256(json.dumps([list(nums), denominator], separators=(',', ':')).encode()).hexdigest())
    except TimeoutError as error:
        metadata.update(status='compute-budget-exceeded', reason=str(error), last_progress=last)
    metadata['elapsed_seconds'] = perf_counter()-started
    return metadata


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded proposal artifact')
    result = run_one(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(result, allow_nan=False))
