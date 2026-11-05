"""Bounded FLOAT g16 h2 gate through two dense-weight coefficients.

Fix source(8,0) at r17 and the cohort of15 two-row and two dense groups.
Each dense row set remains fixed across16 columns. Closed-kernel Holder
weights depend only on output-mask weight, so two polynomial indices sum
all labelled mask pairs without a 4^16 table. The regional matrices still
act in order. Empty/nonempty group-column markers remain separate until
after the complete sparse-hole and dense-mask mixture.

The105 unordered arity pairs partition all65399^2 ordered dense row-set
assignments. Unequal arities get one factor2 outside the column power.
A separate control merges all arities inside one exact Holder moment;
that control gets no arity or row-set multiplicity outside the moment.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_closed_kernel import arity_holder_moments
from hd_two_sum_gf16_arity_dense import arity_weights
from hd_two_sum_gf16_full_dense import Z, LOGP, column_vectors, _check_deadline
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS, ordered_regions
from hd_xor_dense_holder_classes import class_log_weights
from hd_xor_dense_two_feature import positive_bivariate_tail
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_mask_floor import column_log_coefficients
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


def arity_pairs():
    return tuple((a, b, 1 if a == b else 2) for a in range(3, 17) for b in range(a, 17))


def aggregate_moments(group=16, columns=16, dimension=8):
    """Exact all-arity moments with one common denominator; fixed A only."""
    if (group not in (8, 16) or type(columns) is not int or type(dimension) is not int
            or not 0 < dimension < columns):
        raise ValueError('Checked group and positive proper constituent dimension required')
    denominator = 255**((columns-dimension)*group)
    result = [0]*(group+1)
    for arity in range(3, group+1):
        nums, divisor = arity_holder_moments(arity, group, columns, dimension)
        scale, remainder = divmod(denominator, divisor)
        if remainder:
            raise ArithmeticError('The exact arity denominator must divide the common denominator')
        for w, number in enumerate(nums):
            result[w] += number*scale
    return tuple(result), denominator


def two_weight_vectors(regions, log_potential, sparse_groups=15, q=256, *, deadline=None):
    """Category/w1/w2/state vectors: sum all labelled dense mask pairs.

    At each region the four bit pairs00,01,10,11 increment the two weight
    indices and use occupancies t-h,t-h+1,t-h+1,t-h+2. The two middle
    branches share an action but both contribute. Binomial hole factors
    count assignments of labelled sparse holes. No dense weights, maxima,
    or binomial mask multiplicities are inserted in this recurrence.
    """
    tables = np.asarray(regions, dtype=float)
    if (type(sparse_groups) is not int or not 0 <= sparse_groups <= 15
            or type(q) is not int or q < 3
            or tables.ndim != 4 or tables.shape[2:] != (2, 2)
            or not 1 <= len(tables) <= 16 or len(tables) >= q-1
            or tables.shape[1] <= sparse_groups+2 or not np.isfinite(log_potential)
            or np.any(np.isnan(tables)) or np.any(np.isposinf(tables))):
        raise ValueError('Complete ordered log matrices and a sparse count0..15 required')
    kinds = categories(sparse_groups)
    indices = {kind: index for index, kind in enumerate(kinds)}
    result = np.full((len(kinds), len(tables)+1, len(tables)+1, 2), -np.inf)
    for occupied in range(sparse_groups+1):
        _check_deadline(deadline)
        values = np.full((1, 1, occupied+1, 2), -np.inf)
        values[0, 0, 0] = [0., log_potential]
        binomials = [np.array([log2(comb(k, holes)) for k in range(holes, occupied+1)])
                     for holes in range(occupied+1)]
        for table in reversed(tables):
            width = len(values)
            following = np.full((width+1, width+1, occupied+1, 2), -np.inf)
            for holes in range(occupied+1):
                previous = values[:, :, :occupied+1-holes]
                for occupancy in range(3):
                    matrix = table[occupied-holes+occupancy]
                    action = np.logaddexp2(matrix[:, 0]+previous[..., 0, None],
                                          matrix[:, 1]+previous[..., 1, None])
                    action += binomials[holes][None, None, :, None]
                    shifts = ((0, 0),) if occupancy == 0 else (((1, 1),) if occupancy == 2 else ((0, 1), (1, 0)))
                    for left, right in shifts:
                        target = following[left:left+width, right:right+width, holes:]
                        np.logaddexp2(target, action, out=target)
            values = following
        for double in range(occupied+1):
            probabilities = np.array([log2(comb(double, k))
                +(double-k)*log2(q-1-len(tables))-double*log2(q-1) for k in range(double+1)])
            result[indices[occupied-double, double]] = np.logaddexp2.reduce(
                values[:, :, :double+1]+probabilities[None, None, :, None], axis=2)
    _check_deadline(deadline)
    return result


def pair_branch_caps(vectors, log_potential, left_weights, right_weights):
    """Mix weights by0/1/2 nonempty dense masks before one potential quotient."""
    vectors = np.asarray(vectors, dtype=float)
    left, right = np.asarray(left_weights, dtype=float), np.asarray(right_weights, dtype=float)
    if (left.ndim != 1 or right.ndim != 1 or min(len(left), len(right)) < 2
            or vectors.ndim != 4 or vectors.shape[1:] != (len(left), len(right), 2)
            or not np.isfinite(log_potential)
            or any(np.any(np.isnan(a)) or np.any(np.isposinf(a)) for a in (vectors, left, right))):
        raise ValueError('Complete finite-or-zero two-weight vectors and one shared potential required')
    weighted = vectors+left[None, :, None, None]+right[None, None, :, None]
    empty = weighted[:, 0, 0]
    once = np.logaddexp2(np.logaddexp2.reduce(weighted[:, 1:, 0], axis=1),
                        np.logaddexp2.reduce(weighted[:, 0, 1:], axis=1))
    twice = np.logaddexp2.reduce(np.logaddexp2.reduce(weighted[:, 1:, 1:], axis=2), axis=1)
    return np.max(np.stack((empty, once, twice), axis=1)-np.array([0., log_potential]), axis=2)


def run(*, max_seconds=120, progress=None):
    """One frozen witness,105 arity pairs plus one merged-arity control."""
    if (isinstance(max_seconds, bool) or not isinstance(max_seconds, (int, float))
            or not np.isfinite(max_seconds) or not 0 < max_seconds <= 120):
        raise ValueError('A positive budget at most120 seconds is required')
    started = perf_counter()
    deadline = started+max_seconds
    source_count = matching_joint_counts(8192, 17)[8, 0]
    fixed_rows = sum(comb(16, a) for a in range(3, 17))
    pairs = arity_pairs()
    if sum(mult*comb(16, a)*comb(16, b) for a, b, mult in pairs) != fixed_rows**2:
        raise ArithmeticError('Unordered arity classes must count every ordered row-set assignment once')
    result = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        field_polynomial='0x11b', active_groups=17, cohort_counts=[0, 15, 2],
        source_matching_counts=[8, 0], source_supports=source_count,
        source_colorings=comb(17, 2), sparse_pair_assignments=120**15,
        forbidden_weight=LIMIT, total_columns=16, constituent_dimension=8,
        sparse_support_floor=270, dense_nonempty_column_floor=18,
        dense_nonempty_floor_is_aggregate_relaxation=True,
        sparse_normalizer_exponent=240, holder_exponent=16,
        fixed_dense_row_sets_per_group=fixed_rows, ordered_dense_rowset_assignments=fixed_rows**2,
        requested_arity_pairs=105, dense_weight_states=[17, 17], enumerated_dense_masks=0,
        dense_rowsets_charged_once_inside_holder=True, extra_dense_rowset_factor=False,
        unordered_arity_multiplicity_once_outside_power=True, extra_target_mask_binomial=False,
        extra_dense_beta_factor=False, fixed_rows_across_columns=True,
        fixed_source_geometry=True, iid_rowset_resampling_assumed=False, pair_swaps=True,
        labelled_region_order_preserved=True, hole_mixture_exact=True,
        potential_quotient_after_hole_and_dense_mixtures=True,
        projective_normalizations=1, endpoint_potential_factors=1,
        mark=str(Z), log2_potential=str(LOGP), alpha=1, witness_banks=1,
        fitted_new_witnesses=0, raw_regional_banks=0,
        complete_source_coverage=False, directed=False, whole_code_certificate=False,
        old_certificates_transferred=False, max_seconds=max_seconds, entries=[])
    try:
        _check_deadline(deadline)
        raw = pair_regions_screen(float(Z), 17)
        result['raw_regional_banks'] = 1
        tables = ordered_regions(raw, 17, 8, 0)
        vectors = two_weight_vectors(tables, float(LOGP), deadline=deadline)
        result['vector_bank_seconds'] = perf_counter()-started
        weights, pins = {}, {}
        for a in range(3, 17):
            weights[a], pins[a] = arity_weights(a)
        full = pair_branch_caps(vectors, float(LOGP), weights[16], weights[16])
        independent = np.max(column_vectors(tables, float(LOGP), 15, 2, deadline=deadline)
                             -np.array([0., float(LOGP)]), axis=2)-16*log2(255)
        difference = float(np.max(np.abs(full-independent)))
        if difference > 1e-9:
            raise ArithmeticError('The two-weight bank differs from the independent full-arity seen-count DP')
        result['fullarity_maximum_column_cap_difference'] = difference
        prefix = (log2(source_count)+log2(comb(17, 2))+15*log2(120)
                  -log2(255)-LIMIT*log2(float(Z))-float(LOGP)-240*log2(255))
        for a, b, mult in pairs:
            _check_deadline(deadline)
            caps = pair_branch_caps(vectors, float(LOGP), weights[a], weights[b])
            branches = np.array([column_log_coefficients(caps[:, e], 15) for e in range(3)])
            tail, peak = positive_bivariate_tail(branches, 16, 270, 18)
            row = dict(active_constituents=[a, b], ordered_arity_assignments=mult,
                fixed_row_sets=[comb(16, a), comb(16, b)], exact_moments_sha256=[pins[a], pins[b]],
                failure_log2_screen=prefix+log2(mult)+tail, **peak)
            result['entries'].append(row)
            if progress:
                progress(dict(stage='arity-pair-complete', **row))
        _check_deadline(deadline)
        nums, denominator = aggregate_moments()
        merged_weights = class_log_weights(nums, denominator, 16)
        merged_caps = pair_branch_caps(vectors, float(LOGP), merged_weights, merged_weights)
        branches = np.array([column_log_coefficients(merged_caps[:, e], 15) for e in range(3)])
        tail, peak = positive_bivariate_tail(branches, 16, 270, 18)
        result['aggregate_control'] = dict(failure_log2_screen=prefix+tail, **peak,
            holder_fixed_row_sets_per_group=fixed_rows, extra_arity_multiplicity=False,
            exact_moments_sha256=sha256(json.dumps([list(nums), denominator], separators=(',', ':')).encode()).hexdigest())
        _check_deadline(deadline)
        result.update(status='complete-selected-cohort-and-source', complete_selected_cohort_and_source=True,
            failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in result['entries']])))
    except TimeoutError:
        result.update(status='compute-budget-exceeded', complete_selected_cohort_and_source=False)
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=120)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded FLOAT artifact')
    result = run(max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, 'hd_two_sum_closed_kernel.py', 'hd_two_sum_gf16_arity_dense.py',
                     'hd_two_sum_gf16_full_dense.py', 'hd_two_sum_gf16_sparse.py',
                     'hd_xor_dense_holder_classes.py', 'hd_xor_dense_two_feature.py',
                     'hd_xor_mixed_mask_floor.py', 'hd_xor_pair_route.py', 'hd_xor_group_geometry.py')}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False))
