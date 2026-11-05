"""FLOAT g16 h1 gate for every fixed dense row set of arity3..16.

The closed kernel makes the arity-wise Holder weight depend only on the
output mask's cardinality. A positive coefficient DP sums the labelled
regional actions by that cardinality, without commuting the matrices.
Row sets remain fixed across16 columns. Their multiplicity occurs inside
the exact Holder moment, never as an extra binomial factor.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_closed_kernel import arity_holder_moments
from hd_two_sum_gf16_full_dense import Z, LOGP, _check_deadline
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS, ordered_regions
from hd_xor_dense_holder_classes import class_log_weights, mask_branch_caps
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import column_log_coefficients, bivariate_log_tail
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


def weight_vectors(regions, log_potential, sparse_groups=16, q=256, *, deadline=None):
    """Sum exact hole-mixture vectors over all dense masks of each weight.

    The weight variable counts branch choices; its coefficient already
    sums all labelled masks. Do not multiply by C(number_of_regions,w).
    """
    tables = np.asarray(regions, dtype=float)
    if (type(sparse_groups) is not int or not 0 <= sparse_groups <= 16
            or type(q) is not int or q < 3
            or tables.ndim != 4 or tables.shape[2:] != (2, 2)
            or not 1 <= len(tables) <= 16 or len(tables) >= q-1
            or tables.shape[1] <= sparse_groups+1 or not np.isfinite(log_potential)
            or np.any(np.isnan(tables)) or np.any(np.isposinf(tables))):
        raise ValueError('Complete finite-or-zero ordered matrices and a bounded sparse count required')
    kinds = categories(sparse_groups)
    indices = {kind: i for i, kind in enumerate(kinds)}
    result = np.full((len(kinds), len(tables)+1, 2), -np.inf)
    for occupied in range(sparse_groups+1):
        _check_deadline(deadline)
        values = np.full((1, occupied+1, 2), -np.inf)
        values[0, 0] = [0., log_potential]
        binomials = [np.array([log2(comb(k, holes)) for k in range(holes, occupied+1)])
                     for holes in range(occupied+1)]
        for table in reversed(tables):
            following = np.full((len(values)+1, occupied+1, 2), -np.inf)
            for bit in (0, 1):
                for holes in range(occupied+1):
                    previous = values[:, :occupied+1-holes]
                    matrix = table[occupied-holes+bit]
                    action = np.logaddexp2(matrix[:, 0]+previous[..., 0, None],
                                          matrix[:, 1]+previous[..., 1, None])
                    target = following[bit:bit+len(values), holes:]
                    np.logaddexp2(target, action+binomials[holes][None, :, None], out=target)
            values = following
        for double in range(occupied+1):
            probabilities = np.array([log2(comb(double, k))
                +(double-k)*log2(q-1-len(tables))-double*log2(q-1) for k in range(double+1)])
            result[indices[occupied-double, double]] = np.logaddexp2.reduce(
                values[:, :double+1]+probabilities[None, :, None], axis=1)
    _check_deadline(deadline)
    return result


def arity_weights(arity, group=16, columns=16, dimension=8):
    """FLOAT logs of exact per-mask moments; no target-mask multiplicity."""
    nums, denominator = arity_holder_moments(arity, group, columns, dimension)
    weights = class_log_weights(nums, denominator, columns)
    pin = sha256(json.dumps([list(nums), denominator], separators=(',', ':')).encode()).hexdigest()
    return weights, pin


def g8_control(raw, *, deadline=None):
    """Reproduce the existing labelled g8 six-arity control at the same literal."""
    from hd_two_sum_affine_dense import affine_class_moments, dense_mask_vectors
    tables = ordered_region_tables(raw, 17, 8, 0)
    labelled = dense_mask_vectors(tables, float(LOGP))
    reduced = weight_vectors(tables, float(LOGP), deadline=deadline)
    source_count = matching_joint_counts(8192, 17)[8, 0]
    prefix = (log2(source_count)+log2(17)+16*log2(28)-log2(255)
              -LIMIT*log2(float(Z))-float(LOGP)-512*log2(255))
    rows = []
    for arity in range(3, 9):
        _check_deadline(deadline)
        members = tuple(mask for mask in range(256) if mask.bit_count() == arity)
        nums, denominator = affine_class_moments(members)
        old_weights = class_log_weights(nums, denominator)
        new_weights, _ = arity_weights(arity, 8, 32, 16)
        old_caps = mask_branch_caps(labelled, float(LOGP), old_weights)
        new_caps = mask_branch_caps(reduced, float(LOGP), new_weights)
        difference = max(float(np.max(np.abs(a-b))) for a, b in zip(old_caps, new_caps))
        if difference > 1e-9:
            raise ArithmeticError('Cardinality aggregation changed a labelled g8 column cap')
        bounds = []
        for caps in (old_caps, new_caps):
            tail, _ = bivariate_log_tail(*(column_log_coefficients(branch, 16) for branch in caps),
                                        columns=32, x_floor=544, y_floor=17)
            bounds.append(prefix+tail)
        if abs(bounds[0]-bounds[1]) > 1e-8:
            raise ArithmeticError('Cardinality aggregation changed the labelled g8 gate')
        rows.append(dict(active_constituents=arity, labelled_failure_log2_screen=bounds[0],
            cardinality_failure_log2_screen=bounds[1], maximum_column_cap_difference=difference))
    return dict(complete=True, arities=list(range(3, 9)), rows=rows,
        labelled_union_log2_screen=float(np.logaddexp2.reduce([r['labelled_failure_log2_screen'] for r in rows])),
        cardinality_union_log2_screen=float(np.logaddexp2.reduce([r['cardinality_failure_log2_screen'] for r in rows])),
        exact_same_literal=True, directed=False)


def run(*, max_seconds=120, progress=None):
    """One source(8,0), one frozen witness, all14 h1 arities; no fitting."""
    if (isinstance(max_seconds, bool) or not isinstance(max_seconds, (int, float))
            or not np.isfinite(max_seconds) or not 0 < max_seconds <= 120):
        raise ValueError('A positive compute budget at most120 seconds is required')
    started = perf_counter()
    deadline = started+max_seconds
    count = matching_joint_counts(8192, 17)[8, 0]
    result = dict(construction=CONSTRUCTION, labels=list(LABELS), field_order=256,
        field_polynomial='0x11b', active_groups=17, cohort_counts=[0, 16, 1],
        source_matching_counts=[8, 0], source_supports=count, source_colorings=17,
        sparse_pair_assignments=120**16, forbidden_weight=LIMIT, total_columns=16,
        constituent_dimension=8, sparse_support_floor=288, dense_nonempty_column_floor=9,
        sparse_normalizer_exponent=256, holder_exponent=16, arities=list(range(3, 17)),
        fixed_dense_row_sets=sum(comb(16, a) for a in range(3, 17)),
        deterministic_arity_classes=14, dense_weight_states=17, enumerated_dense_masks=0,
        dense_rowsets_charged_once_inside_holder=True, extra_dense_rowset_factor=False,
        extra_target_mask_binomial=False, extra_dense_beta_factor=False,
        fixed_rows_across_columns=True, fixed_source_geometry=True,
        iid_rowset_resampling_assumed=False, pair_swaps=True,
        labelled_region_order_preserved=True, hole_mixture_exact=True,
        potential_quotient_after_hole_and_dense_mixtures=True,
        projective_normalizations=1, endpoint_potential_factors=1,
        mark=str(Z), log2_potential=str(LOGP), alpha=1, witness_banks=1,
        fitted_new_witnesses=0, raw_regional_banks=0,
        complete_source_coverage=False, complete_higher_dense_cohort_coverage=False,
        directed=False, whole_code_certificate=False, old_certificates_transferred=False,
        max_seconds=max_seconds, entries=[])
    try:
        _check_deadline(deadline)
        raw = pair_regions_screen(float(Z), 17)
        result['raw_regional_banks'] = 1
        result['g8_labelled_control'] = g8_control(raw, deadline=deadline)
        if progress:
            progress(dict(stage='g8-control', **{key: value for key, value in result['g8_labelled_control'].items()
                                              if key != 'rows'}))
        tables = ordered_regions(raw, 17, 8, 0)
        vectors = weight_vectors(tables, float(LOGP), deadline=deadline)
        prefix = (log2(count)+log2(17)+16*log2(120)-log2(255)
                  -LIMIT*log2(float(Z))-float(LOGP)-256*log2(255))
        for arity in range(3, 17):
            _check_deadline(deadline)
            weights, moment_pin = arity_weights(arity)
            empty, nonempty = mask_branch_caps(vectors, float(LOGP), weights)
            tail, degrees = bivariate_log_tail(column_log_coefficients(empty, 16),
                column_log_coefficients(nonempty, 16), columns=16, x_floor=288, y_floor=9)
            peak = max(degrees, key=lambda row: row['log2_mass'])
            row = dict(active_constituents=arity, fixed_row_sets=comb(16, arity),
                exact_moments_sha256=moment_pin, failure_log2_screen=prefix+tail,
                largest_dense_degree=peak['nonempty_dense_columns'],
                largest_sparse_degree=peak['largest_sparse_degree'])
            result['entries'].append(row)
            if progress:
                progress(dict(stage='arity-complete', **row))
        _check_deadline(deadline)
        result.update(status='complete-selected-cohort-and-source', complete_selected_cohort_and_source=True,
            failure_log2_screen=float(np.logaddexp2.reduce([r['failure_log2_screen'] for r in result['entries']])))
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
        for name in (Path(__file__).name, 'hd_two_sum_closed_kernel.py', 'hd_two_sum_gf16_full_dense.py',
                     'hd_two_sum_gf16_sparse.py', 'hd_two_sum_affine_dense.py', 'hd_two_sum_affine_outer.py',
                     'hd_xor_dense_holder_classes.py', 'hd_xor_mixed_mask_floor.py',
                     'hd_xor_pair_route.py', 'hd_xor_group_geometry.py')}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('entries', 'g8_labelled_control')}, allow_nan=False))
