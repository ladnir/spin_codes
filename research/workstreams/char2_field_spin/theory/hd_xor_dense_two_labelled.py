"""One bounded FLOAT h2 control retaining all labelled dense occupancies.

The frozen olda4/arity-five type is used for both dense source groups.
Ternary suffix sharing and exact reachability pruning change computation,
not the resource envelope. No directed arithmetic or full type sweep.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_dense_holder_classes import class_moments, class_log_weights
from hd_xor_dense_two_feature import (frozen_type_and_witness, positive_bivariate_tail,
    COMPARATOR_NAME, COMPARATOR_SHA256)
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import column_log_coefficients
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION


DIRECTORY = Path(__file__).resolve().parent


def occupancy_convolution(mask_weights, groups):
    """Positive whole-mask convolution, indexed by (nonempty count, code).

    Code digits are little-endian labelled regional occupancies, in base
    groups+1. Two binary mask codes add without carries in base three.
    Integer/Fraction weights retain exact arithmetic for the toy tests.
    """
    size = len(mask_weights)
    length = size.bit_length()-1
    if (size != 1 << length or length < 1 or type(groups) is not int or groups not in (1, 2)
            or any(value < 0 for value in mask_weights)):
        raise ValueError('Complete nonnegative labelled mask weights and one/two groups required')
    base = groups+1
    codes = [sum(((mask >> i) & 1)*base**i for i in range(length)) for mask in range(size)]
    result = [[0 for _ in range(base**length)] for _ in range(groups+1)]
    if groups == 1:
        for mask, value in enumerate(mask_weights):
            result[int(mask != 0)][codes[mask]] += value
    else:
        for left, lw in enumerate(mask_weights):
            if not lw:
                continue
            for right, rw in enumerate(mask_weights):
                if rw:
                    result[int(left != 0)+int(right != 0)][codes[left]+codes[right]] += lw*rw
    return result


def reachable_sums(length, width, targets):
    """Exact partial-sum down-set for the specified final resource targets."""
    targets = tuple(sorted(set(targets)))
    if (type(length) is not int or length < 1 or type(width) is not int or width < 0
            or not targets or any(type(t) is not int or not 0 <= t <= length*width for t in targets)):
        raise ValueError('Attainable integer final sums are required')
    return tuple(tuple(j for j in range(stage*width+1)
                       if any(j <= target <= j+(length-stage)*width for target in targets))
                 for stage in range(length+1))


def offset_vectors(regions, log_potential, lower, upper, dense_groups, targets, *, check_time=None):
    """Shared binary/ternary suffix DP; output[code,target_index,state]."""
    tables = np.asarray(regions, dtype=float)
    if (tables.ndim != 4 or tables.shape[2:] != (2, 2) or not 1 <= len(tables) <= 8
            or type(dense_groups) is not int or dense_groups not in (1, 2)
            or type(lower) is not int or type(upper) is not int
            or not 0 <= lower <= upper < tables.shape[1]-dense_groups
            or not np.isfinite(log_potential) or np.any(np.isnan(tables))
            or np.any(np.isposinf(tables))):
        raise ValueError('Complete labelled log matrices including dense offsets required')
    width, base = upper-lower, dense_groups+1
    keep = reachable_sums(len(tables), width, targets)
    values = np.array([[[0., log_potential]]])
    candidates, largest_array = 0, values.size
    for stage, table in enumerate(reversed(tables), 1):
        if check_time:
            check_time()
        final = np.asarray(keep[stage], dtype=int)
        lookup = np.full((stage-1)*width+1, -1, dtype=int)
        lookup[list(keep[stage-1])] = np.arange(len(keep[stage-1]))
        following = np.full((base**stage, len(final), 2), -np.inf)
        largest_array = max(largest_array, following.size)
        for delta in range(width+1):
            predecessor = final-delta
            valid = (predecessor >= 0) & (predecessor < len(lookup))
            positions = np.flatnonzero(valid)
            indices = lookup[predecessor[valid]]
            positions, indices = positions[indices >= 0], indices[indices >= 0]
            if not len(indices):
                continue
            suffix = values[:, indices, :]
            for bit in range(base):
                matrix = table[lower+delta+bit]
                action = np.logaddexp2(matrix[None, None, :, 0]+suffix[:, :, 0, None],
                                      matrix[None, None, :, 1]+suffix[:, :, 1, None])
                block = following[bit::base]
                block[:, positions, :] = np.maximum(block[:, positions, :], action)
                candidates += len(values)*len(indices)
        values = following
    return values, dict(two_state_matrix_candidates=candidates,
                        largest_single_array_doubles=largest_array,
                        target_sums=list(keep[-1]))


def labelled_caps(regions, log_potential, mass_logs, sparse_groups, dense_groups,
                  *, q=256, check_time=None, progress=None):
    """Mix dense occupancies and cancellation coins before a single quotient."""
    tables, masses = np.asarray(regions), np.asarray(mass_logs)
    length = len(tables)
    if (length % 2 or not 2 <= length <= 8 or sparse_groups < 1
            or dense_groups not in (1, 2) or masses.shape != (dense_groups+1, (dense_groups+1)**length)
            or np.any(np.isnan(masses)) or np.any(np.isposinf(masses)) or q <= 2):
        raise ValueError('Complete labelled offset masses and even regional count required')
    kinds = categories(sparse_groups)
    requests = {}
    for index, (single, double) in enumerate(kinds):
        for cancelled in range(double+1):
            lower, upper = double-cancelled, single+double
            target = (length-1)*single+2*cancelled
            log_probability = (log2(comb(double, cancelled))
                +(double-cancelled)*log2(q-2)-double*log2(q-1))
            requests.setdefault((lower, upper), []).append((index, target, log_probability))
    vectors = np.full((len(kinds), dense_groups+1, 2), -np.inf)
    candidates, largest = 0, 0
    for done, ((lower, upper), uses) in enumerate(sorted(requests.items()), 1):
        if check_time:
            check_time()
        values, stats = offset_vectors(tables, log_potential, lower, upper, dense_groups,
                                       [target for _, target, _ in uses], check_time=check_time)
        candidates += stats['two_state_matrix_candidates']
        largest = max(largest, stats['largest_single_array_doubles'])
        index_of = {target: i for i, target in enumerate(stats['target_sums'])}
        # Positive summation over complete dense occupancy vectors may precede
        # summation over cancellation coins. Neither is moved through max.
        mixed = np.stack([np.logaddexp2.reduce(values+mass[:, None, None], axis=0)
                          for mass in masses], axis=1)
        for index, target, log_probability in uses:
            vectors[index] = np.logaddexp2(vectors[index], mixed[index_of[target]]+log_probability)
        if progress and (done % 25 == 0 or done == len(requests)):
            progress(dict(stage='labelled-offset-bank', completed=done, total=len(requests)))
    caps = np.max(vectors-np.array([0., log_potential]), axis=2)
    return caps, dict(resource_intervals=len(requests), sparse_categories=len(kinds),
                      two_state_matrix_candidates=candidates, largest_single_array_doubles=largest,
                      labelled_dense_occupancy_vectors=masses.shape[1])


def evaluate(dense_groups, *, seconds=60., progress=None):
    if type(dense_groups) is not int or dense_groups not in (1, 2) or not 0 < seconds <= 60:
        raise ValueError('Only the frozen h1/h2 controls, with a budget at most60s, are allowed')
    started = perf_counter()
    last_progress = {}
    def emit(row):
        last_progress.update(row)
        if progress:
            progress(dict(row, elapsed_seconds=perf_counter()-started))
    def check_time():
        if perf_counter()-started > seconds:
            raise TimeoutError('The bounded FLOAT control exceeded its compute budget')
    sparse_groups = 17-dense_groups
    masks, z, logp = frozen_type_and_witness('a4_mask240_fit')
    metadata = dict(construction=CONSTRUCTION, pair_swaps=True, cohort_counts=[0, sparse_groups, dense_groups],
        active_groups=17, source_matching_counts=[8, 0], field_order=256, forbidden_weight=LIMIT,
        total_columns=32, sparse_pair_assignments=28**sparse_groups,
        source_supports=matching_joint_counts(8192, 17)[8, 0], source_colorings=comb(17, dense_groups),
        fixed_dense_type=dict(witness_class='a4_mask240_fit', arity=5, active_row_masks=list(masks)),
        fixed_row_sets_per_dense_group=36, ordered_row_set_assignments=36**dense_groups,
        dense_rowsets_charged_once_inside_holder=True, extra_dense_beta_factor=False,
        extra_dense_rowset_factor=False, dense_rowsets_fixed_across_columns=True,
        sparse_geometry_fixed_across_columns=True, iid_rowset_resampling_assumed=False,
        one_global_mark_and_potential=True, masks_and_cancellations_mixed_before_potential=True,
        ordered_regional_matrices=True, labelled_dense_occupancies_retained=True,
        dense_nonempty_floor_is_aggregate_relaxation=dense_groups == 2,
        sparse_support_floor=34*sparse_groups, dense_nonempty_column_floor=17*dense_groups,
        sparse_normalizer_exponent=32*sparse_groups, projective_normalizations=1,
        endpoint_potential_factors=1, mark=str(z), log2_potential=str(logp), alpha=1,
        holder_exponent=32, fitted_new_witnesses=0, directed=False, whole_code_certificate=False,
        compute_budget_seconds=seconds)
    try:
        numerators, denominator = class_moments(masks)
        weight_logs = class_log_weights(numerators, denominator)
        weights = np.exp2(weight_logs)
        if np.any(~np.isfinite(weights)) or np.any((weights == 0) & np.isfinite(weight_logs)):
            raise ArithmeticError('FLOAT Holder weight overflow/underflow is not permitted')
        masses = np.asarray(occupancy_convolution(weights, dense_groups), dtype=float)
        if np.any(~np.isfinite(masses)):
            raise ArithmeticError('FLOAT mask convolution overflow is not permitted')
        with np.errstate(divide='ignore'):
            mass_logs = np.log2(masses)
        raw = pair_regions_screen(float(z), 17)
        regions = ordered_region_tables(raw, 17, 8, 0)
        check_time()
        caps, costs = labelled_caps(regions, float(logp), mass_logs, sparse_groups, dense_groups,
                                    check_time=check_time, progress=emit)
        operator_seconds = perf_counter()-started
        branches = np.array([column_log_coefficients(caps[:, e], sparse_groups)
                             for e in range(dense_groups+1)])
        check_time()
        tail, peak = positive_bivariate_tail(branches, 32, 34*sparse_groups, 17*dense_groups)
        check_time()
        prefix = (log2(metadata['source_supports'])+log2(comb(17, dense_groups))
                  +sparse_groups*log2(28)-log2(255)-LIMIT*log2(float(z))-float(logp)
                  -32*sparse_groups*log2(255))
        metadata.update(status='completed', failure_log2_screen=prefix+tail,
            operator_seconds=operator_seconds, **costs, **peak,
            exact_holder_moments_sha256=sha256(json.dumps([list(numerators), denominator],
                                                        separators=(',', ':')).encode()).hexdigest())
        if dense_groups == 1:
            raw = (DIRECTORY/COMPARATOR_NAME).read_bytes()
            if sha256(raw).hexdigest() != COMPARATOR_SHA256:
                raise ValueError('Frozen h1 comparator changed')
            row = next(row for row in json.loads(raw)['entries']
                       if row['witness_label'] == 'a4_mask240_fit' and row['active_constituents'] == 5)
            difference = metadata['failure_log2_screen']-row['failure_log2_screen']
            if abs(difference) > 1e-8:
                raise ArithmeticError('The labelled h1 control does not reproduce its frozen comparator')
            metadata.update(comparator_sha256=COMPARATOR_SHA256,
                labelled_h1_comparator_log2_screen=row['failure_log2_screen'],
                comparator_difference_bits=difference)
    except TimeoutError as error:
        metadata.update(status='compute-budget-exceeded', reason=str(error), last_progress=last_progress)
    metadata['elapsed_seconds'] = perf_counter()-started
    return metadata


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dense-groups', type=int, choices=(1, 2), required=True)
    parser.add_argument('--seconds', type=float, default=60.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded control artifact')
    result = evaluate(args.dense_groups, seconds=args.seconds,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)
