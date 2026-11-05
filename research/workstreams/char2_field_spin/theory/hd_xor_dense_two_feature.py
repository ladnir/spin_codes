"""Bounded FLOAT h2 interface using positive dense feature/resource envelopes.

Only one fixed type B is used: the36 arity-five row sets assigned to the
old a4 witness. Holder counts fixed row sets inside W_B; it does not
resample them. Feature convolution loses labelled mask placement but not
the ordered regional matrix chronology. No directed certificate is claimed.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_dense_holder_classes import class_moments, class_log_weights, frozen_classes
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import column_log_coefficients
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
COMPARATOR_NAME = 'hd_xor_dense_holder_8_0_h1_arity12.json'
COMPARATOR_SHA256 = '4577166d2aefbe1912ac6bc096c8f08a3604439c03f72ff0a2dd7ad4fef74a0e'
CONTROLS = {'h1-olda4': (1, 'a4_mask240_fit'),
            'h2-olda4': (2, 'a4_mask240_fit'), 'h2-a8': (2, 'a8_fit')}


def feature_convolution(mask_weights, groups):
    """Positive exact algebra: (nonempty masks,total weight,full masks).

    Rational/integer inputs give exact rational/integer outputs. Actual
    FLOAT controls supply W_B approximations; no symmetry factor C(8,w)
    is used because the supplied labelled mask weights need not be symmetric.
    """
    size = len(mask_weights)
    width = size.bit_length()-1
    if size != 1 << width or width < 1 or type(groups) is not int or not 0 <= groups <= 2:
        raise ValueError('Complete labelled mask weights and at most two dense groups required')
    one = {}
    for mask, value in enumerate(mask_weights):
        if value < 0:
            raise ValueError('Dense mask weights must be nonnegative')
        if value:
            key = int(mask != 0), mask.bit_count(), int(mask == size-1)
            one[key] = one.get(key, 0)+value
    result = {(0, 0, 0): 1}
    for _ in range(groups):
        following = {}
        for (e, t, f), value in result.items():
            for (de, dt, df), coefficient in one.items():
                key = e+de, t+dt, f+df
                following[key] = following.get(key, 0)+value*coefficient
        result = following
    return result


def ordered_box_vectors(regions, log_potential, lower, upper):
    """Componentwise suffix upper vectors for every shifted occupancy sum."""
    regions = np.asarray(regions, dtype=float)
    if (regions.ndim != 4 or regions.shape[2:] != (2, 2) or not 1 <= len(regions) <= 8
            or not 0 <= lower <= upper < regions.shape[1] or not np.isfinite(log_potential)
            or np.any(np.isnan(regions)) or np.any(np.isposinf(regions))):
        raise ValueError('Complete labelled log-matrix tables and a valid occupancy box required')
    width = upper-lower
    values = np.array([[0., log_potential]])
    for table in reversed(regions):
        following = np.full((len(values)+width, 2), -np.inf)
        for delta in range(width+1):
            matrix = table[lower+delta]
            action = np.logaddexp2(matrix[None, :, 0]+values[:, 0, None],
                                  matrix[None, :, 1]+values[:, 1, None])
            target = following[delta:delta+len(values)]
            np.maximum(target, action, out=target)
        values = following
    return values


def feature_caps(regions, log_potential, features, sparse_groups, dense_groups, q=256):
    """Mix cancellation and positive feature masses before the potential max."""
    length = len(regions)
    if length % 2 or not 0 <= sparse_groups <= 16 or not 1 <= dense_groups <= 2 or q <= 2:
        raise ValueError('Even region count, bounded group counts and q>2 required')
    tables = np.asarray(regions, dtype=float)
    @lru_cache(None)
    def resource(lower, upper):
        return ordered_box_vectors(tables, log_potential, lower, upper)
    result = np.full((len(categories(sparse_groups)), dense_groups+1, 2), -np.inf)
    for index, (single, double) in enumerate(categories(sparse_groups)):
        for cancelled in range(double+1):
            log_probability = (log2(comb(double, cancelled))
                +(double-cancelled)*log2(q-2)-double*log2(q-1))
            for (e, weight, full), mass in features.items():
                if mass <= 0 or not 0 <= full <= e <= dense_groups:
                    raise ValueError('Positive attainable dense feature masses required')
                lower, upper = double-cancelled+full, single+double+e
                target = (length-1)*single+2*cancelled+weight-length*full
                vectors = resource(lower, upper)
                if not 0 <= target < len(vectors):
                    raise ValueError('Dense feature lies outside its occupancy box')
                result[index, e] = np.logaddexp2(result[index, e],
                    vectors[target]+log_probability+log2(mass))
    caps = np.max(result-np.array([0., log_potential]), axis=2)
    return caps, dict(resource_intervals=resource.cache_info().currsize,
                      dense_features=len(features), sparse_categories=len(result))


def positive_bivariate_tail(branches, columns, x_floor, y_floor):
    """Positive ordered-column DP: all coefficients of P(X,Y)^columns.

    Iteration over labelled columns supplies the multinomial coefficients;
    no further binomial/multinomial factor is added. No FFT or subtraction.
    """
    branches = np.asarray(branches, dtype=float)
    if (branches.ndim != 2 or not branches.size or type(columns) is not int or columns < 0
            or type(x_floor) is not int or x_floor < 0 or type(y_floor) is not int or y_floor < 0
            or np.any(np.isnan(branches)) or np.any(np.isposinf(branches))):
        raise ValueError('Positive log-coefficient branches and nonnegative integer floors required')
    dy, dx = np.array(branches.shape)-1
    values = np.array([[0.]])
    terms = [(e, degree, value) for (e, degree), value in np.ndenumerate(branches) if value != -np.inf]
    for _ in range(columns):
        following = np.full((values.shape[0]+dy, values.shape[1]+dx), -np.inf)
        for e, degree, coefficient in terms:
            target = following[e:e+values.shape[0], degree:degree+values.shape[1]]
            np.logaddexp2(target, coefficient+values, out=target)
        values = following
    tail = values[y_floor:, x_floor:]
    if not tail.size:
        return -np.inf, dict(largest_dense_degree=None, largest_sparse_degree=None)
    peak = np.unravel_index(np.argmax(tail), tail.shape)
    return float(np.logaddexp2.reduce(tail.ravel())), dict(
        largest_dense_degree=int(y_floor+peak[0]), largest_sparse_degree=int(x_floor+peak[1]))


def frozen_type_and_witness(label):
    classes, _ = frozen_classes()
    members = next(masks for name, _, _, masks, _ in classes if name == 'a4_mask240_fit')
    members = tuple(mask for mask in members if mask.bit_count() == 5)
    if len(members) != 36:
        raise ValueError('The bounded type must contain the36 frozen olda4 arity-five row sets')
    z, logp = next((z, p) for name, z, p, _, _ in classes if name == label)
    return members, z, logp


def evaluate(control):
    if control not in CONTROLS:
        raise ValueError('Only the three named frozen controls are authorized')
    dense_groups, label = CONTROLS[control]
    sparse_groups = 17-dense_groups
    members, z, logp = frozen_type_and_witness(label)
    started = perf_counter()
    numerators, denominator = class_moments(members)
    logs = class_log_weights(numerators, denominator)
    weights = np.exp2(logs)
    if np.any(~np.isfinite(weights)) or np.any((weights == 0) & np.isfinite(logs)):
        raise ArithmeticError('FLOAT Holder weight overflow or underflow is not permitted')
    features = feature_convolution(weights, dense_groups)
    raw = pair_regions_screen(float(z), 17)
    regions = ordered_region_tables(raw, 17, 8, 0)
    caps, costs = feature_caps(regions, float(logp), features, sparse_groups, dense_groups)
    branches = np.array([column_log_coefficients(caps[:, e], sparse_groups)
                         for e in range(dense_groups+1)])
    operator_seconds = perf_counter()-started
    tail, peak = positive_bivariate_tail(branches, 32, 34*sparse_groups, 17*dense_groups)
    source_count = matching_joint_counts(8192, 17)[8, 0]
    prefix = (log2(source_count)+log2(comb(17, dense_groups))+sparse_groups*log2(28)
              -log2(255)-LIMIT*log2(float(z))-float(logp)-32*sparse_groups*log2(255))
    result = dict(control=control, construction=CONSTRUCTION, pair_swaps=True,
        cohort_counts=[0, sparse_groups, dense_groups], active_groups=17,
        source_matching_counts=[8, 0], field_order=256, forbidden_weight=LIMIT,
        total_columns=32, sparse_pair_assignments=28**sparse_groups,
        source_supports=source_count, source_colorings=comb(17, dense_groups),
        fixed_dense_type=dict(witness_class='a4_mask240_fit', arity=5, active_row_masks=list(members)),
        fixed_row_sets_per_dense_group=36, ordered_row_set_assignments=36**dense_groups,
        dense_rowsets_charged_once_inside_holder=True, extra_dense_beta_factor=False,
        extra_dense_rowset_factor=False, dense_rowsets_fixed_across_columns=True,
        sparse_geometry_fixed_across_columns=True, iid_rowset_resampling_assumed=False,
        one_global_mark_and_potential=True, masks_and_cancellations_mixed_before_potential=True,
        ordered_regional_matrices=True, dense_placement_box_sum_relaxation=True,
        dense_nonempty_floor_is_aggregate_relaxation=dense_groups > 1,
        sparse_support_floor=34*sparse_groups, dense_nonempty_column_floor=17*dense_groups,
        sparse_normalizer_exponent=32*sparse_groups, projective_normalizations=1,
        endpoint_potential_factors=1, mark=str(z), log2_potential=str(logp), witness_label=label,
        exact_holder_moment_numerators=[str(n) for n in numerators],
        exact_holder_moment_denominator=str(denominator), holder_exponent=32,
        failure_log2_screen=prefix+tail, **peak, **costs,
        operator_seconds=operator_seconds, elapsed_seconds=perf_counter()-started,
        fitted_new_witnesses=0, alpha=1, directed=False, whole_code_certificate=False)
    if dense_groups == 1:
        data = (DIRECTORY/COMPARATOR_NAME).read_bytes()
        if sha256(data).hexdigest() != COMPARATOR_SHA256:
            raise ValueError('Frozen h1 comparator changed')
        row = next(row for row in json.loads(data)['entries']
                   if row['witness_label'] == label and row['active_constituents'] == 5)
        result.update(labelled_mask_comparator_log2_screen=row['failure_log2_screen'],
            feature_relaxation_loss_bits=result['failure_log2_screen']-row['failure_log2_screen'],
            comparator_sha256=COMPARATOR_SHA256)
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--control', choices=tuple(CONTROLS), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded diagnostic artifact')
    result = evaluate(args.control)
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('exact_holder_moment_numerators', 'exact_holder_moment_denominator')}, allow_nan=False))
