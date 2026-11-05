"""Bounded FLOAT gate for one fixed-mask dense group and sixteen row pairs.

This gate keeps the dense row set A fixed across32 columns. For each dense
output mask T, sparse occupancy m in labelled region i uses R_i[m+1_T(i)].
The sparse cancellation mixture and then exact N_A(T)/256^a mixture precede
the common-potential quotient. Positive X/Y coefficients count sparse
support and nonempty dense columns, respectively. The retained tail is
X>=544,Y>=17. No dense product domination or geometry resampling occurs.

The initial executable scope is source matching counts(8,0), mu=16 among
the16 two-row groups. Its resource pair constraint is inactive. This is a
FLOAT proposal, not directed arithmetic and not a whole-code certificate.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_fixed_mask_floor import beta_log2, comparison_counts
from hd_xor_geometry_strata import category_resources, geometry_stratum_counts
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import inherited_witnesses, ordered_region_tables, SOURCE
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT


def offset_resource_vectors(regions, log_potential, lower, upper):
    """All dense-mask offset envelopes for a box and a sparse occupancy sum.

    Output[T,j] bounds product_i R_i[m_i+1_T(i)] v componentwise for every
    lower<=m_i<=upper with sum_i(m_i-lower)=j. Sparse geometry constraints
    may be relaxed to this box. Each dense T remains fixed throughout the
    ordered product. The local mask index is little-endian by region.
    """
    regions = np.asarray(regions, dtype=float)
    if (regions.ndim != 4 or regions.shape[2:] != (2, 2)
            or not 1 <= len(regions) <= 8
            or type(lower) is not int or type(upper) is not int
            or not 0 <= lower <= upper < regions.shape[1]-1
            or not np.isfinite(log_potential)
            or np.any(np.isnan(regions)) or np.any(np.isposinf(regions))):
        raise ValueError('Expected labelled log matrices including each dense offset')
    width = upper-lower
    values = np.array([[[0., log_potential]]])
    for table in reversed(regions):
        following = np.full((2*len(values), values.shape[1]+width, 2), -np.inf)
        for bit in (0, 1):
            target = following[bit::2]
            for delta in range(width+1):
                matrix = table[lower+delta+bit]
                action = np.logaddexp2(matrix[None, None, :, 0]+values[:, :, 0, None],
                                      matrix[None, None, :, 1]+values[:, :, 1, None])
                block = target[:, delta:delta+values.shape[1]]
                np.maximum(block, action, out=block)
        values = following
    return values


def mixed_mask_vectors(regions, log_potential, groups=16, mu=16, progress=None):
    """Mix sparse cancellations, leaving dense masks and both states intact.

    This bounded implementation requires mu=groups. Then k<=mu makes the
    pair bound inactive for every category. Using the result for a smaller
    mu would remain an upper bound but would discard that stratum's gain;
    callers must not describe it as the refined small-mu resource cap.
    """
    if type(groups) is not int or not 1 <= groups <= 16 or mu != groups or len(regions) != 8:
        raise ValueError('The bounded offset evaluator currently requires mu=groups<=16')
    kinds = categories(groups)
    requests = {}
    for index, (single, double) in enumerate(kinds):
        for cancelled in range(double+1):
            lower, upper, limit, target = category_resources(single, double, cancelled, mu)
            if limit < 2*upper:
                raise ArithmeticError('This envelope cannot drop an active pair constraint')
            log_weight = (log2(comb(double, cancelled))
                          +(double-cancelled)*log2(Q-2)-double*log2(Q-1))
            requests.setdefault((lower, upper), []).append((index, target, log_weight))
    result = np.full((len(kinds), 1 << len(regions), 2), -np.inf)
    started = perf_counter()
    for finished, ((lower, upper), uses) in enumerate(sorted(requests.items()), 1):
        vectors = offset_resource_vectors(regions, log_potential, lower, upper)
        for index, target, log_weight in uses:
            result[index] = np.logaddexp2(result[index], vectors[:, target]+log_weight)
        if progress and (finished % 25 == 0 or finished == len(requests)):
            progress(dict(stage='dense-mask-resource-bank', completed=finished,
                          total=len(requests), elapsed_seconds=perf_counter()-started))
    return result


def fixed_mask_branch_caps(vectors, active_mask, log_potential, group=8, q=256):
    """Exact dense-mask weights are mixed before either potential quotient."""
    vectors = np.asarray(vectors, dtype=float)
    if (vectors.ndim != 3 or vectors.shape[1:] != (1 << group, 2)
            or np.any(np.isnan(vectors)) or np.any(np.isposinf(vectors))
            or not np.isfinite(log_potential)):
        raise ValueError('Expected category/dense-mask/state log vectors')
    counts = comparison_counts(active_mask, group, q)
    normalizer = active_mask.bit_count()*log2(q)
    weights = np.array([log2(count)-normalizer if count else -np.inf for count in counts[1:]])
    empty_vector = vectors[:, 0]-normalizer
    nonempty_vector = np.logaddexp2.reduce(vectors[:, 1:]+weights[None, :, None], axis=1)
    quotient = np.array([0., log_potential])
    return np.max(empty_vector-quotient, axis=1), np.max(nonempty_vector-quotient, axis=1)


def column_log_coefficients(caps, groups, q=256):
    """Positive source-support polynomial, including (q-1)^X once."""
    caps = np.asarray(caps, dtype=float)
    kinds = categories(groups)
    if caps.shape != (len(kinds),) or q < 2 or np.any(np.isnan(caps)) or np.any(np.isposinf(caps)):
        raise ValueError('Complete finite-or-zero category caps required')
    result = np.full(2*groups+1, -np.inf)
    for (single, double), value in zip(kinds, caps):
        degree = single+2*double
        count = comb(groups, single)*comb(groups-single, double)*2**single
        result[degree] = np.logaddexp2(result[degree], value+log2(count)+degree*log2(q-1))
    return result


def log_convolve(left, right):
    """Positive convolution without FFT cancellation or a coefficient floor."""
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if left.ndim != 1 or right.ndim != 1 or min(len(left), len(right)) < 1:
        raise ValueError('Nonempty one-dimensional coefficient vectors required')
    if len(left) > len(right):
        left, right = right, left
    result = np.full(len(left)+len(right)-1, -np.inf)
    for degree, coefficient in enumerate(left):
        if coefficient != -np.inf:
            target = result[degree:degree+len(right)]
            np.logaddexp2(target, coefficient+right, out=target)
    return result


def log_powers(polynomial, maximum):
    if type(maximum) is not int or maximum < 0:
        raise ValueError('A nonnegative integer power bound is required')
    result = [np.array([0.])]
    for _ in range(maximum):
        result.append(log_convolve(result[-1], polynomial))
    return result


def bivariate_log_tail(empty, nonempty, columns=32, x_floor=544, y_floor=17,
                       *, empty_powers=None, empty_shift=0.):
    """Exact positive X/Y tail in FLOAT log arithmetic; no Chernoff Y bound.

    Scalar branch caps commute, so [Y^h](P0+Y P1)^n=C(n,h)P0^(n-h)P1^h.
    An optional cached power bank represents unshifted P0; empty_shift is
    added for every empty column. It is an algebraic scaling, not a fit.
    """
    if (type(columns) is not int or columns < 0 or type(y_floor) is not int
            or not 0 <= y_floor <= columns or type(x_floor) is not int or x_floor < 0):
        raise ValueError('Invalid column count or coefficient floors')
    if empty_powers is None:
        if empty_shift:
            raise ValueError('An empty shift requires a cached unshifted power bank')
        empty_powers = log_powers(empty, columns-y_floor)
    if len(empty_powers) <= columns-y_floor:
        raise ValueError('Incomplete empty-branch power bank')
    dense_power = np.array([0.])
    rows = []
    for degree in range(columns+1):
        if degree >= y_floor:
            polynomial = log_convolve(empty_powers[columns-degree], dense_power)
            tail = polynomial[x_floor:]
            log_mass = (float(np.logaddexp2.reduce(tail)) if len(tail) else -np.inf)
            log_mass += log2(comb(columns, degree))+(columns-degree)*empty_shift
            rows.append(dict(nonempty_dense_columns=degree, log2_mass=log_mass,
                             largest_sparse_degree=int(x_floor+np.argmax(tail)) if len(tail) else None))
        if degree < columns:
            dense_power = log_convolve(dense_power, nonempty)
    return float(np.logaddexp2.reduce([row['log2_mass'] for row in rows])), rows


def run(*, source_path=SOURCE, maximum_row_sets=219, progress=None):
    """One frozen witness, one ordered source class, and mu16 only."""
    if type(maximum_row_sets) is not int or not 1 <= maximum_row_sets <= 219:
        raise ValueError('Bounded row-set count must lie in1..219')
    started = perf_counter()
    witness_rows, witness_hash = inherited_witnesses(source_path)
    choices = [(z, potential) for low, high, z, potential in witness_rows[8, 0] if low <= 16 <= high]
    if len(choices) != 1:
        raise ValueError('One inherited mu16 witness is required')
    z, logp = choices[0]
    raw = pair_regions_screen(float(z), 17)
    tables = ordered_region_tables(raw, 17, 8, 0)
    vectors = mixed_mask_vectors(tables, float(logp), progress=progress)
    base_empty_caps = np.max(vectors[:, 0]-[0., float(logp)], axis=1)
    base_empty = column_log_coefficients(base_empty_caps, 16)
    empty_powers = log_powers(base_empty, 15)
    source_count = matching_joint_counts(REGION_LENGTH, 17)[8, 0]
    geometry_count = geometry_stratum_counts(16)[16]
    prefix = (log2(source_count)+log2(17)+log2(geometry_count)-log2(Q-1)
              -LIMIT*log2(float(z))-float(logp)-512*log2(Q-1))
    masks = [mask for mask in range(1, 256) if mask.bit_count() >= 3][:maximum_row_sets]
    rows = []
    for number, active_mask in enumerate(masks, 1):
        row_start = perf_counter()
        empty, nonempty = fixed_mask_branch_caps(vectors, active_mask, float(logp))
        empty_poly = column_log_coefficients(empty, 16)
        nonempty_poly = column_log_coefficients(nonempty, 16)
        empty_shift = -8*active_mask.bit_count()
        if not np.allclose(empty_poly, base_empty+empty_shift, rtol=0, atol=1e-11):
            raise ArithmeticError('Empty dense-mask normalization disagrees with cached powers')
        tail, by_dense_degree = bivariate_log_tail(empty_poly, nonempty_poly,
            empty_powers=empty_powers, empty_shift=empty_shift)
        peak = max(by_dense_degree, key=lambda row: row['log2_mass'])
        row = dict(active_row_mask=active_mask, active_constituents=active_mask.bit_count(),
                   failure_log2_screen=prefix+beta_log2(active_mask.bit_count())+tail,
                   largest_dense_degree=peak['nonempty_dense_columns'],
                   largest_sparse_degree_at_peak=peak['largest_sparse_degree'],
                   dense_degree_log2_masses=by_dense_degree,
                   elapsed_seconds=perf_counter()-row_start)
        rows.append(row)
        if progress and (number == 1 or number % 20 == 0 or number == len(masks)):
            progress({key: value for key, value in row.items() if key != 'dense_degree_log2_masses'} |
                     dict(stage='fixed-row-set-tail', completed=number, total=len(masks)))
    arities = []
    for active in range(3, 9):
        selected = [row for row in rows if row['active_constituents'] == active]
        if selected:
            arities.append(dict(active_constituents=active, fixed_row_sets=len(selected),
                failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in selected]))))
    return dict(construction=CONSTRUCTION, active_groups=17, cohort_counts=[0, 16, 1],
        source_matching_counts=[8, 0], source_classes=1, source_supports=source_count,
        dense_source_colorings=17, mu_range=[16, 16], two_row_geometries=geometry_count,
        fixed_row_sets=len(rows), all_dense_row_sets=len(rows) == 219,
        active_constituent_counts=[3, 4, 5, 6, 7, 8], sparse_support_floor=544,
        dense_nonempty_column_floor=17, total_columns=32, sparse_support_normalizer_exponent=512,
        mark=str(z), log2_potential=str(logp), witness_source_sha256=witness_hash,
        field_order=Q, forbidden_weight=LIMIT, pair_swaps=True,
        exact_uniform_dense_mask_mixture=True, dense_input_includes_zero=True,
        dense_beta_formula='q^(32a)/(q-1)^(16a)', conditional_nonzero_amplitude_operator=True,
        dense_row_set_fixed_across_columns=True, sparse_geometry_fixed_across_columns=True,
        masks_and_cancellations_mixed_before_potential=True, same_potential_both_dense_branches=True,
        exact_positive_bivariate_tail=True, chernoff_dense_floor=False,
        product_density_comparison_used=False, projective_normalizations=1,
        endpoint_potential_factors=1, alpha=1, directed=False, fitted_new_witnesses=0,
        proposal_requires_fresh_directed_replay=True, complete_mixed_class=False,
        whole_code_certificate=False, entries=rows, arity_summaries=arities,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows])),
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum-row-sets', type=int, default=219)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen diagnostic artifact')
    result = run(maximum_row_sets=args.maximum_row_sets,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False), flush=True)
