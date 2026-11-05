"""Bounded affine two-sum h1 dense-cohort FLOAT gate, not a certificate.

Scope: r17/source(8,0), sixteen two-row groups plus one dense group, all
219 fixed dense row sets of arity3..8. The affine u8..15 kernel is rebuilt;
no XOR2 support-count kernel is reused. Six deterministic arity classes
use Holder only across fixed row sets, never resampling them by column.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_affine_outer import CONSTRUCTION, LABELS, WITNESSES, uniform_mask_counts, matrix
from hd_xor_dense_holder_classes import class_log_weights, mask_branch_caps
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import ordered_region_tables
from hd_xor_mixed_mask_floor import column_log_coefficients, bivariate_log_tail
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


DIRECTORY = Path(__file__).resolve().parent
ALL_DENSE = tuple(mask for mask in range(256) if mask.bit_count() >= 3)
ARITY_CLASSES = tuple((arity, tuple(mask for mask in ALL_DENSE if mask.bit_count() == arity))
                      for arity in range(3, 9))


@lru_cache(None)
def affine_class_moments(masks, columns=32, dimension=16):
    """Exact affine-kernel W(T)^columns, including dense comparison factors.

    W(T)^32=sum_A N_A(T)^32/255^(16|A|). All N_A are computed from the
    affine matrix, not imported from the previous XOR2 construction.
    """
    masks = tuple(masks)
    if (not masks or len(set(masks)) != len(masks) or tuple(sorted(masks)) != masks
            or any(type(mask) is not int or not 1 <= mask < 256 for mask in masks)
            or type(columns) is not int or type(dimension) is not int
            or not 0 < dimension < columns or LABELS != tuple(range(8, 16))):
        raise ValueError('Distinct fixed row sets and the exact g8 affine kernel required')
    denominator = 255**((columns-dimension)*8)
    numerators = [0]*256
    for active in masks:
        factor = 255**((columns-dimension)*(8-active.bit_count()))
        counts = uniform_mask_counts(active)
        if counts[0] != 1 or sum(counts) != 256**active.bit_count() or min(counts) < 0:
            raise ArithmeticError('Invalid affine uniform-input support counts')
        for target, count in enumerate(counts):
            if count:
                numerators[target] += factor*count**columns
    return tuple(numerators), denominator


def offset_hole_suffix(regions, log_potential, occupied):
    """All labelled dense masks and all sparse hole counts in one suffix DP.

    Output[T,k] sums the ordered action over allocations of k labelled
    holes to the regions. At region i use R_i(occupied-h+1_T(i)) with
    binomial factor C(k,h). No potential quotient is taken here.
    """
    tables = np.asarray(regions, dtype=float)
    if (tables.ndim != 4 or tables.shape[2:] != (2, 2) or not 1 <= len(tables) <= 8
            or type(occupied) is not int or not 0 <= occupied < tables.shape[1]-1
            or not np.isfinite(log_potential) or np.any(np.isnan(tables))
            or np.any(np.isposinf(tables))):
        raise ValueError('Complete labelled log matrices including a dense offset required')
    log_binomials = [np.array([log2(comb(k, h)) for k in range(h, occupied+1)])
                    for h in range(occupied+1)]
    values = np.full((1, occupied+1, 2), -np.inf)
    values[0, 0] = [0., log_potential]
    for table in reversed(tables):
        following = np.full((2*len(values), occupied+1, 2), -np.inf)
        for bit in (0, 1):
            target = following[bit::2]
            for holes in range(occupied+1):
                prior = values[:, :occupied+1-holes]
                local = table[occupied-holes+bit]
                action = np.logaddexp2(local[None, None, :, 0]+prior[:, :, 0, None],
                                      local[None, None, :, 1]+prior[:, :, 1, None])
                block = target[:, holes:]
                np.logaddexp2(block, action+log_binomials[holes][None, :, None], out=block)
        values = following
    return values


def dense_mask_vectors(regions, log_potential, sparse_groups=16, q=256, progress=None):
    """Exactly mix the sparse joint hole law, keeping dense T labelled."""
    if type(sparse_groups) is not int or not 0 <= sparse_groups <= 16 or q-1 <= len(regions):
        raise ValueError('Bounded sparse group count and positive full-support mass required')
    kinds = categories(sparse_groups)
    indices = {kind: index for index, kind in enumerate(kinds)}
    result = np.full((len(kinds), 1 << len(regions), 2), -np.inf)
    for occupied in range(sparse_groups+1):
        # The suffix recurrence depends on t=s+d, not on d separately.
        # Reuse one bank for every d<=t, then mix its own binomial law.
        vectors = offset_hole_suffix(regions, log_potential, occupied)
        for double in range(occupied+1):
            probabilities = np.array([log2(comb(double, k))
                +(double-k)*log2(q-1-len(regions))-double*log2(q-1)
                for k in range(double+1)])
            result[indices[occupied-double, double]] = np.logaddexp2.reduce(
                vectors[:, :double+1]+probabilities[None, :, None], axis=1)
        if progress and (occupied % 4 == 0 or occupied == sparse_groups):
            progress(dict(stage='dense-labelled-hole-bank', completed=occupied+1, total=sparse_groups+1))
    return result


def run(*, mark=WITNESSES[0][1], log_potential=WITNESSES[0][2], witness_label='inherited', progress=None):
    z, logp = F(mark), F(log_potential)
    if not 0 < z < 1 or logp >= 0:
        raise ValueError('An exact rational mark in(0,1) and negative log potential required')
    started = perf_counter()
    weights, moment_pins = {}, {}
    for arity, masks in ARITY_CLASSES:
        nums, denominator = affine_class_moments(masks)
        weights[arity] = class_log_weights(nums, denominator)
        moment_pins[arity] = sha256(json.dumps([list(nums), denominator], separators=(',', ':')).encode()).hexdigest()
    kernel_seconds = perf_counter()-started
    raw = pair_regions_screen(float(z), 17)
    tables = ordered_region_tables(raw, 17, 8, 0)
    vectors = dense_mask_vectors(tables, float(logp), progress=progress)
    operator_seconds = perf_counter()-started-kernel_seconds
    source_count = matching_joint_counts(8192, 17)[8, 0]
    prefix = (log2(source_count)+log2(17)+16*log2(28)-log2(255)
              -LIMIT*log2(float(z))-float(logp)-512*log2(255))
    entries = []
    for arity, masks in ARITY_CLASSES:
        empty, nonempty = mask_branch_caps(vectors, float(logp), weights[arity])
        tail, degrees = bivariate_log_tail(column_log_coefficients(empty, 16),
            column_log_coefficients(nonempty, 16), columns=32, x_floor=544, y_floor=17)
        largest = max(degrees, key=lambda row: row['log2_mass'])
        row = dict(active_constituents=arity, active_row_masks=list(masks), fixed_row_sets=len(masks),
            exact_moments_sha256=moment_pins[arity], failure_log2_screen=prefix+tail,
            largest_dense_degree=largest['nonempty_dense_columns'],
            largest_sparse_degree=largest['largest_sparse_degree'])
        entries.append(row)
        if progress:
            progress(dict(stage='arity-complete', **{key: value for key, value in row.items()
                                                    if key != 'active_row_masks'}))
    return dict(construction=CONSTRUCTION, labels=list(LABELS), field_polynomial='0x11b',
        field_order=256, forbidden_weight=LIMIT, active_groups=17, cohort_counts=[0, 16, 1],
        source_matching_counts=[8, 0], source_classes=1, source_supports=source_count,
        source_colorings=17, sparse_pair_assignments=28**16, pair_swaps=True,
        total_columns=32, sparse_support_floor=544, dense_nonempty_column_floor=17,
        sparse_normalizer_exponent=512, holder_exponent=32, fixed_row_sets=219,
        deterministic_classes=6, dense_type_partition='all-fixed-row-sets-by-arity',
        dense_rowsets_charged_once_inside_holder=True, extra_dense_beta_factor=False,
        extra_dense_rowset_factor=False, fixed_rows_across_columns=True,
        fixed_source_geometry=True, iid_rowset_resampling_assumed=False,
        affine_mask_kernel_rebuilt=True, old_xor_mask_counts_used=False,
        affine_matrix_sha256=sha256(bytes(value for row in matrix() for value in row)).hexdigest(),
        labelled_dense_masks_retained=True, sparse_hole_histogram_exact=True,
        potential_quotient_after_hole_and_dense_mixtures=True, one_global_mark_and_potential=True,
        regional_matrices_are_upper_bounds=True, ordered_regional_matrices=True,
        pair_multiplicity_strata_used=False, resource_maximum_used=False,
        projective_normalizations=1, endpoint_potential_factors=1, alpha=1,
        mark=str(z), log2_potential=str(logp), witness_label=witness_label,
        witness_banks=1, fitted_new_witnesses=0, entries=entries,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])),
        directed=False, old_certificates_transferred=False, whole_code_certificate=False,
        complete_selected_cohort_and_source=True, complete_source_coverage=False,
        kernel_seconds=kernel_seconds, operator_seconds=operator_seconds, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded research artifact')
    result = run(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False))
