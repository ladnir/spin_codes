"""FLOAT full-fiber bounds for all labelled arities10..16 of H*.

Every source set A spans all five affine coordinates, so the actual
outside image is im(V_O), O=A^c. The inside difference code is ker C,
C=H*[O,A], because H* is involutory. Its information sets complement
full-rank column bases of C. Only ranks and counts are cached by affine
hull; each labelled A retains its ordered matrix word and full sum term.

This scaffold has not yet been independently validated or used in a
production control. It is parked while structural proof repairs are tested.
The single branch includes outside zero and uses X1. These raw bounds
contain no source, RS, projective or Chernoff prefix and are not certificates.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, log2

import numpy as np

from hd_bit_syndrome_arity9_gate import gap_geometry, quiet_cover_orbits
from hd_bit_syndrome_counted_general import counted_column_matrices, quiet_group_values_general
from hd_bit_syndrome_flat_geometry import flat_hull_mask, rank_cross_masks
from hd_bit_syndrome_gate import ordered_cover_products, weighted_cover_matrices
from hd_bit_syndrome_outer import CONSTRUCTION, affine_rank
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power

FULL = (1 << 16)-1


@lru_cache(maxsize=307)
def global_bases(outside_flat, rank):
    """Column bases of V_F B V_S^T on the full binary-labelled ground set.

    Overlap is allowed in this auxiliary bilinear rank. Only S subset A
    is used later, where the actual cross map has no identity contribution.
    """
    if type(rank) is not int or not 0 <= rank <= 5 or flat_hull_mask(outside_flat) != outside_flat:
        raise ValueError('A canonical outside flat and exact rank0..5 required')
    masks = []
    for positions in combinations(range(16), rank):
        mask = sum(1 << i for i in positions)
        if rank_cross_masks(outside_flat, mask) == rank:
            masks.append(mask)
    answer = np.array(masks, dtype=np.uint32)
    answer.setflags(write=False)
    return answer


def actual_cover(source):
    source = tuple(source)
    if (not 9 <= len(source) <= 16 or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('Nine to sixteen distinct sorted source labels required')
    mask = sum(1 << i for i in source)
    outside = FULL ^ mask
    rank = affine_rank(outside) if outside else 0
    flat = flat_hull_mask(outside)
    if rank_cross_masks(outside, mask) != rank:
        raise AssertionError('A full-affine-rank source must reach the full outside image')
    bases = global_bases(flat, rank)
    bases = bases[(bases & np.uint32(outside)) == 0]
    if not len(bases):
        raise AssertionError('The actual cross map needs a full-rank column basis')
    cover = tuple(F(len(bases)-int(np.count_nonzero(bases & np.uint32(1 << i))), len(bases))
                  for i in source)
    if sum(cover) != len(source)-rank:
        raise AssertionError('Actual information cover has the wrong dimension')
    return cover, rank


def prepare_geometry(arity, *, check_time=lambda: None):
    if type(arity) is not int or not 10 <= arity <= 16:
        raise ValueError('This production geometry covers arities10..16')
    if arity == 16:
        return dict(rowsets=(tuple(range(16)),), covers=((F(1),)*16,), classes=np.array([0]),
                    outside_caps=((1,),), outside_ranks=(0,), fiber_dimensions=(16,),
                    orbit_sizes=(1,), actual_cover_count=1)
    prepared = quiet_cover_orbits(16, arity, check_time=check_time)
    covers = list(prepared['covers'])
    ranks, dimensions, changed = [], [], 0
    for index, source in enumerate(prepared['rowsets']):
        if not index % 64:
            check_time()
        mask = sum(1 << i for i in source)
        if affine_rank(mask) != 5:
            raise AssertionError('At least nine source labels must have affine rank5')
        rank = affine_rank(FULL ^ mask)
        if rank < 5:
            covers[index], actual_rank = actual_cover(source)
            if actual_rank != rank:
                raise AssertionError('Actual and expected cross ranks disagree')
            changed += 1
        if sum(covers[index]) != arity-rank:
            raise AssertionError('Inside difference dimension mismatch')
        ranks.append(rank)
        dimensions.append(arity-rank)
    if len(covers) != comb(16, arity):
        raise AssertionError('Missing labelled source sets')
    prepared.update(covers=tuple(covers), outside_ranks=tuple(ranks),
                    fiber_dimensions=tuple(dimensions), actual_cover_count=changed)
    return prepared


def contract(prepared, measure, groups, *, columns=16, batch_size=128, check_time=lambda: None):
    rows, covers = prepared['rowsets'], prepared['covers']
    arity = len(rows[0])
    if type(groups) is not int or groups < 1 or type(columns) is not int or columns < 1:
        raise ValueError('Positive group and column counts required')
    if arity == 16:
        values, parents, labels, codes = (F(1),), tuple([-1]+list(range(16))), tuple([-1]+[0]*16), np.array([[16]])
    else:
        values, parents, labels, codes = gap_geometry(rows, covers, 16)
    check_time()
    bank = weighted_cover_matrices(measure, groups, values, q=256)
    products = ordered_cover_products(bank['cover_matrices'], parents, labels)
    terminal = np.stack([groups*quiet_group_values_general(
        [log2(c) if c else -np.inf for c in caps], measure, 0.) for caps in prepared['outside_caps']])
    column, costs = counted_column_matrices(terminal, prepared['classes'], measure, 16-arity,
        products, codes, batch_size=batch_size, check_time=check_time)
    word = log_matrix_power(column, columns)
    values = np.logaddexp2.reduce(word[:, 0, :], axis=1)
    return dict(rowsets=rows, per_rowset_word_log2=values,
                combined_word_log2=combine_row_bounds(values, groups),
                column_log2_matrices=column, costs=costs)


def full_row_bounds(arity, measure, groups=18, *, check_time=lambda: None):
    prepared = prepare_geometry(arity, check_time=check_time)
    answer = contract(prepared, measure, groups, check_time=check_time)
    answer['metadata'] = dict(outer_construction=CONSTRUCTION, arity=arity, groups=groups,
        columns=16, field_order=256, rowset_count=comb(16, arity),
        input_marker=1, outside_zero_included=True, young_split=False,
        outside_image='actual image equals im(V_O)', inside_fiber='actual kernel of H[O,A]',
        actual_cover_count=prepared['actual_cover_count'],
        outside_coefficients=[list(row) for row in prepared['outside_caps']],
        count_orbit_sizes=list(prepared['orbit_sizes']),
        source_rowsets_fixed_across_columns=True, all_labelled_rowsets=True,
        common_positive_atom_measure=True, matrix_chronology_preserved=True,
        rs_normalization_per_group='255^-'+str(8*arity), aggregate_input_floor=9*arity*groups,
        floor_used_numerically=False, global_prefix_included=False,
        orbit_use='outside counts and rank5 quiet covers only; no matrix-orbit collapse',
        directed=False, whole_code_certificate=False)
    return answer
