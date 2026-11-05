"""FLOAT full-fiber contraction for all labelled arity-seven row sets.

Fix A of size7, O=A^c, and the actual outside map C=H*[O,A]. Nine
outside labels make V_O injective, so ker C is the quiet kernel and
the inside difference fiber has dimension7-rank(V_A). Its information
cover may be prepared by affine transport. Outside IMAGE counts cannot
always be transported: when rank(V_A)=4 (240 sets), im C is a proper
subspace of im(V_O). Those240 images are counted separately from actual H.

One full-fiber column upper is C_A(t_O) H_A(t_A), including outside zero.
The input marker is fixed at1. A common positive matrix-atom measure
contracts this bound, then each labelled A retains its Q_A^16 before
whole-word group Minkowski. There is no Young split or global prefix.
No provider or distance certificate is supplied by this FLOAT engine.
"""
from itertools import combinations
from math import comb, log2

import numpy as np

from gf8_packet_maps import INV, MUL
from hd_bit_syndrome_arity9_gate import gap_geometry, quiet_cover_orbits
from hd_bit_syndrome_cross_fibers import cross_fiber
from hd_bit_syndrome_gate import ordered_cover_products, weighted_cover_matrices
from hd_bit_syndrome_outer import CONSTRUCTION, affine_rank
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power


def image_support_counts(generator, *, q=256, max_subsets=512):
    """Exact local-mask counts of an actual linear image, by subset inversion.

    Independent generator vectors are rows. For supported positions S,
    the count is q^(dimension-rank(generator restricted to S^c)).
    All restriction ranks share incremental exact field-echelon bases.
    Non-AES toy fields are allowed only for binary generator matrices.
    """
    generator = tuple(tuple(row) for row in generator)
    if (not generator or not generator[0] or type(q) is not int or q < 2 or q & (q-1)
            or type(max_subsets) is not int or max_subsets < 1):
        raise ValueError('Nonempty independent field generator and explicit subset budget required')
    dimension, length = len(generator), len(generator[0])
    if (dimension > length or (1 << length) > max_subsets
            or any(len(row) != length or any(type(v) is not int or not 0 <= v < 256 for v in row)
                   for row in generator)
            or (q != 256 and any(v not in (0, 1) for row in generator for v in row))):
        raise ValueError('Valid bounded AES generator, or binary toy generator, required')
    columns = tuple(tuple(row[i] for row in generator) for i in range(length))
    bases = [(None,)*dimension] * (1 << length)
    ranks = [0] * (1 << length)
    for mask in range(1, 1 << length):
        bit = mask & -mask
        previous = mask ^ bit
        basis = list(bases[previous])
        vector = list(columns[bit.bit_length()-1])
        rank = ranks[previous]
        for pivot in range(dimension):
            value = vector[pivot]
            if not value:
                continue
            old = basis[pivot]
            if old is None:
                inverse = INV[value]
                basis[pivot] = tuple(MUL[inverse][x] for x in vector)
                rank += 1
                break
            multiples = MUL[value]
            vector = [x ^ multiples[y] for x, y in zip(vector, old)]
        bases[mask], ranks[mask] = tuple(basis), rank
    if ranks[-1] != dimension:
        raise ValueError('Generator vectors must be linearly independent')
    full = (1 << length)-1
    counts = [q**(dimension-ranks[full ^ mask]) for mask in range(1 << length)]
    for bit in range(length):
        step = 1 << bit
        for base in range(0, len(counts), 2*step):
            for offset in range(step):
                counts[base+step+offset] -= counts[base+offset]
    if min(counts) < 0 or counts[0] != 1 or sum(counts) != q**dimension:
        raise AssertionError('Exact actual-image support inversion failed')
    return tuple(counts)


def actual_image_caps(source):
    """Count one labelled C image; no affine covariance of H* is assumed."""
    if len(source) != 7:
        raise ValueError('Seven source positions required')
    fiber = cross_fiber(source)
    counts = image_support_counts(fiber['outside_image_basis'])
    caps = [0]*10
    for mask, count in enumerate(counts):
        caps[mask.bit_count()] = max(caps[mask.bit_count()], count)
    if (fiber['cross_rank'] != fiber['quiet_constraint_rank']
            or not fiber['zero_outside_equals_quiet']):
        raise AssertionError('Arity-seven zero-outside fiber must equal its quiet kernel')
    return tuple(caps), fiber['cross_rank'], counts


def prepare_geometry(*, check_time=lambda: None):
    """Exact generic count orbits plus240 independently counted actual images."""
    prepared = quiet_cover_orbits(16, 7, check_time=check_time)
    rows = prepared['rowsets']
    unique_caps, cap_index, classes, ranks, dimensions, exceptional = [], {}, [], [], [], []
    for index, (source, cover, old_class) in enumerate(zip(rows, prepared['covers'], prepared['classes'])):
        if not index % 64:
            check_time()
        mask = sum(1 << i for i in source)
        rank = affine_rank(mask)
        if affine_rank(((1 << 16)-1) ^ mask) != 5 or rank not in (4, 5):
            raise AssertionError('Seven/nine-position affine ranks changed')
        if rank == 4:
            caps, actual_rank, _ = actual_image_caps(source)
            if actual_rank != rank:
                raise AssertionError('Exceptional actual outside-image dimension changed')
            exceptional.append(index)
        else:
            caps = prepared['outside_caps'][int(old_class)]
        dimension = 7-rank
        if sum(cover) != dimension:
            raise AssertionError('Quiet information-cover dimension changed')
        cls = cap_index.get(caps)
        if cls is None:
            cls = len(unique_caps)
            cap_index[caps] = cls
            unique_caps.append(caps)
        classes.append(cls)
        ranks.append(rank)
        dimensions.append(dimension)
    if len(rows) != comb(16, 7) or len(exceptional) != 240:
        raise AssertionError('Complete arity-seven partition changed')
    prepared.update(outside_caps=tuple(unique_caps), classes=np.asarray(classes, dtype=np.int16),
        outside_ranks=tuple(ranks), fiber_dimensions=tuple(dimensions),
        exceptional_indices=tuple(exceptional),
        cap_class_rule='identical exact weight-cap vectors only; no chronological row merging')
    return prepared


def contract(prepared, measure, groups, *, size=16, columns=16, q=256,
             batch_size=128, check_time=lambda: None):
    """Contract supplied full-fiber bounds; geometry is a mathematical premise."""
    from hd_bit_syndrome_counted_general import counted_column_matrices, quiet_group_values_general
    rows, covers = prepared['rowsets'], prepared['covers']
    caps, classes = prepared['outside_caps'], np.asarray(prepared['classes'])
    if (type(groups) is not int or groups < 1 or type(columns) is not int or columns < 1
            or type(q) is not int or q < 2 or q & (q-1)
            or not rows or classes.shape != (len(rows),)):
        raise ValueError('Positive scope and one class per labelled row required')
    marked = size-len(rows[0])
    if not 1 <= marked <= 9:
        raise ValueError('This full-fiber interface marks one to nine outside positions')
    if (not caps or any(len(row) != marked+1 or row[0] != 1
            or any(type(value) is not int or value < 0 for value in row) for row in caps)):
        raise ValueError('Exact nonnegative complete-image caps with constant1 required')
    check_time()
    values, parents, labels, gap_codes = gap_geometry(rows, covers, size)
    bank = weighted_cover_matrices(measure, groups, values, q=q)
    products = ordered_cover_products(bank['cover_matrices'], parents, labels)
    log_caps = [[log2(value) if value else -np.inf for value in row] for row in caps]
    terminal = np.stack([groups*quiet_group_values_general(row, measure, 0.) for row in log_caps])
    column, costs = counted_column_matrices(terminal, classes, measure, marked,
        products, gap_codes, batch_size=batch_size, check_time=check_time)
    check_time()
    word = log_matrix_power(column, columns)
    per_row = np.logaddexp2.reduce(word[:, 0, :], axis=1)
    return dict(rowsets=tuple(rows), information_cover_exponents=tuple(covers),
        outside_ranks=prepared.get('outside_ranks'), fiber_dimensions=prepared.get('fiber_dimensions'),
        per_rowset_word_log2=per_row, combined_word_log2=combine_row_bounds(per_row, groups),
        column_log2_matrices=column, gap_log2_matrices=bank, costs=costs,
        outside_zero_included=True, young_split=False, extra_rowset_factor=False,
        common_positive_measure=True, terminal_vector=[1, 1], input_marker=1,
        whole_code_certificate=False, directed=False)


def full_row_bounds(measure, groups=18, *, prepared=None, batch_size=128,
                    check_time=lambda: None):
    """Complete raw J_A for a7; caller owns bank authentication and prefix.

    Each homogeneous r-group bound carries255^(-56r); each heterogeneous
    J_A^(1/r) instead carries255^-56. The input floor63r is inert at X=1.
    A supplied prepared object is trusted: construct it here or authenticate
    its complete contents. This engine does not validate external caches.
    """
    check_time()
    prepared = prepare_geometry(check_time=check_time) if prepared is None else prepared
    if prepared['rowsets'] != tuple(combinations(range(16), 7)):
        raise ValueError('All labelled arity-seven rows in lexicographic order required')
    answer = contract(prepared, measure, groups, batch_size=batch_size, check_time=check_time)
    answer['metadata'] = dict(outer_construction=CONSTRUCTION, intended_inner='original-one-state',
        groups=groups, arity=7, columns=16, field_order=256, rowset_count=comb(16, 7),
        outside_marked_positions=9, complete_outside_zero_branch=True,
        exceptional_actual_image_count=len(prepared['exceptional_indices']),
        exceptional_actual_image_dimension=4, generic_image_dimension=5,
        exceptional_fiber_dimension=3, generic_fiber_dimension=2,
        outside_coefficients=[list(row) for row in prepared['outside_caps']],
        count_orbit_sizes=list(prepared['orbit_sizes']),
        orbit_use='quiet covers and generic full-rank outside counts only',
        exceptional_image_counts='actual H[O,A] image, exact support Mobius per labelled exception',
        source_rowsets_fixed_across_columns=True, all_labelled_rowsets=True,
        common_positive_atom_measure=True, atoms_above_one_supported=True,
        input_marker=1, aggregate_input_floor=63*groups, floor_used_numerically=False,
        young_split=False, exact_sparse_restoration=False, global_prefix_included=False,
        rs_normalization_per_group='255^-56', no_extra_rowset_count=True,
        bank_loading='caller-owned; no regional provider or bank loaded',
        directed=False, whole_code_certificate=False)
    return answer
