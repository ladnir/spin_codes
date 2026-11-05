"""FLOAT full-fiber contraction for every labelled arity-nine row set.

Fix A of size9 and O=A^c. The actual cross map C=H*[O,A] has
image im(V_O): nine binary affine labels span all five syndrome rows.
Let D=H*[A,A] ker(C). Every outside value has an affine inside fiber
with difference space D. A common information-set cover bounds every
translate by H_A(t_A), including atoms above1 via max(1,t) when c_i=0.
The complete column polynomial is bounded by C_A(t_O) H_A(t_A).
It includes outside zero; it does NOT identify that fiber with quiet.

This first interface fixes X=1, so it has no Young split or input-marker
factor. Every labelled A retains its ordered gaps, Q_A^16, and position
in the group Minkowski sum. No RS/source/projective prefix is included.
No provider is loaded. All numerical results are FLOAT proposals only.
"""
from collections import deque
from fractions import Fraction as F
from itertools import combinations
from math import comb, log2

import numpy as np

from hd_bit_syndrome_cross_fibers import cross_fiber, information_sets
from hd_bit_syndrome_gate import ordered_cover_products, weighted_cover_matrices
from hd_bit_syndrome_kernel import affine_image_counts, information_cover
from hd_bit_syndrome_orbit_gate import counted_column_matrices
from hd_bit_syndrome_outer import CONSTRUCTION, affine_rank
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power, quiet_group_values


def _permutations(bits):
    size = 1 << bits
    answer = [tuple(i ^ ((((i >> j) ^ (i >> (j+1))) & 1)
                          * ((1 << j) | (1 << (j+1)))) for i in range(size))
              for j in range(bits-1)]
    if bits > 1:
        answer.append(tuple(i ^ ((i >> 1) & 1) for i in range(size)))
    answer.append(tuple(i ^ 1 for i in range(size)))
    return tuple(answer)


def quiet_cover_orbits(size=16, arity=9, *, q=256, check_time=lambda: None):
    """Prepare affine COUNT orbits and quiet covers, never whole-word orbits.

    For production these covers are only provisional: rank4 outside sets
    are overridden by actual cross-fiber covers in prepare_geometry().
    The small binary-label variants are exact combinatorial test fixtures.
    """
    if (type(size) is not int or not 2 <= size <= 16 or size & (size-1)
            or type(arity) is not int or not 1 <= arity < size
            or type(q) is not int or q < 2 or q & (q-1)):
        raise ValueError('Binary label set, nonempty sides, and characteristic-two q required')
    bits = size.bit_length()-1
    rows = tuple(combinations(range(size), arity))
    unseen, covers, classes = set(rows), {}, {}
    representatives, class_sizes, caps = [], [], []
    permutations = _permutations(bits)
    while unseen:
        check_time()
        representative = min(unseen)
        cls = len(representatives)
        representatives.append(representative)
        cover = information_cover(representative, bits=bits)
        covers[representative] = tuple(cover[i] for i in representative)
        classes[representative] = cls
        unseen.remove(representative)
        queue, count = deque([representative]), 0
        while queue:
            if not count % 256:
                check_time()
            source = queue.popleft()
            count += 1
            for permutation in permutations:
                moved = sorted((permutation[i], c) for i, c in zip(source, covers[source]))
                target = tuple(i for i, _ in moved)
                if target in unseen:
                    covers[target] = tuple(c for _, c in moved)
                    classes[target] = cls
                    unseen.remove(target)
                    queue.append(target)
        class_sizes.append(count)
        outside = tuple(i for i in range(size) if i not in representative)
        values = [0] * (len(outside)+1)
        for support, multiplicity in affine_image_counts(outside, q=q, bits=bits).items():
            weight = support.bit_count()
            values[weight] = max(values[weight], multiplicity)
        if values[0] != 1:
            raise AssertionError('The complete outside image must contain its zero word once')
        caps.append(tuple(values))
    return dict(rowsets=rows, covers=tuple(covers[row] for row in rows),
                classes=np.array([classes[row] for row in rows], dtype=np.int16),
                outside_caps=tuple(caps), representatives=tuple(representatives),
                orbit_sizes=tuple(class_sizes))


def actual_fiber_cover(source):
    """Exact uniform information-set cover of the ACTUAL inside difference code."""
    fiber = cross_fiber(source)
    choices = information_sets(fiber)
    cover = tuple(F(sum(i in choice for choice in choices), len(choices)) for i in source)
    if sum(cover) != fiber['inside_zero_fiber_dimension']:
        raise AssertionError('Information-cover dimension changed')
    return cover, fiber['cross_rank'], fiber['inside_zero_fiber_dimension']


def source_cover(source, quiet_cover):
    """Use transported quiet data only when the actual cross kernel is quiet."""
    if (len(source) != 9 or tuple(source) != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < 16 for i in source)
            or len(quiet_cover) != 9 or any(not 0 <= c <= 1 for c in quiet_cover)
            or affine_rank(sum(1 << i for i in source)) != 5):
        raise ValueError('Nine distinct source labels spanning five affine rows required')
    rank = affine_rank(((1 << 16)-1) ^ sum(1 << i for i in source))
    if rank == 4:
        cover, actual_rank, dimension = actual_fiber_cover(source)
        if actual_rank != rank or dimension != 5:
            raise AssertionError('Exceptional actual fiber dimensions changed')
        return cover, rank, dimension, True
    if rank != 5 or sum(quiet_cover) != 4:
        raise AssertionError('Generic arity-nine quiet cover has the wrong dimension')
    return tuple(quiet_cover), rank, 4, False


def gap_geometry(rows, covers, size):
    """Pool only identical ordered cover words; never identify row-set terms."""
    rows, covers = tuple(map(tuple, rows)), tuple(tuple(F(c) for c in row) for row in covers)
    if (not rows or len(rows) != len(covers) or type(size) is not int or size < 2
            or not 1 <= len(rows[0]) < size):
        raise ValueError('Nonempty matched labelled rows and covers required')
    arity = len(rows[0])
    if any(len(row) != arity or row != tuple(sorted(set(row)))
           or any(type(i) is not int or not 0 <= i < size for i in row)
           or len(cover) != arity or any(not 0 <= c <= 1 for c in cover)
           for row, cover in zip(rows, covers)):
        raise ValueError('Sorted row labels and exact cover exponents in[0,1] required')
    values = tuple(sorted({c for cover in covers for c in cover}))
    indices = {c: index for index, c in enumerate(values)}
    words, parents, labels, codes = {(): 0}, [-1], [-1], []
    for source, cover in zip(rows, covers):
        positions = dict(zip(source, cover))
        outside = tuple(i for i in range(size) if i not in positions)
        previous, row_codes = -1, []
        for following in (*outside, size):
            word, parent = (), 0
            for i in range(previous+1, following):
                label = indices[positions[i]]
                word += (label,)
                index = words.get(word)
                if index is None:
                    index = len(words)
                    words[word] = index
                    parents.append(parent)
                    labels.append(label)
                parent = index
            row_codes.append(parent)
            previous = following
        codes.append(row_codes)
    return values, tuple(parents), tuple(labels), np.asarray(codes, dtype=np.int32)


def prepare_geometry(*, check_time=lambda: None):
    """All11,440 A; affine quiet transport plus240 actual exceptional covers.

    For |A|9, rank(V_A)=5, hence im(C)=im(V_O). If rank(V_O)=5,
    D is the quiet code. The remaining240 sets have dimension5 fibers,
    while their quiet dimension is4; no affine transport is used for D.
    """
    prepared = quiet_cover_orbits(check_time=check_time)
    rows = prepared['rowsets']
    covers = list(prepared['covers'])
    ranks, dimensions, exceptional = [], [], []
    for index, source in enumerate(rows):
        if not index % 128:
            check_time()
        covers[index], rank, dimension, is_exceptional = source_cover(source, covers[index])
        if is_exceptional:
            exceptional.append(index)
        ranks.append(rank)
        dimensions.append(dimension)
    if len(rows) != comb(16, 9) or len(exceptional) != 240:
        raise AssertionError('The complete arity-nine partition changed')
    prepared.update(covers=tuple(covers), outside_ranks=tuple(ranks),
                    fiber_dimensions=tuple(dimensions), exceptional_indices=tuple(exceptional))
    return prepared


def contract(prepared, measure, groups, *, size=16, columns=16, q=256,
             batch_size=128, check_time=lambda: None):
    """One positive full-fiber pass; prepared caps/covers are caller premises.

    The small-test interface permits selected rows. full_row_bounds()
    supplies the complete production partition. No Young/prefix is used.
    """
    rows, covers = prepared['rowsets'], prepared['covers']
    caps, classes = prepared['outside_caps'], np.asarray(prepared['classes'])
    if (type(groups) is not int or groups < 1 or type(columns) is not int or columns < 1
            or type(q) is not int or q < 2 or q & (q-1)
            or not rows or classes.shape != (len(rows),)):
        raise ValueError('Positive scope and one class per labelled row required')
    marked = size-len(rows[0])
    if not 1 <= marked <= 8:
        raise ValueError('This count-DP interface marks one to eight outside positions')
    if (not caps or any(len(row) != marked+1 or row[0] != 1
            or any(type(value) is not int or value < 0 for value in row) for row in caps)):
        raise ValueError('Exact nonnegative complete-image caps with constant1 required')
    check_time()
    values, parents, labels, gap_codes = gap_geometry(rows, covers, size)
    bank = weighted_cover_matrices(measure, groups, values, q=q)
    products = ordered_cover_products(bank['cover_matrices'], parents, labels)
    log_caps = [[log2(value) if value else -np.inf for value in row] for row in caps]
    terminal = np.stack([groups*quiet_group_values(row, measure, 0.) for row in log_caps])
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
    """Complete arity9, original-inner raw J_A; no bank load or global prefix.

    A homogeneous r-group bound subsequently multiplies by255^(-72r).
    In a heterogeneous norm, each J_A^(1/r) instead carries255^-72.
    The aggregate input floor81r is inert at this first control's X=1.
    A supplied prepared object is a trusted mathematical input: a replay
    wrapper must construct it here or authenticate the complete geometry.
    """
    check_time()
    prepared = prepare_geometry(check_time=check_time) if prepared is None else prepared
    expected = tuple(combinations(range(16), 9))
    if prepared['rowsets'] != expected:
        raise ValueError('All labelled arity-nine row sets in lexicographic order required')
    answer = contract(prepared, measure, groups, batch_size=batch_size, check_time=check_time)
    answer['metadata'] = dict(outer_construction=CONSTRUCTION, intended_inner='original-one-state',
        groups=groups, arity=9, columns=16, field_order=256, rowset_count=comb(16, 9),
        outside_marked_positions=7, complete_outside_zero_branch=True,
        exceptional_actual_fiber_count=len(prepared['exceptional_indices']),
        exceptional_actual_fiber_dimension=5, generic_fiber_dimension=4,
        outside_coefficients=[list(row) for row in prepared['outside_caps']],
        count_orbit_sizes=list(prepared['orbit_sizes']),
        orbit_use='outside coefficient counts and generic quiet covers only',
        source_rowsets_fixed_across_columns=True, all_labelled_rowsets=True,
        common_positive_atom_measure=True, atoms_above_one_supported=True,
        input_marker=1, aggregate_input_floor=81*groups, floor_used_numerically=False,
        young_split=False, exact_sparse_restoration=False, global_prefix_included=False,
        rs_normalization_per_group='255^-72', no_extra_rowset_count=True,
        bank_loading='caller-owned; no regional provider or bank loaded',
        directed=False, whole_code_certificate=False)
    return answer
