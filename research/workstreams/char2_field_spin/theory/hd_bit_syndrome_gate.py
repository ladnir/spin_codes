"""FLOAT contraction engine for the bit-syndrome outer, with fixed labelled A.

This module loads no regional bank and performs no work on import. The
caller must supply a common positive moment measure for the ORIGINAL
one-state inner and pin its provenance separately. Pure matrix/count-DP
routines are reused from the reviewed quiet helper, not its construction
metadata or post-refresh bank loader.

For each A, quiet masks use universal exact coefficient caps. Nonquiet
outside masks use rank4/5 affine-code coefficient caps. Information-set
inside positions integrate (1+255*t)^groups against the same atom measure;
other inside positions integrate max(1,t)^groups. Thus atoms above1 are
retained safely. Two-way Young combines the COLUMN matrices before the
whole-word power; whole-word group Minkowski then sums all labelled A.
No source/RS/projective/Chernoff-floor prefix is included here.
"""
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import log2

import numpy as np

from hd_bit_syndrome_kernel import information_cover as exact_information_cover, outside_weight_envelopes
from hd_bit_syndrome_outer import CONSTRUCTION, affine_rank, quiet_support_counts
from hd_post_refresh_quiet_rows import (
    combine_row_bounds, log_matmul, log_matrix_power, moment_atoms,
    quiet_column_matrices, quiet_group_values, transition_geometry,
)

ENDPOINTS = ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1))
LAMBDAS = (F(99, 100), F(1, 100))


def _logs(values):
    answer = np.asarray(values, dtype=float)
    if np.any(np.isnan(answer)) or np.any(np.isposinf(answer)):
        raise ValueError('Finite-or-zero logarithms required')
    return answer


def checked_lambdas(values):
    values = tuple(F(value) for value in values)
    if len(values) != 2 or any(value <= 0 for value in values) or sum(values) != 1:
        raise ValueError('Two positive exact Young weights summing to one required')
    return values


def lexicographic_information(source, *, bits=4):
    """Lexicographically smallest quiet information set, by reverse greedy basis."""
    source = tuple(source)
    if (type(bits) is not int or not 1 <= bits <= 4 or not source
            or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < 1 << bits for i in source)):
        raise ValueError('Sorted distinct binary labels required')
    basis_mask, old_rank = 0, 0
    for i in reversed(source):
        trial = basis_mask | (1 << i)
        new_rank = affine_rank(trial, bits)
        if new_rank > old_rank:
            basis_mask, old_rank = trial, new_rank
    return tuple(i for i in source if not basis_mask >> i & 1)


@lru_cache(maxsize=4)
def row_geometry(size=16, arity=8):
    """Every labelled A, its information set, rank, and ordered binary gap words.

    A gap code starts with a sentinel1 and appends0 for constrained inside
    positions or1 for information positions, in chronological order.
    The marked positions are A's complement, not a relabelled orbit.
    """
    if (type(size) is not int or size < 2 or size > 16 or size & (size - 1)
            or type(arity) is not int or not 1 <= arity <= 8 or not 1 <= size - arity <= 8):
        raise ValueError('Power-of-two size2..16 and one to eight positions per side required')
    bits = size.bit_length() - 1
    sources = tuple(combinations(range(size), arity))
    information, ranks, gaps = [], [], []
    for source in sources:
        inside = set(source)
        info = lexicographic_information(source, bits=bits)
        info_set = set(info)
        outside = tuple(i for i in range(size) if i not in inside)
        rank = affine_rank(sum(1 << i for i in outside), bits)
        codes, previous = [], -1
        for following in (*outside, size):
            code = 1
            for i in range(previous + 1, following):
                code = (code << 1) | int(i in info_set)
            codes.append(code)
            previous = following
        information.append(info)
        ranks.append(rank)
        gaps.append(codes)
    ranks = np.asarray(ranks, dtype=np.int8)
    gaps = np.asarray(gaps, dtype=np.int16)
    ranks.setflags(write=False)
    gaps.setflags(write=False)
    return sources, tuple(information), ranks, gaps


@lru_cache(maxsize=4)
def information_covers(size=16, arity=8):
    """Exact c_i for EVERY labelled A, transporting only the quiet matroid.

    Affine index changes preserve information sets of ker(V_A^T). They
    may therefore prepare the inclusion frequencies, but never identify
    different chronological row-set terms in the final word sum.
    """
    rows = row_geometry(size, arity)[0]
    bits = size.bit_length() - 1
    permutations = [tuple(i ^ ((((i >> j) ^ (i >> (j + 1))) & 1)
                                      * ((1 << j) | (1 << (j + 1)))) for i in range(size))
                    for j in range(bits - 1)]
    if bits > 1:
        permutations.append(tuple(i ^ ((i >> 1) & 1) for i in range(size)))
    permutations.append(tuple(i ^ 1 for i in range(size)))
    unseen = set(rows)
    covers = {}
    while unseen:
        source = min(unseen)
        cover = exact_information_cover(source, bits=bits)
        covers[source] = tuple(cover[i] for i in source)
        unseen.remove(source)
        queue = deque([source])
        while queue:
            source = queue.popleft()
            values = covers[source]
            for permutation in permutations:
                moved = sorted((permutation[i], c) for i, c in zip(source, values))
                target = tuple(i for i, _ in moved)
                if target in unseen:
                    covers[target] = tuple(c for _, c in moved)
                    unseen.remove(target)
                    queue.append(target)
    return tuple(covers[source] for source in rows)


@lru_cache(maxsize=4)
def cover_gap_geometry(size=16, arity=8):
    """Pool actual ordered c_i gap words; retain all labelled A and positions."""
    rows = row_geometry(size, arity)[0]
    covers = information_covers(size, arity)
    values = tuple(sorted({c for cover in covers for c in cover}))
    value_index = {c: i for i, c in enumerate(values)}
    words, word_index, parents, labels, codes = [()], {(): 0}, [-1], [-1], []
    for source, cover in zip(rows, covers):
        positions = dict(zip(source, cover))
        outside = tuple(i for i in range(size) if i not in positions)
        previous, row_codes = -1, []
        for following in (*outside, size):
            word, parent = (), 0
            for i in range(previous + 1, following):
                label = value_index[positions[i]]
                word += (label,)
                index = word_index.get(word)
                if index is None:
                    index = len(words)
                    word_index[word] = index
                    words.append(word)
                    parents.append(parent)
                    labels.append(label)
                parent = index
            row_codes.append(parent)
            previous = following
        codes.append(row_codes)
    codes = np.asarray(codes, dtype=np.int32)
    codes.setflags(write=False)
    return covers, values, tuple(parents), tuple(labels), codes


def weighted_cover_matrices(measure, groups, cover_values, *, q=256):
    """Integrate h_c(t)^groups under the SAME common positive atom measure.

    h_c(t)=(1+(q-1)t^(1/c))^c for c>0. For c=0, use max(1,t),
    including when an atom exceeds one. This is a Holder envelope, not
    an assertion that coordinate values in the fiber are independent.
    """
    bank = weighted_gap_matrices(measure, groups, q=q)
    values = tuple(F(c) for c in cover_values)
    if not values or any(c < 0 or c > 1 for c in values):
        raise ValueError('Nonempty exact cover exponents in [0,1] required')
    matrices = np.full((len(values), 2, 2), -np.inf)
    for index, c in enumerate(values):
        atom_weights = (groups * np.maximum(0., bank['log2_atoms']) if not c else
            groups * float(c) * np.logaddexp2(0., log2(q-1) + bank['log2_atoms'] / float(c)))
        for (i, j), mass, weight in zip(ENDPOINTS, bank['log2_masses'], atom_weights):
            matrices[index, i, j] = np.logaddexp2(matrices[index, i, j], mass + weight)
    bank.update(cover_values=values, cover_matrices=matrices)
    return bank


def ordered_cover_products(matrices, parents, labels):
    """Chronological prefix products for pooled gap words; index0 is empty."""
    matrices = _logs(matrices)
    if (matrices.ndim != 3 or matrices.shape[1:] != (2, 2)
            or len(parents) != len(labels) or not parents or parents[0] != -1 or labels[0] != -1):
        raise ValueError('2x2 cover matrices and rooted prefix geometry required')
    result = np.full((len(parents), 2, 2), -np.inf)
    result[0, 0, 0] = result[0, 1, 1] = 0.
    for index in range(1, len(parents)):
        parent, label = parents[index], labels[index]
        if not 0 <= parent < index or not 0 <= label < len(matrices):
            raise ValueError('Invalid chronological prefix dependency')
        result[index] = log_matmul(result[parent], matrices[label])
    return result


def weighted_gap_matrices(measure, groups, *, q=256):
    """EI and ED from the SAME common atoms, including zero and above-one atoms."""
    if type(groups) is not int or groups < 1 or type(q) is not int or q < 2:
        raise ValueError('Positive groups and field-count parameter q>=2 required')
    masses, empty = moment_atoms(measure)
    atoms = np.array([-np.inf, measure['log2_atom00'], measure['log2_atom01'],
                      measure['log2_atom10'], 0.])
    info_weights = groups * np.logaddexp2(0., log2(q - 1) + atoms)
    dropped_weights = groups * np.maximum(0., atoms)
    info, dropped = np.full((2, 2), -np.inf), np.full((2, 2), -np.inf)
    for (i, j), mass, value_i, value_d in zip(ENDPOINTS, masses, info_weights, dropped_weights):
        info[i, j] = np.logaddexp2(info[i, j], mass + value_i)
        dropped[i, j] = np.logaddexp2(dropped[i, j], mass + value_d)
    return dict(information=info, constrained=dropped, empty=empty,
                log2_atoms=atoms, log2_masses=np.asarray(masses))


def ordered_gap_products(information, constrained, maximum=8):
    """All binary gap words of bounded length, retaining noncommuting order."""
    information, constrained = _logs(information), _logs(constrained)
    if (information.shape != (2, 2) or constrained.shape != (2, 2)
            or type(maximum) is not int or not 0 <= maximum <= 8):
        raise ValueError('Two 2x2 matrices and maximum gap0..8 required')
    answer = np.full((1 << (maximum + 1), 2, 2), -np.inf)
    answer[1, 0, 0] = answer[1, 1, 1] = 0.
    for length in range(1, maximum + 1):
        start, stop = 1 << length, 1 << (length + 1)
        parents = answer[np.arange(start, stop) >> 1]
        right = np.stack((constrained, information))[np.arange(start, stop) & 1]
        answer[start:stop] = log_matmul(parents, right)
    return answer


def nonquiet_column_matrices(log_outside_caps, groups, measure, *, size=16, arity=8,
                             q=256, batch_size=128, check_time=lambda: None,
                             information_cover=False):
    """Projected nonquiet columns, BEFORE X^(arity*groups) and Young factors."""
    if (type(groups) is not int or groups < 1 or type(batch_size) is not int
            or not 1 <= batch_size <= 512):
        raise ValueError('Positive groups and explicit batch1..512 required')
    check_time()
    rows, information, ranks, gap_codes = row_geometry(size, arity)
    marked = size - arity
    keys = tuple(sorted(log_outside_caps))
    if not keys or not set(map(int, ranks)).issubset(keys):
        raise ValueError('Missing coefficient envelope for an outside affine rank')
    caps = [_logs(log_outside_caps[key]) for key in keys]
    if any(values.shape != (marked + 1,) or values[0] != -np.inf for values in caps):
        raise ValueError('Nonquiet cap needs one coefficient per weight and zero constant')
    if type(information_cover) is not bool:
        raise ValueError('information_cover must be an explicit bool')
    covers = None
    if information_cover:
        covers, cover_values, parents, labels, gap_codes = cover_gap_geometry(size, arity)
        bank = weighted_cover_matrices(measure, groups, cover_values, q=q)
        products = ordered_cover_products(bank['cover_matrices'], parents, labels)
    else:
        bank = weighted_gap_matrices(measure, groups, q=q)
        products = ordered_gap_products(bank['information'], bank['constrained'], arity)
    counts, layers = transition_geometry(marked)
    values = np.stack([groups * quiet_group_values(row, measure, 0.) for row in caps])
    class_index = np.array([keys.index(int(rank)) for rank in ranks], dtype=np.int8)
    masses = bank['log2_masses']
    result = np.full((len(rows), 2, 2), -np.inf)
    for start in range(0, len(rows), batch_size):
        check_time()
        stop = min(len(rows), start + batch_size)
        codes = gap_codes[start:stop]
        current = products[codes[:, 0]][None, ...]
        for position, (next_size, targets) in enumerate(layers):
            check_time()
            following = np.full((next_size, stop - start, 2, 2), -np.inf)
            after = products[codes[:, position + 1]]
            for label, ((incoming, outgoing), mass) in enumerate(zip(ENDPOINTS, masses)):
                if mass == -np.inf:
                    continue
                term = current[..., incoming, None] + after[None, :, None, outgoing, :] + mass
                target = targets[label]
                following[target] = np.logaddexp2(following[target], term)
            current = following
        group_values = values[class_index[start:stop]].T
        result[start:stop] = np.logaddexp2.reduce(current + group_values[:, :, None, None], axis=0)
    return dict(rowsets=rows, information_sets=information, information_cover_exponents=covers,
        information_cover=information_cover, outside_ranks=ranks,
        column_log2_matrices=result, gap_log2_matrices=bank,
        costs=dict(rowset_count=len(rows), final_atom_count_classes=len(counts),
            coefficient_envelope_ranks=list(keys), batch_size=batch_size,
            maximum_dp_layer_bytes=len(counts) * min(batch_size, len(rows)) * 4 * 8,
            gap_word_slots=len(products)),
        common_positive_measure=True, exact_empty_replacement=False)


def combine_columns(quiet, nonquiet, groups, arity, log2_x, *, lambdas=LAMBDAS):
    quiet, nonquiet = _logs(quiet), _logs(nonquiet)
    if (quiet.shape != nonquiet.shape or quiet.shape[-2:] != (2, 2)
            or type(groups) is not int or groups < 1 or type(arity) is not int or arity < 1
            or not np.isfinite(log2_x) or log2_x < 0):
        raise ValueError('Matched column matrices, positive scope, and finite X>=1 required')
    weights = checked_lambdas(lambdas)
    shifts = np.array([(1 - groups) * log2(weights[0]),
                       (1 - groups) * log2(weights[1]) + arity * groups * log2_x])
    return np.logaddexp2(quiet + shifts[0], nonquiet + shifts[1]), shifts


def contract(log_quiet_caps, log_outside_caps, measure, groups, log2_x, *,
             size=16, columns=16, q=256, lambdas=LAMBDAS, batch_size=128,
             check_time=lambda: None, progress=None, information_cover=False):
    """Generic tested contraction; all labelled row sets are always included."""
    quiet_caps = _logs(log_quiet_caps)
    if (quiet_caps.ndim != 1 or not 2 <= len(quiet_caps) <= 9
            or type(columns) is not int or columns < 1 or not np.isfinite(log2_x) or log2_x < 0):
        raise ValueError('Quiet caps througharity1..8, positive columns and X>=1 required')
    checked_lambdas(lambdas)
    arity = len(quiet_caps) - 1
    check_time()
    quiet = quiet_column_matrices(quiet_caps, groups, measure, log2_x, size=size,
                                  batch_size=batch_size, check_time=check_time)
    if progress:
        progress('quiet_complete')
    nonquiet = nonquiet_column_matrices(log_outside_caps, groups, measure, size=size,
        arity=arity, q=q, batch_size=batch_size, check_time=check_time,
        information_cover=information_cover)
    if progress:
        progress('nonquiet_complete')
    if quiet['rowsets'] != nonquiet['rowsets']:
        raise AssertionError('Branch row-set order differs')
    combined, shifts = combine_columns(quiet['column_log2_matrices'], nonquiet['column_log2_matrices'],
                                        groups, arity, log2_x, lambdas=lambdas)
    check_time()
    word = log_matrix_power(combined, columns)
    values = np.logaddexp2.reduce(word[:, 0, :], axis=1)
    return dict(rowsets=quiet['rowsets'], information_sets=nonquiet['information_sets'],
        information_cover_exponents=nonquiet['information_cover_exponents'],
        information_cover=information_cover,
        outside_ranks=nonquiet['outside_ranks'], per_rowset_word_log2=values,
        combined_word_log2=combine_row_bounds(values, groups), column_log2_matrices=combined,
        branch_column_log2_matrices=dict(quiet=quiet['column_log2_matrices'],
                                        nonquiet=nonquiet['column_log2_matrices']),
        costs=dict(quiet=quiet['costs'], nonquiet=nonquiet['costs']),
        young_log2_scales=shifts, gap_log2_matrices=nonquiet['gap_log2_matrices'],
        terminal_vector=[1, 1], extra_rowset_factor=False,
        common_positive_measure=True, all_labelled_rowsets=True, branch_sum_before_word_power=True)


@lru_cache(maxsize=1)
def coefficient_tables():
    quiet = [0] * 9
    for mask, count in enumerate(quiet_support_counts()):
        weight = mask.bit_count()
        if weight <= 8:
            quiet[weight] = max(quiet[weight], count)
    outside = outside_weight_envelopes()
    return tuple(quiet), outside


def full_row_bounds(measure, groups, log2_x, *, lambdas=LAMBDAS, batch_size=128,
                    check_time=lambda: None, progress=None, information_cover=False):
    """Production-scope coefficients for r groups, arity8, and sixteen columns.

    The caller supplies/pins original-one-state regional curves through r,
    the resulting valid moment measure, and all global prefix factors.
    This raw function does not infer provenance from arbitrary floats.
    """
    check_time()
    quiet, outside = coefficient_tables()
    def logs(values):
        return np.array([log2(value) if value else -np.inf for value in values])
    answer = contract(logs(quiet), {rank: logs(values) for rank, values in outside},
        measure, groups, log2_x, size=16, columns=16, q=256, lambdas=lambdas,
        batch_size=batch_size, check_time=check_time, progress=progress,
        information_cover=information_cover)
    answer['metadata'] = dict(outer_construction=CONSTRUCTION, intended_inner='original-one-state',
        bank_loading='caller-owned; no bank or provider loaded by this engine',
        groups=groups, arity=8, columns=16, field_order=256,
        quiet_coefficients=list(quiet), outside_coefficients={rank: list(values) for rank, values in outside},
        lexicographic_information_sets=not information_cover, information_cover=information_cover,
        cover_preparation='quiet-matroid affine transport only; all labelled A retained',
        source_rowsets_fixed_across_columns=True,
        common_positive_atom_measure=True, atoms_above_one_supported=True,
        young_weights=[str(value) for value in checked_lambdas(lambdas)],
        exact_sparse_restoration=False, global_prefix_included=False,
        directed=False, whole_code_certificate=False)
    return answer
