"""Exact positive column majorants for the new bit-syndrome outer.

For a fixed eight-row set A, retain every quiet word exactly. Nonquiet
fibers are projected to A^c and an information set of the quiet kernel.
The projection polynomial is an upper bound only for output markers in
[0,1], with a common input marker X>=1. evaluate_kernel adds the necessary
constrained-coordinate factors for arbitrary nonnegative finite markers.
Neither interface is a raw regional-matrix inequality.
Singleton and pair input words can be restored exactly by removing their
projected terms with exact integer arithmetic, then adding the true terms.
No probability normalization, RS factor, or row-set multiplicity is here.
"""
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
from itertools import combinations

from hd_bit_syndrome_outer import (
    CONSTRUCTION, affine_rank, matrix, pair_mask_counts, quiet_support_counts,
)


def _source(source, bits):
    source = tuple(source)
    if (type(bits) is not int or not 1 <= bits <= 4 or not source
            or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < 1 << bits for i in source)):
        raise ValueError('Sorted distinct binary-label indices required')
    return source


def _mask(source):
    return sum(1 << i for i in source)


def _field_order(q):
    if type(q) is not int or q < 2 or q & (q - 1):
        raise ValueError('Characteristic-two field order required')


def affine_image_counts(positions, *, q=256, bits=4):
    """Full-support counts of the evaluation code im(V_positions).

    The result uses GLOBAL support masks. For a subset S of the positions,
    the count supported within S is q^(rank(V_positions)-rank(V_positions-S)).
    Exact subset inversion produces nonnegative integer support counts.
    """
    positions = _source(positions, bits)
    _field_order(q)
    rank_all = affine_rank(_mask(positions), bits)
    size = len(positions)
    counts = [q ** (rank_all - affine_rank(_mask(positions[i] for i in range(size)
                                                       if not (mask >> i & 1)), bits))
              for mask in range(1 << size)]
    for i in range(size):
        step = 1 << i
        for base in range(0, len(counts), 2 * step):
            for j in range(step):
                counts[base + step + j] -= counts[base + j]
    if min(counts) < 0 or sum(counts) != q ** rank_all:
        raise AssertionError('Affine evaluation support count identity failed')
    return {_mask(positions[i] for i in range(size) if mask >> i & 1): count
            for mask, count in enumerate(counts) if count}


def information_sets(source, *, bits=4):
    """Information sets of ker(V_A^T), as complements of column bases."""
    source = _source(source, bits)
    r = affine_rank(_mask(source), bits)
    answer = []
    for basis in combinations(source, r):
        if affine_rank(_mask(basis), bits) == r:
            answer.append(tuple(i for i in source if i not in basis))
    return tuple(answer)


def information_cover(source, *, bits=4):
    """Exact inclusion frequencies for uniform Holder over information sets.

    If c_i>0, a translated-fiber bound is
    product_i (1+(q-1)*t_i^(1/c_i))^c_i. Coordinates with c_i=0 are
    constant on the fiber and are dropped using t_i<=1. This function
    returns the rational exponents only; it performs no root arithmetic.
    """
    source = _source(source, bits)
    sets = information_sets(source, bits=bits)
    return {i: Fraction(sum(i in row for row in sets), len(sets)) for i in source}


@lru_cache(maxsize=1)
def affine_eight_orbits():
    """Eight-point subsets under binary affine index changes, for COUNTS only.

    Adjacent bit swaps and one shear generate GL(4,2); conjugating the one
    translation supplies all translations. This symmetry belongs to V's
    evaluation code and quiet ranks. It is NOT a symmetry of chronological
    regional operators or necessarily of H*: do not collapse labelled A
    in a whole-word distance sum using these orbit sizes.
    """
    permutations = []
    for j in range(3):
        permutations.append(tuple(i ^ ((((i >> j) ^ (i >> (j + 1))) & 1) * ((1 << j) | (1 << (j + 1))))
                                  for i in range(16)))
    permutations += [tuple(i ^ ((i >> 1) & 1) for i in range(16)),
                     tuple(i ^ 1 for i in range(16))]
    # Four nibble lookups make the complete finite orbit walk inexpensive.
    tables = tuple(tuple(tuple(_mask(p[4 * block + bit] for bit in range(4) if value >> bit & 1)
                               for value in range(16)) for block in range(4))
                   for p in permutations)
    unseen = {_mask(row) for row in combinations(range(16), 8)}
    result = []
    while unseen:
        representative = min(unseen)
        unseen.remove(representative)
        queue = deque([representative])
        count = 0
        while queue:
            mask = queue.popleft()
            count += 1
            for table in tables:
                target = (table[0][mask & 15] | table[1][mask >> 4 & 15]
                          | table[2][mask >> 8 & 15] | table[3][mask >> 12])
                if target in unseen:
                    unseen.remove(target)
                    queue.append(target)
        result.append((representative, count))
    return tuple(result)


@lru_cache(maxsize=1)
def outside_weight_envelopes():
    """Exact per-support-weight maxima, grouped by outside affine rank4/5.

    If c[r][w] is returned, C_O,nonzero(t) <= sum_w c[r][w] e_w(t_O)
    for EVERY fixed eight-position O of rank r and nonnegative markers.
    Orbit sizes are irrelevant to this universal coefficient envelope;
    they are not inserted as additional codeword counts.
    """
    caps = {4: [0] * 9, 5: [0] * 9}
    for mask, _ in affine_eight_orbits():
        positions = tuple(i for i in range(16) if mask >> i & 1)
        r = affine_rank(mask)
        for support, count in affine_image_counts(positions).items():
            if support:
                weight = support.bit_count()
                caps[r][weight] = max(caps[r][weight], count)
    return tuple((r, tuple(values)) for r, values in sorted(caps.items()))


def _fiber_projection_terms(source, *, information=None, q=256, bits=4):
    """Combinatorial fiber model; caller must establish the cross-rank premise.

    Premise: quiet kernel is ker(V_A^T), and outside image equals im(V_Ac)
    with the same rank. The production A8 wrapper below satisfies this.
    Smaller dimensions are used only by independent exact toy tests.
    """
    source = _source(source, bits)
    _field_order(q)
    n, a = 1 << bits, len(source)
    outside = tuple(i for i in range(n) if i not in source)
    if not outside or affine_rank(_mask(source), bits) != affine_rank(_mask(outside), bits):
        raise ValueError('Equal inside/outside affine ranks required for this model')
    choices = information_sets(source, bits=bits)
    information = choices[0] if information is None else tuple(information)
    if information not in choices:
        raise ValueError('Requested coordinates are not an information set')
    quiet_table = quiet_support_counts(q, bits)
    quiet = {}
    for local in range(1 << a):
        mask = _mask(source[i] for i in range(a) if local >> i & 1)
        if quiet_table[mask]:
            quiet[(local.bit_count(), mask)] = quiet_table[mask]
    outside_counts = affine_image_counts(outside, q=q, bits=bits)
    nonquiet = {}
    for mask, count in outside_counts.items():
        if not mask:
            continue
        for local in range(1 << len(information)):
            projected = mask | _mask(information[j] for j in range(len(information)) if local >> j & 1)
            nonquiet[(a, projected)] = count * (q - 1) ** local.bit_count()
    return quiet, nonquiet, information, outside_counts


def exact_sparse_terms(source):
    """Exact (input weight, complete output mask) counts through input weight2."""
    source = _source(source, 4)
    answer = Counter({(0, 0): 1})
    h = matrix()
    for i in source:
        answer[(1, sum(bool(h[j][i]) << j for j in range(16)))] += 255
    for i, j in combinations(source, 2):
        for mask, count in pair_mask_counts(i, j):
            answer[(2, mask)] += 255 * count
    return dict(answer)


def build_kernel(source, *, information=None, exact_sparse=True):
    """Complete positive A8 majorant, with optional exact sparse restoration.

    Counts include every x in F256^A once at unit markers. Input activity
    is exact for quiet/sparse words, and upper-marked by eight otherwise.
    A remains fixed through every later column; no row-set factor is added.
    """
    source = _source(source, 4)
    if len(source) != 8 or type(exact_sparse) is not bool:
        raise ValueError('This proved production interface requires arity8')
    quiet, nonquiet, information, outside_counts = _fiber_projection_terms(
        source, information=information)
    outside_mask = ((1 << 16) - 1) ^ _mask(source)
    projection_mask = outside_mask | _mask(information)
    sparse = exact_sparse_terms(source) if exact_sparse else {}
    if exact_sparse:
        for (weight, output), count in sparse.items():
            if not weight:
                continue
            if not output & outside_mask:
                raise AssertionError('A nonzero sparse word was unexpectedly quiet')
            key = (8, output & projection_mask)
            remainder = nonquiet.get(key, 0) - count
            if remainder < 0:
                raise AssertionError('Sparse projection exceeds the exact fiber count')
            if remainder:
                nonquiet[key] = remainder
            else:
                nonquiet.pop(key, None)
    terms = Counter(quiet)
    terms.update(nonquiet)
    if exact_sparse:
        terms.update({key: count for key, count in sparse.items() if key[0]})
    if sum(terms.values()) != 256 ** 8 or min(terms.values()) <= 0:
        raise AssertionError('Complete positive column mass changed')
    return dict(terms=dict(terms), quiet_terms=quiet,
        nonquiet_projection_terms=nonquiet, exact_sparse_terms=sparse,
        source=source, information=information,
        outside_evaluation_counts=outside_counts,
        metadata=dict(construction=CONSTRUCTION, field_order=256, arity=8,
            rowset_fixed_across_columns=True, exact_column_count_at_unit_markers=str(256 ** 8),
            output_marker_interval=[0, 1], input_marker_minimum=1,
            nonquiet_input_activity_upper=8, exact_sparse_restored=exact_sparse,
            dropped_inside_positions=tuple(i for i in source if i not in information),
            general_marker_evaluator='evaluate_kernel with constrained-coordinate max(1,t) factors',
            rs_normalization_included=False, rowset_multiplicity_included=False,
            raw_regional_matrix_inequality=False, whole_code_certificate=False))


def _evaluate_nonnegative(terms, markers, x=1):
    markers = tuple(markers)
    if (not markers or any(not value >= 0 or value == float('inf') for value in markers)
            or not x >= 1 or x == float('inf')):
        raise ValueError('Finite nonnegative output markers and common input X>=1 required')
    answer = 0
    for (degree, mask), count in terms.items():
        if count < 0 or degree < 0 or mask < 0 or mask >= 1 << len(markers):
            raise ValueError('Positive term table exceeds the marker scope')
        term = count * x ** degree
        while mask:
            bit = mask & -mask
            term *= markers[bit.bit_length() - 1]
            mask ^= bit
        answer += term
    return answer


def evaluate(terms, markers, x=1):
    """Evaluate the positive projection polynomial in its [0,1] domain."""
    markers = tuple(markers)
    if any(not 0 <= value <= 1 for value in markers):
        raise ValueError('The projection polynomial requires output markers in[0,1]')
    return _evaluate_nonnegative(terms, markers, x)


def evaluate_kernel(kernel, markers, x=1):
    """Sound scalar majorant for ANY nonnegative finite output markers.

    A coordinate omitted from a projected nonquiet mask contributes at
    most max(1,t_i), regardless of its true zero/nonzero status. Keep that
    factor outside the entire remaining nonquiet sum. Quiet and restored
    sparse terms still use their true masks, without that surcharge.
    """
    markers = tuple(markers)
    quiet = _evaluate_nonnegative(kernel['quiet_terms'], markers, x)
    sparse = _evaluate_nonnegative({key: value for key, value in kernel['exact_sparse_terms'].items()
                                   if key[0]}, markers, x)
    projected = _evaluate_nonnegative(kernel['nonquiet_projection_terms'], markers, x)
    factor = 1
    for i in kernel['source']:
        if i not in kernel['information']:
            factor *= max(1, markers[i])
    return quiet + sparse + factor * projected
