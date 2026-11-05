"""Exact input/output support counts for at most four binary source columns.

The caller supplies the ordered binary row forms of its actual matrix.
Source columns, input permutations, and output permutations are not averaged.
For x in GF(q)^a, phi_x(v)=sum_i v_i*x_i is an F2-linear map. Its kernel
determines every original-input and output zero. Kernel multiplicities count
injections of the quotient into GF(q), regarded as an F2 vector space.

This is raw allowed-input counting, with no probability, RS, projective,
row-choice, or setup normalization. It is not a regional provider.
"""
from collections import Counter
from functools import lru_cache

MAX_ARITY = 4


def _arity(arity):
    if type(arity) is not int or not 0 <= arity <= MAX_ARITY:
        raise ValueError('The bounded interface requires arity0..4')


def _field(q):
    if type(q) is not int or q < 2 or q & (q-1):
        raise ValueError('A characteristic-two field order is required')


@lru_cache(maxsize=MAX_ARITY+1)
def kernel_subspaces(arity):
    """Canonical subspace bitsets; bit v is1 iff vector v belongs to K."""
    _arity(arity)
    subspaces = {1}  # Only the zero vector.
    pending = [1]
    while pending:
        kernel = pending.pop()
        members = tuple(v for v in range(1 << arity) if kernel >> v & 1)
        for extra in range(1 << arity):
            if kernel >> extra & 1:
                continue
            enlarged = kernel | sum(1 << (v ^ extra) for v in members)
            if enlarged not in subspaces:
                subspaces.add(enlarged)
                pending.append(enlarged)
    return tuple(sorted(subspaces, key=lambda value: (value.bit_count(), value)))


def injection_count(dimension, q=256):
    """Number of injective F2-linear maps from F2^dimension into GF(q)."""
    _field(q)
    if type(dimension) is not int or dimension < 0:
        raise ValueError('A nonnegative quotient dimension is required')
    if dimension > q.bit_length()-1:
        return 0
    count = 1
    for h in range(dimension):
        count *= q-(1 << h)
    return count


def joint_counts(row_forms, arity, *, q=256):
    """Return exact (original_input_weight, ordered_output_mask)->count.

    row_forms[j] is the a-bit binary linear form for output coordinate j.
    The ordered source i contributes input coordinate e_i. In particular,
    a caller may form these rows after any specified P_j/Q_j permutation,
    while retaining the same original source order across all columns.
    Arbitrary binary matrices, including singular ones, are supported.
    """
    _arity(arity)
    _field(q)
    row_forms = tuple(row_forms)
    if any(type(form) is not int or not 0 <= form < 1 << arity for form in row_forms):
        raise ValueError('Every output row must be an in-range binary form')
    result = Counter()
    for kernel in kernel_subspaces(arity):
        size = kernel.bit_count()
        kernel_dimension = size.bit_length()-1
        count = injection_count(arity-kernel_dimension, q)
        if not count:
            continue
        input_weight = sum(not (kernel >> (1 << i) & 1) for i in range(arity))
        output_mask = sum((not (kernel >> form & 1)) << j for j, form in enumerate(row_forms))
        result[input_weight, output_mask] += count
    if sum(result.values()) != q**arity:
        raise AssertionError('Kernel partition lost input arrays')
    return dict(sorted(result.items()))


def joint_counts_from_matrix(matrix, *, q=256):
    """Convenience wrapper for a nonempty ordered rectangular binary matrix."""
    matrix = tuple(tuple(row) for row in matrix)
    if not matrix:
        raise ValueError('Supply at least one output row, or use joint_counts')
    arity = len(matrix[0])
    _arity(arity)
    if any(len(row) != arity or any(type(x) is not int or x not in (0, 1) for x in row)
           for row in matrix):
        raise ValueError('A rectangular matrix with binary integer entries is required')
    forms = tuple(sum(value << i for i, value in enumerate(row)) for row in matrix)
    return joint_counts(forms, arity, q=q)


def select_columns(matrix, sources):
    """Preserve caller-specified source order when restricting an actual matrix."""
    matrix = tuple(tuple(row) for row in matrix)
    sources = tuple(sources)
    if not matrix or any(len(row) != len(matrix[0]) for row in matrix):
        raise ValueError('A nonempty rectangular matrix is required')
    _arity(len(sources))
    if (len(set(sources)) != len(sources)
            or any(type(i) is not int or not 0 <= i < len(matrix[0]) for i in sources)):
        raise ValueError('Distinct in-range source indices are required')
    return tuple(tuple(row[i] for i in sources) for row in matrix)
