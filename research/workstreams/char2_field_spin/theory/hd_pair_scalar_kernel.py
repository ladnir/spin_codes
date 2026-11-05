"""Exact column support counts and a FLOAT one-layer pair contraction.

Fix a labelled allowed input set A.  Its coordinates range over the whole
field, including zero.  For a supplied generator G with m output rows,
n[S] counts INPUT vectors whose output has exact support S.  For each
zero set Z the count is q**(a-rank(G[Z,:])).  Integer subset inversion
recovers n without enumerating field-valued inputs.  Python integers
retain counts larger than int64, and exact zeros survive the inversion.

The pair layer matches row i of the first column with row i xor delta of
the second.  Independently for every pair, alpha is uniform nonzero and
(u,v) maps to (u+alpha*v, u+2*alpha*v).  Its scalar support kernel is
F=[[1,t**2],[t**2,(2*t+(q-3)*t**2)/(q-1)]].  A positive tensor transform
and permuted dot products contract two columns with the SAME fixed A.

The returned polynomial degree marks the number of nonzero ORIGINAL
input columns, not the number of active mixed output columns.  It has
degrees0,1,2.  No input arity normalization, RS floor, row-set sum,
regional operator, second pair layer, or distance certificate is included.
Only the support counting is exact; contractions are ordinary FLOAT.
"""
from functools import lru_cache

import numpy as np

from hd_xor_state_local import field_tables


MAX_ROWS = 16
MAX_BATCH = 8192


def _field(q):
    if type(q) is not int or q not in (4, 8, 16, 256):
        raise ValueError('Pair kernels require GF4, GF8, GF16, or AES GF256')
    return field_tables(q)


def _matrix(generator, q):
    _field(q)
    matrix = np.asarray(generator)
    if (matrix.ndim != 2 or not 1 <= matrix.shape[0] <= MAX_ROWS
            or not 0 <= matrix.shape[1] <= MAX_ROWS
            or not np.issubdtype(matrix.dtype, np.integer)
            or np.any(matrix < 0) or np.any(matrix >= q)):
        raise ValueError('An integer field matrix with1..16 rows and0..16 columns is required')
    return matrix.astype(np.uint8, copy=True)


def _rank_batch(work, table, inverse):
    """The existing vectorized elimination pattern, bounded to16x16 here."""
    count, rows, columns = work.shape
    ranks = np.zeros(count, dtype=np.int64)
    row_indices = np.arange(rows)[None, :]
    for column in range(columns):
        eligible = (row_indices >= ranks[:, None]) & (work[:, :, column] != 0)
        chosen = np.flatnonzero(eligible.any(axis=1))
        if not len(chosen):
            continue
        pivot = np.argmax(eligible[chosen], axis=1)
        target = ranks[chosen]
        old = work[chosen, target, :].copy()
        work[chosen, target, :] = work[chosen, pivot, :]
        work[chosen, pivot, :] = old
        normalized = table[inverse[work[chosen, target, column]][:, None],
                           work[chosen, target, :]]
        factors = work[chosen, :, column].copy()
        factors[np.arange(len(chosen)), target] = 0
        work[chosen] ^= table[factors[:, :, None], normalized[:, None, :]]
        work[chosen, target, :] = normalized
        ranks[chosen] += 1
    return ranks


def zero_set_ranks(generator, *, q=256, batch_size=2048, check_time=None):
    """Exact ranks indexed by the mask of output rows forced to zero.

    At most batch_size dense matrices exist simultaneously.  The callback,
    if supplied, is invoked between batches; it may raise to stop a run.
    """
    matrix = _matrix(generator, q)
    if type(batch_size) is not int or not 1 <= batch_size <= MAX_BATCH:
        raise ValueError('Rank batch size must be1..8192')
    rows, columns = matrix.shape
    table, inverse = _field(q)
    answer = np.empty(1 << rows, dtype=np.uint8)
    shifts = np.arange(rows, dtype=np.uint32)[None, :]
    for start in range(0, len(answer), batch_size):
        if check_time is not None:
            check_time()
        end = min(start+batch_size, len(answer))
        selected = ((np.arange(start, end, dtype=np.uint32)[:, None] >> shifts) & 1) != 0
        work = np.where(selected[:, :, None], matrix[None, :, :], np.uint8(0))
        answer[start:end] = _rank_batch(work, table, inverse)
    return answer


def support_counts_from_matrix(generator, *, q=256, batch_size=2048, check_time=None):
    """Return exact integer counts n[S] for all2**m labelled output masks.

    Generator columns need not be independent.  Counts include kernel
    multiplicity, so sum(n)==q**a and n[0]==q**nullity(G).
    """
    matrix = _matrix(generator, q)
    rows, columns = matrix.shape
    ranks = zero_set_ranks(matrix, q=q, batch_size=batch_size, check_time=check_time)
    powers = tuple(q**dimension for dimension in range(columns+1))
    full = (1 << rows)-1
    counts = np.array([powers[columns-int(ranks[full ^ support])]
                       for support in range(1 << rows)], dtype=object)
    for bit in range(rows):
        if check_time is not None:
            check_time()
        width = 1 << bit
        blocks = counts.reshape(-1, 2, width)
        blocks[:, 1, :] -= blocks[:, 0, :]
    result = tuple(int(value) for value in counts)
    if (min(result) < 0 or sum(result) != q**columns
            or result[0] != powers[columns-int(ranks[full])]):
        raise ArithmeticError('Integer support inversion failed its mass/kernel checks')
    return result


def support_counts(source, *, h=None, q=256, batch_size=2048, check_time=None):
    """Allowed-input counts for a fixed labelled subset of the columns of h.

    The default h is the actual16x16 AES-field H*.  Supplying another h
    does not transfer an H* construction or distance claim to that matrix.
    """
    if h is None:
        if q != 256:
            raise ValueError('The default H-star matrix is defined over AES GF256')
        from hd_bit_syndrome_outer import matrix
        h = matrix()
    h = _matrix(h, q)
    source = tuple(source)
    if (source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < h.shape[1] for i in source)):
        raise ValueError('A sorted labelled subset of generator columns is required')
    return support_counts_from_matrix(h[:, source], q=q, batch_size=batch_size,
                                      check_time=check_time)


def _count_vector(counts):
    counts = tuple(counts)
    length = len(counts)
    if (length < 2 or length & (length-1) or length > (1 << MAX_ROWS)
            or any(type(value) is not int or value < 0 for value in counts)):
        raise ValueError('Nonnegative exact integer counts on all output masks are required')
    return counts, length.bit_length()-1


def pair_kernel(t, *, q=256):
    """FLOAT2x2 support interaction for independent uniform-NZ pair ratios."""
    _field(q)
    t = float(t)
    if not np.isfinite(t) or not 0 <= t <= 1:
        raise ValueError('Scalar mark must be finite and lie in[0,1]')
    square = t*t
    return np.array(((1.0, square), (square, (2*t+(q-3)*square)/(q-1))))


def tensor_transform(values, interaction):
    """Apply interaction**tensor(rows) using only positive2x2 butterflies.

    Mask bit0 is processed first.  The returned vector is indexed by the
    FIRST support mask; the input vector is indexed by the SECOND mask.
    """
    result = np.asarray(values, dtype=float)
    if (result.ndim != 1 or len(result) < 2 or len(result) & (len(result)-1)
            or len(result) > (1 << MAX_ROWS)
            or not np.all(np.isfinite(result)) or np.any(result < 0)):
        raise ValueError('A finite nonnegative mask vector is required')
    interaction = np.asarray(interaction, dtype=float)
    if (interaction.shape != (2, 2) or not np.all(np.isfinite(interaction))
            or np.any(interaction < 0)):
        raise ValueError('A finite nonnegative2x2 interaction is required')
    result = result.copy()
    for bit in range(len(result).bit_length()-1):
        blocks = result.reshape(-1, 2, 1 << bit)
        low, high = blocks[:, 0, :].copy(), blocks[:, 1, :].copy()
        blocks[:, 0, :] = interaction[0, 0]*low+interaction[0, 1]*high
        blocks[:, 1, :] = interaction[1, 0]*low+interaction[1, 1]*high
    if not np.all(np.isfinite(result)):
        raise FloatingPointError('Tensor contraction overflowed')
    return result


@lru_cache(maxsize=32)
def xor_mask_indices(rows, delta):
    """Map every support mask under the label permutation i -> i xor delta."""
    if (type(rows) is not int or not 1 <= rows <= MAX_ROWS or rows & (rows-1)
            or type(delta) is not int or not 0 <= delta < rows):
        raise ValueError('Power-of-two row count and an in-range xor displacement required')
    masks = np.arange(1 << rows, dtype=np.uint32)
    permuted = np.zeros_like(masks)
    for bit in range(rows):
        permuted |= ((masks >> bit) & 1) << (bit ^ delta)
    permuted.flags.writeable = False
    return permuted


def pair_block_coefficients(counts, t, *, deltas=tuple(range(8, 16)), q=256):
    """Return FLOAT degree0/1/2 coefficients, separately for every delta.

    Split off exactly ONE zero input before any conversion to float.
    Nonzero inputs that map to output zero remain in degree1 or2.  This
    also gives correct input activity for a noninjective supplied matrix.
    One transform is reused across all deltas; each degree2 coefficient
    needs one permuted dot product.  Matching cases are not averaged.
    """
    counts, rows = _count_vector(counts)
    if counts[0] < 1:
        raise ValueError('Allowed-input counts must contain the unique zero input')
    deltas = tuple(deltas)
    if not deltas or len(set(deltas)) != len(deltas):
        raise ValueError('A nonempty list of distinct xor displacements is required')
    permutations = tuple(xor_mask_indices(rows, delta) for delta in deltas)
    nonzero = (counts[0]-1,)+counts[1:]
    vector = np.array(nonzero, dtype=float)
    transformed = tensor_transform(vector, pair_kernel(t, q=q))
    result = {}
    for delta, permutation in zip(deltas, permutations):
        coefficients = np.array((1.0, 2*transformed[0],
                                 np.dot(vector[permutation], transformed)))
        if not np.all(np.isfinite(coefficients)) or np.any(coefficients < 0):
            raise FloatingPointError('Invalid FLOAT pair coefficients')
        result[delta] = coefficients
    return dict(coefficients=result, rows=rows, field_order=q, scalar_mark=float(t),
                original_input_count=sum(counts), metadata=dict(
                    count_stage_exact=True, contraction='FLOAT-not-directed',
                    tensor_transforms=1, permuted_dot_products=len(deltas),
                    same_labelled_source_in_both_columns=True,
                    matching_cases_separate=True, arity_resampled=False,
                    activity_degree='number-of-nonzero-original-input-columns',
                    scalar_mark='total-mixed-output-support-weight',
                    rowset_sum_included=False, rs_normalization_included=False,
                    regional_provider=False, whole_code_certificate=False))
