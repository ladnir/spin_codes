"""Exact GF4 relation categories for at most four allowed source coordinates.

Each packed row uses two bits per coefficient, in caller-supplied source
order. Symbols 0,1,2,3 describe GF4=F2[t]/(t^2+t+1), NOT AES bytes2,3.
The explicit embedding into AES-GF256 is (0,1,0xbc,0xbd).

For x in GF(q)^a, phi_x(v)=sum_i v_i*x_i is GF4-linear. Its kernel K
determines original-input and output zeros. If codim(K)=d, exactly
prod_{h<d}(q-4^h) arrays have that kernel. Counts are raw ambient masses:
there is no RS, projective, row-choice, or setup normalization here.

Categories remain unmerged and identically ordered for every matrix/column.
That stability does not assert that a message has one common K across its
columns. The caller retains its fixed original source set and chronology.
No controller, regional provider, native implementation, or certificate.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product

MAX_ARITY = 4
GF4_MODULUS = 0x7
AES_MODULUS = 0x11b
AES_EMBEDDING = (0, 1, 0xbc, 0xbd)
_AES_DECODE = {value: i for i, value in enumerate(AES_EMBEDDING)}
_MUL = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))


def _arity(arity):
    if type(arity) is not int or not 0 <= arity <= MAX_ARITY:
        raise ValueError('The bounded interface requires arity0..4')


def _field(q):
    if type(q) is not int or q < 4 or q & (q-1) or (q.bit_length()-1) % 2:
        raise ValueError('The ambient field order must be a positive power of4')


def multiply_gf4(a, b):
    if type(a) is not int or type(b) is not int or not 0 <= a < 4 or not 0 <= b < 4:
        raise ValueError('GF4 symbols must lie in0..3')
    return _MUL[a][b]


def embed_aes(value):
    if type(value) is not int or not 0 <= value < 4:
        raise ValueError('A GF4 symbol0..3 is required')
    return AES_EMBEDDING[value]


def pack_row(coefficients):
    coefficients = tuple(coefficients)
    _arity(len(coefficients))
    if any(type(value) is not int or not 0 <= value < 4 for value in coefficients):
        raise ValueError('Rows must contain GF4 symbols0..3')
    return sum(value << (2*i) for i, value in enumerate(coefficients))


def _scale_row(row, scalar, arity):
    return sum(_MUL[scalar][row >> (2*i) & 3] << (2*i) for i in range(arity))


@lru_cache(maxsize=MAX_ARITY+1)
def kernel_subspaces(arity):
    """Canonical RREF bases in increasing (dimension,packed-row-tuple) order.

    Every pivot is1, all other pivot-column entries are0, and free entries
    occur only strictly to the right of the row's pivot. Each subspace is
    therefore generated once, without numerical rank or a deduplication guess.
    """
    _arity(arity)
    result = []
    for dimension in range(arity+1):
        for pivots in combinations(range(arity), dimension):
            free = tuple((i, j) for i, pivot in enumerate(pivots)
                         for j in range(pivot+1, arity) if j not in pivots)
            for entries in product(range(4), repeat=len(free)):
                rows = [1 << (2*pivot) for pivot in pivots]
                for (i, j), value in zip(free, entries):
                    rows[i] |= value << (2*j)
                result.append(tuple(rows))
    return tuple(sorted(result, key=lambda rows: (len(rows), rows)))


@lru_cache(maxsize=MAX_ARITY+1)
def kernel_bitsets(arity):
    """Parallel to kernel_subspaces; bit v records membership of packed v."""
    _arity(arity)
    result = []
    for rows in kernel_subspaces(arity):
        members = (0,)
        for row in rows:
            scaled = tuple(_scale_row(row, scalar, arity) for scalar in range(4))
            members = tuple(value ^ addition for value in members for addition in scaled)
        if len(set(members)) != 4**len(rows):
            raise AssertionError('Generated RREF rows were not independent')
        result.append(sum(1 << value for value in members))
    return tuple(result)


def injection_count(dimension, q=256):
    """Number of GF4-linear quotient injections into the ambient GF(q)."""
    _field(q)
    if type(dimension) is not int or dimension < 0:
        raise ValueError('A nonnegative quotient dimension is required')
    if dimension > (q.bit_length()-1)//2:
        return 0
    count = 1
    for h in range(dimension):
        count *= q-4**h
    return count


@lru_cache(maxsize=32)
def _category_tuples(arity, q):
    _arity(arity)
    _field(q)
    input_masks, weights, counts = [], [], []
    for rows, kernel in zip(kernel_subspaces(arity), kernel_bitsets(arity)):
        mask = sum((not (kernel >> (1 << (2*i)) & 1)) << i for i in range(arity))
        input_masks.append(mask)
        weights.append(mask.bit_count())
        counts.append(injection_count(arity-len(rows), q))
    if sum(counts) != q**arity:
        raise AssertionError('The kernel partition lost ambient input arrays')
    return tuple(input_masks), tuple(weights), tuple(counts)


def category_data(arity, q=256):
    """Stable unmerged category metadata; zero-count categories are retained."""
    masks, weights, counts = _category_tuples(arity, q)
    return dict(arity=arity, q=q, coefficient_field_order=4,
                category_order='increasing-dimension-then-canonical-RREF',
                kernel_rrefs=kernel_subspaces(arity), kernel_bitsets=kernel_bitsets(arity),
                category_input_masks=masks, category_input_weights=weights,
                category_counts=counts, raw_allowed_input_mass=q**arity,
                aes_embedding=AES_EMBEDDING if q == 256 else None)


def column_masks(row_forms, arity):
    """One ordered output mask per stable kernel category, independent of q."""
    _arity(arity)
    row_forms = tuple(row_forms)
    if any(type(form) is not int or not 0 <= form < 4**arity for form in row_forms):
        raise ValueError('Every packed GF4 row form must be in range')
    return tuple(sum((not (kernel >> form & 1)) << j for j, form in enumerate(row_forms))
                 for kernel in kernel_bitsets(arity))


def column_categories(row_forms, arity, *, q=256):
    """Unmerged (original_input_weight,output_mask,count) in stable order."""
    _, weights, counts = _category_tuples(arity, q)
    masks = column_masks(row_forms, arity)
    return tuple(zip(weights, masks, counts))


def joint_counts(row_forms, arity, *, q=256):
    """Merge exact raw counts only after the original input weight is known."""
    result = Counter()
    for weight, mask, count in column_categories(row_forms, arity, q=q):
        if count:
            result[weight, mask] += count
    return dict(sorted(result.items()))


def row_forms_from_matrix(matrix, *, encoding='gf4'):
    """Pack an ordered rectangular matrix; AES decoding must be explicit.

    encoding='gf4' accepts symbols0..3. encoding='aes' accepts ONLY bytes
    0,1,0xbc,0xbd; AES lambda2 is not in this subfield and is rejected.
    Neither input columns nor output rows are reordered.
    """
    matrix = tuple(tuple(row) for row in matrix)
    if not matrix or encoding not in ('gf4', 'aes'):
        raise ValueError('A nonempty matrix and explicit supported encoding are required')
    arity = len(matrix[0])
    _arity(arity)
    if any(len(row) != arity for row in matrix):
        raise ValueError('A rectangular coefficient matrix is required')
    if encoding == 'aes':
        if any(type(value) is not int or value not in _AES_DECODE for row in matrix for value in row):
            raise ValueError('An AES coefficient lies outside the explicit GF4 subfield')
        matrix = tuple(tuple(_AES_DECODE[value] for value in row) for row in matrix)
    return tuple(pack_row(row) for row in matrix)


def joint_counts_from_matrix(matrix, *, q=256, encoding='gf4'):
    matrix = tuple(tuple(row) for row in matrix)
    row_forms = row_forms_from_matrix(matrix, encoding=encoding)
    return joint_counts(row_forms, len(matrix[0]), q=q)
