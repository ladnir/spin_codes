"""Exact support counts and Renyi density moments for dense Hlambda columns.

For mapped complement O of size b<=4, use the FULL parity matrix H[O,:].
Its binary part is M[o,i]=1[o=i]+dot_F2(o,i), with H=M+lambda*J.
The identity contribution is essential; the older cross-rank helper is not
used. For lambda outside F2, rank(M_S+lambda*J)=max(rank M_S,rank(M_S+J)).

For every allowed support S, the code has q**(|S|-rank(H[:,S])) words
supported inside S. Exact integer Mobius inversion yields n_C(S). No field
amplitudes are enumerated and no signed FLOAT transform is used.

After independent nonzero output diagonals, the normalized code measure
has likelihood L(y)=q**b*n_C(supp y)/(q-1)**|supp y| relative to uniform
ambient y. Thus I_p=E L**p equals
q**(p*b-n)*sum_S n_C(S)**p/(q-1)**((p-1)*|S|).
The result is an exact Fraction, grouped into only n+1 rational terms.
There is no RS, source, rowset, regional operator, or certificate factor.
"""
from fractions import Fraction
from functools import lru_cache
from math import log1p, log, log2
from time import perf_counter

from hd_binary_outer_relations import kernel_subspaces

MAX_WIDTH, MAX_CHECKS = 16, 4


def _field(q):
    if type(q) is not int or q not in (4, 256):
        raise ValueError('Only GF4 toy and AES-GF256 counts are bound to this helper')


@lru_cache(maxsize=MAX_CHECKS+1)
def _span_transitions(checks):
    """Positive binary column-span transitions; at most67 states for b=4."""
    spaces = kernel_subspaces(checks)
    lookup = {space: i for i, space in enumerate(spaces)}
    transitions = []
    for vector in range(1 << checks):
        rows = []
        for space in spaces:
            extension = space | sum(1 << (v ^ vector) for v in range(1 << checks) if space >> v & 1)
            rows.append(lookup[extension])
        transitions.append(bytes(rows))
    dimensions = bytes(space.bit_count().bit_length()-1 for space in spaces)
    return tuple(transitions), dimensions, lookup[1]


def subset_ranks(binary_rows, width):
    """All exact ranks of (M+lambda*J)[:,S], including empty S.

    Rows are width-bit packed BINARY M rows, including any identity entries.
    Two <=67-state binary-span DPs share the same subset order. Returned
    bytes are indexed by the actual labelled support mask S.
    """
    rows = tuple(binary_rows)
    if (type(width) is not int or not 1 <= width <= MAX_WIDTH or len(rows) > MAX_CHECKS
            or any(type(row) is not int or not 0 <= row < 1 << width for row in rows)):
        raise ValueError('At most four binary rows and width1..16 are required')
    checks = len(rows)
    transitions, dimensions, zero = _span_transitions(checks)
    columns = tuple(sum(((row >> j) & 1) << i for i, row in enumerate(rows)) for j in range(width))
    ones = (1 << checks)-1
    low = tuple(transitions[column] for column in columns)
    high = tuple(transitions[column ^ ones] for column in columns)
    first, second, ranks = bytearray(1 << width), bytearray(1 << width), bytearray(1 << width)
    first[0] = second[0] = zero
    for mask in range(1, 1 << width):
        bit = mask & -mask
        j, previous = bit.bit_length()-1, mask ^ bit
        a, b = low[j][first[previous]], high[j][second[previous]]
        first[mask], second[mask] = a, b
        ranks[mask] = max(dimensions[a], dimensions[b])
    return bytes(ranks)


def support_counts_from_binary(binary_rows, width, *, q=256):
    """Return labelled exact counts in ker(M+lambda*J), with its true rank."""
    _field(q)
    ranks = subset_ranks(binary_rows, width)
    counts = [q**(mask.bit_count()-rank) for mask, rank in enumerate(ranks)]
    for j in range(width):
        step = 1 << j
        for start in range(0, 1 << width, step*2):
            for offset in range(step):
                counts[start+step+offset] -= counts[start+offset]
    codimension = ranks[-1]
    if counts[0] != 1 or min(counts) < 0 or sum(counts) != q**(width-codimension):
        raise ArithmeticError('Integer support inversion failed its exact mass/positivity checks')
    return dict(width=width, q=q, codimension=codimension,
                subset_ranks=ranks, counts=tuple(counts), mass=q**(width-codimension))


def parity_binary_rows(complement):
    """Full Hlambda parity rows, BEFORE the physical output permutation Q."""
    complement = tuple(complement)
    if (len(complement) > MAX_CHECKS or len(set(complement)) != len(complement)
            or any(type(i) is not int or not 0 <= i < 16 for i in complement)):
        raise ValueError('At most four distinct mapped complement labels are required')
    return tuple(sum((int(o == i) ^ ((o & i).bit_count() & 1)) << i for i in range(16))
                 for o in complement)


def support_counts(complement, *, q=256):
    """Arity12..16 counts. Input complement is already mapped through P_j.

    Q_j permutes the resulting masks and does not change I_p. For a uniform
    maximum over mapped complements there is no additional schedule factor.
    This does not identify individual labelled column distributions.
    """
    complement = tuple(complement)
    result = support_counts_from_binary(parity_binary_rows(complement), 16, q=q)
    if result['codimension'] != len(complement):
        raise ArithmeticError('Full involutive Hlambda rows lost independence')
    return dict(result, complement=complement, arity=16-len(complement),
                lambda_aes_byte=0xbc, includes_diagonal_identity=True)


def density_moment_from_counts(counts, codimension, *, q=256, order=3):
    """Exact E_U[L**order]; raw counts and uniform-space normalization once."""
    _field(q)
    counts = tuple(counts)
    if (not counts or len(counts) & (len(counts)-1) or not 2 <= len(counts) <= 1 << MAX_WIDTH
            or any(type(n) is not int or n < 0 for n in counts)):
        raise ValueError('A power-of-two labelled table of nonnegative Python integer counts is required')
    width = len(counts).bit_length()-1
    if (type(codimension) is not int or not 0 <= codimension <= width
            or type(order) is not int or not 1 <= order <= 4
            or counts[0] != 1 or sum(counts) != q**(width-codimension)):
        raise ValueError('Incompatible exact code mass, codimension, or moment order')
    powers = [0]*(width+1)
    for mask, n in enumerate(counts):
        powers[mask.bit_count()] += n**order
    value = sum((Fraction(n, (q-1)**((order-1)*weight)) for weight, n in enumerate(powers)), Fraction())
    exponent = order*codimension-width
    value *= Fraction(q**exponent) if exponent >= 0 else Fraction(1, q**(-exponent))
    if value < 1 or (order == 1 and value != 1):
        raise ArithmeticError('Normalized density moment violated its exact probability check')
    return value


def renyi_moment(complement, *, q=256, order=3):
    """Exact rational moment plus support diagnostics for one mapped set."""
    start = perf_counter()
    complement = tuple(complement)
    data = support_counts(complement, q=q)
    value = density_moment_from_counts(data['counts'], data['codimension'], q=q, order=order)
    histogram = [0]*17
    nonzero_supports = [0]*17
    for mask, count in enumerate(data['counts']):
        weight = mask.bit_count()
        histogram[weight] += count
        nonzero_supports[weight] += bool(count)
    excess = value-1
    # Convenience FLOAT diagnostics only; the Fraction itself is exact.
    log_value = log1p(float(excess))/log(2.) if excess < 1 else log2(value.numerator)-log2(value.denominator)
    return dict(complement=tuple(complement), arity=data['arity'], order=order, q=q,
        moment=value, excess=excess, log2_moment_float=log_value,
        log2_excess_float=(log2(excess.numerator)-log2(excess.denominator) if excess else None),
        support_weight_histogram=tuple(histogram), nonzero_supports_by_weight=tuple(nonzero_supports),
        raw_code_mass=data['mass'], exact_rational=True, directed_log_claim=False,
        includes_diagonal_identity=True, regional_provider_runs=0, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-start)
