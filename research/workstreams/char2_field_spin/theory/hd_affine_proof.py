"""Exact source-law/rank regressions for affine placement, not a distance certificate.

The numerical examples use GF4 and GRS[4,2]. The source-law and block-rank
identities apply to any MDS[b,k] over GF(q). No native code or timing work.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
from math import comb


def prescribed_mds_count(b, k, weight, q):
    """Codewords with one prescribed exact support, including the zero word."""
    if not 0 <= weight <= b or not 0 < k <= b:
        raise ValueError("Invalid MDS parameters")
    if weight == 0:
        return 1
    parity = b-k
    if weight <= parity:
        return 0
    return (q-1)*sum((-1)**j * comb(weight-1, j) * q**(weight-parity-1-j)
                     for j in range(weight-parity))


def row_probability(b, k, weight, q):
    return F(prescribed_mds_count(b, k, weight, q), (q-1)**weight)


def source_probability(differences, b=4, k=2, q=4):
    if len(differences) % b:
        raise ValueError("Incomplete source row")
    answer = F(1)
    for start in range(0, len(differences), b):
        weight = sum(value != 0 for value in differences[start:start+b])
        answer *= row_probability(b, k, weight, q)
    return answer


def mul(a, b):
    """GF4 with polynomial X^2+X+1."""
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return value


def inverse(a):
    if not a:
        raise ZeroDivisionError
    return mul(a, a)


def rank(matrix, width=None):
    if not matrix:
        return 0
    work = [list(row) for row in matrix]
    width = len(work[0]) if width is None else width
    pivot = 0
    for column in range(width):
        next_row = next((i for i in range(pivot, len(work)) if work[i][column]), None)
        if next_row is None:
            continue
        work[pivot], work[next_row] = work[next_row], work[pivot]
        scale = inverse(work[pivot][column])
        work[pivot] = [mul(scale, value) for value in work[pivot]]
        for i in range(len(work)):
            if i != pivot and work[i][column]:
                scale = work[i][column]
                work[i] = [left ^ mul(scale, right)
                           for left, right in zip(work[i], work[pivot])]
        pivot += 1
        if pivot == len(work):
            break
    return pivot


def matvec(matrix, vector):
    result = []
    for row in matrix:
        value = 0
        for a, b in zip(row, vector):
            value ^= mul(a, b)
        result.append(value)
    return result


def inverse_matrix(matrix):
    n = len(matrix)
    work = [list(row)+[int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = next(i for i in range(column, n) if work[i][column])
        work[column], work[pivot] = work[pivot], work[column]
        scale = inverse(work[column][column])
        work[column] = [mul(scale, value) for value in work[column]]
        for i in range(n):
            if i != column and work[i][column]:
                scale = work[i][column]
                work[i] = [left ^ mul(scale, right)
                           for left, right in zip(work[i], work[column])]
    return [row[n:] for row in work]


def difference_rows(middle, rows, strides, offsets=None):
    """Global prefix inverse difference, then the inverse fixed placement."""
    offsets = (0,)*len(strides) if offsets is None else offsets
    result = [0]*len(middle)
    previous = 0
    for c, (stride, offset) in enumerate(zip(strides, offsets)):
        if len({(stride*j+offset) % rows for j in range(rows)}) != rows:
            raise ValueError("Placement must be bijective")
        for j in range(rows):
            value = middle[c*rows+j]
            r = (stride*j+offset) % rows
            result[r*len(strides)+c] = value ^ previous
            previous = value
    return result


def admissible_tuples():
    return tuple((c, d1, d2, f) for c, d1, d2, f in product(range(1, 4), repeat=4)
                 if c ^ d1 ^ mul(f, d2))


def inverse_post(output, coefficients):
    s1 = s2 = 0
    result = []
    for y, (c, d1, d2, f) in zip(output, coefficients):
        x = mul(inverse(c), y ^ s1 ^ mul(f, s2))
        result.append(x)
        s1 ^= mul(d1, x)
        s2 ^= mul(d2, x)
    return result


def tiny_output_matrix(support=(0, 3, 6), rows=2):
    """A_T with global post/pre carry and a fixed regional permutation."""
    b, n = 4, 4*rows
    tuples = admissible_tuples()
    coefficients = [tuples[(13*i+7) % len(tuples)] for i in range(n)]
    strides = (1, 1, 1, 1) if rows == 2 else (1, 3, 1, 3)
    route = tuple(tuple((j+c) % rows for j in range(rows)) for c in range(b))
    columns = []
    for position in support:
        y = [0]*n
        y[position] = 1
        post_input = inverse_post(y, coefficients)
        middle = [post_input[c*rows+route[c][j]] for c in range(b) for j in range(rows)]
        columns.append(difference_rows(middle, rows, strides))
    return [list(column) for column in zip(*columns)]


def projective_vectors(dimension, q=4):
    """First nonzero coordinate is one; one representative of every line."""
    for first in range(dimension):
        for tail in product(range(q), repeat=dimension-first-1):
            yield (0,)*first+(1,)+tail


def exact_source_first_moment(matrix, b=4, k=2):
    return sum((source_probability(matvec(matrix, word), b, k, 4)
                for word in projective_vectors(len(matrix[0]))), F(0))


def block_rank_moment(matrix, b=4, k=2, q=4):
    """Upper source first moment using only ranks of whole inactive rows."""
    rows, dimension = len(matrix)//b, len(matrix[0])
    z = F(1, (q-1)**(b-k))
    answer = F(0)
    for mask in range(1 << rows):
        selected = [matrix[r*b+c] for r in range(rows) if mask >> r & 1 for c in range(b)]
        nullity = dimension-rank(selected)
        size = mask.bit_count()
        answer += z**(rows-size)*(1-z)**size*(q**nullity-1)
    return answer/(q-1)


def block_rank_weight_moment(matrix, output_mark, b=4, k=2, q=4):
    """Positive two-sided rank contraction, preserving exact output weights.

    Equals sum over projective y of output_mark^wt(y) z^active_rows(Ay).
    The two auxiliary subsets are independent Bernoulli sets, not a route.
    """
    rows, dimension = len(matrix)//b, len(matrix[0])
    mark, z = F(output_mark), F(1, (q-1)**(b-k))
    if not 0 < mark <= 1:
        raise ValueError("Use an output mark in (0,1]")
    answer = F(0)
    for row_mask in range(1 << rows):
        row_size = row_mask.bit_count()
        row_weight = z**(rows-row_size)*(1-z)**row_size
        for column_mask in range(1 << dimension):
            columns = [i for i in range(dimension) if column_mask >> i & 1]
            selected = [[matrix[r*b+c][i] for i in columns]
                        for r in range(rows) if row_mask >> r & 1 for c in range(b)]
            nullity = len(columns)-rank(selected)
            weight = mark**len(columns)*(1-mark)**(dimension-len(columns))
            answer += row_weight*weight*(q**nullity-1)
    return answer/(q-1)


def generic_rank(matrix, b=4, k=2):
    """Matroid-intersection rank; the MDS column matroid is block-uniform."""
    n = len(matrix)
    best = len(matrix[0])
    for mask in range(1 << n):
        left = sum(min(b-k, ((mask >> start) & ((1 << b)-1)).bit_count())
                   for start in range(0, n, b))
        right = rank([row for i, row in enumerate(matrix) if not (mask >> i & 1)])
        best = min(best, left+right)
    return best


def exact_preweight_ranks(matrix):
    """All GF4* diagonals; GRS[4,2] parity rows are (1,1,1,1),(0,1,2,3)."""
    n, dimension = len(matrix), len(matrix[0])
    histogram = Counter()
    for diagonal in product(range(1, 4), repeat=n):
        checks = []
        for start in range(0, n, 4):
            for parity in ((1, 1, 1, 1), (0, 1, 2, 3)):
                row = [0]*dimension
                for c in range(4):
                    scale = mul(parity[c], diagonal[start+c])
                    for j in range(dimension):
                        row[j] ^= mul(scale, matrix[start+c][j])
                checks.append(row)
        histogram[rank(checks)] += 1
    total = 3**n
    failure = F(sum(count for r, count in histogram.items() if r < dimension), total)
    first = sum((F(count*(4**(dimension-r)-1), 3*total)
                 for r, count in histogram.items()), F(0))
    return dict(histogram=dict(sorted(histogram.items())), failure=failure,
                projective_first_moment=first, generic_rank=generic_rank(matrix))


def regional_source_average(rows, strides, multiset=(0, 0, 1, 1)):
    """Exact regional route average for four equal prescribed value histograms.

    Small finite regression only. Distinct arrangements, not distinct values,
    are uniform. Prefix carry is global. No post-ensemble averaging is done.
    """
    if len(multiset) != rows or len(strides) != 4:
        raise ValueError("Four regions of the prescribed size are required")
    arrangements = tuple(sorted(set(permutations(multiset))))
    # Cached local edges retain the incoming last value from the prior region.
    local = {}
    for c, stride in enumerate(strides):
        for incoming in set(multiset) | {0}:
            entries = []
            for values in arrangements:
                previous, mask = incoming, 0
                for j, value in enumerate(values):
                    if value != previous:
                        mask |= 1 << ((stride*j) % rows)
                    previous = value
                entries.append((mask, previous))
            local[c, incoming] = entries
    histogram = Counter()
    for a, carry_a in local[0, 0]:
        for b, carry_b in local[1, carry_a]:
            for c, carry_c in local[2, carry_b]:
                for d, _ in local[3, carry_c]:
                    at_least_three = (a & b & c) | (a & b & d) | (a & c & d) | (b & c & d)
                    if a | b | c | d != at_least_three:
                        continue
                    four = a & b & c & d
                    histogram[((at_least_three ^ four).bit_count(), four.bit_count())] += 1
    count = len(arrangements)**4
    probability = sum((F(number, 3**(2*three+3*four))
                       for (three, four), number in histogram.items()), F(0))/count
    return dict(arrangements=count, feasible_arrangements=sum(histogram.values()),
                row_weight_histogram=dict(sorted(histogram.items())), probability=probability)


def dense_profile():
    """Directed-Arb calculation of a sufficient, but impossible uniform premise."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = 128
    try:
        rows, b, parity, q = 32768, 64, 32, 65536
        n, cutoff = rows*b, (2*rows*b)//5
        logq, logtwo = arb(q-1).log(), arb(2).log()
        logbinomial = arb(n+1).lgamma()-arb(cutoff+1).lgamma()-arb(n-cutoff+1).lgamma()
        numerator = arb(cutoff).log()+logbinomial+(cutoff-1)*logq+40*logtwo
        threshold = numerator/(parity*logq)
        sufficient = int(float(threshold))+1
        def exponent(active):
            return (arb(cutoff).log()+logbinomial+(cutoff-1-parity*active)*logq)/logtwo
        if not exponent(sufficient) < -40 or not exponent(sufficient-1) > -40:
            raise ArithmeticError("Adjacent integer thresholds not certified")
        # At this many active source rows, dimension alone guarantees a
        # nonzero output in EVERY prescribed set of cutoff coordinates.
        forced_sparse_rows = (n-cutoff)//b+1
        return dict(cutoff=cutoff, uniform_active_rows_needed=sufficient,
                    required_fraction=sufficient/rows,
                    upper_log2_at_threshold=str(exponent(sufficient)),
                    impossible_uniform_premise_by_rows=forced_sparse_rows,
                    probability_certificate=False)
    finally:
        ctx.prec = previous


def main():
    matrix = tiny_output_matrix()
    print("tiny full-map ranks:", exact_preweight_ranks(matrix))
    print("exact source first moment:", exact_source_first_moment(matrix))
    print("block-rank upper:", block_rank_moment(matrix))
    print("weight-marked block-rank upper:", block_rank_weight_moment(matrix, F(1, 2)))
    for rows, multiset in ((4, (0, 0, 1, 1)), (4, (0, 0, 1, 2)),
                           (8, (0, 0, 0, 0, 0, 0, 1, 1))):
        for strides in ((1, 1, 1, 1), (1, 3, 5, 7)):
            print("regional source average:", rows, multiset, strides,
                  regional_source_average(rows, strides, multiset))
    print("K20 dense-profile limitation:", dense_profile())


if __name__ == "__main__":
    main()
