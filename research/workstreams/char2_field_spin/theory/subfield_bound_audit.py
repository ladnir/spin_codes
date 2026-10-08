"""Independent algebra for regional GF(q) whole-code bounds.

The functions accepting rational matrices are exact. ``bulk_screen`` uses
floating-point logarithms and is only a search aid, not a certificate.
"""

import argparse
from fractions import Fraction as F
from itertools import product
import json
import math


IDENTITY = ((F(1), F(0)), (F(0), F(1)))
ZERO = ((F(0), F(0)), (F(0), F(0)))


def matadd(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(2)) for i in range(2))


def matscale(a, scale):
    return tuple(tuple(scale * value for value in row) for row in a)


def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def matpower(a, exponent):
    result = IDENTITY
    while exponent:
        if exponent & 1:
            result = matmul(result, a)
        a = matmul(a, a)
        exponent >>= 1
    return result


def terminal_from_zero(a):
    return a[0][0] + a[0][1]


def fixed_support_multiplicity(length, dimension, weight, q):
    """Nonzero MDS words with precisely one prescribed nonempty support."""
    if weight == 0:
        return 0
    return sum((-1) ** j * math.comb(weight, j)
               * q ** max(0, weight - j - (length - dimension))
               for j in range(weight + 1))


def nonzero_word_enumerator(length, dimension, q, x):
    return sum(math.comb(length, weight)
               * fixed_support_multiplicity(length, dimension, weight, q)
               * x ** weight
               for weight in range(length - dimension + 1, length + 1))


def coupled_fixed_rows_bound(regional, active_rows, length, dimension, q):
    """Moment bound for specified active rows, before projective division.

regional[h] already includes its binom(L,h) normalization. The relaxed
row-support count is q**(weight-(length-dimension)), even for impossible
small supports. No extra nonzero-value factor belongs in this function.
"""
    combined = ZERO
    for h in range(active_rows + 1):
        combined = matadd(combined, matscale(
            regional[h], math.comb(active_rows, h) * q ** h))
    return (F(1, q ** ((length - dimension) * active_rows))
            * terminal_from_zero(matpower(combined, length)))


def literal_fixed_rows_moment(regional, active_rows, length, dimension, q):
    """Small exact test oracle; enumerate outer support tuples, not values."""
    supports = []
    for mask in range(1, 1 << length):
        multiplicity = fixed_support_multiplicity(
            length, dimension, mask.bit_count(), q)
        if multiplicity:
            supports.append((mask, multiplicity))
    total = F(0)
    for rows in product(supports, repeat=active_rows):
        matrix = IDENTITY
        multiplicity = math.prod(count for _, count in rows)
        for coordinate in range(length):
            h = sum((mask >> coordinate) & 1 for mask, _ in rows)
            matrix = matmul(matrix, regional[h])
        total += multiplicity * terminal_from_zero(matrix)
    return total


def potential_norm(matrix, eta):
    assert eta > 0
    return max(matrix[0][0] + eta * matrix[0][1],
               matrix[1][0] / eta + matrix[1][1])


def scalar_enumerator_bound(regional, active_rows, length, dimension, q,
                            eta, rho):
    """Check a proposed regional envelope and apply its outer enumerator."""
    assert potential_norm(regional[0], eta) <= 1
    assert all(potential_norm(regional[h], eta) <= rho ** h
               for h in range(1, active_rows + 1))
    return (max(F(1), 1 / eta)
            * nonzero_word_enumerator(length, dimension, q, rho) ** active_rows)


def _log2_add(values):
    peak = max(values)
    return peak + math.log2(sum(2.0 ** (value - peak) for value in values))


def _log2_outer(length, dimension, q, x):
    return _log2_add([
        math.log2(math.comb(length, weight)
                   * fixed_support_multiplicity(length, dimension, weight, q))
        + weight * math.log2(x)
        for weight in range(length - dimension + 1, length + 1)])


def bulk_screen(k=2 ** 16, length=16, dimension=8, q=256, distance=0.1):
    """Search one tilt bounding the entire active-row tail by row sums.

For a fixed tilt, successive binomial terms have decreasing ratio.
When r=(L-a)/(a+1)*P(f)<1, their tail is at most its first term/(1-r).
No claim of outward rounding is made by this numerical search.
"""
    rows = k // dimension
    threshold = math.ceil(distance * rows * length) - 1
    best = None
    for index in range(1001):
        z = math.exp(-12 + 11.999 * index / 1000)
        f = z + (1 - z) / (q - 1)
        lp = _log2_outer(length, dimension, q, f)
        if lp >= 0:
            continue
        p = 2.0 ** lp

        def margin(first):
            ratio = (rows - first) / (first + 1) * p
            if ratio >= 1:
                return -math.inf
            lc = (math.lgamma(rows + 1) - math.lgamma(first + 1)
                  - math.lgamma(rows - first + 1)) / math.log(2)
            return -(lc + first * lp - math.log2(q - 1)
                     - threshold * math.log2(z) - math.log2(1 - ratio))

        if margin(rows) < 40:
            continue
        low, high = 1, rows
        while low < high:
            middle = (low + high) // 2
            if margin(middle) >= 40:
                high = middle
            else:
                low = middle + 1
        record = dict(first_active_rows=low, last_active_rows=rows, z=z,
                      margin_bits=margin(low), log2_outer_moment=lp,
                      rigorous=False)
        if best is None or low < best['first_active_rows']:
            best = record
    return best


def exact_bulk_tail_check(first=2879, z=F(7, 2500), margin_bits=64,
                          k=2 ** 16, length=16, dimension=8, q=256,
                          distance=F(1, 10)):
    """Exact integer proof of the whole geometric active-row tail.

    Comparison is strict, including one global projective division and
    the discarded-final-state terminal sum. No floating-point arithmetic
    or numerical logarithm participates in the claimed result.
    """
    assert k % dimension == 0 and 0 < z < 1 and margin_bits >= 0
    rows = k // dimension
    assert 1 <= first <= rows
    desired = distance * rows * length
    threshold = (desired.numerator + desired.denominator - 1) // desired.denominator - 1
    f = z + (1 - z) / (q - 1)
    p = nonzero_word_enumerator(length, dimension, q, f)
    pn, pd = p.numerator, p.denominator
    tail_denominator = (first + 1) * pd - (rows - first) * pn
    assert tail_denominator > 0
    # Multiply through every positive denominator in
    # binom(L,a) P^a z^-e / ((q-1)(1-r)) < 2^-margin_bits.
    left = (math.comb(rows, first) * pn ** first * z.denominator ** threshold
            * (first + 1) * pd * (1 << margin_bits))
    right = ((q - 1) * pd ** first * z.numerator ** threshold * tail_denominator)
    return dict(first_active_rows=first, last_active_rows=rows,
                threshold=threshold, z=[z.numerator, z.denominator],
                margin_bits=margin_bits, rigorous=left < right,
                comparison_bit_lengths=[left.bit_length(), right.bit_length()])


def exact_window_lower_check(active_rows=32, window=816, state=2,
                             lower_bits=80, k=2 ** 16, length=16,
                             dimension=8, q=256, width=16,
                             distance=F(1, 10)):
    """Exact lower bound on expected bad projective source lines.

    Fix one active-row set and a suffix window in every region. The
    event concerns all coordinates of those rows, regardless of their
    eventual message. Requiring zero state at the first length-1 region
    ends imposes at most (length-1)*state homogeneous equations. This
    lower bound is NOT the same-size lower bound on code-failure probability.
    """
    rows = k // dimension
    assert rows % width == window % width == 0
    assert 1 <= active_rows <= window <= rows
    desired = distance * rows * length
    threshold = (desired.numerator + desired.denominator - 1) // desired.denominator - 1
    assert length * window <= threshold
    surviving_dimension = dimension * active_rows - (length - 1) * state
    assert surviving_dimension > 0
    numerator = (math.comb(window, active_rows) ** length
                 * (q ** surviving_dimension - 1))
    denominator = (math.comb(rows, active_rows) ** length * (q - 1))
    return dict(active_rows=active_rows, window_symbols=window,
                maximum_output_weight=length * window, threshold=threshold,
                surviving_dimension=surviving_dimension,
                expected_bad_lines_greater_than_power_of_two=lower_bits,
                rigorous=numerator > denominator * (1 << lower_bits),
                expectation_not_failure_probability=True,
                comparison_bit_lengths=[numerator.bit_length(), denominator.bit_length()])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bulk-screen', action='store_true')
    parser.add_argument('--verify-bulk', action='store_true')
    parser.add_argument('--verify-bulk-three-percent', action='store_true')
    parser.add_argument('--verify-bulk-gf16', action='store_true')
    parser.add_argument('--verify-window-obstruction', action='store_true')
    args = parser.parse_args()
    if args.bulk_screen:
        print(json.dumps(bulk_screen(), indent=2))
    if args.verify_bulk:
        print(json.dumps(exact_bulk_tail_check(), indent=2))
    if args.verify_bulk_three_percent:
        print(json.dumps(exact_bulk_tail_check(
            first=909, z=F(1, 400), distance=F(3, 100), margin_bits=64), indent=2))
    if args.verify_bulk_gf16:
        print(json.dumps(exact_bulk_tail_check(
            first=5998, z=F(1, 52), q=16, distance=F(1, 10), margin_bits=64), indent=2))
    if args.verify_window_obstruction:
        print(json.dumps(exact_window_lower_check(), indent=2))
