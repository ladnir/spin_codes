"""Exact-arithmetic bound on joint supports from an ordinary weight enumerator.

This is a proof-component check, not a SPIN certificate or benchmark helper.
Run `python joint_support.py` for exhaustive small-code checks.
"""

from fractions import Fraction
from itertools import product
from math import factorial, prod


def gf2_rank(words):
    basis = {}
    for value in words:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def griesmer(d, h):
    return sum((d + (1 << i) - 1) >> i for i in range(h))


def cumulative_basis_counts(spectrum, h):
    """Nonzero h-tuples of codewords, allowing dependence, by total weight."""
    nonzero = [(w, count) for w, count in enumerate(spectrum) if w and count]
    coefficients = [1]
    for _ in range(h):
        next_coefficients = [0] * (len(coefficients) + len(spectrum) - 1)
        for weight, count in nonzero:
            for old_weight, old_count in enumerate(coefficients):
                next_coefficients[old_weight + weight] += old_count * count
        coefficients = next_coefficients
    cumulative, count = [], 0
    for coefficient in coefficients:
        count += coefficient
        cumulative.append(count)
    return cumulative


def joint_bound(spectrum, g, h, u, *, cumulative=None):
    """Upper-bound ordered g-tuples of rank h with union support <= u.

    spectrum must be the exact enumerator or coefficientwise upper bounds;
    its first occupied nonzero shell must be the true minimum distance or a
    valid lower bound. Inputs and outputs use exact integer/rational arithmetic.
    """
    if not 1 <= h <= g or not 0 <= u < len(spectrum):
        raise ValueError("invalid rank/group/support")
    d = next(w for w, count in enumerate(spectrum) if w and count)
    if h == 1:
        return ((1 << g) - 1) * sum(spectrum[1:u + 1]), u
    if u < griesmer(d, h):
        return 0, None
    if cumulative is None:
        cumulative = cumulative_basis_counts(spectrum, h)
    mean = Fraction(h * (1 << (h - 1)) * u, (1 << h) - 1)
    minimum = h * d
    if mean < minimum:
        return 0, None
    bases = prod((1 << h) - (1 << j) for j in range(h))
    spanning = prod((1 << g) - (1 << j) for j in range(h))
    best, best_threshold = None, None
    # At these thresholds the guaranteed fraction of inexpensive bases is >0.
    first = max(minimum, mean.numerator // mean.denominator)
    for threshold in range(first, len(cumulative)):
        # No admissible lower count of good bases exceeds all bases. Once
        # even that denominator cannot improve the result, later thresholds
        # cannot improve it either (cumulative is nondecreasing).
        if best is not None and spanning * cumulative[threshold] // bases >= best:
            break
        numerator = bases * ((threshold + 1) * mean.denominator - mean.numerator)
        denominator = (threshold + 1 - minimum) * mean.denominator
        good_bases = max(factorial(h), -(-numerator // denominator))
        integer_upper = spanning * cumulative[threshold] // good_bases
        if best is None or integer_upper < best:
            best, best_threshold = integer_upper, threshold
    if best is None:
        raise AssertionError("no admissible basis threshold")
    return best, best_threshold


def span(rows):
    if gf2_rank(rows) != len(rows):
        raise ValueError("dependent test generator")
    words = [0]
    for row in rows:
        words += [word ^ row for word in words]
    return words


def check_code(rows, n, g):
    words = span(rows)
    spectrum = [0] * (n + 1)
    for word in words:
        spectrum[word.bit_count()] += 1
    actual = [[0] * (n + 1) for _ in range(g + 1)]
    for word_tuple in product(words, repeat=g):
        union = 0
        for word in word_tuple:
            union |= word
        actual[gf2_rank(word_tuple)][union.bit_count()] += 1
    tested = 0
    for h in range(1, min(g, len(rows)) + 1):
        cumulative = cumulative_basis_counts(spectrum, h)
        count = 0
        for u in range(n + 1):
            count += actual[h][u]
            upper, _ = joint_bound(spectrum, g, h, u, cumulative=cumulative)
            if upper < count:
                raise AssertionError((rows, g, h, u, count, upper))
            # At full length, all tuples are counted and dependence is still
            # conservatively included in the basis polynomial.
            tested += 1
    return tested, len(words) ** g


def self_test():
    cases = [
        ([1, 2, 4], 3, 4),
        ([0b1111000, 0b1100110, 0b1010101], 7, 4),
        ([0b1111, 0b0011], 4, 4),
        ([0b10010111, 0b01001011, 0b00101101, 0b00011110], 8, 3),
        ([0b111000000, 0b000111000, 0b000000111], 9, 4),
    ]
    inequalities = tuples = 0
    for rows, n, g in cases:
        tested, enumerated = check_code(rows, n, g)
        inequalities += tested
        tuples += enumerated
    # An upper enumerator remains safe under coefficientwise inflation.
    exact = [1, 0, 0, 0, 7, 0, 0, 0]
    inflated = [1, 0, 0, 0, 14, 0, 0, 0]
    for u in range(8):
        for h in range(1, 4):
            assert joint_bound(inflated, 4, h, u)[0] >= joint_bound(exact, 4, h, u)[0]
    for u in range(8):
        assert joint_bound(exact, 4, 1, u)[0] == 15 * sum(exact[1:u + 1])
    print(f"passed {inequalities} inequalities against {tuples} exhaustively enumerated tuples")
    print("passed coefficientwise inflation and exact rank-one checks")


if __name__ == "__main__":
    self_test()
