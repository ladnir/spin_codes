"""Exact local moments for the GF(q) diagonal / refreshed-state candidate.

These are local algebra checks, not a whole-code distance certificate.
The model uses A_j=(1,j), C=A^T, independent nonzero input scalars,
and a uniform nonzero refreshed state when the incoming state is nonzero.
"""

from fractions import Fraction
from itertools import combinations


def multiply(a, b, modulus):
    result = 0
    degree = modulus.bit_length() - 1
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        if a >> degree:
            a ^= modulus
        b >>= 1
    return result


def nonzero_sum_count(columns, target, q, modulus):
    """Count x_j != 0 with sum_j (1,j)*x_j = target exactly."""
    n = len(columns)
    h0, h1 = target
    if n == 0:
        return int(h0 == h1 == 0)
    sign = -1 if n & 1 else 1
    if h0 == h1 == 0:
        numerator = (q - 1)**n + sign * (q*q - 1 - n*q*(q - 1))
    else:
        on_line = int(h0 != 0 and any(
            h1 == multiply(j, h0, modulus) for j in columns))
        numerator = (q - 1)**n + sign * (-1 + n*q - on_line*q*q)
    assert numerator >= 0 and numerator % (q*q) == 0
    return numerator // (q*q)


def zero_state_probability(weight, q):
    if weight == 0:
        return Fraction(1)
    sign = -1 if weight & 1 else 1
    numerator = ((q - 1)**weight
                 + sign * (q*q - 1 - weight*q*(q - 1)))
    return Fraction(numerator, q*q * (q - 1)**weight)


def active_output_moment(weight, width, q, z):
    """E[z^wt(Y)] for a nonzero incoming state, before its refresh."""
    z = Fraction(z)
    assert 0 <= weight <= width <= q and 0 <= z <= 1
    f = Fraction(1, q - 1) + Fraction(q - 2, q - 1)*z
    total = Fraction(q + 1 - width, q + 1)*z**(width - weight)*f**weight
    if weight:
        total += Fraction(weight, q + 1)*z**(width - weight + 1)*f**(weight - 1)
    if weight < width:
        total += Fraction(width - weight, q + 1)*z**(width - weight - 1)*f**weight
    return total


def active_return_moment(support, width, q, modulus, z):
    """E[1[Q_next=0] z^wt(Y)] via projective states and subset expansion.

    Cost is O((q+1) * 2^len(support) * len(support)^2). This is a
    transparent exact oracle for short supports, not a fast certificate engine.
    """
    support = tuple(sorted(support))
    assert len(support) == len(set(support))
    assert all(0 <= j < width for j in support)
    assert width <= q
    z = Fraction(z)
    assert 0 <= z <= 1
    w = len(support)
    support_set = set(support)
    total = Fraction(0)
    # Representatives for all q+1 one-dimensional subspaces of F_q^2.
    for r0, r1 in [(1, a) for a in range(q)] + [(0, 1)]:
        values = [r0 ^ multiply(j, r1, modulus) for j in range(width)]
        forced_nonzero = sum(values[j] != 0 for j in range(width) if j not in support_set)
        eligible = [j for j in support if values[j] != 0]
        for size in range(len(eligible) + 1):
            for fixed in combinations(eligible, size):
                fixed_set = set(fixed)
                h0, h1 = r0, r1
                for j in fixed:
                    h0 ^= values[j]
                    h1 ^= multiply(j, values[j], modulus)
                remaining = tuple(j for j in support if j not in fixed_set)
                count = nonzero_sum_count(remaining, (h0, h1), q, modulus)
                total += (z**(forced_nonzero + w - size)
                          * (1-z)**size * count)
    return total / ((q + 1) * (q - 1)**w)


def transition(support, width, q, modulus, z):
    """Exact 2x2 matrix, rows/columns ordered zero-state, nonzero-state."""
    w = len(support)
    z = Fraction(z)
    p0 = zero_state_probability(w, q)
    g = active_output_moment(w, width, q, z)
    j = active_return_moment(support, width, q, modulus, z)
    assert 0 <= j <= g
    return ((z**w * p0, z**w * (1-p0)), (j, g-j))
