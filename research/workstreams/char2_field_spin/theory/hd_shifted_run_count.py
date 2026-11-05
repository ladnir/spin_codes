"""Exact whole-box support counts using weight and cyclic-run markers.

For an indexed binary cycle of length 2*M, let r be its weight and let
a,b count 11 edges in the alternating matchings. Its run statistic is
c=r-a-b, including c=0 for both constant cycles. The positive matrix

    T = [[1, x*y], [1, x]]

satisfies trace(T**(2*M)) = sum_U x**r(U) y**c(U). A 0-to-1 edge starts
a run. The trace closes the indexed cycle; it does not divide by rotations.

For an (r,a,b) box, divide this generating value by the minimum of
(x*y)**r * y**(-a) * y**(-b) on the box. Its integer floor bounds the
TOTAL number of supports in the box. No profile-count factor is needed.
"""
from fractions import Fraction as F
from functools import lru_cache

from hd_shifted_unsigned_geometry import unsigned_box_witness


def _markers(x, y):
    if isinstance(x, float) or isinstance(y, float):
        raise TypeError('Use exact rational run markers, not floats')
    x, y = F(x), F(y)
    if min(x, y) <= 0:
        raise ValueError('Positive weight and run markers are required')
    return x, y


def _product(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


@lru_cache(32)
def _value(pairs, x, y):
    base = (F(1), x*y, F(1), x)
    result = (F(1), F(0), F(0), F(1))
    exponent = 2*pairs
    while exponent:
        if exponent & 1:
            result = _product(result, base)
        exponent >>= 1
        if exponent:
            base = _product(base, base)
    return result[0]+result[3]


def transfer_value(pairs, x, y):
    """Exact rational trace(T**(2*pairs)), including both constants."""
    if not isinstance(pairs, int) or pairs < 1:
        raise ValueError('pairs must be a positive integer')
    x, y = _markers(x, y)
    return _value(pairs, x, y)


def box_count_upper(pairs, ranges, x, y):
    """Exact integer upper bound on all indexed supports in the box.

    Each range is inclusive. Empty physical boxes return zero. If the
    box contains (0,0,0), the count includes the zero support.
    """
    ranges = tuple(tuple(row) for row in ranges)
    x, y = _markers(x, y)
    if unsigned_box_witness(pairs, ranges) is None:
        return 0
    minimum = F(1)
    for marker, (lo, hi) in zip((x*y, 1/y, 1/y), ranges):
        minimum *= marker**(lo if marker >= 1 else hi)
    bound = transfer_value(pairs, x, y)/minimum
    return bound.numerator//bound.denominator


def box_log2_upper(pairs, ranges, x, y, *, precision=256):
    """Directed Arb upper endpoint of log2 of the whole-box count bound."""
    from flint import arb, ctx
    if precision < 128:
        raise ValueError('At least 128 bits of precision are required')
    count = box_count_upper(pairs, ranges, x, y)
    if count == 0:
        raise ValueError('An empty box has no finite logarithmic count')
    with ctx.workprec(precision):
        return (arb(count).log()/arb(2).log()).upper()


def stationary_markers(pairs, weight, runs):
    """Exact marker proposal matching a two-state empirical transition law.

    This chooses a witness only; every count bound is checked independently.
    Strictly interior run counts are required so both markers are positive.
    """
    if (not isinstance(pairs, int) or pairs < 1
            or any(not isinstance(n, int) for n in (weight, runs))
            or not 0 < runs < min(weight, 2*pairs-weight)):
        raise ValueError('An interior integer weight/run pair is required')
    length = 2*pairs
    x = F((weight-runs)*(length-weight), weight*(length-weight-runs))
    y = F(runs*runs, (weight-runs)*(length-weight-runs))
    return x, y
