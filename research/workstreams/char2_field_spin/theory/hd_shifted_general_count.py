"""Exact support counts for two matchings of an indexed binary cycle.

Counts and coefficient bounds use integer or rational arithmetic. The CLI
prints floating logarithms for inspection; it is not a distance certificate.
"""

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import comb, log2

from hd_shifted_pairs_count import circular_run_count


def _validate_pairs(pairs):
    if not isinstance(pairs, int) or pairs < 1:
        raise ValueError("pairs must be a positive integer")


def matching_histograms(pairs, even, odd, even_full, odd_full):
    """Return ordered (00,10,01,11) histograms, even matching first."""
    weight = even + odd
    return ((pairs - weight + even_full, even - even_full,
             odd - even_full, even_full),
            (pairs - weight + odd_full, odd - odd_full,
             even - odd_full, odd_full))


def feasible_witness(pairs, even, odd, even_full, odd_full):
    """Return (even-cycle runs, selected odd bits on 11 edges), or None.

    Constant even cycles have zero runs. Every returned witness constructs
    at least one support in the class; this is an exact feasibility test.
    """
    _validate_pairs(pairs)
    values = (even, odd, even_full, odd_full)
    if any(not isinstance(value, int) or not 0 <= value <= pairs
           for value in values):
        return None
    if any(min(histogram) < 0 for histogram in
           matching_histograms(pairs, *values)):
        return None
    runs_lo = int(0 < even < pairs)
    runs_hi = min(even, pairs - even)
    selected_lo = max(0, even_full + odd_full - odd)
    selected_hi = min(even_full, odd_full)
    total_lo = max(even_full, odd_full, runs_lo + selected_lo)
    total_hi = min(even, pairs - even - odd + even_full + odd_full,
                   runs_hi + selected_hi)
    if total_lo > total_hi:
        return None
    runs = max(runs_lo, total_lo - selected_hi)
    selected = total_lo - runs
    return runs, selected


def joint_is_feasible(pairs, even, odd, even_full, odd_full):
    return feasible_witness(pairs, even, odd, even_full, odd_full) is not None


def oriented_pairing_count(pairs, even, odd, full):
    """Exact marginal count for a matching with both lane counts fixed."""
    _validate_pairs(pairs)
    empty, first, second, both = (pairs - even - odd + full,
                                 even - full, odd - full, full)
    if min(empty, first, second, both) < 0:
        return 0
    return (comb(pairs, both) * comb(pairs - both, first)
            * comb(pairs - both - first, second))


def joint_count_upper(pairs, even, odd, even_full, odd_full):
    """Sound integer upper bound for the specified four-parameter class.

    This is inexpensive for large cycles. Use exact_joint_count for the
    stronger O(pairs**2) reference. All inputs count indexed positions.
    """
    if not joint_is_feasible(pairs, even, odd, even_full, odd_full):
        return 0
    weight = even + odd
    runs = weight - even_full - odd_full
    if weight in (0, 2 * pairs):
        return 1
    if even == odd and even in (even_full, odd_full):
        return circular_run_count(pairs, even, runs)
    return min(comb(pairs, even) * comb(pairs, odd),
               oriented_pairing_count(pairs, even, odd, even_full),
               oriented_pairing_count(pairs, even, odd, odd_full),
               circular_run_count(2 * pairs, weight, runs))


def _cycle_run_range(pairs, even, odd, even_full, odd_full):
    """Feasible even-cycle run interval, assuming a feasible full profile."""
    runs = even + odd - even_full - odd_full
    selected_lo = max(0, even_full + odd_full - odd)
    selected_hi = min(even_full, odd_full)
    return (max(int(0 < even < pairs),
                max(even_full, odd_full) - selected_hi),
            min(even, pairs - even,
                min(even, pairs - runs) - selected_lo))


def _run_count_interval(length, weight, low, high):
    if low > high:
        return 0
    count = circular_run_count(length, weight, low)
    result = count
    for runs in range(low, high):
        numerator = count * (weight - runs) * (length - weight - runs)
        denominator = runs * (runs + 1)
        assert numerator % denominator == 0
        count = numerator // denominator
        result += count
    return result


def joint_run_count_upper(pairs, even, odd, even_full, odd_full):
    """Refine joint_count_upper using runs of each parity's own cycle.

    For fixed even bits, either matching admits an exact conditional
    binomial count for the odd bits. Restricting the even cycles to their
    feasible run interval gives a stronger upper bound in O(pairs) exact
    integer operations. Exchange even and odd for a second bound.
    """
    upper = joint_count_upper(pairs, even, odd, even_full, odd_full)
    if upper <= 1:
        return upper
    for first, second in ((even, odd), (odd, even)):
        low, high = _cycle_run_range(pairs, first, second, even_full, odd_full)
        cycles = _run_count_interval(pairs, first, low, high)
        conditional = min(comb(first, full) * comb(pairs - first, second - full)
                          for full in (even_full, odd_full))
        upper = min(upper, cycles * conditional)
    return upper


def exact_joint_count(pairs, even, odd, even_full, odd_full):
    """Exact class count from even-cycle runs and odd-bit placements.

    At most O(pairs**2) binomial products. Intended for reference classes,
    not an enumeration of every large-cycle class.
    """
    if not joint_is_feasible(pairs, even, odd, even_full, odd_full):
        return 0
    weight = even + odd
    runs = weight - even_full - odd_full
    if weight in (0, 2 * pairs):
        return 1
    if even == odd and even in (even_full, odd_full):
        return circular_run_count(pairs, even, runs)
    # Complement and lane exchange preserve support multiplicities. Choose
    # a parity with the smallest nonconstant cycle to shorten the run sum.
    variants = ((even, odd, even_full, odd_full),
                (odd, even, even_full, odd_full),
                (pairs - even, pairs - odd, pairs - weight + even_full,
                 pairs - weight + odd_full),
                (pairs - odd, pairs - even, pairs - weight + even_full,
                 pairs - weight + odd_full))
    even, odd, even_full, odd_full = min(variants, key=lambda item: item[0])
    weight = even + odd
    runs = weight - even_full - odd_full
    selected_lo = max(0, even_full + odd_full - odd)
    selected_hi = min(even_full, odd_full)
    cycle_runs_lo, cycle_runs_hi = _cycle_run_range(pairs, even, odd, even_full, odd_full)
    result = 0
    for cycle_runs in range(cycle_runs_lo, cycle_runs_hi + 1):
        n00 = pairs - even - cycle_runs
        n11 = even - cycle_runs
        lo = max(selected_lo, even_full - cycle_runs,
                 odd_full - cycle_runs)
        hi = min(selected_hi, n11, pairs - runs - cycle_runs)
        if lo > hi:
            continue
        placements = 0
        for selected in range(lo, hi + 1):
            placements += (comb(n00, odd - even_full - odd_full + selected)
                           * comb(cycle_runs, even_full - selected)
                           * comb(cycle_runs, odd_full - selected)
                           * comb(n11, selected))
        result += circular_run_count(pairs, even, cycle_runs) * placements
    return result


def exact_profiles(pairs):
    """Small exact DP mapping (even,odd,even_full,odd_full) to counts."""
    _validate_pairs(pairs)
    result = defaultdict(int)
    for first in (0, 1):
        states = {(first, 0, 0, 0, first): 1}
        for position in range(1, 2 * pairs):
            following = defaultdict(int)
            for (even, odd, a, b, last), multiplicity in states.items():
                for bit in (0, 1):
                    product = last * bit
                    key = (even + (bit if position % 2 == 0 else 0),
                           odd + (bit if position % 2 else 0),
                           a + (product if position % 2 else 0),
                           b + (product if position % 2 == 0 else 0), bit)
                    following[key] += multiplicity
            states = following
        for (even, odd, a, b, last), multiplicity in states.items():
            result[even, odd, a, b + last * first] += multiplicity
    return dict(result)


def _matrix_multiply(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return (a * e + b * g, a * f + b * h,
            c * e + d * g, c * f + d * h)


def transfer_value(pairs, markers=(1, 1, 1, 1)):
    """Exact positive generating function evaluated at rational markers.

    Markers are ordered (even, odd, even_full, odd_full). The result is
    trace(K**pairs), K[u,v] = x**v * (1 + y * p**u * q**v).
    """
    _validate_pairs(pairs)
    if len(markers) != 4:
        raise ValueError("four positive rational markers required")
    if any(isinstance(value, float) for value in markers):
        raise TypeError("use integers or Fraction markers, not floats")
    x, y, p, q = tuple(map(Fraction, markers))
    if min(x, y, p, q) <= 0:
        raise ValueError("markers must be positive")
    base = (1 + y, x * (1 + y * q), 1 + y * p, x * (1 + y * p * q))
    result = (Fraction(1), Fraction(0), Fraction(0), Fraction(1))
    exponent = pairs
    while exponent:
        if exponent & 1:
            result = _matrix_multiply(result, base)
        exponent >>= 1
        if exponent:
            base = _matrix_multiply(base, base)
    return result[0] + result[3]


@dataclass(frozen=True)
class ProfileBox:
    """Inclusive ranges of (even, odd, even_full, odd_full)."""

    even: tuple
    odd: tuple
    even_full: tuple
    odd_full: tuple

    @property
    def ranges(self):
        return self.even, self.odd, self.even_full, self.odd_full

    def split(self):
        """Bisect a longest coordinate; children partition the parent."""
        widths = tuple(hi - lo for lo, hi in self.ranges)
        coordinate = max(range(4), key=widths.__getitem__)
        if widths[coordinate] == 0:
            raise ValueError("cannot split a singleton box")
        lo, hi = self.ranges[coordinate]
        middle = (lo + hi) // 2
        left, right = list(self.ranges), list(self.ranges)
        left[coordinate] = lo, middle
        right[coordinate] = middle + 1, hi
        return ProfileBox(*left), ProfileBox(*right)


def box_count_upper(pairs, box, markers=(1, 1, 1, 1)):
    """Bound the total support multiplicity in a box without listing it.

    Positive rational markers give an exact coefficient bound. The min of
    this bound and the parity-binomial marginal is returned. A singleton
    additionally uses joint_count_upper, including its feasibility check.
    """
    _validate_pairs(pairs)
    if any(not isinstance(value, int) or not 0 <= value <= pairs
           for interval in box.ranges for value in interval):
        raise ValueError("box endpoints must be integers in [0,pairs]")
    if any(lo > hi for lo, hi in box.ranges):
        return 0
    if len(markers) != 4:
        raise ValueError("four positive rational markers required")
    generating_value = transfer_value(pairs, markers)
    minimum_monomial = Fraction(1)
    for marker, (lo, hi) in zip(map(Fraction, markers), box.ranges):
        minimum_monomial *= marker ** (lo if marker >= 1 else hi)
    tilted = generating_value / minimum_monomial
    upper = tilted.numerator // tilted.denominator
    parity_count = 1
    for lo, hi in (box.even, box.odd):
        parity_count *= sum(comb(pairs, count) for count in range(lo, hi + 1))
    upper = min(upper, parity_count)
    if all(lo == hi for lo, hi in box.ranges):
        upper = min(upper, joint_count_upper(pairs, *(lo for lo, _ in box.ranges)))
    return upper


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=4096)
    parser.add_argument("--even", type=int, default=512)
    parser.add_argument("--odd", type=int, default=512)
    parser.add_argument("--even-full", type=int, default=512)
    parser.add_argument("--odd-full", type=int, default=384)
    parser.add_argument("--exact", action="store_true")
    args = parser.parse_args()
    profile = (args.pairs, args.even, args.odd, args.even_full, args.odd_full)
    count = exact_joint_count(*profile) if args.exact else joint_count_upper(*profile)
    print({"count": count, "log2_count": log2(count) if count else float("-inf"),
           "exact_count": args.exact, "distance_certificate": False})
