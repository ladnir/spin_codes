"""Exact indexed cyclic support counts for the two alternating matchings.

This is a scalar combinatorial helper, not a distance certificate. A positive
last-exit-tree formula counts each fixed phase flow. Consecutive nonzero
four-vertex terms use exact integer division; constant-phase endpoints are
handled directly. Swap, complement, and reflection reduce the work. No
transfer coefficients, support sums, or frozen providers are modified.
"""

from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path


MAX_PAIRS = 4096
NOTE = Path(__file__).resolve().parent.parent/'HD_FIELD128_EXACT_PROFILE_COUNTS.md'


def _checked(pairs, weight, first, second):
    if (any(type(x) is not int for x in (pairs, weight, first, second))
            or not 1 <= pairs <= MAX_PAIRS or not 0 <= weight <= 2*pairs):
        raise ValueError('Expected integer P in 1..4096 and weight in 0..2P')
    low, high = max(0, weight-pairs), weight//2
    return low <= first <= high and low <= second <= high


def _bounds(weight, first, second):
    singles = weight-2*first
    return max(0, second-first), min(singles, weight-first-second)


def _direct(pairs, weight, first, second, u):
    """One feasible fixed-flow count, using four binomials, not factorials."""
    v = weight-2*first-u
    e, f = first+u-second, first+v-second
    e0, e1 = pairs-first-u, first+u
    o0, o1 = pairs-first-v, first+v
    # A missing phase vertex fixes every bit on that phase. The remaining
    # phase is an arbitrary indexed word of its prescribed weight.
    if not e0:
        return comb(pairs, first)
    if not e1:
        return comb(pairs, v)
    if not o0:
        return comb(pairs, first)
    if not o1:
        return comb(pairs, u)
    tau = u*(pairs-weight+second)*o1 + first*f*o0
    if not tau:
        return 0
    numerator = pairs*tau*comb(e0, v)*comb(e1, first)*comb(o0, e)*comb(o1, second)
    value, remainder = divmod(numerator, e0*e1*o0*o1)
    if remainder:
        raise ArithmeticError('Nonintegral last-exit-tree count')
    return value


def fixed_flow_count(pairs, weight, first, second, u):
    """Count supports with exactly u even 10 pairs; invalid flows have count 0."""
    feasible = _checked(pairs, weight, first, second)
    if type(u) is not int:
        raise ValueError('Expected an integer even10 count')
    low, high = _bounds(weight, first, second)
    if not feasible or not low <= u <= high:
        return 0
    return _direct(pairs, weight, first, second, u)


def _evaluate(pairs, weight, first, second, *, record=False):
    """Positive half-sum; record=True exposes terms only for tiny tests."""
    low, high = _bounds(weight, first, second)
    singles = weight-2*first
    last = min(high, singles//2)
    total = previous = previous_tau = 0
    previous_four = False
    recurrences = direct = zeros = 0
    records = [] if record else None
    for u in range(low, last+1):
        v = singles-u
        e, f = first+u-second, first+v-second
        e0, e1 = pairs-first-u, first+u
        o0, o1 = pairs-first-v, first+v
        four = min(e0, e1, o0, o1) > 0
        tau = u*(pairs-weight+second)*o1 + first*f*o0 if four else 0
        if four and not tau:
            value = 0
            zeros += 1
        elif four and previous_four and previous_tau:
            # The previous flow has u-1. Its degrees are
            # e0+1,e1-1,o0-1,o1+1 and its variable edges v+1,f+1,u-1,e-1.
            numerator = previous*tau*(e1-1)*(o0-1)*(v+1)*(f+1)
            denominator = previous_tau*e0*o1*u*e
            value, remainder = divmod(numerator, denominator)
            if remainder:
                raise ArithmeticError('Nonintegral consecutive-flow recurrence')
            recurrences += 1
        else:
            value = _direct(pairs, weight, first, second, u)
            direct += 1
        multiplicity = 1 if 2*u == singles else 2
        total += multiplicity*value
        if record:
            records.append((u, multiplicity, value))
        previous, previous_tau, previous_four = value, tau, four
    diagnostics = (high-low+1, last-low+1, recurrences, direct, zeros)
    return total, diagnostics, tuple(records) if record else None


@lru_cache(maxsize=4096)
def _cached(pairs, weight, first, second):
    return _evaluate(pairs, weight, first, second)[:2]


def _canonical(pairs, weight, first, second):
    if weight > pairs:
        first, second = pairs-weight+first, pairs-weight+second
        weight = 2*pairs-weight
    if first < second:
        first, second = second, first
    return pairs, weight, first, second


def joint_count(pairs, weight, first, second):
    """Exact number of indexed binary supports with profile (r,a,b).

    Invalid matching marginals return 0. The cache retains at most 4096 scalar
    counts/diagnostics, never a profile-pair table. This function does not
    select profiles or silently substitute an upper bound on a timeout.
    """
    if not _checked(pairs, weight, first, second):
        return 0
    return _cached(*_canonical(pairs, weight, first, second))[0]


def count_diagnostics(pairs, weight, first, second):
    if not _checked(pairs, weight, first, second):
        return dict(count=0, feasible_marginals=False, evaluated_flow_terms=0)
    canonical = _canonical(pairs, weight, first, second)
    value, (flows, evaluated, recurrence, direct, zeros) = _cached(*canonical)
    return dict(count=value, feasible_marginals=True, canonical_profile=list(canonical),
                feasible_flow_terms=flows, evaluated_flow_terms=evaluated,
                recurrence_steps=recurrence, direct_initializations=direct,
                disconnected_four_vertex_terms=zeros, reflection_half_sum=True,
                indexed_not_rotation_classes=True, exact_integer_arithmetic=True)


def source_pins():
    path = Path(__file__).resolve()
    return {path.name: sha256(path.read_bytes()).hexdigest(),
            '../'+NOTE.name: sha256(NOTE.read_bytes()).hexdigest()}
