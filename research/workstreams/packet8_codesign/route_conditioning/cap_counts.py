"""Route bounds for H_c=sum_steps min(potential occupancy, c).

Fix q group identities. Independent regional permutations place them in
uniform q-subsets of ``slots`` positions. Each physical step has ``windows``
positions. The bad event is that some q-subset has H_c < h across all regions.
Only group subsets enter this union; message labels do not.

Rational witnesses use exact integer polynomial coefficients. Floating
moments are search aids, not certified probability endpoints.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import hashlib
from math import comb, exp, isfinite, log
from pathlib import Path


def _geometry(q, cap, slots, windows):
    if (any(type(value) is not int for value in (q, cap, slots, windows))
            or slots <= 0 or windows <= 0 or slots % windows
            or not 0 <= q <= slots or not 1 <= cap <= windows):
        raise ValueError('integer q in0..slots, cap in1..windows and complete positive steps required')
    return slots//windows


def _threshold(q, h, cap, regions, slots, windows):
    steps = _geometry(q, cap, slots, windows)
    if (type(regions) is not int or regions <= 0 or type(h) is not int
            or not 0 <= h <= regions*steps*cap+1):
        raise ValueError('positive regions and threshold in0..total_capped_capacity+1 required')
    minimum = regions*((q//windows)*cap+min(q % windows, cap))
    maximum = regions*min(q, cap*steps)
    return steps, minimum, maximum


def _rational(x):
    if isinstance(x, float):
        raise ValueError('use a Fraction or rational string, not a floating tilt')
    x = Fraction(x)
    if not 0 < x <= 1:
        raise ValueError('rational Chernoff variable must be in(0,1]')
    return x


@lru_cache(maxsize=32)
def _weighted_count(q, u, v, cap, slots, windows):
    """[t^q](sum_j C(windows,j)u^min(j,c)v^(c-min(j,c))t^j)^steps."""
    weights = [comb(windows, j)*u**min(j, cap)*v**(cap-min(j, cap))
               for j in range(min(windows, q)+1)]
    coefficients = [0]*(q+1)
    coefficients[0] = 1
    for step in range(slots//windows):
        following = [0]*(q+1)
        for degree in range(min(q, (step+1)*windows)+1):
            following[degree] = sum(coefficients[degree-j]*weights[j]
                for j in range(min(windows, degree, len(weights)-1)+1))
        coefficients = following
    return coefficients[q]


def one_region_moment(q, x, *, cap=2, slots=512, windows=8):
    """Exact expectation of x^H_c for one region and a fixed q-subset."""
    steps = _geometry(q, cap, slots, windows)
    x = _rational(x)
    value = _weighted_count(q, x.numerator, x.denominator, cap, slots, windows)
    return Fraction(value, comb(slots, q)*x.denominator**(cap*steps))


def moment_and_mean(q, eta, *, cap=2, slots=512, windows=8):
    """Floating (log E exp(-eta H_c), tilted E H_c) for one region."""
    steps = _geometry(q, cap, slots, windows)
    eta = float(eta)
    if not isfinite(eta) or eta < 0:
        raise ValueError('finite nonnegative Chernoff parameter required')
    weights = [comb(windows, j)*exp(-eta*min(j, cap))
               for j in range(min(windows, q)+1)]
    coefficients, moments = [0.]*(q+1), [0.]*(q+1)
    coefficients[0] = 1.
    accumulated_scale = 0.
    for step in range(steps):
        following, next_moments = [0.]*(q+1), [0.]*(q+1)
        for degree in range(min(q, (step+1)*windows)+1):
            for j in range(min(windows, degree, len(weights)-1)+1):
                following[degree] += weights[j]*coefficients[degree-j]
                next_moments[degree] += weights[j]*(moments[degree-j]
                    +min(j, cap)*coefficients[degree-j])
        scale = max(following)
        if not isfinite(scale) or scale <= 0:
            raise FloatingPointError('nonpositive capped-occupancy polynomial scale')
        coefficients = [x/scale for x in following]
        moments = [x/scale for x in next_moments]
        accumulated_scale += log(scale)
    if coefficients[q] <= 0 or not isfinite(coefficients[q]):
        raise FloatingPointError('capped-occupancy coefficient underflowed')
    return (log(coefficients[q])+accumulated_scale-log(comb(slots, q)),
            moments[q]/coefficients[q])


def bad_route_margin(q, h, eta, *, cap=2, regions=32, slots=512, windows=8):
    """Floating margin for the subset union Pr[exists S: H_c(S) < h]."""
    _, minimum, maximum = _threshold(q, h, cap, regions, slots, windows)
    eta = float(eta)
    if not isfinite(eta) or eta < 0:
        raise ValueError('finite nonnegative Chernoff parameter required')
    if h <= minimum:
        return float('inf')
    if h > maximum:
        return 0.
    log_moment, _ = moment_and_mean(q, eta, cap=cap, slots=slots, windows=windows)
    return -(log(comb(slots, q))+eta*(h-1)+regions*log_moment)/log(2)


def rational_bad_route_bound(q, h, x, *, cap=2, regions=32, slots=512, windows=8):
    """Exact rational Chernoff upper bound, capped at one.

    Write x=u/v and let T be the integer coefficient in _weighted_count.
    The group-subset union is T^regions /
    (C(slots,q)^(regions-1) v^(cap*steps*regions-h+1) u^(h-1)).
    """
    steps, minimum, maximum = _threshold(q, h, cap, regions, slots, windows)
    x = _rational(x)
    if h <= minimum:
        return Fraction(0)
    if h > maximum:
        return Fraction(1)
    u, v = x.numerator, x.denominator
    value = _weighted_count(q, u, v, cap, slots, windows)
    numerator = value**regions
    denominator = (comb(slots, q)**(regions-1)
        *v**(cap*steps*regions-h+1)*u**(h-1))
    return min(Fraction(1), Fraction(numerator, denominator))


def rational_witness(q, h, x, *, cap=2, target_bits=60, regions=32, slots=512, windows=8):
    """Compact witness with an exact integer pass check and source pin."""
    if type(target_bits) is not int or target_bits < 0:
        raise ValueError('nonnegative integer target bits required')
    bound = rational_bad_route_bound(q, h, x, cap=cap, regions=regions, slots=slots, windows=windows)
    if bound:
        certified_bits = bound.denominator.bit_length()-bound.numerator.bit_length()
        if bound.numerator << certified_bits > bound.denominator:
            certified_bits -= 1
        passes = bound.numerator << target_bits <= bound.denominator
    else:
        certified_bits, passes = 'infinite', True
    source = Path(__file__).resolve()
    return dict(q=q, h=h, cap=cap, statistic='sum_steps min(potential_occupancy,cap)',
        event='some q-group subset has H_cap < h', slots=slots, windows=windows,
        regions=regions, chernoff_x=str(_rational(x)), union_over_group_subsets=True,
        union_over_message_labels=False, target_bits=target_bits,
        certified_margin_bits_floor=certified_bits, exact_integer_check_passed=passes,
        whole_code_certificate=False,
        source_sha256={str(source): hashlib.sha256(source.read_bytes()).hexdigest()})
