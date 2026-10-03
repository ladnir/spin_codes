"""Exact counts and Chernoff bounds for universal regional dispersion.

Fix q group identities among ``slots`` groups. Each independent regional
permutation places those identities in a uniform q-subset of ``slots``
positions. A physical inner step contains ``windows`` consecutive positions.
Let H count occupied physical steps across all regions. The bad routing event
is that some q-subset of group identities has H < h.

The union factor counts group subsets only, not message labels. This module
certifies that routing event; it does not certify the complete code.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import hashlib
from math import comb, exp, isfinite, log
from pathlib import Path


def _geometry(q, slots, windows):
    if (any(type(value) is not int for value in (q, slots, windows))
            or slots <= 0 or windows <= 0 or slots % windows or not 0 <= q <= slots):
        raise ValueError('integer 0 <= q <= slots and complete positive physical steps required')
    return slots // windows


def _threshold(q, h, regions, slots, windows):
    steps = _geometry(q, slots, windows)
    if (type(regions) is not int or regions <= 0 or type(h) is not int
            or not 0 <= h <= regions*steps+1):
        raise ValueError('positive region count and integer threshold in 0..total_steps+1 required')
    return steps


@lru_cache(maxsize=32)
def counts(q, slots=512, windows=8):
    """Number of q-subsets occupying exactly b steps, indexed by b.

    Entry b equals C(steps,b) [X^q] ((1+X)^windows-1)^b.
    Its denominator as a probability is C(slots,q).
    """
    steps = _geometry(q, slots, windows)
    coefficients = [0]*(q+1)
    coefficients[0] = 1
    choose = [comb(windows, k) for k in range(min(windows, q)+1)]
    result = []
    for b in range(steps+1):
        if b:
            following = [0]*(q+1)
            for degree in range(b, min(q, b*windows)+1):
                following[degree] = sum(
                    coefficients[degree-k]*choose[k]
                    for k in range(1, min(windows, degree, len(choose)-1)+1))
            coefficients = following
        result.append(comb(steps, b)*coefficients[q])
    if sum(result) != comb(slots, q):
        raise ArithmeticError('occupied-step counts fail the exact subset partition')
    return tuple(result)


one_region_counts = counts


def bad_route_margin(q, h, eta, *, regions=32, slots=512, windows=8):
    """Floating Chernoff margin, including the union over group subsets.

    The event is H <= h-1. This function is only a numerical proposal;
    ``rational_bad_route_bound`` provides an exact rational alternative.
    """
    _threshold(q, h, regions, slots, windows)
    eta = float(eta)
    if not isfinite(eta) or eta < 0:
        raise ValueError('finite nonnegative Chernoff parameter required')
    distribution = counts(q, slots, windows)
    support = [b for b, value in enumerate(distribution) if value]
    if h <= regions*support[0]:
        return float('inf')
    if h > regions*support[-1]:
        return 0.
    denominator = comb(slots, q)
    terms = [log(distribution[b])-eta*b for b in support]
    top = max(terms)
    log_moment = top+log(sum(exp(value-top) for value in terms))-log(denominator)
    return -(log(denominator)+eta*(h-1)+regions*log_moment)/log(2)


def rational_bad_route_bound(q, h, x, *, regions=32, slots=512, windows=8):
    """Exact upper bound for the bad routing event, using rational 0 < x <= 1.

    Set x=u/v and S=sum_b counts[b] u^b v^(steps-b). Markov's inequality
    and the subset union give S^regions /
    (C(slots,q)^(regions-1) v^(steps*regions-h+1) u^(h-1)).
    The result is capped at one. Integer arithmetic includes every term.
    """
    steps = _threshold(q, h, regions, slots, windows)
    if isinstance(x, float):
        raise ValueError('use a Fraction or rational string, not a floating tilt')
    x = Fraction(x)
    if not 0 < x <= 1:
        raise ValueError('rational Chernoff variable must be in (0,1]')
    distribution = counts(q, slots, windows)
    support = [b for b, value in enumerate(distribution) if value]
    if h <= regions*support[0]:
        return Fraction(0)
    if h > regions*support[-1]:
        return Fraction(1)
    u, v = x.numerator, x.denominator
    total = sum(value*u**b*v**(steps-b) for b, value in enumerate(distribution))
    numerator = total**regions
    denominator = (comb(slots, q)**(regions-1)
                   *v**(steps*regions-h+1)*u**(h-1))
    return min(Fraction(1), Fraction(numerator, denominator))


def rational_witness(q, h, x, *, target_bits=60, regions=32, slots=512, windows=8):
    """Compact metadata whose pass flag comes from an exact integer comparison."""
    if type(target_bits) is not int or target_bits < 0:
        raise ValueError('nonnegative integer target bits required')
    bound = rational_bad_route_bound(q, h, x, regions=regions, slots=slots, windows=windows)
    if bound:
        certified_bits = bound.denominator.bit_length()-bound.numerator.bit_length()
        if bound.numerator << certified_bits > bound.denominator:
            certified_bits -= 1
        passes = bound.numerator << target_bits <= bound.denominator
    else:
        certified_bits, passes = 'infinite', True
    source = Path(__file__).resolve()
    return dict(q=q, h=h, event='some q-group subset has H < h',
        slots=slots, windows=windows, regions=regions, chernoff_x=str(Fraction(x)),
        union_over_group_subsets=True, union_over_message_labels=False,
        target_bits=target_bits, certified_margin_bits_floor=certified_bits,
        exact_integer_check_passed=passes, whole_code_certificate=False,
        source_sha256={str(source): hashlib.sha256(source.read_bytes()).hexdigest()})
