"""Scalar density bounds for independently shuffled BCH rows.

Let A[w] count row messages whose encoded word has weight w. After a uniform
coordinate permutation, each support of that weight has mass A[w]/C(n,w).
For weights in an interval I and a Bernoulli(p) reference measure, this mass
is pointwise at most gamma_I(p) times the reference mass, where

    gamma_I(p) = max_{w in I, Abar[w]>0} Abar[w]/beta(n,w,p).

Here Abar is a coefficientwise spectrum cap and beta is the binomial mass.
The maximum, not a sum, suffices because different weights have disjoint
supports. Consequently the bound applies to every nonnegative downstream
function, including an inner-encoder moment. Excluding weight zero must be
explicit. A row forced to zero or all ones can instead use p=0 or p=1.

For a four-row group with r active rows and disjoint/equal intervals I_i,
the bound for all row-label assignments is

    group_box_multiplicity(I_i) * product_i gamma_{I_i}(p_i)

times the corresponding reference expectation. Inactive rows remain zero.
No coordinate-permutation or lane-permutation multiplicity is added: those
are probability measures already included in the reference experiment.

Heterogeneous groups may be dominated by one common packet reference after
the fresh uniform lane shuffle. packet_reference_factor gives a one-region
density factor d_g against Binomial(R,p*). Multiply by d_g**n for n regions,
in addition to the row factors. This can be expensive. At equal probabilities
p and r<=R, d_g=(1-p)**(r-R), so the extra cost is (1-p)**(-n*(R-r)).
Homogeneous-r analysis avoids that extra comparison. These scalar helpers
do not establish an inner bound or a complete occupancy cover.
"""
from collections import Counter
from fractions import Fraction as Q
from math import comb, factorial, log

import numpy as np
from flint import arb
from scipy.special import gammaln


def _inputs(caps, lo, hi, exclude_zero):
    caps = tuple(Q(x) for x in caps)
    if not caps or any(x < 0 for x in caps):
        raise ValueError('nonempty nonnegative spectrum caps required')
    n = len(caps) - 1
    if not isinstance(lo, int) or not isinstance(hi, int) or not 0 <= lo <= hi <= n:
        raise ValueError('invalid row-weight interval')
    weights = tuple(w for w in range(max(lo, 1 if exclude_zero else 0), hi + 1) if caps[w])
    return caps, n, weights


def _probability(value):
    value = Q(value)
    if not 0 <= value <= 1:
        raise ValueError('probability must lie in [0,1]')
    return value


def _log_exact(value):
    return log(value.numerator) - log(value.denominator) if value else -np.inf


def _arb_exact(value):
    return arb(value.numerator) / value.denominator


def _up(value):
    assert value.is_finite()
    return arb(value.upper())


def forced_reference(caps, lo, hi, *, exclude_zero=False):
    """Return 0/1 for a forced endpoint row, or None otherwise.

    Empty classes return None and have density factor zero for every p.
    """
    _, n, weights = _inputs(caps, lo, hi, exclude_zero)
    if weights == (0,):
        return Q(0)
    if weights == (n,):
        return Q(1)
    return None


def row_gamma_exact(caps, lo, hi, p, *, exclude_zero=False):
    """Exact rational scalar domination for a row-weight interval.

    Degenerate p=0/1 is allowed only when all positive capped mass lies at
    that endpoint. This excludes the undefined positive-mass/zero-reference
    ratio rather than silently dropping the unsupported weight.
    """
    caps, n, weights = _inputs(caps, lo, hi, exclude_zero)
    p = _probability(p)
    ratios = []
    for w in weights:
        mass = comb(n, w) * p**w * (1 - p)**(n - w)
        if not mass:
            raise ValueError('reference gives zero mass to a permitted positive shell')
        ratios.append(caps[w] / mass)
    return max(ratios, default=Q(0))


def row_gamma_arb(caps, lo, hi, p, *, exclude_zero=False):
    """Outward Arb upper endpoint with an exact rational witness p."""
    caps, n, weights = _inputs(caps, lo, hi, exclude_zero)
    p = _probability(p)
    if p in (0, 1):
        return _up(_arb_exact(row_gamma_exact(caps, lo, hi, p, exclude_zero=exclude_zero)))
    probability = _arb_exact(p)
    result = arb(0)
    for w in weights:
        mass = arb(comb(n, w)) * probability**w * (1 - probability)**(n - w)
        lower = arb(mass.lower())
        if not lower > 0:
            raise ArithmeticError('precision is insufficient for a positive binomial-mass lower bound')
        result = max(result, _up(_arb_exact(caps[w]) / lower))
    return result


def row_gamma_function(caps, lo, hi, *, exclude_zero=False):
    """Precompute the binary64 log-gamma proposal objective; not a bound.

    Endpoint references with incompatible support return +infinity, allowing
    optimizers to reject them. Empty classes return -infinity.
    """
    caps, n, weights = _inputs(caps, lo, hi, exclude_zero)
    supports = np.asarray(weights, dtype=np.int64)
    constants = np.asarray([_log_exact(caps[w]) for w in weights])
    constants -= gammaln(n + 1) - gammaln(supports + 1) - gammaln(n - supports + 1)
    def evaluate(p):
        if not 0 <= p <= 1:
            raise ValueError('probability must lie in [0,1]')
        if not weights:
            return -np.inf
        if p == 0:
            return _log_exact(caps[0]) if weights == (0,) else np.inf
        if p == 1:
            return _log_exact(caps[n]) if weights == (n,) else np.inf
        return float(np.max(constants - supports * np.log(p) - (n - supports) * np.log1p(-p)))
    return evaluate


def row_gamma_log(caps, lo, hi, p, *, exclude_zero=False):
    return row_gamma_function(caps, lo, hi, exclude_zero=exclude_zero)(p)


def group_box_multiplicity(intervals, *, rows=4):
    """Count row-label assignments of an unordered active-interval box.

    Intervals describe positive row weights. Equal intervals are identical
    labels; distinct intervals must be disjoint for an exact partition count.
    Inactive rows all receive the single zero label. The empty box represents
    an inactive group, which must not enter an active-group occupancy class.
    """
    intervals = tuple(tuple(interval) for interval in intervals)
    if not isinstance(rows, int) or rows < 1 or len(intervals) > rows:
        raise ValueError('positive row count and at most that many intervals required')
    if any(len(interval) != 2 or not all(isinstance(x, int) for x in interval)
           or not 1 <= interval[0] <= interval[1] for interval in intervals):
        raise ValueError('active intervals must contain positive integer weights')
    distinct = sorted(set(intervals))
    if any(left[1] >= right[0] for left, right in zip(distinct, distinct[1:])):
        raise ValueError('distinct intervals must be disjoint')
    denominator = factorial(rows - len(intervals))
    for count in Counter(intervals).values():
        denominator *= factorial(count)
    return factorial(rows) // denominator


def packet_reference_factor(probabilities, reference_rows, reference_probability, *, lanes=4):
    """Exact one-region domination after independent uniform lane shuffling.

    probabilities lists Bernoulli parameters for the non-fixed-zero rows.
    The target uses reference_rows Bernoulli(reference_probability) rows.
    Remaining lanes are zero in both experiments. Conditional on its weight,
    the lane-shuffled mask is uniform, so the ratio of weight probabilities
    also bounds every mask. Across n independent regions the factor is d**n.
    """
    probabilities = tuple(_probability(p) for p in probabilities)
    p = _probability(reference_probability)
    if (not isinstance(lanes, int) or lanes < 1
            or not isinstance(reference_rows, int) or not 0 <= reference_rows <= lanes
            or len(probabilities) > lanes):
        raise ValueError('invalid lane count or reference row count')
    source = [Q(1)]
    for probability in probabilities:
        next_source = [Q(0)] * (len(source) + 1)
        for w, mass in enumerate(source):
            next_source[w] += mass * (1 - probability)
            next_source[w + 1] += mass * probability
        source = next_source
    result = Q(0)
    for w, mass in enumerate(source):
        if not mass:
            continue
        reference = (comb(reference_rows, w) * p**w * (1 - p)**(reference_rows - w)
                     if w <= reference_rows else Q(0))
        if not reference:
            raise ValueError('packet reference has zero mass on a permitted weight')
        result = max(result, mass / reference)
    assert result >= 1
    return result
