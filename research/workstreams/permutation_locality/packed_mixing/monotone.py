"""Pointwise expected-shell bounds from cumulative canonical-block counts.

This is not Fourier thinning. The result bounds each expected shell itself,
so a positive shell majorant may use different output tilts in later cells.
No whole-code or distance claim is made here.
"""
from fractions import Fraction as Q

from local_models import checked_pmf, convolve, full_block


def _cdf(caps):
    caps = tuple(caps)
    if (len(caps) < 2 or caps[0] != 0
            or any(type(v) not in (int, Q) or v < 0 for v in caps)
            or any(a > b for a, b in zip(caps, caps[1:]))):
        raise ValueError('nonnegative rational nonzero-message CDF caps required')
    return tuple(map(Q, caps))


def suffix_expectation(caps, values):
    """Bound sum_h a[h]*values[h] given prefix bounds on a[1:].

    The caller establishes sum_{j<=h} a[j] <= caps[h], a[0]=0, and
    nonnegative a. Neither the counts nor the caps need to be integral.
    Let m[h]=max(values[h:]). Abel summation gives
    sum a[h]*values[h] <= sum a[h]*m[h]
      <= sum (caps[h]-caps[h-1])*m[h].
    This is sharp over nonnegative real counts with these prefix bounds:
    move each cap increment to any later index attaining its suffix maximum.
    """
    caps, values = _cdf(caps), tuple(values)
    if (len(values) != len(caps)
            or any(type(v) not in (int, Q) or v < 0 for v in values)):
        raise ValueError('matching nonnegative rational kernel values required')
    result, suffix = Q(0), Q(0)
    for h in range(len(caps)-1, 0, -1):
        suffix = max(suffix, values[h])
        result += (caps[h]-caps[h-1])*suffix
    return result


def transport_shells(caps, local=full_block(8)):
    """Return pointwise expected-shell caps after independent local maps.

    caps[h] bounds the number (or expected number) of nonzero messages with
    at most h active blocks. local is the common support PMF of an active
    block; it must assign zero probability to support zero. For output
    support w, use suffix_expectation with values[h]=Pr[sum of h iid W=w].

    Each result is at most transport_cdf(caps, local)[w]: for j>=h, coupling
    by adding positive W gives Pr[sum_j W=w] <= Pr[sum_h W<=w]. The latter
    is decreasing in h, so it also dominates the suffix maximum. The result
    has zero mass at zero, but its sum need not equal the message count.
    A fresh independent uniform coordinate shuffle and conditional iid
    nonzero packet labels are still required for pointwise vector domination.
    """
    caps, local = _cdf(caps), tuple(local)
    if any(type(v) not in (int, Q) for v in local):
        raise ValueError('exact rational local probabilities required')
    local = checked_pmf(local)
    if len(local) < 2 or local[0] != 0:
        raise ValueError('nonzero local support required')
    laws = [(Q(1),)]
    for _ in range(1, len(caps)):
        laws.append(convolve(laws[-1], local))
    size = (len(caps)-1)*(len(local)-1)
    result, suffix = [Q(0)]*(size+1), [Q(0)]*(size+1)
    for h in range(len(caps)-1, 0, -1):
        for w, probability in enumerate(laws[h]):
            suffix[w] = max(suffix[w], probability)
        difference = caps[h]-caps[h-1]
        if difference:
            for w, probability in enumerate(suffix):
                result[w] += difference*probability
    if result[0] != 0:
        raise ArithmeticError('nonzero blocks cannot produce empty support')
    return tuple(result)
