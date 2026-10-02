"""A Fourier-monotone comparison, not actual shared-code shell counts.

For a binary linear output map and 0<z<=1, the Fourier expansion of
z^output_weight has nonnegative coefficients. A uniform nonzero q-ary
packet has nontrivial Fourier coefficient -1/(q-1). Replacing it by zero
with probability 2/q and by a uniform nonzero packet otherwise changes
that coefficient to +1/(q-1), increasing the output-weight moment.

All comparison coefficients are now nonnegative. Adding packet positions
can only decrease the moment. Uniform supports of consecutive sizes can
be coupled by inclusion, so an upper support CDF, with the correct final
mass, may replace the actual CDF by stochastic domination. Differences
of that CDF are used ONLY in this comparison, never as actual shell caps.
The initial state must be zero; arbitrary affine offsets would introduce
Fourier signs. Packet labels must be independent of the remaining maps.
"""
from fractions import Fraction as Q
from math import comb


def thinned_shells(cdf, q=16):
    """Return exact masses of active groups after the monotone comparison.

    The caller must establish the CDF upper bounds and exact final mass.
    Empty outcomes retain the nonzero-message/active-group label.
    """
    if (type(q) is not int or q < 4 or q & (q-1)
            or len(cdf) < 2 or cdf[0] != 0
            or any(type(v) is not int or v < 0 for v in cdf)
            or any(a > b for a, b in zip(cdf, cdf[1:]))):
        raise ValueError('monotone integer nonzero-group CDF and binary-extension alphabet required')
    n = len(cdf)-1
    keep, drop = Q(q-2, q), Q(2, q)
    result = [Q(0)]*(n+1)
    for u in range(1, n+1):
        mass = cdf[u]-cdf[u-1]
        if mass:
            for v in range(u+1):
                result[v] += mass*comb(u, v)*keep**v*drop**(u-v)
    assert sum(result) == cdf[-1]
    return result
