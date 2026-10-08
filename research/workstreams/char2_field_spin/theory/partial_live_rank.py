"""Rank moments when only a designated subset of coordinates must be live.

The caller must prove two hypotheses before exposing the fresh diagonals:
each designated-live vector satisfies the equations with probability at
most ``probability``, and every surviving kernel has dimension at most
``kernel_cap``. Merely knowing numerical matrix rank does not establish
the first hypothesis when some coordinates may be zero.
"""
from fractions import Fraction as F


def partial_live_moment(dimension, live_coordinates, kernel_cap, probability,
                        tau, q=256):
    d, p, cap = dimension, live_coordinates, kernel_cap
    if not (1 <= d and 1 <= p < q and 0 <= cap and 0 <= probability <= 1
            and 1 <= tau <= q-1):
        raise ValueError('invalid partial-live rank-moment parameters')
    if isinstance(tau, int):
        tau = F(tau)
    maximum = min(d, cap)
    if maximum == 0:
        return tau*0
    xs = [0]+[q**j-1-p*(q**(j-1)-1) for j in range(1, maximum+1)]
    ys = [tau*0]+[tau**(j-d) for j in range(1, maximum+1)]
    mean = (q**d-1)*probability
    if mean >= xs[-1]:
        return ys[-1]
    best = max(y for x, y in zip(xs, ys) if x <= mean)
    # One mean constraint and total mass one: an extremizer has at most
    # two nonzero masses. Enumerating pairs also handles a nonconcave
    # first segment near tau=q-1; assuming concavity would be unsound.
    for left in range(maximum):
        if xs[left] > mean:
            break
        for right in range(left+1, maximum+1):
            if xs[right] >= mean:
                fraction = (mean-xs[left])/(xs[right]-xs[left])
                best = max(best, ys[left]+fraction*(ys[right]-ys[left]))
    return best


def mds_return_bound(weight, state, dimension_cap, tau, selected_zeros=0, q=256):
    """Joint state-return/output-zero moment for the full-Cauchy packet.

    Pure returns (selected_zeros=0) only require superregular injection C.
    Selected output zeros require the full-Cauchy joint-rank hypothesis.
    No incoming zero-flat enumeration is used; only input support is fixed.
    """
    w, s, d, z = weight, state, dimension_cap, selected_zeros
    if not (1 <= w < q and 1 <= s and 1 <= d <= w+s and 0 <= z):
        raise ValueError('invalid return geometry')
    if z >= w:
        return tau*0
    probability = F(1, (q-1)**(s+z))
    return max(partial_live_moment(actual, w, w-z, probability, tau, q)
               for actual in range(1, d+1))
