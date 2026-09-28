"""Average iid packet weights before computing distinct-window moments.

This improves the high-occupancy fallback of the homogeneous reference
model, without enumerating every packet-weight shape. It does not refine
feedback cancellation or assert uniform state after output tilting.
"""
from fractions import Fraction as Q
from math import comb

import numpy as np
from flint import arb, ctx

from shape_inner import rational
from window_histogram import lane_factor
from occupancy_memory import M, F, U, rounded
from mature_tail import L48, L56
from group_rank_one_verify import up


def histogram_moments(hist, theta, z, rounding=lambda x: x):
    """All occupancies for one expansion-window histogram.

    At j occupied windows, average over uniform distinct windows and iid
    nonzero packet weights theta. Positive coefficients support exact
    rationals, Arb, or diagnostic floats through the same recurrence.
    """
    if len(hist) != 5 or any(not isinstance(n, int) or n < 0 for n in hist):
        raise ValueError('five nonnegative window counts required')
    if len(theta) != 4:
        raise ValueError('four nonzero packet weights required')
    factors = [rounding(sum((theta[b-1]*lane_factor(r,b,z) for b in range(1,5)), z*0))
               for r in range(5)]
    values = [z*0+1]
    for r, number in enumerate(hist):
        for _ in range(number):
            old = values
            values = [old[0]]+[rounding(old[j]+factors[r]*old[j-1]) for j in range(1,len(old))]
            values.append(rounding(factors[r]*old[-1]))
    weight = sum(r*n for r,n in enumerate(hist))
    base = rounding(z**weight)
    return [rounding(base*value/comb(sum(hist),j)) for j,value in enumerate(values)]


def refine(operators, data, theta, tilt, cutoff, *, exact=True):
    """Refine only coarse_lift fallback slots j>cutoff; two updates.

    The fallback uses M alone for outgoing mature mass, and each low-weight
    tail is bounded by that same mass. This routine is not applicable to
    arbitrary operators whose mature bound is split over tail coordinates.
    """
    histograms, numbers, fresh, spectrum = data
    levels = sorted(spectrum)
    size = (1 << 19)-1
    if exact:
        theta = tuple(Q(p) for p in theta)
        if len(theta) != 4 or min(theta) < 0 or sum(theta) != 1:
            raise ValueError('an exact rational packet probability law is required')
        probabilities = [rational(p) for p in theta]
        z = (-arb(tilt)).exp()
        rounding = up
        alpha, beta, zero = arb(1)/4, arb(3)/4, arb(0)
    else:
        probabilities = np.asarray(theta, dtype=float)
        if probabilities.shape != (4,) or min(probabilities) < 0 or abs(sum(probabilities)-1) > 1e-12:
            raise ValueError('invalid proposal law')
        z = np.exp(-float(tilt))
        rounding = lambda x: x
        alpha, beta, zero = .25, .75, 0.
    moments = {h:histogram_moments(hist, probabilities, z, rounding)
               for h,hist in enumerate(histograms) if sum(r*n for r,n in enumerate(hist))}
    classes = {v:[h for h,hist in enumerate(histograms) if sum(r*n for r,n in enumerate(hist)) == v]
               for v in levels}
    result = [t*1 for t in operators]
    for j in range(cutoff+1,len(result)):
        t = result[j]
        if t[L48,M]!=0 or t[L56,M]!=0:
            raise ValueError('this refinement requires the unsplit coarse mature-mass fallback')
        arbitrary = min(1,max(row[j] for row in moments.values()))
        fresh_bound = min(1,max(rounding(sum((n*moments[h][j] for h,n in distribution.items() if n),zero)
                                             /sum(distribution.values())) for distribution in fresh))
        uniform = {v:min(1,rounding(sum((numbers[h]*moments[h][j] for h in classes[v]),zero)/spectrum[v]))
                   for v in levels}
        values = {M:arbitrary,F:fresh_bound,**{U+k:uniform[v] for k,v in enumerate(levels)}}
        for source,value in values.items():
            # This mass bound includes nonzero targets and also zero targets;
            # retaining the latter only enlarges the mature upper bound.
            t[source,M] = min(t[source,M],rounding(alpha*value))
            for k,v in enumerate(levels):
                t[source,U+k] = min(t[source,U+k],rounding(beta*value*spectrum[v]/size))
            for tail in (L48,L56):
                t[source,tail] = min(t[source,tail],t[source,M])
        t[M,0] = min(t[M,0],rounding(beta*arbitrary/size))
        if exact:
            result[j] = rounded(t)
    return result if exact else np.asarray(result)


class AveragedHighInner:
    """Opt-in wrapper; exact and float paths share the coefficient formula."""
    def __init__(self, base):
        if base.shared.histograms is None:
            raise ValueError('the exact expansion-window census is required')
        self.base = base

    def __getattr__(self, name):
        return getattr(self.base,name)

    def mix(self, theta):
        ctx.prec = self.base.precision
        return refine(self.base.mix(theta),self.base.shared.histograms,theta,
                      self.base.tilt,self.base.cut)

    def float_mix(self, theta):
        return refine(self.base.float_mix(theta),self.base.shared.histograms,theta,
                      self.base.tilt,self.base.cut,exact=False)
