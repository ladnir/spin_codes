"""A count-dependent, pointwise majorant for the shuffle density loss.

For all profiles a summing to n, R(a) <= B * product_j tau_j**a_j.
Only the integer anchor counts are witnesses; replay derives B and tau.
"""
from fractions import Fraction as Q
from math import factorial,log,log1p,lgamma

from flint import arb


def validate(n,anchors):
    if type(n) is not int or n<1 or not anchors or any(type(k) is not int or not 0<=k<=n for k in anchors):
        raise ValueError('positive slot count and bounded integer anchors required')


def exact(n,anchors):
    """Exact constants for small exhaustive checks; avoid huge n here."""
    validate(n,anchors)
    ratios=[Q(k,k+1)**k if k else Q(1) for k in anchors]
    base=Q(n**n,factorial(n))*ratios[0]**n
    for k,r in zip(anchors,ratios):base*=Q(factorial(k),k**k)/r**k
    return base,tuple(r/ratios[0] for r in ratios)


def floating(n,anchors):
    validate(n,anchors)
    slopes=[-k*log1p(1/k) if k else 0. for k in anchors]
    logbase=n*log(n)-lgamma(n+1)+n*slopes[0]
    for k,s in zip(anchors,slopes):logbase+=lgamma(k+1)-(k*log(k) if k else 0.)-k*s
    return logbase,tuple(s-slopes[0] for s in slopes)


def outward(n,anchors):
    """Arb enclosures of the exact logarithms, including signed slopes."""
    validate(n,anchors)
    slopes=[k*(arb(k)/(k+1)).log() if k else arb(0) for k in anchors]
    logbase=n*arb(n).log()-arb(factorial(n)).log()+n*slopes[0]
    for k,s in zip(anchors,slopes):
        logbase+=arb(factorial(k)).log()-(k*arb(k).log() if k else arb(0))-k*s
    return logbase,tuple(s-slopes[0] for s in slopes)
