"""Actual uniform-matching moments, evaluated with directed integer rounding.

For a fixed s-vertex set, the l-th binomial moment counts l disjoint internal
edges. Evaluating its positive factorial-moment series retains the dependence
discarded by binomial domination. No construction or route changes.
"""
from fractions import Fraction as F
from functools import lru_cache
from rational_log_bounds import log2_bounds
from layered_matching_chernoff import matching_envelope,SCALE

TILT=F(3,2)


def moment_bounds(M,s,tilt=TILT,bits=96):
    """Return Q-scaled lower/upper bounds on E[tilt^Z], Q=2^bits."""
    assert M>=2 and M%2==0 and 0<=s<=M and tilt>=1 and bits>=16
    a,d=(tilt-1).numerator,(tilt-1).denominator
    Q=1<<bits;lo=hi=Q;lower=upper=Q
    for j in range(s//2):
        num=(s-2*j)*(s-2*j-1)*a
        den=2*(j+1)*(M-2*j-1)*d
        lo=lo*num//den;hi=(hi*num+den-1)//den
        lower+=lo;upper+=hi
    return lower,upper


@lru_cache(maxsize=8192)
def moment_log_upper(M,s,tilt=TILT):
    _,upper=moment_bounds(M,s,tilt)
    return log2_bounds(F(upper,1<<96),SCALE)[1]


@lru_cache(maxsize=16)
def factorial_envelope(M,R):
    assert M>=2 and M%2==0 and R>=1
    previous=matching_envelope(M,R);small=[];choose=1
    tilt_lower=log2_bounds(TILT,SCALE)[0]
    for s in range(M//2+1):
        if s:choose=choose*(M-s+1)//s
        numerator=R*moment_log_upper(M,s)+(choose.bit_length()+77)*SCALE
        threshold=(numerator+tilt_lower-1)//tilt_lower
        small.append(min(previous[s],threshold))
    return tuple(small[s] if s<=M//2 else R*(s-M//2)+small[M-s] for s in range(M+1))
