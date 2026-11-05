"""Sharper matching envelope from the already-proved binomial MGF bound.

The rational tilt is selected deterministically, before observing the route.
Every candidate has the same per-subset tail budget, so taking the smallest
threshold does not require a union over tilt candidates.
"""
from fractions import Fraction as F
from functools import lru_cache
from layered_route_gate import graph_envelope
from rational_log_bounds import log2_bounds

SCALE=1<<24
TILTS=(F(5,4),F(3,2),F(7,4),F(2))


def threshold(n,p,bits,tilt):
    assert n>=0 and 0<=p<=1 and tilt>1 and bits>=0
    lower,_=log2_bounds(tilt,SCALE)
    _,upper=log2_bounds(1+p*(tilt-1),SCALE)
    return min(n,(n*upper+bits*SCALE+lower-1)//lower)


@lru_cache(maxsize=16)
def matching_envelope(M,R):
    assert M>=2 and M%2==0 and R>=1
    old=graph_envelope(M,R);small=[];choose=1
    for s in range(M//2+1):
        if s:choose=choose*(M-s+1)//s
        n=R*(s//2);p=F(s-1 if s%2==0 else s,M-1)
        value=0 if n==0 else min(old[s],*(threshold(n,p,choose.bit_length()+77,t) for t in TILTS))
        small.append(value)
    return tuple(small[s] if s<=M//2 else R*(s-M//2)+small[M-s] for s in range(M+1))
