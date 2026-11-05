"""Exact inverse-refresh bound for a packet with no potential input positions.

For T=D2 H D1 with a totally nonsingular s-by-s H, T^-1 has the same
two-diagonal form around H^-1. Its output support is uniform conditional on
its weight. The weight transition below follows by shortening [I;H].
It needs neither independent codewords nor a uniform-refresh approximation.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb


@lru_cache(None)
def weight_transition(q=256,state=16):
    if q<2 or q&(q-1) or not 1<=state or 2*state>q:
        raise ValueError('Need a characteristic-two field and disjoint Cauchy nodes')
    result=[]
    for a in range(1,state+1):
        row=[]
        for b in range(1,state+1):
            dimension=a+b-state
            # Full support in a shortening of an MDS [state+a,a] code.
            count=0 if dimension<=0 else (q-1)*sum(
                (-1)**j*comb(a+b-1,j)*q**(dimension-1-j) for j in range(dimension))
            if count<0:raise ArithmeticError('Negative MDS support count')
            row.append(F(comb(state,b)*count,(q-1)**a))
        if sum(row)!=1:raise ArithmeticError('Weight transition does not normalize')
        result.append(tuple(row))
    return tuple(result)


@lru_cache(None)
def block_support_counts(state=16,h=4):
    if h<1 or state%h:raise ValueError('State must partition into equal blocks')
    counts={(0,0):1}
    for _ in range(state//h):
        following={}
        for (weight,active),count in counts.items():
            for j in range(h+1):
                key=(weight+j,active+bool(j))
                following[key]=following.get(key,0)+count*comb(h,j)
        counts=following
    result=tuple(tuple(counts.get((b,a),0) for a in range(state//h+1)) for b in range(state+1))
    for b,row in enumerate(result):assert sum(row)==comb(state,b)
    return result


def quiet_cap(z,*,q=256,state=16,width=64,h=4):
    """Rational cap on E_T[z^wt(A T^-1 y)] for every nonzero y.

    A is a direct sum of [width/(state/h),h] MDS encoders. The cap only
    uses their minimum distances, not a uniform output-weight hypothesis.
    """
    z=F(z)
    if not 0<z<=1 or state%h or width%(state//h) or width<state:
        raise ValueError('Invalid mark or MDS block geometry')
    floor=width//(state//h)-h+1
    local=tuple(sum(count*z**(floor*a) for a,count in enumerate(row))/comb(state,b)
                for b,row in enumerate(block_support_counts(state,h)))
    return max(sum(p*value for p,value in zip(row,local[1:])) for row in weight_transition(q,state))


@lru_cache(None)
def float_transition(q,state):
    import numpy as np
    return np.array([[float(p) for p in row] for row in weight_transition(q,state)])


def quiet_screen(z,*,q=256,state=16,width=64,h=4):
    """FLOAT proposal counterpart; directed replay must use quiet_cap."""
    import numpy as np
    floor=width//(state//h)-h+1
    local=np.array([sum(count*z**(floor*a) for a,count in enumerate(row))/comb(state,b)
                    for b,row in enumerate(block_support_counts(state,h))])
    return float((float_transition(q,state)@local[1:]).max())
