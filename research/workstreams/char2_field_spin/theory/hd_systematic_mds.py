"""Exact maps and positive moments for a systematic global-MDS inner.

This is a construction component, not by itself a whole-code certificate.
The field is the AES quadratic tower used by the native GF65536 kernel.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from hd_blockwise_shared_transfer import multiply,power
from hd_sumrank_maps import Q,transpose


def novel(x,j):
    value=1
    for bit in range(j.bit_length()):
        if j>>bit&1:
            top=bottom=1
            for h in range(1<<bit):
                top=multiply(top,x^h,Q)
                bottom=multiply(bottom,(1<<bit)^h,Q)
            value=multiply(value,multiply(top,power(bottom,Q-2,Q),Q),Q)
    return value


@lru_cache(None)
def maps(state=16,width=64):
    if width not in (4,8,16,32,64,128,256) or not 1<=state<=width//2:
        raise ValueError('Invalid self-orthogonal additive RS dimensions')
    a=tuple(tuple(novel(x,j) for j in range(state)) for x in range(width))
    return a,transpose(a)


def moments(j,z,state=16,width=64,q=Q):
    """Exact K,T,V for j potential coordinates and an MDS emission map."""
    z=F(z)
    if not 0<=j<=width or not 0<z<1 or not 1<=state<=width:
        raise ValueError('Invalid moment parameters')
    k=sum(F(comb(j,c))*z**(j-c)*(1-z)**c/q**(c+min(state,j-c))
          for c in range(j+1))
    t=(z+(1-z)/q)**j
    v=t*sum(F(comb(width-j,d))*z**(width-j-d)*(1-z)**d
            *F(q**(state-d)-1,q**state-1) for d in range(min(state,width-j+1)))
    assert 0<k<=t<=1 and 0<v<=1
    return k,t,v


def matrix(j,z,state=16,width=64,q=Q):
    k,t,v=moments(j,z,state,width,q)
    gamma=min(k,F(z)**max(0,width+1-j-state)/q**min(state,j),(q**state-1)*v)
    return k,t-k,v/(q**state-1) if j else F(0),(v if j else F(q-1,q)*v)+gamma/q


def digest(state=16):
    a,c=maps(state);f=tuple(tuple(int(i==j) for j in range(64)) for i in range(64))
    h=0
    for m in (f,a,c):
        for row in m:
            for x in row:h=((h^x)*0x100000001b3)&((1<<64)-1)
    return h
