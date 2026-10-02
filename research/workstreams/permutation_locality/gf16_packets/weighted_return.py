"""Bound cancellation by the feedback moment tilted toward heavy inputs.

For nonzero a, z^wt(Aa+X) <= z^d_A z^-wt(X). Fourier inversion bounds
the tilted feedback event; no independence of feedback and weight is used.
"""
from math import comb
from fractions import Fraction as Q
import numpy as np
from flint import arb
import gf_refresh as base

prepare,actual,activity=base.prepare,base.actual,base.activity
aq,up=base.base.aq,base.base.up


def floating(data,probabilities,tilt):
    matrix=base.floating(data,probabilities,tilt)
    p=np.asarray(probabilities,dtype=float);z=np.exp(-tilt)
    S=1<<data['bits'];W=data['windows'];alpha=2.**-data['updates']
    kernel=base.base.base
    factors=kernel.KRAW@(p*z**-np.arange(5))
    values=kernel.float_products(factors,data['records'],W)
    mult=data['multiplicities'];empty=p[0]**W
    atom=min(float(mult@np.abs(values)),float(mult@np.abs(values-empty)))/S
    nonzero=max(0.,float(factors[0]**W-mult@values/S))
    floor=z**data['distance']
    levels=data['histograms']@np.arange(5)
    mean=float(data['histogram_multiplicities']@(z**levels))/(S-1)
    matrix[1,0]=min(matrix[1,0],alpha*floor*atom+matrix[1,2]/(S-1))
    matrix[2,0]=min(matrix[2,0],alpha*min(floor*nonzero/(S-1),atom*mean)+matrix[2,2]/(S-1))
    return matrix


def feedback_moments(data,probabilities,h):
    """Bound nonzero feedback atoms and mass weighted by h^wt(X)."""
    activity(probabilities)
    if not h>0:raise ValueError('positive input-weight factor required')
    p=list(map(aq,probabilities));W=data['windows'];S=1<<data['bits'];kernel=base.base.base
    factors=[sum((p[b]*h**b*sum((-1)**j*comb(r,j)*comb(4-r,b-j)
                    for j in range(max(0,b-4+r),min(b,r)+1))/comb(4,b)
                  for b in range(5)),arb(0)) for r in range(5)]
    values=kernel.arb_products(factors,[tuple(map(int,row)) for row in data['records']],W)
    mult=list(map(int,data['multiplicities']));empty=p[0]**W
    atom=min(up(sum((n*abs(v) for n,v in zip(mult,values)),arb(0))/S),
             up(sum((n*abs(v-empty) for n,v in zip(mult,values)),arb(0))/S))
    nonzero=up(factors[0]**W-sum((n*v for n,v in zip(mult,values)),arb(0))/S)
    if nonzero<0:raise ArithmeticError('negative tilted nonzero-feedback mass')
    return atom,nonzero


def outward_at_z(data,probabilities,z):
    matrix=base.outward_at_z(data,probabilities,z)
    atom,nonzero=feedback_moments(data,probabilities,1/z)
    S=1<<data['bits'];alpha=arb(2)**-data['updates']
    floor=z**data['distance']
    mean=up(sum((int(n)*z**sum(w*int(c) for w,c in enumerate(row))
                 for n,row in zip(data['histogram_multiplicities'],data['histograms'])),arb(0))/(S-1))
    matrix[1,0]=min(matrix[1,0],up(alpha*floor*atom+matrix[1,2]/(S-1)))
    matrix[2,0]=min(matrix[2,0],up(alpha*min(up(floor*nonzero/(S-1)),up(atom*mean))+matrix[2,2]/(S-1)))
    return matrix


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
