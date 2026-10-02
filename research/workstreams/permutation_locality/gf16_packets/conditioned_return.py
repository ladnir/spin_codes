"""Condition on packet occupancy before applying Holder to a return event."""
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb
import shape_return as base
import occupancy_kernel as fixed

prepare,actual=base.prepare,base.actual
aq,up=base.aq,base.up
POWERS=(2,3,4,6,8,12,16)


def float_profiles(data,z):
    hist=data['histograms'];W=data['windows']
    # Batch all profiles through their W factors instead of looping over
    # the profiles separately for every proposed Chernoff tilt.
    patterns=np.array([np.repeat(np.arange(5),h) for h in hist])
    inactive=z**patterns;active=((1+z)**4-inactive)/15
    values=np.ones((len(hist),1))
    for j in range(W):
        new=np.zeros((len(hist),j+2))
        new[:,:-1]+=inactive[:,j,None]*values
        new[:,1:]+=active[:,j,None]*values
        values=new
    return values/np.array([comb(W,q) for q in range(W+1)])


def float_moments(data,z):
    values=float_profiles(data,z)
    return values.max(axis=0),data['histogram_multiplicities']@values/((1<<data['bits'])-1)


def float_bounds(data,p,z):
    W=data['windows'];L=(1<<data['bits'])-1
    atom=np.array(list(map(float,data['nonzero_atom_caps'])))
    nu=1-np.array(list(map(float,data['zero_probabilities'])))
    floor=z**np.maximum(0,data['distance']-4*np.arange(W+1))
    maximum=atom*floor;mean=nu*floor/L
    levels=data['histograms']@np.arange(5)
    per_state=atom[None,:]*z**np.maximum(0,levels[:,None]-4*np.arange(W+1)[None,:])
    rows=float_profiles(data,z)
    h=rows.max(axis=0);hbar=data['histogram_multiplicities']@rows/L
    maximum=np.minimum(maximum,h);mean=np.minimum(mean,hbar)
    per_state=np.minimum(per_state,rows)
    for power in POWERS:
        rows=float_profiles(data,z**power)
        h=rows.max(axis=0);hbar=data['histogram_multiplicities']@rows/L
        maximum=np.minimum(maximum,h**(1/power)*atom**(1-1/power))
        mean=np.minimum(mean,hbar**(1/power)*(nu/L)**(1-1/power))
        per_state=np.minimum(per_state,rows**(1/power)*atom[None,:]**(1-1/power))
    maximum=np.minimum(maximum,per_state.max(axis=0))
    mean=np.minimum(mean,data['histogram_multiplicities']@per_state/L)
    weights=np.array([comb(W,q)*p**q*(1-p)**(W-q) for q in range(W+1)])
    return float(weights@maximum),float(weights@mean)


def bounds(data,p,z):
    p=Q(p);W=data['windows'];L=(1<<data['bits'])-1
    if not 0<=p<=1 or not 0<z<=1:raise ValueError('probability and output weight required')
    atoms=list(map(aq,data['nonzero_atom_caps']));nu=[aq(1-x) for x in data['zero_probabilities']]
    floors=[z**max(0,data['distance']-4*q) for q in range(W+1)]
    h,hbar=fixed.moments(data,z)
    maximum=[min(v,up(a*f)) for v,a,f in zip(h,atoms,floors)]
    mean=[min(v,up(n*f/L)) for v,n,f in zip(hbar,nu,floors)]
    per_state=[]
    for histogram in data['histograms']:
        level=sum(w*int(n) for w,n in enumerate(histogram))
        row=fixed.polynomial(histogram,z)
        per_state.append([min(up(v),up(a*z**max(0,level-4*q))) for q,(v,a) in enumerate(zip(row,atoms))])
    for power in POWERS:
        rows=[fixed.polynomial(histogram,z**power) for histogram in data['histograms']]
        h=[max(up(row[q]) for row in rows) for q in range(W+1)]
        hbar=[up(sum((int(n)*row[q] for n,row in zip(data['histogram_multiplicities'],rows)),arb(0))/L) for q in range(W+1)]
        # Arb's general nth-root path can return nan at exact zero.
        maximum=[min(v,arb(0) if a==0 else up((moment*a**(power-1)).root(power))) for v,moment,a in zip(maximum,h,atoms)]
        mean=[min(v,arb(0) if n==0 else up((moment*(n/L)**(power-1)).root(power))) for v,moment,n in zip(mean,hbar,nu)]
        per_state=[[min(v,arb(0) if a==0 else up((moment*a**(power-1)).root(power)))
                    for v,moment,a in zip(previous,row,atoms)] for previous,row in zip(per_state,rows)]
    maximum=[min(v,max(row[q] for row in per_state)) for q,v in enumerate(maximum)]
    mean=[min(v,up(sum((int(n)*row[q] for n,row in zip(data['histogram_multiplicities'],per_state)),arb(0))/L)) for q,v in enumerate(mean)]
    weights=[comb(W,q)*aq(p)**q*aq(1-p)**(W-q) for q in range(W+1)]
    return up(sum((w*v for w,v in zip(weights,maximum)),arb(0))),up(sum((w*v for w,v in zip(weights,mean)),arb(0)))


def floating(data,probabilities,tilt):
    matrix=base.floating(data,probabilities,tilt)
    maximum,mean=float_bounds(data,float(1-probabilities[0]),np.exp(-tilt))
    alpha=2.**-data['updates'];L=(1<<data['bits'])-1
    matrix[1,0]=min(matrix[1,0],alpha*maximum+matrix[1,2]/L)
    matrix[2,0]=min(matrix[2,0],alpha*mean+matrix[2,2]/L)
    return matrix


def outward_at_z(data,probabilities,z):
    matrix=base.outward_at_z(data,probabilities,z)
    maximum,mean=bounds(data,base.base.base.activity(probabilities),z)
    alpha=arb(2)**-data['updates'];L=(1<<data['bits'])-1
    matrix[1,0]=min(matrix[1,0],up(alpha*maximum+matrix[1,2]/L))
    matrix[2,0]=min(matrix[2,0],up(alpha*mean+matrix[2,2]/L))
    return matrix


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
