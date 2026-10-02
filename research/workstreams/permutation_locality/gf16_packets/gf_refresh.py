"""IID GF16 refresh envelope with occupancy-conditioned lazy returns.

For j active packets, feedback cancellation and emitted weight are bounded
together by the expansion-distance floor. The setup and state invariant
are unchanged. Input values must be independently uniform nonzero when active.
"""
from math import comb
from fractions import Fraction as Q
import numpy as np
from flint import arb
import refresh_kernel as base
import occupancy_kernel


def prepare(images,columns,bits,updates=2):
    data=occupancy_kernel.prepare(base.prepare(images,columns,bits,updates),True)
    data['distance']=min(sum(w*int(n) for w,n in enumerate(row)) for row in data['histograms'])
    return data


def actual(updates=2):
    images,columns,_=base.base.maps();data=prepare(images,columns,19,updates)
    print('GF16 REFRESH KERNEL exact feedback through',len(data['exact_feedback'])-1,'packets',flush=True)
    return data


def activity(probabilities):
    p=list(map(Q,probabilities))
    if (len(p)!=5 or min(p)<0 or sum(p)!=1
            or any(p[j]!=(1-p[0])*Q(comb(4,j),15) for j in range(1,5))):
        raise ValueError('uniform nonzero GF16 values conditional on activity required')
    return 1-p[0]


def floating(data,probabilities,tilt):
    matrix=base.floating(data,probabilities,tilt,feedback_aware=True)
    p=float(1-probabilities[0]);z=np.exp(-tilt);W=data['windows'];L=(1<<data['bits'])-1
    alpha=2.**-data['updates']
    # Only proposal arithmetic uses binary64; outward() checks the law exactly.
    weights=np.array([comb(W,j)*p**j*(1-p)**(W-j)*z**max(0,data['distance']-4*j) for j in range(W+1)])
    cancel=float(weights@np.array(list(map(float,data['nonzero_atom_caps']))))
    uniform=float(weights@np.array([float(1-x) for x in data['zero_probabilities']]))/L
    # Column 2 is the existing uniform-refresh mass beta H (or beta Hbar).
    matrix[1,0]=min(matrix[1,0],alpha*cancel+matrix[1,2]/L)
    matrix[2,0]=min(matrix[2,0],alpha*uniform+matrix[2,2]/L)
    return matrix


def outward_at_z(data,probabilities,z):
    p=activity(probabilities);matrix=base.outward_at_z(data,probabilities,z,feedback_aware=True)
    W=data['windows'];L=(1<<data['bits'])-1;alpha=arb(2)**-data['updates']
    weights=[comb(W,j)*base.aq(p)**j*base.aq(1-p)**(W-j)*z**max(0,data['distance']-4*j) for j in range(W+1)]
    cancel=sum((w*base.aq(b) for w,b in zip(weights,data['nonzero_atom_caps'])),arb(0))
    uniform=sum((w*base.aq(1-v) for w,v in zip(weights,data['zero_probabilities'])),arb(0))/L
    matrix[1,0]=min(matrix[1,0],base.up(alpha*cancel+matrix[1,2]/L))
    matrix[2,0]=min(matrix[2,0],base.up(alpha*uniform+matrix[2,2]/L))
    return matrix


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-base.aq(tilt)).exp())
