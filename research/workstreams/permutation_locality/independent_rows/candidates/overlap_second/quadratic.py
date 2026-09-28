"""Positive quadratic upper witnesses for the exact second overlap moment.

Binary64 only proposes a nonnegative rational slope. Every discrete
overlap constraint is then rebuilt outward; the selected slope is never
used as evidence that an inequality holds.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement

import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb

from second import moments  # Also establishes the shared research imports.
from mass_density_screen import baseline
from occupancy_memory import C,Z


def upper_intercepts(levels,weight,tilt,c):
    """Outward a_v satisfying exp(-tilt*(v-d)) <= a_v+c*d^2."""
    if (type(weight) is not int or weight<1 or not levels
            or any(type(v) is not int or v<weight for v in levels)
            or not isinstance(c,Q) or c<0 or Q(tilt)<=0):
        raise ValueError('nonnegative exact slope and nonnegative output weights required')
    slope=arb(c.numerator)/c.denominator
    lam=arb(Q(tilt).numerator)/Q(tilt).denominator
    return {v:max(arb(0),*(baseline.up((-lam*(v-k)).exp()-slope*k*k)
                           for k in range(weight%2,weight+1,2))) for v in levels}


def witness(classes,weight,second,tilt):
    if (type(weight) is not int or weight<1 or not isinstance(second,Q) or not 0<=second<=weight*weight
            or not classes or any(type(v) is not int or v<weight or not isinstance(p,Q) or p<0
                                  for v,p in classes.items()) or sum(classes.values(),Q(0))>1 or Q(tilt)<=0):
        raise ValueError('exact class masses, nonnegative moment bound and output exponents required')
    levels=sorted(classes)
    d=np.arange(weight%2,weight+1,2,dtype=float)
    squares=d*d
    f=np.exp(-float(Q(tilt))*(np.array(levels)[:,None]-d[None,:]))
    positive=squares>0
    limit=float(np.max(f[:,positive]/squares[positive]))
    probabilities=np.array([float(classes[v]) for v in levels])
    def objective(x):
        c=limit*x
        a=np.maximum(0,np.max(f-c*squares[None,:],axis=1))
        return float(probabilities@a+c*float(second))
    fit=minimize_scalar(objective,bounds=(0.,1.),method='bounded',options={'xatol':1e-12})
    candidates=[0.,1.]
    if np.isfinite(fit.x):candidates.append(float(np.clip(fit.x,0,1)))
    proposal=limit*min(candidates,key=objective)
    c=Q(round(proposal*(1<<80)),1<<80)
    slope=arb(c.numerator)/c.denominator
    intercepts=upper_intercepts(levels,weight,tilt,c)
    value=baseline.up(sum((arb(p.numerator)/p.denominator*intercepts[v] for v,p in classes.items()),arb(0))
                      +slope*second.numerator/second.denominator)
    return value,c,intercepts


def zero_refine(base,second,feedback,tilt,penalty,*,minimum=5,maximum=10,rounds=2):
    if (not 1<=minimum<=maximum<len(base) or not 0<Q(penalty)<=1 or Q(tilt)<=0
            or type(rounds) is not int or rounds<1):
        raise ValueError('valid complete shape range, witnesses and update count required')
    required={s for j in range(minimum,maximum+1) for s in combinations_with_replacement(range(1,5),j)}
    if not required<=set(second) or not required<=set(feedback):
        raise ValueError('complete second-moment and feedback censuses required')
    upper={j:arb(0) for j in range(minimum,maximum+1)}
    for shape in sorted(required,key=lambda s:(len(s),s)):
        zero,_,D,counts=feedback[shape]
        if D!=second[shape][0] or zero+sum(counts.values())!=D:
            raise ValueError('second moment and feedback counting measures disagree')
        W=sum(shape)
        # Source zero contributes exactly W^2 to the centered square.
        restricted=second[shape][1]-Q(W*W*zero,D)
        value,_,_=witness({v:Q(n,D) for v,n in counts.items()},W,restricted,tilt)
        scale=Q(penalty)**shape.count(4)/2**rounds
        upper[len(shape)]=max(upper[len(shape)],baseline.up(value*scale.numerator/scale.denominator))
    result=[t*1 for t in base]
    changes=0
    for j,value in upper.items():
        changes+=int(value<result[j][C,Z])
        result[j][C,Z]=min(result[j][C,Z],value)
    print('QUADRATIC OVERLAP tightened zero-return entries',changes,flush=True)
    return result
