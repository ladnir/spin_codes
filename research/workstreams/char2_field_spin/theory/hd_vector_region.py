"""Regional proof diagnostics for an additive vector-state observation inner.

Exact fixed-message transfers and a separately justified rank envelope are
provided. Floating screens cover one active-row class, never a whole code.
No native timing or benchmark functionality is included.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from itertools import combinations
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from cauchy_firstmoment_screen import log_choose
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


# Bounded K20 / GF65536 / state8 / distance2/5 search witnesses. The rank
# witnesses were selected before the sharp active-input observation factor;
# they remain valid inputs to the stronger current envelope.
SCREEN_WITNESSES=(
    (128,64,1,True,(7.2589596999143495,11.090354888959125)),
    (128,64,16,True,(5.749238846616855,2.411027938169151)),
    (128,64,64,True,(5.452787400473079,1.7030810705046269)),
    (128,64,1,False,(7.258959695956655,)),
    (128,64,16,False,(6.528115346092129,)),
    (128,64,64,False,(5.8351656549905835,)),
    (128,64,16,True,(5.74923747235794,2.411023647523447,0.)),
    (128,64,64,True,(5.452787425704088,1.703082235510044,0.)),
    (256,128,64,True,(6.148648223918853,1.7073023580674087,0.)),
)


def local(y, q, state, tau=None):
    """Rows/columns denote whole-vector state zero/live, not coordinate flags.

    tau=None: exact fixed-message expectation of y**output_weight.
    tau given: rank envelope retaining all earlier live scalar inputs and
    live state vectors in its conceptual history; not a transition law.
    """
    if not 0 < y <= 1 or q < 2 or state < 1:
        raise ValueError('invalid local parameters')
    one = y*0+1
    zero = y*0
    if tau is None:
        a = F(1,q**state) if isinstance(y,F) else q**(-state)
        g = (one+(q-1)*y)/q
        return (((one,zero),(zero,g)),
                ((a*y,(one-a)*y),(a*g,(one-a)*g)))
    if not 1 <= tau <= q:
        raise ValueError('rank mark must lie between one and the field size')
    inverse = one/tau
    # A rank-one state image dies after a nonzero zero-output equation.
    # Image rank at least two has rank-zero probability at most q^-2.
    g = inverse+(one-inverse)/(q*q)
    h0 = y+(one-y)*g
    # For live X, every nonzero observation equation proportional to X
    # kills the required input liveness. This recovers a sharp 1/tau charge.
    h1 = y+(one-y)*inverse
    # For the new-zero entry, expose the return equations first. Surviving
    # histories have S=-dX with d nonzero. The output-zero constraint then
    # preserves X only for f.d=c, an event of probability exactly 1/q.
    returned = y+(one-y)/q
    a = inverse**state
    birth_zero = F(1,q**state) if isinstance(y,F) else q**(-state)
    # With old state identically zero, a live scalar input gives new state
    # d*X. This branch has no source-rank restriction: d is zero or live.
    return (((one,zero),(zero,h0)),
            ((birth_zero*y,(one-birth_zero)*y),(a*returned,h1)))


def mm(a,b):
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(2)),F(0))
                       for j in range(2)) for i in range(2))


def mv(v,a):
    return tuple(sum((v[i]*a[i][j] for i in range(2)),F(0)) for j in range(2))


IDENTITY=((F(1),F(0)),(F(0),F(1)))
ZERO=((F(0),F(0)),(F(0),F(0)))


def exact_coefficients(rows,y,q,state,tau=None):
    """Small exact reference, retaining chronological matrix multiplication."""
    if not 0 <= rows <= 24:
        raise ValueError('small reference accepts at most24 rows')
    t0,t1=local(F(y),q,state,None if tau is None else F(tau))
    answer=[IDENTITY]
    for _ in range(rows):
        following=[ZERO]*(len(answer)+1)
        for h,old in enumerate(answer):
            for active,t in enumerate((t0,t1)):
                add=mm(old,t)
                following[h+active]=tuple(tuple(following[h+active][i][j]+add[i][j]
                                                for j in range(2)) for i in range(2))
        answer=following
    return answer


def exact_region(rows,occupancy,y,q,state,tau=None):
    numerator=exact_coefficients(rows,y,q,state,tau)[occupancy]
    return tuple(tuple(value/math.comb(rows,occupancy) for value in row) for row in numerator)


def exact_profile(profile,rows,y,q,state,tau=None):
    coeff=exact_coefficients(rows,y,q,state,tau)
    v=(F(1),F(0))
    for h in profile:
        v=mv(v,tuple(tuple(value/math.comb(rows,h) for value in row) for row in coeff[h]))
    return sum(v)


def exact_active_envelope(active,length,dimension,rows,y,q,state,tau=None,lam=F(1)):
    """Small exact upper bound, with optional total input-weight mark."""
    if not 1<=active<=rows or not 1<=dimension<length or lam<1:
        raise ValueError('invalid source geometry or minimum-weight mark')
    mark=q if tau is None else F(tau)
    coeff=exact_coefficients(rows,y,q,state,tau)
    matrix=tuple(tuple(sum((coeff[h][i][j]*F(math.comb(active,h),math.comb(rows,h))
                               *(mark*lam)**h for h in range(active+1)),F(0))
                       for j in range(2)) for i in range(2))
    v=(F(1),F(0))
    for _ in range(length):
        v=mv(v,matrix)
    divisor=q-1 if tau is None else mark
    return (math.comb(rows,active)*sum(v)/mark**((length-dimension)*active)
            /divisor/lam**((length-dimension+1)*active))


def log_coefficients(rows,active,y,q,state,tau):
    with np.errstate(divide='ignore'):
        factor=np.log(np.asarray(local(y,q,state,tau),dtype=float))
    result=np.asarray([[[0.,-math.inf],[-math.inf,0.]]])
    while rows:
        if rows&1:
            result=log_polynomial_product(result,factor,active)
        rows>>=1
        if rows:
            factor=log_polynomial_product(factor,factor,active)
    return result


def screen_bound(active,point,length=128,dimension=64,k_log=20,q=65536,state=8,
                 distance=F(2,5),rank=True,components=False):
    rows=(2**k_log+dimension-1)//dimension
    n=length*rows
    cutoff=(n*distance.numerator+distance.denominator-1)//distance.denominator-1
    logy=-math.exp(float(point[0]))*active/n
    if logy<-700:
        return math.inf
    logmark=float(point[1]) if rank else math.log(q)
    loglam=float(point[2]) if len(point)>2 else 0.
    mark=math.exp(logmark)
    coefficients=log_coefficients(rows,active,math.exp(logy),q,state,mark if rank else None)
    region=regional_from_coefficients(coefficients,active,rows,logmark+loglam)
    terms=dict(source=log_choose(rows,active)-(length-dimension)*active*logmark
                      -(logmark if rank else math.log(q-1)),
               output=-cutoff*logy, minimum_weight=-(length-dimension+1)*active*loglam,
               regional=log_terminal(region,length))
    value=sum(terms.values())/math.log(2)
    return (value,{key:v/math.log(2) for key,v in terms.items()}) if components else value


def search(active,length=128,dimension=64,k_log=20,q=65536,state=8,
           distance=F(2,5),rank=True,iterations=40,start=None,minimum_weight=False):
    args=(length,dimension,k_log,q,state,distance,rank)
    if rank:
        starts=[(math.log(1024),math.log(q)),(math.log(256),math.log(16))]
        bounds=[(-2,13),(0,math.log(q))]
    else:
        starts=[(math.log(1024),),(math.log(256),)]
        bounds=[(-2,13)]
    if minimum_weight:
        starts=[tuple(p)+(0.,) if rank else (p[0],math.log(q),0.) for p in starts]
        bounds=[(-2,13),(0,math.log(q)) if rank else (math.log(q),math.log(q)),(0,math.log(q))]
    if start is not None:
        starts.insert(0,tuple(start))
    results=[minimize(lambda p:screen_bound(active,p,*args),point,method='L-BFGS-B',bounds=bounds,
                      options={'maxiter':iterations,'ftol':1e-10,'maxls':12}) for point in starts]
    best=min(results,key=lambda r:r.fun)
    value,parts=screen_bound(active,best.x,*args,components=True)
    return dict(outer=[length,dimension],q=q,state=state,k_log=k_log,active=active,
                distance=str(distance),mode='rank-envelope' if rank else 'fixed-message-envelope',
                birth_row='exact-q^-state',
                minimum_weight_mark=minimum_weight,margin_bits=-value,log_parameters=best.x.tolist(),
                components_log2=parts,optimizer_success=bool(best.success),iterations=int(best.nit),
                rigorous=False,whole_code=False,scope='one active-row class; floating diagnostic only')


def replay_witnesses():
    for length,dimension,active,is_rank,point in SCREEN_WITNESSES:
        value,parts=screen_bound(active,point,length,dimension,components=True,rank=is_rank)
        yield dict(outer=[length,dimension],active=active,q=65536,state=8,k_log=20,
                   distance='2/5',mode='rank-envelope' if is_rank else 'fixed-message-envelope',
                   log_parameters=list(point),margin_bits=-value,components_log2=parts,
                   rigorous=False,whole_code=False,optimizer_run=False,
                   scope='fixed floating witness; one active-row class only')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outer',type=int,nargs=2,default=(128,64))
    p.add_argument('--q',type=int,default=65536)
    p.add_argument('--state',type=int,default=8)
    p.add_argument('--active',type=int,nargs='+',default=(1,16,64))
    p.add_argument('--exact-message',action='store_true')
    p.add_argument('--minimum-weight',action='store_true')
    p.add_argument('--iterations',type=int,default=40)
    p.add_argument('--replay',action='store_true')
    args=p.parse_args()
    if args.replay:
        for result in replay_witnesses():print(json.dumps(result),flush=True)
    else:
        for active in args.active:
            print(json.dumps(search(active,*args.outer,q=args.q,state=args.state,
                                    rank=not args.exact_message,iterations=args.iterations,
                                    minimum_weight=args.minimum_weight)),flush=True)
