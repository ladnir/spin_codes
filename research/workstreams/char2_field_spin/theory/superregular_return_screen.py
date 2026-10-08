"""Floating rank cover for symmetric Cauchy injection/emission.

CHANGED inner: C=A^T, CA=0, both maps superregular, independent nonzero
input and state diagonals. Every return branch consumes `state` distinct
diagonal pivots after fixing an actual incoming support. Returned outputs
use a coarse deterministic weight lower bound, not selected-zero branches.
No all-subset eta invertibility assumption is imported. No native timing
or directed whole-code certificate is provided by these floating screens.
"""

import argparse
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_choose
from diagonal_refresh_rank import rank_potential_bound
from reduced_distance_certificate import log_polynomial_product


@lru_cache(maxsize=64,typed=True)
def gamma(tau, state, q=256):
    return tuple(max(rank_potential_bound(r,z,tau,q,state)
                     for r in range(1,state+1)) for z in range(state))


def return_branches(weight, state):
    return sum(math.comb(state,h) for h in range(max(1,state-weight+1),state+1))


def return_flat_factor(weight, tau, state):
    """Charge the deterministic rank of each feasible state zero flat.

    Each rank-k flat has a distinct canonical k-row basis. A return with
    w live inputs requires at most w-1 zero state coordinates. Nonzero
    state excludes the rank-s flat. No current diagonal is used here.
    """
    return sum(math.comb(state,k)/tau**k for k in range(min(weight-1,state-1)+1))


def local_matrix(weight, y, tau, state=8, width=16, flat_cover=True):
    if not (0<=weight<=width and 0<state<width and 0<y<=1 and 1<=tau<=255):
        raise ValueError('invalid matrix parameters')
    zero = y*0
    out = [[zero,zero],[zero,zero]]
    if weight == 0:
        out[0][0] = y**0
    else:
        out[0][1] = y**weight
        if weight > state:
            out[0][0] = y**weight/tau**state
        minimum = max(1,width-state+1-weight)
        if state == 8 and width in (16,32) and weight == 1:
            minimum = width-7  # Exact symmetric-map singleton law, eta=2.
        branches = return_flat_factor(weight,tau,state) if flat_cover else return_branches(weight,state)
        out[1][0] = branches*y**minimum/tau**state
    roots = gamma(tau,state)
    out[1][1] = sum(math.comb(weight,i)*math.comb(width-weight,z)
                    *y**(width-i-z)*(1-y)**(i+z)*roots[z]/tau**i
                    for i in range(weight+1)
                    for z in range(min(state-1,width-weight)+1))
    return out


def region(active,y,tau,rows,state=8,flat_cover=True,width=16):
    if rows%width:
        raise ValueError('regional length must be packet aligned')
    local=np.asarray([local_matrix(w,y,tau,state,width,flat_cover=flat_cover)
                      for w in range(min(width,active)+1)],dtype=float)
    with np.errstate(divide='ignore'):
        coefficients=np.log(local)
    for w in range(len(local)):
        coefficients[w]+=math.log(math.comb(width,w))
    result=np.array([[[0.,-math.inf],[-math.inf,0.]]])
    power=rows//width
    while power:
        if power&1:
            result=log_polynomial_product(result,coefficients,active)
        power>>=1
        if power:
            coefficients=log_polynomial_product(coefficients,coefficients,active)
    for w in range(active+1):
        result[w]+=log_choose(active,w)-log_choose(rows,w)+w*math.log(tau)
    return logsumexp(result,axis=0)


def bound(active,point,k_log=20,length=32,dimension=22,state=8,flat_cover=True,width=16):
    g=geometry(k_log,length,dimension)
    logtilt,logtau=point
    logy=-math.exp(logtilt)*active/g['n']
    table=region(active,math.exp(logy),math.exp(logtau),g['rows'],state,flat_cover,width)
    vector=np.array([0.,-math.inf])
    for _ in range(length):
        vector=logsumexp(vector[:,None]+table,axis=0)
    return (log_choose(g['rows'],active)-((length-dimension)*active+1)*logtau
            -g['threshold']*logy+float(logsumexp(vector)))/math.log(2)


def search(active,k_log=20,length=32,dimension=22,state=8,iterations=50,flat_cover=True,width=16):
    options=[]
    for point in ((5.,math.log(255)),(5.,math.log(16))):
        result=minimize(lambda x:bound(active,x,k_log,length,dimension,state,flat_cover,width),point,
                        method='L-BFGS-B',bounds=[(-2,10),(0,math.log(255))],
                        options={'maxiter':iterations,'ftol':1e-10})
        options.append(result)
    best=min(options,key=lambda result:result.fun)
    return dict(active=active,k_log=k_log,outer=[length,dimension],state=state,width=width,
                margin_bits=-float(best.fun),logtilt=float(best.x[0]),logtau=float(best.x[1]),
                success=bool(best.success),iterations=int(best.nit),rigorous=False,flat_cover=flat_cover,
                construction='symmetric-Cauchy fresh-diagonal state')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active',type=int,nargs='+',default=[4,8,16,32,64])
    parser.add_argument('--outer',type=int,nargs=2,default=(32,22))
    parser.add_argument('--state',type=int,default=8)
    parser.add_argument('--width',type=int,choices=(16,32),default=16)
    parser.add_argument('--k-log',type=int,default=20)
    parser.add_argument('--iterations',type=int,default=50)
    parser.add_argument('--raw-support-sum',action='store_true')
    args=parser.parse_args()
    for active in args.active:
        print(json.dumps(search(active,args.k_log,*args.outer,args.state,args.iterations,
                                not args.raw_support_sum,args.width)),flush=True)
