"""Full-state screen using input-live return moments without zero-flat sums.

The construction and all nonreturn bounds are unchanged. A pure return
is now bounded jointly before exposing either current diagonal. Only the
active input coordinates must remain nonzero. No result is a certificate.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import json
import math
import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
import gf8_fullstate_screen as base
from gf8_fullstate_sparse_tilt import FROZEN
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from fractions import Fraction as F
from partial_live_rank import mds_return_bound, partial_live_moment
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


@lru_cache(maxsize=4096)
def return_factor(weight, state, cap, tau):
    return mds_return_bound(weight, state, cap+weight, tau)


@lru_cache(maxsize=8192)
def joint_factor(weight,state,cap,charge,remaining,tau):
    # For a fixed input-live vector, state support h must satisfy
    # h+w>s and h>z for a nonzero probability. The old joint rank proof
    # therefore gives 'charge' BEFORE any current diagonal is exposed.
    # Extended-map invertibility embeds the joint kernel into the
    # remaining output coordinates, hence nullity <= min(w,remaining).
    probability = F(1,255**charge)
    return max(partial_live_moment(d,weight,min(weight,remaining),probability,tau)
               for d in range(1,cap+weight+1))


def return_bound(weight,y,tau,state=8,cap=None):
    if state < 8:
        raise ValueError('joint partial-live screen requires state >= packet width')
    cap = state if cap is None else cap
    pure = y**base.RETURN_MINIMA[state][weight-1]*return_factor(weight,state,cap,tau)
    joint = sum(count*y**(8-i-z)*(1-y)**(i+z)
        *joint_factor(weight,state,cap,charge,8-i-z,tau)
        for i,z,count,charge in base.joint_terms(weight,state))
    return min(base.return_bound(weight,y,tau,state),pure,joint)


def local_matrix(weight, y, tau, state=8):
    if state < 8:
        raise ValueError('joint partial-live screen requires state >= packet width')
    result = [list(row) for row in base.local_matrix(weight,y,tau,state)]
    if weight:
        result[1][0] = return_bound(weight,y,tau,state)
    return result


def packet(logx,logy,logtau,state=8):
    result = base.packet(logx,logy,logtau,state)
    x,y,tau = math.exp(logx),math.exp(logy),min(255.,math.exp(logtau))
    result[1,0] = sum(math.comb(8,w)*x**w*return_bound(w,y,tau,state) for w in range(1,9))
    return result


def sparse_bound(active,point,state=8,dimension=96,k_log=20,length=128):
    g = geometry(k_log,length,dimension)
    logtilt,logtau = point[:2]
    loglambda = point[2] if len(point)>2 else 0.
    if loglambda < 0:
        raise ValueError('lambda must be at least one')
    logy = -math.exp(logtilt)*active/g['n']
    local = np.asarray([local_matrix(w,math.exp(logy),math.exp(logtau),state)
                        for w in range(min(8,active)+1)],dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    for w in range(len(factor)):
        factor[w] += math.log(math.comb(8,w))
    result = np.array([[[0.,-math.inf],[-math.inf,0.]]])
    power = g['rows']//8
    while power:
        if power&1:
            result = log_polynomial_product(result,factor,active)
        power >>= 1
        if power:
            factor = log_polynomial_product(factor,factor,active)
    region = regional_from_coefficients(result,active,g['rows'],logtau+loglambda)
    return (log_choose(g['rows'],active)-((length-dimension)*active+1)*logtau
            -(length-dimension+1)*active*loglambda-g['threshold']*logy
            +log_terminal(region,length))/math.log(2)


def middle_bound(point,weight,state=8,dimension=96,k_log=20,length=128):
    g = geometry(k_log,length,dimension)
    logx,logy,logtau,logtilt = point
    return (length*modal_log_penalty(g['rows'])
            +g['rows']*source_row_log(weight/g['n'],logtau,logtilt,length,dimension,True)
            +log_terminal_power(packet(logx,logy,logtau,state),g['n']//8)
            -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def search(value,sparse=True,state=8,dimension=96,k_log=20,length=128,iterations=50,start=None):
    g = geometry(k_log,length,dimension)
    target = int(value) if sparse else round(value*g['n'])
    if sparse:
        default = FROZEN.get((state,target),(6.25,1.8))
        starts = [tuple(start) if start is not None else default]
        objective = lambda p:sparse_bound(target,p,state,dimension,k_log,length)
        bounds = [(-2,10),(0,math.log(255))]
    else:
        starts = [tuple(start) if start is not None else (-1.,-.4,1.8,1.)]
        objective = lambda p:middle_bound(p,target,state,dimension,k_log,length)
        bounds = [(-30,8),(-20,0),(0,math.log(255)),(-20,20)]
    choices = [minimize(objective,p,method='L-BFGS-B',bounds=bounds,
                options={'maxiter':iterations,'ftol':1e-11,'maxls':25}) for p in starts]
    best = min(choices,key=lambda v:v.fun)
    return dict(outer=[length,dimension],width=8,state=state,k_log=k_log,
                target=target,mode='sparse' if sparse else 'middle',n=g['n'],
                margin_bits=-float(best.fun),log_parameters=best.x.tolist(),
                optimizer_success=bool(best.success),iterations=int(best.nit),
                rigorous=False,whole_code=False,return_bound='partial-live-pure-and-joint')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state',type=int,nargs='+',choices=(8,12,16),default=[8,16])
    parser.add_argument('--active',type=int,nargs='*',default=[64,128])
    parser.add_argument('--density',type=float,nargs='*',default=[])
    parser.add_argument('--iterations',type=int,default=50)
    args = parser.parse_args()
    for state in args.state:
        for active in args.active:
            print(json.dumps(search(active,True,state,iterations=args.iterations)),flush=True)
        last = None
        for density in args.density:
            receipt = search(density,False,state,iterations=args.iterations,start=last)
            last = receipt['log_parameters']
            print(json.dumps(receipt),flush=True)
