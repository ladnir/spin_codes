"""Direct fixed-word sparse first moment for width8/state8,16 controls.

Fresh input/state nonzero diagonals are exposed for one fixed nonzero
message. Selected zero events cost powers of 1/255, not rank-potential
marks. The outer's exact-support multiplicity is bounded by 256**(W-pA).
Projective normalization divides the resulting first moment by 255 once.
This is a floating diagnostic, not a certificate or a failure lower bound.
"""
from __future__ import annotations
import argparse
import json
import math
import numpy as np
from scipy.optimize import minimize_scalar

from cauchy_firstmoment_screen import geometry, log_choose
from gf8_fullstate_screen import RETURN_MINIMA, joint_terms
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


def local_matrix(w,y,state=8,q=256):
    if state not in (8,16) or not 0<=w<=8 or not 0<y<=1:
        raise ValueError('fixed control requires width8 and state8 or16')
    one=y**0
    zero=y*0
    result=[[one if w==0 else zero,y**w if w else zero],[zero,zero]]
    result[1][1]=sum(math.comb(w,i)*math.comb(8-w,z)*y**(8-i-z)*(one-y)**(i+z)/(q-1)**(i+z)
                    for i in range(w+1) for z in range(min(state-1,8-w)+1))
    if w:
        pure=y**RETURN_MINIMA[state][w-1]/(q-1)**state
        joint=sum(count*y**(8-i-z)*(one-y)**(i+z)/(q-1)**charge
                  for i,z,count,charge in joint_terms(w,state))
        result[1][0]=min(pure,joint)
    return result


def sparse_bound(active,logtilt,state=8,dimension=96,k_log=20,length=128):
    g=geometry(k_log,length,dimension)
    logy=-math.exp(logtilt)*active/g['n']
    local=np.asarray([local_matrix(w,math.exp(logy),state) for w in range(min(8,active)+1)],dtype=float)
    with np.errstate(divide='ignore'):
        factor=np.log(local)
    for w in range(len(factor)):
        factor[w]+=math.log(math.comb(8,w))
    result=np.array([[[0.,-math.inf],[-math.inf,0.]]])
    power=g['rows']//8
    while power:
        if power&1:
            result=log_polynomial_product(result,factor,active)
        power >>= 1
        if power:
            factor=log_polynomial_product(factor,factor,active)
    region=regional_from_coefficients(result,active,g['rows'],math.log(256))
    return (log_choose(g['rows'],active)-(length-dimension)*active*math.log(256)-math.log(255)
            -g['threshold']*logy+log_terminal(region,length))/math.log(2)


def search(active,state=8,dimension=96,k_log=20,length=128):
    objective=lambda x:sparse_bound(active,float(x),state,dimension,k_log,length)
    choices=[minimize_scalar(objective,bounds=bounds,method='bounded',options={'xatol':1e-7})
             for bounds in ((-2.,8.),(8.,12.))]
    best=min(choices,key=lambda v:v.fun)
    return dict(outer=[length,dimension],width=8,state=state,active=active,k_log=k_log,
                margin_bits=-float(best.fun),logtilt=float(best.x),
                rigorous=False,whole_code=False,method='fixed-word-first-moment')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state',type=int,nargs='+',default=[8,16])
    parser.add_argument('--active',type=int,nargs='+',default=[64,128])
    args=parser.parse_args()
    for state in args.state:
        for active in args.active:
            print(json.dumps(search(active,state)),flush=True)
