"""Active-row rank screen keeping exact regional capacity and route factors.

Unlike the homogeneous slot-probability relaxation, each region admits at
most A active slots. Exact-A feasible row supports have dimension W-10*A;
the larger summed family can include infeasible supports conservatively.
All numerical results are FLOATING screens, not directed certificates.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_choose
from reduced_distance_certificate import log_polynomial_product
from sparse_multirow_zero_set import local_rank_transfer
from cauchy_joint_return import rank_return_moment


def log_region(active,y,tau,rows,joint_returns=False,width=16):
    if joint_returns and width != 16:
        raise ValueError('joint return table is currently width16 only')
    local=np.asarray([local_rank_transfer(w,y,tau,width=width) for w in range(min(width,active)+1)],dtype=float)
    if joint_returns:
        for w in range(1,min(3,active)+1):
            local[w,1,0]=min(local[w,1,0],rank_return_moment(w,y,tau,average_inputs=True))
    with np.errstate(divide='ignore'):
        factors=np.log(local)
    for w in range(len(local)):
        factors[w]+=math.log(math.comb(width,w))
    result=np.array([[[0.,-math.inf],[-math.inf,0.]]])
    if rows%width:
        raise ValueError('regional length must be a multiple of packet width')
    power=rows//width
    while power:
        if power&1:
            result=log_polynomial_product(result,factors,active)
        power>>=1
        if power:
            factors=log_polynomial_product(factors,factors,active)
    for w in range(active+1):
        result[w]+=log_choose(active,w)-log_choose(rows,w)+w*math.log(tau)
    return logsumexp(result,axis=0)


def bound(active,y,tau,k_log=20,joint_returns=False,length=32,dimension=22,width=16):
    g=geometry(k_log,length,dimension)
    region=log_region(active,y,tau,g['rows'],joint_returns,width)
    vector=np.array([0.,-math.inf])
    for _ in range(length):
        vector=logsumexp(vector[:,None]+region,axis=0)
    return (log_choose(g['rows'],active)-((length-dimension)*active+1)*math.log(tau)
            -g['threshold']*math.log(y)+float(logsumexp(vector)))/math.log(2)


def search(active,k_log=20,joint_returns=False,length=32,dimension=22,width=16):
    g=geometry(k_log,length,dimension)
    def objective(point):
        logtilt,logtau=point
        logy=-math.exp(logtilt)*active/g['n']
        return bound(active,math.exp(logy),math.exp(logtau),k_log,joint_returns,length,dimension,width) if logy>=-25 else 1e100
    starts=[[math.log(tilt),math.log(tau)] for tilt in (100,1000) for tau in (16,255)]
    choices=[minimize(objective,start,method='L-BFGS-B',bounds=[(-2,13),(0,math.log(255))],
                      options={'maxiter':130}) for start in starts]
    best=min(choices,key=lambda result:result.fun)
    return dict(active=active,k_log=k_log,outer=[length,dimension],width=width,margin_bits=-float(best.fun),
                logtilt=float(best.x[0]),logtau=float(best.x[1]),
                rigorous=False,joint_returns=joint_returns,
                method='exact-regional-count and zero-set rank potential')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active',type=int,nargs='+',default=[1,2,4,8,16,32])
    parser.add_argument('--k-log',type=int,default=20)
    parser.add_argument('--joint-returns',action='store_true')
    parser.add_argument('--length',type=int,default=32)
    parser.add_argument('--dimension',type=int,default=22)
    parser.add_argument('--width',type=int,default=16)
    args=parser.parse_args()
    for active in args.active:
        print(json.dumps(search(active,args.k_log,args.joint_returns,args.length,args.dimension,args.width)),flush=True)
