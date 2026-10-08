"""FLOAT fixed-input-weight screen for the symmetric state-eight candidate.

Uses the existing regional rank/source comparison with a changed packet
cover. Finite route conditioning, zero-start/any-finish, and output threshold
are retained. Isolated optimized weights are NOT interval certificates.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_terminal_power
from state_rank_cap_transfer import return_minimum_weight
from state8_joint_return import joint_return_bound
from superregular_return_screen import local_matrix


def packet(logx, logy, logtau, width=16, joint=False):
    y, tau = math.exp(logy), min(255., math.exp(logtau))
    values = np.asarray([local_matrix(w,y,tau,width=width) for w in range(width+1)],dtype=float)
    if width == 16:
        for w in range(1,17):
            old_minimum = 9 if w == 1 else max(1,9-w)
            values[w,1,0] *= y**(return_minimum_weight(w)-old_minimum)
    if joint:
        for w in range(1,width+1):
            values[w,1,0] = joint_return_bound(w,y,tau,width)
    with np.errstate(divide='ignore'):
        logs = np.log(values)
    for w in range(width+1):
        logs[w] += math.log(math.comb(width,w))+w*logx
    return logsumexp(logs,axis=0)


def bound(point, weight, k_log=20, length=256, dimension=176, width=16,joint=False):
    g = geometry(k_log,length,dimension)
    if not 0 < weight < g['n']:
        raise ValueError('weight must lie strictly inside encoded length')
    alpha = weight/g['n']
    logx, logy, logtau, logtilt = map(float,point)
    row = float(logsumexp([length*math.log1p(-alpha),
                          -(length-dimension)*logtau+length*np.logaddexp(
                              math.log1p(-alpha),logtau+logtilt+math.log(alpha))]))
    logarithms = packet(logx,logy,logtau,width,joint)
    shift = float(logarithms.max())
    power = g['n']//width
    inner = power*shift+log_terminal_power(np.exp(logarithms-shift),power)
    return (length*math.log(g['rows']+1)+g['rows']*row+inner
            -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def search(alpha, length=256, dimension=176, k_log=20, iterations=80, start=None,width=16,joint=False):
    g = geometry(k_log,length,dimension)
    weight = max(1,round(alpha*g['n']))
    starts = [start] if start is not None else [(-2,-.2,math.log(16),0),
                                               (0,-1,math.log(255),-1)]
    candidates = [minimize(bound,point,args=(weight,k_log,length,dimension,width,joint),method='L-BFGS-B',
                           bounds=[(-30,12),(-25,0),(0,math.log(255)),(-20,20)],
                           options={'maxiter':iterations,'ftol':1e-11,'maxls':20})
                  for point in starts]
    best = min(candidates,key=lambda item:item.fun)
    return dict(outer=[length,dimension],k_log=k_log,width=width,joint_return=joint,input_weight=weight,
                input_density=weight/g['n'],margin_bits=-float(best.fun),
                log_parameters=best.x.tolist(),iterations=int(best.nit),
                optimizer_success=bool(best.success),rigorous=False,whole_code=False,
                finite_factors_included=True,scope='one exact input-weight class')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer',type=int,nargs=2,default=[256,176])
    parser.add_argument('--alpha',type=float,nargs='+',default=[.01,.02,.04,.06,.08,.1])
    parser.add_argument('--k-log',type=int,default=20)
    parser.add_argument('--iterations',type=int,default=80)
    parser.add_argument('--width',type=int,choices=(16,32),default=16)
    parser.add_argument('--joint',action='store_true')
    parser.add_argument('--start',type=float,nargs=4)
    args = parser.parse_args()
    previous = args.start
    for alpha in args.alpha:
        value = search(alpha,*args.outer,args.k_log,args.iterations,previous,args.width,args.joint)
        previous = value['log_parameters']
        print(json.dumps(value),flush=True)
