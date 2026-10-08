"""Floating active-row rank bound with a direct regional route inequality.

For a fixed set of A outer rows and h_j selected slots in region j, the
probability that all selected slots originate in those rows is exactly
prod_j binom(A,h_j)/binom(L,h_j) <= (A/L)**W. Exact-A feasible supports
have shortened dimension W-10*A. Thus the unconditioned transfer at
x=tau*A/L bounds the full active-A class; no Bernoulli conditioning cost.
This includes extra infeasible/empty supports conservatively. Screens only.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_terminal_power, log_choose
from regional_refresh_screen import log_packet
from sparse_multirow_zero_set import transfer


def value(active, log_y, log_tau, k_log=20, zero_set=False, log_x=None,
          weight_min=None, weight_max=None):
    g = geometry(k_log)
    rows = g['rows']
    if not 1 <= active <= rows:
        raise ValueError('active count must fit row count')
    natural_x = log_tau+math.log(active/rows)
    if log_x is None:
        log_x = natural_x
    weight_charge=0.
    if weight_min is not None and weight_max is not None:
        if weight_max < weight_min:
            raise ValueError('empty support-weight interval')
        weight_charge=(natural_x-log_x)*(weight_max if natural_x>=log_x else weight_min)
    if zero_set:
        matrix = np.asarray(transfer(math.exp(log_x),math.exp(log_y),math.exp(log_tau)),dtype=float)
        with np.errstate(divide='ignore'):
            table = np.log(matrix)
    else:
        table = log_packet(log_x, log_y, log_tau, True, refresh_law='diagonal')
    shift = float(table.max())
    moment = (log_terminal_power(np.exp(table-shift), g['n']//16)
              +shift*(g['n']//16))
    return (log_choose(rows, active)-(10*active+1)*log_tau
            -g['threshold']*log_y+moment+weight_charge)/math.log(2)


def search(active, k_log=20, zero_set=False, cap_weight=False):
    g = geometry(k_log)
    if cap_weight:
        weight_min=11*active
        weight_max=min(32*active,(11*g['n']+99)//100-1)
        def objective(point):
            scale, log_decay, log_tau=point
            log_y=-math.exp(log_decay)*active/g['n']
            if log_y < -25:
                return 1e100
            return value(active,log_y,log_tau,
                         k_log,zero_set,math.log(active/g['rows'])+scale,
                         weight_min,weight_max)
        starts=[[sx,math.log(decay),math.log(tau)]
                for sx in (0.,1.) for decay in (100,1000) for tau in (16,255)]
        candidates=[minimize(objective,start,method='L-BFGS-B',
                            bounds=[(-15,15),(-4,15),(0,math.log(255))],
                            options={'maxiter':400}) for start in starts]
        best=min(candidates,key=lambda result:result.fun)
        return dict(k_log=k_log,active=active,weight_range=[weight_min,weight_max],
                    margin_bits=-float(best.fun),parameters=best.x.tolist(),
                    zero_set=zero_set,rigorous=False,cap_weight=True)
    starts = [[-scale*active/g['n'], math.log(tau)]
              for scale in (10, 100, 1000) for tau in (16,255)]
    candidates = [minimize(lambda p: value(active, p[0], p[1], k_log,zero_set), start,
                           method='L-BFGS-B', bounds=[(-25,0),(0,math.log(255))],
                           options={'maxiter':300}) for start in starts]
    best = min(candidates, key=lambda result: result.fun)
    return dict(k_log=k_log, rows=g['rows'], n=g['n'], active=active,
                log2_upper=float(best.fun), margin_bits=-float(best.fun),
                log_y=float(best.x[0]), log_tau=float(best.x[1]), rigorous=False,
                zero_set=zero_set,
                method='rank potential and direct selected-row route probability')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active', type=int, nargs='+',default=[1,2,4,16,64,256,1024,4096])
    parser.add_argument('--k-log', type=int,default=20)
    parser.add_argument('--zero-set',action='store_true')
    parser.add_argument('--cap-weight',action='store_true')
    args = parser.parse_args()
    for active in args.active:
        print(json.dumps(search(active,args.k_log,args.zero_set,args.cap_weight)),flush=True)
