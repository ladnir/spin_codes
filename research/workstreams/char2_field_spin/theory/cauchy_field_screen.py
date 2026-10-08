"""Exploratory middle-row screens using the exact Cauchy packet transfer.

Fixed GRS[b,k] outer, iid nonzero F_q input diagonals, independent uniform
regional permutations, and fresh independent GL(s,q) state refresh at
every packet. One global routing stage. This is not a finite certificate.
"""

import argparse
import json
import math
from functools import lru_cache

import numpy as np
from scipy.optimize import minimize

from cauchy_inner_transfer import floating_transfer, validate
from coupled_chernoff import log_terminal_power
from field_rate_tradeoff import log_choose
from q1_batched_screen import fixed_support_multiplicity


@lru_cache(maxsize=None)
def support_corrections(length, dimension, q):
    parity = length-dimension
    return tuple((w, math.log(fixed_support_multiplicity(length, dimension, w, q))
                  -(w-parity)*math.log(q)) for w in range(parity+1, length+1))


def geometry(k_log, length, dimension, width):
    if not 0 < dimension < length or length > 256 or k_log < 1:
        raise ValueError('invalid outer or size')
    rows = (((1 << k_log)+dimension*width-1)//(dimension*width))*width
    return rows, rows*length, (rows*length+9)//10-1


def middle(k_log, active, x, z, rho, length=33, dimension=24, width=16, state=4, q=256):
    validate(width, state, q)
    rows, n, threshold = geometry(k_log, length, dimension, width)
    if not 1 <= active <= rows or x <= 0 or rho <= 0 or not 0 < z <= 1:
        raise ValueError('invalid witness')
    eta = max(value-w*math.log(rho) for w, value in support_corrections(length, dimension, q))
    coefficients = floating_transfer(z, width, state, q)
    matrix = sum(math.comb(width, w)*(q*rho*x)**w*(1+x)**(width-w)*coefficients[w]
                 for w in range(width+1))
    value = (-(length-1)*log_choose(rows, active)-(length-dimension)*active*math.log(q)
             +active*eta-math.log(q-1)-threshold*math.log(z)-length*active*math.log(x)
             +log_terminal_power(matrix, n//width))
    return -value/math.log(2)


def search(k_log, active, **parameters):
    width, dimension = parameters.get('width', 16), parameters.get('dimension', 24)
    rows, n, _ = geometry(k_log, parameters.get('length', 33), dimension, width)
    q = parameters.get('q', 256)
    best = None
    def objective(point):
        logx, logz, logrho = point
        if not -40 < logx < 14 or not -30 < logz < 0 or not -4 < logrho < 3:
            return 1e100
        return -middle(k_log, active, math.exp(logx), math.exp(logz), math.exp(logrho), **parameters)
    for initial_z in (.03, .3, .8, .99):
        result = minimize(objective, [math.log(active/(max(1, rows-active)*q)), math.log(initial_z), 0.],
                          method='Nelder-Mead', options={'maxiter': 500, 'fatol': 1e-6, 'xatol': 1e-9})
        if best is None or result.fun < best.fun:
            best = result
    return dict(k=rows*dimension, n=n, rows=rows, active=active, **parameters,
                margin_bits=-float(best.fun), x=math.exp(best.x[0]), z=math.exp(best.x[1]),
                rho=math.exp(best.x[2]), rigorous=False, scope='one-active-row-count-class-only')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--length', type=int, default=33)
    parser.add_argument('--dimension', type=int, default=24)
    parser.add_argument('--width', type=int, default=16)
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--fractions', type=float, nargs='+', default=[.01, .05, .1, .2, .4, .6, .8])
    args = parser.parse_args()
    parameters = dict(length=args.length, dimension=args.dimension, width=args.width, state=args.state, q=256)
    rows, _, _ = geometry(args.k_log, args.length, args.dimension, args.width)
    for fraction in args.fractions:
        print(json.dumps(search(args.k_log, max(1, int(rows*fraction)), **parameters)), flush=True)
