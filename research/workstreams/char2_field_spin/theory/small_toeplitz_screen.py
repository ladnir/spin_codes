"""Floating screens for smaller Toeplitz outers and eight-lane GF16 inners.

Every outer row uses an independent uniform systematic [2d,d] Toeplitz
matrix. There is no input diagonal. The expected message-value identity
in K20_RANK_CONTINUATION.md supplies the product value-counting measure.
Each coordinate has one independently permuted region; state never resets.

A_j=(1,j,...,j^(s-1)), C=A^T, j=0,...,t-1 in GF16 modulo X^4+X+1.
The model computes CA and I+AC rather than assuming the sixteen-lane law.
All numerical bounds below are floating screens, not certificates.
"""

import argparse
from collections import Counter
from functools import lru_cache
import json
import math

import numpy as np

from large_state_screen import mul, rank, xor_sum
from state_screen import local as generic_local


def geometry(k_log, dimension, width):
    if dimension < 1 or width not in (8, 16):
        raise ValueError('positive outer dimension and width8/16 required')
    multiple = dimension*width
    k = (((1 << k_log)+multiple-1)//multiple)*multiple
    n = 2*k
    return k, n, k//dimension, (n+9)//10-1


@lru_cache(maxsize=None)
def algebra(width, state):
    if width not in (8, 16) or not 1 <= state <= width:
        raise ValueError('unsupported inner geometry')
    a = []
    for point in range(width):
        row = [1]
        for _ in range(1, state):
            row.append(mul(row[-1], point))
        a.append(row)
    gram = [[xor_sum(mul(row[i], row[j]) for row in a) for j in range(state)] for i in range(state)]
    b = [[int(i == j) ^ xor_sum(mul(x, y) for x, y in zip(a[i], a[j]))
          for j in range(width)] for i in range(width)]
    weights = [sum(row[j] != 0 for row in b) for j in range(width)]
    return dict(width=width, state=state, ca_zero=not any(any(row) for row in gram),
                ca=gram, return_rank=rank(b), singleton_weights=weights,
                singleton_histogram=dict(sorted(Counter(weights).items())))


def local(z, width, state):
    info = algebra(width, state)
    if info['return_rank'] != width:
        raise ValueError('invertible-B bound is unavailable for this map')
    result = generic_local(z, width, state, 16)
    denominator = 16**state-1
    joint = sum(z**w for w in info['singleton_weights'])/(width*denominator)
    result[1, 1, 1] += z**width/denominator-joint
    result[1, 1, 0] = joint
    for w in range(2, width+1):
        inverse = sum(math.comb(width, h)*z**(width-h)*(1-z)**h
                      /15**max(0, h-width+w) for h in range(width+1))/denominator
        result[w, 1, 0] = min(result[w, 1, 0], inverse)
    return result


def screen(k_log, active, state, dimension, width, z, exact_q1=True):
    from reduced_distance_certificate import log_polynomial_product
    from scipy.special import logsumexp
    _, _, rows, threshold = geometry(k_log, dimension, width)
    with np.errstate(divide='ignore'):
        factors = np.log(local(z, width, state)[:active+1])
    for w in range(len(factors)):
        factors[w] += math.log(math.comb(width, w))
    regions = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows//width
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, factors, active)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, active)
    for h in range(active+1):
        regions[h] -= math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
    source = np.array([math.log(math.comb(active, h))+h*math.log(15)
                       for h in range(active+1)])
    matrix = logsumexp(regions+source[:, None, None], axis=0)
    if active == 1 and exact_q1:
        # Positive flag DP keeps the nonempty systematic support constraint.
        # Only the first d regions correspond to systematic coordinates.
        empty, seen = np.array([0., -math.inf]), np.array([-math.inf, -math.inf])
        active_matrix = regions[1]+math.log(15)
        for _ in range(dimension):
            seen = np.logaddexp(logsumexp(seen[:, None]+matrix, axis=0),
                                logsumexp(empty[:, None]+active_matrix, axis=0))
            empty = logsumexp(empty[:, None]+regions[0], axis=0)
        vector = seen
        for _ in range(dimension):
            vector = logsumexp(vector[:, None]+matrix, axis=0)
    else:
        vector = np.array([0., -math.inf])
        for _ in range(2*dimension):
            vector = logsumexp(vector[:, None]+matrix, axis=0)
    moment = float(logsumexp(vector))-dimension*active*math.log(16)
    choose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    return -(choose-math.log(15)-threshold*math.log(z)+moment)/math.log(2)


def search(k_log, active, state, dimension, width):
    from scipy.optimize import minimize_scalar
    k, n, rows, _ = geometry(k_log, dimension, width)
    def objective(logtilt):
        z = math.exp(-math.exp(logtilt)*active/n)
        return -screen(k_log, active, state, dimension, width, z)
    optimum = minimize_scalar(objective, bounds=(1, 10), method='bounded', options={'xatol': 2e-6})
    return dict(k=k, n=n, rows=rows, active_rows=active, state=state, width=width,
                outer=[2*dimension, dimension], margin_bits=-float(optimum.fun),
                z=math.exp(-math.exp(optimum.x)*active/n), rigorous=False,
                source='exact-nonempty-systematic' if active == 1 else 'positive-source-envelope')


def middle_screen(k_log, active, state, dimension, width, x, z):
    from coupled_chernoff import log_terminal_power
    _, _, rows, threshold = geometry(k_log, dimension, width)
    transfer = local(z, width, state)
    matrix = sum(math.comb(width, w)*(15*x)**w*(1+x)**(width-w)*transfer[w] for w in range(width+1))
    choose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    value = (-(2*dimension-1)*choose-dimension*active*math.log(16)-math.log(15)
             -threshold*math.log(z)-2*dimension*active*math.log(x)
             + log_terminal_power(matrix, rows*2*dimension//width))
    return -value/math.log(2)


def middle_search(k_log, active, state, dimension, width):
    from scipy.optimize import minimize
    k, n, rows, _ = geometry(k_log, dimension, width)
    def objective(point):
        logx, logz = point
        if not -35 < logx < 10 or not -20 < logz < 0:
            return 1e100
        return -middle_screen(k_log, active, state, dimension, width, math.exp(logx), math.exp(logz))
    best = None
    for initial_z in (.5, .9, .99, .999):
        initial_x = active/(max(1, rows-active)*16)
        optimum = minimize(objective, [math.log(initial_x), math.log(initial_z)], method='Nelder-Mead',
                           options={'maxiter': 300, 'xatol': 1e-9, 'fatol': 1e-7})
        candidate = dict(k=k, n=n, rows=rows, active_rows=active, state=state, width=width,
                         outer=[2*dimension, dimension], margin_bits=-float(optimum.fun),
                         x=math.exp(optimum.x[0]), z=math.exp(optimum.x[1]), rigorous=False)
        if best is None or candidate['margin_bits'] > best['margin_bits']:
            best = candidate
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width', type=int, choices=(8, 16), default=8)
    parser.add_argument('--dimension', type=int, default=16)
    parser.add_argument('--state', nargs='+', type=int, default=[2, 3, 4])
    parser.add_argument('--active', nargs='+', type=int, default=[1])
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--middle', action='store_true')
    args = parser.parse_args()
    for state in args.state:
        print(json.dumps(algebra(args.width, state)), flush=True)
        for active in args.active:
            method = middle_search if args.middle else search
            print(json.dumps(method(args.k_log, active, state, args.dimension, args.width)), flush=True)


if __name__ == '__main__':
    main()
