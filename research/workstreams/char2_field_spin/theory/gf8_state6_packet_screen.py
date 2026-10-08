"""Bounded proof screen for a changed width8/state6 GF256 construction.

A[u,j]=1/(u+(128+j)), eta=2 in the AES field. The default C=A^T;
the asymmetric option uses C[v,j]=1/(v+(64+j)) for v=128..133.
Both variants are superregular and have an invertible extended packet map,
but CA is nonzero.
This changes the construction and has no native implementation or full cert.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from diagonal_refresh_rank import rank_potential_bound
from gf8_packet_maps import INV, MUL, exact_return_minima, matmul, rank
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


WIDTH, STATE = 8, 6
MINIMA = (7, 7, 3, 4, 3, 3, 1, 2)
ASYMMETRIC_MINIMA = (8, 7, 5, 4, 3, 2, 1, 1)


@lru_cache(maxsize=2)
def maps(asymmetric=False):
    a = tuple(tuple(INV[u^v] for v in range(128, 134)) for u in range(8))
    c = (tuple(tuple(INV[v^j] for j in range(64, 72)) for v in range(128, 134))
         if asymmetric else tuple(zip(*a)))
    return a, c


def validate_maps(asymmetric=False):
    a, c = maps(asymmetric)
    ac, ca = matmul(a, c), matmul(c, a)
    returned = tuple(tuple(int(i == j)^MUL[2][ac[i][j]] for j in range(8)) for i in range(8))
    schur = tuple(tuple(int(i == j)^MUL[2][ca[i][j]] for j in range(6)) for i in range(6))
    minima = exact_return_minima(returned)
    expected = ASYMMETRIC_MINIMA if asymmetric else MINIMA
    if rank(returned) != WIDTH or rank(schur) != STATE or minima != expected:
        raise AssertionError('fixed state6 map validation failed')
    return dict(width=WIDTH, state=STATE, eta=2, return_rank=rank(returned),
                state_schur_rank=rank(schur), return_minima=list(minima),
                ca_is_zero=all(value == 0 for row in ca for value in row), asymmetric=asymmetric,
                superregular_reason=('A has disjoint nodes0..7/128..133; C has128..133/64..71'
                                     if asymmetric else 'disjoint Cauchy nodes0..7 and128..133; C is transpose'),
                whole_code=False)


@lru_cache(maxsize=64)
def gamma(tau):
    return tuple(max(rank_potential_bound(r, z, tau, 256, STATE)
                     for r in range(1, STATE+1)) for z in range(STATE))


@lru_cache(maxsize=16)
def terms(weight, minima=MINIMA):
    return tuple((i, z, math.comb(weight, i)*math.comb(WIDTH-weight, z),
                  max(STATE, i+z, weight-WIDTH+i+2*z+1))
                 for i in range(weight+1)
                 for z in range(min(STATE-1, WIDTH-weight)+1)
                 if WIDTH-i-z >= minima[weight-1])


def return_bound(weight, y, tau, minima=MINIMA):
    flats = sum(math.comb(STATE, k)/tau**k for k in range(min(weight-1, STATE-1)+1))
    pure = flats*y**minima[weight-1]/tau**STATE
    joint = sum(count*y**(WIDTH-i-z)*(1-y)**(i+z)/tau**charge
                *sum(math.comb(STATE, k)/tau**k for k in range(min(weight-1, STATE-1-z)+1))
                for i, z, count, charge in terms(weight, minima))
    return min(pure, joint)


def local_matrix(weight, y, tau, minima=MINIMA):
    if not 0 <= weight <= WIDTH or not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('invalid state6 packet marks')
    y, tau = (F(x) if isinstance(x, int) else x for x in (y, tau))
    zero = y*0+tau*0
    result = [[zero, zero], [zero, zero]]
    if weight == 0:
        result[0][0] = y**0
    else:
        result[0][1] = y**weight
        if weight > STATE:
            result[0][0] = y**weight/tau**STATE
        result[1][0] = return_bound(weight, y, tau, minima)
    roots = gamma(tau)
    result[1][1] = sum(math.comb(weight, i)*math.comb(WIDTH-weight, z)
                       *y**(WIDTH-i-z)*(1-y)**(i+z)*roots[z]/tau**i
                       for i in range(weight+1)
                       for z in range(min(STATE-1, WIDTH-weight)+1))
    return result


@lru_cache(maxsize=2)
def return_arrays(minima=MINIMA):
    return tuple(np.asarray([(WIDTH-i-z, i+z, charge, min(weight-1, STATE-1-z), count)
                             for i, z, count, charge in terms(weight, minima)], dtype=float).T
                 for weight in range(1, WIDTH+1))


def packet(logx, logy, logtau, minima=MINIMA):
    x, y, tau = math.exp(logx), math.exp(logy), min(255., math.exp(logtau))
    powers_y = np.asarray([y**i for i in range(WIDTH+1)])
    powers_z = np.asarray([(1-y)**i for i in range(WIDTH+1)])
    powers_tau = np.asarray([tau**(-i) for i in range(2*WIDTH+1)])
    flats = np.cumsum([math.comb(STATE, k)*powers_tau[k] for k in range(STATE)])
    births = [math.comb(WIDTH, w)*(x*y)**w for w in range(WIDTH+1)]
    result = np.zeros((2, 2))
    result[0, 0] = 1+sum(births[STATE+1:])*powers_tau[STATE]
    result[0, 1] = sum(births[1:])
    base = y+x*y+x*(1-y)/tau
    result[1, 1] = sum(math.comb(WIDTH, z)*(1-y)**z*base**(WIDTH-z)*gamma(tau)[z]
                       for z in range(STATE))
    for weight, arrays in enumerate(return_arrays(minima), 1):
        yp, zp, charges, caps, counts = arrays
        joint = np.sum(counts*powers_y[yp.astype(int)]*powers_z[zp.astype(int)]
                       *powers_tau[charges.astype(int)]*flats[caps.astype(int)])
        pure = flats[min(weight-1, STATE-1)]*powers_y[minima[weight-1]]*powers_tau[STATE]
        result[1, 0] += math.comb(WIDTH, weight)*x**weight*min(pure, joint)
    return result


def middle_bound(point, weight, dimension=96, k_log=20, minima=MINIMA):
    g = geometry(k_log, 128, dimension)
    logx, logy, logtau, logtilt = map(float, point)
    return (128*modal_log_penalty(g['rows'])
        + g['rows']*source_row_log(weight/g['n'], logtau, logtilt, 128, dimension, True)
        + log_terminal_power(packet(logx, logy, logtau, minima), g['n']//WIDTH)
        - weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def sparse_bound(active, point, dimension=96, k_log=20, minima=MINIMA):
    g = geometry(k_log, 128, dimension)
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    local = np.asarray([local_matrix(w, math.exp(logy), tau, minima) for w in range(min(WIDTH, active)+1)])
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    for w in range(len(factor)):
        factor[w] += math.log(math.comb(WIDTH, w))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    packets = g['rows']//WIDTH
    while packets:
        if packets & 1:
            result = log_polynomial_product(result, factor, active)
        packets >>= 1
        if packets:
            factor = log_polynomial_product(factor, factor, active)
    regional = regional_from_coefficients(result, active, g['rows'], logtau)
    return (log_choose(g['rows'], active)-((128-dimension)*active+1)*logtau
            -g['threshold']*logy+log_terminal(regional, 128))/math.log(2)


def search(value, mode='middle', dimension=96, k_log=20, iterations=60, start=None, asymmetric=False):
    g = geometry(k_log, 128, dimension)
    minima = ASYMMETRIC_MINIMA if asymmetric else MINIMA
    if mode == 'middle':
        target = max(1, round(value*g['n']))
        objective = lambda p: middle_bound(p, target, dimension, k_log, minima)
        starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
        bounds = [(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)]
    elif mode == 'sparse':
        target = int(value)
        objective = lambda p: sparse_bound(target, p, dimension, k_log, minima)
        starts = [(math.log(640), math.log(255)), (math.log(640), math.log(8))]
        bounds = [(-2, 10), (0, math.log(255))]
    else:
        raise ValueError('unknown screen mode')
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(objective, point, method='L-BFGS-B', bounds=bounds,
                       options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for point in starts]
    best = min(choices, key=lambda r: r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[128, dimension], k_log=k_log, width=WIDTH, state=STATE,
        eta=2, return_minima=list(minima), mode=mode, target=target, asymmetric=asymmetric,
        margin_bits=-float(best.fun), log_parameters=best.x.tolist(),
        optimizer_success=bool(best.success), iterations=int(best.nit),
        rigorous=False, whole_code=False,
        scope='one weight-class floating screen for the changed Cauchy8/6 inner')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('sparse', 'middle'), default='middle')
    parser.add_argument('--active', type=int, nargs='+', default=[64])
    parser.add_argument('--density', type=float, nargs='+', default=[.10, .14])
    parser.add_argument('--dimension', type=int, default=96)
    parser.add_argument('--iterations', type=int, default=60)
    parser.add_argument('--check-maps', action='store_true')
    parser.add_argument('--asymmetric', action='store_true')
    args = parser.parse_args()
    if args.check_maps:
        print(json.dumps(validate_maps(args.asymmetric)), flush=True)
    else:
        previous = None
        for value in args.active if args.mode == 'sparse' else args.density:
            result = search(value, args.mode, args.dimension, iterations=args.iterations, start=previous,
                            asymmetric=args.asymmetric)
            previous = result['log_parameters']
            print(json.dumps(result), flush=True)
