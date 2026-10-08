"""Proof-only Cauchy width8/state7,8,12,16 GF256 screens.

These changed inners use non-self-orthogonal A,C. Invertibility of the
extended packet map replaces CA=0. The full-state transfer deliberately
includes all-output-zero nonreturn terms when there are active inputs.
No result here is a whole-code certificate or a native performance claim.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from fractions import Fraction as F
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


WIDTH = 8
RETURN_MINIMA = {7: (8, 6, 5, 4, 3, 2, 2, 1), 8: (8, 7, 6, 5, 4, 3, 2, 1),
                 12: (5, 5, 5, 4, 3, 3, 2, 1), 16: (8, 7, 6, 5, 4, 3, 2, 1)}


def validate_state(state):
    if state not in RETURN_MINIMA:
        raise ValueError('this screen fixes width8 and state7,8,12,16')


@lru_cache(maxsize=4)
def maps(state=8):
    validate_state(state)
    a = tuple(tuple(INV[u^v] for v in range(128, 128+state)) for u in range(WIDTH))
    c = tuple(tuple(INV[v^j] for j in range(64, 64+WIDTH)) for v in range(128, 128+state))
    return a, c


def validate_maps(state=8):
    a, c = maps(state)
    ac, ca = matmul(a, c), matmul(c, a)
    returned = tuple(tuple(int(i == j)^MUL[2][ac[i][j]] for j in range(WIDTH))
                     for i in range(WIDTH))
    schur = tuple(tuple(int(i == j)^MUL[2][ca[i][j]] for j in range(state))
                  for i in range(state))
    minima = exact_return_minima(returned)
    if rank(returned) != WIDTH or rank(schur) != state or minima != RETURN_MINIMA[state]:
        raise AssertionError('map validation failed')
    return dict(width=WIDTH, state=state, eta=2, return_rank=rank(returned),
                state_schur_rank=rank(schur), return_minima=list(minima),
                a_nodes=[list(range(WIDTH)), list(range(128, 128+state))],
                c_nodes=[list(range(128, 128+state)), list(range(64, 64+WIDTH))],
                ca_is_zero=not any(value for row in ca for value in row),
                whole_code=False)


@lru_cache(maxsize=128)
def gamma(tau, state):
    return tuple(max(rank_potential_bound(r, z, tau, 256, state)
                     for r in range(1, state+1)) for z in range(min(state-1, WIDTH)+1))


@lru_cache(maxsize=32)
def joint_terms(weight, state):
    return tuple((i, z, math.comb(weight, i)*math.comb(WIDTH-weight, z),
                  max(state, i+z, weight-WIDTH+i+2*z+1))
                 for i in range(weight+1)
                 for z in range(min(state-1, WIDTH-weight)+1)
                 if WIDTH-i-z >= RETURN_MINIMA[state][weight-1])


def return_bound(weight, y, tau, state=8):
    minimum = RETURN_MINIMA[state][weight-1]
    flats = sum(math.comb(state, k)/tau**k for k in range(min(weight-1, state-1)+1))
    pure = flats*y**minimum/tau**state
    joint = sum(count*y**(WIDTH-i-z)*(1-y)**(i+z)/tau**charge
                *sum(math.comb(state, k)/tau**k for k in range(min(weight-1, state-1-z)+1))
                for i, z, count, charge in joint_terms(weight, state))
    return min(pure, joint)


def local_matrix(weight, y, tau, state=8):
    validate_state(state)
    if not 0 <= weight <= WIDTH or not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('invalid local marks')
    y, tau = (F(v) if isinstance(v, int) else v for v in (y, tau))
    zero = y*0+tau*0
    result = [[zero, zero], [zero, zero]]
    if weight == 0:
        result[0][0] = y**0
    else:
        result[0][1] = y**weight
        if weight > state:
            result[0][0] = y**weight/tau**state
        result[1][0] = return_bound(weight, y, tau, state)
    roots = gamma(tau, state)
    result[1][1] = sum(math.comb(weight, i)*math.comb(WIDTH-weight, z)
                       *y**(WIDTH-i-z)*(1-y)**(i+z)*roots[z]/tau**i
                       for i in range(weight+1)
                       for z in range(min(state-1, WIDTH-weight)+1))
    return result


@lru_cache(maxsize=4)
def return_arrays(state):
    return tuple(np.asarray([(WIDTH-i-z, i+z, charge, min(weight-1, state-1-z), count)
                             for i, z, count, charge in joint_terms(weight, state)], dtype=float).T
                 for weight in range(1, WIDTH+1))


def packet(logx, logy, logtau, state=8):
    validate_state(state)
    x, y, tau = math.exp(logx), math.exp(logy), min(255., math.exp(logtau))
    powers_y = np.asarray([y**i for i in range(WIDTH+1)])
    powers_z = np.asarray([(1-y)**i for i in range(WIDTH+1)])
    powers_tau = np.asarray([tau**(-i) for i in range(2*WIDTH+1)])
    flats = np.cumsum([math.comb(state, k)*powers_tau[k] for k in range(state)])
    births = [math.comb(WIDTH, w)*(x*y)**w for w in range(WIDTH+1)]
    result = np.zeros((2, 2))
    result[0, 0] = 1+sum(births[state+1:])*powers_tau[state]
    result[0, 1] = sum(births[1:])
    base = y+x*y+x*(1-y)/tau
    result[1, 1] = sum(math.comb(WIDTH, z)*(1-y)**z*base**(WIDTH-z)*gamma(tau, state)[z]
                       for z in range(min(state-1, WIDTH)+1))
    for weight, arrays in enumerate(return_arrays(state), 1):
        yp, zp, charges, caps, counts = arrays
        joint = np.sum(counts*powers_y[yp.astype(int)]*powers_z[zp.astype(int)]
                       *powers_tau[charges.astype(int)]*flats[caps.astype(int)])
        pure = flats[min(weight-1, state-1)]*powers_y[RETURN_MINIMA[state][weight-1]]*powers_tau[state]
        result[1, 0] += math.comb(WIDTH, weight)*x**weight*min(pure, joint)
    return result


def middle_bound(point, weight, state=8, dimension=96, k_log=20, components=False):
    g = geometry(k_log, 128, dimension)
    logx, logy, logtau, logtilt = map(float, point)
    terms = dict(route=128*modal_log_penalty(g['rows']),
        source=g['rows']*source_row_log(weight/g['n'], logtau, logtilt, 128, dimension, True),
        inner=log_terminal_power(packet(logx, logy, logtau, state), g['n']//WIDTH),
        selection=-weight*(logx+logtilt), output=-g['threshold']*logy,
        rank_threshold=-logtau)
    value = sum(terms.values())/math.log(2)
    return (value, {k: v/math.log(2) for k, v in terms.items()}) if components else value


def sparse_bound(active, point, state=8, dimension=96, k_log=20, components=False):
    g = geometry(k_log, 128, dimension)
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    local = np.asarray([local_matrix(w, math.exp(logy), tau, state)
                        for w in range(min(WIDTH, active)+1)])
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
    terms = dict(source=log_choose(g['rows'], active)-((128-dimension)*active+1)*logtau,
                 output=-g['threshold']*logy, regional=log_terminal(regional, 128))
    value = sum(terms.values())/math.log(2)
    return (value, {k: v/math.log(2) for k, v in terms.items()}) if components else value


def search(value, mode='middle', state=8, dimension=96, k_log=20, iterations=60, start=None):
    g = geometry(k_log, 128, dimension)
    if mode == 'middle':
        target = max(1, round(value*g['n']))
        objective = lambda p: middle_bound(p, target, state, dimension, k_log)
        starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
        bounds = [(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)]
    elif mode == 'sparse':
        target = int(value)
        objective = lambda p: sparse_bound(target, p, state, dimension, k_log)
        starts = [(math.log(640), math.log(255)), (math.log(640), math.log(8))]
        bounds = [(-2, 10), (0, math.log(255))]
    else:
        raise ValueError('unknown screen mode')
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(objective, point, method='L-BFGS-B', bounds=bounds,
                       options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for point in starts]
    best = min(choices, key=lambda r: r.fun if math.isfinite(r.fun) else math.inf)
    bound = sparse_bound if mode == 'sparse' else middle_bound
    _, terms = bound(target, best.x, state, dimension, k_log, True) if mode == 'sparse' else bound(
        best.x, target, state, dimension, k_log, True)
    return dict(outer=[128, dimension], k_log=k_log, width=WIDTH, state=state,
        eta=2, return_minima=list(RETURN_MINIMA[state]), mode=mode, target=target,
        rows=g['rows'], n=g['n'], margin_bits=-float(best.fun),
        log_parameters=best.x.tolist(), log2_components=terms,
        optimizer_success=bool(best.success), iterations=int(best.nit),
        rigorous=False, whole_code=False,
        scope='one floating weight-class bound for a changed non-self-orthogonal inner')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=int, choices=(7, 8, 12, 16), nargs='+', default=[8])
    parser.add_argument('--mode', choices=('sparse', 'middle'), default='middle')
    parser.add_argument('--active', type=int, nargs='+', default=[64, 128])
    parser.add_argument('--density', type=float, nargs='+', default=[.03, .10, .14, .20, .24, .28])
    parser.add_argument('--dimension', type=int, default=96)
    parser.add_argument('--iterations', type=int, default=60)
    parser.add_argument('--check-maps', action='store_true')
    args = parser.parse_args()
    for state in args.state:
        checked = validate_maps(state)
        if args.check_maps:
            print(json.dumps(checked), flush=True)
        else:
            previous = None
            for value in args.active if args.mode == 'sparse' else args.density:
                result = search(value, args.mode, state, args.dimension,
                                iterations=args.iterations, start=previous)
                previous = result['log_parameters']
                print(json.dumps(result), flush=True)
