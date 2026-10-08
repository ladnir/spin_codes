"""Two-state full-Cauchy packet geometry screens; not certificates.

Uses a Cauchy matrix of size width+state with invertibly normalized state
rows. Retains the conservative zero-flat return cover as a cheap geometry
gate. No concrete matrix or native implementation is claimed by this file.
"""
from functools import lru_cache
import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from diagonal_refresh_rank import rank_potential_bound
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


@lru_cache(maxsize=256, typed=True)
def roots(tau, width, state):
    return tuple(max(rank_potential_bound(r, z, tau, 256, state)
                     for r in range(1, state+1)) for z in range(min(width, state-1)+1))


def local_matrices(y, tau, width=16, state=8):
    if not (width in (4, 8, 16, 32) and 1 <= state <= 32
            and 0 < y <= 1 and 1 <= tau <= 255):
        raise ValueError('invalid full-Cauchy geometry or marks')
    zero = y*0+tau*0
    result = [[[zero, zero], [zero, zero]] for _ in range(width+1)]
    result[0][0][0] = y**0
    gamma = roots(tau, width, state)
    flats = []
    value = zero
    for k in range(state):
        value += math.comb(state, k)/tau**k
        flats.append(value)
    terms = [math.comb(width, z)*y**(width-z)*(1-y)**z for z in range(width+1)]
    for w in range(width+1):
        for z in range(width+1):
            pivots, remaining = min(w, z), max(0, z-w)
            if remaining < state:
                result[w][1][1] += terms[z]*gamma[remaining]/tau**pivots
            if w and z < w:
                result[w][0][1] += terms[z]/tau**z
            if w > state+z:
                result[w][0][0] += terms[z]/tau**(state+z)
        if w:
            result[w][0][1] = min(result[w][0][1], y**(width+1-w))
            pure = y**(width+1-w)/tau**state*flats[min(w-1, state-1)]
            joint = sum(terms[z]/tau**(state+z)*flats[min(w-z-1, state-1)]
                        for z in range(w))
            result[w][1][0] = min(pure, joint)
    return result


def packet(point, width, state):
    logx, logy, logtau = point[:3]
    tau = min(255., math.exp(logtau))
    local = np.asarray(local_matrices(math.exp(logy), tau, width, state), dtype=float)
    return np.einsum('w,wij->ij',
        [math.comb(width, w)*math.exp(w*logx) for w in range(width+1)], local)


def middle_bound(point, weight, width=16, state=8, length=128, dimension=96, k_log=20):
    g = geometry(k_log, length, dimension)
    logx, logy, logtau, logtilt = map(float, point)
    return (length*modal_log_penalty(g['rows'])
        + g['rows']*source_row_log(weight/g['n'], logtau, logtilt, length, dimension, True)
        + log_terminal_power(packet(point, width, state), g['n']//width)
        - weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def sparse_bound(active, point, width=16, state=8, length=128, dimension=96, k_log=20):
    g = geometry(k_log, length, dimension)
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    local = np.asarray(local_matrices(math.exp(logy), tau, width, state)[:min(width, active)+1], dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    for w in range(len(factor)):
        factor[w] += math.log(math.comb(width, w))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    steps = g['rows']//width
    while steps:
        if steps & 1:
            result = log_polynomial_product(result, factor, active)
        steps >>= 1
        if steps:
            factor = log_polynomial_product(factor, factor, active)
    regional = regional_from_coefficients(result, active, g['rows'], logtau)
    return (log_choose(g['rows'], active)-((length-dimension)*active+1)*logtau
        -g['threshold']*logy+log_terminal(regional, length))/math.log(2)


def search(value, sparse=False, width=16, state=8, length=128, dimension=96,
           iterations=70, start=None):
    g = geometry(20, length, dimension)
    if sparse:
        target = int(value)
        objective = lambda p: sparse_bound(target, p, width, state, length, dimension)
        starts = [(6.5, 2.), (7., 3.)]
        bounds = [(-2, 10), (0, math.log(255))]
    else:
        target = round(value*g['n'])
        objective = lambda p: middle_bound(p, target, width, state, length, dimension)
        starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
        bounds = [(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)]
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(objective, p, method='L-BFGS-B', bounds=bounds,
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for p in starts]
    best = min(choices, key=lambda r: r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[length, dimension], width=width, state=state, k_log=20,
        rows=g['rows'], n=g['n'], mode='sparse' if sparse else 'middle', target=target,
        margin_bits=-float(best.fun), log_parameters=best.x.tolist(),
        optimizer_success=bool(best.success), iterations=int(best.nit),
        rigorous=False, whole_code=False, concrete_matrix_validated=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width', type=int, default=16)
    parser.add_argument('--state', type=int, default=8)
    parser.add_argument('--outer', type=int, nargs=2, default=(128,96))
    parser.add_argument('--sparse', type=int, nargs='*', default=[64,128])
    parser.add_argument('--density', type=float, nargs='*', default=[.10,.14])
    parser.add_argument('--iterations', type=int, default=70)
    args = parser.parse_args()
    for sparse, values in ((True,args.sparse),(False,args.density)):
        previous = None
        for value in values:
            result = search(value,sparse,args.width,args.state,*args.outer,args.iterations,previous)
            previous = result['log_parameters']
            print(json.dumps(result),flush=True)
