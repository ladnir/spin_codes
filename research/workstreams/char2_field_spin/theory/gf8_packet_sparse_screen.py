"""Floating active-row screens for candidate GF256 packet geometries.

Uses the same positive regional averaging as joint_outer_screen, with the
explicit packet dimensions and proved return minima supplied by the caller.
Every result covers one active-row class only, not the whole code. No native
timing is performed. This module never substitutes a new map into a frozen
certificate.
"""
from __future__ import annotations

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose
from gf8_packet_width_screen import local_matrix, validate
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


def log_coefficients(active, y, tau, rows, width, state, minima=()):
    validate(width, state)
    if rows <= 0 or rows % width or not 0 <= active <= rows:
        raise ValueError('invalid packet-aligned region geometry')
    local = np.asarray([local_matrix(w, y, tau, width, state, minima)
                        for w in range(min(width, active)+1)], dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    for weight in range(len(local)):
        factor[weight] += math.log(math.comb(width, weight))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    packets = rows//width
    while packets:
        if packets & 1:
            result = log_polynomial_product(result, factor, active)
        packets >>= 1
        if packets:
            factor = log_polynomial_product(factor, factor, active)
    return result


def bound(active, point, k_log=20, dimension=84, width=8, state=4,
          minima=(), components=False):
    g = geometry(k_log, 128, dimension)
    if not 1 <= active <= g['rows'] or not 1 <= dimension < 128:
        raise ValueError('invalid outer or active-row class')
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    coefficients = log_coefficients(active, math.exp(logy), tau, g['rows'],
                                    width, state, minima)
    regional = regional_from_coefficients(coefficients, active, g['rows'], math.log(tau))
    terms = dict(source=log_choose(g['rows'], active)-((128-dimension)*active+1)*math.log(tau),
                 output=-g['threshold']*logy, regional=log_terminal(regional, 128))
    value = sum(terms.values())/math.log(2)
    return (value, {key: val/math.log(2) for key, val in terms.items()}) if components else value


def search(active, k_log=20, dimension=84, width=8, state=4, minima=(),
           iterations=40, start=None):
    starts = [(math.log(640), math.log(255)), (math.log(640), math.log(8))]
    if start is not None:
        starts.insert(0, start)
    results = [minimize(lambda p: bound(active, p, k_log, dimension, width, state, minima),
                        initial, method='L-BFGS-B', bounds=[(-2, 10), (0, math.log(255))],
                        options={'maxiter': iterations, 'ftol': 1e-10, 'maxls': 16})
               for initial in starts]
    best = min(results, key=lambda result: result.fun if np.isfinite(result.fun) else math.inf)
    value, terms = bound(active, best.x, k_log, dimension, width, state, minima, True)
    g = geometry(k_log, 128, dimension)
    return dict(outer=[128, dimension], width=width, state=state, k_log=k_log,
                requested_k=2**k_log, padded_k=g['k'], rows=g['rows'], n=g['n'],
                output_threshold=g['threshold'], distance='1/10', active=active,
                margin_bits=-value, log_parameters=best.x.tolist(), log2_components=terms,
                return_minima=list(minima), optimizer_success=bool(best.success),
                iterations=int(best.nit), rigorous=False, whole_code=False,
                scope='one active-row class, evaluated in floating point')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width', type=int, default=8)
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--dimension', type=int, nargs='+', default=[84, 85])
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--active', type=int, nargs='+', default=[1, 2, 4, 8, 16, 32, 64, 128, 256])
    parser.add_argument('--return-minima', type=int, nargs='+', default=[])
    parser.add_argument('--iterations', type=int, default=40)
    args = parser.parse_args()
    for dimension in args.dimension:
        previous = None
        for active in args.active:
            result = search(active, args.k_log, dimension, args.width, args.state,
                            tuple(args.return_minima), args.iterations, previous)
            previous = result['log_parameters']
            print(json.dumps(result), flush=True)
