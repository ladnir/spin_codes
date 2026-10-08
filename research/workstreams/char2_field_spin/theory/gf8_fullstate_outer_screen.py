"""Outer-length screens for the non-self-orthogonal width-eight inner.

Uses the separately validated state7/state8 local transfer. Outer changes
remain GF256 MDS codes of length at most256. Results are floating point
class screens, not complete distance certificates.
"""
from __future__ import annotations

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from gf8_fullstate_screen import local_matrix, packet, validate_state
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


def profile(length=256, dimension=192, k_log=20):
    if not 1 <= dimension < length <= 256:
        raise ValueError('require a GF256 MDS outer with dimension<length<=256')
    return geometry(k_log, length, dimension)


def middle_bound(point, weight, length=256, dimension=192, state=8, k_log=20):
    g = profile(length, dimension, k_log)
    if not 0 < weight < g['n']:
        raise ValueError('middle weight outside output length')
    logx, logy, logtau, logtilt = map(float, point)
    return (length*modal_log_penalty(g['rows'])
        + g['rows']*source_row_log(weight/g['n'], logtau, logtilt, length, dimension, True)
        + log_terminal_power(packet(logx, logy, logtau, state), g['n']//8)
        - weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def sparse_bound(active, point, length=256, dimension=192, state=8, k_log=20):
    g = profile(length, dimension, k_log)
    if not 1 <= active <= g['rows']:
        raise ValueError('active rows outside geometry')
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    values = np.asarray([local_matrix(w, math.exp(logy), tau, state)
                         for w in range(min(8, active)+1)], dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(values)
    for w in range(len(factor)):
        factor[w] += math.log(math.comb(8, w))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    packets = g['rows']//8
    while packets:
        if packets & 1:
            result = log_polynomial_product(result, factor, active)
        packets >>= 1
        if packets:
            factor = log_polynomial_product(factor, factor, active)
    regional = regional_from_coefficients(result, active, g['rows'], logtau)
    return (log_choose(g['rows'], active)-((length-dimension)*active+1)*logtau
            -g['threshold']*logy+log_terminal(regional, length))/math.log(2)


def search(value, mode='middle', length=256, dimension=192, state=8,
           k_log=20, iterations=80, start=None):
    validate_state(state)
    g = profile(length, dimension, k_log)
    if mode == 'middle':
        target = round(value*g['n'])
        objective = lambda p: middle_bound(p, target, length, dimension, state, k_log)
        starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
        bounds = [(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)]
    elif mode == 'sparse':
        target = int(value)
        objective = lambda p: sparse_bound(target, p, length, dimension, state, k_log)
        starts = [(math.log(640), math.log(255)), (math.log(640), math.log(8))]
        bounds = [(-2, 10), (0, math.log(255))]
    else:
        raise ValueError('unknown mode')
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(objective, point, method='L-BFGS-B', bounds=bounds,
                       options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25})
               for point in starts]
    best = min(choices, key=lambda r: r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[length, dimension], width=8, state=state, k_log=k_log,
        requested_k=2**k_log, padded_k=g['k'], rows=g['rows'], n=g['n'],
        output_threshold=g['threshold'], mode=mode, target=target,
        margin_bits=-float(best.fun), log_parameters=best.x.tolist(),
        optimizer_success=bool(best.success), iterations=int(best.nit),
        rigorous=False, whole_code=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, default=(256, 192))
    parser.add_argument('--state', type=int, choices=(7, 8), default=8)
    parser.add_argument('--mode', choices=('middle', 'sparse'), default='middle')
    parser.add_argument('--density', type=float, nargs='+', default=[.03, .10, .14])
    parser.add_argument('--active', type=int, nargs='+', default=[32, 64, 128])
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--iterations', type=int, default=80)
    args = parser.parse_args()
    previous = None
    for value in args.density if args.mode == 'middle' else args.active:
        result = search(value, args.mode, *args.outer, args.state, args.k_log,
                        args.iterations, previous)
        previous = result['log_parameters']
        print(json.dumps(result), flush=True)
