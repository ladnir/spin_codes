"""Joint active-row/total-weight positive screens, not certificates.

The source keeps the active-row count A and total support W. A separate
Bernoulli comparison bounds the regional route; it does not make the
actual outer coordinates or input forms independent. The full-Cauchy
local transfer is unchanged from gf8_mds_geometry_screen.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from gf8_mds_geometry_screen import local_matrices, packet
from gf8_rank_history_screen import modal_log_penalty
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


@lru_cache(maxsize=16)
def binomial_logs(length):
    return np.array([math.log(math.comb(length, h)) for h in range(length+1)])


def row_logs(alpha, logtau, loglambda, length=128, dimension=96):
    """Empty/occupied expectations after justified density homogenization.

    For beta=tau*lambda>=2, occupied rows have exact feasible support
    weights p+1..length. Below that threshold the geometric product is
    a valid relaxation. Both row factors separately increase under
    pairwise averaging of unequal regional Bernoulli probabilities.
    """
    if not 0 < alpha < 1 or not 0 < dimension < length:
        raise ValueError('invalid source geometry')
    p = length-dimension
    logzero, logbeta = math.log1p(-alpha), logtau+loglambda
    empty = length*logzero
    if logbeta < math.log(2):
        occupied = -p*logtau+length*np.logaddexp(logzero, logbeta+math.log(alpha))
    else:
        h = np.arange(p+1, length+1)
        occupied = logsumexp(binomial_logs(length)[p+1:]-p*logtau
                            +h*(logbeta+math.log(alpha))+(length-h)*logzero)
    return float(empty), float(occupied)


def source_joint_log(alpha, active, rows, logtau, loglambda, length=128, dimension=96):
    if not 0 <= active <= rows:
        raise ValueError('invalid active-row count')
    empty, occupied = row_logs(alpha, logtau, loglambda, length, dimension)
    return log_choose(rows, active)+(rows-active)*empty+active*occupied


def activation_mark_log(alpha, logmark, rows, logtau, loglambda, length=128, dimension=96):
    """All-A polynomial evaluated at a positive activation mark v."""
    empty, occupied = row_logs(alpha, logtau, loglambda, length, dimension)
    return rows*float(np.logaddexp(empty, logmark+occupied))


def feasible(active, weight, rows, length=128, dimension=96):
    return 1 <= active <= rows and (length-dimension+1)*active <= weight <= length*active


def middle_bound(point, active, weight, width=8, state=32, length=128, dimension=96, k_log=20):
    g = geometry(k_log, length, dimension)
    if not feasible(active, weight, g['rows'], length, dimension) or not 0 < weight < g['n']:
        return -math.inf
    logx, logy, logtau, loglambda = map(float, point)
    return (length*modal_log_penalty(g['rows'])
        +source_joint_log(weight/g['n'], active, g['rows'], logtau, loglambda, length, dimension)
        +log_terminal_power(packet(point, width, state), g['n']//width)
        -weight*(logx+loglambda)-g['threshold']*logy-logtau)/math.log(2)


@lru_cache(maxsize=8)
def regional_coefficients(logy, logtau, active, rows, width, state):
    """Reuse coefficients when only the total-weight mark changes."""
    tau = min(255., math.exp(logtau))
    local = np.asarray(local_matrices(math.exp(logy), tau, width, state)[:min(width, active)+1], dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    for w in range(len(factor)):
        factor[w] += math.log(math.comb(width, w))
    coefficients = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    steps = rows//width
    while steps:
        if steps & 1:
            coefficients = log_polynomial_product(coefficients, factor, active)
        steps >>= 1
        if steps:
            factor = log_polynomial_product(factor, factor, active)
    return coefficients


def sparse_bound(point, active, weight, width=8, state=32, length=128, dimension=96, k_log=20):
    """Fix W before extracting the existing sparse regional coefficients.

    The marker lambda is unrestricted positive; lambda^-W cancels it on
    the selected class. This is not the minimum-weight-only lambda bound.
    """
    g = geometry(k_log, length, dimension)
    if not feasible(active, weight, g['rows'], length, dimension):
        return -math.inf
    logtilt, logtau, loglambda = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    coefficients = regional_coefficients(logy, logtau, active, g['rows'], width, state)
    regional = regional_from_coefficients(coefficients, active, g['rows'], logtau+loglambda)
    return (log_choose(g['rows'], active)-((length-dimension)*active+1)*logtau
        -weight*loglambda-g['threshold']*logy+log_terminal(regional, length))/math.log(2)


def search(active, weight, mode='middle', width=8, state=32, length=128, dimension=96,
           k_log=20, iterations=70, start=None):
    g = geometry(k_log, length, dimension)
    if not feasible(active, weight, g['rows'], length, dimension):
        raise ValueError('empty active-count/weight class')
    if mode == 'middle':
        objective = lambda p: middle_bound(p, active, weight, width, state, length, dimension, k_log)
        starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
        bounds = [(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)]
    elif mode == 'sparse':
        objective = lambda p: sparse_bound(p, active, weight, width, state, length, dimension, k_log)
        starts = [(6.5, 2., 0.), (7., 3., 0.)]
        bounds = [(-2, 10), (0, math.log(255)), (-12, 12)]
    else:
        raise ValueError('unknown cover')
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(objective, p, method='L-BFGS-B', bounds=bounds,
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for p in starts]
    best = min(choices, key=lambda r: r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[length, dimension], width=width, state=state, k_log=k_log,
                rows=g['rows'], n=g['n'], active_rows=active, input_weight=weight,
                margin_bits=-float(best.fun), mode=mode, log_parameters=best.x.tolist(),
                optimizer_success=bool(best.success), iterations=int(best.nit),
                rigorous=False, whole_code=False)


def dominant_active(point, weight, length=128, dimension=96, k_log=20):
    """Binomial source mode at a supplied existing middle witness."""
    g = geometry(k_log, length, dimension)
    empty, occupied = row_logs(weight/g['n'], point[2], point[3], length, dimension)
    probability = math.exp(occupied-float(np.logaddexp(empty, occupied)))
    lower = (weight+length-1)//length
    upper = min(g['rows'], weight//(length-dimension+1))
    return min(upper, max(lower, math.floor((g['rows']+1)*probability)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('middle', 'sparse'), default='middle')
    parser.add_argument('--active', type=int, nargs='+', required=True)
    parser.add_argument('--weight', type=int, nargs='+', required=True)
    parser.add_argument('--width', type=int, default=8)
    parser.add_argument('--state', type=int, default=32)
    parser.add_argument('--outer', type=int, nargs=2, default=(128, 96))
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--iterations', type=int, default=70)
    args = parser.parse_args()
    for active in args.active:
        for weight in args.weight:
            print(json.dumps(search(active, weight, args.mode, args.width, args.state,
                                    *args.outer, args.k_log, args.iterations)), flush=True)
