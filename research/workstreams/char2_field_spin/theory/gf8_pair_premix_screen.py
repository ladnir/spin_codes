"""Conservative positive cover for deterministic 2x2 MDS packet premixing.

Order: regional route, four fixed MDS pairs, fresh nonzero input diagonal,
unchanged repaired cubic 8/4 inner. The original raw source support and
source rank remain unchanged. Mixed zero flats are imposed in the common
source quotient before any current input/state diagonal is exposed.
This module is an exploratory bound, not a whole-code certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import log_choose, log_terminal_power
from gf8_high_rate_outer_screen import RETURN_MINIMA, profile
from gf8_packet_width_screen import local_matrix
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


def pair_patterns(pairs=4):
    """(singletons,doubles,raw_weight,placement_multiplicity)."""
    for singles in range(pairs+1):
        for doubles in range(pairs-singles+1):
            count = math.comb(pairs, doubles)*math.comb(pairs-doubles, singles)*2**singles
            yield singles, doubles, singles+2*doubles, count


def mixed_cover(singles, doubles, y, tau, locals_=None):
    """Uniform fixed-support two-state cover for these raw pair counts.

    Rank-k mixed-zero flats have <=C(d,k)2^k canonical bases. A feasible
    flat has k<=z<=d vanished coordinates, and mixed weight 2a+2d-z.
    Rank-zero flats may already contain cancellations: there is no extra
    1/(q-1) probability for the deterministic mixer.
    """
    if not 0 <= singles <= 4 or not 0 <= doubles <= 4-singles:
        raise ValueError('invalid pair counts')
    if locals_ is None:
        locals_ = [local_matrix(w, y, tau, 8, 4, RETURN_MINIMA) for w in range(9)]
    zero = y*0+tau*0
    result = [[zero, zero], [zero, zero]]
    for rank in range(doubles+1):
        coefficient = math.comb(doubles, rank)*2**rank/tau**rank
        weights = range(2*singles+doubles, 2*singles+2*doubles-rank+1)
        for row in range(2):
            for column in range(2):
                result[row][column] += coefficient*max(locals_[w][row][column] for w in weights)
    return result


@lru_cache(maxsize=16)
def weighted_matrices(y, tau):
    """Sum over all raw supports of each weight, not a pointwise weight bound."""
    local = [local_matrix(w, y, tau, 8, 4, RETURN_MINIMA) for w in range(9)]
    zero = y*0+tau*0
    result = [[[zero, zero], [zero, zero]] for _ in range(9)]
    for singles, doubles, weight, count in pair_patterns():
        cover = mixed_cover(singles, doubles, y, tau, local)
        for row in range(2):
            for column in range(2):
                result[weight][row][column] += count*cover[row][column]
    return result


def packet(logx, logy, logtau):
    x, y, tau = math.exp(logx), math.exp(logy), min(255., math.exp(logtau))
    local = np.asarray(weighted_matrices(y, tau), dtype=float)
    result = np.einsum('w,wij->ij', np.asarray([x**w for w in range(9)]), local)
    with np.errstate(divide='ignore'):
        return np.log(result)


def middle_bound(point, weight, length=128, dimension=96, k_log=20, components=False):
    g = profile(length, dimension, k_log)
    logx, logy, logtau, logtilt = map(float, point)
    logs = packet(logx, logy, logtau)
    shift = float(logs.max())
    steps = g['n']//8
    terms = dict(route=length*modal_log_penalty(g['rows']),
                 source=g['rows']*source_row_log(weight/g['n'], logtau, logtilt, length, dimension, True),
                 inner=steps*shift+log_terminal_power(np.exp(logs-shift), steps),
                 selection=-weight*(logx+logtilt), output=-g['threshold']*logy,
                 rank_threshold=-logtau)
    value = sum(terms.values())/math.log(2)
    return (value, {key: val/math.log(2) for key, val in terms.items()}) if components else value


def regional_coefficients(active, y, tau, rows):
    if rows % 8 or not 1 <= active <= rows:
        raise ValueError('invalid packet-aligned sparse geometry')
    local = np.asarray(weighted_matrices(y, tau)[:min(8, active)+1], dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(local)
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    packets = rows//8
    while packets:
        if packets & 1:
            result = log_polynomial_product(result, factor, active)
        packets >>= 1
        if packets:
            factor = log_polynomial_product(factor, factor, active)
    return result


def sparse_bound(active, point, length=128, dimension=96, k_log=20, components=False):
    g = profile(length, dimension, k_log)
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(255., math.exp(logtau))
    coefficients = regional_coefficients(active, math.exp(logy), tau, g['rows'])
    regional = regional_from_coefficients(coefficients, active, g['rows'], math.log(tau))
    terms = dict(source=log_choose(g['rows'], active)-((length-dimension)*active+1)*math.log(tau),
                 output=-g['threshold']*logy, regional=log_terminal(regional, length))
    value = sum(terms.values())/math.log(2)
    return (value, {key: val/math.log(2) for key, val in terms.items()}) if components else value


def search(value, mode='middle', length=128, dimension=96, k_log=20, iterations=60, start=None):
    g = profile(length, dimension, k_log)
    if mode == 'middle':
        target = max(1, round(value*g['n']))
        objective = lambda p, components=False: middle_bound(p, target, length, dimension, k_log, components)
        starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
        bounds = [(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)]
    elif mode == 'sparse':
        target = int(value)
        objective = lambda p, components=False: sparse_bound(target, p, length, dimension, k_log, components)
        starts = [(math.log(640), math.log(255)), (math.log(640), math.log(8))]
        bounds = [(-2, 10), (0, math.log(255))]
    else:
        raise ValueError('unknown screen mode')
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(objective, point, method='L-BFGS-B', bounds=bounds,
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for point in starts]
    best = min(choices, key=lambda r: r.fun if np.isfinite(r.fun) else math.inf)
    result, terms = objective(best.x, True)
    g.update(mode=mode, target=target, margin_bits=-result, log_parameters=best.x.tolist(),
             log2_components=terms, optimizer_success=bool(best.success), iterations=int(best.nit),
             construction='deterministic-pair-premix-before-input-diagonals-cubic8/4',
             scope='one weight class with mixed-zero-flat cover, floating arithmetic')
    return g


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, default=(128, 96))
    parser.add_argument('--mode', choices=('middle', 'sparse'), default='middle')
    parser.add_argument('--density', type=float, nargs='+', default=[.10, .14])
    parser.add_argument('--active', type=int, nargs='+', default=[64, 128])
    parser.add_argument('--iterations', type=int, default=60)
    parser.add_argument('--k-log', type=int, default=20)
    args = parser.parse_args()
    previous = None
    for value in args.density if args.mode == 'middle' else args.active:
        result = search(value, args.mode, *args.outer, args.k_log, args.iterations, previous)
        previous = result['log_parameters']
        print(json.dumps(result), flush=True)
