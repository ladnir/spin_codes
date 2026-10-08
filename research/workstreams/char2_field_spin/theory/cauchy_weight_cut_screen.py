"""Floating regional rank screen retaining the feasible total input weight.

This uses the same two-state upper transfer as cauchy_regional_rank_screen.
It keeps the polynomial through all outer regions before discarding weights
below (outer parity+1)*active or above the already certified dense cutoff.
It is a diagnostic of lost outer-support constraints, not a certificate.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_choose
from cauchy_joint_return import rank_return_moment
from reduced_distance_certificate import log_polynomial_product
from sparse_multirow_zero_set import local_rank_transfer


def log_power(polynomial, exponent, limit):
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    while exponent:
        if exponent & 1:
            result = log_polynomial_product(result, polynomial, limit)
        exponent >>= 1
        if exponent:
            polynomial = log_polynomial_product(polynomial, polynomial, limit)
    return result


def region_coefficients(active, y, tau, rows):
    if rows % 16:
        raise ValueError('rows must be a multiple of sixteen')
    local = np.asarray([local_rank_transfer(w, y, tau)
                        for w in range(min(16, active)+1)], dtype=float)
    for w in range(1, min(3, active)+1):
        local[w, 1, 0] = min(local[w, 1, 0],
            rank_return_moment(w, y, tau, average_inputs=True))
    with np.errstate(divide='ignore'):
        coefficients = np.log(local)
    for w in range(len(local)):
        coefficients[w] += math.log(math.comb(16, w))
    coefficients = log_power(coefficients, rows//16, active)
    for w in range(active+1):
        coefficients[w] += (log_choose(active, w)-log_choose(rows, w)
                            +w*math.log(tau))
    return coefficients


def evaluate(active, y, tau, k_log=20, length=32, dimension=22):
    g = geometry(k_log, length, dimension)
    minimum = (length-dimension+1)*active
    # W < ceil(.11N), the complement of the directed intermediate class.
    maximum = min(length*active, (11*g['n']+99)//100-1)
    if minimum > maximum:
        return dict(margin_bits=math.inf, removed_bits=0., mean_weight=None,
                    rigorous=False)
    regional = region_coefficients(active, y, tau, g['rows'])
    coefficients = log_power(regional, length, length*active)
    logs = logsumexp(coefficients[:, 0, :], axis=1)
    full = float(logsumexp(logs))
    cut = float(logsumexp(logs[minimum:maximum+1]))
    prefactor = (log_choose(g['rows'], active)
                 -((length-dimension)*active+1)*math.log(tau)
                 -g['threshold']*math.log(y))
    return dict(active=active, k_log=k_log, minimum_weight=minimum,
                maximum_weight=maximum, margin_bits=-(prefactor+cut)/math.log(2),
                uncut_margin_bits=-(prefactor+full)/math.log(2),
                removed_bits=(full-cut)/math.log(2),
                mean_weight=float(np.dot(np.arange(len(logs)), np.exp(logs-full))),
                peak_weight=int(np.argmax(logs)), rigorous=False)


def evaluate_point(active, logtilt, logtau, k_log=20):
    g = geometry(k_log)
    y = math.exp(-math.exp(logtilt)*active/g['n'])
    return evaluate(active, y, math.exp(logtau), k_log)


def search(active, point, k_log=20, iterations=30):
    objective = lambda p: -evaluate_point(active, *p, k_log)['margin_bits']
    result = minimize(objective, point, method='Nelder-Mead',
                      bounds=[(-2, 13), (0, math.log(255))],
                      options={'maxiter': iterations, 'xatol': 1e-5, 'fatol': 1e-5})
    out = evaluate_point(active, *result.x, k_log)
    out.update(logtilt=float(result.x[0]), logtau=float(result.x[1]))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active', type=int, required=True)
    parser.add_argument('--logtilt', type=float, required=True)
    parser.add_argument('--logtau', type=float, required=True)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--iterations', type=int, default=0)
    args = parser.parse_args()
    if args.iterations:
        out = search(args.active, [args.logtilt, args.logtau], args.k_log, args.iterations)
    else:
        out = evaluate_point(args.active, args.logtilt, args.logtau, args.k_log)
    print(json.dumps(out), flush=True)
