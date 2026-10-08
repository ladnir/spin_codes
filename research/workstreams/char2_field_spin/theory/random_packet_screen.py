"""Exact local law and floating whole-class screens for a GF16 packet inner.

Each packet samples independent A,C rectangular Toeplitz matrices and a B
whose action sends each fixed nonzero input to a uniform nonzero vector.
The recurrence is Y=B X+A R, R'=R+C X. Only B must be invertible.
Fresh A and C give the exact zero/nonzero state transfer for every fixed
input and incoming state. No input diagonal or GL state refresh is used.

The outer has independent systematic rectangular Toeplitz rows. One
regional permutation routes each coordinate. This file reuses positive
source counting and coefficient Chernoff normalization; numeric searches
are FLOATING SCREENS, never certificates.
"""

import argparse
from fractions import Fraction
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

from coupled_chernoff import log_terminal_power
from field_rate_tradeoff import log_choose
from reduced_distance_certificate import log_polynomial_product


def geometry(k_log, length, dimension, width):
    if not 0 < dimension < length or width < 1:
        raise ValueError('invalid outer or packet width')
    rows = (((1 << k_log)+dimension*width-1)//(dimension*width))*width
    n = length*rows
    return dict(k=dimension*rows, n=n, rows=rows, threshold=(n+9)//10-1,
                length=length, dimension=dimension, width=width)


def local_matrices(z, width, state, q=16, input_weight=None):
    """Return exact Fraction matrices when z is Fraction, floats otherwise."""
    if not 0 <= z <= 1 or width < 1 or state < 1:
        raise ValueError('invalid packet parameters')
    one = Fraction(1) if isinstance(z, Fraction) else 1.
    zero = one-one
    output_sum = (one+(q-1)*z)**width
    uniform = output_sum/q**width
    if input_weight is None:
        positive_sum = output_sum-one if isinstance(z, Fraction) else math.expm1(width*math.log1p((q-1)*z))
        nonzero = positive_sum/(q**width-1)
    else:
        nonzero = z**input_weight
    returned = one/q**state
    empty = [[one, zero], [zero, uniform]]
    active = [[nonzero*returned, nonzero*(one-returned)],
              [uniform*returned, uniform*(one-returned)]]
    return empty, active


def local_by_weight(z, width, state, identity=False):
    empty, active = local_matrices(z, width, state)
    return np.asarray([empty]+[local_matrices(z, width, state, input_weight=w)[1] if identity else active
                               for w in range(1, width+1)], dtype=float)


def aggregate(x, z, width, state, identity=False, mark=1.):
    empty, active = (np.asarray(matrix) for matrix in local_matrices(z, width, state))
    inactive = (1+x)**width
    # All nonzero-value choices contribute (15*x)^w. Each zero input
    # position can belong to a selected source row, contributing 1+x.
    nonzero = math.expm1(width*math.log1p(15*x*mark/(1+x)))*inactive
    result = inactive*empty+nonzero*active
    if identity:
        born = math.expm1(width*math.log1p(15*x*z*mark/(1+x)))*inactive
        result[0] = [inactive+born/16**state, born*(1-16.**(-state))]
    return result


def log_regions(z, maximum, rows, width, state, identity=False):
    with np.errstate(divide='ignore'):
        factors = np.log(local_by_weight(z, width, state, identity)[:maximum+1])
    for weight in range(len(factors)):
        factors[weight] += math.log(math.comb(width, weight))
    regions = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows//width
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, factors, maximum)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, maximum)
    for weight in range(maximum+1):
        regions[weight] -= log_choose(rows, weight)
    return regions


def positive_source_moment(regions, active, length, dimension):
    """Exact nonempty-systematic row count, with a positive Pascal transform.

    H_seen(new)=sum_old C(seen,old)*15**(old+new)*R_(old+new).
    Therefore H_(seen+1)(new)=H_seen(new)+H_seen(new+1).
    This builds the source transition in quadratic, not cubic, work.
    """
    weighted = regions+np.arange(active+1)[:, None, None]*math.log(15)
    matrix = logsumexp(weighted+np.array([log_choose(active, h) for h in range(active+1)])[:, None, None], axis=0)
    transition = np.full((active+1, 2, active+1, 2), -math.inf)
    pascal = weighted
    for seen in range(active+1):
        counts = np.array([log_choose(active-seen, new) for new in range(active-seen+1)])
        transition[seen, :, seen:, :] = np.transpose(pascal+counts[:, None, None], (1, 0, 2))
        pascal = np.logaddexp(pascal[:-1], pascal[1:])
    current = np.full((active+1, 2), -math.inf)
    current[0, 0] = 0.
    for _ in range(dimension):
        current = logsumexp(current[:, :, None, None]+transition, axis=(0, 1))
    vector = current[active]
    for _ in range(length-dimension):
        vector = logsumexp(vector[:, None]+matrix, axis=0)
    return float(logsumexp(vector))-(length-dimension)*active*math.log(16)


def sparse_bound(active, z, k_log=20, length=48, dimension=32, width=4, state=1, identity=False):
    g = geometry(k_log, length, dimension, width)
    if not 1 <= active <= g['rows']:
        raise ValueError('invalid active row count')
    regions = log_regions(z, active, g['rows'], width, state, identity)
    moment = positive_source_moment(regions, active, length, dimension)
    value = log_choose(g['rows'], active)-math.log(15)-g['threshold']*math.log(z)+moment
    return -value/math.log(2)


def sparse_search(active, **parameters):
    g = geometry(parameters.get('k_log', 20), parameters.get('length', 48), parameters.get('dimension', 32), parameters.get('width', 4))
    def objective(tilt):
        z = math.exp(-math.exp(tilt)*active/g['n'])
        return -sparse_bound(active, z, **parameters)
    # The optimum can migrate from a near-one tilt to a conventional
    # weight Chernoff tilt. Keep both intervals in this exploratory search.
    candidates = [minimize_scalar(objective, bounds=bounds, method='bounded') for bounds in ((-2, 8), (6, 14))]
    best = min(candidates, key=lambda result: result.fun)
    return dict(**g, state=parameters.get('state', 1), active=active, z=math.exp(-math.exp(best.x)*active/g['n']),
                margin_bits=-float(best.fun), rigorous=False, method='positive-exact-systematic-source',
                input_map='identity' if parameters.get('identity', False) else 'uniform-nonzero-transitive')


def middle_bound(active, x, z, k_log=20, length=48, dimension=32, width=4, state=1, identity=False):
    g = geometry(k_log, length, dimension, width)
    value = (-(length-1)*log_choose(g['rows'], active)-(length-dimension)*active*math.log(16)
             -math.log(15)-g['threshold']*math.log(z)-length*active*math.log(x)
             +log_terminal_power(aggregate(x, z, width, state, identity), g['rows']*length//width))
    return -value/math.log(2)


def middle_search(active, **parameters):
    g = geometry(parameters.get('k_log', 20), parameters.get('length', 48), parameters.get('dimension', 32), parameters.get('width', 4))
    def objective(point):
        logx, logz = point
        if not -35 < logx < 15 or not -25 < logz < 0:
            return 1e100
        return -middle_bound(active, math.exp(logx), math.exp(logz), **parameters)
    initial_x = active/(16*max(1, g['rows']-active))
    candidates = [minimize(objective, [math.log(initial_x), math.log(z)], method='Nelder-Mead',
                           options={'maxiter': 400, 'xatol': 1e-9, 'fatol': 1e-7}) for z in (.05, .5, .9, .99, .9999)]
    best = min(candidates, key=lambda result: result.fun)
    return dict(**g, state=parameters.get('state', 1), active=active,
                x=math.exp(best.x[0]), z=math.exp(best.x[1]), margin_bits=-float(best.fun),
                rigorous=False, method='coupled-coefficient-chernoff',
                input_map='identity' if parameters.get('identity', False) else 'uniform-nonzero-transitive')


def log_matrix_power(matrix, exponent):
    scale = float(np.max(matrix))
    matrix = matrix/scale
    logarithm = math.log(scale)
    result = np.eye(2)
    accumulated = 0.
    while exponent:
        if exponent & 1:
            result = result @ matrix
            norm = float(np.max(result))
            result /= norm
            accumulated += logarithm+math.log(norm)
        exponent >>= 1
        if exponent:
            matrix = matrix @ matrix
            norm = float(np.max(matrix))
            matrix /= norm
            logarithm = 2*logarithm+math.log(norm)
    with np.errstate(divide='ignore'):
        return np.log(result)+accumulated


def regional_entries(active, z, rows, width, state, identity=False, mark=1.):
    entries, witnesses = np.zeros((2, 2)), np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            def objective(logx):
                logs = log_matrix_power(aggregate(math.exp(logx), z, width, state, identity, mark), rows//width)
                return logs[i, j]-active*logx-log_choose(rows, active)
            optimum = minimize_scalar(objective, bounds=(-30, 15), method='bounded', options={'xatol': 1e-8})
            entries[i, j], witnesses[i, j] = optimum.fun, math.exp(optimum.x)
    return entries, witnesses


def regional_bound(active, z, k_log=20, length=48, dimension=32, width=4, state=1, identity=False,
                   details=False, mark=1.):
    """Separate regional entry bounds, optionally marking systematic weight.

    Every genuinely active row has a nonzero systematic input, hence total
    systematic weight is at least active. For mark>=1, divide the marked
    all-source moment by mark**active. This is a relaxation, not the exact
    per-row nonempty condition used by sparse_bound.
    """
    if mark < 1:
        raise ValueError('systematic lower-weight marker must be at least one')
    g = geometry(k_log, length, dimension, width)
    entries, witnesses = regional_entries(active, z, g['rows'], width, state, identity)
    systematic, systematic_x = (entries, witnesses) if mark == 1 else regional_entries(active, z, g['rows'], width, state, identity, mark)
    vector = np.array([0., -math.inf])
    for region in range(length):
        vector = logsumexp(vector[:, None]+(systematic if region < dimension else entries), axis=0)
    value = (log_choose(g['rows'], active)-math.log(15)-(length-dimension)*active*math.log(16)
             -g['threshold']*math.log(z)+float(logsumexp(vector))-active*math.log(mark))
    margin = -value/math.log(2)
    return dict(margin_bits=margin, x=witnesses.tolist(), systematic_x=systematic_x.tolist(), mark=mark) if details else margin


def regional_search(active, **parameters):
    g = geometry(parameters.get('k_log', 20), parameters.get('length', 48), parameters.get('dimension', 32), parameters.get('width', 4))
    def objective(tilt):
        z = math.exp(-math.exp(tilt)*active/g['n'])
        return -regional_bound(active, z, **parameters)
    limit = math.log(20*g['n']/active)
    optimum = minimize_scalar(objective, bounds=(-2, min(14, limit)), method='bounded', options={'xatol': 1e-6})
    z = math.exp(-math.exp(optimum.x)*active/g['n'])
    return dict(**g, active=active, state=parameters.get('state', 1), z=z,
                **regional_bound(active, z, details=True, **parameters), rigorous=False,
                method='entrywise-regional-coefficient-chernoff',
                input_map='identity' if parameters.get('identity', False) else 'uniform-nonzero-transitive')


def marked_search(active, **parameters):
    g = geometry(parameters.get('k_log', 20), parameters.get('length', 48), parameters.get('dimension', 32), parameters.get('width', 4))
    def objective(point):
        tilt, logmark = point
        if not -2 < tilt < min(14, math.log(20*g['n']/active)) or not 0 <= logmark < 15:
            return 1e100
        z = math.exp(-math.exp(tilt)*active/g['n'])
        return -regional_bound(active, z, mark=math.exp(logmark), **parameters)
    initial = regional_search(active, **parameters)
    tilt = math.log(-math.log(initial['z'])*g['n']/active)
    optimum = minimize(objective, [tilt, 1.], method='Nelder-Mead', options={'maxiter': 250, 'xatol': 1e-6, 'fatol': 1e-5})
    if -optimum.fun < initial['margin_bits']:
        return initial
    z = math.exp(-math.exp(optimum.x[0])*active/g['n'])
    return dict(**g, active=active, state=parameters.get('state', 1), z=z,
                **regional_bound(active, z, mark=math.exp(optimum.x[1]), details=True, **parameters),
                rigorous=False, method='regional-entrywise-with-total-systematic-weight-marker',
                input_map='identity' if parameters.get('identity', False) else 'uniform-nonzero-transitive')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', nargs=2, type=int, default=(48, 32))
    parser.add_argument('--width', type=int, default=4)
    parser.add_argument('--state', type=int, default=1)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--active', nargs='+', type=int, default=(1, 2, 4, 8, 16))
    parser.add_argument('--middle', action='store_true')
    parser.add_argument('--identity', action='store_true', help='fixed B=I, preserving the exact input weight at births')
    parser.add_argument('--regional', action='store_true', help='optimize each regional boundary-state entry separately')
    parser.add_argument('--marked', action='store_true', help='also enforce the total systematic weight lower bound')
    args = parser.parse_args()
    method = marked_search if args.marked else regional_search if args.regional else middle_search if args.middle else sparse_search
    for active in args.active:
        print(json.dumps(method(active, k_log=args.k_log, length=args.outer[0], dimension=args.outer[1],
                                width=args.width, state=args.state, identity=args.identity)), flush=True)
