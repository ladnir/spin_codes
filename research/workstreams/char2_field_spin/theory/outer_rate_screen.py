"""One-route GF16 screens with a lower-dimensional GRS[16,d] outer.

The inner remains the monomial width-16 construction from gf16_sizes.
Searches use floating point and never certify a complete construction.
Rational witnesses can be replayed for individual sparse or middle classes.
No implementation timing is performed by this module.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math

import numpy as np

from gf16_sizes import (local_model, interval_local_envelope, upper_regions,
                        interval_fraction)
from gf16_certificate import HIGH, LOW, Interval
from q1_batched_screen import fixed_support_multiplicity
from subfield_certificate import integer_margin


def geometry(requested_k=1 << 20, dimension=7):
    """Round up to complete outer rows and sixteen-symbol regional steps."""
    if requested_k < 1 or not 1 <= dimension <= 8:
        raise ValueError('require a positive requested size and dimension 1..8')
    rows = 16*((requested_k+16*dimension-1)//(16*dimension))
    k, n = dimension*rows, 16*rows
    return dict(requested_k=requested_k, k=k, padding=k-requested_k,
                rows=rows, n=n, dimension=dimension, threshold=(n+9)//10-1)


def source_eta(dimension, rho):
    """Exact positive envelope for each possible nonzero row support."""
    rho = Fraction(rho)
    if rho <= 0 or not 1 <= dimension <= 8:
        raise ValueError('invalid dimension or source tilt')
    return max(Fraction(fixed_support_multiplicity(16, dimension, w, 16),
                        16**(w-(16-dimension)))/rho**w
               for w in range(17-dimension, 17))


def log_regions(z, active, rows, state):
    from reduced_distance_certificate import log_polynomial_product
    with np.errstate(divide='ignore'):
        local = np.log(local_model(z, state)[:active+1])
    for w in range(len(local)):
        local[w] += math.log(math.comb(16, w))
    regions = np.array([[[0, -math.inf], [-math.inf, 0]]])
    power = rows//16
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, local, active)
        power >>= 1
        if power:
            local = log_polynomial_product(local, local, active)
    for h in range(active+1):
        regions[h] -= math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
    return regions


def exact_source(regions, active, dimension, directed=False):
    """Enumerate every outer support mask through four active rows."""
    if not 1 <= active <= 4:
        raise ValueError('exact source computation is limited to four rows')
    dtype = object if directed else float
    scalar = Decimal if directed else float
    matrices = np.asarray(regions, dtype=dtype).reshape((-1, 2, 2))
    masks = [tuple((mask >> j) & 1 for j in range(active)) for mask in range(1 << active)]
    with localcontext(HIGH):
        current = np.full((1,)*active+(2,), scalar(0), dtype=dtype)
        current[(0,)*active+(0,)] = scalar(1)
        for region in range(16):
            following = np.full((region+2,)*active+(2,), scalar(0), dtype=dtype)
            for mask in masks:
                target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
                following[target] += current @ matrices[sum(mask)]
            current = following
        counts = [scalar(fixed_support_multiplicity(16, dimension, w, 16))
                  if w >= 17-dimension else scalar(0) for w in range(17)]
        for axis in range(active):
            shape = [1]*(active+1)
            shape[axis] = 17
            current *= np.array(counts, dtype=dtype).reshape(shape)
        return np.sum(current)


def screen_bound(active, z, requested_k=1 << 20, dimension=7, state=4, rho=.871):
    from scipy.special import logsumexp
    geo = geometry(requested_k, dimension)
    rows, threshold = geo['rows'], geo['threshold']
    if not 1 <= active < rows or not 0 < z < 1:
        raise ValueError('invalid active-row class or output tilt')
    regions = log_regions(z, active, rows, state)
    if active <= 4:
        scale = float(np.max(regions))
        moment = exact_source(np.exp(regions-scale), active, dimension)
        moment_log = math.log(float(moment))+16*scale
    else:
        source = np.array([math.log(math.comb(active, h))+h*math.log(16*rho)
                           for h in range(active+1)])
        matrix = logsumexp(regions+source[:, None, None], axis=0)
        vector = np.array([0, -math.inf])
        for _ in range(16):
            vector = logsumexp(vector[:, None]+matrix, axis=0)
        moment_log = (float(logsumexp(vector))+active*math.log(float(source_eta(dimension, rho)))
                      -(16-dimension)*active*math.log(16))
    choose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    return -(choose-math.log(15)-threshold*math.log(z)+moment_log)/math.log(2)


def search(active, requested_k=1 << 20, dimension=7, state=4, rho=.871):
    from scipy.optimize import minimize_scalar
    geo = geometry(requested_k, dimension)
    def objective(logtilt):
        z = math.exp(-math.exp(logtilt)*active/geo['n'])
        return -screen_bound(active, z, requested_k, dimension, state, rho)
    optimum = minimize_scalar(objective, bounds=(1, 9), method='bounded', options=dict(xatol=1e-6))
    return dict(**geo, state=state, active_rows=active,
                z=math.exp(-math.exp(optimum.x)*active/geo['n']),
                rho=rho if active > 4 else None, margin_bits=-float(optimum.fun), rigorous=False)


def search_source(active, requested_k=1 << 20, dimension=7, state=4):
    from scipy.optimize import minimize
    if active <= 4:
        return search(active, requested_k, dimension, state)
    geo = geometry(requested_k, dimension)
    initial = search(active, requested_k, dimension, state)
    initial_logtilt = math.log(-math.log(initial['z'])*geo['n']/active)
    def objective(point):
        lt, lr = point
        if not 1 <= lt <= 9 or not -2 <= lr <= 2:
            return 1e100
        z = math.exp(-math.exp(lt)*active/geo['n'])
        return -screen_bound(active, z, requested_k, dimension, state, math.exp(lr))
    optimum = minimize(objective, [initial_logtilt, math.log(.871)], method='Nelder-Mead',
                       options=dict(maxiter=160, xatol=1e-6, fatol=1e-6))
    return dict(**geo, state=state, active_rows=active,
                z=math.exp(-math.exp(optimum.x[0])*active/geo['n']),
                rho=math.exp(optimum.x[1]), margin_bits=-float(optimum.fun), rigorous=False)


def class_upper(active, z, requested_k=1 << 20, dimension=7, state=4, rho=Fraction(871, 1000)):
    """Directed bound for one active-row class, with explicit natural size."""
    geo = geometry(requested_k, dimension)
    zi, rhoi = interval_fraction(z), interval_fraction(rho)
    regions = upper_regions(zi, active, geo['rows'], state)
    if active <= 4:
        moment = exact_source(regions, active, dimension, directed=True)
    else:
        matrix = [Decimal(0)]*4
        for h, regional in enumerate(regions):
            weight = HIGH.multiply(Decimal(math.comb(active, h)*16**h), (rhoi**h).hi)
            matrix = [HIGH.add(a, HIGH.multiply(weight, b)) for a, b in zip(matrix, regional)]
        vector = [Decimal(1), Decimal(0)]
        for _ in range(16):
            vector = [HIGH.add(HIGH.multiply(vector[0], matrix[0]), HIGH.multiply(vector[1], matrix[2])),
                      HIGH.add(HIGH.multiply(vector[0], matrix[1]), HIGH.multiply(vector[1], matrix[3]))]
        moment = HIGH.multiply(HIGH.add(*vector), (interval_fraction(source_eta(dimension, rho))**active).hi)
        moment = HIGH.divide(moment, Decimal(16**((16-dimension)*active)))
    top = HIGH.multiply(moment, Decimal(math.comb(geo['rows'], active)))
    bottom = LOW.multiply(Decimal(15), (zi**geo['threshold']).lo)
    bound = HIGH.divide(top, bottom)
    return dict(**geo, state=state, active_rows=active, z=str(Fraction(z)),
                rho=str(Fraction(rho)) if active > 4 else None,
                failure_upper=str(bound), integer_bits=integer_margin(bound), rigorous=True, complete=False)


def middle_screen(active, x, z, requested_k=1 << 20, dimension=7, state=4, rho=Fraction(871, 1000)):
    from coupled_chernoff import log_terminal_power
    geo = geometry(requested_k, dimension)
    rows, threshold = geo['rows'], geo['threshold']
    local = local_model(z, state)
    matrix = sum(math.comb(16, w)*(16*float(rho)*x)**w*(1+x)**(16-w)*local[w] for w in range(17))
    choose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    value = (-15*choose-(16-dimension)*active*math.log(16)
             +active*math.log(float(source_eta(dimension, rho)))-math.log(15)
             -threshold*math.log(z)-16*active*math.log(x)+log_terminal_power(matrix, rows))
    return -value/math.log(2)


def middle_search(active, requested_k=1 << 20, dimension=7, state=4, rho=Fraction(871, 1000)):
    from scipy.optimize import minimize
    geo = geometry(requested_k, dimension)
    def objective(point):
        logx, logz = point
        if not -35 < logx < 10 or not -20 < logz < 0:
            return 1e100
        return -middle_screen(active, math.exp(logx), math.exp(logz), requested_k, dimension, state, rho)
    best = None
    for initial_z in [.5, .9, .99, .999]:
        initial_x = active/(max(1, geo['rows']-active)*16)
        optimum = minimize(objective, [math.log(initial_x), math.log(initial_z)], method='Nelder-Mead',
                           options=dict(maxiter=300, xatol=1e-9, fatol=1e-7))
        candidate = dict(**geo, state=state, active_rows=active, x=math.exp(optimum.x[0]),
                         z=math.exp(optimum.x[1]), rho=str(rho), margin_bits=-float(optimum.fun), rigorous=False)
        if best is None or candidate['margin_bits'] > best['margin_bits']:
            best = candidate
    return best


def middle_upper(active, x, z, requested_k=1 << 20, dimension=7, state=4, rho=Fraction(871, 1000)):
    from coupled_chernoff import interval_bound
    geo = geometry(requested_k, dimension)
    rhoi = interval_fraction(rho)
    def local_provider(zi):
        return [[[entry*rhoi**w for entry in row] for row in matrix]
                for w, matrix in enumerate(interval_local_envelope(zi, state))]
    return (interval_bound(active, x, z, local_provider, rows=geo['rows'], dimension=dimension)
            *interval_fraction(source_eta(dimension, rho))**active)


def sparse_summary(maximum=128, requested_k=1 << 20, dimension=7, state=4,
                   rho=.96, optimize_through=16):
    entries = [search_source(a, requested_k, dimension, state)
               if 5 <= a <= optimize_through else search(a, requested_k, dimension, state, rho)
               for a in range(1, maximum+1)]
    total = sum(2**(-entry['margin_bits']) for entry in entries)
    return dict(**geometry(requested_k, dimension), state=state,
                active_rows=[1, maximum], fixed_rho=rho, optimize_through=optimize_through,
                sum_margin_bits=-math.log2(total), worst=sorted(entries, key=lambda r: r['margin_bits'])[:5],
                rigorous=False, complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--requested-k', type=int, default=1 << 20)
    parser.add_argument('--dimension', type=int, choices=[6, 7, 8], default=7)
    parser.add_argument('--state', type=int, choices=range(1, 9), default=4)
    parser.add_argument('--active', nargs='+', type=int, default=[1, 4, 8, 16, 32, 64, 128])
    parser.add_argument('--middle', action='store_true')
    parser.add_argument('--optimize-source', action='store_true')
    parser.add_argument('--source-rho', type=Fraction, default=Fraction(871, 1000))
    parser.add_argument('--summary', action='store_true')
    parser.add_argument('--optimize-through', type=int, default=0)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    if args.summary:
        print(json.dumps(sparse_summary(max(args.active), args.requested_k, args.dimension, args.state,
                                        float(args.source_rho), args.optimize_through)), flush=True)
        raise SystemExit(0)
    for active in args.active:
        if args.optimize_source and not args.middle:
            result = search_source(active, args.requested_k, args.dimension, args.state)
        else:
            method = middle_search if args.middle else search
            result = method(active, args.requested_k, args.dimension, args.state, float(args.source_rho))
        if args.replay:
            z = Fraction(round(result['z']*2**40), 2**40)
            if args.middle:
                x = Fraction(round(result['x']*10**12), 10**12)
                bound = middle_upper(active, x, z, args.requested_k, args.dimension, args.state, args.source_rho)
                result.update(x=str(x), z=str(z), failure_upper=str(bound.hi),
                              integer_bits=integer_margin(bound.hi), rigorous=True, complete=False)
            else:
                rho = Fraction(round(result['rho']*10**6), 10**6) if result['rho'] else Fraction(871, 1000)
                result = class_upper(active, z, args.requested_k, args.dimension, args.state, rho)
        print(json.dumps(result), flush=True)
