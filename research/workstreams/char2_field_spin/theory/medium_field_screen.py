"""Bounded floating screens with genuinely uniform GF(2^32) diagonals.

The fixed GF256 symmetric Cauchy maps are embedded into GF(q); q is 256
or 2^32. Both input and state diagonals sample all of GF(q)^*. This is
not scalar extension of a sampled GF256 encoder. Positive packet covers
generalize the existing rank argument. Samples are never certificates.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from diagonal_refresh_rank import rank_potential_bound
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal
from symmetric_density_screen import joint_terms, minimum_return


FIELDS = (256, 2**32)


def validate(q, state):
    if q not in FIELDS or not 1 <= state <= 8:
        raise ValueError('require q in {256,2^32} and a state1..8 map restriction')


@lru_cache(maxsize=256, typed=True)
def root_factors(tau, q, state):
    validate(q, state)
    if not 1 <= tau <= q-1:
        raise ValueError('rank tilt outside the nonzero diagonal alphabet')
    return tuple(max(rank_potential_bound(r, z, tau, q, state)
                     for r in range(1, state+1)) for z in range(state))


def packet(weight, y, tau, q=2**32, state=4, joint=True):
    """Positive operator uniform over preceding common equation histories.

    Fractions remain exact. State8 uses its checked returned-weight law;
    smaller states use only the generic MDS lower bound. Joint and pure
    return bounds are complete covers, so taking their minimum is valid.
    """
    validate(q, state)
    if not 0 <= weight <= 16 or not 0 < y <= 1 or not 1 <= tau <= q-1:
        raise ValueError('invalid packet weight, output mark, or rank tilt')
    y = Fraction(y) if isinstance(y, int) else y
    tau = Fraction(tau) if isinstance(tau, int) else tau
    zero = y*0
    result = [[zero, zero], [zero, zero]]
    if not weight:
        result[0][0] = y**0
    else:
        result[0][1] = y**weight
        if weight > state:
            result[0][0] = y**weight/tau**state
        flats = sum(math.comb(state, k)/tau**k
                    for k in range(min(weight-1, state-1)+1))
        pure = flats*y**minimum_return(weight, state)/tau**state
        if joint:
            combined = sum(count*y**(16-i-z)*(1-y)**(i+z)/tau**rank
                           *sum(math.comb(state, k)/tau**k
                                for k in range(min(weight-1, state-1-z)+1))
                           for i, z, count, rank in joint_terms(weight, state))
            result[1][0] = min(pure, combined)
        else:
            result[1][0] = pure
    roots = root_factors(tau, q, state)
    # Factor the independent active-input pivot sum. This equals the
    # existing positive double sum and avoids repeating the root terms.
    result[1][1] = (y+(1-y)/tau)**weight*sum(
        math.comb(16-weight, z)*y**(16-weight-z)*(1-y)**z*roots[z]
        for z in range(min(state-1, 16-weight)+1))
    return result


def log_coefficients(active, y, tau, rows, q, state):
    if rows <= 0 or rows % 16 or not 0 <= active <= rows:
        raise ValueError('invalid aligned region geometry')
    values = np.asarray([packet(w, y, tau, q, state)
                         for w in range(min(active, 16)+1)], dtype=float)
    with np.errstate(divide='ignore'):
        factor = np.log(values)
    for w in range(len(factor)):
        factor[w] += math.log(math.comb(16, w))
    out = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    count = rows//16
    while count:
        if count & 1:
            out = log_polynomial_product(out, factor, active)
        count >>= 1
        if count:
            factor = log_polynomial_product(factor, factor, active)
    return out


def active_bound(point, active, q=2**32, state=4, length=64, dimension=48, k_log=20):
    validate(q, state)
    g = geometry(k_log, length, dimension)
    if not 0 < dimension < length <= 256 or not 1 <= active <= g['rows']:
        raise ValueError('invalid outer or active-row count')
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    tau = min(float(q-1), math.exp(logtau))
    coefficients = log_coefficients(active, math.exp(logy), tau, g['rows'], q, state)
    regional = regional_from_coefficients(coefficients, active, g['rows'], math.log(tau))
    return (log_choose(g['rows'], active)-((length-dimension)*active+1)*math.log(tau)
            -g['threshold']*logy+log_terminal(regional, length))/math.log(2)


def aggregate(logx, logy, logtau, q, state):
    tau = min(float(q-1), math.exp(logtau))
    values = np.asarray([packet(w, math.exp(logy), tau, q, state) for w in range(17)])
    with np.errstate(divide='ignore'):
        logs = np.log(values)
    for w in range(17):
        logs[w] += math.log(math.comb(16, w))+w*logx
    return logsumexp(logs, axis=0)


def density_bound(point, weight, q=2**32, state=4, length=64, dimension=48, k_log=20):
    validate(q, state)
    g = geometry(k_log, length, dimension)
    if not 0 < dimension < length <= 256 or not 0 < weight < g['n']:
        raise ValueError('invalid outer or input weight')
    alpha = weight/g['n']
    logx, logy, logtau, logtilt = map(float, point)
    logtau = min(math.log(q-1), logtau)
    row = float(logsumexp([length*math.log1p(-alpha),
                          -(length-dimension)*logtau+length*np.logaddexp(
                              math.log1p(-alpha), logtau+logtilt+math.log(alpha))]))
    logarithms = aggregate(logx, logy, logtau, q, state)
    shift = float(logarithms.max())
    steps = g['n']//16
    inner = steps*shift+log_terminal_power(np.exp(logarithms-shift), steps)
    return (length*math.log(g['rows']+1)+g['rows']*row+inner
            -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def search(kind, value, q=2**32, state=4, length=64, dimension=48,
           k_log=20, iterations=35, start=None):
    """Bounded multistart search; every returned number is floating only."""
    validate(q, state)
    g = geometry(k_log, length, dimension)
    if kind == 'active':
        value = int(value)
        starts = [(math.log(5*length), math.log(16)),
                  (math.log(5*length), math.log(q-1))]
        bounds = [(-2, 10), (0, math.log(q-1))]
        function = active_bound
    elif kind == 'density':
        if not 0 < value < 1:
            raise ValueError('density must lie in (0,1)')
        value = max(1, round(value*g['n']))
        starts = [(-2., -.2, math.log(16), 0.),
                  (0., -1., math.log(q-1), -1.)]
        bounds = [(-30, 30), (-25, 0), (0, math.log(q-1)), (-35, 25)]
        function = density_bound
    else:
        raise ValueError('unknown class kind')
    if q != 256:
        starts.append((math.log(5*length), math.log(255)) if kind == 'active'
                      else (0., -1., math.log(255), -1.))
    if start is not None:
        starts.insert(0, tuple(start))
    options = [minimize(function, point, args=(value, q, state, length, dimension, k_log),
                        method='L-BFGS-B', bounds=bounds,
                        options={'maxiter': iterations, 'ftol': 1e-10, 'maxls': 16})
               for point in starts]
    best = min(options, key=lambda option: option.fun if np.isfinite(option.fun) else math.inf)
    log_parameters = best.x.tolist()
    # A same-witness control only exists when the q32 tilt is admissible
    # for q256. Otherwise explicitly mark its absence, rather than clip it.
    logtau = log_parameters[1 if kind == 'active' else 2]
    control = (-function(best.x, value, 256, state, length, dimension, k_log)
               if logtau <= math.log(255) else None)
    return dict(q=q, outer=[length, dimension], state=state, width=16,
                k_log=k_log, rows=g['rows'], k=g['k'], n=g['n'],
                output_threshold=g['threshold'], kind=kind, class_value=value,
                input_density=value/g['n'] if kind == 'density' else None,
                margin_bits=-float(best.fun), log_parameters=log_parameters,
                tau=math.exp(logtau), q256_same_witness_margin_bits=control,
                optimizer_success=bool(best.success), iterations=int(best.nit),
                rigorous=False, whole_code=False,
                diagonal_law='independent uniform nonzero GF(q), input and state')


def compare(kind, value, state=4, length=64, dimension=48, k_log=20, iterations=35):
    """Same geometry, plus the GF256 witness explicitly seeded into GF32."""
    small = search(kind, value, 256, state, length, dimension, k_log, iterations)
    large = search(kind, value, 2**32, state, length, dimension, k_log, iterations,
                   small['log_parameters'])
    function = active_bound if kind == 'active' else density_bound
    at_small = -function(small['log_parameters'], small['class_value'],
                         2**32, state, length, dimension, k_log)
    return dict(q256=small, q32=large,
                q32_at_q256_witness_margin_bits=at_small,
                same_witness_gain_bits=at_small-small['margin_bits'],
                rigorous=False, whole_code=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, default=(64, 48))
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--q', type=int, choices=FIELDS, default=2**32)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--active', type=int, nargs='+')
    parser.add_argument('--density', type=float, nargs='+')
    parser.add_argument('--iterations', type=int, default=35)
    parser.add_argument('--compare', action='store_true', help='paired q256/q32 screens at identical geometry')
    args = parser.parse_args()
    if not args.active and not args.density:
        parser.error('select explicit bounded --active and/or --density samples')
    for kind, samples in (('active', args.active or ()), ('density', args.density or ())):
        previous = None
        for value in samples:
            if args.compare:
                result = compare(kind, value, args.state, *args.outer, args.k_log, args.iterations)
            else:
                result = search(kind, value, args.q, args.state, *args.outer,
                                args.k_log, args.iterations, previous)
                previous = result['log_parameters']
            print(json.dumps(result), flush=True)
