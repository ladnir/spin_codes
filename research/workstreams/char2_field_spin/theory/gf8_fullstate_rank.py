"""Rank-cap proof refinement for the non-self-orthogonal width8 candidates.

The construction is unchanged from gf8_fullstate_screen. States retain
the cap on incoming nonzero source-form rank. Sparse point evaluations
use rational marks and outward Arb arithmetic; no whole-code claim.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math

from flint import arb, arb_poly, ctx
import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_terminal_power
from field128_sparse_certificate import region
from gf8_fullstate_screen import RETURN_MINIMA, joint_terms, validate_state
from gf8_outer_rate_sparse import rank_terminal_power
from gf8_rank_history_screen import cone_preserved, modal_log_penalty, source_row_log
from live_kernel_rank import live_kernel_bound as exact_live_kernel_bound
from state8_sparse_certificate import ball, integer_bits, rational_marks
from state_rank_cap_transfer import rank_cap_weights


WIDTH = 8


def live_kernel_bound(tau, charge, dimension, ambient, q=256):
    """Rank potential with an actual full-support return witness.

    For U of dimension d, a live j-dimensional kernel contains at least
    q^(j-1)(q-m) full-support vectors. Each fixed full-support vector is
    killed by B D with probability at most (q-1)^(-rank B). The displayed
    bound maximizes d up to dimension, and never reuses a current diagonal.
    """
    if not (1 <= dimension <= ambient < q and 1 <= charge <= ambient and 1 <= tau <= q-1):
        raise ValueError('invalid live-kernel parameters')
    return min(tau**(-charge), exact_live_kernel_bound(dimension, ambient, charge, tau, q))


def return_bound(weight, y, tau, cap, state=8, sharp_return=False):
    minimum = RETURN_MINIMA[state][weight-1]
    def cost(charge, k):
        return (live_kernel_bound(tau, charge, cap-k+weight, state-k+weight)
                if sharp_return else tau**(-charge))
    pure = (sum(math.comb(state, k)*cost(state, k)/tau**k for k in range(min(weight-1, cap-1)+1))
            *y**minimum)
    joint = sum(count*y**(WIDTH-i-z)*(1-y)**(i+z)
                *sum(math.comb(state, k)*cost(charge, k)/tau**k
                     for k in range(min(weight-1, cap-1, state-1-z)+1))
                for i, z, count, charge in joint_terms(weight, state))
    return min(pure, joint)


def local_matrix(weight, y, tau, state=8, sharp_return=False):
    validate_state(state)
    if not 0 <= weight <= WIDTH or not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('invalid local marks')
    y, tau = (F(v) if isinstance(v, int) else v for v in (y, tau))
    zero = y*0+tau*0
    result = [[zero for _ in range(state+1)] for _ in range(state+1)]
    if weight == 0:
        result[0][0] = y**0
    else:
        result[0][min(state, weight)] = y**weight
        if weight > state:
            result[0][0] = y**weight/tau**state
    for cap in range(1, state+1):
        if weight:
            result[cap][0] = return_bound(weight, y, tau, cap, state, sharp_return)
        for i in range(weight+1):
            for z in range(min(state-1, WIDTH-weight)+1):
                factor = (math.comb(weight, i)*math.comb(WIDTH-weight, z)
                          *y**(WIDTH-i-z)*(1-y)**(i+z)/tau**i)
                for h, mass in enumerate(rank_cap_weights(cap, z, tau, 256, state), 1):
                    result[cap][min(state, h+weight-i)] += factor*mass
    return tuple(tuple(row) for row in result)


@lru_cache(maxsize=16)
def _mass_table(tau, state):
    return np.asarray([[rank_cap_weights(cap, z, tau, 256, state)
                        for z in range(min(state-1, WIDTH)+1)] for cap in range(1, state+1)], dtype=float)


def packet(logx, logy, logtau, state=8, sharp_return=False):
    x, y, tau = math.exp(logx), math.exp(logy), min(255., math.exp(logtau))
    result = np.zeros((state+1, state+1))
    result[0, 0] = 1.
    masses = _mass_table(tau, state)
    base = y+x*(1-y)/tau
    for z in range(min(state-1, WIDTH)+1):
        for d in range(WIDTH-z+1):
            coefficient = (math.comb(WIDTH, z)*math.comb(WIDTH-z, d)
                           *(x*y)**d*(1-y)**z*base**(WIDTH-z-d))
            for h in range(1, state+1):
                result[1:, min(state, h+d)] += coefficient*masses[:, z, h-1]
    for weight in range(1, WIDTH+1):
        multiplicity = math.comb(WIDTH, weight)*x**weight
        birth = multiplicity*y**weight
        result[0, min(state, weight)] += birth
        if weight > state:
            result[0, 0] += birth/tau**state
        for cap in range(1, state+1):
            result[cap, 0] += multiplicity*return_bound(weight, y, tau, cap, state, sharp_return)
    scale = float(result.max())
    if not 0 < scale < math.inf or not cone_preserved(result/scale, 2e-12):
        raise ArithmeticError('rank-cap continuation cone failed')
    return result


def middle_bound(point, weight, state=8, dimension=96, k_log=20, sharp_return=False):
    g = geometry(k_log, 128, dimension)
    logx, logy, logtau, logtilt = map(float, point)
    return (128*modal_log_penalty(g['rows'])
            +g['rows']*source_row_log(weight/g['n'], logtau, logtilt, 128, dimension, True)
            +log_terminal_power(packet(logx, logy, logtau, state, sharp_return), g['n']//WIDTH)
            -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def search_middle(density, state=8, dimension=96, k_log=20, iterations=60, start=None, sharp_return=False):
    g = geometry(k_log, 128, dimension)
    weight = max(1, round(density*g['n']))
    starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(lambda p: middle_bound(p, weight, state, dimension, k_log, sharp_return), p,
                       method='L-BFGS-B', bounds=[(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)],
                       options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for p in starts]
    best = min(choices, key=lambda r: r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[128, dimension], width=WIDTH, state=state, k_log=k_log,
                input_weight=weight, margin_bits=-float(best.fun),
                log_parameters=best.x.tolist(), optimizer_success=bool(best.success),
                iterations=int(best.nit), rigorous=False, whole_code=False, sharp_return=sharp_return,
                scope='one floating density class, retaining incoming rank cap')


def polynomial_power(local, packets, degree):
    size = len(local[0])
    result = [[arb_poly([int(i == j)]) for j in range(size)] for i in range(size)]
    factor = [[arb_poly([ball(local[w][i][j])*math.comb(WIDTH, w)
                         for w in range(min(WIDTH, degree)+1)])
               for j in range(size)] for i in range(size)]
    def product(left, right):
        return [[sum((left[i][k]*right[k][j] for k in range(size)), arb_poly([])).truncate(degree+1)
                 for j in range(size)] for i in range(size)]
    while packets:
        if packets & 1:
            result = product(result, factor)
        packets >>= 1
        if packets:
            factor = product(factor, factor)
    return tuple(value for row in result for value in row)


def sparse_point(active, point, state=8, dimension=96, k_log=20, precision=1024, sharp_return=False):
    g = geometry(k_log, 128, dimension)
    y, tau = rational_marks(active, point, g['n'])
    previous = ctx.prec
    try:
        ctx.prec = precision
        local = [local_matrix(w, y, tau, state, sharp_return) for w in range(WIDTH+1)]
        coefficients = polynomial_power(local, g['rows']//WIDTH, active)
        value = (arb(math.comb(g['rows'], active))*ball(tau)**(-((128-dimension)*active+1))
                 *ball(y)**(-g['threshold'])
                 *rank_terminal_power(region(coefficients, active, g['rows'], tau), 128))
        endpoint = value.upper()
        if not endpoint.is_finite() or not endpoint > 0:
            raise ArithmeticError('insufficient Arb precision')
        return dict(outer=[128, dimension], width=WIDTH, state=state, k_log=k_log,
                    active=active, y=str(y), tau=str(tau), precision_bits=precision,
                    margin_bits_lower=integer_bits(endpoint),
                    log2_upper_approx=float(endpoint.log()/arb(2).log()),
                    rigorous=True, whole_code=False, sharp_return=sharp_return,
                    scope='one active-row class, retaining incoming rank cap')
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=int, choices=(7, 8, 12, 16), default=8)
    parser.add_argument('--density', type=float, nargs='+', default=[.03, .10, .14])
    parser.add_argument('--dimension', type=int, default=96)
    parser.add_argument('--iterations', type=int, default=60)
    parser.add_argument('--sharp-return', action='store_true')
    args = parser.parse_args()
    previous = None
    for density in args.density:
        result = search_middle(density, args.state, args.dimension, iterations=args.iterations, start=previous,
                               sharp_return=args.sharp_return)
        previous = result['log_parameters']
        print(json.dumps(result), flush=True)
