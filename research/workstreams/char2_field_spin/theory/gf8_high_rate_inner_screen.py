"""Proof-only high-rate screens for the fixed repaired cubic 8/4 inner.

The five-state transfer retains the incoming state-form rank cap. It acts
on nonnegative continuations that increase on caps 1..4; state zero is
independent. This changes the proof, not the encoder or its randomness.
Sparse point replays use rational marks and outward Arb arithmetic. Middle
searches are floating point and never establish a whole-code certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize
from flint import arb, arb_poly, ctx

from cauchy_firstmoment_screen import geometry, log_terminal_power
from field128_sparse_certificate import region
from gf8_packet_maps import CUBIC_RETURN_MINIMA
from gf8_packet_width_screen import joint_terms
from gf8_rank_history_screen import cone_preserved, modal_log_penalty, source_row_log
from gf8_outer_rate_sparse import rank_terminal_power
from state8_sparse_certificate import ball, integer_bits, rational_marks
from state_rank_cap_transfer import rank_cap_weights


WIDTH, STATE = 8, 4
MINIMA = CUBIC_RETURN_MINIMA


def return_bound(weight, y, tau, cap):
    if not 1 <= weight <= WIDTH or not 1 <= cap <= STATE:
        raise ValueError('positive packet weight and state cap required')
    minimum = MINIMA[weight-1]
    pure = (sum(math.comb(STATE, k)/tau**k
                for k in range(min(weight-1, cap-1)+1))
            * y**minimum/tau**STATE)
    joint = sum(count*y**(WIDTH-i-z)*(1-y)**(i+z)/tau**charge
                *sum(math.comb(STATE, k)/tau**k
                     for k in range(min(weight-1, cap-1, STATE-1-z)+1))
                for i, z, count, charge in joint_terms(weight, WIDTH, STATE, MINIMA))
    return min(pure, joint)


def local_matrix(weight, y, tau):
    if not 0 <= weight <= WIDTH or not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('invalid packet parameters')
    y = F(y) if isinstance(y, int) else y
    tau = F(tau) if isinstance(tau, int) else tau
    zero = y*0+tau*0
    result = [[zero for _ in range(STATE+1)] for _ in range(STATE+1)]
    if weight == 0:
        result[0][0] = y**0
    else:
        result[0][min(STATE, weight)] = y**weight
        if weight > STATE:
            result[0][0] = y**weight/tau**STATE
    for cap in range(1, STATE+1):
        if weight:
            result[cap][0] = return_bound(weight, y, tau, cap)
        for i in range(weight+1):
            for z in range(min(STATE-1, WIDTH-weight)+1):
                factor = (math.comb(weight, i)*math.comb(WIDTH-weight, z)
                          *y**(WIDTH-i-z)*(1-y)**(i+z)/tau**i)
                for h, mass in enumerate(rank_cap_weights(cap, z, tau, 256, STATE), 1):
                    result[cap][min(STATE, h+weight-i)] += factor*mass
    return tuple(tuple(row) for row in result)


@lru_cache(maxsize=16)
def _mass_table(tau):
    return np.asarray([[rank_cap_weights(cap, z, tau, 256, STATE)
                        for z in range(STATE)] for cap in range(1, STATE+1)], dtype=float)


def packet(logx, logy, logtau):
    x, y, tau = math.exp(logx), math.exp(logy), min(255., math.exp(logtau))
    result = np.zeros((STATE+1, STATE+1))
    result[0, 0] = 1.
    masses = _mass_table(tau)
    base = y+x*(1-y)/tau
    for z in range(STATE):
        for d in range(WIDTH-z+1):
            coefficient = (math.comb(WIDTH, z)*math.comb(WIDTH-z, d)
                           *(x*y)**d*(1-y)**z*base**(WIDTH-z-d))
            for h in range(1, STATE+1):
                result[1:, min(STATE, h+d)] += coefficient*masses[:, z, h-1]
    for weight in range(1, WIDTH+1):
        multiplicity = math.comb(WIDTH, weight)*x**weight
        birth = multiplicity*y**weight
        result[0, min(STATE, weight)] += birth
        if weight > STATE:
            result[0, 0] += birth/tau**STATE
        for cap in range(1, STATE+1):
            result[cap, 0] += multiplicity*return_bound(weight, y, tau, cap)
    scale = float(result.max())
    if not 0 < scale < math.inf or not cone_preserved(result/scale, 2e-12):
        raise ArithmeticError('continuation cone or floating scale failed')
    return result


def middle_bound(point, weight, dimension=96, k_log=20, components=False, no_returns=False):
    g = geometry(k_log, 128, dimension)
    if not 0 < weight < g['n']:
        raise ValueError('invalid inverse input weight')
    logx, logy, logtau, logtilt = map(float, point)
    matrix = packet(logx, logy, logtau)
    if no_returns:
        # Deliberately optimistic diagnostic, NOT a proved upper bound.
        matrix[1:, 0] = 0.
        matrix[0, 0] = 1.
    terms = dict(route=128*modal_log_penalty(g['rows']),
        source=g['rows']*source_row_log(weight/g['n'], logtau, logtilt, 128, dimension, True),
        inner=log_terminal_power(matrix, g['n']//WIDTH),
        selection=-weight*(logx+logtilt), output=-g['threshold']*logy,
        rank_threshold=-logtau)
    value = sum(terms.values())/math.log(2)
    return (value, {key: val/math.log(2) for key, val in terms.items()}) if components else value


def search_middle(density, dimension=96, k_log=20, iterations=55, start=None, no_returns=False):
    g = geometry(k_log, 128, dimension)
    weight = max(1, round(density*g['n']))
    starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(lambda p: middle_bound(p, weight, dimension, k_log, no_returns=no_returns), p,
                       method='L-BFGS-B', bounds=[(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)],
                       options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25}) for p in starts]
    best = min(choices, key=lambda r: r.fun if math.isfinite(r.fun) else math.inf)
    value, terms = middle_bound(best.x, weight, dimension, k_log, True, no_returns)
    return dict(outer=[128, dimension], width=WIDTH, state=STATE, k_log=k_log,
        input_weight=weight, input_density=weight/g['n'], margin_bits=-value,
        log_parameters=best.x.tolist(), log2_components=terms,
        optimizer_success=bool(best.success), iterations=int(best.nit),
        rigorous=False, whole_code=False, no_returns_diagnostic=no_returns,
        scope=('optimistic impossible-return diagnostic, not a probability bound' if no_returns
               else 'one floating inverse-input-weight class'))


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


def sparse_point(active, point, dimension=96, k_log=20, precision=1024, no_returns=False):
    g = geometry(k_log, 128, dimension)
    if not 1 <= active <= g['rows']:
        raise ValueError('invalid active-row count')
    y, tau = rational_marks(active, point, g['n'])
    previous = ctx.prec
    try:
        ctx.prec = precision
        local = [local_matrix(w, y, tau) for w in range(WIDTH+1)]
        if no_returns:
            local = [[list(row) for row in matrix] for matrix in local]
            for weight, matrix in enumerate(local):
                for cap in range(1, STATE+1):
                    matrix[cap][0] = F(0)
                if weight:
                    matrix[0][0] = F(0)
        coefficients = polynomial_power(local, g['rows']//WIDTH, active)
        value = (arb(math.comb(g['rows'], active))*ball(tau)**(-((128-dimension)*active+1))
                 *ball(y)**(-g['threshold'])
                 *rank_terminal_power(region(coefficients, active, g['rows'], tau), 128))
        endpoint = value.upper()
        if not endpoint.is_finite() or not endpoint > 0:
            raise ArithmeticError('insufficient Arb precision')
        return dict(outer=[128, dimension], width=WIDTH, state=STATE, k_log=k_log,
            active=active, y=str(y), tau=str(tau), precision_bits=precision,
            margin_bits_lower=integer_bits(endpoint),
            log2_upper=str(endpoint.log()/arb(2).log()),
            log2_upper_approx=float(endpoint.log()/arb(2).log()),
            rigorous=not no_returns, whole_code=False,
            no_returns_diagnostic=no_returns,
            scope=('optimistic impossible-return diagnostic, not a probability bound' if no_returns
                   else 'one active-row class at fixed rational marks'))
    finally:
        ctx.prec = previous


def search_sparse(active, start, dimension=96, k_log=20, iterations=35, no_returns=False):
    """Bounded numerical mark search; the final point is replayed separately."""
    def objective(point):
        if not -2 <= point[0] <= 10 or not 0 <= point[1] <= math.log(255):
            return 1e100
        return sparse_point(active, point, dimension, k_log, precision=512,
                            no_returns=no_returns)['log2_upper_approx']
    choice = minimize(objective, start, method='Nelder-Mead',
                      options={'maxiter': iterations, 'xatol': 1e-5, 'fatol': 1e-4})
    result = sparse_point(active, choice.x, dimension, k_log, no_returns=no_returns)
    result.update(log_parameters=choice.x.tolist(), optimizer_success=bool(choice.success),
                  iterations=int(choice.nit))
    return result


def gv_comparison(rate=F(3, 4), distance=F(1, 10)):
    """Reference exponent for an ideal uniform random GF256 linear code."""
    delta = float(distance)
    entropy_bits = -delta*math.log2(delta)-(1-delta)*math.log2(1-delta)+delta*math.log2(255)
    return dict(rate=str(rate), distance=str(distance), qary_entropy_bits=entropy_bits,
        ideal_random_code_exponent_bits_per_output=8*(1-float(rate))-entropy_bits,
        gv_rate=1-entropy_bits/8, scope='asymptotic ideal-ensemble comparison, not a SPIN bound')


def feedthrough_profile(mu=4, eta=2):
    """Exact local branch distances for a changed cheap feedthrough map.

    Set J=CX, Q'=R+J, and Y=X+A(eta R+mu J). Since CA=0, the new
    feedthrough I+mu AC is invertible. At birth Q=0; at return Q'=0.
    Both maps can therefore be checked by finite support enumeration.
    This does not port the existing output-root/rank-pivot argument.
    """
    from gf8_packet_maps import (MUL, cubic_matrices, exact_return_minima,
                                matmul)
    if not 0 < mu < 256 or not 0 < eta < 256 or mu == eta:
        raise ValueError('distinct nonzero GF256 coefficients required')
    a, c = cubic_matrices()
    ac = matmul(a, c)
    for phase in range(1, 8):
        a_phase, c_phase = cubic_matrices(phase)
        if matmul(a_phase, c_phase) != ac:
            raise AssertionError('feedthrough branch table must match every phase')
    def shifted(value):
        return tuple(tuple(int(i == j)^MUL[value][ac[i][j]] for j in range(WIDTH))
                     for i in range(WIDTH))
    birth, returned = shifted(mu), shifted(mu^eta)
    identity = tuple(tuple(int(i == j) for j in range(WIDTH)) for i in range(WIDTH))
    if matmul(birth, birth) != identity or matmul(returned, returned) != identity:
        raise AssertionError('CA=0 should make both maps involutions')
    return dict(mu=mu, eta=eta, birth_minima=list(exact_return_minima(birth)),
                return_minima=list(exact_return_minima(returned)),
                local_algebra_checked=True, whole_code=False,
                proof_status='requires new joint root/rank transfer; identity feedthrough proof does not apply',
                schedule='C pass precedes A pass; old paired five-GFNI512 schedule is not preserved')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension', type=int, nargs='+', default=[96])
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--density', type=float, nargs='+', default=[.03, .06, .10, .14, .20, .28])
    parser.add_argument('--iterations', type=int, default=55)
    parser.add_argument('--sparse', type=int, nargs='+')
    parser.add_argument('--optimize-control', action='store_true')
    parser.add_argument('--gv', action='store_true')
    parser.add_argument('--no-returns', action='store_true',
                        help='optimistic diagnostic only; removes real return events')
    parser.add_argument('--feedthrough', type=int, nargs='+')
    args = parser.parse_args()
    if args.feedthrough:
        for mu in args.feedthrough:
            print(json.dumps(feedthrough_profile(mu)), flush=True)
    elif args.gv:
        print(json.dumps(gv_comparison()), flush=True)
    elif args.sparse:
        from gf8_packet_sparse_certificate import ANCHORS
        from gf8_packet_sparse_screen import search as control_search
        for dimension in args.dimension:
            for active in args.sparse:
                if args.optimize_control:
                    control = control_search(active, args.k_log, dimension, WIDTH, STATE, MINIMA, args.iterations)
                    marks = control['log_parameters']
                else:
                    marks = ANCHORS[min(ANCHORS, key=lambda a: abs(a-active))]
                    control = None
                result = sparse_point(active, marks, dimension, args.k_log, no_returns=args.no_returns)
                if control:
                    result['two_state_margin_bits'] = control['margin_bits']
                    result['log_parameters'] = marks
                print(json.dumps(result), flush=True)
    else:
        for dimension in args.dimension:
            previous = None
            for density in args.density:
                result = search_middle(density, dimension, args.k_log, args.iterations, previous, args.no_returns)
                previous = result['log_parameters']
                print(json.dumps(result), flush=True)
