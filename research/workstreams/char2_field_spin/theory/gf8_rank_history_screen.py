"""Floating GF256 density screen retaining the incoming state-rank cap.

The positive transfer acts on nonnegative continuations nondecreasing in
the nonzero rank cap. State zero denotes the actual zero witness state.
The fixed source quotient and all fresh-diagonal charges are those of
state_rank_cap_transfer and symmetric_density_screen; there is exactly one
regional route. Individual optimized density classes are not certificates.
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

from cauchy_firstmoment_screen import geometry, log_terminal_power
from higher_distance_screen import output_threshold
from state_rank_cap_transfer import rank_cap_weights, return_flat_factor
from symmetric_density_screen import joint_terms, minimum_return


def return_bound(weight, y, tau, cap, state=6, joint=True, cap_joint=True, return_minima=()):
    """One diagonal charge, with the incoming zero flat fixed beforehand.

    A proper zero flat of an incoming form space of rank r has rank k<r,
    hence k<=cap-1. Selected inactive zeros additionally require k<s-z.
    Both complete covers increase with cap, so their minimum does too.
    cap_joint=False retains the older uniform joint bound as a control.
    """
    if not (state in (4, 6, 8) and 1 <= cap <= state and
            1 <= weight <= 16 and 0 < y <= 1 and 1 <= tau <= 255):
        raise ValueError('invalid return parameters')
    y = Fraction(y) if isinstance(y, int) else y
    tau = Fraction(tau) if isinstance(tau, int) else tau
    minimum = (return_minima[weight-1] if weight <= len(return_minima)
               else minimum_return(weight, state))
    if not minimum_return(weight, state) <= minimum <= 16:
        raise ValueError('return minimum override must strengthen the standard bound')
    pure = return_flat_factor(weight, tau, cap, state)*y**minimum/tau**state
    if not joint:
        return pure
    flat_cap = cap if cap_joint else state
    combined = sum(count*y**(16-i-z)*(1-y)**(i+z)/tau**rank
                   *sum(math.comb(state, k)/tau**k
                        for k in range(min(weight-1, flat_cap-1, state-1-z)+1))
                   for i, z, count, rank in joint_terms(weight, state) if 16-i-z >= minimum)
    return min(pure, combined)


def local_matrix(weight, y, tau, state=6, joint=True, cap_joint=True, return_minima=()):
    """Exact-capable packet transfer for restricted states 4, 6, and 8."""
    if not (state in (4, 6, 8) and 0 <= weight <= 16 and
            0 < y <= 1 and 1 <= tau <= 255):
        raise ValueError('invalid packet parameters')
    y = Fraction(y) if isinstance(y, int) else y
    tau = Fraction(tau) if isinstance(tau, int) else tau
    zero = y*0+tau*0
    result = [[zero for _ in range(state+1)] for _ in range(state+1)]
    if not weight:
        result[0][0] = y**0
    else:
        result[0][min(state, weight)] = y**weight
        if weight > state:
            result[0][0] = y**weight/tau**state
    for cap in range(1, state+1):
        if weight:
            result[cap][0] = return_bound(weight, y, tau, cap, state, joint, cap_joint, return_minima)
        for i in range(weight+1):
            for z in range(min(state-1, 16-weight)+1):
                factor = (math.comb(weight, i)*math.comb(16-weight, z)
                          *y**(16-i-z)*(1-y)**(i+z)/tau**i)
                if factor:
                    for h, mass in enumerate(rank_cap_weights(cap, z, tau, 256, state), 1):
                        result[cap][min(state, h+weight-i)] += factor*mass
    return tuple(tuple(row) for row in result)


def cone_preserved(matrix, tolerance=0):
    """Test all rays of the continuation cone, including independent f(0)."""
    matrix = np.asarray(matrix)
    if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]
            or any(x < -tolerance for x in matrix.flat)):
        return False
    if any(matrix[cap, 0] > matrix[cap+1, 0]+tolerance
           for cap in range(1, len(matrix)-1)):
        return False
    return all(sum(matrix[cap, first:]) <= sum(matrix[cap+1, first:])+tolerance
               for cap in range(1, len(matrix)-1) for first in range(1, len(matrix)))


@lru_cache(maxsize=16)
def _mass_table(tau, state):
    return np.asarray([[rank_cap_weights(cap, z, tau, 256, state)
                        for z in range(state)] for cap in range(1, state+1)], dtype=float)


def packet(logx, logy, logtau, state=6, joint=True, cap_joint=True, return_minima=()):
    """Log aggregate sum_w C(16,w)x^w T_w, using a batched binomial sum.

    For nonreturn terms put d=w-i. Summing the selected active positions
    yields C(16,z)C(16-z,d)(xy)^d(1-y)^z
    (y+x(1-y)/tau)^(16-z-d), avoiding 17 packet matrix constructions.
    """
    if state not in (4, 6, 8):
        raise ValueError('unsupported restricted state')
    x, y, tau = math.exp(logx), math.exp(logy), min(255., math.exp(logtau))
    if not (0 < y <= 1 and 1 <= tau <= 255):
        raise ValueError('invalid packet tilts')
    matrix = np.zeros((state+1, state+1))
    matrix[0, 0] = 1.
    masses = _mass_table(tau, state)
    base = y+x*(1-y)/tau
    for z in range(state):
        for d in range(17-z):
            coefficient = (math.comb(16, z)*math.comb(16-z, d)
                           *(x*y)**d*(1-y)**z*base**(16-z-d))
            if coefficient:
                for h in range(1, state+1):
                    matrix[1:, min(state, h+d)] += coefficient*masses[:, z, h-1]
    for weight in range(1, 17):
        multiplicity = math.comb(16, weight)*x**weight
        birth = multiplicity*y**weight
        matrix[0, min(state, weight)] += birth
        if weight > state:
            matrix[0, 0] += birth/tau**state
        for cap in range(1, state+1):
            matrix[cap, 0] += multiplicity*return_bound(weight, y, tau, cap, state, joint, cap_joint, return_minima)
    scale = float(matrix.max())
    if not 0 < scale < math.inf or not cone_preserved(matrix/scale, 2e-12):
        raise ArithmeticError('packet continuation cone or scaling failed')
    with np.errstate(divide='ignore'):
        return np.log(matrix)


@lru_cache(maxsize=16)
def _log_binomial_row(length):
    return np.asarray([math.log(math.comb(length, h)) for h in range(length+1)])


def source_row_log(alpha, logtau, loglambda, length, dimension, feasible=False):
    """Positive row sum, optionally excluding infeasible nonempty GRS supports.

    Homogenizing unequal regional densities requires beta=tau*lambda>=2:
    then the coefficient sequence 1,0,...,0,tau^-p*beta^(p+1),... is
    discretely convex, so pairwise averaging raises its Bernoulli moment.
    Below beta=2 retain the established relaxed two-product envelope.
    """
    if not 0 < alpha < 1 or not 0 < dimension < length:
        raise ValueError('invalid source geometry')
    parity = length-dimension
    logzero = math.log1p(-alpha)
    logbeta = logtau+loglambda
    empty = length*logzero
    if not feasible or logbeta < math.log(2):
        occupied = -parity*logtau+length*np.logaddexp(logzero, logbeta+math.log(alpha))
        return float(np.logaddexp(empty, occupied))
    h = np.arange(parity+1, length+1)
    terms = (_log_binomial_row(length)[parity+1:]-parity*logtau
             +h*(logbeta+math.log(alpha))+(length-h)*logzero)
    return float(logsumexp(np.concatenate(([empty], terms))))


@lru_cache(maxsize=16)
def modal_log_penalty(rows):
    """Log exact central-mode reciprocal, uniform over each column's mean."""
    if rows <= 0 or rows % 2:
        raise ValueError('positive even regional length required')
    return math.log(float(Fraction(2**rows, math.comb(rows, rows//2))))


def bound(point, weight, k_log=20, length=128, dimension=86, state=6,
          distance=Fraction(1, 10), joint=True, cap_joint=True, components=False, return_minima=(),
          feasible_source=False, modal=False):
    """Finite density-class source and one-route comparison from the control."""
    g = geometry(k_log, length, dimension)
    if not 0 < dimension < length <= 256 or not 0 < weight < g['n']:
        raise ValueError('invalid GRS geometry or inverse input weight')
    alpha = weight/g['n']
    logx, logy, logtau, logtilt = map(float, point)
    row = source_row_log(alpha, logtau, logtilt, length, dimension, feasible_source)
    logarithms = packet(logx, logy, logtau, state, joint, cap_joint, return_minima)
    shift = float(logarithms.max())
    steps = g['n']//16
    inner = steps*shift+log_terminal_power(np.exp(logarithms-shift), steps)
    route = length*(modal_log_penalty(g['rows']) if modal else math.log(g['rows']+1))
    terms = dict(route=route, source=g['rows']*row,
                 inner=inner, selection=-weight*(logx+logtilt),
                 output=-output_threshold(g['n'], distance)*logy, rank_threshold=-logtau)
    result = sum(terms.values())/math.log(2)
    return (result, {name: value/math.log(2) for name, value in terms.items()}) if components else result


def search(alpha, k_log=20, length=128, dimension=86, state=6,
           distance=Fraction(1, 10), iterations=55, start=None, joint=True,
           cap_joint=True, compare=False, return_minima=(), feasible_source=False, modal=False):
    if not 0 < alpha < 1:
        raise ValueError('input density must lie in (0,1)')
    g = geometry(k_log, length, dimension)
    weight = max(1, round(alpha*g['n']))
    starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
    control = None
    if compare:
        from symmetric_density_screen import search as control_search
        control = control_search(alpha, k_log, length, dimension, state, distance, iterations, joint=joint)
        starts.insert(0, control['log_parameters'])
    if start is not None:
        starts.insert(0, start)
    choices = [minimize(lambda p: bound(p, weight, k_log, length, dimension, state, distance, joint, cap_joint,
                                       return_minima=return_minima, feasible_source=feasible_source, modal=modal),
                        p, method='L-BFGS-B',
                        bounds=[(-30, 12), (-25, 0), (0, math.log(255)), (-20, 20)],
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 20}) for p in starts]
    best = min(choices, key=lambda choice: choice.fun if np.isfinite(choice.fun) else math.inf)
    value, terms = bound(best.x, weight, k_log, length, dimension, state, distance, joint, cap_joint, True,
                         return_minima, feasible_source, modal)
    result = dict(outer=[length, dimension], state=state, width=16, field=256, k_log=k_log,
                  rate=str(Fraction(dimension, length)), distance=str(Fraction(distance)),
                  input_weight=weight, input_density=weight/g['n'], margin_bits=-value,
                  log_parameters=best.x.tolist(), log2_components=terms,
                  optimizer_success=bool(best.success), iterations=int(best.nit),
                  rigorous=False, whole_code=False, joint_return=joint, cap_joint=cap_joint,
                  return_minima=list(return_minima),
                  feasible_source=feasible_source, modal_conditioning=modal,
                  source_tail_used=bool(feasible_source and best.x[2]+best.x[3] >= math.log(2)),
                  scope='one sampled inverse-input-weight class')
    if control is not None:
        result['two_state_margin_bits'] = control['margin_bits']
        result['rank_history_at_control_margin_bits'] = -bound(
            control['log_parameters'], weight, k_log, length, dimension, state, distance, joint, cap_joint,
            return_minima=return_minima, feasible_source=feasible_source, modal=modal)
    if feasible_source or modal:
        row_old = source_row_log(weight/g['n'], best.x[2], best.x[3], length, dimension)
        row_new = source_row_log(weight/g['n'], best.x[2], best.x[3], length, dimension, feasible_source)
        result['source_gain_at_parameters_bits'] = g['rows']*(row_old-row_new)/math.log(2)
        result['conditioning_gain_bits'] = (length*(math.log(g['rows']+1)-modal_log_penalty(g['rows']))/math.log(2)
                                             if modal else 0.)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, default=(128, 86))
    parser.add_argument('--state', type=int, choices=(4, 6, 8), default=6)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--distance', type=Fraction, default=Fraction(1, 10))
    parser.add_argument('--density', type=float, nargs='+', default=[.10, .14])
    parser.add_argument('--iterations', type=int, default=55)
    parser.add_argument('--pure', action='store_true')
    parser.add_argument('--uncapped-joint', action='store_true')
    parser.add_argument('--compare', action='store_true')
    parser.add_argument('--return-minima', type=int, nargs='+', default=[],
                        help='independently proved return minima for weights 1, 2, ...')
    parser.add_argument('--feasible-source', action='store_true',
                        help='truncate infeasible row supports when tau*lambda>=2; otherwise use old envelope')
    parser.add_argument('--modal', action='store_true', help='use the sharper central-mode conditioning factor')
    args = parser.parse_args()
    previous = None
    for alpha in args.density:
        value = search(alpha, args.k_log, *args.outer, args.state, args.distance,
                       args.iterations, previous, not args.pure, not args.uncapped_joint, args.compare,
                       tuple(args.return_minima), args.feasible_source, args.modal)
        previous = value['log_parameters']
        print(json.dumps(value), flush=True)
