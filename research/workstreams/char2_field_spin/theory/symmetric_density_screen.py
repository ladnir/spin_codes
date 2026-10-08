"""FLOAT density screen for restricted symmetric Cauchy states1..8.

Restricting the proved sixteen-by-eight A to any s columns and C to the
same transpose rows retains superregularity and CA=0. This screen uses
the generic pure cover or its dimension-independent split-joint refinement.
It does not assume the special eight-state returned-weight law for smaller
states. The GRS source comparison and finite factors are unchanged. Sampled
weights and floating optimization do NOT form a whole-code certificate.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_terminal_power
from superregular_return_screen import local_matrix
from higher_distance_screen import output_threshold
from state_rank_cap_transfer import return_minimum_weight


def minimum_return(weight, state):
    if not 1 <= state <= 8 or not 1 <= weight <= 16:
        raise ValueError('invalid restricted state or live input weight')
    return return_minimum_weight(weight) if state == 8 else max(1, 17-state-weight)


@lru_cache(maxsize=128)
def joint_terms(weight, state):
    minimum = minimum_return(weight, state)
    return tuple((i, z, math.comb(weight, i)*math.comb(16-weight, z),
                  max(state, i+z, weight-16+i+2*z+1))
                 for i in range(weight+1)
                 for z in range(min(state-1, 16-weight)+1)
                 if 16-i-z >= minimum)


def return_bound(weight, y, tau, state):
    """One fresh-diagonal pivot charge for the combined zero/return matrix.

    For incoming live support h and selected inactive zeros z, z<h.
    Return rank is s, selected-output rank at least i+z, and invertibility
    gives joint rank at least h+w-v. Zero-flat ranks are charged before
    current diagonals and satisfy k<=min(w-1,s-1-z). No marginal-root
    factor is multiplied into this combined cover.
    """
    if not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('invalid output mark or rank tilt')
    if isinstance(y, int):
        y = Fraction(y)
    if isinstance(tau, int):
        tau = Fraction(tau)
    flats = sum(math.comb(state,k)/tau**k for k in range(min(weight-1,state-1)+1))
    pure = flats*y**minimum_return(weight,state)/tau**state
    joint = sum(count*y**(16-i-z)*(1-y)**(i+z)/tau**rank
                *sum(math.comb(state,k)/tau**k for k in range(min(weight-1,state-1-z)+1))
                for i,z,count,rank in joint_terms(weight,state))
    return min(pure,joint)


def packet(logx, logy, logtau, state, joint=False):
    if not 1 <= state <= 8:
        raise ValueError('only restrictions of the checked eight-state map are covered')
    y, tau = math.exp(logy), min(255., math.exp(logtau))
    values = np.asarray([local_matrix(w, y, tau, state=state, width=16)
                         for w in range(17)], dtype=float)
    if joint:
        for w in range(1,17):
            values[w,1,0] = return_bound(w,y,tau,state)
    with np.errstate(divide='ignore'):
        logs = np.log(values)
    for w in range(17):
        logs[w] += math.log(math.comb(16, w))+w*logx
    return logsumexp(logs, axis=0)


def bound(point, weight, k_log=20, length=128, dimension=86, state=6,
          distance=Fraction(1, 10), joint=False):
    g = geometry(k_log, length, dimension)
    if not 0 < dimension < length <= 256 or not 0 < weight < g['n']:
        raise ValueError('invalid GRS geometry or inverse input weight')
    alpha = weight/g['n']
    logx, logy, logtau, logtilt = map(float, point)
    row = float(logsumexp([length*math.log1p(-alpha),
                          -(length-dimension)*logtau+length*np.logaddexp(
                              math.log1p(-alpha), logtau+logtilt+math.log(alpha))]))
    logarithms = packet(logx, logy, logtau, state, joint)
    shift = float(logarithms.max())
    steps = g['n']//16
    inner = steps*shift+log_terminal_power(np.exp(logarithms-shift), steps)
    return (length*math.log(g['rows']+1)+g['rows']*row+inner
            -weight*(logx+logtilt)-output_threshold(g['n'], distance)*logy-logtau)/math.log(2)


def search(alpha, k_log=20, length=128, dimension=86, state=6,
           distance=Fraction(1, 10), iterations=55, start=None, joint=False):
    if not 0 < alpha < 1:
        raise ValueError('inverse input density must lie in (0,1)')
    g = geometry(k_log, length, dimension)
    weight = max(1, round(alpha*g['n']))
    starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
    if start is not None:
        starts.insert(0, start)
    options = [minimize(lambda p: bound(p, weight, k_log, length, dimension, state, distance, joint),
                        p, method='L-BFGS-B',
                        bounds=[(-30, 12), (-25, 0), (0, math.log(255)), (-20, 20)],
                        options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 20})
               for p in starts]
    best = min(options, key=lambda v: v.fun if np.isfinite(v.fun) else math.inf)
    return dict(outer=[length, dimension], state=state, width=16, k_log=k_log,
                rate=str(Fraction(dimension, length)), distance=str(Fraction(distance)),
                input_weight=weight, input_density=weight/g['n'], margin_bits=-float(best.fun),
                log_parameters=best.x.tolist(), optimizer_success=bool(best.success),
                iterations=int(best.nit), rigorous=False, whole_code=False,
                joint_return=joint, scope='one sampled inverse-input-weight class')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, default=(128, 86))
    parser.add_argument('--state', type=int, default=6)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--distance', type=Fraction, default=Fraction(1, 10))
    parser.add_argument('--density', type=float, nargs='+', default=[.02, .04, .06, .1, .2, .3])
    parser.add_argument('--iterations', type=int, default=55)
    parser.add_argument('--joint', action='store_true')
    args = parser.parse_args()
    previous = None
    for alpha in args.density:
        value = search(alpha, args.k_log, *args.outer, args.state, args.distance, args.iterations, previous, args.joint)
        previous = value['log_parameters']
        print(json.dumps(value), flush=True)
