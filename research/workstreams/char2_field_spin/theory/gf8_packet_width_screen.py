"""Floating packet-width screen for self-orthogonal superregular GF256 maps.

Assumes A is width-by-state superregular, C is state-by-width
superregular, CA=0, and eta!=0. The matrices need not be transposes.
Independent nonzero input and state diagonals and one regional route are
required. This is a candidate-class screen, never a whole-code certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_terminal_power
from diagonal_refresh_rank import rank_potential_bound
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from higher_distance_screen import output_threshold


def validate(width, state):
    if width not in (8, 16, 32) or not 1 <= state <= width//2:
        raise ValueError('supported self-orthogonal geometries have width8/16/32 and state<=width/2')


def minimum_return(weight, width, state, minima=()):
    validate(width, state)
    if not 1 <= weight <= width:
        raise ValueError('return requires nonempty input')
    generic = max(1, width-state+1-weight)
    supplied = minima[weight-1] if weight <= len(minima) else generic
    if not generic <= supplied <= width:
        raise ValueError('provided return minimum must strengthen the generic bound')
    return supplied


@lru_cache(maxsize=256)
def joint_terms(weight, width, state, minima=()):
    minimum = minimum_return(weight, width, state, minima)
    return tuple((i, z, math.comb(weight, i)*math.comb(width-weight, z),
                  max(state, i+z, weight-width+i+2*z+1))
                 for i in range(weight+1)
                 for z in range(min(state-1, width-weight)+1)
                 if width-i-z >= minimum)


@lru_cache(maxsize=128, typed=True)
def gamma(tau, state):
    return tuple(max(rank_potential_bound(rank, zeros, tau, 256, state)
                     for rank in range(1, state+1)) for zeros in range(state))


def return_bound(weight, y, tau, width, state, minima=()):
    """Min of complete pure/joint covers, using one charge per diagonal.

    A return needs h+w>s for incoming support h. Selected inactive zeros
    satisfy z<h. The joint equation rank is >=s, >=i+z, and >=h+w-v,
    where v=width-i-z. The flat basis has rank <=min(w-1,s-1-z).
    """
    minimum = minimum_return(weight, width, state, minima)
    if not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('invalid packet tilts')
    y = F(y) if isinstance(y, int) else y
    tau = F(tau) if isinstance(tau, int) else tau
    flats = sum(math.comb(state, k)/tau**k for k in range(min(weight-1, state-1)+1))
    pure = flats*y**minimum/tau**state
    joint = sum(count*y**(width-i-z)*(1-y)**(i+z)/tau**rank
                *sum(math.comb(state, k)/tau**k
                     for k in range(min(weight-1, state-1-z)+1))
                for i, z, count, rank in joint_terms(weight, width, state, minima))
    return min(pure, joint)


def local_matrix(weight, y, tau, width, state, minima=()):
    validate(width, state)
    if not 0 <= weight <= width or not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('invalid local matrix parameters')
    zero = y*0+tau*0
    result = [[zero, zero], [zero, zero]]
    if not weight:
        result[0][0] = y**0
    else:
        result[0][1] = y**weight
        if weight > state:
            result[0][0] = y**weight/tau**state
        result[1][0] = return_bound(weight, y, tau, width, state, minima)
    roots = gamma(tau, state)
    result[1][1] = sum(math.comb(weight, i)*math.comb(width-weight, z)
                       *y**(width-i-z)*(1-y)**(i+z)*roots[z]/tau**i
                       for i in range(weight+1)
                       for z in range(min(state-1, width-weight)+1))
    return result


@lru_cache(maxsize=64)
def _return_arrays(width, state, minima):
    arrays = []
    for weight in range(1, width+1):
        rows = [(width-i-z, i+z, rank, min(weight-1, state-1-z), count)
                for i, z, count, rank in joint_terms(weight, width, state, minima)]
        arrays.append(np.asarray(rows, dtype=float).T)
    return arrays


def packet(logx, logy, logtau, width, state, minima=()):
    """Aggregate C(width,w)x^w T_w, with the nonreturn sum batched exactly."""
    validate(width, state)
    x, y, tau = math.exp(logx), math.exp(logy), min(255., math.exp(logtau))
    roots = np.asarray(gamma(tau, state), dtype=float)
    powers_y = np.asarray([y**i for i in range(width+1)])
    powers_z = np.asarray([(1-y)**i for i in range(width+1)])
    powers_tau = np.asarray([tau**(-i) for i in range(2*width+1)])
    flats = np.cumsum([math.comb(state, k)*powers_tau[k] for k in range(state)])
    birth = [math.comb(width, w)*(x*y)**w for w in range(width+1)]
    result = np.zeros((2, 2))
    result[0, 0] = 1+sum(birth[state+1:])*powers_tau[state]
    result[0, 1] = sum(birth[1:])
    base = y+x*y+x*(1-y)/tau
    result[1, 1] = sum(math.comb(width, z)*(1-y)**z*base**(width-z)*roots[z]
                       for z in range(state))
    for weight, rows in enumerate(_return_arrays(width, state, tuple(minima)), 1):
        yp, zp, charges, caps, counts = rows
        joint = np.sum(counts*powers_y[yp.astype(int)]*powers_z[zp.astype(int)]
                       *powers_tau[charges.astype(int)]*flats[caps.astype(int)])
        pure = flats[min(weight-1, state-1)]*powers_y[minimum_return(weight, width, state, minima)]*powers_tau[state]
        result[1, 0] += math.comb(width, weight)*x**weight*min(pure, joint)
    with np.errstate(divide='ignore'):
        return np.log(result)


def bound(point, weight, k_log=20, dimension=88, width=8, state=4,
          distance=F(1, 10), minima=(), components=False):
    g = geometry(k_log, 128, dimension)
    if g['rows'] % width or not 0 < weight < g['n']:
        raise ValueError('unaligned regional length or invalid inverse input weight')
    logx, logy, logtau, logtilt = map(float, point)
    alpha = weight/g['n']
    row = source_row_log(alpha, logtau, logtilt, 128, dimension, True)
    logs = packet(logx, logy, logtau, width, state, minima)
    scale = float(logs.max())
    steps = g['n']//width
    inner = steps*scale+log_terminal_power(np.exp(logs-scale), steps)
    terms = dict(route=128*modal_log_penalty(g['rows']), source=g['rows']*row,
                 inner=inner, selection=-weight*(logx+logtilt),
                 output=-output_threshold(g['n'], distance)*logy, rank_threshold=-logtau)
    value = sum(terms.values())/math.log(2)
    return (value, {key: v/math.log(2) for key, v in terms.items()}) if components else value


def search(density, k_log=20, dimension=88, width=8, state=4,
           iterations=60, start=None, minima=()):
    g = geometry(k_log, 128, dimension)
    weight = max(1, round(density*g['n']))
    starts = [(-2., -.2, math.log(16), 0.), (0., -1., math.log(255), -1.)]
    if start is not None:
        starts.insert(0, start)
    candidates = [minimize(lambda p: bound(p, weight, k_log, dimension, width, state, minima=minima),
                           p, method='L-BFGS-B',
                           bounds=[(-30, 8), (-20, 0), (0, math.log(255)), (-20, 20)],
                           options={'maxiter': iterations, 'ftol': 1e-11, 'maxls': 25})
                  for p in starts]
    best = min(candidates, key=lambda r: r.fun if np.isfinite(r.fun) else math.inf)
    value, terms = bound(best.x, weight, k_log, dimension, width, state,
                         minima=minima, components=True)
    return dict(outer=[128, dimension], width=width, state=state, q=256,
                k_log=k_log, requested_k=2**k_log, padded_k=g['k'], n=g['n'], rows=g['rows'],
                requested_rate=str(F(2**k_log, g['n'])), nominal_rate=str(F(dimension, 128)),
                distance='1/10', input_weight=weight, input_density=weight/g['n'],
                margin_bits=-value, log_parameters=best.x.tolist(), log2_components=terms,
                optimizer_success=bool(best.success), iterations=int(best.nit),
                return_minima=list(minima), rigorous=False, whole_code=False,
                scope='one sampled inverse-input-weight class for the stated map assumptions')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width', type=int, default=8)
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--dimension', type=int, nargs='+', default=[87, 88, 96])
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--density', type=float, nargs='+', default=[.03, .06, .10, .14, .20, .24])
    parser.add_argument('--iterations', type=int, default=60)
    parser.add_argument('--return-minima', type=int, nargs='+', default=[])
    args = parser.parse_args()
    for dimension in args.dimension:
        previous = None
        for density in args.density:
            result = search(density, args.k_log, dimension, args.width, args.state,
                            args.iterations, previous, tuple(args.return_minima))
            previous = result['log_parameters']
            print(json.dumps(result), flush=True)
