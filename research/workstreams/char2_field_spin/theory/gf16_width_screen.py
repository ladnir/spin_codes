"""Floating GF16 screens with a shorter inner group and unchanged outer/route.

GRS[16,8], sixteen independently permuted regions, uniform nonzero diagonal,
and uniform GL(s,16) pre-emission refresh. A uses monomials on the first t
abstract GF16 elements and C=A^T. These screens are not certificates.
"""

import argparse
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

from batched_transfer import multiply
from gf16_sizes import geometry, source_eta
from q1_batched_screen import fixed_support_multiplicity
from reduced_distance_certificate import log_polynomial_product
from state_screen import local as generic_local


@lru_cache(maxsize=None)
def singleton_weights(width, state):
    if width not in (8, 16) or not 1 <= state <= width//2:
        raise ValueError('CA=0 is checked only for these additive-point profiles')
    powers = []
    for j in range(width):
        row = [1]
        for _ in range(state-1):
            row.append(multiply(row[-1], j, 0x13))
        powers.append(row)
    for a in range(state):
        for b in range(state):
            total = 0
            for row in powers:
                total ^= multiply(row[a], row[b], 0x13)
            if total:
                raise ValueError('CA is not zero for these evaluation points')
    weights = []
    for active in range(width):
        weight = 0
        for out in range(width):
            value = int(out == active)
            for p in range(state):
                value ^= multiply(powers[active][p], powers[out][p], 0x13)
            weight += value != 0
        weights.append(weight)
    return tuple(weights)


def local(z, width, state):
    result = generic_local(z, width, state, 16)
    denominator = 16**state-1
    joint = sum(z**v for v in singleton_weights(width, state))/(width*denominator)
    result[1, 1, 1] += z**width/denominator-joint
    result[1, 1, 0] = joint
    for w in range(2, width+1):
        inverse = sum(math.comb(width, h)*z**(width-h)*(1-z)**h
                      /15**max(0, h-width+w) for h in range(width+1))/denominator
        result[w, 1, 0] = min(result[w, 1, 0], inverse)
    return result


def screen(k_log, active, z, width=8, state=4, rho=.871):
    rows, _, threshold = geometry(k_log)
    if rows % width or not 1 <= active < rows:
        raise ValueError('invalid regional geometry')
    with np.errstate(divide='ignore'):
        factors = np.log(local(z, width, state)[:active+1])
    for w in range(len(factors)):
        factors[w] += math.log(math.comb(width, w))
    regions = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows//width
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, factors, active)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, active)
    for h in range(active+1):
        regions[h] -= math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
    if active <= 4:
        # Exact support counts; this is the same outer as the t16 certificate.
        matrices = [np.exp(regions[h])*16**(h-active/2) for h in range(active+1)]
        masks = [tuple((mask >> j) & 1 for j in range(active)) for mask in range(1 << active)]
        current = np.zeros((1,)*active+(2,))
        current[(0,)*active+(0,)] = 1
        scale = 0.
        for region in range(16):
            next_state = np.zeros((region+2,)*active+(2,))
            for mask in masks:
                target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
                next_state[target] += current @ matrices[sum(mask)]
            norm = float(np.max(next_state))
            current = next_state/norm
            scale += math.log(norm)
        weights = np.array([fixed_support_multiplicity(16, 8, w, 16)/16**(w-8) if w >= 9 else 0
                            for w in range(17)])
        for axis in range(active):
            shape = [1]*(active+1)
            shape[axis] = 17
            current *= weights.reshape(shape)
        moment_log = scale+math.log(float(np.sum(current)))
    else:
        source = np.array([math.log(math.comb(active, h))+(h-active/2)*math.log(16)+h*math.log(rho)
                           for h in range(active+1)])
        matrix = logsumexp(regions+source[:, None, None], axis=0)
        vector = np.array([0., -math.inf])
        for _ in range(16):
            vector = logsumexp(vector[:, None]+matrix, axis=0)
        moment_log = float(logsumexp(vector))+active*math.log(float(source_eta(rho)))
    choose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    return -(choose-math.log(15)-threshold*math.log(z)+moment_log)/math.log(2)


def search(k_log, active, width, state, optimize_source=False):
    _, n, _ = geometry(k_log)
    def objective(lt, lr=math.log(.871)):
        z = math.exp(-math.exp(lt)*active/n)
        return -screen(k_log, active, z, width, state, math.exp(lr))
    optimum = minimize_scalar(objective, bounds=(1, 8), method='bounded')
    lt, lr, result = optimum.x, math.log(.871), -optimum.fun
    if optimize_source and active > 4:
        def joint(point):
            if not 1 <= point[0] <= 8 or not -2 <= point[1] <= 2:
                return 1e100
            return objective(*point)
        better = minimize(joint, [lt, lr], method='Nelder-Mead',
                          options=dict(maxiter=150, xatol=1e-6, fatol=1e-6))
        lt, lr, result = *better.x, -better.fun
    return dict(k_log=k_log, active_rows=active, width=width, state=state,
                z=math.exp(-math.exp(lt)*active/n), source_rho=math.exp(lr),
                margin_bits=float(result), rigorous=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, choices=[16, 18, 20], default=20)
    parser.add_argument('--width', type=int, choices=[8, 16], default=8)
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--rows', nargs='+', type=int, default=[1, 4, 8, 16, 32, 64, 128, 256])
    parser.add_argument('--optimize-source', action='store_true')
    args = parser.parse_args()
    for active in args.rows:
        print(json.dumps(search(args.k_log, active, args.width, args.state, args.optimize_source)), flush=True)
