"""Floating screens for higher-rate rectangular Toeplitz and GRS outers.

Each row maps u in GF16^k to (u, D*u), where D has p=b-k rows and k
columns. Every row has an independent uniform Toeplitz matrix. A fixed
nonzero u has uniform D*u in GF16^p, giving the complete averaged value
measure 16^-p times counting measure on u!=0 and arbitrary parity values.
No input diagonal is sampled. One independently permuted region per outer
coordinate feeds the width-16 monomial inner with independent GL(s,16).

These computations are floating screens, never distance certificates.
The optional fixed GRS[16,12] model instead requires independent nonzero
GF16 input multipliers and uses the exact MDS support multiplicities.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

from coupled_chernoff import log_terminal_power
from gf16_sizes import local_model
from reduced_distance_certificate import log_polynomial_product
from q1_batched_screen import fixed_support_multiplicity
from state_screen import birth_zero_count


def geometry(k_log, length, dimension):
    if k_log not in (16, 18, 20) or not 0 < dimension < length or length % 16:
        raise ValueError("require K16/K18/K20 and 0 < dimension < length divisible by 16")
    multiple = dimension*16
    k = (((1 << k_log)+multiple-1)//multiple)*multiple
    rows = k//dimension
    n = rows*length
    return k, n, rows, (n+9)//10-1


def log_choose(n, k):
    if not 0 <= k <= n:
        raise ValueError("invalid binomial arguments")
    k = min(k, n-k)
    if k <= 64:
        return sum(math.log(n-j)-math.log(j+1) for j in range(k))
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)


def source_moment(regions, active, length, dimension, exact_systematic=True):
    """Log of the dominating source-weighted regional matrix product."""
    source = np.array([log_choose(active, h)+h*math.log(15) for h in range(active+1)])
    matrix = logsumexp(regions+source[:, None, None], axis=0)
    if exact_systematic:
        # Exchangeability reduces the set of rows already seen to its size.
        # From t seen rows, choose i old and j new active coordinates.
        # All entries are positive; no inclusion-exclusion is used.
        transition = np.full((active+1, 2, active+1, 2), -math.inf)
        for seen in range(active+1):
            for new in range(active-seen+1):
                terms = np.array([regions[old+new]+log_choose(seen, old)
                                  +(old+new)*math.log(15) for old in range(seen+1)])
                transition[seen, :, seen+new, :] = logsumexp(terms, axis=0)+log_choose(active-seen, new)
        current = np.full((active+1, 2), -math.inf)
        current[0, 0] = 0.
        for _ in range(dimension):
            current = logsumexp(current[:, :, None, None]+transition, axis=(0, 1))
        vector = current[active]
        remaining = length-dimension
    else:
        vector = np.array([0., -math.inf])
        remaining = length
    for _ in range(remaining):
        vector = logsumexp(vector[:, None]+matrix, axis=0)
    return float(logsumexp(vector))-(length-dimension)*active*math.log(16)


def sparse_screen(k_log, active, state, length, dimension, z, exact_systematic=True):
    _, _, rows, threshold = geometry(k_log, length, dimension)
    if not 1 <= active < rows or not 1 <= state <= 8:
        raise ValueError("invalid active class or state size")
    with np.errstate(divide='ignore'):
        factors = np.log(local_model(z, state)[:active+1])
    for w in range(len(factors)):
        factors[w] += math.log(math.comb(16, w))
    regions = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows//16
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, factors, active)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, active)
    for h in range(active+1):
        regions[h] -= log_choose(rows, h)
    moment = source_moment(regions, active, length, dimension, exact_systematic)
    return -(log_choose(rows, active)-math.log(15)-threshold*math.log(z)+moment)/math.log(2)


def sparse_search(k_log, active, state, length, dimension):
    k, n, rows, _ = geometry(k_log, length, dimension)
    def objective(logtilt):
        z = math.exp(-math.exp(logtilt)*active/n)
        return -sparse_screen(k_log, active, state, length, dimension, z)
    optimum = minimize_scalar(objective, bounds=(-2, 12), method='bounded', options={'xatol': 2e-6})
    return dict(k=k, n=n, rows=rows, active_rows=active, state=state, width=16,
                outer=[length, dimension], margin_bits=-float(optimum.fun),
                z=math.exp(-math.exp(optimum.x)*active/n), rigorous=False,
                source='exact-nonempty-systematic-all-active-rows')


def middle_screen(k_log, active, state, length, dimension, x, z, aggregate=None):
    _, _, rows, threshold = geometry(k_log, length, dimension)
    if not 1 <= active <= rows or not math.isfinite(x) or x <= 0:
        raise ValueError("invalid active class or source tilt")
    if aggregate is None:
        transfer = local_model(z, state)
        matrix = sum(math.comb(16, w)*(15*x)**w*(1+x)**(16-w)*transfer[w] for w in range(17))
    else:
        matrix = (1+x)**16*np.asarray(aggregate(x/(1+x), z, state))
    value = (-(length-1)*log_choose(rows, active)-(length-dimension)*active*math.log(16)
             -math.log(15)-threshold*math.log(z)-length*active*math.log(x)
             +log_terminal_power(matrix, rows*length//16))
    return -value/math.log(2)


def middle_search(k_log, active, state, length, dimension, aggregate=None):
    k, n, rows, _ = geometry(k_log, length, dimension)
    def objective(point):
        logx, logz = point
        if not -35 < logx < 10 or not -20 < logz < 0:
            return 1e100
        return -middle_screen(k_log, active, state, length, dimension, math.exp(logx), math.exp(logz), aggregate)
    best = None
    for initial_z in (.05, .5, .9, .99, .999):
        initial_x = active/(max(1, rows-active)*16)
        optimum = minimize(objective, [math.log(initial_x), math.log(initial_z)], method='Nelder-Mead',
                           options={'maxiter': 400, 'xatol': 1e-9, 'fatol': 1e-7})
        candidate = dict(k=k, n=n, rows=rows, active_rows=active, state=state, width=16,
                         outer=[length, dimension], margin_bits=-float(optimum.fun),
                         x=math.exp(optimum.x[0]), z=math.exp(optimum.x[1]), rigorous=False)
        if best is None or candidate['margin_bits'] > best['margin_bits']:
            best = candidate
    return best


def mds_log_eta(rho, dimension=12):
    return max(math.log(fixed_support_multiplicity(16, dimension, w, 16))
               -(w-(16-dimension))*math.log(16)-w*math.log(rho)
               for w in range(17-dimension, 17))


def mds_sparse_screen(k_log, active, state, z, dimension=12):
    """Exact MDS support counts through four rows; requires input diagonal."""
    if not 1 <= active <= 4:
        raise ValueError('exact MDS source screen supports one through four rows')
    _, _, rows, threshold = geometry(k_log, 16, dimension)
    with np.errstate(divide='ignore'):
        factors = np.log(local_model(z, state)[:active+1])
    for w in range(len(factors)):
        factors[w] += math.log(math.comb(16, w))
    regions = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows//16
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, factors, active)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, active)
    for h in range(active+1):
        regions[h] -= log_choose(rows, h)
    matrices = [np.exp(regions[h])*16**(h-active*(16-dimension)/16) for h in range(active+1)]
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
        next_state /= norm
        scale += math.log(norm)
        current = next_state
    counts = np.array([fixed_support_multiplicity(16, dimension, w, 16)/16**(w-(16-dimension))
                       if w >= 17-dimension else 0. for w in range(17)])
    for axis in range(active):
        shape = [1]*(active+1)
        shape[axis] = 17
        current *= counts.reshape(shape)
    moment = scale+math.log(float(np.sum(current)))
    return -(log_choose(rows, active)-math.log(15)-threshold*math.log(z)+moment)/math.log(2)


def mds_sparse_search(k_log, active, state):
    k, n, rows, _ = geometry(k_log, 16, 12)
    def objective(tilt):
        return -mds_sparse_screen(k_log, active, state, math.exp(-math.exp(tilt)*active/n))
    optimum = minimize_scalar(objective, bounds=(-2, 12), method='bounded')
    return dict(k=k, n=n, rows=rows, outer=[16, 12], active_rows=active, state=state,
                margin_bits=-float(optimum.fun), z=math.exp(-math.exp(optimum.x)*active/n),
                input_diagonal='independent-uniform-nonzero-GF16', rigorous=False)


def mds_middle_screen(k_log, active, state, x, z, rho, dimension=12, aggregate=None):
    _, _, rows, threshold = geometry(k_log, 16, dimension)
    if aggregate is None:
        transfer = local_model(z, state)
        matrix = sum(math.comb(16, w)*(16*rho*x)**w*(1+x)**(16-w)*transfer[w] for w in range(17))
    else:
        # aggregate(v,z,s) = sum_w binom(16,w)*15**w*v**w*T_w(z).
        matrix = (1+x)**16*np.asarray(aggregate(16*rho*x/(15*(1+x)), z, state))
    value = (-15*log_choose(rows, active)-(16-dimension)*active*math.log(16)
             +active*mds_log_eta(rho, dimension)-math.log(15)-threshold*math.log(z)
             -16*active*math.log(x)+log_terminal_power(matrix, rows))
    return -value/math.log(2)


def mds_middle_search(k_log, active, state, aggregate=None):
    k, n, rows, _ = geometry(k_log, 16, 12)
    def objective(point):
        logx, logz, logrho = point
        if not -35 < logx < 10 or not -20 < logz < 0 or not -3 < logrho < 2:
            return 1e100
        return -mds_middle_screen(k_log, active, state, math.exp(logx), math.exp(logz), math.exp(logrho), aggregate=aggregate)
    best = None
    for initial_z in (.05, .5, .9, .99, .999):
        initial_x = active/(max(1, rows-active)*16)
        optimum = minimize(objective, [math.log(initial_x), math.log(initial_z), math.log(.9)], method='Nelder-Mead',
                           options={'maxiter': 600, 'xatol': 1e-8, 'fatol': 1e-6})
        candidate = dict(k=k, n=n, rows=rows, outer=[16, 12], active_rows=active, state=state,
                         margin_bits=-float(optimum.fun), x=math.exp(optimum.x[0]), z=math.exp(optimum.x[1]),
                         rho=math.exp(optimum.x[2]), input_diagonal='independent-uniform-nonzero-GF16', rigorous=False)
        if best is None or candidate['margin_bits'] > best['margin_bits']:
            best = candidate
    return best


def uniform_input_check(z, state):
    """Exact aggregate at nonzero-value fugacity one; diagnostic, no new bound."""
    total = ((1+15*z)/16)**16
    kernel = sum(math.comb(16, w)*birth_zero_count(w, state, 16)*z**w for w in range(17))/16**16
    returned = (total-kernel)/(16**state-1)
    exact = np.array([[kernel, total-kernel], [returned, total-returned]])
    transfer = local_model(z, state)
    envelope = sum(math.comb(16, w)*15**w*transfer[w] for w in range(17))/16**16
    return dict(z=z, state=state, exact=exact.tolist(), envelope=envelope.tolist(),
                exact_row_sum=total, envelope_nonzero_row_ratio=float(sum(envelope[1])/total),
                diagnostic_only=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--length', type=int, default=32)
    parser.add_argument('--dimension', type=int, default=24)
    parser.add_argument('--state', nargs='+', type=int, default=[3, 4, 5, 6])
    parser.add_argument('--active', nargs='+', type=int, default=[1])
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--middle', action='store_true')
    parser.add_argument('--mds', action='store_true', help='fixed GRS[16,12] with independent input diagonal')
    parser.add_argument('--exact-local', action='store_true', help='exact character aggregate in floating middle screens')
    args = parser.parse_args()
    aggregate = None
    if args.exact_local:
        if not args.middle or any(state > 3 for state in args.state):
            parser.error('--exact-local requires --middle and state sizes at most three')
        from exact_inner_character import aggregate
    for state in args.state:
        for active in args.active:
            if args.mds:
                method = mds_middle_search if args.middle else mds_sparse_search
                result = method(args.k_log, active, state, aggregate) if args.middle else method(args.k_log, active, state)
            else:
                method = middle_search if args.middle else sparse_search
                result = (method(args.k_log, active, state, args.length, args.dimension, aggregate)
                          if args.middle else method(args.k_log, active, state, args.length, args.dimension))
            if args.exact_local:
                result['local_model'] = 'exact-character-aggregate-floating-evaluation'
            print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
