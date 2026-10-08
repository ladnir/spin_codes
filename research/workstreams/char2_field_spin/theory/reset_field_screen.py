"""Floating screens for fresh full-matrix feedback, with no state refresh.

Fix A as a [16,s] MDS evaluation map over GF16. Independently sample every
entry of C_i uniformly from GF16, including zero. The causal step is
Y=X+A Q; Q'=Q+C_i X. Every nonzero X makes Q' uniform and independent of
the current output. Zero inputs hold the exact state, including its weight.

Only the state weights of A Q are retained: zero and 17-s through 16.
The regional support average is exact for a fixed source word. Whole-code
screens use exact outer multiplicities for up to four active outer rows,
and the positive q^(w-8) envelope thereafter. Results are NOT certificates.
This separate proposal does not change the preserved GL4 construction.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
import json
import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from q1_batched_screen import fixed_support_multiplicity


def local_transfer(z, width=16, q=16, state=4):
    if not 0 < z <= 1 or not 1 <= state <= width <= q:
        raise ValueError('invalid reset-transfer parameters')
    weights = [0]+list(range(width-state+1, width+1))
    counts = [1]+[math.comb(width, v)*fixed_support_multiplicity(width, state, v, q)
                  for v in weights[1:]]
    assert sum(counts) == q**state
    reset = np.array(counts, dtype=float)/q**state
    result = np.zeros((width+1, state+1, state+1))
    f = z+(1-z)/(q-1)
    for category, v in enumerate(weights):
        result[0, category, category] = z**v
        for w in range(1, width+1):
            moment = sum(math.comb(width-v, a)*math.comb(v, w-a)
                         *z**(v-w+2*a)*f**(w-a)
                         for a in range(max(0, w-v), min(w, width-v)+1))
            result[w, category] = moment/math.comb(width, w)*reset
    return weights, result


def log_polynomial_product(left, right, limit):
    size = min(limit+1, len(left)+len(right)-1)
    states = left.shape[1]
    out = np.full((size, states, states), -math.inf)
    for degree, matrix in enumerate(left):
        count = min(len(right), size-degree)
        if count <= 0:
            break
        contribution = logsumexp(matrix[None, :, :, None]+right[:count, None, :, :], axis=2)
        out[degree:degree+count] = np.logaddexp(out[degree:degree+count], contribution)
    return out


def log_regions(z, maximum_h, rows, state=4):
    if rows % 16 or not 1 <= maximum_h < rows:
        raise ValueError('invalid regional geometry')
    _, local = local_transfer(z, state=state)
    with np.errstate(divide='ignore'):
        local = np.log(local[:maximum_h+1])
    for w in range(len(local)):
        local[w] += math.log(math.comb(16, w))
    out = np.full((1, state+1, state+1), -math.inf)
    np.fill_diagonal(out[0], 0)
    power = rows//16
    while power:
        if power & 1:
            out = log_polynomial_product(out, local, maximum_h)
        power >>= 1
        if power:
            local = log_polynomial_product(local, local, maximum_h)
    for h in range(maximum_h+1):
        out[h] -= math.log(math.comb(rows, h))
    return out


def screen_bound(k_log, active, z, state=4, distance=Fraction(1, 10), exact_source=None,
                 source_rho=None):
    if k_log not in (16, 18, 20):
        raise ValueError('supported geometry is K16/K18/K20')
    rows = (1 << k_log)//8
    desired = Fraction(distance)*rows*16
    threshold = (desired.numerator+desired.denominator-1)//desired.denominator-1
    if exact_source is None:
        exact_source = active <= 4
    regional = log_regions(z, active, rows, state)
    if exact_source:
        if not 1 <= active <= 4:
            raise ValueError('exact source is limited to four active rows')
        matrices = [np.exp(regional[h])*16**(h-active/2) for h in range(active+1)]
        masks = list(product((0, 1), repeat=active))
        current = np.zeros((1,)*active+(state+1,))
        current[(0,)*active+(0,)] = 1
        scale = 0.0
        for region in range(16):
            following = np.zeros((region+2,)*active+(state+1,))
            for mask in masks:
                target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
                following[target] += current @ matrices[sum(mask)]
            norm = float(np.max(following))
            following /= norm
            scale += math.log(norm)
            current = following
        counts = np.array([fixed_support_multiplicity(16, 8, w, 16)/16**(w-8)
                           if w >= 9 else 0 for w in range(17)])
        for axis in range(active):
            shape = [1]*(active+1)
            shape[axis] = 17
            current *= counts.reshape(shape)
        moment_log = scale+math.log(float(np.sum(current)))
    else:
        log_rho = math.log(source_rho) if source_rho is not None else 0.0
        log_eta = (max(math.log(fixed_support_multiplicity(16, 8, w, 16))
                       -(w-8)*math.log(16)-w*log_rho for w in range(9, 17))
                   if source_rho is not None else 0.0)
        source = np.array([math.log(math.comb(active, h))+(h-active/2)*math.log(16)+h*log_rho
                           for h in range(active+1)])
        matrix = logsumexp(regional+source[:, None, None], axis=0)
        vector = np.full(state+1, -math.inf)
        vector[0] = 0
        for _ in range(16):
            vector = logsumexp(vector[:, None]+matrix, axis=0)
        moment_log = float(logsumexp(vector))+active*log_eta
    logbound = math.log(math.comb(rows, active))-math.log(15)-threshold*math.log(z)+moment_log
    return -logbound/math.log(2)


def search(k_log, active, state=4, source_rho=None):
    n = 1 << (k_log+1)
    def objective(logtilt):
        z = math.exp(-math.exp(logtilt)*active/n)
        return -screen_bound(k_log, active, z, state, source_rho=source_rho)
    optimum = minimize_scalar(objective, bounds=(1, 8), method='bounded', options=dict(xatol=1e-6))
    return dict(k_log=k_log, active_rows=active, state=state,
                z=math.exp(-math.exp(optimum.x)*active/n), margin_bits=-float(optimum.fun),
                exact_source=active <= 4, rigorous=False,
                source_rho=source_rho,
                construction='fresh-unrestricted-feedback/no-state-refresh')


def self_check():
    """Independent GF4 enumeration checks the averaged emission/reset law."""
    def mul(a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            b >>= 1
            a <<= 1
            if a & 4:
                a ^= 7
        return value
    # A=(1,j), so the four state categories reduce to weights 0,3,4.
    states = list(product(range(4), repeat=2))
    images = [tuple(r[0] ^ mul(j, r[1]) for j in range(4)) for r in states]
    weights, actual = local_transfer(2/3, width=4, q=4, state=2)
    reset = [Fraction(sum(sum(v != 0 for v in image) == w for image in images), 16)
             for w in weights]
    # For every fixed nonzero input, all unrestricted feedback rows give
    # uniform output; translated incoming Q changes no count.
    for x in product(range(4), repeat=4):
        if not any(x):
            continue
        histogram = [0]*4
        for c in product(range(4), repeat=4):
            result = 0
            for a, b in zip(c, x):
                result ^= mul(a, b)
            histogram[result] += 1
        assert histogram == [64]*4
    z = Fraction(2, 3)
    for image in images:
        category = weights.index(sum(v != 0 for v in image))
        for w in range(1, 5):
            total = Fraction(0)
            for support in combinations(range(4), w):
                for values in product(range(1, 4), repeat=w):
                    y = list(image)
                    for lane, value in zip(support, values):
                        y[lane] ^= value
                    total += z**sum(v != 0 for v in y)
            moment = total/(math.comb(4, w)*3**w)
            for next_category, probability in enumerate(reset):
                assert math.isclose(actual[w, category, next_category], float(moment*probability),
                                    rel_tol=2e-14, abs_tol=1e-15)
        expected = [0.0]*len(weights)
        expected[category] = float(z**weights[category])
        assert np.allclose(actual[0, category], expected, rtol=2e-14, atol=1e-15)
    for state in (2, 4, 6):
        _, stochastic = local_transfer(1, state=state)
        assert np.allclose(np.sum(stochastic, axis=2), 1, rtol=2e-14, atol=1e-15)
    return 'GF4 exhaustive feedback/emission and GF16 row sums pass'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', nargs='+', type=int, default=[16])
    parser.add_argument('--rows', nargs='+', type=int, default=[1, 2, 3, 4, 8, 16])
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--source-rho', type=float)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        print(self_check(), flush=True)
    else:
        for k_log in args.k_log:
            for active in args.rows:
                print(json.dumps(search(k_log, active, args.state, args.source_rho)), flush=True)
