"""Reduced-distance investigation for the unchanged regional GF256 code.

At this checkpoint the Perron search is a floating-point screen only.
It uses the binomial chord bound and exact MDS source enumerator.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

from q1_batched_screen import fixed_support_multiplicity
from subfield_certificate import local_envelope
from state_screen import local as state_local


def local_model(z, width=16, q=256, state=2):
    if state == 2 and q == 256:
        return local_envelope(z, width, q)
    result = state_local(z, width, state, q)
    if q == 16 and state == 4 and width == 16:
        # For the monomial Vandermonde rows (1,j,j^2,j^3), C=A^T.
        # This specialization is not claimed for an arbitrary MDS basis.
        ret = 1/(q**state-1)
        joint = (13*z**15+z**16+2*z**14)/(16*(q**state-1))
        result[1, 1, 1] += ret*z**16-joint
        result[1, 1, 0] = joint
        for w in range(2, width+1):
            invertible_bound = sum(math.comb(width, h)*z**(width-h)*(1-z)**h
                                   /(q-1)**max(0, h-width+w) for h in range(width+1))/(q**state-1)
            result[w, 1, 0] = min(result[w, 1, 0], invertible_bound)
    return result


def log_polynomial_product(left, right, limit):
    size = min(limit+1, len(left)+len(right)-1)
    out = np.full((size, 2, 2), -math.inf)
    for degree, matrix in enumerate(left):
        count = min(len(right), size-degree)
        for row in range(2):
            for col in range(2):
                contribution = np.logaddexp(matrix[row, 0]+right[:count, 0, col],
                                            matrix[row, 1]+right[:count, 1, col])
                out[degree:degree+count, row, col] = np.logaddexp(
                    out[degree:degree+count, row, col], contribution)
    return out


def log_regions(z, maximum_h, rows=8192, width=16, q=256, state=2):
    """Logarithmic coefficients avoid the double underflow at a>=256."""
    local = local_model(z, width, q, state)[:maximum_h+1]
    with np.errstate(divide='ignore'):
        local = np.log(local)
    for w in range(len(local)):
        local[w] += math.log(math.comb(width, w))
    out = np.array([[[0, -math.inf], [-math.inf, 0]]])
    power = rows//width
    while power:
        if power & 1:
            out = log_polynomial_product(out, local, maximum_h)
        power >>= 1
        if power:
            local = log_polynomial_product(local, local, maximum_h)
    for h in range(maximum_h+1):
        out[h] -= math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
    return out


def coupled_log_bound(active_rows, z, distance=0.03, q=256, state=2):
    regions = log_regions(z, active_rows, q=q, state=state)
    source = np.array([math.log(math.comb(active_rows, h))+(h-active_rows/2)*math.log(q)
                       for h in range(active_rows+1)])
    matrix = logsumexp(regions+source[:, None, None], axis=0)
    vector = np.array([0, -math.inf])
    for region in range(16):
        vector = logsumexp(vector[:, None]+matrix, axis=0)
    threshold = math.ceil(distance*131072)-1
    choose = math.lgamma(8193)-math.lgamma(active_rows+1)-math.lgamma(8193-active_rows)
    return -(choose-math.log(q-1)-threshold*math.log(z)+float(logsumexp(vector)))/math.log(2)


def coupled_log_screen(active_rows, distance=0.03, search_steps=12, q=256, state=2):
    best = None
    for iz in range(search_steps+1):
        tilt = math.exp(1+6*iz/search_steps)
        z = math.exp(-tilt*active_rows/131072)
        bits = coupled_log_bound(active_rows, z, distance, q, state)
        candidate = dict(active_rows=active_rows, distance=distance, margin_bits=bits,
                         z=z, q=q, state=state, rigorous=False)
        if best is None or bits > best['margin_bits']:
            best = candidate
    return best


def exact_source_log_bound(active_rows, z, distance=0.1, q=16, state=4):
    """Exact MDS source multiplicities, for at most four selected rows."""
    if not 1 <= active_rows <= 4:
        raise ValueError('exact source-type screen is limited to four active rows')
    regions = np.exp(log_regions(z, active_rows, q=q, state=state))
    matrices = [regions[h]*q**(h-active_rows/2) for h in range(active_rows+1)]
    masks = [tuple((mask >> j) & 1 for j in range(active_rows)) for mask in range(1 << active_rows)]
    current = np.zeros((1,)*active_rows+(2,))
    current[(0,)*active_rows+(0,)] = 1
    cumulative = 0.0
    for region in range(16):
        next_state = np.zeros((region+2,)*active_rows+(2,))
        for mask in masks:
            target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
            next_state[target] += current @ matrices[sum(mask)]
        scale = float(np.max(next_state))
        next_state /= scale
        cumulative += math.log(scale)
        current = next_state
    factors = np.array([fixed_support_multiplicity(16, 8, w, q)/q**(w-8) if w >= 9 else 0
                        for w in range(17)])
    for axis in range(active_rows):
        shape = [1]*(active_rows+1)
        shape[axis] = 17
        current *= factors.reshape(shape)
    choose = math.lgamma(8193)-math.lgamma(active_rows+1)-math.lgamma(8193-active_rows)
    threshold = math.ceil(distance*131072)-1
    bound = choose-math.log(q-1)-threshold*math.log(z)+cumulative+math.log(float(np.sum(current)))
    return -bound/math.log(2)


def exact_source_screen(active_rows, distance=0.1, search_steps=24, q=16, state=4):
    best = None
    for iz in range(search_steps+1):
        tilt = math.exp(1+6*iz/search_steps)
        z = math.exp(-tilt*active_rows/131072)
        bits = exact_source_log_bound(active_rows, z, distance, q, state)
        candidate = dict(active_rows=active_rows, distance=distance, margin_bits=bits,
                         z=z, q=q, state=state, exact_source=True, rigorous=False)
        if best is None or bits > best['margin_bits']:
            best = candidate
    return best


def refined_screen(active_rows, distance=0.1, q=16, state=4, exact_source=False):
    fn = exact_source_log_bound if exact_source else coupled_log_bound
    def objective(logtilt):
        z = math.exp(-math.exp(logtilt)*active_rows/131072)
        return -fn(active_rows, z, distance, q, state)
    optimum = minimize_scalar(objective, bounds=(1, 7), method='bounded', options=dict(xatol=1e-6))
    return dict(active_rows=active_rows, distance=distance, q=q, state=state,
                z=math.exp(-math.exp(optimum.x)*active_rows/131072), margin_bits=-optimum.fun,
                exact_source=exact_source, rigorous=False)


def bernoulli_screen(active_rows, distance=0.1, q=16, state=4):
    """Use binom(a,h)/binom(L,h) <= (a/L)^h for every region."""
    rows, width = 8192, 16
    threshold = math.ceil(distance*131072)-1
    x = q*active_rows/rows
    choose = math.lgamma(rows+1)-math.lgamma(active_rows+1)-math.lgamma(rows-active_rows+1)

    def evaluate(logz, record=False):
        z = math.exp(logz)
        local = local_model(z, width, q, state)
        matrix = sum(math.comb(width, w)*x**w*local[w] for w in range(width+1))
        scale = 0
        power = rows
        vector = np.array([1.0, 0.0])
        while power:
            if power & 1:
                vector = vector @ matrix
                norm = float(np.max(vector))
                vector /= norm
                scale += math.log(norm)
            power >>= 1
            if power:
                matrix = matrix @ matrix
                norm = float(np.max(matrix))
                matrix /= norm
                # A matrix rescaling propagates into every remaining use.
                scale += math.log(norm)*power
        bound = (choose-math.log(q-1)-8*active_rows*math.log(q)
                 -threshold*logz+scale+math.log(float(np.sum(vector))))
        if record:
            return dict(active_rows=active_rows, margin_bits=-bound/math.log(2),
                        z=z, x=x, q=q, state=state, method='bernoulli', rigorous=False)
        return bound

    result = minimize_scalar(evaluate, bounds=(-10, -1e-8), method='bounded')
    return evaluate(result.x, True)


def chord_screen(active_rows, distance=0.03, rows=8192, width=16, q=256, state=2):
    n = 16*rows
    threshold = math.ceil(distance*n)-1
    logchoose = math.lgamma(rows+1)-math.lgamma(active_rows+1)-math.lgamma(rows-active_rows+1)
    logroot = logchoose/active_rows
    counts = np.array([math.log(math.comb(16, w)*fixed_support_multiplicity(16, 8, w, q))
                       for w in range(9, 17)])
    weights = np.arange(9, 17)

    def evaluate(parameters, record=False):
        logx, logz = parameters
        if not (-30 < logx < 5 and -15 < logz < 0):
            return 1e100
        x, z = math.exp(logx), math.exp(logz)
        local = local_model(z, width, q, state)
        matrix = sum(math.comb(width, w)*x**w*local[w] for w in range(width+1))
        a, b, c, d = matrix.ravel()
        gap = math.hypot(a-d, 2*math.sqrt(b*c))
        eigenvalue = (a+d+gap)/2
        # Stable right Perron vector (1,eta), avoiding lambda-a cancellation.
        eta = 2*c/(a-d+gap) if a >= d else (d-a+gap)/(2*b)
        logu = -logx-logroot
        terms = counts+weights*logu
        peak = float(np.max(terms))
        logp = peak+math.log(float(np.sum(np.exp(terms-peak))))
        bound = (logchoose-math.log(q-1)+max(0, -math.log(eta))
                 + (n//width)*math.log(eigenvalue)+active_rows*logp-threshold*logz)
        if record:
            return dict(active_rows=active_rows, distance=distance, margin_bits=-bound/math.log(2),
                        x=x, z=z, eta=eta, lambda_upper=eigenvalue, rigorous=False)
        return bound

    best = None
    for start_z in [0.1, 0.5, 0.8, 0.95, 0.99]:
        initial = [math.log(1/start_z-1), math.log(start_z)]
        result = minimize(evaluate, initial, method='Nelder-Mead',
                          options=dict(maxiter=250, xatol=1e-9, fatol=1e-7))
        candidate = evaluate(result.x, True)
        if isinstance(candidate, dict) and (best is None or candidate['margin_bits'] > best['margin_bits']):
            best = candidate
    return best


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', nargs='+', type=int, default=[16, 32, 64, 128, 256, 512, 908])
    parser.add_argument('--distance', type=float, default=0.03)
    parser.add_argument('--coupled', action='store_true')
    parser.add_argument('--steps', type=int, default=12)
    parser.add_argument('--q', type=int, default=256)
    parser.add_argument('--state', type=int, default=2)
    parser.add_argument('--exact-source', action='store_true')
    parser.add_argument('--bernoulli', action='store_true')
    args = parser.parse_args()
    for active in args.rows:
        answer = (bernoulli_screen(active, args.distance, args.q, args.state) if args.bernoulli
                  else exact_source_screen(active, args.distance, args.steps, args.q, args.state) if args.exact_source
                  else coupled_log_screen(active, args.distance, args.steps, args.q, args.state) if args.coupled
                  else chord_screen(active, args.distance, q=args.q, state=args.state))
        print(json.dumps(answer), flush=True)
