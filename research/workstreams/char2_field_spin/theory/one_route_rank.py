"""Rank-potential screens for one regional routing stage.

This is NOT a finite certificate. The 'sum' mode upper-bounds ambiguous
next-rank continuations by their sum. The 'max' mode is HEURISTIC after
source-type aggregation: those types can change the incoming quotient,
so max cannot in general be pulled outside their sum. Both use an MDS
source-support envelope, not a count of field-valued messages. A uniform global permutation
with GRS[16,8] has a separate late-activation obstruction and is not tested.
"""

import argparse
import json
import math
from functools import lru_cache

import numpy as np

from batched_transfer import multiply


def return_rank_lower(weight, canceled, inactive_zeros):
    r"""Lower-bound rank((I+A C)_(Z,S)) for |Z\S| <= 1.

    A row-kernel relation is determined by a two-vector annihilating all
    C columns in S\Z. Distinct columns are pairwise independent.
    """
    uncanceled = weight-canceled
    return max(0, canceled+inactive_zeros-max(0, 2-uncanceled))


@lru_cache(maxsize=None)
def return_patterns(width=16, q=256):
    """Exact counts (w, outputweight, diagonalrank, kappam, forced, count).

    kappa is the number of M-only state constraints before exposing D.
    Full rank kills the required nonzero incoming state and is omitted.
    For p<=1 uncanceled lanes, subset XOR aggregates all cancellation sets.
    """
    if width < 4 or width & (width-1) or width > q:
        raise ValueError('the lane labels must form a binary subspace of size at least four')
    modulus = {4: 0b111, 8: 0b1011, 16: 0b10011, 256: 0x11b}[q]
    products = [[multiply(a, b, modulus) for b in range(q)] for a in range(q)]
    result = {}

    def add(w, out, rankd, rankm, forced, count):
        if count:
            key = (w, out, rankd, rankm, forced)
            result[key] = result.get(key, 0)+count

    # With at least two uncanceled inputs, their C columns span F_q^2.
    for p in range(2, width+1):
        for r in range(width-p+1):
            for u in range(min(1, width-p-r)+1):
                count = math.comb(width, r)*math.comb(width-r, p)
                if u:
                    count *= width-r-p
                add(r+p, width-r-u, r+2, u, width-p-u, count)

    for p in [0, 1]:
        choices_p = [()] if p == 0 else [(j,) for j in range(width)]
        for uncanceled in choices_p:
            for zero_inactive in [()]+[(j,) for j in range(width) if j not in uncanceled]:
                u = len(zero_inactive)
                excluded = set(uncanceled+zero_inactive)
                available = [j for j in range(width) if j not in excluded]
                counts = np.zeros((len(available)+1, width), dtype=np.int64)
                counts[0, 0] = 1
                for index, j in enumerate(available):
                    for size in range(index+1, 0, -1):
                        for xor in range(width):
                            counts[size, xor ^ j] += counts[size-1, xor]
                for r in range(len(available)+1):
                    if r+p == 0:
                        continue
                    for xor in range(width):
                        count = int(counts[r, xor])
                        if not count:
                            continue
                        d = [[1 ^ (r & 1), xor], [xor, 1 ^ products[xor][xor]]]
                        if p:
                            j = uncanceled[0]
                            equations = [(products[j][d[0][0]] ^ d[1][0],
                                          products[j][d[0][1]] ^ d[1][1])]
                        else:
                            equations = [tuple(d[0]), tuple(d[1])]
                        if u:
                            equations.append((1, zero_inactive[0]))
                        nonzero = [row for row in equations if row != (0, 0)]
                        rankm = 0
                        if nonzero:
                            a, b = nonzero[0]
                            rankm = 1
                            if any(products[a][d2] ^ products[b][c2] for c2, d2 in nonzero[1:]):
                                rankm = 2
                            else:
                                # ker(a,b) consists of multiples of (b,a).
                                # Required nonzero A_j R exclude roots except P or U.
                                if any((b ^ products[j][a]) == 0 for j in available):
                                    continue
                                if p:
                                    j = uncanceled[0]
                                    input_value = products[d[0][0]][b] ^ products[d[0][1]][a]
                                    output_value = input_value ^ b ^ products[j][a]
                                    if input_value == 0 or output_value == 0:
                                        continue
                        if rankm < 2:
                            add(r+p, width-r-u, r+p, rankm, width-p-u, count)
    return tuple((*key, value) for key, value in sorted(result.items()))


def local_operator(z, tau, width=16, q=256, sharpen_returns=True):
    """L_w(v)=fixed[w]@v + ambiguous[w]*max(v[1],v[2]).

    Counts output supports and rank potential averaged over input masks of
    weight w. binom(width,w) is inserted separately. Sharpened return counts
    first impose M-only consistency equations, then use every independent D
    pivot; the older complete-state target envelope remains as a control.
    """
    if not (0 < z <= 1 and tau >= 1):
        raise ValueError('require 0 < z <= 1 and tau >= 1')
    if width < 4 or width & (width-1) or width > q:
        raise ValueError('the lane labels must form a binary subspace of size at least four')
    fixed = np.zeros((width+1, 3, 3))
    ambiguous = np.zeros((width+1, 3))
    inverse_tau = 1/tau
    g = 1/(q-1)+(1-1/(q-1))*inverse_tau
    b1 = inverse_tau+(1-inverse_tau)/(q*q-1)
    b2 = (inverse_tau**2+(inverse_tau-inverse_tau**2)/(q-1)
          +(1-inverse_tau)/(q*(q*q-1)))
    fixed[0, 0, 0] = 1
    fixed[0, 1, 1] = (1-width/(q+1))*z**width+width/(q+1)*z**(width-1)
    fixed[0, 2, 1] = width/tau*z**(width-1)
    fixed[0, 2, 2] = z**width
    for w in range(1, width+1):
        if w == 1:
            fixed[w, 0, 1] = z
        else:
            ambiguous[w, 0] = z**w
            if w >= 3:
                fixed[w, 0, 0] = g*g*z**w
        for canceled in range(w+1):
            uncanceled = w-canceled
            for u in range(min(1, width-w)+1):
                multiplicity = math.comb(w, canceled)*((width-w) if u else 1)
                output_weight = width-canceled-u
                monomial = multiplicity*z**output_weight
                for h in [1, 2]:
                    if h == 1:
                        alpha = 1/(q+1) if u else 1-(width-uncanceled)/(q+1)
                        intermediate_rank = 1
                    else:
                        alpha = inverse_tau if u else 1.0
                        intermediate_rank = 1 if u else 2
                    factor = monomial*alpha*g**canceled
                    if intermediate_rank+uncanceled == 1:
                        fixed[w, h, 1] += factor
                    else:
                        ambiguous[w, h] += factor
                    # A zero output block cannot return a nonzero state:
                    # CA=0 implies C Y=C X=R on a return.
                    if output_weight and not sharpen_returns:
                        rank = return_rank_lower(w, canceled, u)
                        return_factor = g**rank*b1
                        if h == 2:
                            return_factor = min(return_factor, b2)
                        fixed[w, h, 0] += monomial*return_factor
    if sharpen_returns:
        for w, out, rankd, rankm, forced, count in return_patterns(width, q):
            factor = count/math.comb(width, w)*z**out*g**rankd
            fixed[w, 1, 0] += factor*(1/(q+1) if rankm else 1-forced/(q+1))
            fixed[w, 2, 0] += factor*(inverse_tau if rankm else 1)
    return fixed, ambiguous


def polynomial_product(left, right, limit):
    size = min(limit+1, len(left)+len(right)-1)
    out = np.zeros((size, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                out[:, i, j] += np.convolve(left[:, i, k], right[:, k, j])[:size]
    return out


def source_weights(active_rows, tau):
    return np.array([math.log(math.comb(active_rows, h))+(h-active_rows/2)*math.log(tau)
                     for h in range(active_rows+1)])


def regional_matrices(z, tau, active_rows, rows=8192, width=16):
    fixed, ambiguous = local_operator(z, tau, width)
    fixed[:, :, 1] += ambiguous
    fixed[:, :, 2] += ambiguous
    p = active_rows/rows
    local = fixed[:min(width, active_rows)+1].copy()
    for w in range(len(local)):
        local[w] *= math.comb(width, w)*p**w*(1-p)**(width-w)
    out = np.eye(3)[None, :, :]
    power = rows//width
    while power:
        if power & 1:
            out = polynomial_product(out, local, active_rows)
        power >>= 1
        if power:
            local = polynomial_product(local, local, active_rows)
    for h in range(active_rows+1):
        logprob = (math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
                   + h*math.log(p)+(rows-h)*math.log1p(-p))
        out[h] /= math.exp(logprob)
    return out


def regional_sum_matrix(z, tau, active_rows, rows=8192, width=16):
    out = regional_matrices(z, tau, active_rows, rows, width)
    weights = source_weights(active_rows, tau)
    scale = float(np.max(weights))
    return np.sum(out*np.exp(weights-scale)[:, None, None], axis=0), scale


def exact_source_bound(active_rows, tau, z, k=2**16, outer_length=16,
                       outer_dimension=8, width=16, delta=0.1):
    """Exact MDS support weights; linear sum-envelope for rank continuations."""
    assert 1 <= active_rows <= 4
    rows = k//outer_dimension
    n = rows*outer_length
    threshold = math.ceil(delta*n)-1
    matrices = regional_matrices(z, tau, active_rows, rows, width)
    for h in range(active_rows+1):
        matrices[h] *= tau**(h-active_rows/2)
    masks = [tuple((mask >> j) & 1 for j in range(active_rows)) for mask in range(1 << active_rows)]
    # Forward transfer of support weights; each local rank continuation
    # is already replaced by a sound linear sum-envelope.
    current = np.zeros((1,)*active_rows+(3,))
    current[(0,)*active_rows+(0,)] = 1
    cumulative = 0.0
    for region in range(outer_length):
        next_state = np.zeros((region+2,)*active_rows+(3,))
        for mask in masks:
            target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
            next_state[target] += current @ matrices[sum(mask)]
        scale = float(np.max(next_state))
        next_state /= scale
        cumulative += math.log(scale)
        current = next_state
    valid = tuple(slice(outer_length-outer_dimension+1, None) for _ in range(active_rows))+(slice(None),)
    total = float(np.sum(current[valid]))
    logchoose = math.lgamma(rows+1)-math.lgamma(active_rows+1)-math.lgamma(rows-active_rows+1)
    return -(logchoose-math.log(tau)-threshold*math.log(z)+cumulative+math.log(total))/math.log(2)


def regional_max_apply(continuation, z, tau, active_rows, rows=8192, width=16):
    """Fixed-future-continuation recursion; inter-region source sums need care.

    Applying this after aggregating future source types is not a proved
    whole-code bound. class_bound(mode='max') is a diagnostic only.
    """
    fixed, ambiguous = local_operator(z, tau, width)
    p = active_rows/rows
    maximum_w = min(width, active_rows)
    for w in range(maximum_w+1):
        scalar = math.comb(width, w)*p**w*(1-p)**(width-w)
        fixed[w] *= scalar
        ambiguous[w] *= scalar
    current = np.zeros((active_rows+1, 3))
    current[0] = continuation
    log_scale = 0.0
    for step in range(rows//width):
        next_state = np.zeros_like(current)
        for w in range(maximum_w+1):
            subset = current[:active_rows+1-w]
            next_state[w:] += (subset @ fixed[w].T
                              + np.maximum(subset[:, 1], subset[:, 2])[:, None]*ambiguous[w])
        scale = float(np.max(next_state))
        if scale == 0:
            return np.zeros(3), -math.inf
        next_state /= scale
        log_scale += math.log(scale)
        current = next_state
    weights = source_weights(active_rows, tau)
    for h in range(active_rows+1):
        weights[h] -= (math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
                       + h*math.log(p)+(rows-h)*math.log1p(-p))
    scale = float(np.max(weights))
    return np.sum(current*np.exp(weights-scale)[:, None], axis=0), log_scale+scale


def class_bound(active_rows, tau, z, k=2**16, outer_length=16,
                outer_dimension=8, width=16, delta=0.1, mode='sum'):
    rows = k//outer_dimension
    n = rows*outer_length
    threshold = math.ceil(delta*n)-1
    assert rows % width == 0 and 0 < active_rows < rows
    assert outer_length == 2*outer_dimension
    vector = np.ones(3)
    total_scale = 0.0
    if mode == 'sum':
        matrix, scale = regional_sum_matrix(z, tau, active_rows, rows, width)
        for region in range(outer_length):
            vector = matrix @ vector
            normalization = float(np.max(vector))
            vector /= normalization
            total_scale += scale+math.log(normalization)
    elif mode == 'max':
        for region in range(outer_length):
            vector, scale = regional_max_apply(vector, z, tau, active_rows, rows, width)
            normalization = float(np.max(vector))
            vector /= normalization
            total_scale += scale+math.log(normalization)
    else:
        raise ValueError(mode)
    logchoose = math.lgamma(rows+1)-math.lgamma(active_rows+1)-math.lgamma(rows-active_rows+1)
    logbound = logchoose-math.log(tau)-threshold*math.log(z)+total_scale+math.log(float(vector[0]))
    return -logbound/math.log(2)


def screen(active_rows, width=16, outer_length=16, mode='sum', search_steps=24, exact_source=False):
    n = 2**17
    best = None
    for tau in [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        for iz in range(search_steps+1):
            tilt = math.exp(1+6*iz/search_steps)
            z = math.exp(-tilt*active_rows/n)
            if exact_source:
                bits = exact_source_bound(active_rows, tau, z, width=width,
                                          outer_length=outer_length, outer_dimension=outer_length//2)
            else:
                bits = class_bound(active_rows, tau, z, width=width,
                                   outer_length=outer_length, outer_dimension=outer_length//2, mode=mode)
            if best is None or bits > best['margin_bits']:
                best = dict(active_rows=active_rows, width=width, outer_length=outer_length,
                            mode=mode, margin_bits=bits, tau=tau, z=z,
                            bound_semantics='linear-envelope' if mode == 'sum' or exact_source else 'heuristic')
    return best


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=int, nargs='+', default=[1, 2, 4, 8, 16])
    parser.add_argument('--width', type=int, default=16)
    parser.add_argument('--outer', type=int, default=16)
    parser.add_argument('--steps', type=int, default=24)
    parser.add_argument('--mode', choices=['sum', 'max'], default='sum')
    parser.add_argument('--tau', type=float)
    parser.add_argument('--z', type=float)
    parser.add_argument('--exact-source', action='store_true')
    args = parser.parse_args()
    for a in args.rows:
        if args.tau is not None and args.z is not None:
            kwargs = dict(width=args.width, outer_length=args.outer, outer_dimension=args.outer//2)
            if args.exact_source:
                margin = exact_source_bound(a, args.tau, args.z, **kwargs)
            else:
                margin = class_bound(a, args.tau, args.z, mode=args.mode, **kwargs)
            answer = dict(active_rows=a, width=args.width, outer_length=args.outer, mode=args.mode,
                          tau=args.tau, z=args.z,
                          margin_bits=margin, exact_source=args.exact_source,
                          bound_semantics='linear-envelope' if args.mode == 'sum' or args.exact_source else 'heuristic')
        else:
            answer = screen(a, width=args.width, outer_length=args.outer, mode=args.mode,
                            search_steps=args.steps, exact_source=args.exact_source)
        print(json.dumps(answer), flush=True)
