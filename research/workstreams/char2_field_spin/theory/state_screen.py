"""Floating whole-code screens for uniform-GL(s,256) batched inners.

These are not certificates. Keep the exact law: A is a [t,s] MDS generator,
C=A^T, uniform nonzero input scalars, and uniform GL(s,256) pre-emission
refresh. A/C entries are in GF256, never full-width field multipliers.
"""

import argparse
from functools import lru_cache
import json
import math

import numpy as np

from q1_batched_screen import fixed_support_multiplicity
from subfield_certificate import polynomial_product


@lru_cache(maxsize=None)
def zero_patterns(width, state, q):
    return tuple(fixed_support_multiplicity(width, state, width-j, q)
                 for j in range(state))


@lru_cache(maxsize=None)
def birth_zero_count(weight, state, q):
    return sum((-1)**i * math.comb(weight, i) * q**max(0, weight-state-i)
               for i in range(weight+1))


def local(z, width=16, state=4, q=256):
    """Return entrywise marked-transition upper bounds for finite 0 < z <= 1.

    The return-event weight bounds use decreasing powers of z. Values outside
    this domain do not give the stated envelope. Results use floating arithmetic
    and are exploratory screens, not directed-rounding certificates.
    """
    if not math.isfinite(z) or not 0 < z <= 1:
        raise ValueError("the marked-transition envelope requires finite 0 < z <= 1")
    assert 1 <= state <= width <= q
    result = np.zeros((width+1, 2, 2))
    f = 1/(q-1) + (q-2)/(q-1)*z
    patterns = zero_patterns(width, state, q)
    denominator = q**state-1
    for w in range(width+1):
        p0 = birth_zero_count(w, state, q)/(q-1)**w
        g = 0.0
        for zeros, count in enumerate(patterns):
            for active_zeros in range(max(0, zeros-(width-w)), min(w, zeros)+1):
                inactive_zeros = zeros-active_zeros
                g += (count/denominator * math.comb(w, active_zeros)
                      * math.comb(width-w, inactive_zeros)
                      * z**(width-w-inactive_zeros+active_zeros)
                      * f**(w-active_zeros))
        ret = (1-p0)/denominator
        result[w] = [[z**w*p0, z**w*(1-p0)],
                     [ret*z**max(0, width-state+1-w), g-ret*z**width]]
    return result


def regions(z, active, rows, width, state):
    assert 0 < active < rows and rows % width == 0
    p = active/rows
    factors = local(z, width, state)
    for w in range(width+1):
        factors[w] *= math.comb(width, w)*p**w*(1-p)**(width-w)
    result = np.eye(2)[None, :, :]
    exponent = rows//width
    while exponent:
        if exponent & 1:
            result = polynomial_product(result, factors, active)
        exponent >>= 1
        if exponent:
            factors = polynomial_product(factors, factors, active)
    for h in range(active+1):
        logprob = (math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
                   + h*math.log(p)+(rows-h)*math.log1p(-p))
        result[h] /= math.exp(logprob)
    return result


def coupled(active, k=2**16, width=16, state=4, b=16, dim=8, delta=0.1, samples=50):
    rows = k//dim
    n = rows*b
    e = math.ceil(delta*n)-1
    logchoose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    best = None
    for iz in range(samples+1):
        tilt = math.exp(1+7*iz/samples)
        logz = -tilt*active/n
        z = math.exp(logz)
        r = regions(z, active, rows, width, state)
        logs = np.array([math.log(math.comb(active, h))+(h-active*(b-dim)/b)*math.log(256)
                         for h in range(active+1)])
        scale = float(np.max(logs))
        matrix = np.sum(r*np.exp(logs-scale)[:, None, None], axis=0)
        vector = np.array([1.0, 0.0])
        total = b*scale
        for _ in range(b):
            vector = vector @ matrix
            norm = float(np.max(vector))
            if norm == 0:
                break
            vector /= norm
            total += math.log(norm)
        if norm == 0:
            continue
        bits = -(logchoose-math.log(255)-e*logz+total+math.log(float(np.sum(vector))))/math.log(2)
        if best is None or bits > best['margin_bits']:
            best = dict(active_rows=active, state=state, width=width, outer=[b, dim],
                        k=k, delta=delta, z=z, margin_bits=bits, certificate=False)
    return best


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active', nargs='+', type=int, default=[1, 2, 3, 4, 8, 16, 32, 64, 128])
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--width', type=int, default=16)
    parser.add_argument('--outer', type=int, default=16)
    parser.add_argument('--dimension', type=int, default=8)
    parser.add_argument('--k-log', type=int, default=16)
    parser.add_argument('--distance', type=float, default=0.1)
    parser.add_argument('--samples', type=int, default=50)
    args = parser.parse_args()
    for active in args.active:
        print(json.dumps(coupled(active, 2**args.k_log, args.width, args.state,
                                args.outer, args.dimension, args.distance, args.samples)), flush=True)
