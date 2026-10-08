"""Floating-point one-active-outer-row screen; NOT a finite certificate.

Exact local probabilities and exact MDS multiplicities feed a scalar
Chernoff search. Each active outer coordinate lands uniformly in its own
region. Independent nonzero GF(q) input multipliers and a uniform SL(2,q)
pre-output refresh are required. This does not model the native-field
small-coefficient / three-shear implementation.
"""

import argparse
import json
import math

from batched_transfer import multiply


def fixed_support_multiplicity(length, dimension, weight, q):
    return sum((-1)**j * math.comb(weight, j)
               * q**max(0, weight-j-(length-dimension))
               for j in range(weight+1))


def return_weights(width, q=256, modulus=0x11b):
    # On a return from one active lane j, R=(x,j*x).
    return [sum((int(l == j) ^ 1 ^ multiply(l, j, modulus)) != 0
                for l in range(width)) for j in range(width)]


def q1_moment(z, k, width, outer_length=8, outer_dimension=4, q=256):
    rows = k // outer_dimension
    assert k % outer_dimension == 0 and rows % width == 0
    steps = rows // width
    f = 1/(q-1) + (q-2)/(q-1)*z
    # No active input: a nonzero state never disappears.
    a = (1-width/(q+1))*z**width + width/(q+1)*z**(width-1)
    region_survival = a**steps
    geometric_mean = (1-region_survival)/(steps*(1-a)) if a != 1 else 1.0
    g1 = ((1-width/(q+1))*z**(width-1)*f
          + z**width/(q+1)
          + (width-1)/(q+1)*z**(width-2)*f)
    j1 = sum(z**w for w in return_weights(width, q))/(width*(q*q-1))
    # One region, averaged over all positions of its single active symbol.
    h01 = z*geometric_mean
    h10 = j1*geometric_mean
    h11 = (g1-j1)*a**(steps-1)
    states = [(1.0, 0.0)] + [(0.0, 0.0)]*outer_length
    for region in range(outer_length):
        next_states = [(0.0, 0.0)]*(outer_length+1)
        for weight in range(region+1):
            v0, v1 = states[weight]
            x0, x1 = next_states[weight]
            next_states[weight] = (x0+v0, x1+v1*region_survival)
            x0, x1 = next_states[weight+1]
            next_states[weight+1] = (x0+v1*h10, x1+v0*h01+v1*h11)
        states = next_states
    # Count projective source lines, since a nonzero scalar does not
    # change a realized output support. The support choices are already
    # summed by the state polynomial, so no extra binomial factor occurs.
    terms = {}
    for weight in range(outer_length-outer_dimension+1, outer_length+1):
        lines = fixed_support_multiplicity(outer_length, outer_dimension, weight, q)//(q-1)
        terms[weight] = rows * lines * sum(states[weight])
    return sum(terms.values()), terms


def screen(k=2**16, width=16, outer_length=8, outer_dimension=4, delta=0.1):
    n = k*outer_length//outer_dimension
    e = math.ceil(delta*n)-1
    best = None
    # z=e^(-a/N), useful tilts from roughly one to a few hundred
    # expected active-region lifetimes. This is only a numerical screen.
    for index in range(1001):
        a = math.exp(-2 + 14*index/1000)
        log_z = -a/n
        z = math.exp(log_z)
        moment, terms = q1_moment(z, k, width, outer_length, outer_dimension)
        if moment == 0:
            continue
        bits = -(math.log(moment)-e*log_z)/math.log(2)
        if best is None or bits > best['margin_bits']:
            best = dict(k=k, n=n, width=width, outer=[outer_length, outer_dimension],
                        delta=delta, threshold=e, z=z, margin_bits=bits,
                        terms_at_best={str(w): -(math.log(value)-e*log_z)/math.log(2)
                                       for w, value in terms.items() if value})
    return best


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, default=16)
    parser.add_argument('--width', type=int, default=16)
    parser.add_argument('--outer-length', type=int, default=8)
    parser.add_argument('--outer-dimension', type=int, default=4)
    parser.add_argument('--distance', type=float, default=0.1)
    args = parser.parse_args()
    print(json.dumps(screen(2**args.k_log, args.width, args.outer_length,
                            args.outer_dimension, args.distance), indent=2))
