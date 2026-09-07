"""Directed batched evaluation of every (occupancy, all-one count) case.

This evaluates the same V_d(h) as certify_constant_split, sharing the fold
across occupancies. Positive 4x4 powers use explicit outward binary64 sums,
with an exact binary exponent. No BLAS or floating logarithm enters a bound.
"""
import math
from fractions import Fraction as F

import numpy as np
from flint import arb

import certify_constant_split as single

density, core, base, sparse = single.density, single.core, single.base, single.sparse
POLICY = 'zero-excluded;all-one-count-exact;ordinary-nonconstant-bands-only-v1'


def ordinary_bands():
    bands = sparse.BANDS[:-1]
    core.require(sparse.BANDS[-1] == (256,) and len(bands) == 12 and
                 all(b and all(0 < w < 256 for w in b) for b in bands),
                 'Zero/all-one row leaked into ordinary envelope')
    return bands


def up(x):
    return np.nextafter(x, np.inf)


def folds(mantissas, exponents, left, right):
    """Yield all offsets h at each ordinary-row count d, not only Q-d."""
    current, powers = mantissas.copy(), exponents.copy()
    yield 0, current, powers
    for d in range(1, len(current)):
        common = np.maximum(powers[:-1], powers[1:])
        a = up(np.ldexp(current[:-1], (powers[:-1]-common)[:, None, None]))
        b = up(np.ldexp(current[1:], (powers[1:]-common)[:, None, None]))
        updated = up(up(left[0]*a)+up(right[0]*b))
        for x, y in zip(left[1:], right[1:]):
            np.maximum(updated, up(up(x*a)+up(y*b)), out=updated)
        maximum = updated.max(axis=(1, 2))
        core.require(bool((maximum > 0).all() and np.isfinite(maximum).all()), 'Invalid fold')
        _, shift = np.frexp(maximum)
        current = up(np.ldexp(updated, -shift[:, None, None]))
        powers = common+shift
        yield d, current, powers


def terminal_256(matrices, exponents):
    """Upper bounds on e_Z^T (matrix * 2^exponent)^256 1.

Every term is nonnegative. Round each multiply/add upward separately and
normalize by exact powers of two after each square. Int64 exponents are
ample for the supported instances; reject unreasonable inputs before use.
"""
    core.require(matrices.ndim == 3 and matrices.shape[1:] == (4, 4) and
                 exponents.shape == (len(matrices),) and
                 bool(np.isfinite(matrices).all() and (matrices >= 0).all() and
                      (matrices <= 2).all() and (np.abs(exponents) < 10**9).all()),
                 'Invalid normalized positive matrices')
    a, powers = matrices.copy(), exponents.copy()
    for _ in range(8):
        # Fixed-width 4-term dot products: never delegate rounding to BLAS.
        out = up(a[:, :, 0, None]*a[:, None, 0, :])
        out = up(out+up(a[:, :, 1, None]*a[:, None, 1, :]))
        out = up(out+up(a[:, :, 2, None]*a[:, None, 2, :]))
        out = up(out+up(a[:, :, 3, None]*a[:, None, 3, :]))
        maximum = out.max(axis=(1, 2))
        core.require(bool((maximum > 0).all() and np.isfinite(maximum).all()), 'Invalid square')
        _, shift = np.frexp(maximum)
        a = up(np.ldexp(out, -shift[:, None, None]))
        powers = 2*powers+shift
    total = up(up(up(a[:, 0, 0]+a[:, 0, 1])+a[:, 0, 2])+a[:, 0, 3])
    return total, powers


def dyadic_ceiling(mantissa, exponent, count):
    """Exact ceil(log2(count * binary64_mantissa * 2^exponent)), capped at -80."""
    core.require(math.isfinite(mantissa) and mantissa > 0 and count > 0, 'Invalid upper bound')
    numerator, denominator = float(mantissa).as_integer_ratio()
    n = numerator*count
    ceiling = n.bit_length()-int(n & (n-1) == 0)
    return max(-80, ceiling+int(exponent)-(denominator.bit_length()-1))


def evaluate(region, ps, spec, lo, hi, tilt):
    bands = ordinary_bands()
    core.require(len(region) == hi+1 and len(ps) == len(bands) and 1 <= lo <= hi <= spec['rows'],
                 'Wrong grid shape')
    probs = [base.decode(p) for p in ps]
    core.require(all(0 < p < 1 for p in probs), 'Invalid ordinary probability')
    costs = sparse.costs_for(bands, probs, core.inputs.caps_module.caps())
    roots = [((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probabilities = [arb(p.numerator)/p.denominator for p in probs]
    left = up(np.array([float((r*(1-p)).upper()) for r, p in zip(roots, probabilities)]))
    right = up(np.array([float((r*p).upper()) for r, p in zip(roots, probabilities)]))
    keep = sparse.hull.indices(left, right)
    correction = (spec['cutoff']*(arb(tilt)/10).exp()).exp().upper()
    man, exp = correction.man_exp()
    correction_exp = int(exp)+int(man).bit_length()
    correction_man = float(np.nextafter(float(correction*arb(2)**(-correction_exp)), np.inf))
    mantissas, exponents = density.initial(region)
    result = [[None]*(q+1) for q in range(lo, hi+1)]
    locations = [math.comb(spec['rows'], q) for q in range(hi+1)]
    ordinary_count = 1
    for d, matrices, powers in folds(mantissas, exponents, left[keep], right[keep]):
        first_h = max(0, lo-d)
        terminal, terminal_exp = terminal_256(matrices[first_h:], powers[first_h:])
        terminal = up(terminal*correction_man)
        for offset, (value, exponent) in enumerate(zip(terminal, terminal_exp)):
            h = first_h+offset
            q = d+h
            count = locations[q]*math.comb(q, h)*ordinary_count
            result[q-lo][h] = dyadic_ceiling(value, int(exponent)+correction_exp, count)
        ordinary_count *= len(bands)
    core.require(all(all(type(p) is int for p in row) for row in result), 'Unvisited case')
    return result


def merge(best, values):
    core.require(len(best) == len(values) and all(len(a) == len(b) for a, b in zip(best, values)),
                 'Case shape mismatch')
    for a, b in zip(best, values):
        for h, power in enumerate(b):
            a[h] = power if a[h] is None else min(a[h], power)


def bounds(powers, lo, hi):
    core.require(len(powers) == hi-lo+1 and all(len(row) == q+1 and
                 all(type(p) is int and -80 <= p <= 10**8 for p in row)
                 for q, row in zip(range(lo, hi+1), powers)), 'Incomplete or invalid case table')
    return [sum((F(2)**p for p in row), F(0)) for row in powers]
