"""Sparse support screen with periodic independent GL4 state orientation.

This is a CHANGED construction, not an optimization of the fixed map.
All reported exponents omit finite boundary/conditioning factors and are
exploratory. No native benchmarks are run. The exact local intersection
probabilities are separately checked over every GF2 subspace of dimension4.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from return_cover import screen as return_screen


def gaussian(n, k, q):
    if k < 0 or k > n:
        return 0
    result = Fraction(1)
    for i in range(k):
        result *= Fraction(q**(n-i)-1, q**(k-i)-1)
    assert result.denominator == 1
    return int(result)


def intersection_probability(rank, roots, following, q=256, state=4):
    """Pr[dim(U intersect fixed codimension-roots space)=following]."""
    h, r, z, s = following, rank, roots, state
    if not max(0, r-z) <= h <= min(r, s-z):
        return Fraction(0)
    return Fraction(gaussian(s-z, h, q) * gaussian(z, r-h, q)
                    * q**((s-z-h)*(r-h)), gaussian(s, r, q))


@lru_cache(maxsize=None)
def probabilities(q=256):
    return np.array([[[float(intersection_probability(r, z, h, q))
                       for h in range(5)] for z in range(4)]
                     for r in range(1, 5)])


def root_factors(tau, q=256):
    """Worst-rank alive contraction for each root count z=0..3."""
    contraction = np.zeros((4, 4))
    for r in range(1, 5):
        for h in range(1, r+1):
            contraction[r-1] += probabilities(q)[r-1, :, h] / tau**(r-h)
    return contraction.max(axis=0)


@lru_cache(maxsize=None)
def support_types():
    result = []
    for z in range(4):
        for i in range(17-z):
            for a in range(17-z-i):
                o = 16-z-i-a
                c = math.comb(16, z)*math.comb(16-z, i)*math.comb(16-z-i, a)
                result.append((math.log(c), z, i, a, i+a, o+a))
    return np.array(result).T


def log_packet(log_x, log_y, log_tau, refreshed, one_output_filter=True,
               refresh_law='gl'):
    """Two-state positive upper transfer; states mean zero/nonzero witness."""
    coeff, z, only_i, both, w, v = support_types()
    if refresh_law == 'diagonal':
        from diagonal_refresh_rank import root_factors as factors_for
    elif refresh_law == 'gl':
        factors_for = root_factors
    else:
        raise ValueError('unknown refresh law')
    factors = factors_for(math.exp(log_tau)) if refreshed else np.ones(4)
    terms = (coeff + w*log_x+v*log_y-only_i*log_tau
             + np.log(factors)[z.astype(int)])
    active = float(logsumexp(terms))
    valid = (w > 0) & (v > 0) & (both > 0)
    if one_output_filter:
        # At Q'=0 and wt(Y)=1, X=(I+eta AC)Y has weight16.
        # eta is outside GF16; all off-diagonal coefficients are nonzero.
        valid &= (v != 1) | (w == 16)
    returns = float(logsumexp((terms-np.minimum(both, 4)*log_tau)[valid]))
    xy = log_x+log_y
    zero = float(logsumexp([0]+[math.log(math.comb(16, k))+k*xy-4*log_tau
                               for k in range(5, 17)]))
    births = float(logsumexp([math.log(math.comb(16, k))+k*xy for k in range(1, 17)]))
    return np.array([[zero, births], [returns, active]])


def log_period_radius(point, period, refresh_law='gl'):
    log_x, log_y, log_tau = point[:3]
    # period=0 means no additional refresh, preserved eta reference.
    phases = [False] if period == 0 else [True]+[False]*(period-1)
    product = np.eye(2)
    scale = 0.
    for fresh in phases:
        logarithms = log_packet(log_x, log_y, log_tau, fresh, refresh_law=refresh_law)
        shift = float(logarithms.max())
        product = product @ np.exp(logarithms-shift)
        maximum = float(product.max())
        product /= maximum
        scale += shift+math.log(maximum)
    a, b, c, d = product.ravel()
    eigenvalue = (a+d+math.sqrt((a-d)**2+4*b*c))/2
    return (scale+math.log(eigenvalue))/len(phases)


def objective(point, alpha, period, length=32, parity=10, refresh_law='gl'):
    log_x, log_y, log_tau, log_tilt = point
    log_row = logsumexp([length*math.log1p(-alpha),
        -parity*log_tau+length*np.logaddexp(math.log1p(-alpha),
                                         log_tau+log_tilt+math.log(alpha))])
    return (log_period_radius(point, period, refresh_law)/16-alpha*(log_x+log_tilt)
            -.1*log_y+log_row/length)/math.log(2)


def search(alpha, period, length=32, parity=10, refresh_law='gl'):
    base = return_screen(alpha)['logarithmic_parameters'][:4]
    bounds = [(-30, 12), (-25, 0), (0, math.log(255)), (-20, 20)]
    starts = [base, [-5, -.02, math.log(255), 0], [-1, -1, math.log(255), 0]]
    candidates = [minimize(objective, start, args=(alpha, period, length, parity, refresh_law),
                           method='L-BFGS-B', bounds=bounds,
                           options={'maxiter': 250}) for start in starts]
    best = min(candidates, key=lambda result: result.fun)
    return dict(alpha=alpha, refresh_every_symbols=16*period if period else None,
                refresh_law=refresh_law,
                outer=[length, length-parity], log2_upper_per_N=float(best.fun),
                log_parameters=best.x.tolist(), rigorous=False,
                whole_code=False, finite_factors_included=False)


def exact_checks():
    count = 0
    for r in range(1, 5):
        spaces = set()
        for basis in combinations(range(1, 16), r):
            span = {0}
            for vector in basis:
                span |= {x ^ vector for x in tuple(span)}
            if len(span) == 2**r:
                spaces.add(tuple(sorted(span)))
        assert len(spaces) == gaussian(4, r, 2)
        for z in range(5):
            histogram = [0]*5
            for space in spaces:
                size = sum((vector & ((1 << z)-1)) == 0 for vector in space)
                h = size.bit_length()-1
                histogram[h] += 1
            for h in range(5):
                assert Fraction(histogram[h], len(spaces)) == intersection_probability(r, z, h, 2)
                count += 1
    # A strict correction to the tempting but false bound Gamma_z<=tau^-z.
    p = intersection_probability(2, 1, 2)
    assert (1-p)/255+p > Fraction(1, 255)
    return count


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--screen', action='store_true')
    parser.add_argument('--alpha', nargs='+', type=float, default=[.001, .01, .04, .1, .15, .18])
    parser.add_argument('--period', nargs='+', type=int, default=[0, 1, 4])
    parser.add_argument('--law', choices=['gl', 'diagonal'], default='gl')
    args = parser.parse_args()
    if args.check:
        print(json.dumps({'exact_intersection_cases': exact_checks()}), flush=True)
    if args.screen:
        for period in args.period:
            for alpha in args.alpha:
                print(json.dumps(search(alpha, period, refresh_law=args.law)), flush=True)
