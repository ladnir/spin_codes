"""FLOAT fixed-word first-moment screen for the Cauchy diagonal-state inner.

The exact local helper accepts Fraction.  The fixed [b,k] MDS source has
per-exact-support multiplicity <=q**(w-(b-k)); this is a deterministic
count bound, not a second use of the random input diagonal.  No result
here is a directed certificate or a construction-failure lower bound.
"""

import argparse
from fractions import Fraction
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp


def geometry(k_log=20, length=32, dimension=22, period=64):
    rows = ((2**k_log+dimension*period-1)//(dimension*period))*period
    n = rows*length
    return dict(rows=rows, n=n, k=rows*dimension,
                threshold=(n+9)//10-1, length=length, dimension=dimension)


def local_by_weight(y, width=16, state=4, q=256, one_input_upgrade=False):
    """Two-state upper matrices, uniformly over fixed input values/state.

    T11 is the full marginal moment, so returns are conservatively counted
    again in T10.  Every return term omits i=w (no remaining overlap).
    """
    if not 0 < y <= 1 or not 1 <= state <= width or q < 2:
        raise ValueError('invalid local parameters')
    one = Fraction(1) if isinstance(y, Fraction) else 1.0
    zero = one-one
    result = []
    for w in range(width+1):
        marginal = zero
        returned = zero
        for i in range(w+1):
            for z in range(min(state-1, width-w)+1):
                term = (math.comb(w, i)*math.comb(width-w, z)
                        *y**(width-i-z)*(one-y)**(i+z)/(q-1)**(i+z))
                marginal += term
                if w and i < w:
                    returned += term/(q-1)**min(w-i, state)
        zero_zero = one if w == 0 else y**w/(q-1)**state if w > state else zero
        zero_nonzero = y**w if w else zero
        if one_input_upgrade and w == 1:
            # Every injection column has four nonzero entries.  A return
            # therefore either is impossible or fixes all state diagonals.
            # The returned output is a nonzero-state emission plus one
            # input coordinate, hence has weight at least width-state.
            returned = min(returned, y**(width-state)/(q-1)**state)
        result.append([[zero_zero, zero_nonzero], [returned, marginal]])
    return result


def aggregate(x, y, width=16, state=4, q=256, nonzero_only=False, one_input_upgrade=False):
    local = np.asarray(local_by_weight(y, width, state, q, one_input_upgrade), dtype=float)
    weights = np.asarray([math.comb(width, w)*((q-1)*x)**w*(1+x)**(width-w)
                          for w in range(width+1)])
    matrix = np.einsum('w,wij->ij', weights, local)
    if not nonzero_only:
        return matrix
    # Track whether any nonzero input was selected.  This removes the
    # allzero relaxed source by a positive construction, not subtraction.
    following = np.zeros((3, 3))
    following[0, 0] = weights[0]
    following[0, 1] = sum(weights[w]*local[w, 0, 0] for w in range(1, width+1))
    following[0, 2] = matrix[0, 1]
    following[1:, 1:] = matrix
    return following


def log_terminal_power(matrix, exponent, nonzero_only=False):
    vector = np.zeros(len(matrix))
    vector[0] = 1
    logarithm = 0.0
    scale = float(np.max(matrix))
    if not 0 < scale < math.inf:
        return math.inf
    matrix = matrix/scale
    logarithm += exponent*math.log(scale)
    while exponent:
        if exponent & 1:
            vector = vector@matrix
            scale_v = float(np.max(vector))
            if not 0 < scale_v < math.inf:
                return math.inf
            vector /= scale_v
            logarithm += math.log(scale_v)
        exponent >>= 1
        if exponent:
            matrix = matrix@matrix
            scale_m = float(np.max(matrix))
            if not 0 < scale_m < math.inf:
                return math.inf
            matrix /= scale_m
            logarithm += exponent*math.log(scale_m)
    terminal = float(sum(vector[1:] if nonzero_only else vector))
    # A finite positive x,y has positive nonzero-source mass.  Its floating
    # underflow is an unusable bound, never a zero failure probability.
    return logarithm+math.log(terminal) if terminal > 0 else math.inf


def log_choose(n, k):
    # Exact integer conversion avoids cancellation when k is very small.
    return math.log(math.comb(n, k))


def middle_bound(active, x, y, *, k_log=20, length=32, dimension=22,
                 width=16, state=4, q=256, nonzero_only=False, components=False,
                 one_input_upgrade=False):
    g = geometry(k_log, length, dimension)
    rows = g['rows']
    if not 1 <= active <= rows or x <= 0 or not 0 < y <= 1:
        raise ValueError('invalid geometry or tilts')
    matrix = aggregate(x, y, width, state, q, nonzero_only, one_input_upgrade)
    terms = dict(route=-(length-1)*log_choose(rows, active),
        source=-(length-dimension)*active*math.log(q)
               +length*active*math.log(q/(q-1)),
        projective=-math.log(q-1), output=-g['threshold']*math.log(y),
        selection=-length*active*math.log(x),
        packets=log_terminal_power(matrix, rows*length//width, nonzero_only))
    margin = -sum(terms.values())/math.log(2)
    if components:
        return margin, {name: value/math.log(2) for name, value in terms.items()}
    return margin


def search(active, *, nonzero_only=False, k_log=20, one_input_upgrade=False):
    g = geometry(k_log)
    def objective(point):
        logx, logdecay = point
        if not -35 < logx < 6 or not -18 < logdecay < 5:
            return 1e100
        x, y = math.exp(logx), math.exp(-math.exp(logdecay))
        return -middle_bound(active, x, y, k_log=k_log, nonzero_only=nonzero_only,
                             one_input_upgrade=one_input_upgrade)
    initial_x = active/(256*max(1, g['rows']-active))
    starts = [(math.log(initial_x), math.log(decay))
              for decay in (active/g['n'], 30*active/g['n'], .02, 1.0)]
    candidates = [minimize(objective, start, method='Nelder-Mead',
                           options={'maxiter':300, 'xatol':1e-9, 'fatol':1e-7})
                  for start in starts]
    best = min(candidates, key=lambda result: result.fun)
    x, y = math.exp(best.x[0]), math.exp(-math.exp(best.x[1]))
    margin, terms = middle_bound(active, x, y, k_log=k_log,
                                nonzero_only=nonzero_only, components=True,
                                one_input_upgrade=one_input_upgrade)
    return dict(**g, active=active, x=x, y=y, margin_bits=margin,
                log2_components=terms, nonzero_only=nonzero_only,
                rigorous=False, one_input_upgrade=one_input_upgrade,
                method='coupled-selection-coefficient')


def exact_mds_projective_support(weight, parity=10, q=256):
    """Number of projective MDS words on one particular exact support."""
    if weight <= parity:
        return 0
    return sum((-1)**j*math.comb(weight-1, j)*q**(weight-parity-1-j)
               for j in range(weight-parity))


def one_row_regional_bound(y, *, exact_source=False, k_log=20,
                           length=32, dimension=22, width=16, state=4, q=256,
                           one_input_upgrade=False):
    g = geometry(k_log, length, dimension)
    packets = g['rows']//width
    local = np.asarray(local_by_weight(y, width, state, q, one_input_upgrade), dtype=float)
    empty_active = local[0, 1, 1]
    if empty_active == 1:
        average = 1.0
    else:
        average = -math.expm1(packets*math.log(empty_active))/(packets*(1-empty_active))
    empty = np.diag([1.0, empty_active**packets])
    active = np.array([[0.0, local[1, 0, 1]*average],
                       [local[1, 1, 0]*average,
                        local[1, 1, 1]*empty_active**(packets-1)]])
    if not exact_source:
        moment = log_terminal_power(empty+(q-1)*active, length)
        value = (log_choose(g['rows'], 1)-(length-dimension)*math.log(q)
                 +length*math.log(q/(q-1))-math.log(q-1)
                 -g['threshold']*math.log(y)+moment)
    else:
        # Coefficient w sums the noncommuting regional products over all
        # support subsets of size w.  Scalar exact-support counts then apply.
        with np.errstate(divide='ignore'):
            log_empty, log_active = np.log(empty), np.log(active)
        current = np.full((length+1, 2), -math.inf)
        current[0, 0] = 0.0
        for region in range(length):
            following = np.full_like(current, -math.inf)
            for w in range(region+1):
                following[w] = np.logaddexp(following[w],
                    logsumexp(current[w, :, None]+log_empty, axis=0))
                following[w+1] = np.logaddexp(following[w+1],
                    logsumexp(current[w, :, None]+log_active, axis=0))
            current = following
        terms = [float(logsumexp(current[w]))
                 +math.log(exact_mds_projective_support(w, length-dimension, q))
                 for w in range(length-dimension+1, length+1)]
        value = log_choose(g['rows'], 1)-g['threshold']*math.log(y)+float(logsumexp(terms))
    return -value/math.log(2)


def one_row_regional_search(exact_source=False, k_log=20, one_input_upgrade=False):
    def objective(logdecay):
        return -one_row_regional_bound(math.exp(-math.exp(logdecay)),
                                      exact_source=exact_source, k_log=k_log,
                                      one_input_upgrade=one_input_upgrade)
    choices = [minimize_scalar(objective, bounds=interval, method='bounded')
               for interval in ((-16, -6), (-6, 0))]
    best = min(choices, key=lambda result: result.fun)
    return dict(active=1, exact_source=exact_source, y=math.exp(-math.exp(best.x)),
                margin_bits=-float(best.fun), rigorous=False,
                one_input_upgrade=one_input_upgrade,
                method='exact-one-row-regional-support-transfer')


def regional_bound(active, y, *, k_log=20, one_input_upgrade=False):
    """Exact per-region selected-count coefficients; dominated source only."""
    from reduced_distance_certificate import log_polynomial_product
    g = geometry(k_log)
    local = np.asarray(local_by_weight(y, one_input_upgrade=one_input_upgrade))
    factors = []
    for selected in range(min(16, active)+1):
        coefficient = math.comb(16, selected)*sum(
            math.comb(selected, w)*255.0**w*local[w] for w in range(selected+1))
        factors.append(coefficient)
    with np.errstate(divide='ignore'):
        factors = np.log(np.asarray(factors))
    regions = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = g['rows']//16
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, factors, active)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, active)
    region = regions[active]-log_choose(g['rows'], active)
    vector = np.array([0., -math.inf])
    for _ in range(32):
        vector = logsumexp(vector[:, None]+region, axis=0)
    value = (log_choose(g['rows'], active)-80*active*math.log(2)
             +32*active*math.log(256/255)-math.log(255)
             -g['threshold']*math.log(y)+float(logsumexp(vector)))
    return -value/math.log(2)


def regional_search(active, k_log=20, one_input_upgrade=False):
    g = geometry(k_log)
    def objective(logtilt):
        y = math.exp(-math.exp(logtilt)*active/g['n'])
        return -regional_bound(active, y, k_log=k_log,
                               one_input_upgrade=one_input_upgrade)
    choices = [minimize_scalar(objective, bounds=interval, method='bounded',
                               options={'xatol':1e-5}) for interval in ((0, 8), (8, 12))]
    best = min(choices, key=lambda result: result.fun)
    return dict(active=active, one_input_upgrade=one_input_upgrade,
                y=math.exp(-math.exp(best.x)*active/g['n']),
                margin_bits=-float(best.fun), rigorous=False,
                method='exact-regional-selected-count-dominated-source')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active', type=int, nargs='+', default=[1,16,64,256,1024,4096])
    parser.add_argument('--nonzero-only', action='store_true')
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--regional-one', action='store_true')
    parser.add_argument('--regional', action='store_true')
    parser.add_argument('--one-input-upgrade', action='store_true')
    args = parser.parse_args()
    if args.regional_one:
        for exact_source in (False, True):
            print(json.dumps(one_row_regional_search(exact_source, args.k_log,
                                                    args.one_input_upgrade)), flush=True)
    elif args.regional:
        for active in args.active:
            print(json.dumps(regional_search(active, args.k_log, args.one_input_upgrade)), flush=True)
    else:
        for active in args.active:
            print(json.dumps(search(active, nonzero_only=args.nonzero_only, k_log=args.k_log,
                                    one_input_upgrade=args.one_input_upgrade)), flush=True)
