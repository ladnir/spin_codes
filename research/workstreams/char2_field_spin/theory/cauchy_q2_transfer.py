"""Two-active-row regional transfer with exact per-row MDS support counts.

Floating screens only until a directed witness is explicitly replayed.
Distinct route slots, same-packet collisions, both nonempty source rows, and
one GLOBAL projective quotient are retained.  No iid row-support assumption.
"""

import argparse
import json
import math
from fractions import Fraction

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

from cauchy_firstmoment_screen import (exact_mds_projective_support, geometry,
                                      local_by_weight)
from reduced_distance_certificate import log_polynomial_product


WITNESSES = {16: Fraction(4979,5000), 18: Fraction(49947,50000),
             20: Fraction(24993,25000)}


def packet_matrices(y, return2_upper=None):
    matrices = local_by_weight(y)
    matrices[1][1][0] = y**16/255**4
    if return2_upper is not None:
        if return2_upper < 0:
            raise ValueError('negative return upper bound')
        matrices[2][1][0] = min(matrices[2][1][0], return2_upper)
    return matrices[:3]


def matrix_product(left, right):
    zero = left[0][0]*0
    return [[sum((left[i][k]*right[k][j] for k in range(2)), zero)
             for j in range(2)] for i in range(2)]


def polynomial_product(left, right, limit=2):
    zero = left[0][0][0]*0
    count = min(limit+1, len(left)+len(right)-1)
    result = [[[zero, zero], [zero, zero]] for _ in range(count)]
    for i, a in enumerate(left):
        for j, b in enumerate(right[:max(0, count-i)]):
            product = matrix_product(a, b)
            for r in range(2):
                for k in range(2):
                    result[i+j][r][k] += product[r][k]
    return result


def regional_coefficients(local, rows, width=16):
    """Positive generic arithmetic reference for c=0,1,2 active positions."""
    if rows % width or rows < max(width, 2) or len(local) < 3:
        raise ValueError('invalid regional geometry')
    zero, one = local[0][0][0]*0, local[0][0][0]*0+1
    factors = [[[entry*math.comb(width, c) for entry in row] for row in local[c]]
               for c in range(3)]
    result = [[[one, zero], [zero, one]]]
    power = rows//width
    while power:
        if power & 1:
            result = polynomial_product(result, factors)
        power >>= 1
        if power:
            factors = polynomial_product(factors, factors)
    return [[[entry/math.comb(rows, c) for entry in row] for row in result[c]]
            for c in range(3)]


def log_regions(y, rows, return2_upper=None):
    local = np.asarray(packet_matrices(y, return2_upper), dtype=float)
    with np.errstate(divide='ignore'):
        factors = np.log(local)
    for c in range(3):
        factors[c] += math.log(math.comb(16, c))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows//16
    while power:
        if power & 1:
            result = log_polynomial_product(result, factors, 2)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, 2)
    for c in range(3):
        result[c] -= math.log(math.comb(rows, c))
    return result


def two_row_polynomial(regions, length=32):
    """Positive generic reference; coefficient[a][b][state]."""
    zero, one = regions[0][0][0]*0, regions[0][0][0]*0+1
    current = [[[zero, zero] for _ in range(length+1)] for _ in range(length+1)]
    current[0][0][0] = one
    for region in range(length):
        following = [[[zero, zero] for _ in range(length+1)] for _ in range(length+1)]
        for a in range(region+1):
            for b in range(region+1):
                for da, db in ((0,0), (1,0), (0,1), (1,1)):
                    matrix = regions[da+db]
                    for end in range(2):
                        following[a+da][b+db][end] += sum(
                            (current[a][b][start]*matrix[start][end] for start in range(2)), zero)
        current = following
    return current


def log_two_row_polynomial(regions, length=32):
    current = np.full((length+1, length+1, 2), -math.inf)
    current[0, 0, 0] = 0.
    for region in range(length):
        following = np.full_like(current, -math.inf)
        for da, db in ((0,0), (1,0), (0,1), (1,1)):
            matrix = regions[da+db]
            for end in range(2):
                contribution = np.logaddexp(current[:region+1, :region+1, 0]+matrix[0, end],
                                             current[:region+1, :region+1, 1]+matrix[1, end])
                target = following[da:da+region+1, db:db+region+1, end]
                np.logaddexp(target, contribution, out=target)
        current = following
    return current


def bound(y, k_log=20, return2_upper=None, details=False):
    g = geometry(k_log)
    regions = log_regions(y, g['rows'], return2_upper)
    polynomial = log_two_row_polynomial(regions)
    # n_w counts projective words on ONE specified support.  Two nonzero
    # source rows contribute (q-1)*n_a*n_b after one global quotient.
    weights = np.full(33, -math.inf)
    for w in range(11, 33):
        weights[w] = math.log(exact_mds_projective_support(w))
    terms = polynomial+weights[:, None, None]+weights[None, :, None]
    moment = float(logsumexp(terms))
    value = (math.log(math.comb(g['rows'], 2))+math.log(255)
             -g['threshold']*math.log(y)+moment)
    margin = -value/math.log(2)
    if details:
        index = np.unravel_index(np.argmax(terms), terms.shape)
        return margin, dict(dominant_weights=[int(index[0]), int(index[1])],
                            dominant_terminal_state=int(index[2]),
                            collision_probability=15/(g['rows']-1))
    return margin


def search(k_log=20):
    g = geometry(k_log)
    def objective(logtilt):
        y = math.exp(-math.exp(logtilt)*2/g['n'])
        return -bound(y, k_log)
    choices = [minimize_scalar(objective, bounds=interval, method='bounded',
                               options={'xatol':1e-6}) for interval in ((0,8),(8,12))]
    best = min(choices, key=lambda result: result.fun)
    y = math.exp(-math.exp(best.x)*2/g['n'])
    margin, detail = bound(y, k_log, details=True)
    return dict(**g, k_log=k_log, active_rows=2, y=y, margin_bits=margin,
                local_return='originalR2; exactR1=y^16/255^4',
                rigorous=False, whole_code=False, **detail)


def directed_packets(mark):
    from cauchy_q1_certificate import interval, marginal
    from subfield_certificate import ONE, ZERO
    y, omitted = interval(mark), interval(1-mark)
    m2, r2 = ZERO, ZERO
    for i in range(3):
        for z in range(4):
            term = (math.comb(2, i)*math.comb(14, z)
                    *y**(16-i-z)*omitted**(i+z)/255**(i+z))
            m2 += term
            if i < 2:
                r2 += term/255**(2-i)
    return [[[ONE, ZERO], [ZERO, marginal(mark, 0)]],
            [[ZERO, y], [y**16/255**4, marginal(mark, 1)]],
            [[ZERO, y**2], [r2, m2]]]


def directed_expression(rows, mark, length=32, parity=10):
    from cauchy_q1_certificate import interval, projective_support_count
    from subfield_certificate import ZERO
    regions = regional_coefficients(directed_packets(mark), rows)
    polynomial = two_row_polynomial(regions, length)
    counts = [projective_support_count(w, parity) for w in range(length+1)]
    total = sum(((polynomial[a][b][0]+polynomial[a][b][1])*counts[a]*counts[b]
                 for a in range(parity+1, length+1)
                 for b in range(parity+1, length+1)), ZERO)
    threshold = (length*rows+9)//10-1
    return math.comb(rows, 2)*255*total/interval(mark)**threshold


def certify(k_log):
    from subfield_certificate import integer_margin
    if k_log not in WITNESSES:
        raise ValueError('only frozen K16/K18/K20 witnesses are certified')
    g, mark = geometry(k_log), WITNESSES[k_log]
    result = directed_expression(g['rows'], mark)
    if not 0 < result.lo <= result.hi < 1:
        raise ArithmeticError('invalid directed Q2 upper expression')
    margin = integer_margin(result.hi)
    if margin < 40:
        raise ArithmeticError('Q2 margin below40bits')
    return dict(**g, k_log=k_log, active_rows=2, y=str(mark),
                failure_upper=str(result.hi), class_margin_bits_lower=margin,
                local_return='originalR2; exactR1=y^16/255^4',
                arithmetic='80-digit outward-rounded Decimal intervals',
                rigorous=True, whole_code=False,
                scope='exactly two nonzero GRS outer rows; all remaining rows zero')


def log_all_regions(y, rows, active):
    """FLOAT conditional route matrices through count ``active``."""
    if not 1 <= active <= rows or rows % 16:
        raise ValueError('invalid route geometry')
    local = local_by_weight(y)
    local[1][1][0] = y**16/255**4
    with np.errstate(divide='ignore'):
        factors = np.log(np.asarray(local[:min(16, active)+1], dtype=float))
    for c in range(len(factors)):
        factors[c] += math.log(math.comb(16, c))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows//16
    while power:
        if power & 1:
            result = log_polynomial_product(result, factors, active)
        power >>= 1
        if power:
            factors = log_polynomial_product(factors, factors, active)
    for c in range(len(result)):
        result[c] -= math.log(math.comb(rows, c))
    return result


def log_union_source(active, logtheta, length=32, parity=10):
    """Exact source union-count recurrence, evaluated in floating logs.

    Result[j] is the marked product of projective per-support counts over
    A labeled rows whose support union equals one specified j-element set.
    Nonempty/minimum-weight constraints hold separately for EVERY row.
    """
    if active < 0:
        raise ValueError('negative row count')
    counts = [exact_mds_projective_support(w, parity) for w in range(length+1)]
    step = np.full((length+1, length+1), -math.inf)
    for target in range(length+1):
        for previous in range(target+1):
            added = target-previous
            terms = [math.log(math.comb(previous, overlap))+math.log(counts[added+overlap])
                     +(added+overlap)*logtheta
                     for overlap in range(previous+1) if counts[added+overlap]]
            if terms:
                step[target, previous] = math.log(math.comb(target, previous))+logsumexp(terms)
    vector = np.full(length+1, -math.inf)
    vector[0] = 0.
    while active:
        if active & 1:
            vector = logsumexp(step+vector[None, :], axis=1)
        active >>= 1
        if active:
            step = logsumexp(step[:, :, None]+step[None, :, :], axis=1)
    return vector


def log_region_union_polynomial(empty, nonempty, length=32):
    """e0 (empty+z nonempty)^length, with chronological matrix products."""
    current = np.full((length+1, 2), -math.inf)
    current[0, 0] = 0.
    for region in range(length):
        following = np.full_like(current, -math.inf)
        for delta, matrix in ((0, empty), (1, nonempty)):
            for end in range(2):
                contribution = np.logaddexp(current[:region+1, 0]+matrix[0, end],
                                             current[:region+1, 1]+matrix[1, end])
                target = following[delta:delta+region+1, end]
                np.logaddexp(target, contribution, out=target)
        current = following
    return logsumexp(current, axis=1)


def union_envelope_bound(active, y, logtheta, k_log=20, regions=None, details=False):
    """Proved fixed-word envelope; FLOAT evaluation is not a certificate."""
    g = geometry(k_log)
    if regions is None:
        regions = log_all_regions(y, g['rows'], active)
    scaled = regions[1:active+1]-np.arange(1, active+1)[:, None, None]*logtheta
    nonempty = np.max(scaled, axis=0)
    source = log_union_source(active, logtheta)
    routed = log_region_union_polynomial(regions[0], nonempty)
    terms = source+routed
    logarithm = (math.log(math.comb(g['rows'], active))+(active-1)*math.log(255)
                 -g['threshold']*math.log(y)+float(logsumexp(terms)))
    margin = -logarithm/math.log(2)
    if details:
        return margin, dict(dominant_union_size=int(np.argmax(terms)),
                            envelope_maximizers=(np.argmax(scaled, axis=0)+1).tolist())
    return margin


def union_envelope_search(active, k_log=20):
    """Bounded two-tilt exploratory search, preserving full region chronology."""
    g = geometry(k_log)
    cached_y, cached_regions = None, None

    def objective(point):
        nonlocal cached_y, cached_regions
        logtheta, logtilt = point
        if not -20 < logtheta < 8 or not -2 < logtilt < 14:
            return 1e100
        y = math.exp(-math.exp(logtilt)*active/g['n'])
        if y != cached_y:
            cached_y = y
            cached_regions = log_all_regions(y, g['rows'], active)
        return -union_envelope_bound(active, y, logtheta, k_log, cached_regions)

    starts = [(-6., 5.), (-3., 5.), (0., 5.)]
    choices = [minimize(objective, start, method='Nelder-Mead',
                        options={'maxiter':120, 'xatol':2e-5, 'fatol':1e-5})
               for start in starts]
    best = min(choices, key=lambda result: result.fun)
    logtheta, logtilt = best.x
    y = math.exp(-math.exp(logtilt)*active/g['n'])
    margin, detail = union_envelope_bound(active, y, logtheta, k_log, details=True)
    return dict(**g, active_rows=active, y=y, theta=math.exp(logtheta),
                margin_bits=margin, rigorous=False, whole_code=False,
                method='exact-source-union; regional-positive-count-monomial-envelope',
                **detail)


def word_moment_obstruction(k_log=20, active=16):
    """Exact lower bound on expected bad words, NOT on bad-code probability.

    Fix ``active`` row labels.  Route all their coordinates into aligned
    suffix windows in every region.  Zeroing the state after the first 31
    regions leaves at least 22*active-4*31 source dimensions.  This witness
    can contain words with fewer active rows and does not certify a failure
    lower bound for any one exact row-activity class.
    """
    g = geometry(k_log)
    window = (g['threshold']//(32*16))*16
    dimension = 22*active-4*31
    if active > window or dimension <= 0:
        raise ValueError('window or source dimension too small')
    probability = Fraction(math.comb(window, active), math.comb(g['rows'], active))**32
    expectation = probability*Fraction(256**dimension-1, 255)
    return dict(rows=g['rows'], active_rows_at_most=active, window=window,
                dimension=dimension, expected_projective_lines_lower=expectation)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, nargs='+', default=[16,18,20])
    parser.add_argument('--certify', action='store_true')
    parser.add_argument('--union-active', type=int, nargs='+')
    args = parser.parse_args()
    for k_log in args.k_log:
        if args.union_active:
            for active in args.union_active:
                print(json.dumps(union_envelope_search(active, k_log)), flush=True)
        else:
            print(json.dumps(certify(k_log) if args.certify else search(k_log)), flush=True)
