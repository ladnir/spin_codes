"""Exploratory GF16 screens with scaled feedback C = gamma A^T.

The construction is GRS[16,8], one independent individual-symbol permutation
in each of sixteen regions, independent nonzero GF16 input diagonals, and
uniform GL(s,16) refresh before emission. A evaluates monomials of degree
below s at all sixteen field elements. Only 1 <= s <= 15 is considered.

Algebraic checks are exact. Floating searches are not certificates. The
optional directed-rounding middle-class replay certifies that class only,
not a whole-code distance claim. This file does not alter preserved models.
"""

import argparse
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
import json
import math

import numpy as np

from gf16_sizes import geometry, source_eta, interval_fraction
from gf16_certificate import Interval, HIGH, LOW, upper_polynomial_product, upper_exact_source
from gf16_sizes import upper_tilted_source
from q1_batched_screen import fixed_support_multiplicity
from state_screen import birth_zero_count, zero_patterns, local as generic_local
from subfield_certificate import ZERO, integer_margin


def mul(a, b):
    """Independent polynomial arithmetic in GF(2)[X]/(X^4+X+1)."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 0x13
    return result


def rank(matrix):
    rows = [list(row) for row in matrix]
    if not rows:
        return 0
    pivot = 0
    for column in range(len(rows[0])):
        found = next((i for i in range(pivot, len(rows)) if rows[i][column]), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        inverse = next(x for x in range(1, 16) if mul(rows[pivot][column], x) == 1)
        rows[pivot] = [mul(x, inverse) for x in rows[pivot]]
        for i in range(len(rows)):
            if i != pivot and rows[i][column]:
                factor = rows[i][column]
                rows[i] = [a ^ mul(factor, b) for a, b in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def xor_sum(values):
    result = 0
    for value in values:
        result ^= value
    return result


@lru_cache(maxsize=None)
def matrices(state, gamma=2):
    if not isinstance(state, int) or not 1 <= state <= 15:
        raise ValueError('state must be an integer from one through fifteen; state16 is excluded')
    if not isinstance(gamma, int) or not 1 <= gamma < 16:
        raise ValueError('gamma must be a nonzero abstract GF16 element')
    a = []
    for point in range(16):
        row = [1]
        for _ in range(1, state):
            row.append(mul(row[-1], point))
        a.append(tuple(row))
    gram = tuple(tuple(xor_sum(mul(row[i], row[j]) for row in a)
                       for j in range(state)) for i in range(state))
    small = tuple(tuple(int(i == j) ^ mul(gamma, gram[i][j])
                        for j in range(state)) for i in range(state))
    b = tuple(tuple(int(i == j) ^ mul(gamma, xor_sum(mul(x, y) for x, y in zip(a[i], a[j])))
                    for j in range(16)) for i in range(16))
    return tuple(a), gram, small, b


@lru_cache(maxsize=None)
def algebra(state, gamma=2):
    _, gram, small, b = matrices(state, gamma)
    pattern = all(gram[i][j] == int(i+j == 15)
                  for i in range(state) for j in range(state))
    return dict(state=state, gamma=gamma, gram_anti_diagonal=pattern,
                state_feedback_rank=rank(small), output_feedback_rank=rank(b),
                singleton_histogram=dict(sorted(Counter(sum(row[j] != 0 for row in b)
                                                         for j in range(16)).items())))


@lru_cache(maxsize=None)
def singleton_weights(state, gamma=2):
    report = algebra(state, gamma)
    if not report['gram_anti_diagonal'] or report['output_feedback_rank'] != 16:
        raise ValueError('the invertible-output return bound is not valid for this parameter set')
    b = matrices(state, gamma)[3]
    return tuple(sum(row[j] != 0 for row in b) for j in range(16))


def local_model(z, state=10, gamma=2):
    """Entrywise upper envelope at finite 0 < z <= 1, with exact singleton return."""
    weights = singleton_weights(state, gamma)
    result = generic_local(z, 16, state, 16)
    denominator = 16**state-1
    joint = sum(z**w for w in weights)/(16*denominator)
    result[1, 1, 1] += z**16/denominator-joint
    result[1, 1, 0] = joint
    for w in range(2, 17):
        inverse = sum(math.comb(16, h)*z**(16-h)*(1-z)**h
                      /15**max(0, h-16+w) for h in range(17))/denominator
        result[w, 1, 0] = min(result[w, 1, 0], inverse)
    return result


def interval_local_envelope(z, state=10, gamma=2):
    """Directed-rounding enclosure of the same entrywise envelope."""
    weights = singleton_weights(state, gamma)
    if z.lo <= 0 or z.hi > 1:
        raise ValueError('require 0 < z <= 1')
    f = Interval.rational(1, 15)+Interval.rational(14, 15)*z
    zs, fs, cs = [z**i for i in range(17)], [f**i for i in range(17)], [(1-z)**i for i in range(17)]
    patterns, denominator = zero_patterns(16, state, 16), 16**state-1
    result = []
    for w in range(17):
        p0 = Interval.rational(birth_zero_count(w, state, 16), 15**w)
        g = ZERO
        for zeros, count in enumerate(patterns):
            for active_zeros in range(max(0, zeros-(16-w)), min(w, zeros)+1):
                inactive_zeros = zeros-active_zeros
                coefficient = count*math.comb(w, active_zeros)*math.comb(16-w, inactive_zeros)
                g += (Interval.rational(coefficient, denominator)
                      *zs[16-w-inactive_zeros+active_zeros]*fs[w-active_zeros])
        ret = (1-p0)/denominator
        jlo, jhi = ret*zs[16], ret*zs[max(0, 17-state-w)]
        if w == 1:
            jlo = jhi = sum((zs[v] for v in weights), ZERO)/(16*denominator)
        elif w >= 2:
            inverse = sum((Interval.rational(math.comb(16, h), denominator*15**max(0, h-16+w))
                           *zs[16-h]*cs[h] for h in range(17)), ZERO)
            jhi = Interval(min(jhi.lo, inverse.lo), min(jhi.hi, inverse.hi))
        result.append([[zs[w]*p0, zs[w]*(1-p0)], [jhi, g-jlo]])
    return result


def screen_bound(k_log, active, z, state=10, gamma=2, exact_source=None, rho=Fraction(871, 1000)):
    """Log-domain regional coefficient screen; source normalization matches gf16_sizes."""
    from reduced_distance_certificate import log_polynomial_product
    from scipy.special import logsumexp
    rows, _, threshold = geometry(k_log)
    if not 1 <= active < rows:
        raise ValueError('invalid active-row count')
    exact_source = active <= 4 if exact_source is None else exact_source
    with np.errstate(divide='ignore'):
        local = np.log(local_model(z, state, gamma)[:active+1])
    for w in range(len(local)):
        local[w] += math.log(math.comb(16, w))
    regions = np.array([[[0, -math.inf], [-math.inf, 0]]])
    exponent = rows//16
    while exponent:
        if exponent & 1:
            regions = log_polynomial_product(regions, local, active)
        exponent >>= 1
        if exponent:
            local = log_polynomial_product(local, local, active)
    for h in range(active+1):
        regions[h] -= math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
    if exact_source:
        if not 1 <= active <= 4:
            raise ValueError('exact source screen is limited to four rows')
        matrices_ = [np.exp(regions[h])*16**(h-active/2) for h in range(active+1)]
        masks = [tuple((mask >> j) & 1 for j in range(active)) for mask in range(1 << active)]
        current = np.zeros((1,)*active+(2,))
        current[(0,)*active+(0,)] = 1
        scale = 0.0
        for region in range(16):
            next_state = np.zeros((region+2,)*active+(2,))
            for mask in masks:
                target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
                next_state[target] += current @ matrices_[sum(mask)]
            norm = float(np.max(next_state))
            next_state /= norm
            scale += math.log(norm)
            current = next_state
        factors = np.array([fixed_support_multiplicity(16, 8, w, 16)/16**(w-8) if w >= 9 else 0
                            for w in range(17)])
        for axis in range(active):
            shape = [1]*(active+1)
            shape[axis] = 17
            current *= factors.reshape(shape)
        moment_log = scale+math.log(float(np.sum(current)))
    else:
        log_rho, log_eta = math.log(float(rho)), math.log(float(source_eta(rho)))
        source = np.array([math.log(math.comb(active, h))+(h-active/2)*math.log(16)+h*log_rho
                           for h in range(active+1)])
        matrix = logsumexp(regions+source[:, None, None], axis=0)
        vector = np.array([0, -math.inf])
        for _ in range(16):
            vector = logsumexp(vector[:, None]+matrix, axis=0)
        moment_log = float(logsumexp(vector))+active*log_eta
    choose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    return -(choose-math.log(15)-threshold*math.log(z)+moment_log)/math.log(2)


def search(k_log, active, state=10, gamma=2, rho=Fraction(871, 1000)):
    from scipy.optimize import minimize_scalar
    _, n, _ = geometry(k_log)
    def objective(logtilt):
        return -screen_bound(k_log, active, math.exp(-math.exp(logtilt)*active/n), state, gamma, rho=rho)
    optimum = minimize_scalar(objective, bounds=(1, 8), method='bounded', options=dict(xatol=1e-6))
    return dict(k_log=k_log, state=state, gamma=gamma, active_rows=active,
                z=math.exp(-math.exp(optimum.x)*active/n), margin_bits=-float(optimum.fun),
                source_rho=str(rho), method='regional-coefficient', rigorous=False)


def middle_screen(k_log, state, active, x, z, gamma=2, rho=Fraction(871, 1000)):
    from coupled_chernoff import log_terminal_power
    rows, _, threshold = geometry(k_log)
    local = local_model(z, state, gamma)
    matrix = sum(math.comb(16, w)*(16*float(rho)*x)**w*(1+x)**(16-w)*local[w] for w in range(17))
    choose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    value = (-15*choose-8*active*math.log(16)+active*math.log(float(source_eta(rho)))-math.log(15)
             -threshold*math.log(z)-16*active*math.log(x)+log_terminal_power(matrix, rows))
    return -value/math.log(2)


def middle_search(k_log, state, active, gamma=2, rho=Fraction(871, 1000)):
    from scipy.optimize import minimize
    rows, _, _ = geometry(k_log)
    def objective(point):
        logx, logz = point
        if not -35 < logx < 10 or not -20 < logz < 0:
            return 1e100
        return -middle_screen(k_log, state, active, math.exp(logx), math.exp(logz), gamma, rho)
    best = None
    for initial_z in [.5, .9, .99, .999]:
        initial_x = active/(max(1, rows-active)*16)
        optimum = minimize(objective, [math.log(initial_x), math.log(initial_z)], method='Nelder-Mead',
                           options=dict(maxiter=300, xatol=1e-9, fatol=1e-7))
        candidate = dict(k_log=k_log, state=state, gamma=gamma, active_rows=active,
                         x=math.exp(optimum.x[0]), z=math.exp(optimum.x[1]),
                         source_rho=str(rho), margin_bits=-float(optimum.fun),
                         method='coupled-coefficient-chernoff', rigorous=False)
        if best is None or candidate['margin_bits'] > best['margin_bits']:
            best = candidate
    return best


def middle_upper(k_log, state, active, x, z, gamma=2, rho=Fraction(871, 1000)):
    from coupled_chernoff import interval_bound
    rows, _, _ = geometry(k_log)
    rhoi = interval_fraction(rho)
    def local_provider(zi):
        return [[[entry*rhoi**w for entry in row] for row in matrix]
                for w, matrix in enumerate(interval_local_envelope(zi, state, gamma))]
    return interval_bound(active, x, z, local_provider, rows=rows)*interval_fraction(source_eta(rho))**active


def class_upper(k_log, active, z, state=10, gamma=2, rho=Fraction(871, 1000)):
    """Rational regional replay for one active-row class, not a full certificate."""
    z = Fraction(z)
    zi = interval_fraction(z)
    rows, _, threshold = geometry(k_log)
    local = interval_local_envelope(zi, state, gamma)
    factors = [[HIGH.multiply(entry.hi, Decimal(math.comb(16, w))) for row in matrix for entry in row]
               for w, matrix in enumerate(local[:active+1])]
    out = [[Decimal(1), Decimal(0), Decimal(0), Decimal(1)]]
    power = rows//16
    while power:
        if power & 1:
            out = upper_polynomial_product(out, factors, active)
        power >>= 1
        if power:
            factors = upper_polynomial_product(factors, factors, active)
    regions = [[HIGH.divide(entry, Decimal(math.comb(rows, h))) for entry in matrix]
               for h, matrix in enumerate(out)]
    moment = upper_exact_source(regions, active) if active <= 4 else upper_tilted_source(regions, active, rho)
    bound = HIGH.divide(HIGH.multiply(moment, Decimal(math.comb(rows, active))),
                        LOW.multiply(Decimal(15), (zi**threshold).lo))
    return dict(k_log=k_log, state=state, gamma=gamma, active_rows=active, z=str(z),
                source_rho=str(rho), failure_upper=str(bound), integer_bits=integer_margin(bound),
                rigorous=True, scope='one active-row class only')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', nargs='+', type=int, default=[8, 9, 10, 12])
    parser.add_argument('--gamma', type=int, default=2)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--active', nargs='+', type=int, default=[1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
    parser.add_argument('--mode', choices=['algebra', 'sparse', 'middle'], default='middle')
    args = parser.parse_args()
    for state in args.state:
        print(json.dumps(algebra(state, args.gamma)), flush=True)
        if args.mode != 'algebra':
            for active in args.active:
                result = (search(args.k_log, active, state, args.gamma) if args.mode == 'sparse'
                          else middle_search(args.k_log, state, active, args.gamma))
                print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
