"""Complete support-cover relaxation for the regional state-four inner.

The cover and diagonal rank lemma are valid; the numerical envelope below
does NOT close. Optimization is exploratory, not directed certification.
This file performs mathematical checks only, never native timing.
"""

import argparse
from fractions import Fraction
import json
import math
from random import Random

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import gammaln, logsumexp

from region64_reference import (IDENTITY, global_step, phase_geometry,
                                translation_matrix, transpose, matvec)
from rate_design_coset import mul, power, xor_sum


def gf4_mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def diagonal_rank_checks():
    """Exhaust all nonzero ell and all offsets in a two-dimensional GF4 quotient."""
    count = 0
    for ell in ((a, b) for a in range(4) for b in range(4) if a or b):
        line = {tuple(gf4_mul(d, x) for x in ell) for d in range(4)}
        for offset in ((a, b) for a in range(4) for b in range(4)):
            for tau in (Fraction(1), Fraction(2), Fraction(3)):
                potential = Fraction(0)
                for diagonal in (1, 2, 3):
                    equation = tuple(gf4_mul(diagonal, x)^c
                                     for x, c in zip(ell, offset))
                    if equation == (0, 0):
                        potential += 1
                    elif equation not in line:
                        potential += 1/tau
                    # Otherwise imposing the equation kills required ell.
                assert potential/3 <= 1/tau
                count += 1
    return count


def cover_checks():
    rng = Random(20261005)
    _, phases = phase_geometry()
    count = 0
    for trial in range(80):
        state = (0, 0, 0, 0) if trial % 3 == 0 else tuple(rng.randrange(256) for _ in range(4))
        for step in range(8):
            values = tuple(rng.randrange(1, 256) if rng.random() < (trial % 9)/16 else 0
                           for _ in range(16))
            output, following = global_step(values, state, IDENTITY, phases[step % 4])
            left, right = bool(any(state)), bool(any(following))
            wi, wo = sum(bool(v) for v in values), sum(bool(v) for v in output)
            union = sum(bool(x or y) for x, y in zip(values, output))
            if not left:
                assert values == output
            else:
                assert union >= 13
            if not left and not right and wi:
                assert wi == wo and wi >= 5
            if left != right:
                assert wi and wo
            if not wo:
                assert following == state
            state = following
            count += 1
    return count


def emission_checks():
    points, phases = phase_geometry()
    evaluation_points, scales = [], []
    for beta in phases:
        for u in points:
            denominator = 1 ^ mul(beta, u)
            assert denominator
            z = mul(u, power(denominator, 254))
            scale = power(denominator, 3)
            row = matvec(transpose(translation_matrix(beta)),
                         (1, u, mul(u, u), power(u, 3)))
            assert row == tuple(mul(scale, x) for x in (1, z, mul(z, z), power(z, 3)))
            evaluation_points.append(z)
            scales.append(scale)
    assert evaluation_points.count(0) == 4
    nonzero = sorted(set(evaluation_points)-{0})
    assert len(nonzero) == 60
    # Cubics with root zero and any two other available points attain weight58.
    checked = 0
    for index, a in enumerate(nonzero):
        for b in nonzero[index+1:]:
            emitted = [mul(scale, mul(z, mul(z^a, z^b)))
                       for z, scale in zip(evaluation_points, scales)]
            assert sum(bool(x) for x in emitted) == 58
            checked += 1
    return checked


def repaired_packet_checks():
    """Check the separate anti-diagonal-J repair and its joint packet identity."""
    points, phases = phase_geometry()
    labels = tuple(u^beta for beta in phases for u in points)
    vector = lambda z: (1, z, mul(z, z), power(z, 3))
    for z in labels:
        for w in labels:
            assert xor_sum(mul(a, b) for a, b in zip(reversed(vector(z)), vector(w))) == power(z^w, 3)
    rng, checked = Random(20261006), 0
    for beta in phases:
        columns = tuple(vector(u^beta) for u in points)
        for trial in range(32):
            state = tuple(rng.randrange(256) for _ in range(4))
            emitted = tuple(xor_sum(mul(a, b) for a, b in zip(reversed(c), state)) for c in columns)
            values = (emitted if trial == 0 else
                      tuple(rng.randrange(256) if rng.randrange(4) == 0 else 0 for _ in range(16)))
            output = tuple(x^a for x, a in zip(values, emitted))
            following = tuple(state[r]^xor_sum(mul(c[r], x) for c, x in zip(columns, values))
                              for r in range(4))
            # F_J S plus fresh nonzero-output input contributions.
            factored = tuple(state[r]^xor_sum(mul(columns[j][r], emitted[j]^values[j])
                                for j in range(16) if output[j]) for r in range(4))
            assert following == factored
            if not any(output):
                assert following == state
            checked += 1
    return checked


def gf4_rank(matrix):
    work = [list(row) for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((j for j in range(rank, len(work)) if work[j][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = gf4_mul(work[rank][column], work[rank][column])
        work[rank] = [gf4_mul(inverse, a) for a in work[rank]]
        for j in range(len(work)):
            if j != rank and work[j][column]:
                value = work[j][column]
                work[j] = [a^gf4_mul(value, b) for a, b in zip(work[j], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def one_output_rank_checks():
    """Exhaust GF4 dependent forms and all still-fresh diagonals for v=1."""
    from itertools import product
    rng, checked = Random(20261007), 0
    orientations = [[list(row) for row in IDENTITY]]
    while len(orientations) < 3:
        candidate = [[rng.randrange(4) for _ in range(4)] for _ in range(4)]
        if gf4_rank(candidate) == 4:
            orientations.append(candidate)
    for orientation in orientations:
        for rank in range(1, 5):
            state = [row[:rank] for row in orientation]
            for label in range(4):
                c = (1, label, gf4_mul(label, label), gf4_mul(gf4_mul(label, label), label))
                a = tuple(reversed(c))
                functional = [xor_sum(gf4_mul(a[i], state[i][j]) for i in range(4))
                              for j in range(rank)]
                transformed = [[state[i][j]^gf4_mul(c[i], functional[j])
                                for j in range(rank)] for i in range(4)]
                assert gf4_rank(transformed) == rank
                enlarged_rank = gf4_rank([row+[c[i]] for i, row in enumerate(state)])
                for diagonal in (1, 2, 3):
                    independent = [row+[gf4_mul(diagonal, c[i])]
                                   for i, row in enumerate(transformed)]
                    assert gf4_rank(independent) == enlarged_rank
                for ell in product(range(4), repeat=rank):
                    if not any(ell):
                        continue
                    ranks = [gf4_rank([[transformed[i][j]^gf4_mul(diagonal, gf4_mul(c[i], ell[j]))
                                        for j in range(rank)] for i in range(4)])
                             for diagonal in (1, 2, 3)]
                    assert min(ranks) >= rank-1
                    assert sum(value < rank for value in ranks) <= 1
                    if enlarged_rank > rank:
                        assert min(ranks) == rank
                    checked += 1
    return checked


LOG_BINOMIAL = np.array([gammaln(17)-gammaln(w+1)-gammaln(17-w) for w in range(17)])
RETURN_TYPES = []
for union in range(13, 17):
    for only_input in range(union+1):
        for only_output in range(union-only_input+1):
            both = union-only_input-only_output
            if only_input+both and only_output+both:
                RETURN_TYPES.append((LOG_BINOMIAL[union]+gammaln(union+1)
                    -gammaln(only_input+1)-gammaln(only_output+1)-gammaln(both+1),
                    only_input, only_input+both, only_output+both))
RETURN_TYPES = np.array(RETURN_TYPES)


def screen_objective(variables, alpha, delta=.1, length=32, parity=10):
    log_x, log_y, log_tau, log_tilt, log_eta = variables
    log_g = -log_tau  # Valid on 1 <= tau <= q-1, keeping required input forms alive.
    log_union = logsumexp([log_g+log_x, log_y, log_x+log_y])
    zero_zero = logsumexp([0]+[LOG_BINOMIAL[w]+w*(log_x+log_y) for w in range(5, 17)])
    zero_nonzero = logsumexp([LOG_BINOMIAL[w]+w*(log_x+log_y) for w in range(1, 17)])
    nonzero_nonzero = logsumexp([LOG_BINOMIAL[w]+w*log_union for w in range(13, 17)])
    nonzero_zero = logsumexp(RETURN_TYPES[:, 0]+RETURN_TYPES[:, 1]*log_g
                            +RETURN_TYPES[:, 2]*log_x+RETURN_TYPES[:, 3]*log_y)
    shift = max(zero_zero, zero_nonzero, nonzero_zero, nonzero_nonzero)
    a, b, c, d = [math.exp(v-shift) for v in
                  (zero_zero, zero_nonzero, nonzero_zero, nonzero_nonzero)]
    log_radius = shift+math.log((a+d+math.sqrt((a-d)**2+4*b*c))/2)
    log_row = logsumexp([length*math.log1p(-alpha),
        -parity*log_tau-(parity+1)*log_eta
        +length*np.logaddexp(math.log1p(-alpha),
                             log_tau+log_tilt+log_eta+math.log(alpha))])
    return (log_radius/16-alpha*log_x-delta*log_y
            -alpha*log_tilt+log_row/length)/math.log(2)


def screen(alpha):
    candidates = [minimize(screen_objective, [-1, -2, lt, 0, 0], args=(alpha,),
                    method='L-BFGS-B', bounds=[(-25, 10), (-20, 0),
                       (0, math.log(255)), (-15, 15), (0, 15)])
                  for lt in (0, 1, 2, 4)]
    best = min(candidates, key=lambda value: value.fun)
    return dict(alpha=alpha, optimized_log2_envelope_per_N=float(best.fun),
                logarithmic_parameters=best.x.tolist(), certified=False,
                interpretation='positive means this evaluated envelope does not prove decay')


EMPTY_RUN_PATTERNS = []
for root_zero in (0, 1):
    for a in range(4-root_zero):
        for b in range(4-root_zero-a):
            for c in range(4-root_zero-a-b):
                for d in range(4-root_zero-a-b-c):
                    counts = (a, b, c, d)
                    multiplicity = math.prod(math.comb(15, j) for j in counts)
                    prefixes = [16*j-root_zero*j-sum(counts[:j]) for j in range(4)]
                    EMPTY_RUN_PATTERNS.append((multiplicity,
                                               64-4*root_zero-sum(counts), *prefixes))
EMPTY_RUN_PATTERNS = np.array(EMPTY_RUN_PATTERNS)


def empty_run_series(z, y):
    """Count actual periodic emission support masks for a fixed-state empty run."""
    total = 1.
    for packets in (1, 2, 3):
        total += z**packets*sum(math.comb(15*packets, roots)
            * y**(16*packets-roots-root_zero*packets)
            for root_zero in (0, 1) for roots in range(4-root_zero))
    cycle = z**4*y**EMPTY_RUN_PATTERNS[:, 1]
    prefix = sum(z**j*y**EMPTY_RUN_PATTERNS[:, 2+j] for j in range(4))
    return total+float(sum(EMPTY_RUN_PATTERNS[:, 0]*cycle*prefix/(1-cycle)))


def run_screen_objective(variables, alpha):
    log_x, log_y, log_tau, log_tilt, log_eta = variables
    x, y, g = math.exp(log_x), math.exp(log_y), math.exp(-log_tau)
    zero_zero = sum(math.comb(16, w)*(x*y)**w for w in range(5, 17))
    zero_nonzero = math.expm1(16*math.log1p(x*y))
    union = g*x+y+x*y
    nonzero_nonzero = sum(math.comb(16, w)*(union**w-y**w) for w in range(13, 17))
    nonzero_zero = math.exp(logsumexp(RETURN_TYPES[:, 0]-RETURN_TYPES[:, 1]*log_tau
        +RETURN_TYPES[:, 2]*log_x+RETURN_TYPES[:, 3]*log_y))
    def pole(log_z):
        z = math.exp(log_z)
        empty = empty_run_series(z, y)
        a, b = z*zero_zero/(1-z), z*zero_nonzero*empty
        c, d = z*nonzero_zero/(1-z), z*nonzero_nonzero*empty
        return (a+d+math.sqrt((a-d)**2+4*b*c))/2-1
    log_z = (-1e-14 if pole(-1e-14) <= 0
             else brentq(pole, -700, -1e-14, xtol=1e-10))
    row = logsumexp([32*math.log1p(-alpha), -10*log_tau-11*log_eta
        +32*np.logaddexp(math.log1p(-alpha), log_tau+log_tilt+log_eta+math.log(alpha))])
    return (-log_z/16-alpha*log_x-.1*log_y-alpha*log_tilt+row/32)/math.log(2)


def run_screen(alpha):
    candidates = [minimize(run_screen_objective, [-2, -1, lt, 0, 0], args=(alpha,),
                    method='L-BFGS-B', bounds=[(-18, 5), (-10, -1e-5),
                    (0, math.log(255)), (-15, 15), (0, 15)], options={'maxiter': 250})
                  for lt in (1, 3, 5)]
    best = min(candidates, key=lambda value: value.fun)
    return dict(alpha=alpha, fixed_state_run_log2_envelope_per_N=float(best.fun),
                logarithmic_parameters=best.x.tolist(), certified=False,
                scope='asymptotic rate of the complete cover; boundary prefactors omitted')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--screen', action='store_true')
    parser.add_argument('--run-screen', action='store_true')
    args = parser.parse_args()
    if args.check:
        print(json.dumps(dict(diagonal_rank_cases=diagonal_rank_checks(),
                              packet_cover_cases=cover_checks(),
                              minimum_emission_words=emission_checks(),
                              repaired_packet_cases=repaired_packet_checks(),
                              one_output_rank_cases=one_output_rank_checks())), flush=True)
    if args.screen:
        for density in (.0001, .001, .01, .04, .1, .15, .2, .2399):
            print(json.dumps(screen(density)), flush=True)
    if args.run_screen:
        for density in (.001, .01, .1, .2399):
            print(json.dumps(run_screen(density)), flush=True)
