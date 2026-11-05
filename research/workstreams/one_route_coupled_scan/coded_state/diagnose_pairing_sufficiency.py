"""FLOAT design diagnostic with one directed restricted quartet replay.

Joint support/zero-selection count with superregular local constraints.
Size-two groups model the current mixer; size-four groups are a stronger
comparison construction, not two successive pair layers. Never a whole-code
certificate. No files written.
"""

import argparse
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
from fractions import Fraction

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


LN2 = math.log(2)


def choose_bits(n, k):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / LN2


def choices(g, before, after, target, strict):
    """Necessary conditions from superregular A and A+delta I.

    A target coordinate has zero input. Polynomial coefficient j counts
    permitted selections of j zero outputs on the support intersection.
    """
    full = (1 << g) - 1
    active = before.bit_count()
    if active == 0 and after & target:
        return (0,) * (g + 1)
    if strict and active and (target & ~after).bit_count() >= active:
        return (0,) * (g + 1)
    if not strict and active == 1 and target & ~before & ~after:
        return (0,) * (g + 1)
    overlap = before & after
    counts = [0] * (g + 1)
    for zeros in range(1 << g):
        if zeros & ~overlap:
            continue
        forced = full & ~(before | after)
        if strict and active and (target & (forced | zeros)).bit_count() >= active:
            continue
        counts[zeros.bit_count()] += 1
    return tuple(counts)


def local_records(g, strict):
    groups = []
    for a in range(g + 1):
        before = (1 << a) - 1
        records = []
        for after, target in itertools.product(range(1 << g), repeat=2):
            coeff = choices(g, before, after, target, strict)
            if not any(coeff):
                continue
            coverage = target.bit_count() if a else (target & after).bit_count()
            records.append((target.bit_count(), coverage, after.bit_count(),
                            (before & ~after).bit_count(), *coeff))
        groups.append(np.array(records, dtype=float))
    return groups


def coefficient_power_log(coeff_logs, power, degree):
    """Coefficient of a positive polynomial power, with scale at each product."""
    # An exact positive change of variable centers the requested coefficient.
    # Without it, optimizer proposals can underflow before normalization.
    coeff_logs = np.array(coeff_logs, dtype=float)
    indices = np.arange(len(coeff_logs))
    tilt = 0.0
    if 0 < degree < power * (len(coeff_logs) - 1):
        lo, hi = -256.0, 256.0
        for _ in range(40):
            tilt = (lo + hi) / 2
            shifted = coeff_logs + indices * tilt
            terms = np.exp2(shifted - max(shifted))
            mean = float(terms @ indices / terms.sum())
            if mean < degree / power:
                lo = tilt
            else:
                hi = tilt
        coeff_logs = coeff_logs + indices * tilt
    shift = float(max(coeff_logs))
    base = np.exp2(np.array(coeff_logs) - shift)
    acc = np.array([1.0])
    scale = 0.0
    while power:
        if power & 1:
            acc = np.convolve(acc, base)[:degree + 1]
            peak = float(acc.max())
            acc /= peak
            scale += shift + math.log2(peak)
        power >>= 1
        if power:
            base = np.convolve(base, base)[:degree + 1]
            peak = float(base.max())
            base /= peak
            shift = 2 * shift + math.log2(peak)
    if len(acc) <= degree or acc[degree] <= 0:
        return float('inf')  # Reject underflow; never interpret it as a small bound.
    return math.log2(float(acc[degree])) + scale - degree * tilt


def objective(g, strict, target_rows=1974, swaps=32, quartet_fraction=None):
    n, w, L, M, regions = 512, 256, 4096, 8192, 256
    W = 46 * (1 << 21) // 100
    A, D = L * w, (L - 1) * swaps
    z = max(0, A + D - W)
    row_bits, initial_bits = choose_bits(M, target_rows), choose_bits(n, w)
    schedule = [(g, L - 1)]
    if quartet_fraction is not None:
        quartets = round((L - 1) * quartet_fraction)
        schedule = [(2,L - 1 - quartets),(4,quartets)]
    groups = [(size, count, local_records(size, strict)) for size, count in schedule if count]

    def evaluate(values):
        m, lp_odds, lzero, lnext, ldeath = values
        theta = m / 32
        lp = -float(np.logaddexp2(0, -lp_odds))
        lc = -float(np.logaddexp2(0, lp_odds))
        pi = row_bits + target_rows * lp + (M - target_rows) * lc
        continuation = 0.0
        for size, count, records in groups:
            terms = []
            powers = np.exp2(lzero * np.arange(size + 1))
            for a, data in enumerate(records):
                polynomial = data[:, 4:] @ powers
                weights = (data[:, 0] * lp + (size - data[:, 0]) * lc
                           - m * data[:, 1] + lnext * data[:, 2]
                           + ldeath * data[:, 3] + theta * np.log2(polynomial))
                terms.append(float(logsumexp(weights * LN2) / LN2) + math.log2(math.comb(size,a)))
            kernel = coefficient_power_log(terms, n // size, w) - initial_bits
            continuation += count * (kernel - w * lnext - swaps * ldeath)
        return (row_bits + initial_bits + w * lc - regions * pi + theta
                + m * (A - (1 << 20) + 128 * target_rows - z)
                - theta * z * lzero
                + continuation)
    return evaluate


def optimize(g, strict, target_rows=1974, swaps=32, quartet_fraction=None):
    f = objective(g, strict, target_rows, swaps, quartet_fraction)
    bounds = ((0.001,32),(-12,12),(-16,16),(-24,24),(-24,24))
    best = None
    for initial in ((4.,-1.7,-1.,0.,-3.),(16.,-1.7,-1.,0.,-3.)):
        result = minimize(lambda v: f(v) / (1 << 20), initial, method='L-BFGS-B', bounds=bounds,
                          options={'maxiter':180,'ftol':1e-13,'gtol':1e-8,'maxls':30})
        value = float(f(result.x))
        record = dict(group=g, strict=strict, target_rows=target_rows, swaps=swaps,
                      quartet_fraction=quartet_fraction, upper_bits=value, weights=result.x.tolist(),
                      optimizer_success=bool(result.success),
                      evaluations=int(result.nfev), float_only=True)
        if best is None or value < best['upper_bits']:
            best = record
    return best


def rational_quartet_replay():
    """Directed replay for ONE quartet profile, using rational positive markers.

    Uses a coefficient upper bound instead of the floating exact convolution.
    Algebraic roots are rounded upward by exact integer comparisons.
    """
    g, n, w, L, M, r = 4, 512, 256, 4096, 8192, 1974
    m, theta = 6, Fraction(3,16)
    p, zero, nxt, death = Fraction(8,11), Fraction(5,16), Fraction(5,18), Fraction(1,24)
    z = L*w + (L-1)*32 - 46*(1<<21)//100
    roots = {}

    def root_upper(value):
        if value in roots:
            return roots[value]
        cube = value**3
        precision = 40
        numerator = cube.numerator << (16*precision)
        denominator = cube.denominator
        lo, hi = 0, 1 << (precision + 2)
        assert hi**16 * denominator >= numerator
        while hi-lo > 1:
            mid = (lo+hi)//2
            if mid**16 * denominator >= numerator:
                hi = mid
            else:
                lo = mid
        result = Fraction(hi,1<<precision)
        assert result**16 >= cube
        assert Fraction(hi-1,1<<precision)**16 < cube
        roots[value] = result
        return result

    def logs(value):
        value = Fraction(value)
        assert value > 0
        power = 16
        d = (value.numerator**power).bit_length() - (value.denominator**power).bit_length()
        return Fraction(d-1,power), Fraction(d+1,power)

    def log_term_upper(coefficient, value):
        lo, hi = logs(value)
        return coefficient * (hi if coefficient >= 0 else lo)

    local = []
    for a in range(g+1):
        before = (1<<a)-1
        value = Fraction(0)
        for after, target in itertools.product(range(1<<g),repeat=2):
            coeff = choices(g,before,after,target,True)
            if not any(coeff):
                continue
            poly = sum(c*zero**j for j,c in enumerate(coeff))
            count = target.bit_count()
            covered = count if a else (target & after).bit_count()
            value += (p**count * (1-p)**(g-count) * Fraction(1,1 << (m*covered))
                      * nxt**after.bit_count() * death**(before & ~after).bit_count()
                      * root_upper(poly))
        local.append(value * math.comb(g,a))
    pi = math.comb(M,r) * p**r * (1-p)**(M-r)
    fixed = (log_term_upper(1,math.comb(M,r)) + log_term_upper(1,math.comb(n,w))
             + log_term_upper(w,1-p) + log_term_upper(-256,pi)
             + theta + m*(128*r-z) + log_term_upper(-theta*z,zero))
    # All-matrix rejection density for the comparison construction is <2.
    minor_degree = 2*sum(j*math.comb(g,j)**2 for j in range(1,g+1))
    assert minor_degree == 280
    assert Fraction((1<<21)//g * minor_degree,1<<32) < Fraction(1,2)
    best = None
    for incoming in (4,8,12,16,24,32):
        polynomial = sum(value*incoming**a for a,value in enumerate(local))
        kernel = polynomial**(n//g) / (math.comb(n,w)*incoming**w*nxt**w*death**32)
        upper = fixed + log_term_upper(L-1,kernel)
        record = dict(incoming_marker=incoming, log2_upper_rational=str(upper),
                      log2_upper_decimal=float(upper), directed=True,
                      scope='quartet; w=256; 32 swaps at each transition; 1974 inactive rows')
        print(record,flush=True)
        best = upper if best is None else min(best,upper)
    assert best < -41
    print(dict(replay_pass=True, whole_code=False, algebraic_roots_checked=len(roots)),flush=True)


def tests():
    checks = 0
    # Exhaustive GF8 accepted pair maps, all vectors and target masks.
    def mul(a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & 8:
                a ^= 11
        return result
    maps = []
    delta = 2
    for a, b, c, d in itertools.product(range(1, 8), repeat=4):
        if not (a ^ delta) or not (d ^ delta):
            continue
        if not (mul(a, d) ^ mul(b, c)) or not (mul(a ^ delta, d ^ delta) ^ mul(b, c)):
            continue
        maps.append((a, b, c, d))
    # Bound the test size by testing all matrices with eight representative
    # state/input pairs, rather than sampling or enumerating large instances.
    for a, b, c, d in maps:
        for q0, q1, x0, x1 in ((1,0,0,0),(0,1,0,0),(1,1,0,0),(1,2,0,0),
                              (1,0,1,0),(1,1,0,1),(0,0,1,0),(0,0,0,0)):
            nxt = (mul(a,q0)^mul(b,q1)^x0, mul(c,q0)^mul(d,q1)^x1)
            out = (nxt[0]^mul(delta,q0), nxt[1]^mul(delta,q1))
            before = int(q0 != 0) + 2 * int(q1 != 0)
            after = int(nxt[0] != 0) + 2 * int(nxt[1] != 0)
            target = int(x0 == 0) + 2 * int(x1 == 0)
            overlap = before & after
            zero_count = sum(bool(overlap & (1 << i)) and out[i] == 0 for i in range(2))
            assert choices(2,before,after,target,True)[zero_count] > 0
            checks += 1
    for g in (2, 4):
        for p in (1, 2, 4):
            coeff = [math.comb(g, a) for a in range(g + 1)]
            for degree in range(g * p + 1):
                got = coefficient_power_log([math.log2(x) for x in coeff], p, degree)
                assert abs(got - math.log2(math.comb(g * p, degree))) < 1e-10
                checks += 1
    # Degree-zero incoming states can only create nonzero values on allowed inputs.
    for g in (2, 4):
        for after, target in itertools.product(range(1 << g), repeat=2):
            assert bool(any(choices(g,0,after,target,True))) == (after & target == 0)
            checks += 1
    # Exact random-partition coefficient identity, independent of float evaluation.
    for g in (2,4):
        n, w = 2*g, g
        for local in ([Fraction(a+1, a+2) for a in range(g+1)],
                      [Fraction(2**a) for a in range(g+1)]):
            averaged = Fraction(0)
            partitions = 0
            support = set(range(w))
            for rest in itertools.combinations(range(1,n),g-1):
                block = {0,*rest}
                a = len(block & support)
                averaged += local[a] * local[w-a]
                partitions += 1
            averaged /= partitions
            coeff = sum(math.comb(g,a)*math.comb(g,w-a)*local[a]*local[w-a]
                        for a in range(g+1)) / math.comb(n,w)
            assert averaged == coeff
            checks += 1
    print(dict(checks=checks, accepted_GF8_maps=len(maps)), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--group', type=int, choices=(2,4), default=2)
    parser.add_argument('--coarse', action='store_true')
    parser.add_argument('--tests-only', action='store_true')
    parser.add_argument('--screen', action='store_true')
    parser.add_argument('--mixed', action='store_true')
    parser.add_argument('--replay-quartet', action='store_true')
    args = parser.parse_args()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free >= 20 * (1 << 30)
    timer = threading.Timer(60, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.tests_only:
            tests()
            return
        if args.replay_quartet:
            rational_quartet_replay()
            return
        if args.mixed:
            for fraction in (0.25,0.5,0.75):
                print(optimize(2,True,quartet_fraction=fraction),flush=True)
        elif args.screen:
            for swaps in (16,32,64):
                for target_rows in (1024,1974,2560,4096):
                    record = optimize(args.group, not args.coarse, target_rows, swaps)
                    record.pop('weights')
                    print(record, flush=True)
        else:
            print(optimize(args.group, not args.coarse), flush=True)
        print(dict(seconds=time.monotonic()-started, whole_code=False),flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    main()
