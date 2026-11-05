"""Exact selector-capacity thresholds for the unchanged quartet construction.

For fixed eligible rows H, each region contains a uniform h-subset of its
M=8192 slots.  For a quartet with j eligible slots put W=(j-r)_+.
Negative association of disjoint slot sets gives

  E beta^sum W <= B(h,r,beta)^J,
  B = sum_j C(h,j)C(M-h,4-j) beta^(j-r)_+ / C(M,4),

with J=2^19 quartets over all regions.  There is NO Bernoulli quota cost.
For actual support profiles G=F-f, Cout=R-128h, every matching footprint
has sum W >= (4-r)*(G/4)-Cout.  The H-set union costs C(M,h) once.

An independent exact-quota coefficient bound uses positive rational x:
Q(x,beta)=sum_j C(4,j)x^j beta^(j-r)_+, and bounds the same H-unioned
tail by Q^J/(x^(256h)*C(M,h)^255*beta^z).  Rational x is selected by
20 integer derivative bisections; its optimality is not needed for
validity.  Thresholds take the better of these two bounds by default.

Each table row certifies one H-unioned tail <=2^-92.  Choosing the best
beta before stating that tail adds no beta union.  All 8191*3 rows cost
less than2^-77 together; there is no support/profile/death-count union.
The table does not assert that every remaining profile is realizable.

All certificate arithmetic uses integers and rational logarithm intervals;
there is no floating-point logarithm.  S=2^48,20 atanh terms are fixed.
The profile function uses actual-state statistics.  Its L argument is
retained for the common profile interface but does not enter this cut.
Complement reflection is available mathematically via
Z_(H,r)=256h-rJ+Z_(Hc,4-r), but is not used by this implementation.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb
from operator import index
import os
import shutil
import threading
import time


SC_M = 8192
SC_N = 1 << 21
SC_J = SC_N // 4
SC_REGIONS = 256
SC_SCALE = 1 << 48
SC_TERMS = 20
SC_TAIL_BITS = 92
SC_BETAS = (Fraction(9, 8), Fraction(5, 4), Fraction(11, 8),
            Fraction(3, 2), Fraction(2), Fraction(4), Fraction(8))
SC_DEFAULT_H = tuple(range(256, SC_M, 256))


def _sc_ceildiv(a, b):
    if b <= 0:
        raise ValueError('Positive divisor required')
    return -((-a) // b)


def _sc_choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def _sc_ln_mantissa(n, d, scale=SC_SCALE, terms=SC_TERMS):
    """Integer bounds for scale*ln(n/d), where1<=n/d<=2."""
    if not (0 < d <= n <= 2 * d and scale > 0 and terms > 0):
        raise ValueError('Require1<=mantissa<=2 and positive precision')
    if n == d:
        return 0, 0
    a, b = n - d, n + d
    a2, b2 = a * a, b * b
    apow, bpow = a, b
    lower = 0
    for j in range(terms):
        lower += (2 * scale * apow) // ((2 * j + 1) * bpow)
        apow *= a2
        bpow *= b2
    # Tail<=2*z^(2L+1)/((2L+1)*(1-z^2)); each floor loses<1.
    tail = _sc_ceildiv(2 * scale * apow * b2,
                       (2 * terms + 1) * bpow * (b2 - a2))
    return lower, lower + terms + tail


_SC_LN2 = _sc_ln_mantissa(2, 1)


@lru_cache(maxsize=131072)
def selector_log2_interval(n, d=1):
    """Return integer lo,hi with lo<=2^48*log2(n/d)<=hi."""
    n, d = index(n), index(d)
    if n <= 0 or d <= 0:
        raise ValueError('Positive rational required')
    exponent = n.bit_length() - d.bit_length()
    if exponent >= 0:
        if n < (d << exponent):
            exponent -= 1
    elif (n << (-exponent)) < d:
        exponent -= 1
    if exponent >= 0:
        mn, md = n, d << exponent
    else:
        mn, md = n << (-exponent), d
    if mn == md:
        value = exponent * SC_SCALE
        return value, value
    lower, upper = _sc_ln_mantissa(mn, md)
    return (exponent * SC_SCALE + SC_SCALE * lower // _SC_LN2[1],
            exponent * SC_SCALE + _sc_ceildiv(SC_SCALE * upper, _SC_LN2[0]))


@lru_cache(maxsize=8192)
def selector_entropy_bits(h, slots=SC_M):
    """Integer e with C(slots,h)<=2^e, without a large binomial log."""
    h, slots = index(h), index(slots)
    if not 0 <= h <= slots or slots <= 0:
        raise ValueError('Require0<=h<=slots')
    if h in (0, slots):
        return 0
    lp = selector_log2_interval(h, slots)[0]
    lq = selector_log2_interval(slots - h, slots)[0]
    return _sc_ceildiv(-h * lp - (slots - h) * lq, SC_SCALE)


def selector_quartet_mgf(h, r, beta, slots=SC_M):
    """Exact hypergeometric quartet mgf; beta is a rational greater than1."""
    h, r, slots = index(h), index(r), index(slots)
    beta = Fraction(beta)
    if slots < 4 or not 0 <= h <= slots or not 0 <= r <= 3 or beta <= 1:
        raise ValueError('Require slots>=4,0<=h<=slots,0<=r<=3,beta>1')
    a, b = beta.numerator, beta.denominator
    degree = 4 - r
    numerator = sum(_sc_choose(h, j) * _sc_choose(slots - h, 4 - j)
                    * a ** max(0, j - r) * b ** (degree - max(0, j - r))
                    for j in range(5))
    return Fraction(numerator, comb(slots, 4) * b ** degree)


def _sc_compact_log2_interval(n, d=1):
    """Directed rational log after64-bit prefix truncation of large inputs."""
    n, d = index(n), index(d)
    if n <= 0 or d <= 0:
        raise ValueError('Positive rational required')
    sn, sd = max(0, n.bit_length() - 64), max(0, d.bit_length() - 64)
    if not (sn or sd):
        return selector_log2_interval(n, d)
    hn, hd = n >> sn, d >> sd
    un, ud = hn + int(bool(sn)), hd + int(bool(sd))
    offset = (sn - sd) * SC_SCALE
    return (offset + selector_log2_interval(hn, ud)[0],
            offset + selector_log2_interval(un, hd)[1])


@lru_cache(maxsize=8192)
def selector_binomial_log_interval(h, slots=SC_M):
    """Directed log of an exact binomial, using only its64-bit prefix."""
    h, slots = index(h), index(slots)
    if not 0 <= h <= slots or slots < 1:
        raise ValueError('Invalid binomial dimensions')
    return _sc_compact_log2_interval(comb(slots, h))


def _sc_label_coefficients(r, beta):
    beta = Fraction(beta)
    if not 0 <= r <= 3 or beta <= 1:
        raise ValueError('Require0<=r<=3 and beta>1')
    a, b = beta.numerator, beta.denominator
    degree = 4 - r
    return tuple(comb(4, j) * a ** max(0, j - r)
                 * b ** (degree - max(0, j - r)) for j in range(5))


def _sc_label_numerators(n, d, coefficients):
    """Fixed degree-four positive kernel: numerator of Q and x*dQ/dx."""
    n2, d2 = n * n, d * d
    q0 = coefficients[0] * d2 * d2
    q1 = coefficients[1] * n * d * d2
    q2 = coefficients[2] * n2 * d2
    q3 = coefficients[3] * n2 * n * d
    q4 = coefficients[4] * n2 * n2
    return q0 + q1 + q2 + q3 + q4, q1 + 2 * q2 + 3 * q3 + 4 * q4


def selector_label_polynomial(r, beta, x):
    """Exact Q_r(x,beta)=sum_j C(4,j)x^j beta^(j-r)_+."""
    r, beta, x = index(r), Fraction(beta), Fraction(x)
    if x <= 0:
        raise ValueError('Positive x required')
    numerator, _ = _sc_label_numerators(x.numerator, x.denominator,
                                       _sc_label_coefficients(r, beta))
    return Fraction(numerator, x.denominator ** 4 * beta.denominator ** (4 - r))


def selector_quota_tilt(h, r, beta, steps=20, slots=SC_M):
    """Choose rational x near xQ'/Q=4h/M; every returned x is valid."""
    h, r, steps, slots = map(index, (h, r, steps, slots))
    beta = Fraction(beta)
    if not 0 < h < slots or not 1 <= r <= 3 or not 1 <= steps <= 40:
        raise ValueError('Invalid quota-tilt parameters')
    coefficients = _sc_label_coefficients(r, beta)
    a, b = beta.numerator, beta.denominator
    grid = 1 << steps
    low, high = 0, grid
    # x=x0*(1/beta+(1-1/beta)*m/grid), x0=h/(M-h).
    den = (slots - h) * a * grid
    while low < high:
        middle = (low + high) // 2
        num = h * (b * grid + (a - b) * middle)
        value, derivative = _sc_label_numerators(num, den, coefficients)
        if slots * derivative < 4 * h * value:
            low = middle + 1
        else:
            high = middle
    return Fraction(h * (b * grid + (a - b) * low), den)


def selector_maximum(h, r, slots=SC_M, regions=SC_REGIONS):
    """Exact deterministic maximum, by packing eligible slots in quartets."""
    h, r, slots, regions = map(index, (h, r, slots, regions))
    if slots % 4 or not 0 <= h <= slots or not 0 <= r <= 3 or regions < 1:
        raise ValueError('Invalid selector dimensions')
    return regions * ((4 - r) * (h // 4) + max(0, h % 4 - r))


@lru_cache(maxsize=32768)
def selector_threshold_row(h, r, betas=SC_BETAS, tail_bits=SC_TAIL_BITS,
                           exact_quota=True):
    """Certify the selected threshold, or the deterministic impossible tail."""
    h, r, tail_bits = map(index, (h, r, tail_bits))
    if not 1 <= h < SC_M or not 1 <= r <= 3 or tail_bits < 1:
        raise ValueError('Require1<=h<8192,1<=r<=3,tail_bits>=1')
    entropy = selector_entropy_bits(h)
    binomial_lower = selector_binomial_log_interval(h)[0] if exact_quota else None
    maximum = selector_maximum(h, r)
    best = dict(h=h, r=r, threshold=maximum + 1,
                kind='deterministically_impossible', beta=None,
                maximum=maximum, tail_bits=tail_bits,
                entropy_bits=entropy, certificate=True)
    for supplied in betas:
        beta = Fraction(supplied)
        mgf = selector_quartet_mgf(h, r, beta)
        upper = selector_log2_interval(mgf.numerator, mgf.denominator)[1]
        lower = selector_log2_interval(beta.numerator, beta.denominator)[0]
        if lower <= 0:
            raise ArithmeticError('The beta log lower bound must be positive')
        tail_numerator = SC_SCALE * entropy + SC_J * upper
        threshold = _sc_ceildiv(tail_numerator + SC_SCALE * tail_bits, lower)
        if threshold < best['threshold']:
            assert entropy * SC_SCALE + SC_J * upper - threshold * lower <= -tail_bits * SC_SCALE
            best = dict(h=h, r=r, threshold=threshold,
                        kind='hypergeometric_NA', beta=(beta.numerator, beta.denominator),
                        maximum=maximum, tail_bits=tail_bits,
                        entropy_bits=entropy, log_mgf_upper=upper,
                        log_beta_lower=lower, tail_numerator_upper=tail_numerator,
                        certificate=True)
        if exact_quota:
            x = selector_quota_tilt(h, r, beta)
            polynomial = selector_label_polynomial(r, beta, x)
            polynomial_upper = _sc_compact_log2_interval(polynomial.numerator, polynomial.denominator)[1]
            x_lower = selector_log2_interval(x.numerator, x.denominator)[0]
            tail_numerator = (SC_J * polynomial_upper - SC_REGIONS * h * x_lower
                              - (SC_REGIONS - 1) * binomial_lower)
            threshold = _sc_ceildiv(tail_numerator + SC_SCALE * tail_bits, lower)
            if threshold < best['threshold']:
                assert tail_numerator - threshold * lower <= -tail_bits * SC_SCALE
                best = dict(h=h, r=r, threshold=threshold, kind='exact_quota_coefficient',
                            beta=(beta.numerator, beta.denominator),
                            x=(x.numerator, x.denominator), maximum=maximum,
                            tail_bits=tail_bits, entropy_bits=entropy,
                            log_Q_upper=polynomial_upper, log_x_lower=x_lower,
                            log_binomial_lower=binomial_lower, log_beta_lower=lower,
                            tail_numerator_upper=tail_numerator, certificate=True)
    return best


def selector_threshold_table(hvalues=SC_DEFAULT_H, rvalues=(1, 2, 3),
                             betas=SC_BETAS, tail_bits=SC_TAIL_BITS,
                             exact_quota=True):
    """Return certified rows keyed by(h,r); beta ties keep the first entry."""
    betas = tuple(Fraction(beta) for beta in betas)
    return {(index(h), index(r)): selector_threshold_row(h, r, betas, tail_bits, exact_quota)
            for h in hvalues for r in rvalues}


def capacity_witness(h, R, L, f, F, table=None):
    """Return a covering selector row, or None; assumes actual-state statistics."""
    h, R, L, f, F = map(index, (h, R, L, f, F))
    if not (1 <= h < SC_M and 0 <= f <= L <= F <= SC_N and R >= 128 * h):
        raise ValueError('Invalid actual-state profile for selector interface')
    G, outside = F - f, R - 128 * h
    if G % 4:
        raise ValueError('Actual G=F-f must be divisible by four')
    for r in (1, 2, 3):
        row = selector_threshold_row(h, r) if table is None else table.get((h, r))
        if row is None or not row.get('certificate', False):
            continue
        z = (4 - r) * (G // 4) - outside
        if z >= row['threshold']:
            return dict(h=h, r=r, selector_statistic=z, **{
                key: value for key, value in row.items() if key not in ('h', 'r')})
    return None


def capacity_excluded(h, R, L, f, F, table=None):
    """Boolean wrapper for a certified selected-profile class exclusion."""
    return capacity_witness(h, R, L, f, F, table) is not None


def selector_hard_certificate():
    """Exact separate hard-h class: z_2>=J/2-2 has probability<2^-4094."""
    mgf = selector_quartet_mgf(4096, 2, Fraction(3, 2))
    if not mgf <= Fraction(77, 64):
        raise ArithmeticError('Hard marginal bound failed')
    if not 8 * 5929 ** 64 < 6144 ** 64:
        raise ArithmeticError('Hard rational power comparison failed')
    # (77/64)^J*(3/2)^(-J/2+2)
    # = (5929/6144)^(J/2)*(9/4) <2^(-3J/128+2).
    if SC_J % 128 or not Fraction(9, 4) < 4:
        raise ArithmeticError('Hard exponent decomposition failed')
    exponent = SC_M - 3 * (SC_J // 128) + 2
    return dict(h=4096, r=2, threshold=SC_J // 2 - 2,
                strict_log2_probability_upper=exponent,
                global_H_union_bits=SC_M, no_quota_cost=True,
                all_death_counts=True, full_code_certificate=False)


def selector_capacity_checks():
    """Independent exact toys; this routine is run only by the parent."""
    logarithm_checks = entropy_checks = na_checks = quota_checks = 0
    for exponent in range(-40, 41):
        n, d = (1 << exponent, 1) if exponent >= 0 else (1, 1 << (-exponent))
        assert selector_log2_interval(n, d) == (exponent * SC_SCALE,) * 2
        logarithm_checks += 1
    for n in range(1, 18):
        for d in range(1, 18):
            lower, upper = selector_log2_interval(n, d)
            assert lower <= upper and upper - lower < 256
            # Independent coarse logarithmic brackets from exact16th powers.
            for k in range(-80, 81):
                lhs, rhs = n ** 16, d ** 16
                if k >= 0:
                    rhs <<= k
                else:
                    lhs <<= -k
                if lhs >= rhs:
                    assert 16 * upper >= k * SC_SCALE
                if lhs <= rhs:
                    assert 16 * lower <= k * SC_SCALE
            logarithm_checks += 1
    # Check the series floor/tail bounds independently using Fractions.
    for n, d in ((2, 1), (3, 2), (5, 4), (17, 16), (31, 17)):
        z = Fraction(n - d, n + d)
        partial = sum(2 * z ** (2 * j + 1) / (2 * j + 1) for j in range(SC_TERMS))
        tail = 2 * z ** (2 * SC_TERMS + 1) / ((2 * SC_TERMS + 1) * (1 - z * z))
        lower, upper = _sc_ln_mantissa(n, d)
        assert Fraction(lower, SC_SCALE) <= partial
        assert partial + tail <= Fraction(upper, SC_SCALE)
        logarithm_checks += 1
    for slots in range(4, 33):
        for h in range(slots + 1):
            assert comb(slots, h) <= 1 << selector_entropy_bits(h, slots)
            assert selector_binomial_log_interval(h, slots) == selector_log2_interval(comb(slots, h))
            entropy_checks += 1
    for shift in (1, 17, 128, 4096):
        integer = ((1 << 63) + 12345) << shift
        lower, upper = _sc_compact_log2_interval(integer + 1)
        prefix = selector_log2_interval((1 << 63) + 12345)
        assert lower == shift * SC_SCALE + prefix[0]
        assert lower <= upper and upper - lower < 256
        logarithm_checks += 1
    for slots in (4, 8, 12):
        groups = slots // 4
        for h in range(slots + 1):
            counts = [[0] * (slots + 1) for _ in range(4)]
            for selected in combinations(range(slots), h):
                populations = [0] * groups
                for slot in selected:
                    populations[slot // 4] += 1
                for r in range(4):
                    score = sum(max(0, value - r) for value in populations)
                    counts[r][score] += 1
            for r in range(4):
                maximum = max(i for i, count in enumerate(counts[r]) if count)
                assert maximum == selector_maximum(h, r, slots, 1)
                for beta in SC_BETAS:
                    exact = sum(count * beta ** score for score, count in enumerate(counts[r])) / comb(slots, h)
                    marginal = selector_quartet_mgf(h, r, beta, slots)
                    assert exact <= marginal ** groups
                    na_checks += 1
                    if r and 0 < h < slots:
                        x = selector_quota_tilt(h, r, beta, steps=8, slots=slots)
                        qpoly = selector_label_polynomial(r, beta, x)
                        assert exact <= qpoly ** groups / (x ** h * comb(slots, h))
                        # Independent eligible-count coefficient expansion.
                        coefficients = [Fraction(1)]
                        local = [comb(4, j) * beta ** max(0, j - r) for j in range(5)]
                        for _ in range(groups):
                            product = [Fraction(0)] * (len(coefficients) + 4)
                            for a, first in enumerate(coefficients):
                                for b, second in enumerate(local):
                                    product[a + b] += first * second
                            coefficients = product
                        assert exact == coefficients[h] / comb(slots, h)
                        quota_checks += 1
    hard = selector_hard_certificate()
    assert hard['strict_log2_probability_upper'] == -4094
    table = selector_threshold_table((1, 256, 4096, 8191))
    for row in table.values():
        if row['kind'] == 'deterministically_impossible':
            assert row['threshold'] == row['maximum'] + 1
        else:
            assert row['tail_numerator_upper'] - row['threshold'] * row['log_beta_lower'] <= -SC_TAIL_BITS * SC_SCALE
        na_row = selector_threshold_row(row['h'], row['r'], exact_quota=False)
        assert row['threshold'] <= na_row['threshold']
    assert capacity_excluded(4096, 786432, 786432, 4, 1048576, table)
    assert not capacity_excluded(4096, 786432, 786432, 4, 1048576, {})
    assert (SC_M - 1) * 3 < (1 << 15)
    print(dict(logarithm_checks=logarithm_checks, entropy_checks=entropy_checks,
               exact_NA_toys=na_checks, exact_quota_toys=quota_checks,
               threshold_rows=len(table),
               hard_certificate=hard, full_code_certificate=False), flush=True)


def selector_capacity_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--hvalues', type=int, nargs='*', help='Empty list uses256-step coarse grid')
    mode.add_argument('--all-h', action='store_true')
    parser.add_argument('--na-only', action='store_true', help='Keep the earlier NA control')
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            selector_capacity_checks()
        else:
            hs = range(1, SC_M) if args.all_h else (args.hvalues or SC_DEFAULT_H)
            print(dict(hard_certificate=selector_hard_certificate()), flush=True)
            table = selector_threshold_table(hs, exact_quota=not args.na_only)
            for row in table.values():
                print(row, flush=True)
            print(dict(rows=len(table), pair_union_log2_upper=-77,
                       seconds=round(time.monotonic() - started, 3),
                       exact_integer_certificate=True, full_code_certificate=False), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    selector_capacity_main()
