"""Closed geometry-only availability/occupancy count (closure marker t=1).

For A_0=empty and a free A_T, let L=sum|A_t| and
G=sum|Gamma_t(A_(t-1))|.  For width n and T>=1 updates,

  sum_A s^L z^G = (1+s)^n [1+z^4((1+s)^4-1)]^(n(T-1)/4).

The identity holds for EVERY fixed sequence of quartet partitions.  The
terminal mask has no subsequent quartet charge.  For canonical histories
G=F-f; arbitrary masks counted here need not be canonical or realizable.
No field, routing, or distance certificate is included.

occupancy_geometry_logbound(L,G) returns a FLOAT upper-bound diagnostic
for the expected number of histories.  It optimizes the z marker in
closed form and s by scalar bisection.  It must never be capped at one.
FLOAT evaluation is not directed arithmetic or an exact certificate.
The --checks mode verifies the closed polynomial against independent
physical-mask enumerations with exact integers and rational moments.
Standalone; no local imports, subprocesses or data files.  Run jobs serially.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, exp, isfinite, log, log1p
from operator import index
import os
import shutil
import threading
import time


OG_LN2 = log(2.0)
OG_INV_LN2 = 1.0 / OG_LN2


def _og_dimensions(width, updates):
    width, updates = index(width), index(updates)
    if width < 4 or width % 4:
        raise ValueError('Width must be a positive integer multiple of four')
    if updates < 0:
        raise ValueError('Updates must be nonnegative')
    return width, updates


def _og_entropy(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * log(p) + (1.0 - p) * log1p(-p)) * OG_INV_LN2


def _og_logadd2(left, right):
    if left < right:
        left, right = right, left
    return left + log1p(exp((right - left) * OG_LN2)) * OG_INV_LN2


def _og_log_quartet(log_s):
    # (1+s)^4-1 = s*(4+6s+4s^2+s^3).  This avoids subtraction
    # when s is small, and no large power is materialized when s is large.
    terms = (2.0, log(6.0) * OG_INV_LN2 + log_s,
             2.0 + 2.0 * log_s, 3.0 * log_s)
    maximum = max(terms)
    return log_s + maximum + log(sum(exp((v - maximum) * OG_LN2)
                                    for v in terms)) * OG_INV_LN2


def occupancy_geometry_logmoment(width=512, updates=4096, log_s=0.0, log_z=0.0):
    """FLOAT log2 of the t=1 moment, with arbitrary finite log2 markers."""
    width, updates = _og_dimensions(width, updates)
    log_s, log_z = float(log_s), float(log_z)
    if not isfinite(log_s) or not isfinite(log_z):
        raise ValueError('Logarithmic markers must be finite')
    if updates == 0:
        return 0.0
    free_endpoint = width * _og_logadd2(0.0, log_s)
    groups = (width // 4) * (updates - 1)
    if groups == 0:
        return free_endpoint
    return free_endpoint + groups * _og_logadd2(0.0, 4.0 * log_z + _og_log_quartet(log_s))


def _og_optimum(L, G, width, updates):
    """Return (log2 bound, p=s/(1+s), k=G/(n(T-1)), empty).

    The hot scalar loop uses only arithmetic.  Its monotone equation is
      L = n*p + G/(4-6p+4p^2-p^3).
    The original form is
      L = n*s/(1+s) + G*s*(1+s)^3/((1+s)^4-1).
    """
    if updates == 0:
        return (0.0, 0.0, 0.0, False) if L == 0 and G == 0 else (float('-inf'), None, None, True)
    maximum_G = width * (updates - 1)
    if G < 0 or G > maximum_G or G % 4 or L < G // 4 or L > G + width:
        return float('-inf'), None, None, True
    if G == 0:
        # All masks before the free endpoint are empty.  The exact count is
        # C(width,L); the optimized Chernoff bound is 2^(width*H2(L/width)).
        p = L / width
        return width * _og_entropy(p), p, 0.0, False

    hits = G // 4
    groups = maximum_G // 4
    k = G / maximum_G
    occupancy_entropy = groups * _og_entropy(k)
    if L == hits:
        # s -> 0: the exact coefficient is C(groups,hits)*4^hits.
        return occupancy_entropy + 2.0 * hits, 0.0, k, False
    if L == G + width:
        # s -> infinity: the exact coefficient is C(groups,hits).
        return occupancy_entropy, 1.0, k, False

    lower, upper = 0.0, 1.0
    # 64 rounds exceed binary64's significand precision.  Adjacent
    # representable endpoints stop early; no allocation occurs in the loop.
    for _ in range(64):
        p = (lower + upper) * 0.5
        if p == lower or p == upper:
            break
        denominator = 4.0 - p * (6.0 - p * (4.0 - p))
        if width * p + G / denominator < L:
            lower = p
        else:
            upper = p
    p = (lower + upper) * 0.5
    denominator = 4.0 - p * (6.0 - p * (4.0 - p))
    # All three coefficients below are nonnegative on the supported range.
    # This form avoids cancellation of large log(1+s) and log(s) terms.
    value = occupancy_entropy + (-(L - hits) * log(p)
             - (G + width - L) * log1p(-p)
             + hits * log(denominator)) * OG_INV_LN2
    return value, p, k, False


def occupancy_geometry_logbound(L, G, width=512, updates=4096):
    """FLOAT log2 expected-count upper bound for one integer (L,G) profile.

    The exact supported range is G in 4Z, 0<=G<=width*(updates-1),
    G/4<=L<=G+width.  Empty profiles return -infinity.  With no updates,
    only (0,0) exists, with count one.  The all-empty history is included;
    callers imposing nonzero codewords must exclude it separately.

    Writing J=width*(updates-1)/4 and k=G/(4J), the bound minimizes
      J*H2(k)+(G/4)*log2((1+s)^4-1)
        +width*log2(1+s)-L*log2(s).
    The J*H2(k) term is the optimized z-marker bound, not log2 C(J,G/4).
    """
    width, updates = _og_dimensions(width, updates)
    L, G = index(L), index(G)
    return _og_optimum(L, G, width, updates)[0]


def occupancy_geometry_profile(L, G, width=512, updates=4096):
    """The same FLOAT bound with optimizer and endpoint metadata."""
    width, updates = _og_dimensions(width, updates)
    L, G = index(L), index(G)
    value, p, k, empty = _og_optimum(L, G, width, updates)
    log_s = log_z = None
    if not empty:
        if p == 0.0:
            log_s = float('-inf')
        elif p == 1.0:
            log_s = float('inf')
        else:
            log_s = (log(p) - log1p(-p)) * OG_INV_LN2
            if k == 0.0:
                log_z = float('-inf')
            elif k == 1.0:
                log_z = float('inf')
            else:
                log_z = ((log(k) - log1p(-k)) * OG_INV_LN2 - _og_log_quartet(log_s)) / 4.0
    return dict(width=width, updates=updates, L=L, G=G, empty=empty,
                log2_expected_count_upper=value, optimized_p=p, hit_fraction=k,
                optimized_log_s=log_s, optimized_log_z=log_z,
                boundary_limit=(not empty and (p in (0.0, 1.0) or k in (0.0, 1.0))),
                all_empty_profile=(not empty and L == 0 and G == 0),
                free_endpoint=True, closure_marker=1,
                floating_diagnostic_only=True, full_code_certificate=False)


def _og_exact_table(width, updates):
    """Exact polynomial coefficients for toy checks only; keys are (L,G)."""
    if updates == 0:
        return {(0, 0): 1}
    groups = (width // 4) * (updates - 1)
    table = {(0, 0): 1}
    for _ in range(groups):
        after = dict(table)
        for (L, G), value in table.items():
            for count in range(1, 5):
                key = L + count, G + 4
                after[key] = after.get(key, 0) + value * comb(4, count)
        table = after
    result = {}
    for (L, G), value in table.items():
        for terminal in range(width + 1):
            key = L + terminal, G
            result[key] = result.get(key, 0) + value * comb(width, terminal)
    return result


def _og_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        return [(first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                for first in [1 + sum(1 << j for j in rest)]]
    raise ValueError('Independent physical checks use widths four and eight')


def _og_physical_table(width, updates, partition):
    # Enumerate actual masks.  Only masks before the terminal contribute
    # hit quartets.  With t=1 there is no transition-intersection weight.
    populations = [mask.bit_count() for mask in range(1 << width)]
    covered = [4 * sum(bool(mask & quartet) for quartet in partition)
               for mask in range(1 << width)]
    result = {}
    for masks in product(range(1 << width), repeat=updates):
        L = sum(populations[mask] for mask in masks)
        G = sum(covered[mask] for mask in masks[:-1])
        key = L, G
        result[key] = result.get(key, 0) + 1
    return result


def occupancy_geometry_checks():
    coefficient_checks = moment_checks = endpoint_checks = float_checks = 0
    profiles = []
    for width, updates in ((4, 0), (4, 1), (4, 2), (4, 3), (8, 0), (8, 1), (8, 2)):
        expected = _og_exact_table(width, updates)
        for partition in _og_partitions(width):
            actual = _og_physical_table(width, updates, partition)
            assert actual == expected
            coefficient_checks += len(actual)
        assert sum(expected.values()) == 1 << (width * updates)
        for s, z in ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(2)),
                     (Fraction(1), Fraction(1)), (Fraction(3, 4), Fraction(2, 3)),
                     (Fraction(5, 4), Fraction(7, 5))):
            enumerated = sum(value * s ** L * z ** G for (L, G), value in expected.items())
            formula = (Fraction(1) if updates == 0 else
                       (1 + s) ** width * (1 + z ** 4 * ((1 + s) ** 4 - 1)) ** ((width // 4) * (updates - 1)))
            assert enumerated == formula
            moment_checks += 1
        if updates:
            groups = (width // 4) * (updates - 1)
            for hits in range(groups + 1):
                G = 4 * hits
                assert expected[hits, G] == comb(groups, hits) * 4 ** hits
                assert expected[G + width, G] == comb(groups, hits)
                endpoint_checks += 2
            for L in range(width + 1):
                assert expected[L, 0] == comb(width, L)
                endpoint_checks += 1
        for L in range(width * updates + 2):
            for G in range(width * max(0, updates - 1) + 5):
                result = occupancy_geometry_profile(L, G, width, updates)
                assert result['empty'] == ((L, G) not in expected)
                if not result['empty']:
                    # Floating optimizer diagnostics are checked separately
                    # from the exact identities above, with a stated tolerance.
                    exact_log = log(expected[L, G]) * OG_INV_LN2
                    assert result['log2_expected_count_upper'] >= exact_log - 1e-10
                    p, k = result['optimized_p'], result['hit_fraction']
                    if 0 < p < 1 and 0 < k < 1:
                        log_s, log_z = result['optimized_log_s'], result['optimized_log_z']
                        evaluated = occupancy_geometry_logmoment(width, updates, log_s, log_z) - L * log_s - G * log_z
                        assert abs(evaluated - result['log2_expected_count_upper']) < 1e-9
                        denominator = 4.0 - p * (6.0 - p * (4.0 - p))
                        assert abs(width * p + G / denominator - L) < 1e-10
                    float_checks += 1
        profiles.append((width, updates))
    assert occupancy_geometry_logbound(-1, 0, 4, 2) == float('-inf')
    assert occupancy_geometry_logbound(0, -4, 4, 2) == float('-inf')
    assert occupancy_geometry_logbound(0, 0, 4, 0) == 0.0
    print(dict(exact_physical_coefficient_checks=coefficient_checks,
               exact_rational_moment_checks=moment_checks,
               exact_free_endpoint_checks=endpoint_checks,
               floating_optimizer_checks=float_checks, checked_toys=profiles,
               full_code_certificate=False), flush=True)


def occupancy_geometry_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--width', type=int, default=512)
    parser.add_argument('--updates', type=int, default=4096)
    parser.add_argument('--L', type=int, default=900000)
    parser.add_argument('--G', type=int, default=1130000)
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
            occupancy_geometry_checks()
        else:
            print(dict(seconds=round(time.monotonic() - started, 3),
                       **occupancy_geometry_profile(args.L, args.G, args.width, args.updates)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    occupancy_geometry_main()
