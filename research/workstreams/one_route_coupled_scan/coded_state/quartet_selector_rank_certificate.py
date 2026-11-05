"""Exact actual-L profile classes from live selectors and grouped amplitudes.

Inject quartet_selector_capacity.py directed logarithm helpers.  No matrix
witness or floating-point arithmetic enters this certificate.

The actual-L quotient bound is applied before summing selectors.  The
theta=1 specialization beta4,volume1/4,v_i=16^i gives
M<=(17/16)^(128*(4096-1))*2^512 from five local rational comparisons.
The theta=0 geometry specialization beta2,volume1/2,v_i=4^i gives
M<=(5/4)^(128*(4096-1))*2^256.  Both retain the formal free terminal
without adding a selector layer or a terminal volume marker.

Births cost C(N-G,f).  The capped forced-birth route factor already
includes the H-set union; there is no second H union.  These are fixed
actual-profile bounds, so the strict2^96 profile union remains necessary.
An integer exponent<=-172 certifies its profile class below2^-76.
The full distance goal remains open.  Parent-only, serial and file-free.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time


SRC_N = 1 << 21
SRC_M = 8192
SRC_Q = 1 << 32
SRC_WIDTH = 512
SRC_UPDATES = 4096
SRC_POWER = 128 * (SRC_UPDATES - 1)
SRC_SCALE = 1 << 48
SRC_THRESHOLD = -172
SRC_DEFAULT_PROFILE = (1536, 753664, 730726, 4, 851968)


def _src_ceildiv(a, b):
    if b <= 0:
        raise ValueError('Positive divisor required')
    return -((-a) // b)


def _src_dependencies():
    if globals().get('SC_SCALE') != SRC_SCALE:
        raise RuntimeError('Inject the2^48 selector-capacity directed log helpers')
    for name in ('selector_log2_interval', 'selector_binomial_log_interval'):
        if name not in globals():
            raise RuntimeError('Inject quartet_selector_capacity.py definitions')


def selector_rank_logbinomial_upper(n, k):
    """Directed entropy upper in scaled bits; never constructs a large binomial."""
    n, k = index(n), index(k)
    if n < 0 or not 0 <= k <= n:
        raise ValueError('Invalid binomial dimensions')
    if k in (0, n):
        return 0
    lp = selector_log2_interval(k, n)[0]
    lq = selector_log2_interval(n - k, n)[0]
    return -k * lp - (n - k) * lq


def selector_rank_constants():
    """Exact local supersolutions, projective factor and free-endpoint bounds."""
    for c in range(5):
        assert 1 + Fraction(4) ** (c - 2) <= Fraction(17, 16) * 2 ** c
        assert (1 + Fraction(2) ** (c - 2)) ** 2 <= Fraction(25, 16) * 2 ** c
    assert Fraction(1024, 989) / (SRC_Q - 1) < Fraction(1, 1 << 30)
    assert Fraction(4) ** (4 * 128) / Fraction(16) ** 128 == 1 << 512
    assert Fraction(2) ** (4 * 128) / Fraction(4) ** 128 == 1 << 256
    assert SRC_POWER == 524160
    assert 8192 * (1 << 20) * SRC_N ** 3 == 1 << 96
    return dict(exact_local_comparisons=10, exact_terminal_comparisons=2,
                profile_union_bits=96, qualifying_threshold=SRC_THRESHOLD,
                full_code_certificate=False)


def selector_rank_profile(h, R, Ls, factual, F):
    """Two exact profile exponents, with the routing factor capped before use."""
    _src_dependencies()
    h, R, Ls, factual, F = map(index, (h, R, Ls, factual, F))
    if not 1 <= h <= SRC_M or not 1 <= factual <= Ls <= F <= SRC_N or R < 128 * h:
        raise ValueError('Invalid actual-state profile')
    G = F - factual
    C_H = F - R + 128 * h
    holes_H = 256 * h - C_H
    if G % 4:
        raise ValueError('Actual G=F-f must be divisible by four')
    if not (G <= SRC_N - SRC_WIDTH and 129 * h <= C_H <= 256 * h
            and 0 <= C_H - factual <= G and 0 <= holes_H <= SRC_N - F):
        return dict(h=h, R=R, Ls=Ls, factual=factual, F=F,
                    status='empty', failed_necessary_route_or_occupancy_condition=True,
                    full_code_certificate=False)
    birth_upper = selector_rank_logbinomial_upper(SRC_N - G, factual)
    route_uncapped = (selector_rank_logbinomial_upper(G, C_H - factual)
                      + selector_rank_logbinomial_upper(SRC_N - F, holes_H)
                      - 255 * selector_binomial_log_interval(h)[0])
    route_upper = min(0, route_uncapped)
    log3_lower = selector_log2_interval(3)[0]
    log17_upper = selector_log2_interval(17, 16)[1]
    log5_upper = selector_log2_interval(5, 4)[1]
    actual_free = Ls - factual
    rank_upper = (birth_upper + route_upper
                  + (-30 + 32 * (Ls - R) + 2 * G + 512) * SRC_SCALE
                  - actual_free * log3_lower + SRC_POWER * log17_upper)
    geometry_upper = (birth_upper + route_upper
                      + (G - actual_free + 256) * SRC_SCALE
                      + SRC_POWER * log5_upper)
    values = dict(rank_theta1=_src_ceildiv(rank_upper, SRC_SCALE),
                  geometry_theta0=_src_ceildiv(geometry_upper, SRC_SCALE))
    selected = min(values, key=values.get)
    exponent = values[selected]
    return dict(h=h, R=R, Ls=Ls, factual=factual, F=F, G=G,
                C_H=C_H, C_outside=R - 128 * h, status='nonempty_necessary_conditions',
                log_births_upper=birth_upper, log_route_uncapped_upper=route_uncapped,
                log_route_capped_upper=route_upper, directed_log_scale=SRC_SCALE,
                profile_exponent_uppers=values, selected_branch=selected,
                exact_profile_exponent_upper=exponent,
                qualifies_global_profile_class=exponent <= SRC_THRESHOLD,
                qualifying_threshold=SRC_THRESHOLD, profile_union_bits=96,
                class_failure_exponent=-76, H_union_already_in_route=True,
                actual_L_retained_before_amplitude_enlargement=True,
                free_endpoint=True, production_matrix_constructed=False,
                exact_integer_certificate=True, full_code_certificate=False)


def _src_exact_overlap(width, u, v, beta):
    return sum(Fraction(comb(u, ell) * comb(width - u, v - ell), comb(width, v))
               * beta ** ell for ell in range(max(0, v - width + u), min(u, v) + 1))


def selector_rank_checks():
    _src_dependencies()
    constants = selector_rank_constants()
    row_checks = endpoint_checks = entropy_checks = 0
    for width in (4, 8):
        m = width // 4
        for beta, volume, potential, lam, endpoint in (
                (Fraction(4), Fraction(1, 4), Fraction(16), Fraction(17, 16), 1 << width),
                (Fraction(2), Fraction(1, 2), Fraction(4), Fraction(5, 4), 1 << (width // 2))):
            matrix = [[comb(m, b) * volume ** (4 * b)
                       * _src_exact_overlap(width, 4 * a, 4 * b, beta)
                       for b in range(m + 1)] for a in range(m + 1)]
            vector = [potential ** a for a in range(m + 1)]
            terminal = [beta ** (4 * a) for a in range(m + 1)]
            for a in range(m + 1):
                assert sum(matrix[a][b] * vector[b] for b in range(m + 1)) <= lam ** m * vector[a]
                assert terminal[a] <= endpoint * vector[a]
                row_checks += 1
            row = [Fraction(int(a == 0)) for a in range(m + 1)]
            for updates in range(1, 5):
                exact = sum(row[a] * terminal[a] for a in range(m + 1))
                assert exact <= endpoint * lam ** (m * (updates - 1))
                row = [sum(row[a] * matrix[a][b] for a in range(m + 1)) for b in range(m + 1)]
                endpoint_checks += 1
    for n in range(2, 65):
        for k in range(n + 1):
            upper = selector_rank_logbinomial_upper(n, k)
            assert selector_log2_interval(comb(n, k))[1] <= upper + 1024
            entropy_checks += 1
    print(dict(exact_small_rows=row_checks, exact_free_endpoint_checks=endpoint_checks,
               directed_entropy_checks=entropy_checks, **constants), flush=True)


def selector_rank_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, action='append',
                        metavar=('H', 'R', 'LS', 'F_BIRTH', 'F'))
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
            selector_rank_checks()
        else:
            selector_rank_constants()
            for profile in args.profile or (SRC_DEFAULT_PROFILE,):
                print(selector_rank_profile(*profile), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    selector_rank_main()
