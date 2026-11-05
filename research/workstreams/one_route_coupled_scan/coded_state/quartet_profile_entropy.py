"""Availability-run census times exact-quota routing capacity, by profile.

Fix a candidate availability path before exposing partitions.  Its
canonical footprint is determined by the partitions.  For fixed footprint
size F and a fixed eligible-row set H of size h, the H-labeled positions
are independent uniform h-subsets of each 8192-position region.  Writing
C_H=F-R+128h, dropping the regional restrictions in the favorable numerator
gives probability at most

    C(F-f,C_H-f) C(N-F,256h-C_H) / C(8192,h)^256,

because the f route-independent canonical forced births must all lie in
H.  The original unforced numerator C(F,C_H) remains an explicit control.

Union over H and the availability census B(L,D), WITHOUT a partition
union, gives the geometry expression used below.  It is an expected-count
upper bound; it must not be capped at one before multiplying a conditional
per-footprint field bound.  The routing probability, including the H union,
is capped at one BEFORE multiplication by the availability census.
FLOAT logarithms and a bounded field-marker
grid are diagnostics, never directed certificates or a full profile cover.
All setup, construction, rank grouping, quotas and free endpoints are fixed.
"""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb, lgamma, log2
import os
import shutil
import threading
import time

import numpy as np


PE_LD = np.longdouble
PE_N = 1 << 21
PE_WIDTH = 512
PE_T = 4096
PE_ROWS = 8192
PE_REGIONS = 256
PE_W = 943718
PE_Q = 1 << 32


def profile_logchoose(n, k):
    """Direct natural-log gamma difference, converted once to base two."""
    if n < 0 or k < 0 or k > n:
        return PE_LD('-inf')
    if k == 0 or k == n:
        return PE_LD(0)
    k = min(k, n - k)
    return PE_LD(lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / np.log(PE_LD(2))
def profile_log_census(L, D, width=PE_WIDTH, updates=PE_T):
    N = width * updates
    if not 0 <= L <= N or D < 0 or D > min(L, N - L):
        return PE_LD('-inf')
    if L == 0:
        return PE_LD(0) if D == 0 else PE_LD('-inf')
    terms = []
    for terminal in range(width + 1):
        runs = D + terminal
        if not 1 <= runs <= L or runs > N - L + width:
            continue
        terms.append(profile_logchoose(L - 1, runs - 1)
                     + profile_logchoose(N - L + width, runs))
    if not terms:
        return PE_LD('-inf')
    largest = max(terms)
    run_bound = largest + np.log2(np.exp2(np.array(terms, dtype=PE_LD) - largest).sum())
    return min(profile_logchoose(N, L), run_bound)


def profile_census_constant(j):
    if not isinstance(j, int) or not 0 <= j <= 16:
        raise ValueError('The recorded weighted-census family uses integer j in 0..16')
    if j == 0:
        return PE_N
    if j == 3:
        return 361984
    denominator = 1 << (j + 1)
    return 512 * j + (3 * PE_N + denominator - 1) // denominator


def profile_feasibility(h, R, L, f, F, D):
    """Necessary integer screens and already-proved class exclusions.

    Passing these tests asserts only aggregate feasibility, never an actual
    availability path, routing realization, bad geometry or codeword.
    """
    values = (h, R, L, f, F, D)
    if any(not isinstance(value, (int, np.integer)) for value in values):
        raise ValueError('Profiles must have exact integer coordinates')
    empty = []
    if not 1 <= h <= PE_ROWS:
        empty.append('eligible-row count outside 1..8192')
    if not 1 <= f <= L <= F <= PE_N:
        empty.append('require 1<=f<=L<=F<=N')
    if f > PE_W:
        empty.append('forced nonzero outputs already exceed W')
    C_H = F - R + 128 * h
    C_out = R - 128 * h
    if not 129 * h <= C_H <= 256 * h:
        empty.append('eligible-row degree capacity: h<=F-R<=128h')
    if f > C_H:
        empty.append('forced births must be covered positions in eligible rows')
    if not 0 <= C_out <= 128 * (PE_ROWS - h):
        empty.append('outside-row degree capacity')
    if 2 * R < F or R > F:
        empty.append('require F/2<=R<=F')
    if (F - f) % 4:
        empty.append('full-quartet count requires F-f divisible by four')
    if not max(0, f - PE_WIDTH) <= D <= min(F - L, L, PE_N - L):
        empty.append('availability starts/closures: max(0,f-512)<=D<=min(F-L,L,N-L)')
    if not L - PE_WIDTH <= F - f <= 4 * L:
        empty.append('quartet expansion capacity')
    if 2 * F - L - f > PE_N - PE_WIDTH:
        empty.append('canonical old-hole / hit-quartet exclusion')
    if F - f + 4 * max(0, C_out - L + f) > 4 * L:
        empty.append('accepted-MDS aggregate after eliminating available-H count')
    if (F - f) % 4 == 0:
        terminal_lower = max(0, f - D)
        terminal_upper = min(PE_WIDTH, L - D,
                             L - (F - f) // 4 - max(0, C_out - L + f))
        if terminal_lower > terminal_upper:
            empty.append('incompatible free-terminal range from runs and MDS capacity')
    if empty:
        return dict(status='empty', reasons=empty, C_H=C_H, C_out=C_out)
    certified = []
    if F >= 1736704:
        certified.append('previous dense-footprint class')
    if R >= 1009256:
        certified.append('previous absolute-rank class')
    if R - L >= 65537:
        certified.append('previous rank-minus-availability class')
    if 2 * R + max(0, f - PE_WIDTH) >= 2 * 1009257:
        certified.append('previous rank-plus-half-state-deaths class via forced runs')
    for j in range(17):
        constant = profile_census_constant(j)
        if 32 * (R - L) - 2 * j * D >= constant + 45:
            certified.append(f'weighted rank census j={j}')
        if 32 * R - 16 * L + (15 - 2 * j) * D >= 16 * PE_W + constant + 92:
            certified.append(f'weighted closure census j={j}')
    return dict(status='certified' if certified else 'remaining', reasons=certified,
                C_H=C_H, C_out=C_out, eligible_holes=256 * h - C_H)


def profile_field_grid(step=.25, maximum=32., q=PE_Q):
    if step <= 0 or maximum < 0:
        raise ValueError('Require positive field-grid step and nonnegative range')
    count = int(np.floor(maximum / step))
    logx = -PE_LD(step) * np.arange(count + 1, dtype=PE_LD)
    if count * step < maximum:
        logx = np.append(logx, -PE_LD(maximum))
    x = np.exp2(logx)
    c = (PE_LD(q) + 1) * x
    lam = 1 + c
    closed = np.maximum(np.minimum(PE_LD(1), lam / (PE_LD(q) + 1)),
                        (1 + PE_LD(q) * x * x) / lam)
    return dict(q=q, logx=logx, logc=np.log2(c), loglambda=np.log2(lam), logd=np.log2(closed))


def profile_field_bound(R, L, f, D, grid=None, W=PE_W):
    if grid is None:
        grid = profile_field_grid()
    logq = np.log2(PE_LD(grid['q']))
    values = (1 - logq - logq * R - W * grid['logx'] + f * grid['logc']
              + (L - f) * grid['loglambda'] + D * grid['logd'])
    index = int(np.argmin(values))
    return min(PE_LD(0), values[index]), float(grid['logx'][index])


def profile_entropy_value(h, R, L, f, F, D, grid=None, forced=True):
    screening = profile_feasibility(h, R, L, f, F, D)
    result = dict(h=h, R=R, L=L, f=f, F=F, D=D, forced_birth_routing=bool(forced), **screening)
    if screening['status'] != 'remaining':
        return result
    C_H = screening['C_H']
    census = profile_log_census(L, D)
    routing_raw = (profile_logchoose(F - f, C_H - f) if forced else profile_logchoose(F, C_H))
    routing_raw += (profile_logchoose(PE_N - F, PE_REGIONS * h - C_H)
                    - 255 * profile_logchoose(PE_ROWS, h))
    routing = min(PE_LD(0), routing_raw)
    field, logx = profile_field_bound(R, L, f, D, grid)
    result.update(log2_availability_census=float(census), log2_route_union_before_probability_cap=float(routing_raw),
                  log2_route_factor_with_H_union=float(routing),
                  log2_expected_geometry_count=float(census + routing), log2_conditional_field=float(field),
                  selected_log2_x=logx, log2_diagnostic=float(census + routing + field), certificate=False)
    return result


def _profile_exact_census_bound(L, D, width, updates):
    N = width * updates
    if L == 0:
        return int(D == 0)
    total = 0
    for terminal in range(width + 1):
        runs = D + terminal
        if 1 <= runs <= L and runs <= N - L + width:
            total += comb(L - 1, runs - 1) * comb(N - L + width, runs)
    return min(comb(N, L), total)


def profile_entropy_checks():
    census_checks = quota_checks = forced_quota_checks = field_checks = screening_checks = 0
    for width, updates in ((1, 4), (2, 3), (3, 3), (4, 3)):
        counts = Counter()
        full = (1 << width) - 1
        for history in range(1 << (width * updates)):
            masks = [(history >> (t * width)) & full for t in range(updates)]
            L = sum(mask.bit_count() for mask in masks)
            D = sum((masks[t] & ~masks[t + 1]).bit_count() for t in range(updates - 1))
            counts[L, D] += 1
        for (L, D), count in counts.items():
            bound = _profile_exact_census_bound(L, D, width, updates)
            assert count <= bound
            assert abs(profile_log_census(L, D, width, updates) - PE_LD(log2(bound))) <= PE_LD('2e-12')
            census_checks += 1
    for regions, selected in ((2, 1), (2, 2), (3, 1), (3, 2)):
        width = 4
        choices = [sum(1 << j for j in subset) for subset in combinations(range(width), selected)]
        labels = [sum(mask << (width * region) for region, mask in enumerate(picked))
                  for picked in product(choices, repeat=regions)]
        denominator = comb(width, selected) ** regions
        assert len(labels) == denominator
        for regional_coverage in product(range(width + 1), repeat=regions):
            covered = sum(((1 << v) - 1) << (width * region) for region, v in enumerate(regional_coverage))
            F = sum(regional_coverage)
            distribution = Counter((covered & mask).bit_count() for mask in labels)
            for C_H, count in distribution.items():
                numerator = comb(F, C_H) * comb(width * regions - F, regions * selected - C_H)
                assert Fraction(count, denominator) <= Fraction(numerator, denominator)
                quota_checks += 1
            positions = [j for j in range(width * regions) if covered >> j & 1]
            for f in sorted({0, min(1, F), F // 2, F}):
                births = sum(1 << j for j in positions[:f])
                forced_distribution = Counter((covered & mask).bit_count() for mask in labels if births & ~mask == 0)
                for C_H, count in forced_distribution.items():
                    numerator = comb(F - f, C_H - f) * comb(width * regions - F, regions * selected - C_H)
                    unforced = comb(F, C_H) * comb(width * regions - F, regions * selected - C_H)
                    assert count <= numerator <= unforced
                    forced_quota_checks += 1
    for q in (2, 4, 16):
        for exponent in range(0, 2 * int(log2(q)) + 1):
            x = Fraction(1, 1 << exponent)
            c = (q + 1) * x
            lam = 1 + c
            d = max(min(Fraction(1), lam / (q + 1)), (1 + q * x * x) / lam)
            for R, L, f, D, W in ((2, 3, 1, 1, 3), (4, 5, 2, 2, 4), (1, 2, 1, 0, 2)):
                exact = Fraction(2, q) * Fraction(1, q) ** R * x ** (-W) * c ** f * lam ** (L - f) * d ** D
                grid = profile_field_grid(step=1., maximum=exponent, q=q)
                index = len(grid['logx']) - 1
                value = (1 - np.log2(PE_LD(q)) - R * np.log2(PE_LD(q)) - W * grid['logx'][index]
                         + f * grid['logc'][index] + (L - f) * grid['loglambda'][index] + D * grid['logd'][index])
                expected = log2(exact.numerator) - log2(exact.denominator)
                assert abs(value - expected) <= PE_LD('2e-12')
                field_checks += 1
    old = profile_feasibility(2048, 900000, 900000, 20000, 1200000, 30000)
    assert old['status'] == 'empty' and any('eligible-row' in reason for reason in old['reasons'])
    screening_checks += 1
    new = profile_feasibility(2048, 900000, 900000, 20000, 1150000, 30000)
    assert new['status'] == 'remaining' and new['eligible_holes'] == 12144
    screening_checks += 1
    impossible_run = profile_feasibility(2048, 900000, 900000, 20000, 1150000, 10000)
    assert impossible_run['status'] == 'empty'
    screening_checks += 1
    print(dict(exact_availability_census_checks=census_checks, exact_quota_numerator_checks=quota_checks,
               exact_forced_birth_quota_numerator_checks=forced_quota_checks,
               rational_field_grid_comparisons=field_checks, integer_profile_screen_checks=screening_checks,
               full_code_certificate=False), flush=True)


def profile_entropy_scan(step=1024, field_step=.25, top=5, only=None, forced=True):
    if step < 1 or top < 1:
        raise ValueError('Require positive integer D step and report count')
    grid = profile_field_grid(field_step)
    profiles = [("new_capacity_feasible", 2048, 900000, 900000, 20000, 1150000)]
    # The recorded sensitivities are rounded real diagnostics, not integer paths.
    # Increase f by one to satisfy the necessary F-f=0 mod4 relation.
    profiles.append(("observed_sensitivity_f_adjusted_plus_one", 2048, 903775, 891035, 71427, 1132763))
    if only is not None:
        profiles = [profiles[only]]
    for label, h, R, L, f, F in profiles:
        lower, upper = max(0, f - PE_WIDTH), min(F - L, L, PE_N - L)
        candidates = set(range(lower, upper + 1, step))
        candidates.update((lower, upper))
        candidates.update(range(45000, 60001, 1000))
        critical = PE_W - L
        candidates.update(critical + offset for offset in (-1024, -512, 0, 512, 1024))
        candidates = sorted(D for D in candidates if lower <= D <= upper)
        evaluated = []
        exclusions = Counter()
        for D in candidates:
            result = profile_entropy_value(h, R, L, f, F, D, grid, forced)
            if result['status'] == 'remaining':
                evaluated.append(result)
            else:
                exclusions[result['status']] += 1
        evaluated.sort(key=lambda result: result['log2_diagnostic'], reverse=True)
        print(dict(profile=label, h=h, R=R, L=L, f=f, F=F, D_step=step,
                   evaluated=len(evaluated), exclusions=dict(exclusions), field_grid_step=field_step,
                   forced_birth_routing=bool(forced),
                   finite_scan_only=True, certificate=False), flush=True)
        for result in evaluated[:top]:
            print(dict(profile=label, **result), flush=True)


def profile_entropy_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--D-step', type=int, default=1024)
    parser.add_argument('--field-step', type=float, default=.25)
    parser.add_argument('--top', type=int, default=5)
    parser.add_argument('--only', type=int, choices=(0, 1))
    parser.add_argument('--unforced', action='store_true', help='Use the original routing numerator as a control')
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
            profile_entropy_checks()
        else:
            profile_entropy_scan(args.D_step, args.field_step, args.top, args.only, not args.unforced)
            print(dict(seconds=round(time.monotonic() - started, 3), certificate=False), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    profile_entropy_main()
