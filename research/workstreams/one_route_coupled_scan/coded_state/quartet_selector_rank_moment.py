"""FLOAT fixed-actual-L quotient bound over consistent live selectors.

For G=F-f and m=L-f, let A0 be the sum of adjacent live-set overlaps,
including the formal free endpoint U_(T+1)=[width]. Grouping amplitudes
at fixed(U,B) gives a field bound

 min(1, 2^(32*(L-R)-30)*C(A0,m)).

Apply theta only after this grouping. For0<theta<=1 choose
y=beta^(1/theta)-1, yielding beta^A0*y^(-theta*m). At theta0 retain only
A0>=m, giving beta^(A0-m). These are ACTUAL L and f, not allowed volume.

H_ij(beta)=E beta^|U intersect V| is supplied by temporal_selector_overlap.
K_ij=H_ij*C(width/4,j)*z^(4j). The exact selector moment is
e0*K^(T-1)*t with t_i=beta^(4i). Initially we upper-bound t by
beta^width*1 and reuse the injected positive _st_power_row helper.

The remaining factors are C(N-G,f), z^-G, and min(1,Q), where
Q=C(G,CH-f)*C(N-F,256h-CH)/C(8192,h)^255 and CH=F-R+128h.
Q already contains the H union. No further H or state-support union is
used. The output-weight constraint is not needed by this broader bound.
The production result is FLOAT, not a directed certificate.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, lgamma, log
from operator import index
import os
import shutil
import threading
import time

import numpy as np


SR_WIDTH = 512
SR_UPDATES = 4096
SR_N = SR_WIDTH * SR_UPDATES
SR_M = 8192
SR_LD = np.longdouble
SR_LN2 = np.log(SR_LD(2))
SR_PROFILES = ((1536, 753664, 730726, 4, 851968),
               (2048, 819200, 819200, 4, 950272),
               (2560, 819200, 829030, 4, 983040))


def _sr_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_selector_temporal
    return getattr(quartet_selector_temporal, name)


def _sr_normalize(array):
    maximum = np.max(array)
    if not np.isfinite(maximum) or maximum <= 0 or not np.all(np.isfinite(array)) or np.any(array <= 0):
        raise ArithmeticError('A positive selector-rank action became invalid')
    return array / maximum, np.log2(maximum)


def _sr_logchoose(n, k):
    if not 0 <= k <= n:
        return SR_LD('-inf')
    if k in (0, n):
        return SR_LD(0)
    return SR_LD((lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / log(2))


def selector_rank_prepare(log_beta=1, width=SR_WIDTH):
    """Build overlap once; reuse this object across volume markers and theta."""
    width, log_beta = index(width), SR_LD(log_beta)
    if width < 4 or width % 4 or not np.isfinite(log_beta) or log_beta <= 0:
        raise ValueError('Require quartet width and finite positive log_beta')
    overlap, scale = _sr_dependency('temporal_selector_overlap')(log_beta, width)
    return dict(width=width, log_beta=log_beta, overlap=overlap, overlap_scale=SR_LD(scale))


def selector_rank_moment(log_beta=1, log_z=-1, width=SR_WIDTH,
                         updates=SR_UPDATES, prepared=None):
    """FLOAT moment with beta^width endpoint upper; independent of profile/theta."""
    width, updates = index(width), index(updates)
    log_beta, log_z = SR_LD(log_beta), SR_LD(log_z)
    if updates < 1 or not np.isfinite(log_z):
        raise ValueError('Require updates>=1 and finite log_z')
    prepared = selector_rank_prepare(log_beta, width) if prepared is None else prepared
    if prepared['width'] != width or prepared['log_beta'] != log_beta:
        raise ValueError('Prepared overlap parameters do not match')
    groups = width // 4
    logs = np.asarray([np.log2(SR_LD(comb(groups, j))) + 4 * j * log_z
                       for j in range(groups + 1)], dtype=SR_LD)
    weight_scale = logs.max()
    weights = np.exp2(logs - weight_scale)
    if np.any(weights <= 0) or not np.all(np.isfinite(weights)):
        raise ArithmeticError('Positive selector weight underflow or overflow')
    kernel, extra = _sr_normalize(prepared['overlap'] * weights[None, :])
    scale = prepared['overlap_scale'] + weight_scale + extra
    bulk = SR_LD(_sr_dependency('_st_power_row')(kernel, scale, updates - 1))
    endpoint = width * log_beta
    value = bulk + endpoint
    if not np.isfinite(value):
        raise ArithmeticError('Nonfinite selector moment')
    return dict(log2_moment=value, log2_free_endpoint_sum=bulk,
                log2_endpoint_upper=endpoint, log_beta=log_beta, log_z=log_z,
                width=width, updates=updates, states=groups + 1,
                endpoint='beta^width upper bound', floating_diagnostic_only=True)


def selector_rank_bound(profile=SR_PROFILES[1], log_beta=1, log_z=-1,
                        theta=0, prepared=None, moment=None):
    """One fixed-profile diagnostic; a supplied moment can be reused across theta."""
    profile = tuple(map(index, profile))
    h, rank, volume, births, footprint = profile
    log_beta, log_z, theta = map(SR_LD, (log_beta, log_z, theta))
    if not (1 <= h <= SR_M and 128 * h <= rank <= (1 << 20)
            and 1 <= births <= volume <= footprint <= SR_N):
        raise ValueError('Invalid nonzero actual-state profile')
    if not all(np.isfinite(value) for value in (log_beta, log_z, theta)) or log_beta <= 0 or not 0 <= theta <= 1:
        raise ValueError('Require log_beta>0, finite log_z, and0<=theta<=1')
    G, m = footprint - births, volume - births
    if G % 4 or G > SR_N - SR_WIDTH:
        raise ValueError('Actual G is quartet-aligned and has an empty first layer')
    covered_H = footprint - rank + 128 * h
    route = (_sr_logchoose(G, covered_H - births)
             + _sr_logchoose(SR_N - footprint, 256 * h - covered_H)
             - 255 * _sr_logchoose(SR_M, h))
    if np.isneginf(route):
        return dict(profile=profile, geometrically_empty=True,
                    raw_log2_probability_bound=float('-inf'),
                    floating_diagnostic_only=True, full_code_certificate=False)
    if not np.isfinite(route):
        raise ArithmeticError('Nonfinite routing bound')
    if moment is None:
        moment = selector_rank_moment(log_beta, log_z, prepared=prepared)
    if any(moment[key] != value for key, value in
           (('log_beta', log_beta), ('log_z', log_z), ('width', SR_WIDTH), ('updates', SR_UPDATES))):
        raise ValueError('Supplied moment parameters do not match')
    birth_count = _sr_logchoose(SR_N - G, births)
    route_used = min(SR_LD(0), route)
    field = theta * (32 * (volume - rank) - 30)
    if theta == 0:
        log_y = None
        support = -m * log_beta
    else:
        exponent = log_beta / theta
        log_y = exponent + np.log2(-np.expm1(-exponent * SR_LN2))
        support = -theta * m * log_y
    volume_marker = -G * log_z
    bound = SR_LD(moment['log2_moment']) + birth_count + route_used + field + support + volume_marker
    if not np.isfinite(bound):
        raise ArithmeticError('Nonfinite fixed-profile bound')
    return dict(profile=profile, G=G, actual_nonbirth_volume=m,
                log_beta=float(log_beta), log_z=float(log_z), theta=float(theta),
                log_y=None if log_y is None else float(log_y),
                geometry_log2_moment=float(moment['log2_moment']),
                free_endpoint_log2_upper=float(moment['log2_endpoint_upper']),
                birth_mask_log2_count=float(birth_count),
                routing_log2_uncapped=float(route), routing_log2_used=float(route_used),
                field_log2_correction=float(field), actual_volume_log2_correction=float(support),
                G_log2_marker_correction=float(volume_marker),
                raw_log2_probability_bound=float(bound),
                capped_log2_probability_bound=float(min(SR_LD(0), bound)),
                grouped_amplitudes_before_theta=True, routing_contains_H_union=True,
                additional_H_union=False, additional_state_support_union=False,
                output_weight_constraint_used=False, floating_diagnostic_only=True,
                full_code_certificate=False)


def _sr_exact_overlap(width, old_groups, new_groups, beta):
    old = (1 << (4 * old_groups)) - 1
    return sum((beta ** (old & sum(1 << i for i in subset)).bit_count()
                for subset in combinations(range(width), 4 * new_groups)), Fraction(0)) / comb(width, 4 * new_groups)


def _sr_exact_selector_row(width, updates, beta, z):
    groups = width // 4
    kernel = [[_sr_exact_overlap(width, i, j, beta) * comb(groups, j) * z ** (4 * j)
               for j in range(groups + 1)] for i in range(groups + 1)]
    row = [Fraction(int(i == 0)) for i in range(groups + 1)]
    for _ in range(updates - 1):
        row = [sum(row[i] * kernel[i][j] for i in range(groups + 1))
               for j in range(groups + 1)]
    return row, kernel


def selector_rank_checks():
    fibre_checks = binomial_checks = selector_checks = float_checks = 0
    for q in (2, 4):
        for allowed in range(6):
            counts = [0] * (allowed + 1)
            for values in product(range(q), repeat=allowed):
                counts[sum(value != 0 for value in values)] += 1
            for m, count in enumerate(counts):
                assert count == comb(allowed, m) * (q - 1) ** m
                fibre_checks += 1
    for allowed in range(9):
        for m in range(allowed + 1):
            for y in (Fraction(1, 2), Fraction(1), Fraction(3)):
                assert comb(allowed, m) * y ** m <= (1 + y) ** allowed
                assert Fraction(2) ** (allowed - m) >= 1  # theta0 indicator bound
                for constant in (Fraction(1, 16), Fraction(1), Fraction(16)):
                    value = constant * comb(allowed, m)
                    assert min(Fraction(1), value) ** 2 <= value  # theta1/2, after grouping
                binomial_checks += 1
    beta, z = Fraction(3, 2), Fraction(2, 3)
    for width in (4, 8):
        groups = width // 4
        first_partition = (15,) if width == 4 else (15, 240)
        partitions = [(15,)] if width == 4 else [
            (first, 255 ^ first) for rest in combinations(range(1, 8), 3)
            for first in (1 + sum(1 << i for i in rest),)]
        choices = lambda partition: [sum(block for j, block in enumerate(partition) if selected >> j & 1)
                                      for selected in range(1 << groups)]
        physical = Fraction(0)
        for next_partition in partitions:
            for old, new in product(choices(first_partition), choices(next_partition)):
                G = old.bit_count() + new.bit_count()
                A0 = (old & new).bit_count() + new.bit_count()
                physical += z ** G * beta ** A0 / len(partitions)
                allowed_births = old.bit_count() + (new & ~old).bit_count() + (width - new.bit_count())
                for f in range(4):
                    assert (comb(allowed_births, f) if f <= allowed_births else 0) <= comb(3 * width - G, f)
                selector_checks += 1
        row, _ = _sr_exact_selector_row(width, 3, beta, z)
        exact = sum(value * beta ** (4 * i) for i, value in enumerate(row))
        assert exact == physical
        for updates in (1, 2, 3):
            row, _ = _sr_exact_selector_row(width, updates, beta, z)
            exact = sum(value * beta ** (4 * i) for i, value in enumerate(row))
            upper = beta ** width * sum(row)
            assert exact <= upper
            if updates == 1:
                assert exact == 1
            if updates == 2:
                assert exact == (1 + (beta * z) ** 4) ** groups
            computed = selector_rank_moment(np.log2(SR_LD(3) / 2), np.log2(SR_LD(2) / 3), width, updates)
            reference = np.log2(SR_LD(upper.numerator) / upper.denominator)
            assert abs(computed['log2_moment'] - reference) < SR_LD('1e-10')
            float_checks += 1
        # Independent exact scalar witness for beta4,z1/4.
        _, kernel = _sr_exact_selector_row(width, 1, Fraction(4), Fraction(1, 4))
        witness = [Fraction(16) ** i for i in range(groups + 1)]
        for i in range(groups + 1):
            assert sum(kernel[i][j] * witness[j] for j in range(groups + 1)) <= Fraction(17, 16) ** groups * witness[i]
            selector_checks += 1
    print(dict(exact_fixed_volume_fibres=fibre_checks, exact_binomial_grouping_checks=binomial_checks,
               exact_consistent_selector_checks=selector_checks, floating_moment_checks=float_checks,
               free_endpoint_checked=True, no_extra_H_union=True,
               full_code_certificate=False), flush=True)


def selector_rank_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=SR_PROFILES[1])
    parser.add_argument('--log-beta', type=float, default=1)
    parser.add_argument('--log-z', type=float, default=-1)
    parser.add_argument('--theta', type=float, default=0)
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
            print(selector_rank_bound(args.profile, args.log_beta, args.log_z, args.theta), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    selector_rank_main()
