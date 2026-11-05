"""Positive selector-hit and regional-capacity quotient diagnostics.

For selected quartet unions U_t, let J0 count selected quartets at t
disjoint from U_(t-1).  Every such actual live quartet needs a distinct
forced birth at t-1, so J0<=f.  With0<u<=1, insert u^J0 and correct by
u^-f.  U_1 is empty; the formal full endpoint has no quartet-hit charge.

Given an old union of4i positions, the exact weighted selector kernel is

 K_ij=C(m,j) z^(4j)/C(n,4i) * [X^(4i)]
       (u+4 beta X+6 beta²X²+4 beta³X³+beta⁴X⁴)^j
       (1+X)^(n-4j),  m=n/4.

The prepared positive coefficient matrix is compatible with the frozen
quartet_selector_rank_moment APIs.  Production arithmetic is FLOAT and
does not certify rounding, even when all positive checks finish.

For regional live volumes G_r, matching routes necessarily satisfy
C_H-f<=sum min(h,G_r).  A c>=1 marker is bounded entrywise perregion by
min(K(z*c)^s,c^h*K(z)^s), with s=15 first and16 thereafter.  Compose the
first block and255 later blocks, retain beta^512 as a free-endpoint upper,
and correct by c^(-(C_H-f)).  The route factor is still the frozen one;
it already includes the H union.  No added route or state-support union.

Inject quartet_selector_temporal.py and quartet_selector_rank_moment.py
definitions, or import this helper beside them.  No frozen source is
modified.  Parent-only jobs, with no files, subprocesses or witness dumps.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np


SH_WIDTH = 512
SH_UPDATES = 4096
SH_REGION_UPDATES = 16
SH_LD = np.longdouble
SH_PROFILE = (2560, 786432, 779879, 4, 950272)


def _sh_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_selector_rank_moment
    return getattr(quartet_selector_rank_moment, name)


def _sh_normalize(array):
    maximum = np.max(array)
    if maximum <= 0 or not np.isfinite(maximum) or np.any(array <= 0) or not np.all(np.isfinite(array)):
        raise ArithmeticError('Positive selector-hit action underflow or overflow')
    return array / maximum, np.log2(maximum)


def selector_hit_prepare(log_beta=1, log_u=-8, width=SH_WIDTH):
    """Build the129-population positive coefficient overlap once."""
    width = index(width)
    log_beta, log_u = SH_LD(log_beta), SH_LD(log_u)
    if width < 4 or width % 4 or not np.isfinite(log_beta) or log_beta <= 0:
        raise ValueError('Require quartet width and positive finite log_beta')
    if not np.isfinite(log_u) or log_u > 0:
        raise ValueError('Require finite log_u<=0, hence0<u<=1')
    groups = width // 4
    local = np.array([np.exp2(log_u), 4 * np.exp2(log_beta),
                      6 * np.exp2(2 * log_beta), 4 * np.exp2(3 * log_beta),
                      np.exp2(4 * log_beta)], dtype=SH_LD)
    local, local_scale = _sh_normalize(local)
    old_indices = np.arange(0, width + 1, 4)
    old_log_denominators = np.array([np.log2(SH_LD(comb(width, int(i)))) for i in old_indices], dtype=SH_LD)
    log_matrix = np.empty((groups + 1, groups + 1), dtype=SH_LD)
    selected = np.ones(1, dtype=SH_LD)
    selected_scale = SH_LD(0)
    for j in range(groups + 1):
        remaining = width - 4 * j
        unselected = np.array([SH_LD(comb(remaining, k)) for k in range(remaining + 1)], dtype=SH_LD)
        unselected, unselected_scale = _sh_normalize(unselected)
        coefficients = np.convolve(selected, unselected)[old_indices]
        if np.any(coefficients <= 0) or not np.all(np.isfinite(coefficients)):
            raise ArithmeticError('Positive selector-hit coefficient became invalid')
        log_matrix[:, j] = (np.log2(coefficients) + selected_scale + unselected_scale
                            - old_log_denominators)
        if j < groups:
            selected, extra = _sh_normalize(np.convolve(selected, local))
            selected_scale += local_scale + extra
    scale = np.max(log_matrix)
    overlap = np.exp2(log_matrix - scale)
    if np.any(overlap <= 0) or not np.all(np.isfinite(overlap)):
        raise ArithmeticError('Positive selector-hit overlap became invalid')
    return dict(width=width, log_beta=log_beta, log_u=log_u,
                overlap=overlap, overlap_scale=scale,
                prepared_method='positive selected-quartet coefficient convolutions')


def _sh_validate_prepared(prepared, log_beta, log_u, width):
    if any(prepared[key] != value for key, value in
           (('width', width), ('log_beta', log_beta), ('log_u', log_u))):
        raise ValueError('Prepared selector-hit parameters do not match')


def selector_hit_moment(log_beta=1, log_z=-1, log_u=-8, width=SH_WIDTH,
                        updates=SH_UPDATES, prepared=None):
    """Reuse frozen selector volume weights and positive powers, with the hit overlap."""
    width, updates = index(width), index(updates)
    log_beta, log_z, log_u = map(SH_LD, (log_beta, log_z, log_u))
    prepared = selector_hit_prepare(log_beta, log_u, width) if prepared is None else prepared
    _sh_validate_prepared(prepared, log_beta, log_u, width)
    moment = _sh_dependency('selector_rank_moment')(log_beta, log_z, width, updates, prepared)
    return dict(**moment, log_u=log_u, log_c=SH_LD(0), hit_marker=True,
                regional_capacity_marker=False)


def _sh_kernel(prepared, log_z):
    width = prepared['width']
    groups = width // 4
    logs = np.array([np.log2(SH_LD(comb(groups, j))) + 4 * j * log_z for j in range(groups + 1)], dtype=SH_LD)
    weight_scale = np.max(logs)
    weights = np.exp2(logs - weight_scale)
    if np.any(weights <= 0) or not np.all(np.isfinite(weights)):
        raise ArithmeticError('Positive selector volume weight became invalid')
    kernel, extra = _sh_normalize(prepared['overlap'] * weights[None, :])
    return kernel, prepared['overlap_scale'] + weight_scale + extra


def _sh_matrix_power(kernel, scale, count):
    count = index(count)
    if count < 1:
        raise ValueError('Regional matrix powers must be positive')
    result, result_scale = None, SH_LD(0)
    power, power_scale = kernel.copy(), SH_LD(scale)
    while count:
        if count & 1:
            if result is None:
                result, result_scale = power.copy(), power_scale
            else:
                result, extra = _sh_normalize(result @ power)
                result_scale += power_scale + extra
        count >>= 1
        if count:
            power, extra = _sh_normalize(power @ power)
            power_scale = 2 * power_scale + extra
    return result, result_scale


def _sh_minimum(first, first_scale, second, second_scale):
    logs = np.minimum(np.log2(first) + first_scale, np.log2(second) + second_scale)
    scale = np.max(logs)
    values = np.exp2(logs - scale)
    if np.any(values <= 0) or not np.all(np.isfinite(values)):
        raise ArithmeticError('Positive regional minimum became invalid')
    return values, scale


def _sh_row_power(row, row_scale, kernel, kernel_scale, count):
    power, power_scale = kernel.copy(), SH_LD(kernel_scale)
    count = index(count)
    while count:
        if count & 1:
            row, extra = _sh_normalize(row @ power)
            row_scale += power_scale + extra
        count >>= 1
        if count:
            power, extra = _sh_normalize(power @ power)
            power_scale = 2 * power_scale + extra
    return row_scale + np.log2(np.sum(row))


def selector_hit_regional_moment(h, log_beta=1, log_z=-1, log_u=-8, log_c=1 / 16,
                                 width=SH_WIDTH, updates=SH_UPDATES,
                                 region_updates=SH_REGION_UPDATES, prepared=None):
    """First region has one fewer random selector layer; all later regions are full."""
    h, width, updates, region_updates = map(index, (h, width, updates, region_updates))
    log_beta, log_z, log_u, log_c = map(SH_LD, (log_beta, log_z, log_u, log_c))
    if region_updates < 2 or updates < region_updates or updates % region_updates:
        raise ValueError('Require whole regions of at least two updates')
    if not 0 <= h <= width * region_updates or not np.isfinite(log_c) or log_c < 0 or not np.isfinite(log_z):
        raise ValueError('Invalid regional capacity or marker')
    prepared = selector_hit_prepare(log_beta, log_u, width) if prepared is None else prepared
    _sh_validate_prepared(prepared, log_beta, log_u, width)
    plain, plain_scale = _sh_kernel(prepared, log_z)
    marked, marked_scale = _sh_kernel(prepared, log_z + log_c)
    blocks = []
    for count in (region_updates - 1, region_updates):
        first, first_scale = _sh_matrix_power(marked, marked_scale, count)
        second, second_scale = _sh_matrix_power(plain, plain_scale, count)
        blocks.append(_sh_minimum(first, first_scale, second, second_scale + h * log_c))
    initial, initial_scale = blocks[0]
    regional, regional_scale = blocks[1]
    bulk = _sh_row_power(initial[0].copy(), initial_scale, regional, regional_scale,
                         updates // region_updates - 1)
    endpoint = width * log_beta
    value = bulk + endpoint
    if not np.isfinite(value):
        raise ArithmeticError('Nonfinite regional selector-hit moment')
    return dict(log2_moment=value, log2_free_endpoint_sum=bulk,
                log2_endpoint_upper=endpoint, log_beta=log_beta, log_z=log_z,
                log_u=log_u, log_c=log_c, width=width, updates=updates,
                states=width // 4 + 1, endpoint='beta^width upper bound',
                hit_marker=True, regional_capacity_marker=True, capacity_h=h,
                first_region_selector_layers=region_updates - 1,
                later_region_selector_layers=region_updates,
                regions=updates // region_updates, floating_diagnostic_only=True)


def selector_hit_bound(profile=SH_PROFILE, log_beta=1, log_z=-1, log_u=-8,
                       theta=0, log_c=0, prepared=None, moment=None):
    """Apply both necessary-geometry markers to the frozen actual-L profile bound."""
    profile = tuple(map(index, profile))
    h, rank, volume, births, footprint = profile
    log_beta, log_z, log_u, theta, log_c = map(SH_LD, (log_beta, log_z, log_u, theta, log_c))
    if not np.isfinite(log_u) or log_u > 0 or not np.isfinite(log_c) or log_c < 0:
        raise ValueError('Require log_u<=0 and log_c>=0, both finite')
    if moment is None:
        moment = (selector_hit_regional_moment(h, log_beta, log_z, log_u, log_c, prepared=prepared)
                  if log_c > 0 else selector_hit_moment(log_beta, log_z, log_u, prepared=prepared))
    if moment.get('log_u') != log_u or moment.get('log_c') != log_c:
        raise ValueError('Supplied marker moment parameters do not match')
    if moment.get('regional_capacity_marker') and moment.get('capacity_h') != h:
        raise ValueError('Supplied regional moment has a different capacity')
    base = _sh_dependency('selector_rank_bound')(profile, log_beta, log_z, theta, prepared, moment)
    if base.get('geometrically_empty'):
        return base
    covered_H = footprint - rank + 128 * h
    hit_correction = -births * log_u
    capacity_correction = -(covered_H - births) * log_c
    value = SH_LD(base['raw_log2_probability_bound']) + hit_correction + capacity_correction
    if not np.isfinite(value):
        raise ArithmeticError('Nonfinite selector-hit profile bound')
    base.update(log_u=float(log_u), log_c=float(log_c),
                hit_marker_log2_correction=float(hit_correction),
                regional_capacity_log2_correction=float(capacity_correction),
                regional_capacity_marker=bool(log_c > 0),
                raw_log2_probability_bound=float(value),
                capped_log2_probability_bound=float(min(SH_LD(0), value)),
                selected_empty_quartets_charged_to_births=True,
                full_code_certificate=False)
    return base


def _sh_partitions(width):
    if width == 4:
        return ((15,),)
    if width == 8:
        return tuple((1 + sum(1 << i for i in rest), 255 ^ (1 + sum(1 << i for i in rest)))
                     for rest in combinations(range(1, 8), 3))
    raise ValueError('Exact toys use width4 or8')


def _sh_convolve(first, second):
    result = [Fraction(0)] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return result


def _sh_exact_kernel(width, beta, z, u):
    groups = width // 4
    local = [u] + [comb(4, c) * beta ** c for c in range(1, 5)]
    selected = [Fraction(1)]
    kernel = [[Fraction(0)] * (groups + 1) for _ in range(groups + 1)]
    for j in range(groups + 1):
        polynomial = _sh_convolve(selected, [Fraction(comb(width - 4 * j, k)) for k in range(width - 4 * j + 1)])
        for i in range(groups + 1):
            kernel[i][j] = comb(groups, j) * z ** (4 * j) * polynomial[4 * i] / comb(width, 4 * i)
        selected = _sh_convolve(selected, local)
    return kernel


def _sh_exact_power(matrix, exponent):
    size = len(matrix)
    result = [[Fraction(int(i == j)) for j in range(size)] for i in range(size)]
    for _ in range(exponent):
        result = [[sum(result[i][k] * matrix[k][j] for k in range(size))
                   for j in range(size)] for i in range(size)]
    return result


def selector_hit_checks():
    physical_checks = birth_checks = float_checks = regional_checks = 0
    beta, z, u, c = Fraction(3, 2), Fraction(2, 3), Fraction(1, 4), Fraction(3, 2)
    for width in (4, 8):
        groups = width // 4
        partitions = _sh_partitions(width)
        kernel = _sh_exact_kernel(width, beta, z, u)
        prepared = selector_hit_prepare(log2(float(beta)), log2(float(u)), width)
        reconstructed = prepared['overlap'] * np.exp2(prepared['overlap_scale'])
        for i in range(groups + 1):
            old = (1 << (4 * i)) - 1
            physical = [Fraction(0)] * (groups + 1)
            for partition in partitions:
                for selector in range(1 << groups):
                    blocks = [quartet for j, quartet in enumerate(partition) if (selector >> j) & 1]
                    new = sum(blocks)
                    empty = sum(not (quartet & old) for quartet in blocks)
                    physical[len(blocks)] += z ** new.bit_count() * beta ** (new & old).bit_count() * u ** empty / len(partitions)
                for state in range(1 << width):
                    live = [quartet for quartet in partition if quartet & state]
                    empty = sum(not (quartet & old) for quartet in live)
                    births = (state & ~old).bit_count()
                    assert empty <= births and u ** empty >= u ** births
                    birth_checks += 1
            for j in range(groups + 1):
                assert physical[j] == kernel[i][j]
                reference = kernel[i][j] / (comb(groups, j) * z ** (4 * j))
                assert abs(reconstructed[i, j] - SH_LD(reference.numerator) / reference.denominator) < SH_LD('1e-10')
                physical_checks += 1
        for updates in (1, 2, 3):
            matrix = _sh_exact_power(kernel, updates - 1)
            row = matrix[0]
            exact = sum(row[i] * beta ** (4 * i) for i in range(groups + 1))
            upper = beta ** width * sum(row)
            assert exact <= upper
            computed = selector_hit_moment(log2(float(beta)), log2(float(z)), log2(float(u)), width, updates, prepared)
            assert abs(computed['log2_moment'] - np.log2(SH_LD(upper.numerator) / upper.denominator)) < SH_LD('1e-10')
            float_checks += 1
        # Exact regional matrix minima, plus the genuinely capped path sum.
        capacity = 3
        marked = _sh_exact_kernel(width, beta, z * c, u)
        comparisons = []
        for layers in (1, 2):
            plain_power, marked_power = _sh_exact_power(kernel, layers), _sh_exact_power(marked, layers)
            comparisons.append([[min(marked_power[i][j], c ** capacity * plain_power[i][j])
                                 for j in range(groups + 1)] for i in range(groups + 1)])
        upper_row = [sum(comparisons[0][0][i] * comparisons[1][i][j] for i in range(groups + 1))
                     for j in range(groups + 1)]
        matrix_upper = beta ** width * sum(upper_row)
        exact_cap = Fraction(0)
        for first, second, third in product(range(groups + 1), repeat=3):
            path = kernel[0][first] * kernel[first][second] * kernel[second][third]
            regional_volume = min(capacity, 4 * first) + min(capacity, 4 * (second + third))
            exact_cap += path * c ** regional_volume * beta ** (4 * third)
        assert exact_cap <= matrix_upper
        computed = selector_hit_regional_moment(capacity, log2(float(beta)), log2(float(z)),
                                                log2(float(u)), log2(float(c)), width, 4, 2, prepared)
        assert abs(computed['log2_moment'] - np.log2(SH_LD(matrix_upper.numerator) / matrix_upper.denominator)) < SH_LD('1e-10')
        regional_checks += 1
        neutral = selector_hit_prepare(log2(float(beta)), 0, width)
        neutral_kernel = _sh_exact_kernel(width, beta, z, Fraction(1))
        for i, j in product(range(groups + 1), repeat=2):
            reference = neutral_kernel[i][j] / (comb(groups, j) * z ** (4 * j))
            value = neutral['overlap'][i, j] * np.exp2(neutral['overlap_scale'])
            assert abs(value - SH_LD(reference.numerator) / reference.denominator) < SH_LD('1e-10')
            float_checks += 1
    print(dict(exact_physical_kernel_checks=physical_checks, exact_birth_hit_checks=birth_checks,
               floating_moment_and_neutral_checks=float_checks,
               exact_and_float_regional_cap_checks=regional_checks,
               terminal_has_no_hit_or_volume_charge=True,
               full_code_certificate=False), flush=True)


def selector_hit_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=SH_PROFILE)
    parser.add_argument('--log-beta', type=float, default=1)
    parser.add_argument('--log-z', type=float, default=-1)
    parser.add_argument('--log-u', type=float, default=-8)
    parser.add_argument('--theta', type=float, default=0)
    parser.add_argument('--log-c', type=float, default=0, help='Positive value enables the regional capacity option')
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
            selector_hit_checks()
        else:
            print(selector_hit_bound(args.profile, args.log_beta, args.log_z, args.log_u,
                                     args.theta, args.log_c), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    selector_hit_main()
