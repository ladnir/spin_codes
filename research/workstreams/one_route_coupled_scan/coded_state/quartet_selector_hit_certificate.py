"""Exact dyadic regional certificate for consistent, birth-supported selectors.

Selected new quartets disjoint from the previous live set each require a
distinct forced birth. Their total is at most f. For0<u<=1, marking each
such quartet by u and paying u^-f therefore gives an upper bound.

The selected-quartet polynomial is positive:
P(x)=u+4*beta*x+6*beta^2*x^2+4*beta^3*x^3+beta^4*x^4.
In the fixed similarity basis D_ii=16^i,
K'_ij=C(m,j) z^(4j) [x^(4i)]P(x)^j(1+x)^(n-4j)
      /C(n,4i) *16^(j-i), m=n/4.

Each exact entry is rounded UP to a multiple2^-bits. The production default
is512 bits; coarser controls remained valid but added too much positive
mass to tiny transitions. Object-array integer products are normalized UP
to the same global mantissa precision. Dyadic operands
are aligned exactly before an entrywise regional minimum. Structural
zeros remain zero; positive(0,0) excursions are never removed.

The regional capacity gives the upper operator
min(K'(zc)^ell,c^h*K'(z)^ell), ell15 or16, in this SAME basis.
The first region has15 selector transitions; the other255 have16.
The free endpoint is bounded by max_i(beta^4/16)^i. All H-set accounting
is already in the capped route factor. This certifies selected whole
actual-profile classes, not the full distance goal. No FLOAT arithmetic
is used in the certificate path. Parent-only, serial, guarded jobs.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb
from operator import index
import os
import shutil
import threading
import time

import numpy as np


SH_BITS = 512
SH_SCALE = 1 << 48
SH_WIDTH = 512
SH_REGION_UPDATES = 16
SH_REGIONS = 256
SH_BETA = Fraction(35, 16)
SH_Z = Fraction(459, 1024)
SH_C = Fraction(33, 32)
SH_U = Fraction(1, 256)
SH_PROFILE = (2048, 819200, 819200, 4, 950272)


def _sh_ceildiv(a, b):
    if b <= 0:
        raise ValueError('Positive denominator required')
    return -((-a) // b)


def _sh_dyadic_denominator(value):
    value = Fraction(value)
    denominator = value.denominator
    if denominator & (denominator - 1):
        raise ValueError('Certificate markers must be exact dyadic rationals')
    return value, denominator.bit_length() - 1


def _sh_normalize(array, exponent, bits=SH_BITS):
    """Array*2^exponent, rounded upward to at most bits global mantissa bits."""
    maximum = int(np.max(array))
    if maximum < 0 or np.any(array < 0):
        raise ArithmeticError('A positive dyadic action became negative')
    while maximum.bit_length() > bits:
        shift = maximum.bit_length() - bits
        divisor = 1 << shift
        array = (array + divisor - 1) // divisor
        exponent += shift
        maximum = int(np.max(array))
    return array, exponent


def _sh_multiply(first, second, bits=SH_BITS):
    return _sh_normalize(first[0] @ second[0], first[1] + second[1], bits)


def _sh_scalar_upper(value, bits=SH_BITS):
    """Positive rational -> positive dyadic upper with at most bits bits."""
    value = Fraction(value)
    if value <= 0:
        raise ValueError('Positive rational scalar required')
    n, d = value.numerator, value.denominator
    floor_log = n.bit_length() - d.bit_length()
    if floor_log >= 0:
        if n < d << floor_log:
            floor_log -= 1
    elif n << (-floor_log) < d:
        floor_log -= 1
    exponent = floor_log - bits + 1
    mantissa = (_sh_ceildiv(n, d << exponent) if exponent >= 0
                else _sh_ceildiv(n << (-exponent), d))
    while mantissa.bit_length() > bits:
        mantissa = (mantissa + 1) // 2
        exponent += 1
    return mantissa, exponent


def _sh_scale(first, scalar, bits=SH_BITS):
    return _sh_normalize(first[0] * scalar[0], first[1] + scalar[1], bits)


def _sh_minimum(first, second, bits=SH_BITS):
    exponent = min(first[1], second[1])
    left = first[0] * (1 << (first[1] - exponent))
    right = second[0] * (1 << (second[1] - exponent))
    return _sh_normalize(np.minimum(left, right), exponent, bits)


def _sh_power(matrix, power, bits=SH_BITS):
    power = index(power)
    if power < 0:
        raise ValueError('Nonnegative power required')
    result = (np.eye(matrix[0].shape[0], dtype=object), 0)
    current = matrix
    while power:
        if power & 1:
            result = _sh_multiply(result, current, bits)
        power >>= 1
        if power:
            current = _sh_multiply(current, current, bits)
    return result


def _sh_adjacent_powers(matrix, length, bits=SH_BITS):
    if length != 16:
        return _sh_power(matrix, length - 1, bits), _sh_power(matrix, length, bits)
    square = _sh_multiply(matrix, matrix, bits)
    fourth = _sh_multiply(square, square, bits)
    eighth = _sh_multiply(fourth, fourth, bits)
    sixteenth = _sh_multiply(eighth, eighth, bits)
    fifteenth = _sh_multiply(_sh_multiply(_sh_multiply(eighth, fourth, bits), square, bits), matrix, bits)
    return fifteenth, sixteenth


def selector_hit_kernels(width, beta, volumes, u, bits=SH_BITS, progress=False):
    """One exact positive convolution pass builds every requested volume tilt."""
    width, bits = index(width), index(bits)
    beta, _ = _sh_dyadic_denominator(beta)
    u, _ = _sh_dyadic_denominator(u)
    volumes = tuple(_sh_dyadic_denominator(value) for value in volumes)
    if width < 4 or width % 4 or bits < 2 or beta < 1 or not 0 < u <= 1 or any(value <= 0 for value, _ in volumes):
        raise ValueError('Invalid positive dyadic kernel parameters')
    groups = width // 4
    coefficients = (u, 4 * beta, 6 * beta ** 2, 4 * beta ** 3, beta ** 4)
    denominator = max(value.denominator for value in coefficients)
    denominator_bits = denominator.bit_length() - 1
    polynomial = np.asarray([int(value * denominator) for value in coefficients], dtype=object)
    assert all(Fraction(value, denominator) == exact for value, exact in zip(polynomial, coefficients))
    matrices = [np.zeros((groups + 1, groups + 1), dtype=object) for _ in volumes]
    choose_old = [comb(width, 4 * i) for i in range(groups + 1)]
    polynomial_power = np.asarray([1], dtype=object)
    for j in range(groups + 1):
        remainder = width - 4 * j
        binomial = np.asarray([comb(remainder, k) for k in range(remainder + 1)], dtype=object)
        convolution = np.convolve(polynomial_power, binomial)
        group_count = comb(groups, j)
        for matrix, (volume, volume_den_bits) in zip(matrices, volumes):
            column_factor = group_count * volume.numerator ** (4 * j)
            common_shift = bits + 4 * j - denominator_bits * j - volume_den_bits * 4 * j
            for i in range(groups + 1):
                numerator = int(convolution[4 * i]) * column_factor
                denominator_i = choose_old[i]
                shift = common_shift - 4 * i
                if shift >= 0:
                    numerator <<= shift
                else:
                    denominator_i <<= -shift
                matrix[i, j] = _sh_ceildiv(numerator, denominator_i)
        if progress and (j % 32 == 0 or j == groups):
            print(dict(phase='exact_kernel_columns', completed=j, total=groups), flush=True)
        if j != groups:
            # Fixed quartic shifts; each coefficient and intermediate is positive.
            following = np.zeros(len(polynomial_power) + 4, dtype=object)
            for shift in range(5):
                following[shift:shift + len(polynomial_power)] += polynomial[shift] * polynomial_power
            polynomial_power = following
    return tuple((matrix, -bits) for matrix in matrices)


def selector_hit_prepare(width=SH_WIDTH, region_updates=SH_REGION_UPDATES,
                         beta=SH_BETA, z=SH_Z, c=SH_C, u=SH_U,
                         bits=SH_BITS, progress=False):
    """Reuse exact regional powers across profiles with different eligible counts."""
    width, region_updates, bits = map(index, (width, region_updates, bits))
    beta, z, c, u = map(Fraction, (beta, z, c, u))
    if region_updates < 1 or c < 1:
        raise ValueError('Require region_updates>=1 and c>=1')
    _sh_dyadic_denominator(c)
    plain, marked = selector_hit_kernels(width, beta, (z, z * c), u, bits, progress)
    plain_powers = _sh_adjacent_powers(plain, region_updates, bits)
    if progress:
        print(dict(phase='plain_regional_powers_complete'), flush=True)
    marked_powers = _sh_adjacent_powers(marked, region_updates, bits)
    if progress:
        print(dict(phase='marked_regional_powers_complete'), flush=True)
    return dict(width=width, region_updates=region_updates, bits=bits,
                beta=beta, z=z, c=c, u=u,
                plain=plain_powers, marked=marked_powers)


def selector_hit_regions(prepared, h):
    h = index(h)
    if h < 0:
        raise ValueError('Nonnegative regional capacity required')
    scalar = _sh_scalar_upper(prepared['c'] ** h, prepared['bits'])
    return tuple(_sh_minimum(marked, _sh_scale(plain, scalar, prepared['bits']), prepared['bits'])
                 for plain, marked in zip(prepared['plain'], prepared['marked']))


def selector_hit_moment(first, regular, regions=SH_REGIONS, bits=SH_BITS, progress=False):
    """Dyadic upper on e0*B_first*B_regular^(regions-1)*1."""
    regions = index(regions)
    if regions < 1:
        raise ValueError('Positive region count required')
    row = (first[0][0].copy(), first[1])
    power, count = regular, regions - 1
    while count:
        if count & 1:
            row = _sh_multiply(row, power, bits)
        count >>= 1
        if count:
            power = _sh_multiply(power, power, bits)
    total = sum(int(value) for value in row[0])
    if total <= 0:
        raise ArithmeticError('The selector moment must contain the empty history')
    if progress:
        print(dict(phase='global_region_powers_complete'), flush=True)
    return total, row[1]


def _sh_dependencies():
    names = ('selector_log2_interval', 'selector_rank_profile')
    if all(name in globals() for name in names):
        if globals().get('SC_SCALE') != SH_SCALE:
            raise RuntimeError('Directed log scale must be2^48')
        return globals()['selector_log2_interval'], globals()['selector_rank_profile']
    import quartet_selector_capacity as capacity
    import quartet_selector_rank_certificate as rank
    rank.SC_SCALE = capacity.SC_SCALE
    rank.selector_log2_interval = capacity.selector_log2_interval
    rank.selector_binomial_log_interval = capacity.selector_binomial_log_interval
    return capacity.selector_log2_interval, rank.selector_rank_profile


def selector_hit_profile(profile=SH_PROFILE, prepared=None, progress=False):
    """Exact upper exponent for one whole actual-profile class, all deaths."""
    log_interval, rank_profile = _sh_dependencies()
    profile = tuple(map(index, profile))
    base = rank_profile(*profile)
    if base['status'] == 'empty':
        return dict(profile=profile, status='empty', exact_integer_certificate=True,
                    full_code_certificate=False)
    prepared = selector_hit_prepare(progress=progress) if prepared is None else prepared
    if prepared['width'] != SH_WIDTH or prepared['region_updates'] != SH_REGION_UPDATES:
        raise ValueError('Production profile requires width512 and16-update regions')
    h, _, volume, births, _ = profile
    first, regular = selector_hit_regions(prepared, h)
    mantissa, exponent = selector_hit_moment(first, regular, SH_REGIONS, prepared['bits'], progress)
    moment_upper = log_interval(mantissa)[1] + exponent * SH_SCALE
    beta, z, c, u = (prepared[key] for key in ('beta', 'z', 'c', 'u'))
    endpoint_upper = 0
    if beta > 2:
        endpoint_upper = SH_WIDTH * log_interval(beta.numerator, 2 * beta.denominator)[1]
    moment_upper += endpoint_upper
    G, CH = base['G'], base['C_H']
    upper = (moment_upper + base['log_births_upper'] + base['log_route_capped_upper']
             - G * log_interval(z.numerator, z.denominator)[0]
             - (volume - births) * log_interval(beta.numerator, beta.denominator)[0]
             - births * log_interval(u.numerator, u.denominator)[0]
             - (CH - births) * log_interval(c.numerator, c.denominator)[0])
    integer_exponent = _sh_ceildiv(upper, SH_SCALE)
    return dict(profile=profile, status='whole_actual_profile',
                beta=str(beta), z=str(z), c=str(c), u=str(u),
                mantissa_bits=prepared['bits'], states=SH_WIDTH // 4 + 1,
                moment_log2_upper_scaled=moment_upper, endpoint_log2_upper_scaled=endpoint_upper,
                birth_log2_upper_scaled=base['log_births_upper'],
                route_log2_upper_scaled=base['log_route_capped_upper'],
                profile_log2_upper_scaled=upper, directed_log_scale=SH_SCALE,
                exact_profile_exponent_upper=integer_exponent,
                qualifying_threshold=-172, qualifies_global_profile_class=integer_exponent <= -172,
                global_profile_union_bits=96, class_failure_exponent=-76,
                all_death_counts=True, H_union_already_in_route=True,
                first_region_transitions=15, later_region_transitions=16,
                positive_return_excursions_preserved=True, free_terminal=True,
                exact_integer_certificate=True, full_code_certificate=False)


def _sh_fraction_array(pair):
    factor = Fraction(2) ** pair[1]
    return [[int(value) * factor for value in row] for row in pair[0]]


def _sh_exact_product(first, second):
    size = len(first)
    return [[sum(first[i][k] * second[k][j] for k in range(size))
             for j in range(size)] for i in range(size)]


def _sh_exact_power(matrix, power):
    result = [[Fraction(int(i == j)) for j in range(len(matrix))] for i in range(len(matrix))]
    for _ in range(power):
        result = _sh_exact_product(result, matrix)
    return result


def _sh_physical_kernel(width, beta, z, u):
    groups = width // 4
    partitions = [(15,)] if width == 4 else [
        (first, 255 ^ first) for rest in combinations(range(1, 8), 3)
        for first in (1 + sum(1 << i for i in rest),)]
    matrix = [[Fraction(0) for _ in range(groups + 1)] for _ in range(groups + 1)]
    for i in range(groups + 1):
        old = (1 << (4 * i)) - 1
        for partition in partitions:
            for selected in range(1 << groups):
                chosen = [block for j, block in enumerate(partition) if selected >> j & 1]
                new = sum(chosen)
                j = len(chosen)
                empty = sum(not (block & old) for block in chosen)
                matrix[i][j] += (z ** (4 * j) * beta ** (old & new).bit_count()
                                  * u ** empty * Fraction(16) ** (j - i) / len(partitions))
    return matrix


def selector_hit_checks():
    kernel_entries = product_entries = region_entries = moment_checks = 0
    beta, z, c, u = Fraction(3, 2), Fraction(3, 4), Fraction(5, 4), Fraction(1, 4)
    for width in (4, 8):
        size = width // 4 + 1
        exact = _sh_physical_kernel(width, beta, z, u)
        initial = selector_hit_kernels(width, beta, (z,), u)[0]
        upper = _sh_fraction_array(initial)
        for i, j in product(range(size), repeat=2):
            assert exact[i][j] <= upper[i][j] < exact[i][j] + Fraction(1, 1 << SH_BITS)
            kernel_entries += 1
        for bits in (8, SH_BITS):
            power = _sh_power(initial, 5, bits)
            exact_power = _sh_exact_power(exact, 5)
            power_upper = _sh_fraction_array(power)
            for i, j in product(range(size), repeat=2):
                assert exact_power[i][j] <= power_upper[i][j]
                product_entries += 1
            prepared = selector_hit_prepare(width, 2, beta, z, c, u, bits)
            blocks = selector_hit_regions(prepared, 3)
            exact_marked = _sh_physical_kernel(width, beta, z * c, u)
            exact_blocks = []
            for length, block in zip((1, 2), blocks):
                plain_power = _sh_exact_power(exact, length)
                marked_power = _sh_exact_power(exact_marked, length)
                reference = [[min(marked_power[i][j], c ** 3 * plain_power[i][j])
                              for j in range(size)] for i in range(size)]
                observed = _sh_fraction_array(block)
                for i, j in product(range(size), repeat=2):
                    assert reference[i][j] <= observed[i][j]
                    if length == 1:
                        capped = exact[i][j] * c ** min(3, 4 * j)
                    else:
                        capped = sum(exact[i][k] * exact[k][j] * c ** min(3, 4 * (k + j))
                                     for k in range(size))
                    assert capped <= reference[i][j]
                    region_entries += 1
                exact_blocks.append(reference)
            combined = _sh_exact_product(exact_blocks[0], _sh_exact_power(exact_blocks[1], 2))
            value, scale = selector_hit_moment(*blocks, regions=3, bits=bits)
            assert sum(combined[0]) <= value * Fraction(2) ** scale
            moment_checks += 1
            zero = np.asarray([[0, 1], [0, 0]], dtype=object)
            normalized, _ = _sh_normalize(zero * (1 << 100), -100, bits)
            assert np.array_equal(normalized == 0, zero == 0)
    for value in (Fraction(1, 3), Fraction(33, 32) ** 2048, Fraction(1 << 100)):
        mantissa, exponent = _sh_scalar_upper(value)
        assert mantissa.bit_length() <= SH_BITS and value <= mantissa * Fraction(2) ** exponent
    assert max((beta ** 4 / 16) ** i for i in range(3)) == 1
    assert (SH_BETA ** 4 / 16) ** 128 == (SH_BETA / 2) ** SH_WIDTH
    print(dict(exact_kernel_entries=kernel_entries, dyadic_power_entries=product_entries,
               exact_regional_comparisons=region_entries, dyadic_global_moments=moment_checks,
               upward_rounding_and_structural_zeros=True, free_endpoint_checked=True,
               full_code_certificate=False), flush=True)


def selector_hit_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, action='append')
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
            prepared = selector_hit_prepare(progress=True)
            for profile in args.profile or (SH_PROFILE,):
                print(selector_hit_profile(profile, prepared, progress=True), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    selector_hit_main()
