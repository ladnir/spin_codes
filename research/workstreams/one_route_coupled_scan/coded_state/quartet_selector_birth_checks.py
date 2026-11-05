"""Exact small checks for integrating forced births into selector transitions.

For a selected new quartet containing c old positions, sum all birth
subsets of its 4-c new-only positions.  If c=0, require a nonempty birth
subset.  Its bivariate coefficient is therefore

    (1+w+beta*X)^4 - 1.

The subtraction removes only the constant coefficient; the polynomial
has nonnegative coefficients.  Unselected quartets contribute (1+X)^4.
The formal full endpoint contributes beta^|old| (1+w)^(n-|old|), without
a selected-quartet or volume charge.  This file checks these identities
by independent physical enumeration at widths4 and8.  It does not run
a production moment, certify a floating bound, or write any artifacts.

For fixed w>0, set b=beta/(1+w), zprime=z*(1+w), and
u=1-(1+w)^-4.  The full integrated moment equals (1+w)^n times the raw
selector-hit moment at (b,zprime,u), with its exact free endpoint b^|old|.
Do not use selector_hit_bound: its separate birth-binomial and u^-f
corrections do not belong to this integrated count.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb
import os
import shutil
import threading
import time


def _sb_partitions(width):
    if width == 4:
        return ((15,),)
    if width == 8:
        return tuple((1 + sum(1 << i for i in rest),
                      255 ^ (1 + sum(1 << i for i in rest)))
                     for rest in combinations(range(1, 8), 3))
    raise ValueError('Exact physical checks use width4 or8')


def _sb_subsets(mask):
    subset = mask
    while True:
        yield subset
        if subset == 0:
            break
        subset = (subset - 1) & mask


def _sb_convolve(first, second):
    result = [Fraction(0)] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return result


def _sb_bivariate_product(first, second):
    result = {}
    for (x, w), a in first.items():
        for (y, v), b in second.items():
            key = (x + y, w + v)
            result[key] = result.get(key, Fraction(0)) + a * b
    return result


def selector_birth_exact_kernel(width, beta, z):
    """Return K[i][j][f], the exact coefficient of w^f, for a small width."""
    if width not in (4, 8):
        raise ValueError('This exact checker intentionally supports width4/8 only')
    beta, z = Fraction(beta), Fraction(z)
    if beta <= 0 or z <= 0:
        raise ValueError('Positive rational markers required')
    groups = width // 4
    local = {(c, d): Fraction(comb(4, c) * comb(4 - c, d)) * beta ** c
             for c in range(5) for d in range(5 - c) if c or d}
    selected = {(0, 0): Fraction(1)}
    result = [[[Fraction(0)] * (width + 1) for _ in range(groups + 1)]
              for _ in range(groups + 1)]
    for j in range(groups + 1):
        remaining = width - 4 * j
        unselected = {(c, 0): Fraction(comb(remaining, c))
                      for c in range(remaining + 1)}
        polynomial = _sb_bivariate_product(selected, unselected)
        for i, d in product(range(groups + 1), range(width + 1)):
            result[i][j][d] = (comb(groups, j) * z ** (4 * j)
                               * polynomial.get((4 * i, d), Fraction(0))
                               / comb(width, 4 * i))
        selected = _sb_bivariate_product(selected, local)
    return result


def _sb_hit_kernel(width, beta, z, u):
    """Independent univariate coefficient formula for the raw hit kernel."""
    groups = width // 4
    local = [u] + [Fraction(comb(4, c)) * beta ** c for c in range(1, 5)]
    selected = [Fraction(1)]
    result = [[Fraction(0)] * (groups + 1) for _ in range(groups + 1)]
    for j in range(groups + 1):
        polynomial = _sb_convolve(selected, [Fraction(comb(width - 4 * j, c))
                                             for c in range(width - 4 * j + 1)])
        for i in range(groups + 1):
            result[i][j] = (comb(groups, j) * z ** (4 * j)
                            * polynomial[4 * i] / comb(width, 4 * i))
        selected = _sb_convolve(selected, local)
    return result


def _sb_evaluate(polynomial, marker):
    value = Fraction(0)
    for coefficient in reversed(polynomial):
        value = value * marker + coefficient
    return value


def _sb_birth_moment(kernel, beta, width, updates):
    """Exact birth-degree polynomial, including all free-endpoint births."""
    groups = width // 4
    row = [[Fraction(int(i == 0))] for i in range(groups + 1)]
    for _ in range(updates - 1):
        new = [[Fraction(0)] * (len(row[0]) + width) for _ in range(groups + 1)]
        for i, j in product(range(groups + 1), repeat=2):
            contribution = _sb_convolve(row[i], kernel[i][j])
            for degree, coefficient in enumerate(contribution):
                new[j][degree] += coefficient
        row = new
    result = [Fraction(0)] * (len(row[0]) + width)
    for i in range(groups + 1):
        endpoint = [beta ** (4 * i) * comb(width - 4 * i, d)
                    for d in range(width - 4 * i + 1)]
        for degree, coefficient in enumerate(_sb_convolve(row[i], endpoint)):
            result[degree] += coefficient
    return result


def selector_birth_checks():
    beta, z = Fraction(3), Fraction(2, 3)
    physical_coefficients = physical_birth_choices = reparameterized_entries = 0
    endpoint_checks = moment_checks = coefficient_bounds = 0
    for width in (4, 8):
        groups = width // 4
        partitions = _sb_partitions(width)
        kernel = selector_birth_exact_kernel(width, beta, z)
        for i in range(groups + 1):
            old = (1 << (4 * i)) - 1
            physical = [[Fraction(0)] * (width + 1) for _ in range(groups + 1)]
            for partition in partitions:
                for selector in range(1 << groups):
                    blocks = [block for j, block in enumerate(partition)
                              if (selector >> j) & 1]
                    new = sum(blocks)
                    for births in _sb_subsets(new & ~old):
                        if any(not (block & (old | births)) for block in blocks):
                            continue
                        physical[len(blocks)][births.bit_count()] += (
                            z ** new.bit_count() * beta ** (new & old).bit_count()
                            / len(partitions))
                        physical_birth_choices += 1
            for j, degree in product(range(groups + 1), range(width + 1)):
                assert physical[j][degree] == kernel[i][j][degree]
                physical_coefficients += 1
        for w in (Fraction(1, 4), Fraction(1, 2), Fraction(1)):
            b = beta / (1 + w)
            zprime = z * (1 + w)
            u = 1 - (1 + w) ** -4
            hit = _sb_hit_kernel(width, b, zprime, u)
            for i, j in product(range(groups + 1), repeat=2):
                assert _sb_evaluate(kernel[i][j], w) == hit[i][j]
                reparameterized_entries += 1
            for i in range(groups + 1):
                old = (1 << (4 * i)) - 1
                physical_endpoint = sum(beta ** (4 * i) * w ** births.bit_count()
                                        for births in _sb_subsets(((1 << width) - 1) & ~old))
                assert physical_endpoint == (1 + w) ** width * b ** (4 * i)
                endpoint_checks += 1
            row = [Fraction(int(i == 0)) for i in range(groups + 1)]
            for updates in (1, 2, 3):
                polynomial = _sb_birth_moment(kernel, beta, width, updates)
                exact = _sb_evaluate(polynomial, w)
                transformed = (1 + w) ** width * sum(row[i] * b ** (4 * i)
                                                    for i in range(groups + 1))
                assert exact == transformed
                # The production helper uses the safe b^width endpoint upper.
                assert transformed <= (1 + w) ** width * b ** width * sum(row)
                moment_checks += 1
                for degree, coefficient in enumerate(polynomial):
                    assert coefficient * w ** degree <= exact
                    coefficient_bounds += 1
                row = [sum(row[i] * hit[i][j] for i in range(groups + 1))
                       for j in range(groups + 1)]
    report = dict(exact_physical_birth_coefficients=physical_coefficients,
                  enumerated_valid_birth_choices=physical_birth_choices,
                  exact_reparameterized_kernel_entries=reparameterized_entries,
                  exact_free_endpoint_checks=endpoint_checks,
                  exact_moment_checks=moment_checks,
                  exact_coefficient_upper_checks=coefficient_bounds,
                  separate_birth_binomial=False, separate_hit_correction=False,
                  production_numerics_run=False, full_code_certificate=False)
    print(report, flush=True)
    return report


def selector_birth_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', required=True, action='store_true')
    parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        selector_birth_checks()
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    selector_birth_main()
