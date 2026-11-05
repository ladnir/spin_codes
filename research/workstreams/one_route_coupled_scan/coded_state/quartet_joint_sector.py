"""Four-sector positive comparison of the exact joint availability kernel.

Standalone: no injected helpers or local imports are required.  The four
potentials have disjoint supports, so entrywise quota minima are valid in
the common boundary basis.  All polynomial transforms and regional powers
are positive.  --checks includes independent rational physical-mask checks;
--probe is a FLOAT diagnostic, never a full-code or directed certificate.
Run numerical jobs serially.  This helper never starts a subprocess.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import comb, lgamma
import os
import shutil
import threading
import time

import numpy as np


JS_LD = np.longdouble
JS_N = 1 << 21
JS_W = 943718


def _js_ld(value):
    if isinstance(value, Fraction):
        return JS_LD(value.numerator) / JS_LD(value.denominator)
    return JS_LD(value)


def joint_sector_index(u, e):
    return 2 * int(u > 0) + int(e > 0)


def joint_sector_zero(dtype=JS_LD):
    zero = np.zeros((4, 4), dtype=dtype)
    zero[0, 0] = zero[1, 0] = 1
    return zero


def joint_sector_prepare(width, exact=False):
    if width < 4 or width % 4:
        raise ValueError('Width must be a positive multiple of four')
    dtype = object if exact else JS_LD
    cast = Fraction if exact else _js_ld
    denominator = np.ones((width + 1, width + 1), dtype=dtype)
    for u in range(width + 1):
        denominator[u, :width - u + 1] = [
            cast(comb(width, u) * comb(width - u, e))
            for e in range(width - u + 1)]
    return dict(width=width, exact=exact, dtype=dtype, cast=cast,
                denominator=denominator)


def _js_powers(local, groups, one, dtype):
    powers = [np.array([one], dtype=dtype)]
    for _ in range(groups):
        powers.append(np.convolve(powers[-1], local))
    return powers


def _js_separated_product(left, right, groups, width, dtype, first=0):
    """Sum C(groups,k) left[k](x) right[groups-k](z), positively."""
    lhs = np.zeros((width + 1, groups + 1), dtype=dtype)
    rhs = np.zeros_like(lhs)
    for k in range(first, groups + 1):
        lhs[:len(left[k]), k] = left[k] * comb(groups, k)
        rhs[:len(right[groups - k]), k] = right[groups - k]
    return lhs @ rhs.T


def joint_sector_actions(width, m, b, p, r, alpha, t=1, s=1,
                         coverage=1, exact=False, prepared=None):
    """Return K applied to all four disjoint potential columns.

    Sectors are 00, 0+, +0, ++.  Their potentials are respectively
    1, 1, t**u, t**u*s**e on their supports and zero elsewhere.
    Exact mode uses Fraction coefficients and is intended only for toys.
    """
    if prepared is None:
        prepared = joint_sector_prepare(width, exact)
    if prepared['width'] != width or prepared['exact'] != exact:
        raise ValueError('Prepared data does not match width/exact mode')
    dtype, cast = prepared['dtype'], prepared['cast']
    m, b, p, r, alpha, t, s, coverage = map(cast, (m, b, p, r, alpha, t, s, coverage))
    if min(m, b, p, r, alpha, coverage) < 0 or min(t, s) <= 0:
        raise ValueError('Require nonnegative weights and positive t,s')
    groups = width // 4
    one, zero = cast(1), cast(0)
    birth, route4 = m * p * coverage, (r * coverage) ** 4
    available = np.zeros(5, dtype=dtype)
    excess = np.zeros(5, dtype=dtype)
    excess_unit = np.zeros(5, dtype=dtype)
    both = np.zeros(5, dtype=dtype)
    for i in range(1, 5):
        prefactor = route4 * comb(4, i)
        available[i] = prefactor * (b * t) ** 4
        excess[i] = prefactor * alpha ** i * s ** 4
        excess_unit[i] = prefactor * alpha ** i
        # Select j new available positions, k of them among the i old A.
        both[i] = prefactor * sum(
            comb(i, k) * comb(4 - i, j - k) * alpha ** (i - k)
            * s ** (4 - j) * (b * t) ** j
            for j in range(1, 4)
            for k in range(max(0, j - (4 - i)), min(i, j) + 1))
    total = available + excess + both
    available_p = _js_powers(available, groups, one, dtype)
    excess_p = _js_powers(excess, groups, one, dtype)
    excess_unit_p = _js_powers(excess_unit, groups, one, dtype)
    # Presence of both outgoing labels among the old-hit quartets.
    both_p = [np.array([zero], dtype=dtype), both.copy()]
    for k in range(1, groups):
        both_p.append(np.convolve(both_p[k], total)
                      + np.convolve(available_p[k], excess + both)
                      + np.convolve(excess_p[k], available + both))
    unhit_zero = np.array([cast(comb(4, j)) for j in range(5)], dtype=dtype)
    partial_positive = sum(cast(comb(4, j)) * (birth * t) ** j for j in range(1, 4))
    # Contains at least one new birth; this is B_t-B_0 without subtraction.
    unhit_positive = unhit_zero * partial_positive
    unhit_positive[0] += (birth * t) ** 4
    unhit_total = unhit_zero + unhit_positive
    unhit_zero_p = _js_powers(unhit_zero, groups, one, dtype)
    unhit_total_p = _js_powers(unhit_total, groups, one, dtype)
    unhit_positive_p = [np.array([zero], dtype=dtype)]
    for k in range(groups):
        unhit_positive_p.append(np.convolve(unhit_positive_p[k], unhit_total)
                               + np.convolve(unhit_zero_p[k], unhit_positive))
    actions = np.zeros((4, width + 1, width + 1), dtype=dtype)
    actions[0, 0, :] = one  # K 1_(00) = 1_[old u=0], the exact Z column.
    actions[1] = _js_separated_product(excess_unit_p, unhit_zero_p,
                                       groups, width, dtype, first=1)
    actions[2] = _js_separated_product(available_p, unhit_total_p,
                                       groups, width, dtype, first=1)
    actions[2, 0, :] += unhit_positive_p[groups]
    actions[3] = (_js_separated_product(both_p, unhit_total_p,
                                        groups, width, dtype, first=1)
                  + _js_separated_product(excess_p, unhit_positive_p,
                                           groups, width, dtype, first=1))
    for j in range(1, 4):
        for u in range(width + 1):
            actions[j, u, :width - u + 1] /= prepared['denominator'][u, :width - u + 1]
            actions[j, u, width - u + 1:] = zero
    if not exact and (not np.all(np.isfinite(actions)) or np.any(actions < 0)):
        raise ArithmeticError('Nonfinite/negative positive sector transform')
    return actions


def joint_sector_comparison(width, m, b, p, r, alpha, t=1, s=1,
                            coverage=1, exact=False, prepared=None):
    actions = joint_sector_actions(width, m, b, p, r, alpha, t, s,
                                   coverage, exact, prepared)
    cast = Fraction if exact else _js_ld
    t, s = cast(t), cast(s)
    result = np.zeros((4, 4), dtype=object if exact else JS_LD)
    witnesses = [[None] * 4 for _ in range(4)]
    if exact:
        for u in range(width + 1):
            for e in range(width - u + 1):
                sector = joint_sector_index(u, e)
                potential = cast(1) if u == 0 else t ** u * (s ** e if e else 1)
                for j in range(4):
                    ratio = actions[j, u, e] / potential
                    if witnesses[sector][j] is None or ratio > result[sector, j]:
                        result[sector, j] = ratio
                        witnesses[sector][j] = (u, e)
    else:
        # Inspect every physical state, with four vectorized reductions per
        # row rather than repeated scalar powers in the 131841-state loop.
        t_powers = t ** np.arange(width + 1)
        s_powers = s ** np.arange(width + 1)
        for u in range(width + 1):
            length = width - u + 1
            potential = np.ones(length, dtype=JS_LD) if u == 0 else t_powers[u] * s_powers[:length]
            if not np.all(np.isfinite(potential)) or np.any(potential <= 0):
                raise ArithmeticError('Potential overflow/underflow')
            ratios = actions[:, u, :length] / potential[None, :]
            for start, end, sector in ((0, 1, 0 if u == 0 else 2),
                                       (1, length, 1 if u == 0 else 3)):
                if start == end:
                    continue
                locations = np.argmax(ratios[:, start:end], axis=1) + start
                for j, e in enumerate(locations):
                    ratio = ratios[j, e]
                    if witnesses[sector][j] is None or ratio > result[sector, j]:
                        result[sector, j] = ratio
                        witnesses[sector][j] = (u, int(e))
    zero = joint_sector_zero(result.dtype)
    # The first column is analytically exact, not a rounded subtraction.
    if not np.array_equal(result[:, 0], zero[:, 0]):
        raise AssertionError('The exact general rank-one zero column changed')
    if not exact and not np.all(np.isfinite(result)):
        raise ArithmeticError('Nonfinite sector row bound')
    return result, witnesses


def _js_normalize_positive(matrix, scale=0):
    norm = matrix.max()
    if norm == 0:
        return matrix.copy(), JS_LD('-inf')
    if not np.isfinite(norm) or np.any(matrix < 0):
        raise ArithmeticError('Invalid positive regional matrix')
    return matrix / norm, JS_LD(scale) + np.log2(norm)


def joint_sector_positive_compose(left, left_scale, right, right_scale, zero):
    """(Z+L)(Z+R)-Z = LR + LZ + ZR, with exact idempotent Z."""
    terms = []
    if np.isfinite(left_scale) and np.isfinite(right_scale):
        terms.append((left @ right, left_scale + right_scale))
    if np.isfinite(left_scale):
        terms.append((left @ zero, left_scale))
    if np.isfinite(right_scale):
        terms.append((zero @ right, right_scale))
    if not terms:
        return np.zeros_like(zero), JS_LD('-inf')
    base = max(scale for _, scale in terms)
    result = np.zeros_like(zero)
    for matrix, scale in terms:
        factor = np.exp2(scale - base)
        if factor == 0 and np.any(matrix > 0):
            raise ArithmeticError('Regional scaling would discard a positive term')
        result += matrix * factor
    return _js_normalize_positive(result, base)


def joint_sector_positive_region(comparison, steps=16, zero=None):
    if steps < 1:
        raise ValueError('A region must contain at least one update')
    if zero is None:
        zero = joint_sector_zero()
    positive = comparison.copy()
    # B-Z has an identically zero first column; no near cancellation occurs.
    positive[:, 0] = 0
    if not np.array_equal(comparison[:, 0], zero[:, 0]):
        raise ValueError('Comparison has the wrong exact zero column')
    power, power_scale = _js_normalize_positive(positive)
    accumulated = None
    remaining = steps
    while remaining:
        if remaining & 1:
            if accumulated is None:
                accumulated = (power.copy(), power_scale)
            else:
                accumulated = joint_sector_positive_compose(*accumulated, power, power_scale, zero)
        remaining >>= 1
        if remaining:
            power, power_scale = joint_sector_positive_compose(power, power_scale, power, power_scale, zero)
    return accumulated


def _js_logpositive(values):
    result = np.full(values.shape, JS_LD('-inf'), dtype=JS_LD)
    mask = values > 0
    result[mask] = np.log2(values[mask])
    return result


def joint_sector_occupied_log(log_positive, terminal, regions=256, minimum=129, zero=None):
    """Log of the exact positive comparison tail, using the general Z."""
    if zero is None:
        zero = joint_sector_zero()
    if not 0 <= minimum <= regions:
        raise ValueError('Invalid occupied-region threshold')
    log_zero = _js_logpositive(zero)
    paths = np.full((regions + 1, 4), JS_LD('-inf'), dtype=JS_LD)
    paths[0, 0] = 0
    for step in range(regions):
        after = np.full_like(paths, JS_LD('-inf'))
        current = paths[:step + 1]
        empty = np.logaddexp2.reduce(current[:, :, None] + log_zero[None, :, :], axis=1)
        occupied = np.logaddexp2.reduce(current[:, :, None] + log_positive[None, :, :], axis=1)
        after[:step + 1] = empty
        after[1:step + 2] = np.logaddexp2(after[1:step + 2], occupied)
        paths = after
    return np.logaddexp2.reduce((paths[minimum:] + _js_logpositive(terminal)).ravel())


def _js_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        return [(first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                for first in [1 + sum(1 << j for j in rest)]]
    raise ValueError('Exact physical checks use width four or eight')


def _js_physical_actions(width, m, b, p, r, alpha, t, s, coverage):
    """Independent rational enumeration of partitions and new physical masks."""
    F = Fraction
    result = np.full((4, width + 1, width + 1), F(0), dtype=object)
    partitions = _js_partitions(width)
    full = (1 << width) - 1
    for u in range(width + 1):
        old_a = (1 << u) - 1
        for e in range(width - u + 1):
            old_e = ((1 << e) - 1) << u
            for partition in partitions:
                hit = sum(quad for quad in partition if quad & old_a)
                if hit & old_e:
                    continue
                for new_a in range(full + 1):
                    new_v = hit | new_a
                    if any((quad & old_e) and (quad & new_v) == quad for quad in partition):
                        continue
                    new_e = hit & ~new_a
                    born = new_a & ~hit
                    v, f = new_a.bit_count(), new_e.bit_count()
                    sector = joint_sector_index(v, f)
                    potential = F(1) if v == 0 else t ** v * (s ** f if f else 1)
                    weight = ((r * coverage) ** hit.bit_count()
                              * (m * p * coverage) ** born.bit_count()
                              * b ** (new_a & hit).bit_count()
                              * alpha ** (old_a & ~new_a).bit_count())
                    result[sector, u, e] += weight * potential / len(partitions)
    return result


def joint_sector_checks():
    F = Fraction
    weights = tuple(map(F, (F(3, 4), F(5, 4), F(1, 3), F(2, 5), F(2, 3), F(1, 2), F(2), F(4, 5))))
    entries = composed = tails = inequalities = 0
    for width in (4, 8):
        truth = _js_physical_actions(width, *weights)
        exact = joint_sector_actions(width, *weights, exact=True)
        assert np.array_equal(exact, truth)
        floating = joint_sector_actions(width, *weights)
        for u in range(width + 1):
            for e in range(width - u + 1):
                for j in range(4):
                    expected = _js_ld(truth[j, u, e])
                    assert abs(floating[j, u, e] - expected) <= JS_LD('2e-15') * max(1, abs(expected))
                    entries += 1
        bound, _ = joint_sector_comparison(width, *weights, exact=True)
        approximate, _ = joint_sector_comparison(width, *weights)
        for u in range(width + 1):
            for e in range(width - u + 1):
                i = joint_sector_index(u, e)
                potential = F(1) if u == 0 else weights[5] ** u * (weights[6] ** e if e else 1)
                for j in range(4):
                    assert truth[j, u, e] <= potential * bound[i, j]
                    inequalities += 1
        z_exact = joint_sector_zero(object)
        z_float = joint_sector_zero()
        assert np.array_equal(z_exact @ z_exact, z_exact)
        for steps in (1, 2, 3, 16):
            actual = np.linalg.matrix_power(bound, steps) - z_exact
            scaled, scale = joint_sector_positive_region(approximate, steps)
            for i in range(4):
                for j in range(4):
                    expected = _js_ld(actual[i, j])
                    assert abs(scaled[i, j] * np.exp2(scale) - expected) <= JS_LD('2e-13') * max(1, abs(expected))
                    composed += 1
        q = bound - z_exact
        end = np.array([F(1), F(0), max(F(1), weights[5] ** (-width)), F(0)], dtype=object)
        paths = {(0, 0): F(1)}
        regions = 5
        for _ in range(regions):
            after = {}
            for (k, i), value in paths.items():
                for j in range(4):
                    for increment, matrix in ((0, z_exact), (1, q)):
                        key = (k + increment, j)
                        after[key] = after.get(key, F(0)) + value * matrix[i, j]
            paths = after
        q_float = np.array([[_js_ld(value) for value in row] for row in q])
        end_float = np.array([_js_ld(value) for value in end])
        for minimum in (0, 1, 3, 5):
            expected = sum(value * end[j] for (k, j), value in paths.items() if k >= minimum)
            answer = joint_sector_occupied_log(_js_logpositive(q_float), end_float, regions, minimum, z_float)
            assert abs(answer - np.log2(_js_ld(expected))) < JS_LD('2e-13')
            tails += 1
    print(dict(exact_joint_sector_entries=entries, exact_row_inequalities=inequalities,
               positive_region_entries=composed, occupied_tail_checks=tails,
               no_dense_physical_kernel=True, full_code_certificate=False), flush=True)


def joint_sector_value(h=2048, theta=1/16, logz=-.4375, logt=0, logs=0,
                       tilts=9, progress=False):
    if not 0 < h < 8192 or not 0 < theta <= 1 or logz > 0 or tilts < 1:
        raise ValueError('Invalid field/row/tilt parameters')
    theta, logz, logt, logs = map(JS_LD, (theta, logz, logt, logs))
    q = JS_LD(2) ** 32
    x = np.exp2(logz / theta)
    c = (q + 1) * x
    lam = 1 + c
    d = max(min(JS_LD(1), lam / (q + 1)), (1 + q * x * x) / lam)
    m, b, alpha, field = c ** theta, lam ** theta, d ** theta, q ** theta
    t, s = np.exp2(logt), np.exp2(logs)
    actual = JS_LD(h) / 8192
    low = actual / (field - (field - 1) * actual)
    choose = JS_LD((lgamma(8193) - lgamma(h + 1) - lgamma(8193 - h)) / float(np.log(JS_LD(2))))
    prepared = joint_sector_prepare(512)
    region_logs = None
    for index, fraction in enumerate(np.linspace(0, 1, tilts)):
        started = time.monotonic()
        odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
        p = odds / (1 + odds)
        logpi = choose + h * np.log2(p) + (8192 - h) * np.log2(1 - p)
        comparison, witnesses = joint_sector_comparison(512, m, b, p,
                    p + (1 - p) / field, alpha, t, s, prepared=prepared)
        positive, scale = joint_sector_positive_region(comparison)
        candidate = _js_logpositive(positive) + scale - logpi
        region_logs = candidate if region_logs is None else np.minimum(region_logs, candidate)
        if progress:
            print(dict(joint_sector_tilt=index + 1, logt=float(logt), logs=float(logs),
                       p=float(p), seconds=round(time.monotonic() - started, 3),
                       row_witnesses=witnesses), flush=True)
    terminal = np.array([1, 0, np.exp2(max(JS_LD(0), -512 * logt)), 0], dtype=JS_LD)
    tail = joint_sector_occupied_log(region_logs, terminal)
    bound = theta - 32 * theta - JS_W * logz + choose - 32 * theta * 128 * h + tail
    return float(bound), dict(dimension=4, checked_physical_states=131841,
                             logt=float(logt), logs=float(logs), tilts=tilts,
                             closure_d=float(d), alpha=float(alpha), certificate=False)


def joint_sector_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--logz', type=float, default=-.4375)
    parser.add_argument('--logt', type=float, nargs='+', default=(0.,))
    parser.add_argument('--logs', type=float, default=0.)
    parser.add_argument('--tilts', type=int, default=9)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            joint_sector_checks()
        else:
            for logt in args.logt:
                started = time.monotonic()
                bound, detail = joint_sector_value(args.h, args.theta, args.logz,
                                                   logt, args.logs, args.tilts, progress=True)
                print(dict(joint_sector_rows=args.h, theta=args.theta, logz=args.logz,
                           log2_bound=bound, seconds=round(time.monotonic() - started, 3),
                           **detail), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    joint_sector_main()
