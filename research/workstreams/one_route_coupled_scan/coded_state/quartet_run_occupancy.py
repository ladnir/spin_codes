"""Joint availability-run and quartet-occupancy moment, without binning.

For arbitrary availability histories A_0=empty,A_1,...,A_T, let
L=sum|A_t|, D=sum|A_(t-1) minus A_t| and
G=sum|Gamma_t(A_(t-1))|.  The free endpoint A_T is unrestricted and
Gamma_(T+1)(A_T) is not charged.  This helper computes

  sum_A E_partitions[s^L t^D z^G]
    = e0 [diag(M_u(z)) T(s,t)]^T 1,
  M_u(z)=sum_k P(u,k) z^(4k),
  T_uv=s^v [y^v](t+y)^u(1+y)^(width-u).

For canonical footprints, recovery identifies G=F-f.  The unrestricted
availability sum is an upper count after imposing canonical conditions;
the identity does not assert that every counted history is realizable.
No routing, field, occupied-region or minimum-distance certificate is
implicit.  Probes are FLOAT diagnostics; small checks are rational.
Standalone, no local imports/subprocesses.  Run all numerical jobs serially.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import comb
import os
import shutil
import threading
import time

import numpy as np


RO_LD = np.longdouble


def _ro_cast(value):
    if isinstance(value, Fraction):
        return RO_LD(value.numerator) / RO_LD(value.denominator)
    return RO_LD(value)


def _ro_exact_convolve(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def run_occupancy_probabilities(width, exact=False):
    """P(u,k): a uniform quartet partition has k groups meeting a fixed u-set."""
    if not isinstance(width, int) or width < 4 or width % 4:
        raise ValueError('Width must be a positive integer multiple of four')
    groups = width // 4
    result = np.zeros((width + 1, groups + 1), dtype=object if exact else RO_LD)
    polynomial = [1]
    for k in range(groups + 1):
        multiplier = comb(groups, k)
        for u, coefficient in enumerate(polynomial):
            if coefficient:
                numerator, denominator = multiplier * coefficient, comb(width, u)
                result[u, k] = (Fraction(numerator, denominator) if exact else
                                RO_LD(numerator) / RO_LD(denominator))
        if k < groups:
            polynomial = _ro_exact_convolve(polynomial, (0, 4, 6, 4, 1))
    if not exact and (np.any(result < 0) or not np.all(np.isfinite(result))):
        raise ArithmeticError('Invalid occupancy probabilities')
    return result


def run_occupancy_prepare(width):
    probability = run_occupancy_probabilities(width)
    log_probability = np.full(probability.shape, RO_LD('-inf'), dtype=RO_LD)
    positive = probability > 0
    log_probability[positive] = np.log2(probability[positive])
    choose = tuple(np.array([RO_LD(comb(n, j)) for j in range(n + 1)])
                   for n in range(width + 1))
    return dict(width=width, probabilities=probability,
                log_probabilities=log_probability, choose=choose)


def _ro_normalize(matrix, scale=0):
    norm = matrix.sum(axis=-1).max()
    if norm <= 0 or not np.isfinite(norm) or not np.all(np.isfinite(matrix)) or np.any(matrix < 0):
        raise ArithmeticError('Invalid positive matrix normalization')
    return matrix / norm, RO_LD(scale) + np.log2(norm)


def run_occupancy_matrix(width, log_s, log_t, log_z, prepared=None):
    """Return a positive normalized matrix and its base-two scalar scale.

    T rows are built as convolutions of two binomial distributions.  Their
    exact row masses are kept in logarithmic form, avoiding large powers
    before normalization.  Positive entries must remain representable.
    """
    log_s, log_t, log_z = map(RO_LD, (log_s, log_t, log_z))
    if not all(np.isfinite(value) for value in (log_s, log_t, log_z)):
        raise ValueError('Logarithmic markers must be finite')
    if prepared is None:
        prepared = run_occupancy_prepare(width)
    if prepared['width'] != width:
        raise ValueError('Prepared occupancy data has a different width')
    log_old = np.logaddexp2(log_t, log_s)
    log_new = np.logaddexp2(RO_LD(0), log_s)
    old_keep, old_drop = np.exp2(log_s - log_old), np.exp2(log_t - log_old)
    new_keep, new_drop = np.exp2(log_s - log_new), np.exp2(-log_new)
    if min(old_keep, old_drop, new_keep, new_drop) <= 0:
        raise ArithmeticError('Marker range underflowed a positive local transition')
    old_binomial, new_binomial = [], []
    for n in range(width + 1):
        j = np.arange(n + 1)
        old = prepared['choose'][n] * old_keep ** j * old_drop ** (n - j)
        new = prepared['choose'][n] * new_keep ** j * new_drop ** (n - j)
        if np.any(old <= 0) or np.any(new <= 0):
            raise ArithmeticError('Marker range underflowed a positive binomial coefficient')
        old_binomial.append(old)
        new_binomial.append(new)
    occupancy_logs = np.logaddexp2.reduce(
        prepared['log_probabilities'] + 4 * np.arange(width // 4 + 1) * log_z, axis=1)
    indices = np.arange(width + 1)
    row_logs = occupancy_logs + indices * log_old + (width - indices) * log_new
    common = row_logs.max()
    matrix = np.empty((width + 1, width + 1), dtype=RO_LD)
    for u in range(width + 1):
        matrix[u] = np.convolve(old_binomial[u], new_binomial[width - u]) * np.exp2(row_logs[u] - common)
    if np.any(matrix <= 0):
        raise ArithmeticError('Global normalization underflowed a positive matrix entry')
    return _ro_normalize(matrix, common)


def _ro_power_log(matrix, scale, updates):
    """log2(e0*(2^scale matrix)^updates*1), by positive binary powers."""
    if not isinstance(updates, int) or updates < 0:
        raise ValueError('The update count must be a nonnegative integer')
    if updates == 0:
        return RO_LD(0)
    row = np.zeros(len(matrix), dtype=RO_LD)
    row[0] = 1
    row_scale = RO_LD(0)
    power, power_scale = matrix, RO_LD(scale)
    remaining = updates
    while remaining:
        if remaining & 1:
            row = row @ power
            if np.any(row <= 0):
                raise ArithmeticError('Positive row product underflowed')
            norm = row.sum()
            if not np.isfinite(norm) or norm <= 0:
                raise ArithmeticError('Invalid positive row product')
            row /= norm
            row_scale += power_scale + np.log2(norm)
        remaining >>= 1
        if remaining:
            squared = power @ power
            if np.any(squared <= 0):
                raise ArithmeticError('Positive matrix square underflowed')
            power, power_scale = _ro_normalize(squared, 2 * power_scale)
    return row_scale + np.log2(row.sum())


def run_occupancy_logmoment(width=512, updates=4096, log_s=0, log_t=0, log_z=0,
                            prepared=None):
    """The complete geometry moment; no field/routing prefactor is included."""
    if not isinstance(width, int) or width < 4 or width % 4:
        raise ValueError('Width must be a positive integer multiple of four')
    if not isinstance(updates, int) or updates < 0:
        raise ValueError('Updates must be a nonnegative integer')
    if not all(np.isfinite(value) for value in (log_s, log_t, log_z)):
        raise ValueError('Logarithmic markers must be finite')
    if updates == 0:
        return RO_LD(0)
    if log_z == 0:
        # No occupancy marker: coordinates are independent two-state chains.
        # Avoid an unnecessary 513-square power in this exact special case.
        logarithms = np.array([[0, log_s], [log_t, log_s]], dtype=RO_LD)
        common = logarithms.max()
        small = np.exp2(logarithms - common)
        if np.any(small <= 0):
            raise ArithmeticError('Positive two-state transition underflowed')
        small, scale = _ro_normalize(small, common)
        return width * _ro_power_log(small, scale, updates)
    matrix, scale = run_occupancy_matrix(width, log_s, log_t, log_z, prepared)
    return _ro_power_log(matrix, scale, updates)


def run_occupancy_profile_bound(moment, markers, box):
    """Expected number of histories in an integer (L,D,G) box; never cap at 1."""
    if len(markers) != 3 or len(box) != 3 or any(lo > hi for lo, hi in box):
        raise ValueError('Require three markers and nonempty L,D,G intervals')
    minimum = sum(RO_LD(coefficient) * interval[0 if coefficient >= 0 else 1]
                  for coefficient, interval in zip(markers, box))
    return RO_LD(moment) - minimum


def _ro_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        return [(first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                for first in [1 + sum(1 << j for j in rest)]]
    raise ValueError('Exact partition checks use width four or eight')


def _ro_exact_matrix(width, s, t, z, probability):
    matrix = np.full((width + 1, width + 1), Fraction(0), dtype=object)
    for u in range(width + 1):
        occupancy = sum(probability[u, k] * z ** (4 * k) for k in range(width // 4 + 1))
        old = [comb(u, j) * t ** (u - j) for j in range(u + 1)]
        new = [comb(width - u, j) for j in range(width - u + 1)]
        for v, coefficient in enumerate(_ro_exact_convolve(old, new)):
            matrix[u, v] = occupancy * coefficient * s ** v
    return matrix


def run_occupancy_checks():
    F = Fraction
    s, t, z = F(3, 4), F(2, 3), F(5, 4)
    probabilities = entries = endpoint_checks = no_occupancy_checks = 0
    for width in (4, 8):
        partitions = _ro_partitions(width)
        exact_probability = run_occupancy_probabilities(width, exact=True)
        exact = _ro_exact_matrix(width, s, t, z, exact_probability)
        truth = np.full_like(exact, F(0))
        s_power = [s ** j for j in range(width + 1)]
        t_power = [t ** j for j in range(width + 1)]
        z_power = [z ** j for j in range(width + 1)]
        for u in range(width + 1):
            old = (1 << u) - 1
            histogram = [0] * (width // 4 + 1)
            for partition in partitions:
                hit_count = sum(4 for quartet in partition if quartet & old)
                histogram[hit_count // 4] += 1
                for new in range(1 << width):
                    truth[u, new.bit_count()] += (s_power[new.bit_count()]
                        * t_power[(old & ~new).bit_count()] * z_power[hit_count] / len(partitions))
            for k, count in enumerate(histogram):
                assert exact_probability[u, k] == F(count, len(partitions))
                probabilities += 1
            assert sum(exact_probability[u]) == 1
        assert np.array_equal(exact, truth)
        prepared = run_occupancy_prepare(width)
        floating, scale = run_occupancy_matrix(width, np.log2(_ro_cast(s)),
            np.log2(_ro_cast(t)), np.log2(_ro_cast(z)), prepared)
        for u in range(width + 1):
            for v in range(width + 1):
                expected = _ro_cast(exact[u, v])
                assert abs(floating[u, v] * np.exp2(scale) - expected) < RO_LD('3e-15') * max(1, abs(expected))
                entries += 1
        row = np.array([F(int(u == 0)) for u in range(width + 1)], dtype=object)
        for updates in range(5):
            expected = np.log2(_ro_cast(sum(row)))
            actual = _ro_power_log(floating, scale, updates)
            assert abs(actual - expected) < RO_LD('3e-13')
            if updates == 1:
                # No Gamma(A_T) factor: the first-step moment is (1+s)^width.
                assert sum(row) == (1 + s) ** width
            endpoint_checks += 1
            row = row @ exact
        no_occupancy = _ro_exact_matrix(width, s, t, F(1), exact_probability)
        row = np.array([F(int(u == 0)) for u in range(width + 1)], dtype=object)
        for updates in range(5):
            actual = run_occupancy_logmoment(width, updates, np.log2(_ro_cast(s)), np.log2(_ro_cast(t)), 0)
            assert abs(actual - np.log2(_ro_cast(sum(row)))) < RO_LD('3e-13')
            no_occupancy_checks += 1
            row = row @ no_occupancy
    print(dict(exact_occupancy_probabilities=probabilities, exact_mask_partition_entries=entries,
               free_endpoint_checks=endpoint_checks, independent_coordinate_checks=no_occupancy_checks,
               routing_or_field_factors=False, full_code_certificate=False), flush=True)


def run_occupancy_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--width', type=int, default=512)
    parser.add_argument('--updates', type=int, default=4096)
    parser.add_argument('--log-s', type=float, default=0.)
    parser.add_argument('--log-t', type=float, default=0.)
    parser.add_argument('--log-z', type=float, default=0.)
    parser.add_argument('--profile', type=int, nargs=3, default=(900000, 45000, 1130000), metavar=('L', 'D', 'G'))
    parser.add_argument('--halfwidth', type=int, default=0)
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
            run_occupancy_checks()
            return
        if args.halfwidth < 0:
            raise ValueError('Profile half-width must be nonnegative')
        started = time.monotonic()
        markers = (args.log_s, args.log_t, args.log_z)
        moment = run_occupancy_logmoment(args.width, args.updates, *markers)
        box = tuple((value - args.halfwidth, value + args.halfwidth) for value in args.profile)
        bound = run_occupancy_profile_bound(moment, markers, box)
        print(dict(width=args.width, updates=args.updates, markers=markers,
                   log2_moment=float(moment), profile_L_D_G=tuple(args.profile),
                   halfwidth=args.halfwidth, log2_history_count_bound=float(bound),
                   seconds=round(time.monotonic() - started, 3),
                   marker_optimization_performed=False, routing_or_field_factors=False,
                   certificate=False), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    run_occupancy_main()
