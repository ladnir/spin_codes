"""Two-update separation bound for CANONICAL AVAILABILITY histories.

Write G_i=Gamma_i(A_(i-1)), E_i=G_i minus A_i, and
W=s^L t^D z^G. Canonical availability requires G_i disjoint E_(i-1).
This helper retains every other such requirement. Dropping either parity
of constraints gives

  sum_A E[W product_i separation_i] <= min(Z_even,Z_odd).

Cauchy--Schwarz with the SAME positive measure W also bounds the square
by Z_even Z_odd, but is weaker here. Its geometric mean is diagnostic only.

For an even update count 2m, K is the unrestricted one-update kernel and
B is the exact two-update kernel imposing separation at its second step:

  Z_even=e0 B^m 1,   Z_odd=e0 K B^(m-1) K 1.

Both kernels have width+1 states. The intermediate excess set is generated
by the first update and tested by the second before it is forgotten.
This is an availability relaxation: terminal E_T=empty, routing, forced
birth counts, and occupied-region requirements are not imposed. In
particular, these separation conditions must not be applied to actual-state
support profiles without a separate argument.

Factorization B=C T, with a the number of second-hit coordinates outside
the intermediate availability set:

  J_ua=sum_k P(u,k) z^(4k) choose(n-4k,a),
  C_uv=s^v sum_a P(v,(v+a)/4) z^(v+a) J_ua / choose(n-v,a)
                    * [y^v](t+y)^u(1+y)^(n-a-u).

Only a with v+a divisible by four and u+a<=n contribute. Vandermonde
eliminates the intermediate excess count. The coefficient arrays for
N=n-a obey H_N[0,v]=choose(N,v) and
H_N[u,v]=t H_(N-1)[u-1,v]+H_(N-1)[u-1,v-1], so they can be streamed
in O(n^3) positive operations and O(n^2) storage.

FLOAT probes are diagnostics, not interval certificates. The --checks
mode compares rational formulas to independent width-four/eight mask and
partition enumeration. Standalone; no local imports or subprocesses.
All numerical jobs must be run serially.
"""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb
import os
import shutil
import threading
import time

import numpy as np


TUS_LD = np.longdouble


def _tus_cast(value):
    if isinstance(value, Fraction):
        return TUS_LD(value.numerator) / TUS_LD(value.denominator)
    return TUS_LD(value)


def _tus_convolve(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, value in enumerate(left):
        for j, other in enumerate(right):
            result[i + j] += value * other
    return result


def _tus_probabilities(width, exact=False):
    if not isinstance(width, int) or width < 4 or width % 4:
        raise ValueError('Width must be a positive integer multiple of four')
    groups = width // 4
    result = np.zeros((width + 1, groups + 1), dtype=object if exact else TUS_LD)
    polynomial = [1]
    for k in range(groups + 1):
        for u, coefficient in enumerate(polynomial):
            if coefficient:
                numerator, denominator = comb(groups, k) * coefficient, comb(width, u)
                result[u, k] = (Fraction(numerator, denominator) if exact else
                                TUS_LD(numerator) / TUS_LD(denominator))
        if k < groups:
            polynomial = _tus_convolve(polynomial, (0, 4, 6, 4, 1))
    return result


def two_update_prepare(width=512):
    """Marker-independent occupancy and uniform-subset avoidance data."""
    probability = _tus_probabilities(width)
    choose = np.zeros((width + 1, width + 1), dtype=TUS_LD)
    log_choose = np.full_like(choose, TUS_LD('-inf'))
    for n in range(width + 1):
        for j in range(n + 1):
            choose[n, j] = TUS_LD(comb(n, j))
            log_choose[n, j] = np.log2(choose[n, j])
    avoidance = np.zeros((width // 4 + 1, width + 1), dtype=TUS_LD)
    for k in range(width // 4 + 1):
        remaining = width - 4 * k
        avoidance[k, :remaining + 1] = choose[remaining, :remaining + 1] / choose[width, :remaining + 1]
        if np.any(avoidance[k, :remaining + 1] <= 0):
            raise ArithmeticError('Positive avoidance probability underflowed')
    return dict(width=width, probability=probability, choose=choose,
                log_choose=log_choose, avoidance=avoidance)


def _tus_normalize(matrix, scale=0):
    if np.any(matrix < 0) or not np.all(np.isfinite(matrix)):
        raise ArithmeticError('Invalid nonnegative matrix')
    norm = matrix.sum(axis=-1).max()
    if not np.isfinite(norm) or norm <= 0:
        raise ArithmeticError('Invalid positive normalization')
    result = matrix / norm
    if np.any((matrix > 0) & (result == 0)):
        raise ArithmeticError('Matrix normalization underflowed a positive entry')
    return result, TUS_LD(scale) + np.log2(norm)


def _tus_scaled_block(base, row_logs, column_logs):
    """Positive diagonal scalings, with exact structural zeros retained."""
    rows, columns = np.isfinite(row_logs), np.isfinite(column_logs)
    if not np.any(rows) or not np.any(columns):
        return None
    row_max, column_max = row_logs[rows].max(), column_logs[columns].max()
    row_factor, column_factor = np.exp2(row_logs - row_max), np.exp2(column_logs - column_max)
    result = base * row_factor[:, None] * column_factor[None, :]
    positive = rows[:, None] & columns[None, :] & (base > 0)
    if np.any(positive & (result <= 0)) or not np.all(np.isfinite(result)):
        raise ArithmeticError('Positive block entry underflowed or overflowed')
    return result, row_max + column_max


def _tus_add_block(matrix, scale, block, block_scale):
    common = max(scale, block_scale)
    old_multiplier = np.exp2(scale - common) if np.isfinite(scale) else TUS_LD(0)
    new_multiplier = np.exp2(block_scale - common)
    old = matrix * old_multiplier
    contribution = block * new_multiplier
    if np.any((matrix > 0) & (old <= 0)) or np.any((block > 0) & (contribution <= 0)):
        raise ArithmeticError('Positive summand underflowed during common scaling')
    matrix = old
    matrix[:len(block), :len(block)] += contribution
    if not np.all(np.isfinite(matrix)):
        raise ArithmeticError('Positive block sum overflowed')
    return matrix, common


def two_update_matrices(width=512, log_s=0, log_t=0, log_z=0, prepared=None):
    """Return normalized C, T, B=C*T and unrestricted K, with log2 scales.

    Array recurrences are positive. The coefficient rows are normalized
    binomial convolutions; their masses are restored through logarithmic
    diagonal factors. Finite-range failures raise instead of discarding
    terms. Rounded floating arithmetic is not an exact upper certificate.
    """
    log_s, log_t, log_z = map(TUS_LD, (log_s, log_t, log_z))
    if not all(np.isfinite(value) for value in (log_s, log_t, log_z)):
        raise ValueError('Logarithmic markers must be finite')
    if prepared is None:
        prepared = two_update_prepare(width)
    if prepared['width'] != width:
        raise ValueError('Prepared data has a different width')
    probability = prepared['probability']
    log_choose = prepared['log_choose']
    occupancy_scale = max(TUS_LD(0), width * log_z)
    weighted = probability * np.exp2(4 * np.arange(width // 4 + 1) * log_z - occupancy_scale)
    if np.any((probability > 0) & (weighted <= 0)) or not np.all(np.isfinite(weighted)):
        raise ArithmeticError('Positive weighted occupancy probability underflowed')
    # Q_ua = 2^-occupancy_scale J_ua / choose(width,a).
    q_values = weighted @ prepared['avoidance']
    u_index = np.arange(width + 1)
    a_index = np.arange(width + 1)
    expected_q = a_index[None, :] <= width - 4 * ((u_index[:, None] + 3) // 4)
    if np.any(expected_q & (q_values <= 0)) or not np.all(np.isfinite(q_values)):
        raise ArithmeticError('Positive averaged avoidance factor underflowed')
    log_q = np.full_like(q_values, TUS_LD('-inf'))
    log_q[q_values > 0] = np.log2(q_values[q_values > 0])
    log_second = np.full((width + 1, width + 1), TUS_LD('-inf'), dtype=TUS_LD)
    for v in range(width + 1):
        for k in range(width // 4 + 1):
            a = 4 * k - v
            if 0 <= a <= width and weighted[v, k] > 0:
                log_second[v, a] = np.log2(weighted[v, k])
    log_old_mass = np.logaddexp2(TUS_LD(0), log_t)
    old_drop, old_keep = np.exp2(log_t - log_old_mass), np.exp2(-log_old_mass)
    if min(old_drop, old_keep) <= 0:
        raise ArithmeticError('A positive local coefficient underflowed')
    c_matrix = np.zeros((width + 1, width + 1), dtype=TUS_LD)
    c_scale = TUS_LD('-inf')
    h_previous = None
    for reduced_width in range(width + 1):
        # H rows divided by (1+t)^u 2^(reduced_width-u).
        h_current = np.zeros((reduced_width + 1, reduced_width + 1), dtype=TUS_LD)
        h_current[0] = prepared['choose'][reduced_width, :reduced_width + 1] * np.exp2(TUS_LD(-reduced_width))
        if reduced_width:
            h_current[1:, :-1] = old_drop * h_previous
            h_current[1:, 1:] += old_keep * h_previous
        if np.any(h_current <= 0) or not np.all(np.isfinite(h_current)):
            raise ArithmeticError('Positive transition coefficient underflowed')
        a = width - reduced_width
        indices = np.arange(reduced_width + 1)
        row_logs = log_q[:reduced_width + 1, a] + indices * log_old_mass + reduced_width - indices
        column_logs = (log_second[:reduced_width + 1, a] + indices * log_s
                       + log_choose[width, a] - log_choose[width - indices, a])
        scaled = _tus_scaled_block(h_current, row_logs, column_logs)
        if scaled is not None:
            block, block_scale = scaled
            c_matrix, c_scale = _tus_add_block(c_matrix, c_scale, block, block_scale + 2 * occupancy_scale)
        h_previous = h_current
    c_matrix, c_scale = _tus_normalize(c_matrix, c_scale)
    indices = np.arange(width + 1)
    row_logs = indices * log_old_mass + width - indices
    t_matrix, t_scale = _tus_scaled_block(h_current, row_logs, indices * log_s)
    t_matrix, t_scale = _tus_normalize(t_matrix, t_scale)
    b_product = c_matrix @ t_matrix
    if np.any(b_product <= 0):
        raise ArithmeticError('Positive two-update matrix product underflowed')
    b_matrix, b_scale = _tus_normalize(b_product, c_scale + t_scale)
    occupancy = weighted.sum(axis=1)
    k_product = occupancy[:, None] * t_matrix
    if np.any(k_product <= 0):
        raise ArithmeticError('Positive unrestricted transition underflowed')
    k_matrix, k_scale = _tus_normalize(k_product, occupancy_scale + t_scale)
    return dict(C=c_matrix, C_scale=c_scale, T=t_matrix, T_scale=t_scale,
                B=b_matrix, B_scale=b_scale, K=k_matrix, K_scale=k_scale)


def _tus_row_normalize(row, scale=0):
    if np.any(row < 0) or not np.all(np.isfinite(row)):
        raise ArithmeticError('Invalid positive row')
    norm = row.sum()
    if norm <= 0 or not np.isfinite(norm):
        raise ArithmeticError('Invalid row normalization')
    normalized = row / norm
    if np.any((row > 0) & (normalized <= 0)):
        raise ArithmeticError('Positive row normalization underflowed')
    return normalized, TUS_LD(scale) + np.log2(norm)


def _tus_apply_powers(row, row_scale, matrix, matrix_scale, exponent):
    if not isinstance(exponent, int) or exponent < 0:
        raise ValueError('Matrix exponent must be a nonnegative integer')
    while exponent:
        if exponent & 1:
            product = row @ matrix
            if np.any(product <= 0):
                raise ArithmeticError('Positive row product underflowed')
            row, row_scale = _tus_row_normalize(product, row_scale + matrix_scale)
        exponent >>= 1
        if exponent:
            product = matrix @ matrix
            if np.any(product <= 0):
                raise ArithmeticError('Positive matrix square underflowed')
            matrix, matrix_scale = _tus_normalize(product, 2 * matrix_scale)
    return row, row_scale


def two_update_logmoment(width=512, updates=4096, log_s=0, log_t=0, log_z=0,
                         prepared=None, matrices=None):
    """Best parity upper moment; the two parity values share binary powers."""
    if not isinstance(updates, int) or updates < 0 or updates % 2:
        raise ValueError('This helper requires a nonnegative even update count')
    if updates == 0:
        return dict(log2_bound=TUS_LD(0), log2_even=TUS_LD(0), log2_odd=TUS_LD(0), log2_cauchy=TUS_LD(0))
    if matrices is None:
        matrices = two_update_matrices(width, log_s, log_t, log_z, prepared)
    row = np.zeros(width + 1, dtype=TUS_LD)
    row[0] = 1
    odd_row, odd_scale = _tus_row_normalize(row @ matrices['K'], matrices['K_scale'])
    rows = np.stack((row, odd_row))
    scales = np.array((TUS_LD(0), odd_scale), dtype=TUS_LD)
    power, power_scale = matrices['B'], matrices['B_scale']
    exponent = updates // 2 - 1
    while exponent:
        if exponent & 1:
            products = rows @ power
            if np.any(products <= 0) or not np.all(np.isfinite(products)):
                raise ArithmeticError('Positive parity row product underflowed')
            norms = products.sum(axis=1)
            rows = products / norms[:, None]
            if np.any(rows <= 0) or not np.all(np.isfinite(rows)):
                raise ArithmeticError('Positive parity row normalization failed')
            scales += power_scale + np.log2(norms)
        exponent >>= 1
        if exponent:
            product = power @ power
            if np.any(product <= 0):
                raise ArithmeticError('Positive shared matrix square underflowed')
            power, power_scale = _tus_normalize(product, 2 * power_scale)
    even_product, odd_product = rows[0] @ matrices['B'], rows[1] @ matrices['K']
    if np.any(even_product <= 0) or np.any(odd_product <= 0):
        raise ArithmeticError('Positive final parity product underflowed')
    even_row, even_scale = _tus_row_normalize(even_product, scales[0] + matrices['B_scale'])
    odd_row, odd_scale = _tus_row_normalize(odd_product, scales[1] + matrices['K_scale'])
    even = even_scale + np.log2(even_row.sum())
    odd = odd_scale + np.log2(odd_row.sum())
    return dict(log2_bound=min(even, odd), log2_even=even, log2_odd=odd, log2_cauchy=(even + odd) / 2)


def _tus_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        return [(first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                for first in [1 + sum(1 << bit for bit in rest)]]
    raise ValueError('Exact enumeration uses widths four and eight')


def _tus_gamma(partition, mask):
    result = 0
    for quartet in partition:
        if quartet & mask:
            result |= quartet
    return result


def _tus_exact_matrices(width, s, t, z):
    probability = _tus_probabilities(width, exact=True)
    c_matrix = np.full((width + 1, width + 1), Fraction(0), dtype=object)
    t_matrix = np.full_like(c_matrix, Fraction(0))
    occupancy = [sum(probability[u, k] * z ** (4 * k) for k in range(width // 4 + 1))
                 for u in range(width + 1)]
    for u in range(width + 1):
        coefficients = _tus_convolve([comb(u, j) * t ** (u - j) for j in range(u + 1)],
                                    [comb(width - u, j) for j in range(width - u + 1)])
        for v in range(width + 1):
            t_matrix[u, v] = s ** v * coefficients[v]
        for a in range(width - u + 1):
            j_value = sum(probability[u, k] * z ** (4 * k) * comb(width - 4 * k, a)
                          for k in range(width // 4 + 1) if width - 4 * k >= a)
            coefficients = _tus_convolve([comb(u, j) * t ** (u - j) for j in range(u + 1)],
                                        [comb(width - a - u, j) for j in range(width - a - u + 1)])
            for v, coefficient in enumerate(coefficients):
                if (v + a) % 4 == 0:
                    c_matrix[u, v] += (s ** v * probability[v, (v + a) // 4]
                        * z ** (v + a) * j_value * coefficient / comb(width - v, a))
    k_matrix = np.array([[occupancy[u] * t_matrix[u, v] for v in range(width + 1)]
                         for u in range(width + 1)], dtype=object)
    return dict(C=c_matrix, T=t_matrix, B=c_matrix @ t_matrix, K=k_matrix)


def _tus_enumerated_c(width, s, t, z):
    """Independent labeled-mask/partition census; no occupancy formulas."""
    partitions = _tus_partitions(width)
    gamma = [[_tus_gamma(partition, mask) for mask in range(1 << width)]
             for partition in partitions]
    second_histograms = {}
    truth = np.full((width + 1, width + 1), Fraction(0), dtype=object)
    histogram_entries = 0
    for u in range(width + 1):
        old = (1 << u) - 1
        histogram = Counter()
        for first in range(len(partitions)):
            hit = gamma[first][old]
            g = hit.bit_count()
            for new in range(1 << width):
                excess = hit & ~new
                key = (new, excess)
                if key not in second_histograms:
                    counts = Counter(gamma[second][new].bit_count() for second in range(len(partitions))
                                     if not (gamma[second][new] & excess))
                    second_histograms[key] = counts
                v, drops = new.bit_count(), (old & ~new).bit_count()
                for h, count in second_histograms[key].items():
                    histogram[v, drops, g, h] += count
        for (v, drops, g, h), count in histogram.items():
            truth[u, v] += Fraction(count, len(partitions) ** 2) * s ** v * t ** drops * z ** (g + h)
        histogram_entries += len(histogram)
    return truth, histogram_entries


def _tus_exact_path_moment(width, updates, s, t, z, constrained_parities):
    """Tiny full mask/excess DP, retaining the actual previous excess set."""
    partitions = _tus_partitions(width)
    distribution = {(0, 0): Fraction(1)}
    for step in range(1, updates + 1):
        following = {}
        for (old, old_excess), weight in distribution.items():
            for partition in partitions:
                hit = _tus_gamma(partition, old)
                if step % 2 in constrained_parities and hit & old_excess:
                    continue
                for new in range(1 << width):
                    key = (new, hit & ~new)
                    value = weight * s ** new.bit_count() * t ** (old & ~new).bit_count() * z ** hit.bit_count() / len(partitions)
                    following[key] = following.get(key, Fraction(0)) + value
        distribution = following
    return sum(distribution.values())


def two_update_checks():
    s, t, z = Fraction(3, 4), Fraction(2, 3), Fraction(5, 4)
    entries = histograms = floating_entries = parity_checks = 0
    for width in (4, 8):
        exact = _tus_exact_matrices(width, s, t, z)
        truth, histogram_count = _tus_enumerated_c(width, s, t, z)
        unrestricted_two = exact['K'] @ exact['K']
        assert np.array_equal(exact['C'], truth)
        assert np.array_equal(exact['B'][0], unrestricted_two[0])
        assert all(exact['B'][u, v] <= unrestricted_two[u, v]
                   for u in range(width + 1) for v in range(width + 1))
        entries += (width + 1) ** 2
        histograms += histogram_count
        floating = two_update_matrices(width, np.log2(_tus_cast(s)), np.log2(_tus_cast(t)), np.log2(_tus_cast(z)))
        for name in ('C', 'T', 'B', 'K'):
            for u in range(width + 1):
                for v in range(width + 1):
                    expected = _tus_cast(exact[name][u, v])
                    actual = floating[name][u, v] * np.exp2(floating[name + '_scale'])
                    assert abs(actual - expected) <= TUS_LD('5e-14') * max(TUS_LD(1), abs(expected))
                    floating_entries += 1
        if width == 4:
            for updates in (2, 4):
                even = _tus_exact_path_moment(width, updates, s, t, z, {0})
                odd = _tus_exact_path_moment(width, updates, s, t, z, {1})
                all_constraints = _tus_exact_path_moment(width, updates, s, t, z, {0, 1})
                b_power = np.linalg.matrix_power(exact['B'], updates // 2)
                odd_power = exact['K'] @ np.linalg.matrix_power(exact['B'], updates // 2 - 1) @ exact['K']
                assert even == sum(b_power[0])
                assert odd == sum(odd_power[0])
                assert all_constraints ** 2 <= even * odd
                assert all_constraints <= min(even, odd)
                result = two_update_logmoment(width, updates, matrices=floating)
                assert abs(result['log2_even'] - np.log2(_tus_cast(even))) < TUS_LD('5e-13')
                assert abs(result['log2_odd'] - np.log2(_tus_cast(odd))) < TUS_LD('5e-13')
                parity_checks += 1
    print(dict(exact_mask_partition_C_entries=entries, enumeration_histogram_terms=histograms,
               floating_matrix_entries=floating_entries, exact_parity_endpoint_checks=parity_checks,
               terminal_excess_enforced=False, routing_or_field_factors=False,
               full_code_certificate=False), flush=True)


def two_update_main():
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
            two_update_checks()
            return
        started = time.monotonic()
        markers = (args.log_s, args.log_t, args.log_z)
        result = two_update_logmoment(args.width, args.updates, *markers)
        count_bound = result['log2_bound'] - sum(TUS_LD(marker) * value for marker, value in zip(markers, args.profile))
        print(dict(width=args.width, updates=args.updates, markers=markers,
                   **{key: float(value) for key, value in result.items()},
                   profile_L_D_G=tuple(args.profile), log2_history_count_bound=float(count_bound),
                   seconds=round(time.monotonic() - started, 3),
                   terminal_excess_enforced=False, occupied_region_constraint=False,
                   routing_or_field_factors=False, certificate=False), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    two_update_main()
