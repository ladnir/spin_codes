"""Direct availability-size transfer for the unchanged quartet construction.

For every canonical footprint, V_t=Gamma_t(A_(t-1)) union A_t and its
forced births are A_t minus Gamma_t(A_(t-1)).  Here Gamma is the union of
quartets meeting the old available set.  Thus the availability path
determines the footprint.  Summing arbitrary available paths enlarges the
canonical count without introducing additional footprint multiplicity.
Fix a whole candidate A path before exposing partitions, reconstruct its
V path, and apply the field lemma only on the canonical fixed-point event.
Only then drop that indicator.  In particular, do not condition a fresh
partition on a realized Aprev, which can depend on that partition.

The transfer charges alpha**|Aprev minus A| and the joint field weights
m,b, so it does not use an additional forced-run closure bonus.  At the
free terminal endpoint A_T=V_T: every quartet hit by A_(T-1) must lie in
A_T.  The last-step terminal vector enforces exactly that requirement.

Load after quartet_lowrank_footprint.py with its main block removed.
complete_force_region and complete_force_tail are injected, and no local
files are imported.  --checks is exact at widths 4/8; --probe reports a
floating diagnostic, never a certificate.  Run all numerical jobs serially.
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


DA_LD = np.longdouble
DA_WIDTH = 512
DA_W = 943718


def direct_available_census_constant(j):
    if not isinstance(j, int) or not 0 <= j <= 16:
        raise ValueError('Implemented census family uses integer j in 0..16')
    if j == 0:
        return 1 << 21
    if j == 3:
        return 361984
    denominator = 1 << (j + 1)
    return 512 * j + (3 * (1 << 21) + denominator - 1) // denominator


def direct_available_prepare(width):
    """Exact combinatorial probabilities converted to positive longdouble."""
    if width < 4 or width % 4:
        raise ValueError('Width must be a positive multiple of four')
    groups = width // 4
    factorial = np.ones(width + 1, dtype=DA_LD)
    for j in range(1, width + 1):
        factorial[j] = factorial[j - 1] * j
    if not np.all(np.isfinite(factorial)):
        raise ArithmeticError('Positive transforms require extended longdouble exponent range')
    choose = tuple(np.array([DA_LD(comb(n, j)) for j in range(n + 1)])
                   for n in range(width + 1))
    probabilities = np.zeros((width + 1, groups + 1), dtype=DA_LD)
    power = [1]
    for k in range(groups + 1):
        multiplier = comb(groups, k)
        for u, coefficient in enumerate(power):
            if coefficient:
                probabilities[u, k] = DA_LD(multiplier * coefficient) / DA_LD(comb(width, u))
        if k < groups:
            after = [0] * (len(power) + 4)
            for j in range(1, 5):
                coefficient = comb(4, j)
                for u, value in enumerate(power):
                    after[u + j] += coefficient * value
            power = after
    return dict(width=width, factorial=factorial, choose=choose, probabilities=probabilities)


def direct_available_matrix(width, m, b, p, r, alpha, coverage=1., prepared=None):
    """Return K and the free terminal vector minus e0 by positive arithmetic.

    Rows and columns are actual available-position counts 0..width, so
    this physical count basis is independent of routing-quota tilts.
    The two positive Taylor convolutions require a=m*p*coverage<=b.
    """
    if prepared is None:
        prepared = direct_available_prepare(width)
    elif prepared['width'] != width:
        raise ValueError('Prepared geometry has a different width')
    m, b, p, r, alpha, coverage = map(DA_LD, (m, b, p, r, alpha, coverage))
    if min(m, b, p, r, alpha, coverage) < 0:
        raise ValueError('All transfer weights must be nonnegative')
    a = m * p * coverage
    if a > b:
        raise ValueError('The positive homogenization requires m*p*coverage<=b')
    factorial = prepared['factorial']
    choose = prepared['choose']
    probabilities = prepared['probabilities']
    groups = width // 4
    indices = np.arange(width + 1, dtype=int)
    apowers = a ** indices
    bpowers = b ** indices
    alphapowers = alpha ** indices
    difference_powers = (b - a) ** indices
    route_powers = (r * coverage) ** (4 * np.arange(groups + 1, dtype=int))
    terminal_full = (r * coverage * b) ** (4 * np.arange(groups + 1, dtype=int))
    terminal_full *= (1 + a) ** (width - 4 * np.arange(groups + 1, dtype=int))
    matrix = np.empty((width + 1, width + 1), dtype=DA_LD)
    terminal = probabilities @ terminal_full
    for u in range(width + 1):
        degree = width - u
        coefficients = np.zeros(degree + 1, dtype=DA_LD)
        first_k = (u + 3) // 4
        positions = 4 * np.arange(first_k, groups + 1, dtype=int) - u
        coefficients[positions] = probabilities[u, first_k:] * route_powers[first_k:]
        # F_u(1+z): positive Taylor shift, no polynomial subtraction.
        shifted = np.convolve((coefficients * factorial[:degree + 1])[::-1],
                              1 / factorial[:degree + 1])[:degree + 1][::-1]
        shifted /= factorial[:degree + 1]
        # (1+a*y)^degree F_u(1+(b-a)y/(1+a*y)).
        left = shifted * factorial[degree::-1] * difference_powers[:degree + 1]
        right = apowers[:degree + 1] / factorial[:degree + 1]
        homogenized = np.convolve(left, right)[:degree + 1] / factorial[degree::-1]
        old = choose[u] * alphapowers[u::-1] * bpowers[:u + 1]
        matrix[u, :] = np.convolve(old, homogenized)
    terminal[0] = np.expm1(width * np.log1p(a))
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(terminal)):
        raise ArithmeticError('Nonfinite direct availability transfer; diagnostic invalid')
    return matrix, terminal


def _direct_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        partitions = []
        for other in combinations(range(1, 8), 3):
            first = 1 + sum(1 << j for j in other)
            partitions.append((first, 255 ^ first))
        return partitions
    raise ValueError('Independent exhaustive partition checks use widths 4 and 8')


def _direct_exact_convolve(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def _direct_exact_matrix(width, m, b, p, r, alpha, coverage):
    """Independent enumeration of every fixed old mask and every new mask."""
    partitions = _direct_partitions(width)
    mask_limit = 1 << width
    matrix = [[Fraction(0)] * (width + 1) for _ in range(width + 1)]
    end = [Fraction(0)] * (width + 1)
    seen = set()
    route_powers = [(r * coverage) ** j for j in range(width + 1)]
    birth_powers = [(m * p * coverage) ** j for j in range(width + 1)]
    bulk_powers = [b ** j for j in range(width + 1)]
    close_powers = [alpha ** j for j in range(width + 1)]
    weights = {}
    for old in range(mask_limit):
        u = old.bit_count()
        row = [Fraction(0)] * (width + 1)
        terminal = Fraction(0)
        census = {}
        terminal_census = {}
        for partition in partitions:
            hit = 0
            for quad in partition:
                if quad & old:
                    hit |= quad
            hit_count = hit.bit_count()
            for new in range(mask_limit):
                v = new.bit_count()
                born = (new & ~hit).bit_count()
                closed = (old & ~new).bit_count()
                key = (v, hit_count, born, closed)
                census[key] = census.get(key, 0) + 1
                if not (hit & ~new):
                    terminal_census[key] = terminal_census.get(key, 0) + 1
        for key, count in census.items():
            if key not in weights:
                v, hit_count, born, closed = key
                weights[key] = (route_powers[hit_count] * birth_powers[born]
                                * bulk_powers[v - born] * close_powers[closed] / len(partitions))
            row[key[0]] += count * weights[key]
        terminal = sum(count * weights[key] for key, count in terminal_census.items())
        if u in seen:
            assert row == matrix[u] and terminal == end[u]
        else:
            matrix[u], end[u] = row, terminal
            seen.add(u)
    assert len(seen) == width + 1
    return matrix, end


def _direct_exact_formula(width, m, b, p, r, alpha, coverage):
    """Rational coefficient formula, separate from mask enumeration and float transforms."""
    groups = width // 4
    powers = [[Fraction(1)]]
    for _ in range(groups):
        powers.append(_direct_exact_convolve(powers[-1], [0, 4, 6, 4, 1]))
    result = [[Fraction(0)] * (width + 1) for _ in range(width + 1)]
    terminal = [Fraction(0)] * (width + 1)
    a = m * p * coverage
    for u in range(width + 1):
        total_probability = Fraction(0)
        for k in range(groups + 1):
            if u >= len(powers[k]) or 4 * k < u:
                continue
            probability = Fraction(comb(groups, k), comb(width, u)) * powers[k][u]
            if not probability:
                continue
            total_probability += probability
            old = [comb(u, j) * alpha ** (u - j) * b ** j for j in range(u + 1)]
            hit = [comb(4 * k - u, j) * b ** j for j in range(4 * k - u + 1)]
            outside = [comb(width - 4 * k, j) * a ** j for j in range(width - 4 * k + 1)]
            polynomial = _direct_exact_convolve(_direct_exact_convolve(old, hit), outside)
            for v, coefficient in enumerate(polynomial):
                result[u][v] += probability * (r * coverage) ** (4 * k) * coefficient
            terminal[u] += probability * (r * coverage * b) ** (4 * k) * (1 + a) ** (width - 4 * k)
        assert total_probability == 1
    return result, terminal


def direct_available_checks():
    F = Fraction
    m, b, p, r, alpha, coverage = F(3, 4), F(5, 4), F(1, 3), F(2, 5), F(2, 3), F(4, 5)
    matrix_entries = terminal_checks = probability_checks = zero_checks = composed = 0
    for width in (4, 8):
        exact, exact_end = _direct_exact_matrix(width, m, b, p, r, alpha, coverage)
        rational_formula, formula_end = _direct_exact_formula(width, m, b, p, r, alpha, coverage)
        assert exact == rational_formula and exact_end == formula_end
        prepared = direct_available_prepare(width)
        for u in range(width + 1):
            assert abs(prepared['probabilities'][u].sum() - 1) <= DA_LD('2e-17')
            probability_checks += 1
        C, terminal_positive = direct_available_matrix(width, float(m), float(b), float(p), float(r),
                                                       float(alpha), float(coverage), prepared)
        for u in range(width + 1):
            truth = DA_LD(exact_end[u].numerator) / DA_LD(exact_end[u].denominator)
            assert abs(terminal_positive[u] + int(u == 0) - truth) <= DA_LD('2e-13') * max(1, abs(truth))
            terminal_checks += 1
            for v in range(width + 1):
                truth = DA_LD(exact[u][v].numerator) / DA_LD(exact[u][v].denominator)
                assert abs(C[u, v] - truth) <= DA_LD('2e-13') * max(1, abs(truth))
                matrix_entries += 1
        for special_alpha in (0., 1.):
            special, special_end = direct_available_matrix(width, float(m), float(b), float(p), float(r),
                                                           special_alpha, float(coverage), prepared)
            truth, truth_end = _direct_exact_formula(width, m, b, p, r, F(int(special_alpha)), coverage)
            for u in range(width + 1):
                for v in range(width + 1):
                    expected = DA_LD(truth[u][v].numerator) / DA_LD(truth[u][v].denominator)
                    assert abs(special[u, v] - expected) <= DA_LD('2e-13') * max(1, abs(expected))
                    matrix_entries += 1
        # b-a=0 and a=0 exercise the two positive-transform boundaries.
        for test_m, test_p, test_coverage in ((b, F(1), F(1)), (m, p, F(0))):
            special, special_end = direct_available_matrix(width, float(test_m), float(b), float(test_p),
                                                           float(r), float(alpha), float(test_coverage), prepared)
            truth, truth_end = _direct_exact_formula(width, test_m, b, test_p, r, alpha, test_coverage)
            for u in range(width + 1):
                for v in range(width + 1):
                    expected = DA_LD(truth[u][v].numerator) / DA_LD(truth[u][v].denominator)
                    assert abs(special[u, v] - expected) <= DA_LD('2e-13') * max(1, abs(expected))
                    matrix_entries += 1
            if test_coverage == 0:
                zero = np.zeros_like(special)
                zero[0, 0] = 1
                assert np.array_equal(special, zero) and np.all(special_end == 0)
                zero_checks += 1
        # Free endpoints compose as e0*K^(T-1)*terminal, including T=1.
        row = [F(int(u == 0)) for u in range(width + 1)]
        floating_row = np.zeros(width + 1, dtype=DA_LD)
        floating_row[0] = 1
        full_terminal = terminal_positive.copy()
        full_terminal[0] += 1
        for _ in range(4):
            truth = sum(row[u] * exact_end[u] for u in range(width + 1))
            expected = DA_LD(truth.numerator) / DA_LD(truth.denominator)
            assert abs(floating_row @ full_terminal - expected) <= DA_LD('2e-12') * max(1, abs(expected))
            row = [sum(row[u] * exact[u][v] for u in range(width + 1)) for v in range(width + 1)]
            floating_row = floating_row @ C
            composed += 1
        if 'complete_force_region' in globals():
            Q, scale, end, end_scale = complete_force_region(C, terminal_positive)
            direct = np.linalg.matrix_power(C, 16)
            direct[0, 0] -= 1
            assert np.max(abs(Q * np.exp2(scale) - direct)) <= DA_LD('2e-12') * max(1, direct.max())
            direct_end = np.linalg.matrix_power(C, 15) @ full_terminal
            direct_end[0] -= 1
            assert np.max(abs(end * np.exp2(end_scale) - direct_end)) <= DA_LD('2e-12') * max(1, direct_end.max())
    print(dict(direct_available_matrix_entries=matrix_entries, independent_all_old_and_new_masks=True,
               terminal_checks=terminal_checks, partition_probability_checks=probability_checks,
               zero_operator_checks=zero_checks, composed_free_terminal_checks=composed,
               optional_region_checks=True, full_code_certificate=False), flush=True)


def direct_available_field_weights(theta, logz):
    theta, logz = DA_LD(theta), DA_LD(logz)
    if not 0 < theta <= 1 or logz > 0:
        raise ValueError('Require 0<theta<=1 and log2(z)<=0')
    q = DA_LD(2) ** 32
    x = np.exp2(logz / theta)
    c = (q + 1) * x
    lam = 1 + c
    d = max(min(DA_LD(1), lam / (q + 1)), (1 + q * x * x) / lam)
    return c ** theta, lam ** theta, d ** theta, d


def direct_available_value(h=2048, theta=1/16, logz=-.4375, tilts=9, logcoverage=0.,
                           logrankgap=0., logavailable=0., logoutside=0., progress=False,
                           profile=False, closuregap=0., closure_census_j=0,
                           rank_census_gap=0., rank_census_j=0, logholes=0., logexclusion=0.):
    if 'complete_force_region' not in globals() or 'complete_force_tail' not in globals():
        raise RuntimeError('Inject complete_force_region and complete_force_tail before running a probe')
    if tilts < 1 or not 0 < h < 8192:
        raise ValueError('Require a positive tilt count and 0<h<8192')
    if logcoverage > 0 or logoutside > 0 or logrankgap < 0 or logavailable < 0:
        raise ValueError('Require coverage/outside log markers<=0 and rankgap/availability log markers>=0')
    if closuregap < 0 or rank_census_gap < 0 or logholes < 0 or logexclusion < 0:
        raise ValueError('Census complement markers must be nonnegative')
    m, b, alpha, d = direct_available_field_weights(theta, logz)
    closure_cap = 16 * DA_W + direct_available_census_constant(closure_census_j) + 91
    rank_census_cap = direct_available_census_constant(rank_census_j) + 44
    # Apply the canonical class-complement markers before enlarging paths.
    # Unlike survival tags, this direct positive kernel also permits alpha>1.
    alpha *= np.exp2(DA_LD(-(15 - 2 * closure_census_j) * closuregap
                          + 2 * rank_census_j * rank_census_gap))
    m *= np.exp2(DA_LD(logrankgap + logavailable + 16 * closuregap + 32 * rank_census_gap+logholes+2*logexclusion))
    b *= np.exp2(DA_LD(logrankgap + logavailable + 16 * closuregap + 32 * rank_census_gap+logholes+logexclusion))
    field = np.exp2(DA_LD(32 * (theta + closuregap + rank_census_gap)))
    prepared = direct_available_prepare(DA_WIDTH)
    rows = 8192
    actual = DA_LD(h) / rows
    low = actual / (field - (field - 1) * actual)
    choose = (lgamma(rows + 1) - lgamma(h + 1) - lgamma(rows - h + 1)) / float(np.log(DA_LD(2)))
    logs = endlogs = None
    for tilt_index, fraction in enumerate(np.linspace(0., 1., tilts)):
        started = time.monotonic()
        odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
        p = odds / (1 + odds)
        logpi = choose + h * np.log2(p) + (rows - h) * np.log2(1 - p)
        r = p + (1 - p) * np.exp2(DA_LD(logoutside - logrankgap)) / field
        C, terminal = direct_available_matrix(DA_WIDTH, m, b, p, r, alpha,
                                              np.exp2(DA_LD(logcoverage-5*logholes-2*logexclusion)), prepared)
        if progress:
            print(dict(direct_available_stage='matrix', tilt=tilt_index + 1, dimension=len(C),
                       seconds=round(time.monotonic() - started, 3)), flush=True)
        Q, scale, end, end_scale = complete_force_region(C, terminal)
        if not np.all(Q > 0) or not np.all(end > 0):
            raise ArithmeticError('Unexpected regional zero/underflow; do not use it as an upper bound')
        candidate = np.log2(Q) + scale - logpi
        candidate_end = np.log2(end) + end_scale - logpi
        logs = candidate if logs is None else np.minimum(logs, candidate)
        endlogs = candidate_end if endlogs is None else np.minimum(endlogs, candidate_end)
        if progress:
            print(dict(direct_available_stage='region', tilt=tilt_index + 1,
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    scale = logs.max()
    Q = np.exp2(logs - scale)
    norm = Q.sum(axis=1).max()
    Q /= norm
    scale += np.log2(norm)
    end_scale = endlogs.max()
    end = np.exp2(endlogs - end_scale)
    value = complete_force_tail(Q, scale, end, end_scale)
    bound = float(theta - 32 * theta - DA_W * logz + choose - 32 * theta * 128 * h + value
                  - (1009255 - 128 * h) * logoutside - 1736703 * logcoverage
                  + (65536 - 128 * h) * logrankgap
                  + (closure_cap - 32 * 128 * h) * closuregap
                  + (rank_census_cap - 32 * 128 * h) * rank_census_gap+4*(1<<21)*logholes
                  + ((1<<21)-512)*logexclusion)
    detail = dict(dimension=DA_WIDTH + 1, alpha=float(alpha), closure_d=float(d), tilts=tilts,
                  closuregap=closuregap, closure_census_j=closure_census_j,
                  rank_census_gap=rank_census_gap, rank_census_j=rank_census_j,logholes=logholes,logexclusion=logexclusion)
    if profile:
        forward = np.zeros(len(Q), dtype=DA_LD)
        forward[0] = 1
        backward = end.copy()
        backward /= backward.max()
        for _ in range(128):
            forward = forward @ Q
            forward /= forward.sum()
        for _ in range(127):
            backward = Q @ backward
            backward /= backward.max()
        marginal = forward * backward
        marginal /= marginal.sum()
        top = np.argsort(marginal)[-10:][::-1]
        detail['relaxed_midpoint_top_available'] = [(int(u), float(marginal[u])) for u in top]
        detail['relaxed_midpoint_mean_available'] = float(np.dot(marginal, np.arange(len(Q))))
    return bound, detail


def direct_available_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--tilts', type=int, default=9)
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--logz', type=float, default=-.4375)
    parser.add_argument('--rankgap', type=float, default=0.)
    parser.add_argument('--coverage', type=float, default=0.)
    parser.add_argument('--outside', type=float, default=0.)
    parser.add_argument('--closuregap', type=float, default=0.)
    parser.add_argument('--closure-census-j', type=int, default=0)
    parser.add_argument('--rank-census-gap', type=float, default=0.)
    parser.add_argument('--rank-census-j', type=int, default=0)
    parser.add_argument('--profile', action='store_true')
    parser.add_argument('--holes', type=float, default=0.)
    parser.add_argument('--exclusion', type=float, default=0.)
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
            direct_available_checks()
            return
        started = time.monotonic()
        bound, detail = direct_available_value(args.h, args.theta, args.logz, args.tilts,
                                               logrankgap=args.rankgap, logcoverage=args.coverage,
                                               logoutside=args.outside, profile=args.profile, progress=True,
                                               closuregap=args.closuregap, closure_census_j=args.closure_census_j,
                                               rank_census_gap=args.rank_census_gap, rank_census_j=args.rank_census_j,
                                               logholes=args.holes,logexclusion=args.exclusion)
        print(dict(direct_available_rows=args.h, theta=args.theta, logz=args.logz,
                   log2_diagnostic=bound, seconds=round(time.monotonic() - started, 3),
                   certificate=False, **detail), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    direct_available_main()
