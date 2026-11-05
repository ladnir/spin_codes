"""Positive global availability-tag transfer for the unchanged quartet code.

Load after quartet_lowrank_footprint.py with its main block removed.  The
regional helpers are injected through globals; this file imports no local
helper.  --checks uses exact rational partition enumeration at widths 4/8.
--probe is a FLOAT upper-bound diagnostic, never a certificate.

For a canonical footprint, independent tags on A_(t-1), with probability
rho=1-alpha, survive in A_t with probability alpha**D_A.  Floor rounding
retains floor_J(number of tags) uniformly selected tags.  Barycentric
rounding chooses adjacent tag levels around that number, then samples a
uniform subset of the entire available set; it need not retain original
Bernoulli tags.  Both enlarge the exact survival probability.  Barycentric
rounding improves the floor bound, by convexity of C(ell-g,D)/C(ell,D).
Only after this canonical inequality do we add the positive fully-forced
quartet relaxation.  Extra annotations need not describe geometric A_t.

The all-run closure factor is d=max(e,(1+q*x*x)/lambda), where
e=min(1,lambda/(q+1)).  It pays every closed availability run once.
Birth m=c**theta is UNMODIFIED: do not also apply the old forced-run bonus
or its terminal correction.  The terminal state remains free.
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


AT_LD = np.longdouble
AT_WIDTH = 512
AT_W = 943718


def available_census_constant(j):
    if not isinstance(j,int) or not 0<=j<=16:
        raise ValueError('Implemented census family uses integer j in 0..16')
    if j==0:return 1<<21
    if j==3:return 361984
    denominator=1<<(j+1)
    return 512*j+(3*(1<<21)+denominator-1)//denominator


def available_tag_levels(width, levels):
    levels = tuple(int(value) for value in levels)
    if width < 4 or width % 4 or not levels or levels[0] != 0:
        raise ValueError('Width must be a positive quartet width; levels must start at zero')
    if tuple(sorted(set(levels))) != levels or levels[-1] > width:
        raise ValueError('Tag levels must be strictly increasing and at most the width')
    return levels


def available_tag_states(width, levels):
    levels = available_tag_levels(width, levels)
    return tuple((i, g) for i in range(width // 4 + 1) for g in levels if g <= 4 * i)


def available_tag_probabilities(width, alpha, levels, rounding='floor'):
    """Positive binomial sums; no differences of nearly equal CDF values."""
    levels = available_tag_levels(width, levels)
    alpha = AT_LD(alpha)
    if not 0 <= alpha <= 1:
        raise ValueError('alpha must be a probability')
    if rounding not in ('floor', 'barycentric'):
        raise ValueError('Tag rounding must be floor or barycentric')
    result = np.zeros((len(levels), width + 1), dtype=AT_LD)
    distribution = np.ones(1, dtype=AT_LD)
    local = np.array([alpha, 1 - alpha], dtype=AT_LD)
    intervals = []
    for j, lo in enumerate(levels):
        hi = levels[j + 1] if j + 1 < len(levels) else width + 1
        if rounding == 'barycentric' and j + 1 < len(levels):
            fraction = np.arange(hi - lo, dtype=AT_LD) / (hi - lo)
            intervals.append((lo, hi, 1 - fraction, fraction))
        else:
            intervals.append((lo, hi, None, None))
    for ell in range(width + 1):
        for j, (lo, hi, down, up) in enumerate(intervals):
            if down is not None and hi <= ell:
                # An upper level larger than ell cannot be sampled: use floor there.
                result[j, ell] += np.dot(distribution[lo:hi], down)
                result[j + 1, ell] += np.dot(distribution[lo:hi], up)
            else:
                result[j, ell] += distribution[lo:hi].sum()
        if ell < width:
            distribution = np.convolve(distribution, local)
    return result


def available_tag_prepare(width, alpha, levels, rounding='floor'):
    """Geometry-only data, shared across every routing-quota tilt."""
    levels = available_tag_levels(width, levels)
    states = available_tag_states(width, levels)
    factorial = np.ones(width + 1, dtype=AT_LD)
    for j in range(1, width + 1):
        factorial[j] = factorial[j - 1] * j
    if not np.all(np.isfinite(factorial)):
        raise ArithmeticError('This positive Taylor transform requires extended longdouble range')
    powers = [np.ones(1, dtype=AT_LD)]
    local = np.array([0, 4, 6, 4, 1], dtype=AT_LD)
    for _ in range(width // 4):
        powers.append(np.convolve(powers[-1], local))
    binomial = tuple(np.array([AT_LD(comb(n, j)) for j in range(n + 1)])
                     for n in range(width + 1))
    by_level = {}
    for g in levels:
        rows = np.array([j for j, (_, old_g) in enumerate(states) if old_g == g], dtype=int)
        exponents = np.array([4 * states[j][0] - g for j in rows], dtype=int)
        derivative = np.array([AT_LD(comb(j + g, g)) for j in range(width - g + 1)])
        by_level[g] = rows, exponents, derivative
    return dict(width=width, levels=levels, alpha=AT_LD(alpha), rounding=rounding,
                states=states, factorial=factorial,
                binomial=binomial, powers=powers, by_level=by_level,
                probabilities=available_tag_probabilities(width, alpha, levels, rounding))


def available_tag_matrix(width, m, b, p, r, alpha, levels=(0, 32), coverage=1., prepared=None,
                         rounding='floor'):
    """Return reduced C=VU and the terminal vector minus e_0, all positively.

    The physical U depends on m,b,alpha and levels, but not p,r,coverage.
    Thus regional minima over p use one common physical factor and basis.
    """
    if prepared is None:
        prepared = available_tag_prepare(width, alpha, levels, rounding)
    elif (prepared['width'] != width or prepared['levels'] != tuple(levels)
          or prepared['alpha'] != AT_LD(alpha) or prepared['rounding'] != rounding):
        raise ValueError('Prepared geometry does not match width/levels/alpha/rounding')
    m, b, p, r, coverage = map(AT_LD, (m, b, p, r, coverage))
    if min(m, b, p, r, coverage) < 0:
        raise ValueError('The positive transfer requires nonnegative weights')
    groups = width // 4
    states = prepared['states']
    factorial = prepared['factorial']
    binomial = prepared['binomial']
    probabilities = prepared['probabilities']
    a = m * p * coverage
    apowers = a ** np.arange(width + 1, dtype=int)
    bpowers = b ** np.arange(width + 1, dtype=int)
    prefactor = np.array([
        AT_LD(comb(groups, i)) * (r * coverage) ** (4 * i)
        / (AT_LD(comb(width, 4 * i)) * AT_LD(comb(4 * i, g)))
        for i, g in states], dtype=AT_LD)
    result = np.zeros((len(states), len(states)), dtype=AT_LD)
    level_index = {g: j for j, g in enumerate(prepared['levels'])}
    for column, (k, h) in enumerate(states):
        degree = 4 * k
        coefficients = prepared['powers'][k] * probabilities[level_index[h], :degree + 1]
        # [z^j] P(a+b*z), via a positive reversed factorial convolution.
        shifted = np.convolve((coefficients * factorial[:degree + 1])[::-1],
                              apowers[:degree + 1] / factorial[:degree + 1])[:degree + 1][::-1]
        shifted *= bpowers[:degree + 1] / factorial[:degree + 1]
        outside = binomial[width - degree]
        for g, (rows, exponents, derivative) in prepared['by_level'].items():
            if g > degree:
                continue
            poly = np.convolve(shifted[g:] * derivative[:degree - g + 1], outside)
            result[rows, column] = prefactor[rows] * poly[exponents]
    terminal = np.array([
        AT_LD(comb(groups, i)) * (r * coverage * b) ** (4 * i)
        * (1 + a) ** (width - 4 * i) for i, _ in states], dtype=AT_LD)
    terminal[0] = np.expm1(width * np.log1p(a))
    if not np.all(np.isfinite(result)) or not np.all(np.isfinite(terminal)):
        raise ArithmeticError('Nonfinite tagged transfer; this is not an upper-bound certificate')
    return result, terminal


def _available_exact_weights(ell, alpha, levels, rounding='floor'):
    result = {g: Fraction(0) for g in levels}
    rho = 1 - alpha
    if rounding not in ('floor', 'barycentric'):
        raise ValueError('Unknown rounding rule')
    for count in range(ell + 1):
        kept = max(g for g in levels if g <= count)
        probability = comb(ell, count) * rho ** count * alpha ** (ell - count)
        upper = next((g for g in levels if g > count), None)
        if rounding == 'barycentric' and upper is not None and upper <= ell:
            result[kept] += probability * Fraction(upper - count, upper - kept)
            result[upper] += probability * Fraction(count - kept, upper - kept)
        else:
            result[kept] += probability
    assert sum(result.values()) == 1
    return result


def _available_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        result = []
        for rest in combinations(range(1, 8), 3):
            first = 1 + sum(1 << j for j in rest)
            result.append((first, 255 ^ first))
        return result
    raise ValueError('Independent partition enumeration is restricted to widths 4 and 8')


def _available_exact_factorization(width, m, b, p, r, alpha, levels, coverage, rounding='floor'):
    """Independent rational U: enumerate physical quartet partitions and unions."""
    latent = available_tag_states(width, levels)
    physical = [(f, g, 4 * i - g) for i, g in latent for f in range(width - 4 * i + 1)]
    parts = _available_partitions(width)
    groups = width // 4
    tag_weights = [_available_exact_weights(ell, alpha, levels, rounding) for ell in range(width + 1)]
    U = [[Fraction(0)] * len(latent) for _ in physical]
    V = [[Fraction(0)] * len(physical) for _ in latent]
    for si, (f, g, o) in enumerate(physical):
        forced = (1 << f) - 1
        tags = ((1 << g) - 1) << f
        covered = (1 << (f + g + o)) - 1
        bulk = covered ^ forced
        for li, (k, h) in enumerate(latent):
            for partition in parts:
                for selected in combinations(partition, k):
                    union = 0
                    for quad in selected:
                        union |= quad
                    if (forced | tags) & ~union or any(not (covered & quad) for quad in selected):
                        continue
                    ell = (covered & union).bit_count()
                    U[si][li] += (m ** f * b ** (bulk & union).bit_count() * tag_weights[ell][h]
                                  / (len(parts) * comb(groups, k)))
            if g == h and g + o == 4 * k:
                V[li][si] = comb(groups, k) * comb(width - 4 * k, f) * (p * coverage) ** f * (r * coverage) ** (4 * k)
    C = [[sum(V[i][s] * U[s][j] for s in range(len(physical))) for j in range(len(latent))]
         for i in range(len(latent))]
    terminal = [m ** f * b ** (g + o) for f, g, o in physical]
    end = [sum(V[i][s] * terminal[s] for s in range(len(physical))) for i in range(len(latent))]
    return latent, physical, U, V, C, terminal, end


def _available_exact_convolve(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def _available_exact_old129(width, m, b, p, r, coverage):
    """The old untagged coefficient formula, computed entirely rationally."""
    groups = width // 4
    a = m * p * coverage
    local = [(1 + a) ** 4 - 1] + [comb(4, j) * (1 + a) ** (4 - j) * b ** j for j in range(1, 5)]
    power = [Fraction(1)]
    C = [[Fraction(0)] * (groups + 1) for _ in range(groups + 1)]
    for k in range(groups + 1):
        rest = width - 4 * k
        poly = _available_exact_convolve(power, [Fraction(comb(rest, j)) for j in range(rest + 1)])
        for i in range(groups + 1):
            C[i][k] = Fraction(comb(groups, i), comb(width, 4 * i)) * (r * coverage) ** (4 * i) * poly[4 * i]
        power = _available_exact_convolve(power, local)
    return C


def available_tag_checks():
    F = Fraction
    m, b, p, r, alpha, coverage = F(3, 4), F(5, 4), F(1, 3), F(2, 5), F(2, 3), F(4, 5)
    entries = paths = terminal_checks = thinning_checks = uniform_checks = old_checks = 0
    for width, rounding in ((4, 'floor'), (8, 'floor'), (4, 'barycentric'), (8, 'barycentric')):
        for levels in ((0,), (0, 2), tuple(range(width + 1))):
            latent, physical, U, V, exact, terminal, end = _available_exact_factorization(
                width, m, b, p, r, alpha, levels, coverage, rounding)
            C, gpos = available_tag_matrix(width, float(m), float(b), float(p), float(r),
                                           float(alpha), levels, float(coverage), rounding=rounding)
            for i in range(len(latent)):
                expected_end = comb(width // 4, latent[i][0]) * (r * coverage * b) ** (4 * latent[i][0]) * (1 + m * p * coverage) ** (width - 4 * latent[i][0])
                assert end[i] == expected_end
                truth = AT_LD(end[i].numerator) / AT_LD(end[i].denominator)
                assert abs(gpos[i] + int(i == 0) - truth) <= AT_LD('2e-13') * max(1, abs(truth))
                terminal_checks += 1
                for j in range(len(latent)):
                    truth = AT_LD(exact[i][j].numerator) / AT_LD(exact[i][j].denominator)
                    assert abs(C[i, j] - truth) <= AT_LD('2e-13') * max(1, abs(truth))
                    entries += 1
            if levels == (0,):
                assert exact == _available_exact_old129(width, m, b, p, r, coverage)
                if 'complete_force_matrix' in globals():
                    prior, prior_end = complete_force_matrix(width, float(m), float(b), float(p), float(r), float(coverage))
                    assert np.max(abs(prior - C)) <= AT_LD('2e-13') * max(1, C.max())
                    assert np.max(abs(prior_end - gpos)) <= AT_LD('2e-13') * max(1, gpos.max())
                old_checks += len(latent) ** 2
            physical_transfer = [[sum(U[s][i] * V[i][d] for i in range(len(latent)))
                                  for d in range(len(physical))] for s in range(len(physical))]
            row = [F(int(s == 0)) for s in range(len(physical))]
            reduced = [F(int(i == 0)) for i in range(len(latent))]
            for _ in range(4):
                row = [sum(row[s] * physical_transfer[s][d] for s in range(len(physical)))
                       for d in range(len(physical))]
                assert sum(row[s] * terminal[s] for s in range(len(physical))) == sum(reduced[i] * end[i] for i in range(len(latent)))
                reduced = [sum(reduced[i] * exact[i][j] for i in range(len(latent))) for j in range(len(latent))]
                paths += 1
            zero, zero_end = available_tag_matrix(width, float(m), float(b), float(p), float(r),
                                                 float(alpha), levels, 0., rounding=rounding)
            expected_zero = np.zeros_like(zero)
            expected_zero[0, 0] = 1
            assert np.array_equal(zero, expected_zero) and np.all(zero_end == 0)
            if 'complete_force_region' in globals():
                Q, scale, end_positive, end_scale = complete_force_region(C, gpos)
                direct = np.linalg.matrix_power(C, 16)
                direct[0, 0] -= 1
                assert np.max(abs(Q * np.exp2(scale) - direct)) <= AT_LD('2e-12') * max(1, direct.max())
                gfull = gpos.copy()
                gfull[0] += 1
                direct_end = np.linalg.matrix_power(C, 15) @ gfull
                direct_end[0] -= 1
                assert np.max(abs(end_positive * np.exp2(end_scale) - direct_end)) <= AT_LD('2e-12') * max(1, direct_end.max())
    # Enumerate tag identities and retained subsets, independently of the transfer.
    for ell in range(9):
        for requested in ((0,), (0, 2, 4), tuple(range(9))):
            levels = tuple(g for g in requested if g <= max(ell, 0))
            weights = _available_exact_weights(ell, alpha, levels)
            retained = {}
            for tags in range(1 << ell):
                positions = [j for j in range(ell) if tags >> j & 1]
                count = len(positions)
                kept = max(g for g in levels if g <= count)
                probability = (1 - alpha) ** count * alpha ** (ell - count) / comb(count, kept)
                for subset in combinations(positions, kept):
                    mask = sum(1 << j for j in subset)
                    retained[mask] = retained.get(mask, F(0)) + probability
            assert sum(retained.values()) == 1
            for mask, probability in retained.items():
                kept = mask.bit_count()
                assert probability == weights[kept] / comb(ell, kept)
                uniform_checks += 1
            for holes in range(ell + 1):
                hole_mask = (1 << holes) - 1
                direct = sum(probability for mask, probability in retained.items() if not mask & hole_mask)
                formula = sum(weight * F(comb(ell - holes, g), comb(ell, g))
                              for g, weight in weights.items() if g <= ell - holes)
                assert direct == formula and direct >= alpha ** holes
                if requested == tuple(range(9)):
                    assert direct == alpha ** holes
                thinning_checks += 1
    rounded_gate_checks = rounded_weight_checks = 0
    for alpha_test in (F(0), F(1, 4), F(1, 2), F(2, 3), F(3, 4), F(1)):
        for levels in ((0,), (0, 2, 4, 8, 12), (0, 3, 7, 12), tuple(range(13))):
            floating = available_tag_probabilities(12, float(alpha_test), levels, 'barycentric')
            for ell in range(13):
                rounded = _available_exact_weights(ell, alpha_test, levels, 'barycentric')
                floored = _available_exact_weights(ell, alpha_test, levels)
                assert all(weight >= 0 for weight in rounded.values())
                for j, g in enumerate(levels):
                    truth = AT_LD(rounded[g].numerator) / AT_LD(rounded[g].denominator)
                    assert abs(floating[j, ell] - truth) <= AT_LD('2e-13')
                    rounded_weight_checks += 1
                for holes in range(ell + 1):
                    # P(g)=C(ell-g,D)/C(ell,D), decreasing and discretely convex.
                    rounded_survival = sum(weight * F(comb(ell - g, holes), comb(ell, holes))
                                           for g, weight in rounded.items() if g <= ell - holes)
                    floor_survival = sum(weight * F(comb(ell - g, holes), comb(ell, holes))
                                         for g, weight in floored.items() if g <= ell - holes)
                    assert alpha_test ** holes <= rounded_survival <= floor_survival
                    if levels == tuple(range(13)):
                        assert rounded_survival == alpha_test ** holes
                    rounded_gate_checks += 1
    run_checks = 0
    for q in (2, 4, 16, 1 << 32):
        for x in (F(1), F(1, 2), F(1, 16), F(1, q)):
            c = (q + 1) * x
            lam = 1 + c
            d = max(min(F(1), lam / (q + 1)), (1 + q * x * x) / lam)
            unforced = (F(1), F(0))
            forced = (F(0), q * x)
            for ell in range(1, 25):
                left, right = unforced
                unforced = (left + x * right, q * x * left + (1 + q * x) * right)
                assert unforced[0] + x * unforced[1] <= d * lam ** ell
                assert forced[0] + x * forced[1] <= c * lam ** (ell - 1) * d
                assert sum(unforced) <= lam ** ell
                assert sum(forced) <= c * lam ** (ell - 1)
                left, right = forced
                forced = (left + x * right, q * x * left + (1 + q * x) * right)
                run_checks += 1
    print(dict(available_tag_factorization_entries=entries, exact_composed_paths=paths,
               terminal_checks=terminal_checks, tag_thinning_checks=thinning_checks,
               uniform_retained_subset_checks=uniform_checks, exact_old129_entries=old_checks,
               barycentric_probability_checks=rounded_weight_checks,
               exact_barycentric_gate_sandwich_checks=rounded_gate_checks,
               exact_joint_closed_and_free_run_checks=run_checks,
               zero_operator_and_optional_region_checks=True, full_code_certificate=False), flush=True)


def available_tag_field_weights(theta, logz):
    """One joint run bound; alpha is not scaled by rank/availability markers.

    For H=[[1,qx],[x,1+qx]], the left row e0*H**ell has second/first
    ratio at most q+1, an invariant of its increasing Mobius update.
    With w=(1,x), this implies u*(H*w-lambda*w)<=-u[0]*x*x.
    Thus lambda**(-ell)*e0*H**ell*w decreases and is at most its
    length-one value (1+q*x*x)/lambda.  A forced closed run uses e;
    taking their maximum is a single valid factor for every closure.
    """
    theta, logz = AT_LD(theta), AT_LD(logz)
    if not 0 < theta <= 1 or logz > 0:
        raise ValueError('Require 0 < theta <= 1 and log2(z) <= 0')
    q = AT_LD(2) ** 32
    x = np.exp2(logz / theta)
    c = (q + 1) * x
    lam = 1 + c
    forced_closed = min(AT_LD(1), lam / (q + 1))
    unforced_closed = (1 + q * x * x) / lam
    d = max(forced_closed, unforced_closed)
    return c ** theta, lam ** theta, d ** theta, d


def available_tag_value(h=2048, theta=1/16, logz=-.4375, tilts=3, levels=(0, 32),
                        logcoverage=0., logrankgap=0., logavailable=0., logoutside=0., progress=False,
                        rounding='floor', profile=False, closuregap=0., closure_census_j=0,
                        rank_census_gap=0., rank_census_j=0,logexclusion=0.):
    """Complete-history FLOAT diagnostic, retaining all terminal endpoints."""
    if 'complete_force_region' not in globals() or 'complete_force_tail' not in globals():
        raise RuntimeError('Inject complete_force_region and complete_force_tail from the prior helper first')
    if tilts < 1 or not 0 < h < 8192:
        raise ValueError('Require positive tilt count and 0<h<8192')
    if logcoverage > 0 or logoutside > 0:
        raise ValueError('The upper coverage/outside caps require logcoverage<=0 and logoutside<=0')
    if logrankgap < 0 or logavailable < 0:
        raise ValueError('Require logrankgap>=0 and uncompensated logavailable>=0; use the matrix API for raw moments')
    m, b, alpha, d = available_tag_field_weights(theta, logz)
    if closuregap < 0 or rank_census_gap < 0 or logexclusion < 0:
        raise ValueError('Closure-remainder marker must be nonnegative')
    # Census class complements, applied canonically before tag rounding.
    # Apply it before the positive relaxation and before thinning tags.
    closure_cap = 16 * AT_W + available_census_constant(closure_census_j) + 91
    rank_census_cap = available_census_constant(rank_census_j)+44
    alpha *= np.exp2(AT_LD(-(15-2*closure_census_j)*closuregap
                          +2*rank_census_j*rank_census_gap))
    if alpha>1:
        raise ValueError('Census markers make alpha>1; survival tags no longer apply')
    prepared = available_tag_prepare(AT_WIDTH, alpha, levels, rounding)
    m *= np.exp2(AT_LD(logrankgap + logavailable + 16 * closuregap + 32*rank_census_gap+2*logexclusion))
    b *= np.exp2(AT_LD(logrankgap + logavailable + 16 * closuregap + 32*rank_census_gap+logexclusion))
    field = np.exp2(AT_LD(32 * (theta + closuregap + rank_census_gap)))
    rows = 8192
    actual = AT_LD(h) / rows
    low = actual / (field - (field - 1) * actual)
    choose = (lgamma(rows + 1) - lgamma(h + 1) - lgamma(rows - h + 1)) / float(np.log(AT_LD(2)))
    logs = endlogs = None
    for tilt_index, fraction in enumerate(np.linspace(0., 1., tilts)):
        started = time.monotonic()
        odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
        p = odds / (1 + odds)
        logpi = choose + h * np.log2(p) + (rows - h) * np.log2(1 - p)
        r = p + (1 - p) * np.exp2(AT_LD(logoutside - logrankgap)) / field
        C, terminal_positive = available_tag_matrix(AT_WIDTH, m, b, p, r, alpha, levels,
                                                    np.exp2(AT_LD(logcoverage-2*logexclusion)), prepared, rounding)
        if progress:
            print(dict(available_tag_stage='matrix', tilt=tilt_index + 1, dimension=len(C),
                       seconds=round(time.monotonic() - started, 3)), flush=True)
        Q, scale, end, end_scale = complete_force_region(C, terminal_positive)
        if not np.all(Q > 0) or not np.all(end > 0):
            raise ArithmeticError('Unexpected regional zero/underflow; do not treat it as an upper bound')
        candidate = np.log2(Q) + scale - logpi
        candidate_end = np.log2(end) + end_scale - logpi
        logs = candidate if logs is None else np.minimum(logs, candidate)
        endlogs = candidate_end if endlogs is None else np.minimum(endlogs, candidate_end)
        if progress:
            print(dict(available_tag_stage='region', tilt=tilt_index + 1,
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    scale = logs.max()
    Q = np.exp2(logs - scale)
    norm = Q.sum(axis=1).max()
    Q /= norm
    scale += np.log2(norm)
    end_scale = endlogs.max()
    end = np.exp2(endlogs - end_scale)
    value = complete_force_tail(Q, scale, end, end_scale)
    bound = float(theta - 32 * theta - AT_W * logz + choose - 32 * theta * 128 * h + value
                  - (1009255 - 128 * h) * logoutside - 1736703 * logcoverage
                  + (65536 - 128 * h) * logrankgap
                  + (closure_cap - 32 * 128 * h) * closuregap
                  + (rank_census_cap-32*128*h)*rank_census_gap+((1<<21)-512)*logexclusion)
    detail = dict(dimension=len(prepared['states']), alpha=float(alpha), closure_d=float(d),
                  levels=tuple(levels), tilts=tilts, rounding=rounding, closuregap=closuregap,
                  closure_census_j=closure_census_j,rank_census_gap=rank_census_gap,rank_census_j=rank_census_j,logexclusion=logexclusion)
    if profile:
        # Diagnostic for the all-positive regional term of the upper bound.
        # This is NOT a sampled-code trajectory or a canonical profile.
        forward = np.zeros(len(Q), dtype=AT_LD); forward[0] = 1
        backward = end.copy(); backward /= backward.max()
        for _ in range(128):
            forward = forward @ Q; forward /= forward.sum()
        for _ in range(127):
            backward = Q @ backward; backward /= backward.max()
        marginal = forward * backward; marginal /= marginal.sum()
        top = np.argsort(marginal)[-10:][::-1]
        detail['relaxed_midpoint_top_i_g'] = [(prepared['states'][j], float(marginal[j])) for j in top]
        detail['relaxed_midpoint_mean_optional'] = float(sum(4*i*marginal[j] for j,(i,g) in enumerate(prepared['states'])))
        detail['relaxed_midpoint_tag_mass'] = [(g,float(sum(marginal[j] for j,(_,gg) in enumerate(prepared['states']) if gg==g))) for g in levels]
    return bound, detail


def available_tags_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--tilts', type=int, default=3)
    parser.add_argument('--levels', default='0,32')
    parser.add_argument('--rounding', choices=('floor', 'barycentric'), default='floor')
    parser.add_argument('--profile', action='store_true')
    parser.add_argument('--exclusion', type=float, default=0.)
    parser.add_argument('--closuregap', type=float, default=0.)
    parser.add_argument('--closure-census-j', type=int, default=0)
    parser.add_argument('--rank-census-gap', type=float, default=0.)
    parser.add_argument('--rank-census-j', type=int, default=0)
    parser.add_argument('--rankgap', type=float, default=0.)
    parser.add_argument('--coverage', type=float, default=0.)
    parser.add_argument('--outside', type=float, default=0.)
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--logz', type=float, default=-.4375)
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
            available_tag_checks()
            return
        levels = tuple(int(value) for value in args.levels.split(','))
        started = time.monotonic()
        bound, detail = available_tag_value(args.h, args.theta, args.logz, args.tilts, levels,
                                            progress=True, rounding=args.rounding, profile=args.profile,
                                            closuregap=args.closuregap, logrankgap=args.rankgap,
                                            logcoverage=args.coverage, logoutside=args.outside,
                                            closure_census_j=args.closure_census_j,
                                            rank_census_gap=args.rank_census_gap,rank_census_j=args.rank_census_j,
                                            logexclusion=args.exclusion)
        print(dict(available_tag_rows=args.h, theta=args.theta, logz=args.logz,
                   log2_diagnostic=bound, seconds=round(time.monotonic() - started, 3),
                   certificate=False, **detail), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    available_tags_main()
