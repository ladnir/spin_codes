"""Bounded FLOAT partition-averaged nested-selector potential diagnostic.

For |U|=u and |Aprev|=a, the one-step potential ratio is exactly

 [X^a Y^(u-a)] P(X,Y)^g / (multinomial(n;a,u-a,n-u) A^u B^a),
 P=(1+dX+Y)^4+z^4 A^4[(1+wB+beta B(X+Y))^4-1], g=n/4.

All fifteen coefficients are built positively. Each row uses its own
coefficient Chernoff variables. The three boundary faces are extracted
exactly; no finite tilt substitutes for zero or infinity. Bounded batched
Newton only selects positive tilts, so stopping early preserves the
mathematical upper-bound direction. FLOAT arithmetic is not directed and
the reported values are diagnostics, not probability certificates.

By default, rows satisfy u/4<=a<=u: every selected current quartet meets
the previous allowed set. This invariant is preserved by the required
birth in an old-empty selected quartet. The next partition is fresh, so
its averaged action depends only on (u,a). --all-rows retains the larger
0<=a<=u control. No transition matrix or outer marker search is built.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import comb, lgamma, log, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np


PA_WIDTH, PA_UPDATES, PA_M, PA_W = 512, 4096, 8192, 943718
PA_PROFILE = (4096, 819200, 838861, 65536, 1081344)
PA_LOGS = (-0.9301908045807551, -1.343939327491652,
           0.17125048238006038, -0.45762697461224744,
           0.0, 0.9488661799848431)
_PA_LN2 = log(2)
_PA_EXP = np.asarray([(i, j) for i in range(5) for j in range(5-i)], dtype=np.float64)
_PA_MULT = np.asarray([comb(4, int(i))*comb(4-int(i), int(j))
                      for i, j in _PA_EXP], dtype=np.float64)


def partition_averaged_polynomial(log_beta, log_d, log_z, log_w, log_A, log_B):
    """Return the fifteen positive coefficients in log2 form, without subtraction."""
    markers = np.asarray((log_beta, log_d, log_z, log_w, log_A, log_B), dtype=np.float64)
    if not np.all(np.isfinite(markers)) or np.max(np.abs(markers)) > 512:
        raise ValueError('Require six finite marker logs of absolute value at most512')
    wb = log_w + log_B
    constant_births = np.asarray([log2(comb(4, ell)) + ell*wb for ell in range(1, 5)])
    maximum = float(constant_births.max())
    positive_constant = maximum + log2(float(np.exp2(constant_births-maximum).sum()))
    degree = _PA_EXP.sum(axis=1)
    selected = (np.log2(_PA_MULT) + 4*(log_z+log_A)
                + degree*(log_beta+log_B) + (4-degree)*np.logaddexp2(0, wb))
    selected[0] = 4*(log_z+log_A) + positive_constant
    unselected = np.log2(_PA_MULT) + _PA_EXP[:, 0]*log_d
    return np.logaddexp2(unselected, selected)


def _pa_evaluate(tilts, targets, exponents, coefficients, derivatives=True):
    affine = tilts @ exponents.T + coefficients
    maximum = affine.max(axis=1)
    weights = np.exp2(affine-maximum[:, None])
    totals = weights.sum(axis=1)
    weights /= totals[:, None]
    value = maximum + np.log2(totals) - np.sum(tilts*targets, axis=1)
    if not derivatives:
        return value
    means = weights @ exponents
    gradient = means-targets
    e0 = exponents[:, 0]
    h00 = _PA_LN2*np.maximum(0.0, weights @ (e0*e0)-means[:, 0]**2)
    if exponents.shape[1] == 1:
        return value, gradient, (h00,)
    e1 = exponents[:, 1]
    h11 = _PA_LN2*np.maximum(0.0, weights @ (e1*e1)-means[:, 1]**2)
    h01 = _PA_LN2*(weights @ (e0*e1)-means[:, 0]*means[:, 1])
    return value, gradient, (h00, h01, h11)


def _pa_newton(targets, exponents, coefficients, initial, max_iterations):
    """Independent convex row optimizations, batched; every stored iterate is valid."""
    tilts = np.clip(initial.copy(), -256.0, 256.0)
    best = _pa_evaluate(tilts, targets, exponents, coefficients, False)
    best_tilts = tilts.copy()
    active = np.arange(len(targets))
    iterations = evaluations = 0
    for iteration in range(max_iterations):
        if not len(active):
            break
        value, gradient, hessian = _pa_evaluate(tilts[active], targets[active], exponents, coefficients)
        evaluations += len(active)
        keep = np.max(np.abs(gradient), axis=1) > 2e-10
        active, value, gradient = active[keep], value[keep], gradient[keep]
        hessian = tuple(entry[keep] for entry in hessian)
        if not len(active):
            break
        if exponents.shape[1] == 1:
            direction = -gradient/(hessian[0][:, None]+1e-11)
        else:
            h00, h01, h11 = hessian
            ridge = 1e-11+1e-10*(h00+h11)
            h00, h11 = h00+ridge, h11+ridge
            determinant = np.maximum(h00*h11-h01*h01, 1e-22)
            direction = np.column_stack(((-h11*gradient[:, 0]+h01*gradient[:, 1])/determinant,
                                         (h01*gradient[:, 0]-h00*gradient[:, 1])/determinant))
        direction *= np.minimum(1.0, 8.0/np.maximum(1e-300, np.max(np.abs(direction), axis=1)))[:, None]
        direction = np.clip(tilts[active]+direction, -256.0, 256.0)-tilts[active]
        descent = np.sum(gradient*direction, axis=1)
        bad = ~np.isfinite(descent) | (descent >= -1e-16)
        if np.any(bad):
            fallback = -gradient[bad]/np.maximum(1.0, np.max(np.abs(gradient[bad]), axis=1))[:, None]
            direction[bad] = np.clip(tilts[active[bad]]+fallback, -256.0, 256.0)-tilts[active[bad]]
            descent[bad] = np.sum(gradient[bad]*direction[bad], axis=1)
        movable = descent < -1e-16
        active, value, direction, descent = active[movable], value[movable], direction[movable], descent[movable]
        pending = np.arange(len(active))
        accepted = np.zeros(len(active), dtype=bool)
        alpha = np.ones(len(active))
        for _ in range(18):
            if not len(pending):
                break
            ids = active[pending]
            candidate = tilts[ids]+alpha[pending, None]*direction[pending]
            trial = _pa_evaluate(candidate, targets[ids], exponents, coefficients, False)
            evaluations += len(ids)
            improved = trial < best[ids]
            best[ids[improved]], best_tilts[ids[improved]] = trial[improved], candidate[improved]
            good = trial <= value[pending]+1e-4*alpha[pending]*descent[pending]
            tilts[ids[good]] = candidate[good]
            accepted[pending[good]] = True
            pending = pending[~good]
            alpha[pending] *= 0.5
        active = active[accepted]
        iterations = iteration+1
    _, residual, _ = _pa_evaluate(best_tilts, targets, exponents, coefficients)
    return best, best_tilts, dict(iterations=iterations, row_evaluations=evaluations,
                                 maximum_mean_residual=float(np.max(np.abs(residual))) if len(residual) else 0.0)


def _pa_faces(width, rows):
    """Minimal exposed polynomial faces; each description has a finite optimizer."""
    u, a = rows[:, 0], rows[:, 1]
    interior = (a > 0) & (a < u) & (u < width)
    yield 'interior', np.flatnonzero(interior), np.arange(15), _PA_EXP, np.column_stack((a, u-a))
    yield 'a_zero', np.flatnonzero((a == 0) & (u > 0) & (u < width)), np.flatnonzero(_PA_EXP[:, 0] == 0), None, (u-a)[:, None]
    yield 'a_u', np.flatnonzero((a == u) & (u > 0) & (u < width)), np.flatnonzero(_PA_EXP[:, 1] == 0), None, a[:, None]
    yield 'u_width', np.flatnonzero((u == width) & (a > 0) & (a < width)), np.flatnonzero(_PA_EXP.sum(axis=1) == 4), None, a[:, None]


def partition_averaged_rows(log_beta, log_d, log_z, log_w, log_A, log_B,
                            width=512, physical_rows=True, max_iterations=40, chunk_size=2048):
    """Reusable row arrays and worst ratios; no transition matrix is constructed.

    row_log2_upper is the raw coefficient-Chernoff comparison. The global
    log2_growth_upper also takes the valid old pointwise comparison as a
    fallback. Both are FLOAT diagnostics. Faces/vertices are handled exactly
    at the polynomial level, before floating coefficient evaluation.
    """
    width, max_iterations, chunk_size = map(index, (width, max_iterations, chunk_size))
    if width < 4 or width % 4 or width > 512 or not 1 <= max_iterations <= 96 or not 1 <= chunk_size <= 8192:
        raise ValueError('Require width4..512 divisible4, iterations1..96, chunk1..8192')
    coefficients = partition_averaged_polynomial(log_beta, log_d, log_z, log_w, log_A, log_B)
    groups = width//4
    rows = np.asarray([(u, a) for u in range(0, width+1, 4)
                       for a in range(u//4 if physical_rows else 0, u+1)], dtype=np.int64)
    count = len(rows)
    raw = np.full(count, np.nan)
    tilts = np.full((count, 2), np.nan)
    face_names = np.full(count, '', dtype='<U12')
    statistics = dict(row_evaluations=0, maximum_mean_residual=0.0, maximum_iterations=0)
    for name, indices, terms, exponents, target_counts in _pa_faces(width, rows):
        if not len(indices):
            continue
        if exponents is None:
            coordinate = 1 if name == 'a_zero' else 0
            exponents = _PA_EXP[terms, coordinate, None]
        selected_coefficients = coefficients[terms]
        for start in range(0, len(indices), chunk_size):
            subset = indices[start:start+chunk_size]
            target = target_counts[subset]/groups
            if exponents.shape[1] == 2:
                missing = width-rows[subset, 0]
                initial = np.log2(target_counts[subset]/missing[:, None])
            else:
                initial = np.log2(target_counts[subset]/(width-target_counts[subset]))
            values, selected_tilts, info = _pa_newton(target, exponents, selected_coefficients, initial, max_iterations)
            raw[subset] = groups*values
            tilts[subset, :exponents.shape[1]] = selected_tilts
            face_names[subset] = name
            statistics['row_evaluations'] += info['row_evaluations']
            statistics['maximum_mean_residual'] = max(statistics['maximum_mean_residual'], info['maximum_mean_residual'])
            statistics['maximum_iterations'] = max(statistics['maximum_iterations'], info['iterations'])
    for u, a, exponent in ((0, 0, (0, 0)), (width, 0, (0, 4)), (width, width, (4, 0))):
        positions = np.flatnonzero((rows[:, 0] == u) & (rows[:, 1] == a))
        if len(positions):
            term = int(np.flatnonzero(np.all(_PA_EXP == exponent, axis=1))[0])
            raw[positions] = groups*coefficients[term]
            face_names[positions] = 'vertex'
    if not np.all(np.isfinite(raw)) or np.any(face_names == ''):
        raise ArithmeticError('Uncovered or non-finite coefficient row')
    log_factorial = np.zeros(width+1)
    log_factorial[1:] = np.cumsum(np.log2(np.arange(1, width+1, dtype=np.float64)))
    u, a = rows[:, 0], rows[:, 1]
    multinomial = log_factorial[width]-log_factorial[a]-log_factorial[u-a]-log_factorial[width-u]
    raw -= multinomial+u*log_A+a*log_B
    local = coefficients-np.log2(_PA_MULT)-_PA_EXP.sum(axis=1)*log_A-_PA_EXP[:, 0]*log_B
    pointwise = groups*float(local.max())
    order = np.argsort(raw)[-8:][::-1]
    top = [dict(u=int(u[i]), a=int(a[i]), face=str(face_names[i]), log2_ratio_upper=float(raw[i]),
                finite_log_tilts=[float(x) for x in tilts[i] if np.isfinite(x)]) for i in order]
    return dict(width=width, physical_rows=bool(physical_rows), row_count=count,
                log2_growth_upper=min(float(raw.max()), pointwise),
                raw_chernoff_log2_growth_upper=float(raw.max()), pointwise_log2_growth_upper=pointwise,
                pointwise_fallback_used=bool(pointwise < raw.max()), largest_rows=top,
                face_counts={name: int(np.sum(face_names == name)) for name in np.unique(face_names)},
                optimizer=statistics, rows=rows, row_log2_upper=raw, log_tilts=tilts,
                floating_diagnostic_only=True, full_code_certificate=False)


def _pa_choose(n, k):
    if not 0 <= k <= n:
        return float('-inf')
    if k in (0, n):
        return 0.0
    return (lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/_PA_LN2


def partition_averaged_bound(profile=PA_PROFILE, death_upper=None, logs=PA_LOGS, **row_options):
    """Low-death wrapper: logs=(logz,logw,logA,logB,logb,logk), b,k>=1."""
    h, rank, volume, births, footprint = map(index, profile)
    logs = np.asarray(logs, dtype=np.float64)
    if logs.shape != (6,) or not np.all(np.isfinite(logs)) or np.any(logs[4:] < 0):
        raise ValueError('Require six finite logs, with logb/logk nonnegative')
    N = PA_WIDTH*PA_UPDATES
    if not (1 <= h <= PA_M and 128*h <= rank <= 1 << 20 and 1 <= births <= volume <= footprint <= N):
        raise ValueError('Invalid actual-state profile')
    G, m, CH = footprint-births, volume-births, footprint-rank+128*h
    if G % 4 or G > N-PA_WIDTH:
        raise ValueError('Invalid live-selector volume')
    lower, upper = max(0, births-PA_WIDTH), min(footprint-volume, PA_W//2)
    D0 = min(PA_W-volume if death_upper is None else index(death_upper), upper)
    if D0 < lower or not (129*h <= CH <= 256*h and 0 <= CH-births <= G and 0 <= 256*h-CH <= N-footprint):
        return dict(profile=tuple(profile), status='empty', full_code_certificate=False)
    log_z, log_w, log_A, log_B, log_b, log_k = logs
    if 'width' in row_options and row_options['width'] != PA_WIDTH:
        raise ValueError('Profile wrapper requires width512')
    data = partition_averaged_rows(log_b+log_k, -log_k, log_z, log_w, log_A, log_B, **row_options)
    route = min(0.0, _pa_choose(G, CH-births)+_pa_choose(N-footprint, 256*h-CH)-255*_pa_choose(PA_M, h))
    endpoint = max(float(np.logaddexp2(0, log_w)), log_b+log_k-log_A+max(0.0, -log_B))
    correction = route-G*log_z-births*log_w-m*log_b-(m-D0)*log_k
    raw = (PA_UPDATES-1)*data['log2_growth_upper']+PA_WIDTH*endpoint+correction
    summary = {key: value for key, value in data.items() if key not in ('rows', 'row_log2_upper', 'log_tilts')}
    summary.update(profile=tuple(map(int, profile)), status='low_actual_death_class',
                   death_interval=(lower, D0), logs=list(map(float, logs)), route_log2=route,
                   log2_endpoint_per_coordinate=float(endpoint), external_log2_correction=float(correction),
                   raw_log2_probability_bound=float(raw), births_integrated=True,
                   H_union_already_in_route=True, free_endpoint=True)
    return summary


def _pa_fraction_polynomial(beta, d, z, w, A, B):
    polynomial = {}
    for i in range(5):
        for j in range(5-i):
            r = i+j
            births = sum(Fraction(comb(4-r, ell))*(w*B)**ell for ell in range(int(r == 0), 5-r))
            polynomial[i, j] = comb(4, i)*comb(4-i, j)*(d**i+z**4*A**4*(beta*B)**r*births)
    return polynomial


def _pa_power(polynomial, exponent):
    result = {(0, 0): Fraction(1)}
    for _ in range(exponent):
        following = {}
        for (i, j), value in result.items():
            for (k, ell), coefficient in polynomial.items():
                key = (i+k, j+ell)
                following[key] = following.get(key, Fraction(0))+value*coefficient
        result = following
    return result


def _pa_subsets(mask):
    subset = mask
    while True:
        yield subset
        if subset == 0:
            break
        subset = (subset-1)&mask


def _pa_physical(old, previous, partition, beta, d, z, w, A, B):
    total, transitions = Fraction(0), 0
    for selected in range(1 << len(partition)):
        blocks = [block for j, block in enumerate(partition) if selected >> j & 1]
        new = sum(blocks)
        for births in _pa_subsets(new & ~old):
            if any(not (block&old) and not (block&births) for block in blocks):
                continue
            allowed = (old&new)|births
            assert 4*allowed.bit_count() >= new.bit_count()
            assert all(block&allowed for block in blocks)
            total += (d**(previous&~allowed).bit_count()*z**new.bit_count()
                      *beta**(old&new).bit_count()*w**births.bit_count()
                      *A**new.bit_count()*B**allowed.bit_count())
            transitions += 1
    return total, transitions


def partition_averaged_checks():
    physical = transitions = faces = rational_chernoff = floating = 0
    settings = ((Fraction(2), Fraction(1, 2), Fraction(1, 2), Fraction(1, 8), Fraction(2), Fraction(1, 2)),
                (Fraction(3, 2), Fraction(2, 3), Fraction(3, 4), Fraction(1, 4), Fraction(3, 2), Fraction(3, 4)))
    for beta, d, z, w, A, B in settings:
        polynomial = _pa_fraction_polynomial(beta, d, z, w, A, B)
        logs = tuple(log2(float(x)) for x in (beta, d, z, w, A, B))
        computed = partition_averaged_polynomial(*logs)
        for (i, j), coefficient in polynomial.items():
            term = int(np.flatnonzero(np.all(_PA_EXP == (i, j), axis=1))[0])
            assert abs(computed[term]-log2(float(coefficient))) < 1e-12
        for width in (4, 8):
            groups = width//4
            powered = _pa_power(polynomial, groups)
            partitions = [(15,)] if width == 4 else [
                (first, 255^first) for rest in combinations(range(1, 8), 3)
                for first in (1+sum(1 << i for i in rest),)]
            numeric = partition_averaged_rows(*logs, width=width, physical_rows=False, max_iterations=48, chunk_size=7)
            for row_index, (u, a) in enumerate(numeric['rows']):
                u, a = int(u), int(a)
                old, previous = (1 << u)-1, (1 << a)-1
                averaged = Fraction(0)
                for partition in partitions:
                    value, checked = _pa_physical(old, previous, partition, beta, d, z, w, A, B)
                    averaged += value
                    transitions += checked
                averaged /= len(partitions)
                normalization = comb(width, a)*comb(width-a, u-a)
                assert averaged == powered[a, u-a]/normalization
                physical += 1
                face = {key: value for key, value in polynomial.items()
                        if (a != 0 or key[0] == 0) and (a != u or key[1] == 0)
                        and (u != width or sum(key) == 4)
                        and (u != 0 or sum(key) == 0)}
                assert _pa_power(face, groups).get((a, u-a), 0) == powered[a, u-a]
                faces += 1
                X, Y = Fraction(3, 2), Fraction(2, 3)
                bound = sum(value*X**i*Y**j for (i, j), value in face.items())**groups/(X**a*Y**(u-a))
                assert bound >= powered[a, u-a]
                rational_chernoff += 1
                actual_log = log2(float(averaged/(A**u*B**a)))
                assert numeric['row_log2_upper'][row_index] >= actual_log-2e-10
                assert numeric['log2_growth_upper'] >= actual_log-2e-10
                floating += 1
            filtered = partition_averaged_rows(*logs, width=width, physical_rows=True, max_iterations=1)
            assert np.all(4*filtered['rows'][:, 1] >= filtered['rows'][:, 0])
            assert filtered['row_count'] == sum(3*i+1 for i in range(groups+1))
    # A stopped optimizer still evaluates an ordinary coefficient upper bound.
    target = np.asarray(((1.0, 2.0), (0.25, 0.5)))
    tilt = np.asarray(((0.2, -0.3), (-0.7, 0.4)))
    coefficients = partition_averaged_polynomial(*logs)
    values, gradient, _ = _pa_evaluate(tilt, target, _PA_EXP, coefficients)
    for coordinate in range(2):
        delta = np.zeros_like(tilt)
        delta[:, coordinate] = 1e-5
        plus = _pa_evaluate(tilt+delta, target, _PA_EXP, coefficients, False)
        minus = _pa_evaluate(tilt-delta, target, _PA_EXP, coefficients, False)
        assert np.max(np.abs((plus-minus)/(2e-5)-gradient[:, coordinate])) < 1e-7
    print(dict(exact_physical_population_rows=physical, exact_birth_transition_checks=transitions,
               exact_face_extractions=faces, exact_rational_chernoff_checks=rational_chernoff,
               float_row_comparisons=floating, gradient_coordinates=2,
               physical_row_count_512=sum(3*i+1 for i in range(129)),
               all_row_count_512=sum(4*i+1 for i in range(129)),
               full_code_certificate=False), flush=True)


def partition_averaged_main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--checks', action='store_true')
    modes.add_argument('--probe', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=PA_PROFILE)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--logs', type=float, nargs=6, default=PA_LOGS)
    parser.add_argument('--iterations', type=int, default=40)
    parser.add_argument('--chunk-size', type=int, default=2048)
    parser.add_argument('--all-rows', action='store_true')
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            partition_averaged_checks()
        else:
            print(partition_averaged_bound(args.profile, args.death_upper, args.logs,
                                          physical_rows=not args.all_rows,
                                          max_iterations=args.iterations, chunk_size=args.chunk_size), flush=True)
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    partition_averaged_main()
