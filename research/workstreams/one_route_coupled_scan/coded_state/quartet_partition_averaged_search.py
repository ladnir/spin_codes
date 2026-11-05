"""Bounded FLOAT outer search for the partition-averaged nested potential.

This helper leaves the averaged-row implementation unchanged. It expands
its polynomial into49 positive monomials:15 unselected and34 selected.
For each active population row, Chernoff tilts are minimized on the exact
polynomial face. Envelope gradients are conditional monomial means;
Hessians use the covariance Schur complement in the tilt coordinates.
They are FLOAT derivatives at approximately optimized inner tilts.

An active set only selects candidate markers. Periodic and final complete
physical-row scans supply every reported global bound. The search is
bounded in a stated six-dimensional box and makes no optimality claim.
No SciPy, large transition matrix, or numerical probability certificate.
"""

import argparse
from fractions import Fraction
from math import comb, log, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np

try:
    partition_averaged_rows
except NameError:
    from quartet_partition_averaged_potential import (
        partition_averaged_rows, partition_averaged_bound, partition_averaged_polynomial,
        _pa_newton, _pa_choose, PA_PROFILE, PA_LOGS, PA_WIDTH, PA_UPDATES, PA_M, PA_W)


PS_BOUNDS = ((-8.0, 4.0), (-24.0, 8.0), (-8.0, 8.0),
             (-12.0, 8.0), (0.0, 12.0), (0.0, 12.0))
_PS_LN2 = log(2)


def _ps_model():
    tilts, outer, multiplicities, identities = [], [], [], []
    for i in range(5):
        for j in range(5-i):
            tilts.append((i, j))
            outer.append((0, 0, 0, 0, 0, -i))
            multiplicities.append(comb(4, i)*comb(4-i, j))
            identities.append((False, i, j, 0))
    for i in range(5):
        for j in range(5-i):
            for ell in range(5-i-j):
                if i+j+ell == 0:
                    continue
                tilts.append((i, j))
                outer.append((4, ell, 4, i+j+ell, i+j, i+j))
                multiplicities.append(comb(4, i)*comb(4-i, j)*comb(4-i-j, ell))
                identities.append((True, i, j, ell))
    return (np.asarray(tilts, dtype=np.float64), np.asarray(outer, dtype=np.float64),
            np.asarray(multiplicities, dtype=np.float64), tuple(identities))


_PS_TILT, _PS_OUTER, _PS_MULT, _PS_IDENTITIES = _ps_model()
_PS_COEFFICIENTS = np.log2(_PS_MULT)


def _ps_face(width, u, a):
    keep = np.ones(49, dtype=bool)
    if a == 0:
        keep &= _PS_TILT[:, 0] == 0
    if a == u:
        keep &= _PS_TILT[:, 1] == 0
    if u == width:
        keep &= _PS_TILT.sum(axis=1) == 4
    if u == 0:
        keep &= _PS_TILT.sum(axis=1) == 0
    if u == 0 or (u == width and a in (0, width)):
        coordinates = ()
    elif a == 0:
        coordinates = (1,)
    elif a == u or u == width:
        coordinates = (0,)
    else:
        coordinates = (0, 1)
    return np.flatnonzero(keep), coordinates


def partition_averaged_envelopes(logs, rows, width=512, cache=None, inner_iterations=32):
    """Per-quartet optimized coefficient values, gradients and Schur Hessians.

    Values are positive-coefficient bounds at finite inner tilts. Envelope
    derivatives require an optimized tilt and are therefore approximate;
    they guide a FLOAT search only. Vertices have no tilt variables.
    """
    logs = np.asarray(logs, dtype=np.float64)
    rows = np.asarray(rows, dtype=np.int64)
    width, inner_iterations = index(width), index(inner_iterations)
    if logs.shape != (6,) or not np.all(np.isfinite(logs)) or rows.ndim != 2 or rows.shape[1] != 2:
        raise ValueError('Require six finite logs and(u,a) rows')
    if width < 4 or width % 4 or width > 512 or not 1 <= inner_iterations <= 96:
        raise ValueError('Invalid width or inner iteration budget')
    if np.any(rows[:, 0] % 4) or np.any(rows < 0) or np.any(rows[:, 0] > width) or np.any(rows[:, 1] > rows[:, 0]):
        raise ValueError('Invalid population row')
    if cache is None:
        cache = {}
    groups = width//4
    coefficients = _PS_COEFFICIENTS+_PS_OUTER @ logs
    values = np.empty(len(rows))
    gradients = np.empty((len(rows), 6))
    hessians = np.empty((len(rows), 6, 6))
    factorial = np.zeros(width+1)
    factorial[1:] = np.cumsum(np.log2(np.arange(1, width+1, dtype=np.float64)))
    faces = {}
    for row_index, (u, a) in enumerate(rows):
        terms, coordinates = _ps_face(width, int(u), int(a))
        key = (tuple(terms), coordinates)
        faces.setdefault(key, []).append(row_index)
    maximum_residual = 0.0
    for (term_tuple, coordinates), indices in faces.items():
        terms = np.asarray(term_tuple)
        ids = np.asarray(indices)
        selected = rows[ids]
        u, a = selected[:, 0], selected[:, 1]
        full_targets = np.column_stack((a, u-a))/groups
        outer = _PS_OUTER[terms]
        coeff = coefficients[terms]
        if coordinates:
            exponents = _PS_TILT[terms][:, coordinates]
            targets = full_targets[:, coordinates]
            initial = np.empty_like(targets)
            for position, (old_u, old_a) in enumerate(selected):
                key = (width, int(old_u), int(old_a))
                if key in cache:
                    initial[position] = cache[key]
                elif len(coordinates) == 2:
                    initial[position] = np.log2(np.asarray((old_a, old_u-old_a))/(width-old_u))
                else:
                    count = old_a if coordinates[0] == 0 else old_u-old_a
                    initial[position, 0] = log2(count/(width-count))
            optimized, tilts, info = _pa_newton(targets, exponents, coeff, initial, inner_iterations)
            maximum_residual = max(maximum_residual, info['maximum_mean_residual'])
            for position, (old_u, old_a) in enumerate(selected):
                cache[width, int(old_u), int(old_a)] = tilts[position].copy()
            affine = tilts @ exponents.T+coeff
        else:
            affine = np.broadcast_to(coeff, (len(ids), len(coeff)))
        maximum = affine.max(axis=1)
        weights = np.exp2(affine-maximum[:, None])
        weights /= weights.sum(axis=1)[:, None]
        if not coordinates:
            optimized = maximum+np.log2(np.exp2(affine-maximum[:, None]).sum(axis=1))
        means = weights @ outer
        covariance = np.einsum('nt,ti,tj->nij', weights, outer, outer)-means[:, :, None]*means[:, None, :]
        if coordinates:
            mean_tilt = weights @ exponents
            cross = np.einsum('nt,ti,tj->nij', weights, outer, exponents)-means[:, :, None]*mean_tilt[:, None, :]
            tilt_cov = np.einsum('nt,ti,tj->nij', weights, exponents, exponents)-mean_tilt[:, :, None]*mean_tilt[:, None, :]
            if len(coordinates) == 1:
                covariance -= cross @ np.swapaxes(cross, 1, 2)/np.maximum(tilt_cov[:, 0, 0], 1e-14)[:, None, None]
            else:
                h00 = np.maximum(tilt_cov[:, 0, 0], 0.0)+1e-13
                h11 = np.maximum(tilt_cov[:, 1, 1], 0.0)+1e-13
                h01 = tilt_cov[:, 0, 1]
                det = np.maximum(h00*h11-h01*h01, 1e-26)
                inverse = np.empty_like(tilt_cov)
                inverse[:, 0, 0], inverse[:, 1, 1] = h11/det, h00/det
                inverse[:, 0, 1] = inverse[:, 1, 0] = -h01/det
                covariance -= cross @ inverse @ np.swapaxes(cross, 1, 2)
        covariance = (covariance+np.swapaxes(covariance, 1, 2))*0.5
        # Small FLOAT negative eigenvalues have no proof role; clipping only
        # stabilizes candidate Newton steps, never changes a row bound.
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        covariance = (eigenvectors*np.maximum(eigenvalues, 0.0)[:, None, :]) @ np.swapaxes(eigenvectors, 1, 2)
        normalization = (factorial[width]-factorial[a]-factorial[u-a]-factorial[width-u])/groups
        values[ids] = optimized-normalization-(u*logs[2]+a*logs[3])/groups
        gradients[ids] = means
        gradients[ids, 2] -= u/groups
        gradients[ids, 3] -= a/groups
        hessians[ids] = _PS_LN2*covariance
    return values, gradients, hessians, maximum_residual


def _ps_smooth(values, gradients, hessians, temperature):
    maximum = float(np.max(values))
    weights = np.exp2((values-maximum)/temperature)
    total = float(weights.sum())
    weights /= total
    gradient = weights @ gradients
    centered = gradients-gradient
    hessian = np.einsum('i,ijk->jk', weights, hessians)
    hessian += (_PS_LN2/temperature)*(centered.T @ (weights[:, None]*centered))
    return maximum+temperature*log2(total), gradient, hessian


def _ps_endpoint(logs):
    probability = float(np.exp2(-np.logaddexp2(0.0, -logs[1])))
    first_gradient = np.zeros(6)
    first_gradient[1] = probability
    first_hessian = np.zeros((6, 6))
    first_hessian[1, 1] = _PS_LN2*probability*(1-probability)
    gradients = np.asarray((first_gradient, (0, 0, -1, 0, 1, 1), (0, 0, -1, -1, 1, 1)))
    values = np.asarray((np.logaddexp2(0.0, logs[1]), gradients[1] @ logs, gradients[2] @ logs))
    hessians = np.zeros((3, 6, 6))
    hessians[0] = first_hessian
    return values, gradients, hessians


def partition_averaged_outer_search(profile=PA_PROFILE, death_upper=None, initial=PA_LOGS,
                                    bounds=PS_BOUNDS, max_evaluations=512, max_active_rows=96,
                                    inner_iterations=32, scan_iterations=40):
    """Bounded active-row candidate search; complete scans validate all saved bounds."""
    max_evaluations, max_active_rows = index(max_evaluations), index(max_active_rows)
    if not 32 <= max_evaluations <= 1000 or not 16 <= max_active_rows <= 256:
        raise ValueError('Require32..1000 evaluations and16..256 active rows')
    bounds = np.asarray(bounds, dtype=np.float64)
    initial = np.asarray(initial, dtype=np.float64)
    if (bounds.shape != (6, 2) or not np.all(np.isfinite(bounds)) or np.any(bounds[:, 0] >= bounds[:, 1])
            or np.any(bounds[4:, 0] < 0) or initial.shape != (6,) or not np.all(np.isfinite(initial))):
        raise ValueError('Invalid finite search box or initial logs')
    point = np.clip(initial, bounds[:, 0], bounds[:, 1])
    first = partition_averaged_bound(profile, death_upper, point, max_iterations=scan_iterations)
    if first['status'] == 'empty':
        return first
    h, rank, volume, births, footprint = map(index, profile)
    D0 = first['death_interval'][1]
    G, m = footprint-births, volume-births
    power = (PA_WIDTH//4)*(PA_UPDATES-1)
    linear = np.asarray((-G, -births, 0, 0, -m, -(m-D0)), dtype=np.float64)/power
    route = first['route_log2']/power
    endpoint_weight = PA_WIDTH/power
    active = {(0, 0), (4, 1), (4, 4), (PA_WIDTH, PA_WIDTH//4), (PA_WIDTH, PA_WIDTH)}
    active.update((row['u'], row['a']) for row in first['largest_rows'])
    cache = {}
    best, best_logs = first, point.copy()
    evaluations, scans, accepted_steps, attempted_steps = 0, 1, 0, 0
    maximum_inner_residual = 0.0
    scan_history = [dict(scan=1, raw_log2_bound=first['raw_log2_probability_bound'],
                         active_rows=len(active), row=first['largest_rows'][0])]

    def full_scan(candidate):
        nonlocal best, best_logs, scans
        result = partition_averaged_bound(profile, D0, candidate, max_iterations=scan_iterations)
        scans += 1
        if result['raw_log2_probability_bound'] < best['raw_log2_probability_bound']:
            best, best_logs = result, candidate.copy()
        for row in result['largest_rows']:
            if len(active) >= max_active_rows:
                break
            active.add((row['u'], row['a']))
        scan_history.append(dict(scan=scans, raw_log2_bound=result['raw_log2_probability_bound'],
                                 active_rows=len(active), row=result['largest_rows'][0]))
        return result

    def evaluate(candidate, temperature):
        nonlocal evaluations, maximum_inner_residual
        if evaluations >= max_evaluations:
            return None
        evaluations += 1
        selected_rows = np.asarray(sorted(active), dtype=np.int64)
        values, gradients, hessians, residual = partition_averaged_envelopes(
            candidate, selected_rows, cache=cache, inner_iterations=inner_iterations)
        maximum_inner_residual = max(maximum_inner_residual, residual)
        local = _ps_smooth(values, gradients, hessians, temperature)
        endpoint = _ps_smooth(*_ps_endpoint(candidate), temperature)
        return (local[0]+endpoint_weight*endpoint[0]+linear @ candidate+route,
                local[1]+endpoint_weight*endpoint[1]+linear,
                local[2]+endpoint_weight*endpoint[2])

    for temperature in (1/16, 1/64, 1/256, 1/1024, 1/4096, 1/16384):
        for _ in range(24):
            current = evaluate(point, temperature)
            if current is None:
                break
            value, gradient, hessian = current
            fixed = (((point <= bounds[:, 0]+1e-10) & (gradient >= 0))
                     | ((point >= bounds[:, 1]-1e-10) & (gradient <= 0)))
            free = np.flatnonzero(~fixed)
            if not len(free) or np.max(np.abs(gradient[free])) < 1e-8:
                break
            reduced = hessian[np.ix_(free, free)]
            ridge = max(1e-10, 1e-8*float(np.trace(reduced))/len(free))
            direction = np.zeros(6)
            try:
                direction[free] = np.linalg.solve(reduced+ridge*np.eye(len(free)), -gradient[free])
            except np.linalg.LinAlgError:
                direction[free] = -gradient[free]
            direction *= min(1.0, 1.5/max(1e-300, float(np.max(np.abs(direction)))))
            direction = np.clip(point+direction, bounds[:, 0], bounds[:, 1])-point
            if float(gradient @ direction) >= -1e-15:
                direction = np.clip(point-gradient/max(1.0, float(np.max(np.abs(gradient)))),
                                    bounds[:, 0], bounds[:, 1])-point
            descent = float(gradient @ direction)
            if descent >= -1e-15:
                break
            accepted, alpha = False, 1.0
            for _ in range(16):
                candidate = point+alpha*direction
                trial = evaluate(candidate, temperature)
                if trial is None:
                    break
                if trial[0] <= value+1e-4*alpha*descent:
                    point, accepted = candidate, True
                    accepted_steps += 1
                    break
                alpha *= 0.5
            attempted_steps += 1
            if accepted and accepted_steps % 4 == 0:
                full_scan(point)
            if not accepted or evaluations >= max_evaluations or scans >= 32:
                break
        if evaluations >= max_evaluations or scans >= 32:
            break
        full_scan(point)
    full_scan(point)
    result = dict(best)
    result.update(search_box=bounds.tolist(), optimizer_evaluations=evaluations,
                  attempted_steps=attempted_steps, accepted_steps=accepted_steps,
                  complete_physical_row_scans=scans, final_active_rows=len(active),
                  maximum_inner_mean_residual=maximum_inner_residual,
                  scan_history=scan_history, returned_logs=best_logs.tolist(),
                  every_reported_bound_has_full_row_scan=True,
                  active_objective_only_is_not_a_bound=True,
                  optimizer_optimality_proved=False, scalar_family_impossibility_proved=False)
    return result


def partition_averaged_search_checks():
    coefficient_checks = derivative_checks = 0
    settings = ((Fraction(1, 2), Fraction(1, 8), Fraction(2), Fraction(1, 2), Fraction(1), Fraction(2)),
                (Fraction(3, 4), Fraction(1, 4), Fraction(3, 2), Fraction(3, 4), Fraction(1), Fraction(3, 2)))
    for markers in settings:
        z, w, A, B, b, k = markers
        combined = {(i, j): Fraction(0) for i in range(5) for j in range(5-i)}
        for exponent, coefficient, (_, i, j, ell) in zip(_PS_OUTER, _PS_MULT, _PS_IDENTITIES):
            value = Fraction(int(coefficient))
            for marker, power in zip(markers, exponent):
                value *= marker**int(power)
            combined[i, j] += value
        beta, d = b*k, 1/k
        for (i, j), value in combined.items():
            r = i+j
            births = sum(Fraction(comb(4-r, ell))*(w*B)**ell for ell in range(int(r == 0), 5-r))
            expected = comb(4, i)*comb(4-i, j)*(d**i+z**4*A**4*(beta*B)**r*births)
            assert value == expected
            coefficient_checks += 1
    point = np.asarray(PA_LOGS)
    rows = np.asarray(((0, 0), (4, 1), (4, 2), (4, 4), (8, 2), (8, 8)))
    values, gradients, hessians, residual = partition_averaged_envelopes(point, rows, width=8, inner_iterations=64)
    assert residual < 1e-6
    full = partition_averaged_rows(point[4]+point[5], -point[5], point[0], point[1], point[2], point[3],
                                    width=8, physical_rows=True, max_iterations=64)
    lookup = {tuple(map(int, row)): i for i, row in enumerate(full['rows'])}
    for i, row in enumerate(rows):
        assert abs(2*values[i]-full['row_log2_upper'][lookup[tuple(row)]]) < 1e-8
        assert np.linalg.eigvalsh(hessians[i]).min() > -1e-8
    for coordinate in range(6):
        delta = np.zeros(6)
        delta[coordinate] = 2e-4
        plus = partition_averaged_envelopes(point+delta, rows, width=8, inner_iterations=64)
        minus = partition_averaged_envelopes(point-delta, rows, width=8, inner_iterations=64)
        assert np.max(np.abs((plus[0]-minus[0])/(4e-4)-gradients[:, coordinate])) < 2e-5
        assert np.max(np.abs((plus[1]-minus[1])/(4e-4)-hessians[:, :, coordinate])) < 2e-4
        derivative_checks += 1
    print(dict(positive_monomials=len(_PS_MULT), exact_collapsed_coefficients=coefficient_checks,
               face_and_row_envelopes=len(rows), envelope_gradient_hessian_coordinates=derivative_checks,
               full_code_certificate=False), flush=True)


def partition_averaged_search_main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--checks', action='store_true')
    modes.add_argument('--search', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=PA_PROFILE)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--logs', type=float, nargs=6, default=PA_LOGS)
    parser.add_argument('--evaluations', type=int, default=512)
    parser.add_argument('--active-rows', type=int, default=96)
    parser.add_argument('--inner-iterations', type=int, default=32)
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
            partition_averaged_search_checks()
        else:
            print(partition_averaged_outer_search(args.profile, args.death_upper, args.logs,
                                                 max_evaluations=args.evaluations,
                                                 max_active_rows=args.active_rows,
                                                 inner_iterations=args.inner_iterations), flush=True)
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    partition_averaged_search_main()
