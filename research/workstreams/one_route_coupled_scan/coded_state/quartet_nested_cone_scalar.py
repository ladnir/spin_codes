"""Seven-variable FLOAT scalar control on reachable nested selector states.

Every selected quartet contains at least one allowed position: an old-live
intersection, or a mandatory birth if that intersection is empty.  Thus
all intermediate states(U,Aprev) satisfy |Aprev|>=|U|/4.  The formal full
endpoint is not an intermediate selected state.

For0<t<=1, multiply each local normalized action by t^(r/4-a).
If every adjusted local action is<=lambda, their product is bounded by
lambda^(n/4)*V*t^(|Aprev|-|U|/4)<=lambda^(n/4)*V on the reachable cone.
The free endpoint and all external profile corrections stay unchanged.

This helper wraps the frozen six-variable local actions and profile
accounting.  It adds only the cone variable and a dimension-generic bounded
Newton search.  Exact small tests use t=v^4, avoiding fractional powers.
No production probability or optimizer certificate is claimed.  Parent
alone runs the serial checks/probes; no frozen source is modified.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np


NC_WIDTH = 512
NC_UPDATES = 4096
NC_POWER = (NC_WIDTH // 4) * (NC_UPDATES - 1)
NC_PROFILE = (4096, 819200, 838861, 65536, 1081344)
NC_INITIAL = (-1.25, -4.0, 1.0, -0.5, 0.5, 0.5, -1.0)
NC_BOUNDS = ((-8.0, 4.0), (-24.0, 8.0), (-8.0, 8.0),
             (-12.0, 8.0), (0.0, 12.0), (0.0, 12.0), (-16.0, 0.0))
NC_LOG_NAMES = ('log_z', 'log_w', 'log_A', 'log_B', 'log_b', 'log_k', 'log_t')


def _nc_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_nested_selector_scalar
    return getattr(quartet_nested_selector_scalar, name)


def _nc_actions(logs):
    """Lift the old affine actions; the endpoint has zero cone derivative."""
    rows, endpoints = _nc_dependency('_ns_actions')(logs[:6])
    cases = _nc_dependency('_NS_CASES')
    adjusted, lifted_endpoints = [], []
    for (value, gradient, hessian), (a, r, _, _) in zip(rows, cases):
        slope = r / 4 - a
        new_gradient = np.empty(7)
        new_gradient[:6], new_gradient[6] = gradient, slope
        new_hessian = np.zeros((7, 7))
        new_hessian[:6, :6] = hessian
        adjusted.append((value + slope * logs[6], new_gradient, new_hessian))
    for value, gradient, hessian in endpoints:
        new_gradient = np.zeros(7)
        new_gradient[:6] = gradient
        new_hessian = np.zeros((7, 7))
        new_hessian[:6, :6] = hessian
        lifted_endpoints.append((value, new_gradient, new_hessian))
    return adjusted, lifted_endpoints


def _nc_profile(profile, death_upper):
    data = dict(_nc_dependency('nested_selector_profile')(profile, death_upper))
    data['linear'] = np.append(data['linear'], 0.0)
    return data


def _nc_objective(logs, temperature, data):
    rows, endpoints = _nc_actions(logs)
    rvalues = np.asarray([row[0] for row in rows])
    evalues = np.asarray([row[0] for row in endpoints])
    smooth_max = _nc_dependency('_ns_smooth_max')
    lam = smooth_max(rvalues, np.asarray([row[1] for row in rows]),
                     np.asarray([row[2] for row in rows]), temperature)
    end = smooth_max(evalues, np.asarray([row[1] for row in endpoints]),
                     np.asarray([row[2] for row in endpoints]), temperature)
    endpoint_weight = NC_WIDTH / NC_POWER
    correction = float(data['linear'] @ logs) + data['route_log2'] / NC_POWER
    value = lam[0] + endpoint_weight * end[0] + correction
    gradient = lam[1] + endpoint_weight * end[1] + data['linear']
    hessian = lam[2] + endpoint_weight * end[2]
    raw = float(rvalues.max()) + endpoint_weight * float(evalues.max()) + correction
    return value, gradient, hessian, raw


def _nc_validate_logs(logs):
    logs = np.asarray(logs, dtype=np.float64)
    if (logs.shape != (7,) or not np.all(np.isfinite(logs))
            or np.any(logs[4:6] < 0) or logs[6] > 0):
        raise ValueError('Seven finite logs required: log_b,log_k>=0 and log_t<=0')
    return logs


def nested_cone_bound(profile=NC_PROFILE, death_upper=None, logs=NC_INITIAL):
    logs = _nc_validate_logs(logs)
    data = _nc_profile(profile, death_upper)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    rows, endpoints = _nc_actions(logs)
    values = [row[0] for row in rows]
    log_lam, log_endpoint = max(values), max(row[0] for row in endpoints)
    raw = (NC_POWER * log_lam + NC_WIDTH * log_endpoint + data['route_log2']
           + NC_POWER * float(data['linear'] @ logs))
    if not np.isfinite(raw):
        raise ArithmeticError('Nonfinite cone scalar bound')
    cases = _nc_dependency('_NS_CASES')
    order = sorted(range(len(values)), key=values.__getitem__, reverse=True)
    return dict(profile=data['profile'], status='low_actual_death_class',
                death_interval=(data['death_lower'], data['death_upper']),
                logs=dict(zip(NC_LOG_NAMES, map(float, logs))),
                log_beta=float(logs[4] + logs[5]), log_d=float(-logs[5]),
                log2_lambda=log_lam, log2_endpoint_per_coordinate=log_endpoint,
                route_log2=data['route_log2'], raw_log2_probability_bound=raw,
                tightest_local_cases=[dict(a=cases[i][0], r=cases[i][1],
                                          log2_cone_adjusted_ratio=values[i])
                                     for i in order[:4]],
                reachable_cone='4*|Aprev|>=|U|', cone_marker_has_no_external_cost=True,
                births_integrated=True, separate_birth_mask_factor=False,
                H_union_already_in_route=True, free_endpoint=True,
                floating_diagnostic_only=True, global_obstruction_proved=False,
                full_code_certificate=False)


def nested_cone_search(profile=NC_PROFILE, death_upper=None, initial=NC_INITIAL,
                       bounds=NC_BOUNDS, max_evaluations=1024):
    """Bounded convex Newton scout; all production bounds remain FLOAT."""
    data = _nc_profile(profile, death_upper)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    maximum = index(max_evaluations)
    bounds = np.asarray(bounds, dtype=np.float64)
    if (bounds.shape != (7, 2) or not np.all(np.isfinite(bounds))
            or np.any(bounds[:, 0] >= bounds[:, 1])
            or np.any(bounds[4:6, 0] < 0) or bounds[6, 1] > 0
            or not 32 <= maximum <= 4096):
        raise ValueError('Finite seven-variable box and32..4096 evaluations required')
    initial = _nc_validate_logs(initial)
    point = np.clip(initial, bounds[:, 0], bounds[:, 1])
    dimension = len(point)
    evaluations = iterations = 0
    best_point, best_value = point.copy(), float('inf')
    last_temperature = 0.25

    def evaluate(candidate, temperature):
        nonlocal evaluations, best_point, best_value
        if evaluations >= maximum:
            return None
        evaluations += 1
        result = _nc_objective(candidate, temperature, data)
        if not all(np.all(np.isfinite(value)) for value in result):
            raise ArithmeticError('Nonfinite convex objective or derivative')
        if result[3] < best_value:
            best_value, best_point = result[3], candidate.copy()
        return result

    for temperature in (1/4, 1/16, 1/64, 1/256, 1/1024, 1/4096, 1/16384, 1/65536):
        last_temperature = temperature
        for _ in range(28):
            current = evaluate(point, temperature)
            if current is None:
                break
            value, gradient, hessian, _ = current
            active = (((point <= bounds[:, 0] + 1e-11) & (gradient >= 0))
                      | ((point >= bounds[:, 1] - 1e-11) & (gradient <= 0)))
            free = np.flatnonzero(~active)
            if not len(free) or np.max(np.abs(gradient[free])) < 1e-9:
                break
            reduced = hessian[np.ix_(free, free)]
            ridge = max(1e-10, 1e-8 * float(np.trace(reduced)) / len(free))
            step = np.zeros(dimension)
            try:
                step[free] = np.linalg.solve(reduced + ridge * np.eye(len(free)), -gradient[free])
            except np.linalg.LinAlgError:
                step[free] = -gradient[free]
            length = float(np.max(np.abs(step)))
            if length > 2:
                step *= 2 / length
            direction = np.clip(point + step, bounds[:, 0], bounds[:, 1]) - point
            if float(gradient @ direction) >= -1e-15:
                scale = max(1.0, float(np.max(np.abs(gradient))))
                direction = np.clip(point - gradient / scale, bounds[:, 0], bounds[:, 1]) - point
            descent = float(gradient @ direction)
            if descent >= -1e-15:
                break
            accepted, alpha = False, 1.0
            for _ in range(24):
                candidate = point + alpha * direction
                trial = evaluate(candidate, temperature)
                if trial is None:
                    break
                if trial[0] <= value + 1e-4 * alpha * descent:
                    point, accepted = candidate, True
                    break
                alpha *= 0.5
            iterations += 1
            if not accepted:
                break
        if evaluations >= maximum:
            break
    smooth, gradient, _, _ = _nc_objective(best_point, last_temperature, data)
    corner = np.where(gradient >= 0, bounds[:, 0], bounds[:, 1])
    smoothing_error = last_temperature * (log2(9) + (NC_WIDTH / NC_POWER) * log2(3))
    lower_box = NC_POWER * (smooth + float(gradient @ (corner - best_point)) - smoothing_error)
    result = nested_cone_bound(profile, death_upper, best_point)
    result.update(optimizer_evaluations=evaluations, optimizer_iterations=iterations,
                  search_box=bounds.tolist(), final_temperature=last_temperature,
                  finite_box_lower_diagnostic=lower_box,
                  finite_box_gap_diagnostic=result['raw_log2_probability_bound'] - lower_box,
                  smoothing_error_bits=NC_POWER * smoothing_error,
                  touches_box_boundary=bool(np.any(np.isclose(best_point, bounds[:, 0], atol=1e-7))
                                             or np.any(np.isclose(best_point, bounds[:, 1], atol=1e-7))),
                  lower_diagnostic_is_exact_certificate=False,
                  unbounded_parameter_domain_certified=False)
    return result


def nested_cone_checks():
    local = _nc_dependency('_ns_local')
    transitions = _nc_dependency('_ns_transitions')
    local_checks = physical_checks = reachable_checks = endpoint_checks = moment_checks = 0
    settings = ((Fraction(2), Fraction(1, 2), Fraction(1, 2), Fraction(1, 8), Fraction(2), Fraction(1, 2)),
                (Fraction(3, 2), Fraction(2, 3), Fraction(3, 4), Fraction(1, 4), Fraction(3, 2), Fraction(3, 4)))
    for beta, d, z, w, A, B in settings:
        for v in (Fraction(1), Fraction(3, 4)):
            ratios = {(a, r): local(a, r, beta, d, z, w, A, B) / (A ** r * B ** a)
                                  * v ** (r - 4 * a)
                      for r in range(5) for a in range(r + 1)}
            lam = max(ratios.values())
            assert lam == max(ratios[(a, r)] for r in range(5) for a in (0, r))
            for (a, r), ratio in ratios.items():
                assert ratio <= lam
                assert local(a, r, beta, d, z, w, A, B) <= lam * A ** r * B ** a * v ** (4 * a - r)
                local_checks += 1
            endpoint = max(1 + w, (beta / A) * max(Fraction(1), 1 / B))
            for width in (4, 8):
                partitions = [(15,)] if width == 4 else [
                    (first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                    for first in (1 + sum(1 << i for i in rest),)]
                for population in range(0, width + 1, 4):
                    old = (1 << population) - 1
                    for previous_size in range(population // 4, population + 1):
                        previous = (1 << previous_size) - 1
                        potential = A ** population * B ** previous_size
                        assert beta ** population * (1 + w) ** (width - population) <= endpoint ** width * potential
                        endpoint_checks += 1
                        for partition in partitions:
                            exact = Fraction(0)
                            for new, allowed, weight in transitions(old, previous, partition, beta, d, z, w):
                                assert allowed & ~new == 0
                                assert 4 * allowed.bit_count() >= new.bit_count()
                                exact += weight * A ** new.bit_count() * B ** allowed.bit_count()
                                reachable_checks += 1
                            assert exact <= lam ** (width // 4) * potential * v ** (4 * previous_size - population)
                            assert exact <= lam ** (width // 4) * potential
                            physical_checks += 1
            states = {(0, 0): Fraction(1)}
            for updates in range(1, 5):
                exact = sum(weight * beta ** old.bit_count() * (1 + w) ** (4 - old.bit_count())
                            for (old, previous), weight in states.items())
                assert exact <= lam ** (updates - 1) * endpoint ** 4
                moment_checks += 1
                following = {}
                for (old, previous), weight in states.items():
                    for new, allowed, step in transitions(old, previous, (15,), beta, d, z, w):
                        following[(new, allowed)] = following.get((new, allowed), Fraction(0)) + weight * step
                states = following
    data = _nc_profile(NC_PROFILE, None)
    point = np.asarray(NC_INITIAL)
    value, gradient, hessian, raw = _nc_objective(point, 1/16, data)
    smoothing_error = (1/16) * (log2(9) + NC_WIDTH / NC_POWER * log2(3))
    assert raw <= value + 1e-12 and value <= raw + smoothing_error + 1e-12
    assert np.linalg.eigvalsh(hessian).min() > -1e-9
    for coordinate in range(7):
        delta = np.zeros(7)
        delta[coordinate] = 1e-5
        plus, minus = _nc_objective(point + delta, 1/16, data), _nc_objective(point - delta, 1/16, data)
        assert abs((plus[0] - minus[0]) / 2e-5 - gradient[coordinate]) < 1e-6
        assert np.max(np.abs((plus[1] - minus[1]) / 2e-5 - hessian[:, coordinate])) < 1e-5
    neutral = point.copy()
    neutral[6] = 0
    cone = nested_cone_bound(logs=neutral)
    old = _nc_dependency('nested_selector_bound')(logs=neutral[:6])
    assert abs(cone['raw_log2_probability_bound'] - old['raw_log2_probability_bound']) < 1e-8
    report = dict(exact_local_cone_comparisons=local_checks,
                  exact_reachable_transition_checks=reachable_checks,
                  exact_cone_partition_actions=physical_checks,
                  exact_free_endpoints=endpoint_checks, exact_temporal_moments=moment_checks,
                  float_gradient_hessian_coordinates=7, neutral_cone_matches_frozen_control=True,
                  exact_test_markers_are_fourth_powers=True, full_code_certificate=False)
    print(report, flush=True)
    return report


def nested_cone_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    mode.add_argument('--search', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=NC_PROFILE)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--logs', type=float, nargs=7, default=NC_INITIAL)
    parser.add_argument('--evaluations', type=int, default=1024)
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
            nested_cone_checks()
        elif args.search:
            print(nested_cone_search(args.profile, args.death_upper, args.logs,
                                    max_evaluations=args.evaluations), flush=True)
        else:
            print(nested_cone_bound(args.profile, args.death_upper, args.logs), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    nested_cone_main()
