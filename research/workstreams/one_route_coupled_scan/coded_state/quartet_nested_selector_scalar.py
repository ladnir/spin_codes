"""Bounded FLOAT scalar control for nested live-selector pairs and low deaths.

State(U_t,A_(t-1)) is nested. A new quartet contains r old live positions
and a<=r old allowed positions. Its positive weighted action on the
potential A^|U| B^|Aprev| is

 h(a,r)=d^a+z^4*A^4*(beta*B)^r*((1+w*B)^(4-r)-1[r=0]).

The subtraction removes only an empty birth choice; code evaluates its
positive polynomial. If h(a,r)<=lambda*A^r*B^a for all15 cases, the
integrated selector/birth moment is at most

 lambda^(128*(4096-1))*max(1+w,(beta/A)*max(1,B^-1))^512.

At fixed actual profile m=L-f and D<=D0, necessary conditions are
A0>=m and A0-C_A>=m-D0. Set beta=b*k,d=1/k with b,k>=1.
The exponent adds the capped route factor and subtracts
G*logz+f*logw+m*logb+(m-D0)*logk. No separate birth-mask count,
H union, or old empty-quartet marker correction is added.

For fixed r the normalized action is convex in a, so a=0,r suffice
(nine unique cases). Its log is a positive log-sum-exp of affine forms;
the whole objective is convex in(logz,logw,logA,logB,logb,logk).
The optional bounded search uses smoothed maxima and damped Newton steps
in a stated finite box. It is a FLOAT diagnostic, not an optimizer proof
or a probability certificate. No large state DP or SciPy is used.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, lgamma, log, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np


NS_WIDTH = 512
NS_UPDATES = 4096
NS_N = NS_WIDTH * NS_UPDATES
NS_M = 8192
NS_POWER = (NS_WIDTH // 4) * (NS_UPDATES - 1)
NS_W = 943718
NS_PROFILE = (4096, 819200, 838861, 65536, 1081344)
NS_INITIAL = (-1.25, -4.0, 1.0, -0.5, 0.5, 0.5)
NS_BOUNDS = ((-8.0, 4.0), (-24.0, 8.0), (-8.0, 8.0),
             (-12.0, 8.0), (0.0, 12.0), (0.0, 12.0))
NS_LOG_NAMES = ('log_z', 'log_w', 'log_A', 'log_B', 'log_b', 'log_k')
_NS_LN2 = log(2)


def _ns_model_cases():
    cases = []
    for r in range(5):
        for a in ((0,) if r == 0 else (0, r)):
            exponents = [(0, 0, -r, -a, 0, -a)]
            coefficients = [0.0]
            for ell in range(int(r == 0), 5 - r):
                exponents.append((4, ell, 4 - r, r + ell - a, r, r))
                coefficients.append(log2(comb(4 - r, ell)))
            cases.append((a, r, np.asarray(exponents, dtype=np.float64),
                          np.asarray(coefficients, dtype=np.float64)))
    return tuple(cases)


_NS_CASES = _ns_model_cases()


def _ns_logsum(affine, exponents):
    maximum = float(np.max(affine))
    weights = np.exp2(affine - maximum)
    total = float(weights.sum())
    weights /= total
    gradient = weights @ exponents
    centered = exponents - gradient
    hessian = _NS_LN2 * (centered.T @ (weights[:, None] * centered))
    return maximum + log2(total), gradient, hessian


def _ns_smooth_max(values, gradients, hessians, temperature):
    maximum = float(np.max(values))
    weights = np.exp2((values - maximum) / temperature)
    total = float(weights.sum())
    weights /= total
    gradient = weights @ gradients
    centered = gradients - gradient
    hessian = np.einsum('i,ijk->jk', weights, hessians)
    hessian += (_NS_LN2 / temperature) * (centered.T @ (weights[:, None] * centered))
    return maximum + temperature * log2(total), gradient, hessian


def _ns_actions(logs):
    rows = [_ns_logsum(exponents @ logs + coefficients, exponents)
            for _, _, exponents, coefficients in _NS_CASES]
    endpoint_exponents = np.asarray(((0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0)), dtype=np.float64)
    first = _ns_logsum(endpoint_exponents @ logs, endpoint_exponents)
    linear1 = np.asarray((0, 0, -1, 0, 1, 1), dtype=np.float64)
    linear2 = np.asarray((0, 0, -1, -1, 1, 1), dtype=np.float64)
    endpoints = (first, (float(linear1 @ logs), linear1, np.zeros((6, 6))),
                 (float(linear2 @ logs), linear2, np.zeros((6, 6))))
    return rows, endpoints


def _ns_choose(n, k):
    if not 0 <= k <= n:
        return float('-inf')
    if k in (0, n):
        return 0.0
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / _NS_LN2


def nested_selector_profile(profile=NS_PROFILE, death_upper=None):
    h, rank, volume, births, footprint = map(index, profile)
    if not (1 <= h <= NS_M and 128 * h <= rank <= (1 << 20)
            and 1 <= births <= volume <= footprint <= NS_N):
        raise ValueError('Invalid actual-state profile')
    G, m = footprint - births, volume - births
    if G % 4 or G > NS_N - NS_WIDTH:
        raise ValueError('Invalid actual live-selector volume')
    CH = footprint - rank + 128 * h
    lower, upper = max(0, births - NS_WIDTH), min(footprint - volume, NS_W // 2)
    D0 = NS_W - volume if death_upper is None else index(death_upper)
    D0 = min(D0, upper)
    empty = D0 < lower or not (129 * h <= CH <= 256 * h
                              and 0 <= CH - births <= G
                              and 0 <= 256 * h - CH <= NS_N - footprint)
    route = float('-inf') if empty else min(0.0,
        _ns_choose(G, CH - births) + _ns_choose(NS_N - footprint, 256 * h - CH)
        - 255 * _ns_choose(NS_M, h))
    linear = np.asarray((-G, -births, 0, 0, -m, -(m - D0)), dtype=np.float64) / NS_POWER
    return dict(profile=(h, rank, volume, births, footprint), G=G, m=m,
                C_H=CH, death_lower=lower, death_upper=D0, empty=empty,
                route_log2=route, linear=linear)


def _ns_objective(logs, temperature, profile_data):
    rows, endpoints = _ns_actions(logs)
    rvalues = np.asarray([row[0] for row in rows])
    evalues = np.asarray([row[0] for row in endpoints])
    lam = _ns_smooth_max(rvalues, np.asarray([row[1] for row in rows]),
                         np.asarray([row[2] for row in rows]), temperature)
    end = _ns_smooth_max(evalues, np.asarray([row[1] for row in endpoints]),
                         np.asarray([row[2] for row in endpoints]), temperature)
    endpoint_weight = NS_WIDTH / NS_POWER
    correction = float(profile_data['linear'] @ logs) + profile_data['route_log2'] / NS_POWER
    value = lam[0] + endpoint_weight * end[0] + correction
    gradient = lam[1] + endpoint_weight * end[1] + profile_data['linear']
    hessian = lam[2] + endpoint_weight * end[2]
    raw = float(rvalues.max()) + endpoint_weight * float(evalues.max()) + correction
    return value, gradient, hessian, raw


def nested_selector_bound(profile=NS_PROFILE, death_upper=None, logs=NS_INITIAL):
    """One unsmoothed nine-case FLOAT bound, with all15 exact cases audited."""
    logs = np.asarray(logs, dtype=np.float64)
    if logs.shape != (6,) or not np.all(np.isfinite(logs)) or np.any(logs[4:] < 0):
        raise ValueError('Six finite logs required; log_b and log_k must be nonnegative')
    data = nested_selector_profile(profile, death_upper)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    rows, endpoints = _ns_actions(logs)
    values = [row[0] for row in rows]
    log_lam, log_endpoint = max(values), max(row[0] for row in endpoints)
    raw = (NS_POWER * log_lam + NS_WIDTH * log_endpoint + data['route_log2']
           + NS_POWER * float(data['linear'] @ logs))
    order = sorted(range(len(values)), key=values.__getitem__, reverse=True)
    return dict(profile=data['profile'], status='low_actual_death_class',
                death_interval=(data['death_lower'], data['death_upper']),
                logs=dict(zip(NS_LOG_NAMES, map(float, logs))),
                log_beta=float(logs[4] + logs[5]), log_d=float(-logs[5]),
                log2_lambda=log_lam, log2_endpoint_per_coordinate=log_endpoint,
                route_log2=data['route_log2'], raw_log2_probability_bound=raw,
                tightest_local_cases=[dict(a=_NS_CASES[i][0], r=_NS_CASES[i][1],
                                           log2_ratio=values[i]) for i in order[:4]],
                births_integrated=True, separate_birth_mask_factor=False,
                H_union_already_in_route=True, free_endpoint=True,
                floating_diagnostic_only=True, global_obstruction_proved=False,
                full_code_certificate=False)


def nested_selector_search(profile=NS_PROFILE, death_upper=None, initial=NS_INITIAL,
                            bounds=NS_BOUNDS, max_evaluations=1024):
    """Bounded smoothed-convex Newton scout in six log parameters; no SciPy."""
    data = nested_selector_profile(profile, death_upper)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    max_evaluations = index(max_evaluations)
    bounds = np.asarray(bounds, dtype=np.float64)
    if bounds.shape != (6, 2) or np.any(bounds[:, 0] >= bounds[:, 1]) or np.any(bounds[4:, 0] < 0):
        raise ValueError('Invalid finite convex-search box')
    if not np.all(np.isfinite(bounds)) or not 32 <= max_evaluations <= 4096:
        raise ValueError('Require finite bounds and32..4096 evaluations')
    point = np.clip(np.asarray(initial, dtype=np.float64), bounds[:, 0], bounds[:, 1])
    if point.shape != (6,) or not np.all(np.isfinite(point)):
        raise ValueError('Invalid initial logs')
    evaluations, iterations = 0, 0
    best_point, best_value = point.copy(), float('inf')
    last_temperature = 0.25
    def evaluate(candidate, temperature):
        nonlocal evaluations, best_point, best_value
        if evaluations >= max_evaluations:
            return None
        evaluations += 1
        result = _ns_objective(candidate, temperature, data)
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
            active = ((point <= bounds[:, 0] + 1e-11) & (gradient >= 0)
                      | (point >= bounds[:, 1] - 1e-11) & (gradient <= 0))
            free = np.flatnonzero(~active)
            if not len(free) or np.max(np.abs(gradient[free])) < 1e-9:
                break
            reduced = hessian[np.ix_(free, free)]
            ridge = max(1e-10, 1e-8 * float(np.trace(reduced)) / len(free))
            step = np.zeros(6)
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
            accepted = False
            alpha = 1.0
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
        if evaluations >= max_evaluations:
            break
    # First-order support of the smooth objective, minus its uniform smoothing
    # error, gives a FLOAT lower diagnostic only on this stated finite box.
    smooth, gradient, _, _ = _ns_objective(best_point, last_temperature, data)
    minimizing_corner = np.where(gradient >= 0, bounds[:, 0], bounds[:, 1])
    smoothing_error = last_temperature * (log2(9) + (NS_WIDTH / NS_POWER) * log2(3))
    lower_box = NS_POWER * (smooth + float(gradient @ (minimizing_corner - best_point)) - smoothing_error)
    result = nested_selector_bound(profile, death_upper, best_point)
    result.update(optimizer_evaluations=evaluations, optimizer_iterations=iterations,
                  search_box=bounds.tolist(), final_temperature=last_temperature,
                  finite_box_lower_diagnostic=lower_box,
                  finite_box_gap_diagnostic=result['raw_log2_probability_bound'] - lower_box,
                  smoothing_error_bits=NS_POWER * smoothing_error,
                  touches_box_boundary=bool(np.any(np.isclose(best_point, bounds[:, 0], atol=1e-7))
                                             or np.any(np.isclose(best_point, bounds[:, 1], atol=1e-7))),
                  lower_diagnostic_is_exact_certificate=False,
                  unbounded_parameter_domain_certified=False)
    return result


def _ns_subsets(mask):
    current = mask
    while True:
        yield current
        if current == 0:
            break
        current = (current - 1) & mask


def _ns_local(a, r, beta, d, z, w, A, B):
    selected = sum(comb(4 - r, ell) * (w * B) ** ell
                   for ell in range(int(r == 0), 5 - r))
    return d ** a + z ** 4 * A ** 4 * (beta * B) ** r * selected


def _ns_transitions(old, previous_allowed, partition, beta, d, z, w):
    for selector in range(1 << len(partition)):
        selected = [block for j, block in enumerate(partition) if selector >> j & 1]
        new = sum(selected)
        for births in _ns_subsets(new & ~old):
            if any(not (block & old) and not (block & births) for block in selected):
                continue
            allowed = (old & new) | births
            overlap = (old & new).bit_count()
            closures = (previous_allowed & ~allowed).bit_count()
            weight = z ** new.bit_count() * beta ** overlap * d ** closures * w ** births.bit_count()
            yield new, allowed, weight


def nested_selector_checks():
    local_checks = physical_checks = endpoint_checks = moment_checks = 0
    settings = ((Fraction(2), Fraction(1, 2), Fraction(1, 2), Fraction(1, 8), Fraction(2), Fraction(1, 2)),
                (Fraction(3, 2), Fraction(2, 3), Fraction(3, 4), Fraction(1, 4), Fraction(3, 2), Fraction(3, 4)),
                (Fraction(4), Fraction(1, 4), Fraction(1, 4), Fraction(1, 16), Fraction(2), Fraction(1)))
    for beta, d, z, w, A, B in settings:
        ratios = {}
        for r in range(5):
            old = (1 << r) - 1
            for a in range(r + 1):
                previous = (1 << a) - 1
                exact = sum(weight * A ** new.bit_count() * B ** allowed.bit_count()
                            for new, allowed, weight in _ns_transitions(old, previous, (15,), beta, d, z, w))
                assert exact == _ns_local(a, r, beta, d, z, w, A, B)
                ratios[(a, r)] = exact / (A ** r * B ** a)
                local_checks += 1
        lam = max(ratios.values())
        assert lam == max(ratios[(a, r)] for r in range(5) for a in (0, r))
        end = max(1 + w, (beta / A) * max(Fraction(1), 1 / B))
        for width in (4, 8):
            partitions = [(15,)] if width == 4 else [
                (first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                for first in (1 + sum(1 << i for i in rest),)]
            for population in range(0, width + 1, 4):
                old = (1 << population) - 1
                for previous_size in range(population + 1):
                    previous = (1 << previous_size) - 1
                    terminal = beta ** population * (1 + w) ** (width - population)
                    assert terminal <= end ** width * A ** population * B ** previous_size
                    endpoint_checks += 1
                    for partition in partitions:
                        exact = sum(weight * A ** new.bit_count() * B ** allowed.bit_count()
                                    for new, allowed, weight in _ns_transitions(old, previous, partition, beta, d, z, w))
                        local_product = Fraction(1)
                        for block in partition:
                            local_product *= _ns_local((previous & block).bit_count(), (old & block).bit_count(),
                                                       beta, d, z, w, A, B)
                        assert exact == local_product <= lam ** (width // 4) * A ** population * B ** previous_size
                        physical_checks += 1
        states = {(0, 0): Fraction(1)}
        for updates in range(1, 5):
            exact = sum(weight * beta ** old.bit_count() * (1 + w) ** (4 - old.bit_count())
                        for (old, previous), weight in states.items())
            assert exact <= lam ** (updates - 1) * end ** 4
            if updates == 1:
                assert exact == (1 + w) ** 4
            moment_checks += 1
            following = {}
            for (old, previous), weight in states.items():
                for new, allowed, step in _ns_transitions(old, previous, (15,), beta, d, z, w):
                    following[(new, allowed)] = following.get((new, allowed), Fraction(0)) + weight * step
            states = following
    data = nested_selector_profile()
    point = np.asarray(NS_INITIAL)
    value, gradient, hessian, raw = _ns_objective(point, 1/16, data)
    epsilon = (1/16) * (log2(9) + NS_WIDTH / NS_POWER * log2(3))
    assert raw <= value + 1e-12 and value <= raw + epsilon + 1e-12
    assert np.linalg.eigvalsh(hessian).min() > -1e-9
    for coordinate in range(6):
        delta = np.zeros(6)
        delta[coordinate] = 1e-5
        plus = _ns_objective(point + delta, 1/16, data)
        minus = _ns_objective(point - delta, 1/16, data)
        assert abs((plus[0] - minus[0]) / (2e-5) - gradient[coordinate]) < 1e-6
        assert np.max(np.abs((plus[1] - minus[1]) / (2e-5) - hessian[:, coordinate])) < 1e-5
    print(dict(exact_local_cases=local_checks, exact_nested_partition_cases=physical_checks,
               exact_free_endpoint_cases=endpoint_checks, exact_temporal_moments=moment_checks,
               float_gradient_hessian_coordinates=6, nine_cases_equal_fifteen=True,
               births_integrated=True, full_code_certificate=False), flush=True)


def nested_selector_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    mode.add_argument('--search', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=NS_PROFILE)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--logs', type=float, nargs=6, default=NS_INITIAL)
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
            nested_selector_checks()
        elif args.search:
            print(nested_selector_search(args.profile, args.death_upper, args.logs,
                                         max_evaluations=args.evaluations), flush=True)
        else:
            print(nested_selector_bound(args.profile, args.death_upper, args.logs), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    nested_selector_main()
