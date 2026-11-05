"""Bounded FLOAT epigraph control for the frozen forced-run scalar family.

The variables are the seven/eight marker logs x and epigraphs l,e.  All35
local constraints f_i(x)<l, all4 endpoint constraints g_j(x)<e, and the
strict marker box are retained.  The normalized objective is

 l + (512/FS_POWER)*e + linear.x + route_log2/FS_POWER.

The barrier uses natural logarithms of positive slacks; the objective and
the scalar-helper actions use base-two logarithms.  This convention gives
the displayed derivatives without an extra ln(2) in the barrier terms.
Newton steps use diagonal scaling and a feasibility/Armijo line search.

For any strictly feasible iterate, positive reciprocal-slack weights are
normalized to sum1 on local constraints and512/FS_POWER on endpoints.
These satisfy the epigraph dual equalities in real arithmetic.  Their
weighted convex action has a tangent whose minimum over the marker box
is a dual lower bound.  Thus the reported finite-box gap need not assume
centrality.  All arithmetic here is FLOAT: neither the lower diagnostic
nor a small gap is an exact certificate or an unbounded-domain result.

The existing helper owns the model, profile corrections, compulsory-hit
option, and probability-bound metadata.  This file changes no model.

Parent verification (2026-10-31):38 barrier derivative coordinates,
20 finite-box dual candidate comparisons,8 strict epigraph checks,and
1 independent quadratic check passed.  The combined compulsory-hit
control at the default profile and D0=104857 used259 evaluations,
117 Newton iterations,and14 centered stages.  Its primal was
177152.546568 bits; the best finite-box dual diagnostic was
177140.016189 bits (gap12.530380).  The best dual occurred at mu=6.4e-7,
with free-coordinate residual2.83e-6 and only the log_b lower face near
active.  Later smaller-mu cancellation worsened the dual; retaining the
best earlier dual was intentional.  This bounded FLOAT result supports
moving beyond scalar tuning.  It proves neither an unbounded-family
statement nor a construction obstruction.
"""

import argparse
from math import isfinite
from operator import index
import os
import shutil
import threading
import time

import numpy as np

try:
    _fs_actions
except NameError:
    from quartet_forced_run_scalar import (
        _fs_actions, _fs_profile, forced_run_scalar_bound,
        FS_WIDTH, FS_POWER, FS_PROFILE, FS_INITIAL, FS_BOUNDS, FS_NAMES,
        FC_INITIAL, FC_BOUNDS, FC_NAMES)


FB_RATIO = FS_WIDTH/FS_POWER
FB_LOCAL_COUNT = 35
FB_ENDPOINT_COUNT = 4


def _fb_model(logs, combined, compulsory_hit):
    local, endpoint = _fs_actions(logs, combined, compulsory_hit)
    if len(local) != FB_LOCAL_COUNT or len(endpoint) != FB_ENDPOINT_COUNT:
        raise ValueError('Expected all35 local and4 endpoint constraints')
    rows = local+endpoint
    return (np.asarray([row[0] for row in rows]),
            np.asarray([row[1] for row in rows]),
            np.asarray([row[2] for row in rows]))


def _fb_raw(logs, model, data):
    values = model[0]
    return (float(values[:FB_LOCAL_COUNT].max())
            + FB_RATIO*float(values[FB_LOCAL_COUNT:].max())
            + float(data['linear'] @ logs)+data['route_log2']/FS_POWER)


def _fb_barrier(point, mu, bounds, data, model):
    """Evaluate the barrier from cached model derivatives; None is infeasible."""
    dimension = len(data['linear'])
    logs, values = point[:dimension], model[0]
    lower, upper = logs-bounds[:, 0], bounds[:, 1]-logs
    slack = np.concatenate((point[dimension]-values[:FB_LOCAL_COUNT],
                            point[dimension+1]-values[FB_LOCAL_COUNT:]))
    if (not np.all(np.isfinite(point)) or not np.all(np.isfinite(slack))
            or np.any(lower <= 0) or np.any(upper <= 0) or np.any(slack <= 0)):
        return None
    objective_vector = np.concatenate((data['linear'], (1.0, FB_RATIO)))
    objective = float(objective_vector @ point)+data['route_log2']/FS_POWER
    jacobian = np.zeros((FB_LOCAL_COUNT+FB_ENDPOINT_COUNT, dimension+2))
    jacobian[:, :dimension] = model[1]
    jacobian[:FB_LOCAL_COUNT, dimension] = -1
    jacobian[FB_LOCAL_COUNT:, dimension+1] = -1
    multipliers = mu/slack
    gradient = objective_vector+jacobian.T @ multipliers
    gradient[:dimension] += mu/upper-mu/lower
    scaled_jacobian = jacobian/slack[:, None]
    hessian = mu*(scaled_jacobian.T @ scaled_jacobian)
    hessian[:dimension, :dimension] += np.einsum('i,ijk->jk', multipliers, model[2])
    diagonal = np.arange(dimension)
    hessian[diagonal, diagonal] += mu*(lower**-2+upper**-2)
    value = objective-mu*(float(np.log(slack).sum())
                          +float(np.log(lower).sum())+float(np.log(upper).sum()))
    if (not isfinite(value) or not np.all(np.isfinite(gradient))
            or not np.all(np.isfinite(hessian))):
        return None
    return dict(value=value, gradient=gradient, hessian=hessian, slack=slack,
                jacobian=jacobian, lower=lower, upper=upper, model=model,
                objective=objective, raw=_fb_raw(logs, model, data), mu=mu)


def _fb_dual(point, bounds, data, evaluation, names):
    """A finite-box dual tangent, with normalized epigraph multiplier masses."""
    dimension, slack = len(data['linear']), evaluation['slack']
    logs, model, mu = point[:dimension], evaluation['model'], evaluation['mu']
    local = slack[:FB_LOCAL_COUNT]
    endpoint = slack[FB_LOCAL_COUNT:]
    # Ratios avoid an unnecessary reciprocal of an extremely small slack.
    alpha = local.min()/local
    alpha /= alpha.sum()
    eta = endpoint.min()/endpoint
    eta *= FB_RATIO/eta.sum()
    weights = np.concatenate((alpha, eta))
    gradient = data['linear']+weights @ model[1]
    value = (float(data['linear'] @ logs)+float(weights @ model[0])
             +data['route_log2']/FS_POWER)
    corner = np.where(gradient >= 0, bounds[:, 0], bounds[:, 1])
    lower = value+float(gradient @ (corner-logs))
    face_tolerance = 1e-5*(bounds[:, 1]-bounds[:, 0])
    near_lower = evaluation['lower'] <= face_tolerance
    near_upper = evaluation['upper'] <= face_tolerance
    free = ~(near_lower | near_upper)
    lower_weights = mu/evaluation['lower']
    upper_weights = mu/evaluation['upper']
    stationarity = gradient-lower_weights+upper_weights
    faces = []
    for i, name in enumerate(names):
        if near_lower[i]:
            faces.append(dict(variable=name, side='lower', value=float(logs[i]),
                              bound=float(bounds[i, 0]), distance=float(evaluation['lower'][i])))
        if near_upper[i]:
            faces.append(dict(variable=name, side='upper', value=float(logs[i]),
                              bound=float(bounds[i, 1]), distance=float(evaluation['upper'][i])))
    diagnostics = dict(
        logs=dict(zip(names, map(float, logs))), barrier_mu=float(mu),
        finite_box_dual_lower_bits=FS_POWER*lower,
        local_multiplier_sum=float(alpha.sum()), endpoint_multiplier_sum=float(eta.sum()),
        endpoint_multiplier_target=FB_RATIO,
        raw_local_multiplier_sum=float((mu/local).sum()),
        raw_endpoint_multiplier_sum=float((mu/endpoint).sum()),
        normalized_dual_gradient=dict(zip(names, map(float, gradient))),
        box_lower_multipliers=dict(zip(names, map(float, lower_weights))),
        box_upper_multipliers=dict(zip(names, map(float, upper_weights))),
        free_coordinate_stationarity_inf=float(np.max(np.abs(gradient[free]))) if np.any(free) else 0.0,
        normalized_box_stationarity_inf=float(np.max(np.abs(stationarity))),
        raw_barrier_stationarity_inf=float(np.max(np.abs(evaluation['gradient']))),
        active_box_faces=faces, active_face_relative_tolerance=1e-5,
        local_epigraph_slack_min=float(local.min()), endpoint_epigraph_slack_min=float(endpoint.min()),
        epigraph_objective_bits=FS_POWER*evaluation['objective'],
        raw_primal_bits=FS_POWER*evaluation['raw'],
        nominal_exact_center_gap_bits=FS_POWER*(FB_LOCAL_COUNT+FB_ENDPOINT_COUNT+2*dimension)*mu,
        dual_feasibility_is_float_diagnostic=True)
    return lower, diagnostics


def _fb_newton_direction(gradient, hessian):
    """Solve a tiny diagonally scaled Newton system, with bounded ridge fallback."""
    scaling = 1/np.sqrt(np.maximum(np.diag(hessian), 1e-300))
    matrix = (hessian*scaling[:, None])*scaling[None, :]
    rhs = -scaling*gradient
    eye = np.eye(len(gradient))
    for ridge in (0.0, 1e-12, 1e-10, 1e-8, 1e-6):
        try:
            direction = scaling*np.linalg.solve(matrix+ridge*eye, rhs)
        except np.linalg.LinAlgError:
            continue
        if np.all(np.isfinite(direction)) and float(gradient @ direction) < 0:
            return direction, ridge
    direction = -gradient/np.maximum(np.diag(hessian), 1e-16)
    return direction, None


def forced_run_barrier_search(profile=FS_PROFILE, death_upper=None, initial=None,
                              bounds=None, combined=True, compulsory_hit=False,
                              max_evaluations=1024, mu_initial=0.01):
    """One bounded barrier schedule; no optimum or probability certification."""
    dimension, names = (8, FC_NAMES) if combined else (7, FS_NAMES)
    if initial is None:
        initial = FC_INITIAL if combined else FS_INITIAL
    if bounds is None:
        bounds = FC_BOUNDS if combined else FS_BOUNDS
    initial, bounds = np.asarray(initial, dtype=np.float64), np.asarray(bounds, dtype=np.float64)
    max_evaluations = index(max_evaluations)
    if (initial.shape != (dimension,) or bounds.shape != (dimension, 2)
            or not np.all(np.isfinite(initial)) or not np.all(np.isfinite(bounds))
            or np.any(bounds[:, 0] >= bounds[:, 1]) or np.any(bounds[5:, 0] < 0)
            or not 64 <= max_evaluations <= 4096 or not isfinite(mu_initial)
            or not 1e-6 <= mu_initial <= 0.1):
        raise ValueError('Require a finite valid model box,64..4096 evaluations,and1e-6<=mu<=.1')
    data = _fs_profile(profile, death_upper, combined)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    baseline = np.clip(initial, bounds[:, 0], bounds[:, 1])
    baseline_model = _fb_model(baseline, combined, compulsory_hit)
    best_logs, best_value = baseline.copy(), _fb_raw(baseline, baseline_model, data)
    best_lower, best_dual = float('-inf'), None
    margin = 1e-3*(bounds[:, 1]-bounds[:, 0])
    logs = np.clip(baseline, bounds[:, 0]+margin, bounds[:, 1]-margin)
    model = _fb_model(logs, combined, compulsory_hit)
    point = np.concatenate((logs, (float(model[0][:FB_LOCAL_COUNT].max())+FB_LOCAL_COUNT*mu_initial,
                                   float(model[0][FB_LOCAL_COUNT:].max())+FB_ENDPOINT_COUNT*mu_initial/FB_RATIO)))
    evaluations = iterations = backtracks = completed_stages = 0
    max_ridge = 0.0
    last_decrement = float('inf')
    termination = 'barrier_schedule_complete'
    final_evaluation = None
    last_mu = mu_initial

    def evaluate(candidate, mu):
        nonlocal evaluations, best_logs, best_value, best_lower, best_dual
        if evaluations >= max_evaluations:
            return None
        evaluations += 1
        candidate_logs = candidate[:dimension]
        if (not np.all(np.isfinite(candidate)) or np.any(candidate_logs <= bounds[:, 0])
                or np.any(candidate_logs >= bounds[:, 1])):
            return None
        candidate_model = _fb_model(candidate_logs, combined, compulsory_hit)
        value = _fb_barrier(candidate, mu, bounds, data, candidate_model)
        if value is not None:
            if value['raw'] < best_value:
                best_logs, best_value = candidate_logs.copy(), value['raw']
            lower, diagnostic = _fb_dual(candidate, bounds, data, value, names)
            if lower > best_lower:
                best_lower, best_dual = lower, diagnostic
        return value

    for stage in range(14):
        mu = mu_initial*(0.2**stage)
        last_mu = mu
        centered = False
        for _ in range(64):
            current = evaluate(point, mu)
            if current is None:
                if evaluations >= max_evaluations:
                    termination = 'evaluation_budget'
                    break
                raise ArithmeticError('Accepted barrier iterate lost strict feasibility')
            final_evaluation = current
            direction, ridge = _fb_newton_direction(current['gradient'], current['hessian'])
            if ridge is not None:
                max_ridge = max(max_ridge, ridge)
            descent = float(current['gradient'] @ direction)
            last_decrement = -descent
            if last_decrement <= max(1e-16, mu*1e-5):
                centered = True
                break
            if not isfinite(descent) or descent >= 0 or not np.all(np.isfinite(direction)):
                termination = 'numerical_direction_stall'
                break
            step = 1.0
            negative, positive = direction[:dimension] < 0, direction[:dimension] > 0
            if np.any(negative):
                step = min(step, float(np.min(-0.99*current['lower'][negative]/direction[:dimension][negative])))
            if np.any(positive):
                step = min(step, float(np.min(0.99*current['upper'][positive]/direction[:dimension][positive])))
            constraint_change = current['jacobian'] @ direction
            toward_boundary = constraint_change > 0
            if np.any(toward_boundary):
                step = min(step, float(np.min(0.99*current['slack'][toward_boundary]/constraint_change[toward_boundary])))
            accepted = False
            for _ in range(40):
                trial_point = point+step*direction
                trial = evaluate(trial_point, mu)
                if trial is not None and trial['value'] <= current['value']+1e-4*step*descent:
                    point, final_evaluation, accepted = trial_point, trial, True
                    break
                if evaluations >= max_evaluations:
                    termination = 'evaluation_budget'
                    break
                step *= 0.5
                backtracks += 1
            iterations += 1
            if not accepted:
                if termination != 'evaluation_budget':
                    termination = 'line_search_stall'
                break
        if centered:
            completed_stages += 1
        if evaluations >= max_evaluations:
            termination = 'evaluation_budget'
            break
        if termination in ('numerical_direction_stall', 'line_search_stall'):
            break
    if final_evaluation is None or best_dual is None:
        raise ArithmeticError('No feasible barrier diagnostic was produced')
    # Reporting calls are separate from the explicitly bounded optimizer budget.
    final_model = _fb_model(point[:dimension], combined, compulsory_hit)
    final_evaluation = _fb_barrier(point, last_mu, bounds, data, final_model)
    if final_evaluation is None:
        raise ArithmeticError('Final barrier report lost strict feasibility')
    final_lower, final_dual = _fb_dual(point, bounds, data, final_evaluation, names)
    if final_lower > best_lower:
        best_lower, best_dual = final_lower, final_dual
    result = forced_run_scalar_bound(profile, death_upper, best_logs,
                                     _combined=combined, compulsory_hit=compulsory_hit)
    result.update(
        optimizer='finite_box_epigraph_log_barrier', optimizer_evaluations=evaluations,
        optimizer_iterations=iterations, optimizer_backtracks=backtracks,
        optimizer_centered_stages=completed_stages, optimizer_termination=termination,
        maximum_scaled_newton_ridge=max_ridge, last_newton_decrement_squared=last_decrement,
        search_box=bounds.tolist(), all_local_constraints=FB_LOCAL_COUNT,
        all_endpoint_constraints=FB_ENDPOINT_COUNT, initial_raw_bound_bits=FS_POWER*_fb_raw(baseline, baseline_model, data),
        finite_box_dual_lower_diagnostic=FS_POWER*best_lower,
        finite_box_primal_dual_gap_diagnostic=result['raw_log2_probability_bound']-FS_POWER*best_lower,
        best_dual_diagnostics=best_dual, final_dual_diagnostics=final_dual,
        lower_diagnostic_is_exact_certificate=False, unbounded_parameter_domain_certified=False,
        optimizer_optimality_proved=False, scalar_family_impossibility_proved=False)
    return result


def forced_run_barrier_checks():
    """Small derivative and dual-feasibility checks; no production search."""
    derivatives = duals = epigraphs = 0
    for combined in (False, True):
        dimension, names = (8, FC_NAMES) if combined else (7, FS_NAMES)
        bounds = np.asarray(FC_BOUNDS if combined else FS_BOUNDS, dtype=np.float64)
        data = _fs_profile(FS_PROFILE, None, combined)
        logs = np.asarray(FC_INITIAL if combined else FS_INITIAL, dtype=np.float64).copy()
        logs[5:] = 0.4
        for compulsory_hit in (False, True):
            model = _fb_model(logs, combined, compulsory_hit)
            point = np.concatenate((logs, (model[0][:FB_LOCAL_COUNT].max()+0.7,
                                           model[0][FB_LOCAL_COUNT:].max()+1.3)))
            mu = 0.003
            current = _fb_barrier(point, mu, bounds, data, model)
            assert current is not None
            assert abs(current['raw']-_fb_raw(logs, model, data)) < 1e-14
            assert current['raw'] <= current['objective']
            assert np.linalg.eigvalsh(current['hessian']).min() > -1e-10
            for coordinate in range(dimension+2):
                delta = np.zeros(dimension+2)
                delta[coordinate] = 1e-5
                plus_point, minus_point = point+delta, point-delta
                plus = _fb_barrier(plus_point, mu, bounds, data,
                                   _fb_model(plus_point[:dimension], combined, compulsory_hit))
                minus = _fb_barrier(minus_point, mu, bounds, data,
                                    _fb_model(minus_point[:dimension], combined, compulsory_hit))
                assert plus is not None and minus is not None
                assert abs((plus['value']-minus['value'])/(2e-5)-current['gradient'][coordinate]) < 2e-6
                assert np.max(np.abs((plus['gradient']-minus['gradient'])/(2e-5)
                                     -current['hessian'][:, coordinate])) < 2e-5
                derivatives += 1
            lower, diagnostic = _fb_dual(point, bounds, data, current, names)
            assert abs(diagnostic['local_multiplier_sum']-1) < 1e-14
            assert abs(diagnostic['endpoint_multiplier_sum']-FB_RATIO) < 1e-16
            assert lower <= current['raw']+1e-12
            for candidate in (logs, bounds.mean(axis=1), bounds[:, 0], bounds[:, 1],
                              0.75*bounds[:, 0]+0.25*bounds[:, 1]):
                candidate_model = _fb_model(candidate, combined, compulsory_hit)
                assert lower <= _fb_raw(candidate, candidate_model, data)+1e-9
                duals += 1
            direction, _ = _fb_newton_direction(current['gradient'], current['hessian'])
            assert float(current['gradient'] @ direction) < 0
            invalid = point.copy()
            invalid[dimension] = model[0][:FB_LOCAL_COUNT].max()
            assert _fb_barrier(invalid, mu, bounds, data, model) is None
            invalid = point.copy()
            invalid[dimension+1] = model[0][FB_LOCAL_COUNT:].max()
            assert _fb_barrier(invalid, mu, bounds, data, model) is None
            epigraphs += 2
    # An independent convex quadratic model checks signs and the affine dual.
    bounds = np.asarray(((-2.0, 3.0),))
    data = dict(linear=np.asarray((0.0,)), route_log2=0.0)
    point = np.asarray((0.3, 0.8, 1.2))
    local_values = np.full(FB_LOCAL_COUNT, point[0]**2)
    endpoint_values = np.full(FB_ENDPOINT_COUNT, (point[0]-1)**2)
    gradients = np.concatenate((np.full(FB_LOCAL_COUNT, 2*point[0]),
                                np.full(FB_ENDPOINT_COUNT, 2*(point[0]-1))))[:, None]
    model = (np.concatenate((local_values, endpoint_values)), gradients,
             np.full((FB_LOCAL_COUNT+FB_ENDPOINT_COUNT, 1, 1), 2.0))
    current = _fb_barrier(point, 0.002, bounds, data, model)
    lower, _ = _fb_dual(point, bounds, data, current, ('x',))
    exact_optimum = FB_RATIO/(1+FB_RATIO)
    assert lower <= exact_optimum <= current['raw']
    print(dict(barrier_gradient_hessian_coordinates=derivatives,
               finite_box_dual_candidate_checks=duals, strict_epigraph_checks=epigraphs,
               independent_quadratic_dual_checks=1, full_code_certificate=False), flush=True)


def forced_run_barrier_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--search', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=FS_PROFILE)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--second-only', action='store_true', help='Use the seven-log second-cap model')
    parser.add_argument('--compulsory-hit', action='store_true')
    parser.add_argument('--logs', type=float, nargs='+')
    parser.add_argument('--evaluations', type=int, default=1024)
    parser.add_argument('--mu', type=float, default=0.01)
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
            forced_run_barrier_checks()
        else:
            print(forced_run_barrier_search(
                args.profile, args.death_upper, initial=args.logs, combined=not args.second_only,
                compulsory_hit=args.compulsory_hit, max_evaluations=args.evaluations,
                mu_initial=args.mu), flush=True)
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    forced_run_barrier_main()
