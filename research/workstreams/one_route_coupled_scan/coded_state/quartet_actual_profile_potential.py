"""Bounded FLOAT actual-support/MDS exponential-potential control.

This is the existing actual-S kernel, not a new code construction.  An old
nonempty quartet with r active coordinates permits at most r-1 new zeros
outside H.  An empty old quartet permits new support only in H; these are
exactly its forced births.  The local polynomial P(X,Y) marks old/new
support counts and includes the old-mask binomial C(4,r).

For phi(u)=alpha**u, all physical row actions are

 [X**u] P(X,alpha)**(n/4) / (C(n,u)*alpha**u).

A positive quartic convolution computes these coefficients together with
first/second exponent moments.  Their means and covariances give analytic
logarithmic gradients and Hessians.  This does not build an n-by-n kernel
or replace its coefficients by a Chernoff bound.  Local s and alpha occur
only as s*alpha; five-variable jets are lifted to six variables afterwards.

If every row is at most mu, the free-end moment is at most
mu**T * max(1,alpha**(-n)).  Start phi(0)=1.  At theta=0, D<=Dhi and
0<d<=1, the external correction is

 logC(M,h)-256*log(pi_h(p))-L*logs-Dhi*logd-G*logz
                                    -f*lognu-Cout*logeta.

There is one H union; no external route Q, field factor, or birth-mask
factor.  The death marker is independent and this first control is purely
geometric.  Actual live hits, deaths, births and the MDS exclusion remain
inside the kernel.  Regional exact quotas are upper bounded by the usual
conditioning cost.  No prior census complement is inferred automatically.

The optional convex epigraph search uses all n+1 physical rows, both
endpoint branches, and a stated finite box.  A normalized positive dual
mixture plus its affine tangent gives a finite-box lower diagnostic in
real arithmetic.  All production arithmetic here is FLOAT; neither its
upper nor lower output is an exact probability or optimum certificate.
No other helper globals are changed.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, isfinite, log, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np


AP_WIDTH, AP_UPDATES, AP_ROWS, AP_REGIONS = 512, 4096, 8192, 256
AP_N, AP_W = AP_WIDTH*AP_UPDATES, 943718
AP_PROFILE = (4096, 819200, 838861, 65536, 1081344)
AP_NAMES = ('log_s', 'log_d', 'log_z', 'log_nu', 'log_eta', 'log_alpha')
AP_INITIAL = (1.0, -1.0, -0.5, -3.0, -2.0, 0.0)
AP_BOUNDS = ((-8., 8.), (-16., 0.), (-8., 4.), (-20., 8.), (-12., 8.), (-8., 8.))
AP_LD = np.longdouble
_AP_LN2 = AP_LD(log(2))
_AP_PAIRS = tuple((i, j) for i in range(5) for j in range(i, 5))
_AP_I = np.asarray([i for i, _ in _AP_PAIRS])
_AP_J = np.asarray([j for _, j in _AP_PAIRS])
_AP_LIFT = np.asarray((0, 1, 2, 3, 4, 0))


@lru_cache(maxsize=8)
def _ap_terms(p):
    """Positive monomials (new,death,live,birth,outside), grouped by old r."""
    if not isfinite(p) or not 0 < p < 1:
        raise ValueError('Require 0<p<1')
    lp, lo = log2(p), log2(1-p)
    rows = []
    for r in range(5):
        exponents, constants = [], []
        for v in range(5):
            if r == 0:
                exponents.append((v, 0, 0, v, 0))
                constants.append(log2(comb(4, v))+v*lp)
                continue
            for kept in range(max(0, v-(4-r)), min(r, v)+1):
                death_count = comb(4, r)*comb(r, kept)*comb(4-r, v-kept)
                for outside_new in range(v+1):
                    for outside_zero in range(min(4-v, r-1)+1):
                        outside = outside_new+outside_zero
                        coefficient = death_count*comb(v, outside_new)*comb(4-v, outside_zero)
                        exponents.append((v, r-kept, 4, 0, outside))
                        constants.append(log2(coefficient)+(4-outside)*lp+outside*lo)
        rows.append((np.asarray(exponents, dtype=AP_LD), np.asarray(constants, dtype=AP_LD)))
    return tuple(rows)


@lru_cache(maxsize=8)
def _ap_logchoose_rows(width):
    return np.asarray([log2(comb(width, u)) for u in range(width+1)], dtype=AP_LD)


def _ap_local_jets(logs, p):
    reduced = np.asarray(logs[:5], dtype=AP_LD).copy()
    reduced[0] += logs[5]
    terms = _ap_terms(float(p))
    affine = [exponents @ reduced+constant for exponents, constant in terms]
    scale = max(float(row.max()) for row in affine)
    jets = np.zeros((5, 21), dtype=AP_LD)
    for r, ((exponents, _), values) in enumerate(zip(terms, affine)):
        weights = np.exp2(values-AP_LD(scale))
        if np.any(weights <= 0) or not np.all(np.isfinite(weights)):
            raise ArithmeticError('Local positive monomial underflow/overflow')
        jets[r, 0] = weights.sum()
        jets[r, 1:6] = weights @ exponents
        jets[r, 6:] = weights @ (exponents[:, _AP_I]*exponents[:, _AP_J])
    normalizer = jets[:, 0].max()
    jets /= normalizer
    return jets, AP_LD(scale)+np.log2(normalizer)


def actual_profile_rows(logs=AP_INITIAL, p=.5, width=AP_WIDTH):
    """Exact-coefficient FLOAT row logs and analytic derivatives, no matrix."""
    width = index(width)
    logs = np.asarray(logs, dtype=AP_LD)
    if width < 4 or width % 4 or width > 512 or logs.shape != (6,) or not np.all(np.isfinite(logs)):
        raise ValueError('Require width divisible by4 in4..512 and six finite logs')
    local, local_scale = _ap_local_jets(logs, p)
    buffers = (np.zeros((width+1, 21), dtype=AP_LD), np.zeros((width+1, 21), dtype=AP_LD))
    current, target = buffers
    current[0, 0] = 1
    degree, scale = 0, AP_LD(0)
    for _ in range(width//4):
        target[:degree+5] = 0
        old = current[:degree+1]
        for r in range(5):
            out = target[r:r+degree+1]
            scalar, first, second = local[r, 0], local[r, 1:6], local[r, 6:]
            out[:, 0] += old[:, 0]*scalar
            out[:, 1:6] += old[:, 1:6]*scalar+old[:, :1]*first
            out[:, 6:] += (old[:, 6:]*scalar+old[:, :1]*second
                           +old[:, 1+_AP_I]*first[_AP_J]
                           +old[:, 1+_AP_J]*first[_AP_I])
        degree += 4
        if np.any(target[:degree+1, 0] <= 0) or not np.all(np.isfinite(target[:degree+1])):
            raise ArithmeticError('Positive coefficient convolution lost numerical range')
        normalizer = target[:degree+1, 0].max()
        target[:degree+1] /= normalizer
        scale += local_scale+np.log2(normalizer)
        current, target = target, current
    coefficients = current[:, 0]
    mean = current[:, 1:6]/coefficients[:, None]
    second = current[:, 6:]/coefficients[:, None]
    covariance = np.zeros((width+1, 5, 5), dtype=AP_LD)
    for coordinate, (i, j) in enumerate(_AP_PAIRS):
        covariance[:, i, j] = second[:, coordinate]-mean[:, i]*mean[:, j]
        covariance[:, j, i] = covariance[:, i, j]
    populations = np.arange(width+1, dtype=AP_LD)
    values = np.log2(coefficients)+scale-_ap_logchoose_rows(width)-populations*logs[5]
    gradients = mean[:, _AP_LIFT].copy()
    gradients[:, 5] -= populations
    hessians = _AP_LN2*covariance[:, _AP_LIFT][:, :, _AP_LIFT]
    if not (np.all(np.isfinite(values)) and np.all(np.isfinite(gradients)) and np.all(np.isfinite(hessians))):
        raise ArithmeticError('Nonfinite coefficient row derivatives')
    # Newton is tiny; coefficient propagation above retains extended precision.
    return np.asarray(values, dtype=float), np.asarray(gradients, dtype=float), np.asarray(hessians, dtype=float)


def actual_profile_data(profile=AP_PROFILE, death_upper=None, p=.5):
    h, rank, volume, births, footprint = map(index, profile)
    if not (1 <= h <= AP_ROWS and 128*h <= rank <= (1 << 20)
            and 1 <= births <= volume <= footprint <= AP_N and 0 < p < 1):
        raise ValueError('Invalid actual-state profile or p')
    G, outside = footprint-births, rank-128*h
    CH = footprint-outside
    if G % 4 or G > AP_N-AP_WIDTH:
        raise ValueError('Invalid actual live volume')
    lower = max(0, births-AP_WIDTH)
    upper = min(footprint-volume, AP_W//2)
    Dhi = AP_W-volume if death_upper is None else index(death_upper)
    Dhi = min(Dhi, upper)
    empty = Dhi < lower or not (129*h <= CH <= 256*h and births <= CH
                                and 0 <= CH-births <= G
                                and 0 <= 256*h-CH <= AP_N-footprint
                                and outside <= 128*(AP_ROWS-h))
    choose = log2(comb(AP_ROWS, h))
    quota = choose+h*log2(p)+(AP_ROWS-h)*log2(1-p)
    constant = choose-AP_REGIONS*quota
    linear = np.asarray((-volume, -Dhi, -G, -births, -outside, 0.), dtype=float)/AP_UPDATES
    return dict(profile=(h, rank, volume, births, footprint), G=G, Cout=outside, C_H=CH,
                death_lower=lower, death_upper=Dhi, empty=empty, p=float(p),
                eligible_union=choose, regional_log2_quota=quota,
                constant=constant/AP_UPDATES, linear=linear)


def _ap_model(logs, p):
    values, gradients, hessians = actual_profile_rows(logs, p)
    endpoints = np.asarray((0., -logs[5]))
    endpoint_gradients = np.zeros((2, 6)); endpoint_gradients[1, 5] = -1
    return (np.concatenate((values, endpoints)), np.vstack((gradients, endpoint_gradients)),
            np.concatenate((hessians, np.zeros((2, 6, 6)))))


def _ap_raw(logs, model, data):
    return float(model[0][:513].max()+(AP_WIDTH/AP_UPDATES)*model[0][513:].max()
                 +data['linear'] @ logs+data['constant'])


def actual_profile_bound(profile=AP_PROFILE, death_upper=None, logs=AP_INITIAL, p=.5):
    logs = np.asarray(logs, dtype=float)
    if logs.shape != (6,) or not np.all(np.isfinite(logs)) or logs[1] > 0:
        raise ValueError('Six finite logs with log_d<=0 required')
    data = actual_profile_data(profile, death_upper, p)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    rows = actual_profile_rows(logs, p)
    growth = float(rows[0].max())
    endpoint = max(0., -AP_WIDTH*logs[5])
    correction = AP_UPDATES*float(data['linear'] @ logs+data['constant'])
    order = np.argsort(rows[0])[-8:][::-1]
    return dict(profile=data['profile'], death_interval=(data['death_lower'], data['death_upper']),
                logs=dict(zip(AP_NAMES, map(float, logs))), p=float(p), theta=0,
                log2_row_max=growth, terminal_log2_factor=endpoint,
                external_log2_factor=correction,
                raw_log2_probability_bound=AP_UPDATES*growth+endpoint+correction,
                largest_rows=[dict(population=int(u), log2_ratio=float(rows[0][u])) for u in order],
                H_union_log2=data['eligible_union'], log2_regional_quota=data['regional_log2_quota'],
                actual_live_hits=True, MDS_geometry_retained=True, forced_births_retained=True,
                field_factor_used=False, separate_routing_factor_used=False,
                census_complement_applied=False, positive_support_union_outside_fraction=True,
                row_count=AP_WIDTH+1, transitions=AP_UPDATES, free_endpoint=True,
                coefficient_chernoff_relaxation=False, floating_diagnostic_only=True,
                full_code_certificate=False)


def _ap_barrier(point, mu, bounds, data, model):
    logs = point[:6]
    lower, upper = logs-bounds[:, 0], bounds[:, 1]-logs
    slack = np.concatenate((point[6]-model[0][:513], point[7]-model[0][513:]))
    if (not np.all(np.isfinite(point)) or np.any(lower <= 0) or np.any(upper <= 0)
            or np.any(slack <= 0) or not np.all(np.isfinite(slack))):
        return None
    objective_vector = np.concatenate((data['linear'], (1., AP_WIDTH/AP_UPDATES)))
    objective = float(objective_vector @ point+data['constant'])
    jacobian = np.zeros((515, 8)); jacobian[:, :6] = model[1]
    jacobian[:513, 6] = -1; jacobian[513:, 7] = -1
    multipliers = mu/slack
    gradient = objective_vector+jacobian.T @ multipliers
    gradient[:6] += mu/upper-mu/lower
    scaled = jacobian/slack[:, None]
    hessian = mu*(scaled.T @ scaled)
    hessian[:6, :6] += np.einsum('i,ijk->jk', multipliers, model[2])
    diagonal = np.arange(6)
    hessian[diagonal, diagonal] += mu*(lower**-2+upper**-2)
    value = objective-mu*float(np.log(slack).sum()+np.log(lower).sum()+np.log(upper).sum())
    if not isfinite(value) or not np.all(np.isfinite(gradient)) or not np.all(np.isfinite(hessian)):
        return None
    return dict(value=value, gradient=gradient, hessian=hessian, lower=lower, upper=upper,
                slack=slack, jacobian=jacobian, model=model, raw=_ap_raw(logs, model, data))


def _ap_dual(point, bounds, data, evaluation):
    local, endpoint = evaluation['slack'][:513], evaluation['slack'][513:]
    alpha = local.min()/local; alpha /= alpha.sum()
    eta = endpoint.min()/endpoint; eta *= (AP_WIDTH/AP_UPDATES)/eta.sum()
    weights = np.concatenate((alpha, eta))
    gradient = data['linear']+weights @ evaluation['model'][1]
    logs = point[:6]
    value = float(data['linear'] @ logs+weights @ evaluation['model'][0]+data['constant'])
    corner = np.where(gradient >= 0, bounds[:, 0], bounds[:, 1])
    return value+float(gradient @ (corner-logs))


def actual_profile_search(profile=AP_PROFILE, death_upper=None, initial=AP_INITIAL,
                          bounds=AP_BOUNDS, p=.5, max_evaluations=256, mu_initial=.001):
    """Bounded convex epigraph diagnostic. Every saved bound uses all513 rows."""
    initial, bounds = np.asarray(initial, dtype=float), np.asarray(bounds, dtype=float)
    max_evaluations = index(max_evaluations)
    if (initial.shape != (6,) or bounds.shape != (6, 2) or not np.all(np.isfinite(initial))
            or not np.all(np.isfinite(bounds)) or np.any(bounds[:, 0] >= bounds[:, 1])
            or bounds[1, 1] > 0 or not 32 <= max_evaluations <= 1024
            or not isfinite(mu_initial) or not 1e-6 <= mu_initial <= .1):
        raise ValueError('Require finite six-variable box with logd<=0 and32..1024 evaluations')
    data = actual_profile_data(profile, death_upper, p)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    baseline = np.clip(initial, bounds[:, 0], bounds[:, 1])
    baseline_model = _ap_model(baseline, p)
    best_logs, best_raw = baseline.copy(), _ap_raw(baseline, baseline_model, data)
    initial_raw = best_raw
    best_lower = float('-inf')
    margin = .001*(bounds[:, 1]-bounds[:, 0])
    logs = np.clip(baseline, bounds[:, 0]+margin, bounds[:, 1]-margin)
    model = _ap_model(logs, p)
    point = np.concatenate((logs, (model[0][:513].max()+513*mu_initial,
                                   model[0][513:].max()+2*mu_initial/(AP_WIDTH/AP_UPDATES))))
    evaluations = iterations = stages = backtracks = 0
    termination = 'schedule_complete'

    def evaluate(candidate, mu):
        nonlocal evaluations, best_logs, best_raw, best_lower
        if evaluations >= max_evaluations:
            return None
        evaluations += 1
        if (not np.all(np.isfinite(candidate)) or np.any(candidate[:6] <= bounds[:, 0])
                or np.any(candidate[:6] >= bounds[:, 1])):
            return None
        value = _ap_barrier(candidate, mu, bounds, data, _ap_model(candidate[:6], p))
        if value is not None:
            if value['raw'] < best_raw:
                best_logs, best_raw = candidate[:6].copy(), value['raw']
            best_lower = max(best_lower, _ap_dual(candidate, bounds, data, value))
        return value

    for stage in range(12):
        mu = mu_initial*(.2**stage)
        centered = False
        for _ in range(48):
            current = evaluate(point, mu)
            if current is None:
                if evaluations >= max_evaluations:
                    termination = 'evaluation_budget'; break
                raise ArithmeticError('Accepted epigraph point lost feasibility')
            scaling = 1/np.sqrt(np.maximum(np.diag(current['hessian']), 1e-300))
            matrix = current['hessian']*scaling[:, None]*scaling[None, :]
            direction = None
            for ridge in (0., 1e-12, 1e-10, 1e-8, 1e-6):
                try:
                    candidate = scaling*np.linalg.solve(matrix+ridge*np.eye(8), -scaling*current['gradient'])
                except np.linalg.LinAlgError:
                    continue
                if np.all(np.isfinite(candidate)) and current['gradient'] @ candidate < 0:
                    direction = candidate; break
            if direction is None:
                direction = -current['gradient']/np.maximum(np.diag(current['hessian']), 1e-16)
            descent = float(current['gradient'] @ direction)
            if -descent <= max(1e-16, mu*1e-5):
                centered = True; break
            if not isfinite(descent) or descent >= 0:
                termination = 'numerical_direction_stall'; break
            step = 1.
            negative, positive = direction[:6] < 0, direction[:6] > 0
            if np.any(negative):
                step = min(step, float(np.min(-.99*current['lower'][negative]/direction[:6][negative])))
            if np.any(positive):
                step = min(step, float(np.min(.99*current['upper'][positive]/direction[:6][positive])))
            change = current['jacobian'] @ direction
            toward = change > 0
            if np.any(toward):
                step = min(step, float(np.min(.99*current['slack'][toward]/change[toward])))
            accepted = False
            for _ in range(32):
                candidate = point+step*direction
                trial = evaluate(candidate, mu)
                if trial is not None and trial['value'] <= current['value']+1e-4*step*descent:
                    point, accepted = candidate, True; break
                if evaluations >= max_evaluations:
                    termination = 'evaluation_budget'; break
                step *= .5; backtracks += 1
            iterations += 1
            if not accepted:
                if termination != 'evaluation_budget':
                    termination = 'line_search_stall'
                break
        if centered:
            stages += 1
        if evaluations >= max_evaluations:
            termination = 'evaluation_budget'; break
        if termination in ('line_search_stall', 'numerical_direction_stall'):
            break
    if not isfinite(best_lower):
        raise ArithmeticError('No feasible finite-box dual diagnostic')
    result = actual_profile_bound(profile, death_upper, best_logs, p)
    result.update(optimizer='six_log_actual_state_epigraph', optimizer_evaluations=evaluations,
                  optimizer_iterations=iterations, optimizer_centered_stages=stages,
                  optimizer_backtracks=backtracks, optimizer_termination=termination,
                  initial_raw_bound_bits=AP_UPDATES*initial_raw, search_box=bounds.tolist(),
                  finite_box_dual_lower_diagnostic=AP_UPDATES*best_lower,
                  finite_box_primal_dual_gap_diagnostic=result['raw_log2_probability_bound']-AP_UPDATES*best_lower,
                  lower_diagnostic_is_exact_certificate=False, optimizer_optimality_proved=False,
                  unbounded_parameter_domain_certified=False, scalar_family_impossibility_proved=False)
    return result


def _ap_physical_quartet(s, d, z, nu, eta, alpha, p):
    """Independent labelled old/new/H-mask enumeration, summed over new/H."""
    fixed_old = [Fraction(0) for _ in range(16)]
    for old in range(16):
        r = old.bit_count()
        for new in range(16):
            v, deaths = new.bit_count(), (old & ~new).bit_count()
            for eligible in range(16):
                if (r == 0 and new & ~eligible) or (r and ((~eligible & ~new) & 15).bit_count() >= r):
                    continue
                eligible_count = eligible.bit_count()
                outside = ((~eligible) & (15 if r else new) & 15).bit_count()
                fixed_old[old] += (p**eligible_count*(1-p)**(4-eligible_count)
                                  *(s*alpha)**v*d**deaths*z**(4 if r else 0)
                                  *nu**(v if r == 0 else 0)*eta**outside)
    return fixed_old


def actual_profile_checks():
    """Exact tiny physical comparisons and derivative checks; no large search."""
    physical = rows_checked = derivatives = endpoints = duals = 0
    parameter_sets = ((Fraction(3, 2), Fraction(1, 2), Fraction(2, 3), Fraction(1, 3),
                       Fraction(3, 4), Fraction(5, 4), Fraction(1, 2)),
                      (Fraction(1),)*6+(Fraction(1, 3),))
    for params in parameter_sets:
        s, d, z, nu, eta, alpha, p = params
        fixed = _ap_physical_quartet(*params)
        polynomial = [sum(fixed[old] for old in range(16) if old.bit_count() == r) for r in range(5)]
        for width in (4, 8):
            exact = [Fraction(1)]
            for _ in range(width//4):
                new = [Fraction(0)]*(len(exact)+4)
                for i, value in enumerate(exact):
                    for r, coefficient in enumerate(polynomial):
                        new[i+r] += value*coefficient
                exact = new
            logs = tuple(log2(float(value)) for value in (s, d, z, nu, eta, alpha))
            values, gradients, hessians = actual_profile_rows(logs, float(p), width)
            for u in range(width+1):
                target = exact[u]/comb(width, u)/alpha**u
                assert abs(values[u]-log2(float(target))) < 2e-11
                rows_checked += 1
            if width == 8:
                partitions = []
                for other in combinations(range(1, 8), 3):
                    first = (0,)+other
                    second = tuple(i for i in range(8) if i not in first)
                    partitions.append((first, second))
                for u in range(9):
                    old = (1 << u)-1
                    total = Fraction(0)
                    for partition in partitions:
                        weight = Fraction(1)
                        for quartet in partition:
                            local_old = sum(((old >> coordinate)&1) << i for i, coordinate in enumerate(quartet))
                            weight *= fixed[local_old]
                        total += weight
                    assert total/len(partitions) == exact[u]/comb(8, u)
                    physical += 1
            epsilon = 1e-5
            for coordinate in range(6):
                direction = np.zeros(6); direction[coordinate] = epsilon
                plus = actual_profile_rows(np.asarray(logs)+direction, float(p), width)
                minus = actual_profile_rows(np.asarray(logs)-direction, float(p), width)
                assert np.max(np.abs((plus[0]-minus[0])/(2*epsilon)-gradients[:, coordinate])) < 2e-7
                assert np.max(np.abs((plus[1]-minus[1])/(2*epsilon)-hessians[:, :, coordinate])) < 2e-6
                derivatives += 1
            for u in range(width+1):
                assert max(Fraction(1), alpha**(-width))*alpha**u >= 1
                endpoints += 1
    # A tiny synthetic convex row model audits normalization and dual signs.
    data = dict(linear=np.asarray((.1, -.2, 0., 0., 0., 0.)), constant=.3)
    bounds = np.asarray(((-2., 2.),)*6)
    logs = np.asarray((.2, -.1, .1, 0., -.2, .3))
    values = np.full(513, float(logs @ logs))
    gradients = np.tile(2*logs, (513, 1))
    hessians = np.tile(2*np.eye(6), (513, 1, 1))
    end_grad = np.zeros((2, 6)); end_grad[1, 5] = -1
    model = (np.concatenate((values, (0., -logs[5]))), np.vstack((gradients, end_grad)),
             np.concatenate((hessians, np.zeros((2, 6, 6)))))
    point = np.concatenate((logs, (values.max()+1, 1.)))
    evaluation = _ap_barrier(point, .001, bounds, data, model)
    lower = _ap_dual(point, bounds, data, evaluation)
    for candidate in (logs, np.zeros(6), bounds[:, 0], bounds[:, 1]):
        objective = float(candidate @ candidate+(AP_WIDTH/AP_UPDATES)*max(0., -candidate[5])
                          +data['linear'] @ candidate+data['constant'])
        assert lower <= objective+1e-12
        duals += 1
    print(dict(exact_partition_rows=physical, exact_coefficient_rows=rows_checked,
               gradient_hessian_coordinates=derivatives, free_endpoint_checks=endpoints,
               finite_box_dual_checks=duals, full_code_certificate=False), flush=True)


def actual_profile_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    mode.add_argument('--search', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=AP_PROFILE)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--logs', type=float, nargs=6, default=AP_INITIAL)
    parser.add_argument('--p', type=float, default=.5)
    parser.add_argument('--evaluations', type=int, default=256)
    parser.add_argument('--mu', type=float, default=.001)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124)); timer.daemon = True; timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            actual_profile_checks()
        elif args.probe:
            print(actual_profile_bound(args.profile, args.death_upper, args.logs, args.p), flush=True)
        else:
            print(actual_profile_search(args.profile, args.death_upper, args.logs, p=args.p,
                                        max_evaluations=args.evaluations, mu_initial=args.mu), flush=True)
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    actual_profile_main()
