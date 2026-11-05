"""Bounded FLOAT scalar control for the whole-run marker Omega_2.

E2=2*(closed forced runs)+sum_closed_unforced min(2,length).
The four coordinate labels are Z outside U, O in U outside Aprev,
J an unforced singleton in Aprev, and C every other active position.
For local counts(o,j,c), r=o+j+c<=4, potential
V=alpha^|U| rho^|Aprev| sigma^|J| has positive one-quartet action

 kappa^(-j-2c)
 +z^4 alpha^4(beta*rho)^r sigma^o
       *[(1+w*rho)^(4-r)-1[r=0]].

Divide by alpha^r rho^(j+c) sigma^j and maximize over35 compositions.
The endpoint factor is
 max(1+w,beta/alpha,beta/(alpha*rho),beta/(alpha*rho*sigma))^512.
No run closes at the endpoint. The local labels implement a single
global whole-run statistic in a scalar supersolution, not a state DP.

For D<=D0, m=L-f and beta=b*kappa with b,kappa>=1, the profile correction
is logQ-Glogz-flogw-m logb-(m-2D0)logkappa. Births are integrated;
there is no extra support union, birth binomial, or death penalty.
The bounded smoothed convex search is FLOAT only, with no optimality
claim. Exact small tests separately verify the whole-run label semantics.

The optional eight-log combined interface retains the first cap too:
beta=b*kappa1*kappa2, J closes with(kappa1*kappa2)^-1 and C with
(kappa1*kappa2^2)^-1. Its external correction contains
-(m-D0)logkappa1-(m-2D0)logkappa2. Setting kappa2=sigma=1 exactly
embeds the previous first-cap scalar control. The seven-log interface
and its defaults remain unchanged; use --combined to opt in.

With --compulsory-hit, a genuine selected next quartet containing exactly
one current allowed position forces that position to be nonzero. Only
the selected(o,j,c)=(1,0,0), birth-degree0 term changes: its new label
is C rather than J, so its sigma factor is removed. This promotes a
closed compulsory unforced singleton from charge1 to charge2. It does
not alter f, prescribe an amplitude, or change birth normalization.
The formal full endpoint supplies no such hitting condition.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import comb, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np

try:
    _ns_logsum
except NameError:
    from quartet_nested_selector_scalar import (
        _ns_logsum, _ns_smooth_max, nested_selector_profile, nested_selector_bound)


FS_WIDTH, FS_UPDATES = 512, 4096
FS_POWER = 128*(FS_UPDATES-1)
FS_PROFILE = (4096, 819200, 838861, 65536, 1081344)
FS_INITIAL = (-0.9301908045807551, -1.343939327491652,
              0.17125048238006038, -0.45762697461224744, 0.0, 0.0, 0.9488661799848431)
FS_BOUNDS = ((-8.0, 4.0), (-24.0, 8.0), (-8.0, 8.0), (-12.0, 8.0),
             (-12.0, 12.0), (0.0, 12.0), (0.0, 12.0))
FS_NAMES = ('log_z', 'log_w', 'log_alpha', 'log_rho', 'log_sigma', 'log_b', 'log_kappa')
FC_INITIAL = FS_INITIAL[:6]+(FS_INITIAL[6], 0.0)
FC_BOUNDS = FS_BOUNDS[:6]+(FS_BOUNDS[6], FS_BOUNDS[6])
FC_NAMES = FS_NAMES[:6]+('log_kappa1', 'log_kappa2')


def _fs_cases(combined=False, compulsory_hit=False):
    cases = []
    for o in range(5):
        for j in range(5-o):
            for c in range(5-o-j):
                r = o+j+c
                exponents = [(0, 0, -r, -j-c, -j, 0, -j-2*c)]
                coefficients = [0.0]
                for ell in range(int(r == 0), 5-r):
                    sigma_power = o-j-int(compulsory_hit and o == r == 1 and ell == 0)
                    exponents.append((4, ell, 4-r, o+ell, sigma_power, r, r))
                    coefficients.append(log2(comb(4-r, ell)))
                if combined:
                    exponents = ([exponents[0][:6]+(-j-c,)+exponents[0][6:]]
                                 +[exponent[:6]+(r,)+exponent[6:] for exponent in exponents[1:]])
                cases.append((o, j, c, np.asarray(exponents, dtype=np.float64),
                              np.asarray(coefficients, dtype=np.float64)))
    return tuple(cases)


_FS_CASES = _fs_cases()
_FC_CASES = _fs_cases(True)
_FS_HIT_CASES = _fs_cases(False, True)
_FC_HIT_CASES = _fs_cases(True, True)


def _fs_actions(logs, combined=False, compulsory_hit=False):
    cases = (_FC_HIT_CASES if combined else _FS_HIT_CASES) if compulsory_hit else (_FC_CASES if combined else _FS_CASES)
    dimension = 8 if combined else 7
    rows = [_ns_logsum(exponents @ logs+coefficients, exponents)
            for _, _, _, exponents, coefficients in cases]
    exponents = np.zeros((2, dimension))
    exponents[1, 1] = 1
    endpoints = [_ns_logsum(exponents @ logs, exponents)]
    for exponent in ((0, 0, -1, 0, 0, 1, 1),
                     (0, 0, -1, -1, 0, 1, 1),
                     (0, 0, -1, -1, -1, 1, 1)):
        if combined:
            exponent = exponent[:6]+(1,)+exponent[6:]
        vector = np.asarray(exponent, dtype=np.float64)
        endpoints.append((float(vector @ logs), vector, np.zeros((dimension, dimension))))
    return rows, endpoints


def _fs_profile(profile, death_upper, combined=False):
    data = nested_selector_profile(profile, death_upper)
    m, D0 = data['m'], data['death_upper']
    data['linear'] = np.asarray((-data['G'], -data['profile'][3], 0, 0, 0, -m, -(m-2*D0)), dtype=np.float64)/FS_POWER
    if combined:
        data['linear'] = np.insert(data['linear'], 6, -(m-D0)/FS_POWER)
    return data


def _fs_objective(logs, temperature, data, combined=False, compulsory_hit=False):
    rows, endpoints = _fs_actions(logs, combined, compulsory_hit)
    rvalues = np.asarray([row[0] for row in rows])
    evalues = np.asarray([row[0] for row in endpoints])
    local = _ns_smooth_max(rvalues, np.asarray([row[1] for row in rows]),
                           np.asarray([row[2] for row in rows]), temperature)
    endpoint = _ns_smooth_max(evalues, np.asarray([row[1] for row in endpoints]),
                              np.asarray([row[2] for row in endpoints]), temperature)
    ratio = FS_WIDTH/FS_POWER
    correction = float(data['linear'] @ logs)+data['route_log2']/FS_POWER
    value = local[0]+ratio*endpoint[0]+correction
    gradient = local[1]+ratio*endpoint[1]+data['linear']
    hessian = local[2]+ratio*endpoint[2]
    raw = float(rvalues.max())+ratio*float(evalues.max())+correction
    return value, gradient, hessian, raw


def forced_run_scalar_bound(profile=FS_PROFILE, death_upper=None, logs=FS_INITIAL, _combined=False,
                            compulsory_hit=False):
    logs = np.asarray(logs, dtype=np.float64)
    dimension, cases, names = (8, _FC_CASES, FC_NAMES) if _combined else (7, _FS_CASES, FS_NAMES)
    if compulsory_hit:
        cases = _FC_HIT_CASES if _combined else _FS_HIT_CASES
    if logs.shape != (dimension,) or not np.all(np.isfinite(logs)) or np.any(logs[5:] < 0):
        raise ValueError('Require the model log count, with logb/cap logs nonnegative')
    data = _fs_profile(profile, death_upper, _combined)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    rows, endpoints = _fs_actions(logs, _combined, compulsory_hit)
    values = [row[0] for row in rows]
    log_lambda, log_endpoint = max(values), max(row[0] for row in endpoints)
    exponent = FS_POWER*log_lambda+FS_WIDTH*log_endpoint+data['route_log2']+FS_POWER*float(data['linear'] @ logs)
    order = sorted(range(35), key=values.__getitem__, reverse=True)
    result = dict(profile=data['profile'], death_interval=(data['death_lower'], data['death_upper']),
                status='combined_first_second_cap_control' if _combined else 'forced_run_low_death_control',
                logs=dict(zip(names, map(float, logs))), log_beta=float(np.sum(logs[5:])),
                log2_lambda=log_lambda, log2_endpoint_per_coordinate=log_endpoint,
                route_log2=data['route_log2'], raw_log2_probability_bound=exponent,
                tightest_compositions=[dict(o=cases[i][0], j=cases[i][1], c=cases[i][2],
                                            log2_ratio=values[i]) for i in order[:5]],
                closed_forced_singleton_charge=2, closed_unforced_singleton_charge=1,
                closed_compulsory_unforced_singleton_charge=2 if compulsory_hit else 1,
                no_artificial_terminal_closure=True, births_integrated=True,
                H_union_already_in_route=True, floating_diagnostic_only=True,
                optimizer_optimality_proved=False, scalar_family_impossibility_proved=False,
                first_cap_retained=bool(_combined), compulsory_hit_singletons_retained=bool(compulsory_hit),
                amplitude_birth_normalization_unchanged=True, full_code_certificate=False)
    if _combined:
        result.update(log_kappa1=float(logs[6]), log_kappa2=float(logs[7]))
    else:
        result['log_kappa'] = float(logs[6])
    return result


def forced_run_scalar_search(profile=FS_PROFILE, death_upper=None, initial=FS_INITIAL,
                             bounds=FS_BOUNDS, max_evaluations=1024, _combined=False, compulsory_hit=False):
    """Bounded convex smoothing/Newton control; seven-variable default is unchanged."""
    dimension = 8 if _combined else 7
    data = _fs_profile(profile, death_upper, _combined)
    if data['empty']:
        return dict(profile=data['profile'], status='empty', full_code_certificate=False)
    max_evaluations = index(max_evaluations)
    initial, bounds = np.asarray(initial, dtype=np.float64), np.asarray(bounds, dtype=np.float64)
    if (initial.shape != (dimension,) or bounds.shape != (dimension, 2) or not np.all(np.isfinite(initial))
            or not np.all(np.isfinite(bounds)) or np.any(bounds[:, 0] >= bounds[:, 1])
            or np.any(bounds[5:, 0] < 0) or not 32 <= max_evaluations <= 2048):
        raise ValueError('Invalid finite model box or32..2048 evaluation budget')
    point = np.clip(initial, bounds[:, 0], bounds[:, 1])
    best_point, best_value = point.copy(), float('inf')
    evaluations = iterations = 0
    last_temperature = 0.25

    def evaluate(candidate, temperature):
        nonlocal evaluations, best_point, best_value
        if evaluations >= max_evaluations:
            return None
        evaluations += 1
        result = _fs_objective(candidate, temperature, data, _combined, compulsory_hit)
        if result[3] < best_value:
            best_point, best_value = candidate.copy(), result[3]
        return result

    for temperature in (1/4, 1/16, 1/64, 1/256, 1/1024, 1/4096, 1/16384, 1/65536):
        last_temperature = temperature
        for _ in range(28):
            current = evaluate(point, temperature)
            if current is None:
                break
            value, gradient, hessian, _ = current
            fixed = (((point <= bounds[:, 0]+1e-11) & (gradient >= 0))
                     | ((point >= bounds[:, 1]-1e-11) & (gradient <= 0)))
            free = np.flatnonzero(~fixed)
            if not len(free) or np.max(np.abs(gradient[free])) < 1e-9:
                break
            reduced = hessian[np.ix_(free, free)]
            ridge = max(1e-10, 1e-8*float(np.trace(reduced))/len(free))
            direction = np.zeros(dimension)
            try:
                direction[free] = np.linalg.solve(reduced+ridge*np.eye(len(free)), -gradient[free])
            except np.linalg.LinAlgError:
                direction[free] = -gradient[free]
            direction *= min(1.0, 2/max(1e-300, float(np.max(np.abs(direction)))))
            direction = np.clip(point+direction, bounds[:, 0], bounds[:, 1])-point
            if float(gradient @ direction) >= -1e-15:
                direction = np.clip(point-gradient/max(1.0, float(np.max(np.abs(gradient)))),
                                    bounds[:, 0], bounds[:, 1])-point
            descent = float(gradient @ direction)
            if descent >= -1e-15:
                break
            accepted, step = False, 1.0
            for _ in range(24):
                candidate = point+step*direction
                trial = evaluate(candidate, temperature)
                if trial is None:
                    break
                if trial[0] <= value+1e-4*step*descent:
                    point, accepted = candidate, True
                    break
                step *= 0.5
            iterations += 1
            if not accepted:
                break
        if evaluations >= max_evaluations:
            break
    smooth, gradient, _, _ = _fs_objective(best_point, last_temperature, data, _combined, compulsory_hit)
    corner = np.where(gradient >= 0, bounds[:, 0], bounds[:, 1])
    error = last_temperature*(log2(35)+(FS_WIDTH/FS_POWER)*2)
    lower = FS_POWER*(smooth+float(gradient @ (corner-best_point))-error)
    result = forced_run_scalar_bound(profile, death_upper, best_point, _combined, compulsory_hit)
    result.update(optimizer_evaluations=evaluations, optimizer_iterations=iterations,
                  search_box=bounds.tolist(), final_temperature=last_temperature,
                  smoothing_error_bits=FS_POWER*error, finite_box_lower_diagnostic=lower,
                  finite_box_gap_diagnostic=result['raw_log2_probability_bound']-lower,
                  lower_diagnostic_is_exact_certificate=False,
                  unbounded_parameter_domain_certified=False)
    return result


def forced_run_combined_bound(profile=FS_PROFILE, death_upper=None, logs=FC_INITIAL, compulsory_hit=False):
    """Eight logs(z,w,alpha,rho,sigma,b,kappa1,kappa2), with b/kappa1/kappa2>=1."""
    return forced_run_scalar_bound(profile, death_upper, logs, _combined=True, compulsory_hit=compulsory_hit)


def forced_run_combined_search(profile=FS_PROFILE, death_upper=None, initial=FC_INITIAL,
                               bounds=FC_BOUNDS, max_evaluations=1024, compulsory_hit=False):
    return forced_run_scalar_search(profile, death_upper, initial, bounds, max_evaluations,
                                    _combined=True, compulsory_hit=compulsory_hit)


def _fs_subsets(mask):
    subset = mask
    while True:
        yield subset
        if subset == 0:
            break
        subset = (subset-1)&mask


def _fs_transitions(old, active, singleton, partition, beta, kappa, z, w, kappa1=Fraction(1),
                    compulsory_hit=False):
    for choice in range(1 << len(partition)):
        blocks = [block for j, block in enumerate(partition) if choice >> j & 1]
        new = sum(blocks)
        for births in _fs_subsets(new&~old):
            if any(not (block&(old|births)) for block in blocks):
                continue
            new_active = (old&new)|births
            new_singleton = (old&~active)&new
            if compulsory_hit:
                for block in blocks:
                    if (block&new_active).bit_count() == 1:
                        new_singleton &= ~block
            charge = (singleton&~new).bit_count()+2*(active&~singleton&~new).bit_count()
            weight = (beta**(old&new).bit_count()*kappa**(-charge)*kappa1**(-(active&~new).bit_count())
                      *z**new.bit_count()*w**births.bit_count())
            yield new, new_active, new_singleton, births, charge, weight


def _fs_run_charge(allowed, births, width, compulsory=None):
    charge = 0
    for coordinate in range(width):
        start = None
        for t in range(len(allowed)+1):
            occupied = t < len(allowed) and bool(allowed[t] >> coordinate & 1)
            if occupied and start is None:
                start = t
            if not occupied and start is not None:
                if t < len(allowed):
                    forced = bool(births[start] >> coordinate & 1)
                    if compulsory is not None:
                        forced |= any((compulsory[s] >> coordinate) & 1 for s in range(start, t))
                    charge += 2 if forced else min(2, t-start)
                start = None
    return charge


def _fs_compulsory_hit_checks():
    local = classifications = changed = 0
    z, w, alpha, rho, sigma = Fraction(1, 3), Fraction(1, 4), Fraction(3, 2), Fraction(2, 3), Fraction(2)
    for combined in (False, True):
        cases = _FC_HIT_CASES if combined else _FS_HIT_CASES
        original = _FC_CASES if combined else _FS_CASES
        kappa1, kappa2, beta = Fraction(3, 2) if combined else Fraction(1), Fraction(2), Fraction(3)
        markers = ((z, w, alpha, rho, sigma, beta/(kappa1*kappa2), kappa1, kappa2) if combined
                   else (z, w, alpha, rho, sigma, beta/kappa2, kappa2))
        for case, old_case in zip(cases, original):
            o, j, c, exponents, _ = case
            r = o+j+c
            old, active, singleton = (1 << r)-1, ((1 << (j+c))-1) << o, ((1 << j)-1) << o
            total = Fraction(0)
            for new, new_active, new_singleton, births, _, weight in _fs_transitions(
                    old, active, singleton, (15,), beta, kappa2, z, w, kappa1, True):
                independently_classified = 0
                for coordinate in range(4):
                    bit = 1 << coordinate
                    if not new_active&bit:
                        continue
                    must_be_nonzero = new_active.bit_count() == 1
                    if not (active&bit or births&bit or must_be_nonzero):
                        independently_classified |= bit
                assert new_singleton == independently_classified
                classifications += 1
                total += weight*alpha**new.bit_count()*rho**new_active.bit_count()*sigma**new_singleton.bit_count()
            multiplicities = [1]+[comb(4-r, ell) for ell in range(int(r == 0), 5-r)]
            model = Fraction(0)
            for multiplicity, powers in zip(multiplicities, exponents):
                term = Fraction(multiplicity)
                for marker, power in zip(markers, powers):
                    term *= marker**int(power)
                model += term
            assert model == total/(alpha**r*rho**(j+c)*sigma**j)
            difference = exponents-old_case[3]
            if o == 1 and j == c == 0:
                expected = np.zeros_like(difference)
                expected[1, 4] = -1
                assert np.array_equal(difference, expected)
                changed += 1
            else:
                assert not np.any(difference)
            local += 1
    def mask(*positions):
        return sum(1 << (position-1) for position in positions)
    core = mask(1, 5, 6, 7)
    live = (0, mask(1, 2, 3, 4), 255, core, core, core)
    born = (mask(1), 0, 0, 0, 0, 0)
    partitions = ((mask(1, 2, 3, 4), mask(5, 6, 7, 8)),
                  (mask(1, 2, 3, 5), mask(4, 6, 7, 8)))+((core, 255^core),)*3
    allowed = tuple((live[t]&live[t+1])|born[t] for t in range(5))+(live[-1]|born[-1],)
    assert allowed == (mask(1), mask(1, 2, 3, 4), core, core, core, core)
    compulsory = []
    for t, partition in enumerate(partitions):
        forced = 0
        for block in partition:
            if block&live[t+1] and (block&allowed[t]).bit_count() == 1:
                forced |= block&allowed[t]
        compulsory.append(forced)
    compulsory.append(0)  # The formal full endpoint supplies no quartet-hit test.
    assert _fs_run_charge(allowed, born, 8) == 3
    assert _fs_run_charge(allowed, born, 8, compulsory) == 4
    for retain in (False, True):
        active = singleton = charge = 0
        for t, partition in enumerate(partitions):
            matches = [entry for entry in _fs_transitions(live[t], active, singleton, partition,
                       Fraction(2), Fraction(2), Fraction(1, 2), Fraction(1, 4), compulsory_hit=retain)
                       if entry[0] == live[t+1] and entry[3] == born[t]]
            assert len(matches) == 1
            _, active, singleton, _, extra, _ = matches[0]
            assert active == allowed[t]
            charge += extra
        assert charge == (4 if retain else 3)
    support = (mask(1),)*4+(mask(1, 5),)*2
    assert all(not (support[t]&~allowed[t] or born[t]&~support[t]) for t in range(6))
    deaths = sum((support[t]&~support[t+1]).bit_count() for t in range(5))
    assert sum(value.bit_count() for value in support) == 8 and deaths == 0
    assert sum(value.bit_count() for value in live)//4+support[-1].bit_count() == 8
    assert compulsory[1]&mask(4) and not allowed[2]&mask(4)
    assert not support[1]&mask(4)  # This packing witness fails the required live hit.
    holes = sum(value.bit_count() for value in allowed)-8
    assert holes == 13 and holes+2*deaths >= 4
    # Thus the hit requirement forces a death, but the E2 inequality alone
    # does not exclude this example: its other holes can pay the charge.
    old_logs = (FC_INITIAL[0], FC_INITIAL[1], FC_INITIAL[2], FC_INITIAL[3], FC_INITIAL[5], FC_INITIAL[6])
    old_bound = nested_selector_bound(FS_PROFILE, None, old_logs)
    embedded = forced_run_combined_bound(compulsory_hit=True)
    assert abs(old_bound['raw_log2_probability_bound']-embedded['raw_log2_probability_bound']) < 1e-7
    point, data = np.asarray(FC_INITIAL), _fs_profile(FS_PROFILE, None, True)
    _, gradient, hessian, _ = _fs_objective(point, 1/16, data, True, True)
    for coordinate in range(8):
        delta = np.zeros(8)
        delta[coordinate] = 1e-5
        plus = _fs_objective(point+delta, 1/16, data, True, True)
        minus = _fs_objective(point-delta, 1/16, data, True, True)
        assert abs((plus[0]-minus[0])/(2e-5)-gradient[coordinate]) < 1e-6
        assert np.max(np.abs((plus[1]-minus[1])/(2e-5)-hessian[:, coordinate])) < 1e-5
    return dict(compulsory_hit_local_checks=local, independent_singleton_classifications=classifications,
                changed_monomials_checked=changed, strict_width8_T6_example=True,
                strict_example_old_E2=3, strict_example_new_E2=4,
                strict_example_not_excluded_by_E2_alone=True, compulsory_gradient_hessian_coordinates=8)


def forced_run_scalar_checks():
    local = physical = endpoint = histories = packing = 0
    settings = ((Fraction(2), Fraction(2), Fraction(1, 2), Fraction(1, 8),
                 Fraction(2), Fraction(1, 2), Fraction(3, 2)),
                (Fraction(3), Fraction(3, 2), Fraction(1, 3), Fraction(1, 4),
                 Fraction(3, 2), Fraction(2, 3), Fraction(1, 2)))
    for beta, kappa, z, w, alpha, rho, sigma in settings:
        ratios = []
        for o, j, c, _, _ in _FS_CASES:
            r = o+j+c
            old, active = (1 << r)-1, ((1 << (j+c))-1) << o
            singleton = ((1 << j)-1) << o
            exact = sum(weight*alpha**new.bit_count()*rho**new_active.bit_count()*sigma**new_singleton.bit_count()
                        for new, new_active, new_singleton, _, _, weight in
                        _fs_transitions(old, active, singleton, (15,), beta, kappa, z, w))
            choices = sum(Fraction(comb(4-r, ell))*(w*rho)**ell for ell in range(int(r == 0), 5-r))
            formula = kappa**(-j-2*c)+z**4*alpha**4*(beta*rho)**r*sigma**o*choices
            assert exact == formula
            ratios.append(exact/(alpha**r*rho**(j+c)*sigma**j))
            local += 1
        lam = max(ratios)
        terminal = max(1+w, beta/alpha, beta/(alpha*rho), beta/(alpha*rho*sigma))
        partitions = [(block, 255^block) for rest in combinations(range(1, 8), 3)
                      for block in (1+sum(1 << i for i in rest),)]
        for o, j, c in ((0, 0, 0), (4, 0, 0), (3, 1, 0), (0, 4, 0), (0, 0, 4),
                         (2, 1, 1), (8, 0, 0), (0, 8, 0), (0, 0, 8), (4, 2, 2)):
            r = o+j+c
            old, active, singleton = (1 << r)-1, ((1 << (j+c))-1) << o, ((1 << j)-1) << o
            potential = alpha**r*rho**(j+c)*sigma**j
            assert beta**r*(1+w)**(8-r) <= terminal**8*potential
            endpoint += 1
            for partition in partitions:
                exact = sum(weight*alpha**new.bit_count()*rho**new_active.bit_count()*sigma**new_singleton.bit_count()
                            for new, new_active, new_singleton, _, _, weight in
                            _fs_transitions(old, active, singleton, partition, beta, kappa, z, w))
                assert exact <= lam**2*potential
                physical += 1
        states = [(0, 0, 0, (), (), 0, Fraction(1))]
        for updates in range(1, 4):
            moment = Fraction(0)
            for old, active, singleton, allowed, born, charge, weight in states:
                for final_births in _fs_subsets(15&~old):
                    complete_allowed = allowed+((old|final_births),)
                    complete_births = born+(final_births,)
                    assert charge == _fs_run_charge(complete_allowed, complete_births, 4)
                    moment += weight*beta**old.bit_count()*w**final_births.bit_count()
                    histories += 1
            assert moment <= lam**(updates-1)*terminal**4
            following = []
            for old, active, singleton, allowed, born, charge, weight in states:
                for new, new_active, new_singleton, births, extra, step in _fs_transitions(
                        old, active, singleton, (15,), beta, kappa, z, w):
                    following.append((new, new_active, new_singleton, allowed+(new_active,),
                                      born+(births,), charge+extra, weight*step))
            states = following
    for length in range(1, 6):
        for mask in range(1 << length):
            starts = mask & ~(mask << 1)
            for forced in _fs_subsets(starts):
                allowed = tuple((mask >> t)&1 for t in range(length))
                births = tuple((forced >> t)&1 for t in range(length))
                charge = _fs_run_charge(allowed, births, 1)
                for support in _fs_subsets(mask):
                    if forced & ~support:
                        continue
                    deaths = (support & ~(support >> 1) & ((1 << (length-1))-1)).bit_count()
                    holes = mask.bit_count()-support.bit_count()
                    assert holes+2*deaths >= charge
                    packing += 1
    data, point = _fs_profile(FS_PROFILE, None), np.asarray(FS_INITIAL)
    value, gradient, hessian, raw = _fs_objective(point, 1/16, data)
    error = (1/16)*(log2(35)+(FS_WIDTH/FS_POWER)*2)
    assert raw <= value+1e-12 and value <= raw+error+1e-12
    assert np.linalg.eigvalsh(hessian).min() > -1e-9
    for coordinate in range(7):
        delta = np.zeros(7)
        delta[coordinate] = 1e-5
        plus, minus = _fs_objective(point+delta, 1/16, data), _fs_objective(point-delta, 1/16, data)
        assert abs((plus[0]-minus[0])/(2e-5)-gradient[coordinate]) < 1e-6
        assert np.max(np.abs((plus[1]-minus[1])/(2e-5)-hessian[:, coordinate])) < 1e-5
    combined_local = embedding = 0
    beta, kappa1, kappa2 = Fraction(3), Fraction(3, 2), Fraction(2)
    z, w, alpha, rho, sigma = Fraction(1, 3), Fraction(1, 4), Fraction(3, 2), Fraction(2, 3), Fraction(1, 2)
    for o, j, c, exponents, _ in _FC_CASES:
        r = o+j+c
        old, active, singleton = (1 << r)-1, ((1 << (j+c))-1) << o, ((1 << j)-1) << o
        physical_sum = sum(weight*alpha**new.bit_count()*rho**new_active.bit_count()*sigma**new_singleton.bit_count()
                           for new, new_active, new_singleton, _, _, weight in
                           _fs_transitions(old, active, singleton, (15,), beta, kappa2, z, w, kappa1))
        births = sum(Fraction(comb(4-r, ell))*(w*rho)**ell for ell in range(int(r == 0), 5-r))
        selected = z**4*alpha**4*(beta*rho)**r*sigma**o*births
        assert physical_sum == kappa1**(-j-c)*kappa2**(-j-2*c)+selected
        multiplicities = [1]+[comb(4-r, ell) for ell in range(int(r == 0), 5-r)]
        markers = (z, w, alpha, rho, sigma, beta/(kappa1*kappa2), kappa1, kappa2)
        model = Fraction(0)
        for multiplicity, powers in zip(multiplicities, exponents):
            term = Fraction(multiplicity)
            for marker, power in zip(markers, powers):
                term *= marker**int(power)
            model += term
        assert model == physical_sum/(alpha**r*rho**(j+c)*sigma**j)
        combined_local += 1
        # With kappa2=sigma=1, each composition is exactly the old first-cap
        # local action, irrespective of the split of active positions into J/C.
        combined_first_only = sum(weight*alpha**new.bit_count()*rho**new_active.bit_count()
                                  for new, new_active, _, _, _, weight in
                                  _fs_transitions(old, active, singleton, (15,), beta, Fraction(1), z, w, kappa1))
        old_first_only = kappa1**(-(j+c))+z**4*alpha**4*(beta*rho)**r*births
        assert combined_first_only == old_first_only
        embedding += 1
    old_logs = np.asarray((FC_INITIAL[0], FC_INITIAL[1], FC_INITIAL[2], FC_INITIAL[3], FC_INITIAL[5], FC_INITIAL[6]))
    old_bound = nested_selector_bound(FS_PROFILE, None, old_logs)
    combined_bound = forced_run_combined_bound()
    assert abs(old_bound['raw_log2_probability_bound']-combined_bound['raw_log2_probability_bound']) < 1e-7
    assert abs(old_bound['log2_lambda']-combined_bound['log2_lambda']) < 1e-12
    embedding += 1
    second_only = np.asarray(FS_INITIAL[:6]+(0.0, FS_INITIAL[6]))
    original = forced_run_scalar_bound()
    retained = forced_run_combined_bound(logs=second_only)
    assert abs(original['raw_log2_probability_bound']-retained['raw_log2_probability_bound']) < 1e-7
    embedding += 1
    combined_point = np.asarray(FC_INITIAL)
    combined_data = _fs_profile(FS_PROFILE, None, True)
    _, gradient, hessian, _ = _fs_objective(combined_point, 1/16, combined_data, True)
    for coordinate in range(8):
        delta = np.zeros(8)
        delta[coordinate] = 1e-5
        plus = _fs_objective(combined_point+delta, 1/16, combined_data, True)
        minus = _fs_objective(combined_point-delta, 1/16, combined_data, True)
        assert abs((plus[0]-minus[0])/(2e-5)-gradient[coordinate]) < 1e-6
        assert np.max(np.abs((plus[1]-minus[1])/(2e-5)-hessian[:, coordinate])) < 1e-5
    print(dict(exact_local_compositions=local, exact_width8_partition_controls=physical,
               exact_free_endpoint_cases=endpoint, exact_whole_run_history_checks=histories,
               exact_forced_run_packing_checks=packing, gradient_hessian_coordinates=7,
               exact_combined_local_compositions=combined_local, first_second_cap_embeddings=embedding,
               combined_gradient_hessian_coordinates=8,
               full_code_certificate=False, **_fs_compulsory_hit_checks()), flush=True)


def forced_run_scalar_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    mode.add_argument('--search', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=FS_PROFILE)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--combined', action='store_true')
    parser.add_argument('--compulsory-hit', action='store_true')
    parser.add_argument('--logs', type=float, nargs='+')
    parser.add_argument('--evaluations', type=int, default=1024)
    args = parser.parse_args()
    logs = args.logs if args.logs is not None else (FC_INITIAL if args.combined else FS_INITIAL)
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
            forced_run_scalar_checks()
        elif args.search:
            search = forced_run_combined_search if args.combined else forced_run_scalar_search
            print(search(args.profile, args.death_upper, logs, max_evaluations=args.evaluations,
                         compulsory_hit=args.compulsory_hit), flush=True)
        else:
            bound = forced_run_combined_bound if args.combined else forced_run_scalar_bound
            print(bound(args.profile, args.death_upper, logs, compulsory_hit=args.compulsory_hit), flush=True)
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    forced_run_scalar_main()
