r"""Bounded FLOAT global count of saturated support histories.

Saturate each occupied allowed coordinate run and every open allowed run.
The audited map preserves canonical U, births B, footprint/rank/eligible
row class, raises support volume L, and lowers deaths D. Thus it suffices
to count saturated histories with L>=R+1,D<=W-R-1 for the dimension
obstruction class, at fixed(h,R,f,F). No amplitude count is introduced.

Labels are Z outside U, O=U\Aprev, E=unused Aprev, S=used Aprev. For
counts(o,e,s), r=o+e+s<=4, the potential alpha^u beta^a gamma^s has
unselected action d^s and positive selected action

 z^4 alpha^4 beta^r * sum_{k=0..o,b=0..4-r,s+k+b>=1}
     C(o,k)C(4-r,b) (gamma*v)^(s+k) (w*beta*gamma*v)^b.

Old used runs continue used; old unused runs continue unused. Newly
starting unforced runs may be used or unused, while births are used.
The hit condition s+k+b>=1 is retained on every genuine selected quartet.
At the formal endpoint E must be empty and all remaining open positions
are filled. The endpoint ratio is bounded by
 max(1+w*v,v/alpha,v/(alpha*beta*gamma))^512.

The route factor already includes the H union and requires B subset H.
For v>=1,0<d<=1 the external correction is
 logQ-Glogz-flogw-L0logv-D0logd.
No field factor, birth binomial, or separate support union is added.
The35-case pointwise comparison is a global positive-count control, not
a local-state transition DP. Its finite-box epigraph search is FLOAT;
upper/lower diagnostics are not exact probability/optimality certificates.

Optional --averaged-cone uses the reachable inequality3*S>=O+E. The
quartic IID-label polynomial on that cone has35 positive Bernstein
coefficients at vertices Z,S,(3O+S)/4,(3E+S)/4. They mix the normalized
local actions above. Exact partition conditioning costs at most
C_mode=4^512*(128!)^4/512! per genuine update. The total extra correction
is4095*log2(C_mode), not128 times this quantity. The pointwise and averaged
objectives are optimized separately; their minimum is not treated as convex.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, isfinite, lgamma, log, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np

try:
    _ns_logsum
except NameError:
    from quartet_nested_selector_scalar import _ns_logsum, _ns_smooth_max
try:
    _fb_barrier
except NameError:
    from quartet_forced_run_barrier import (
        _fb_barrier, _fb_dual, _fb_newton_direction, _fb_raw,
        FB_RATIO, FB_LOCAL_COUNT, FB_ENDPOINT_COUNT)


SAT_WIDTH, SAT_UPDATES, SAT_M, SAT_W = 512, 4096, 8192, 943718
SAT_N, SAT_POWER = SAT_WIDTH*SAT_UPDATES, 128*(SAT_UPDATES-1)
SAT_RATIO = SAT_WIDTH/SAT_POWER
SAT_PROFILE = (4096, 819200, 65536, 1081344)
SAT_NAMES = ('log_z', 'log_w', 'log_alpha', 'log_beta', 'log_gamma', 'log_v', 'log_d')
SAT_INITIAL = (-1.0, -1.0, 0.0, 0.0, 0.0, 1.0, -1.0)
SAT_BOUNDS = ((-8., 4.), (-20., 8.), (-8., 8.), (-12., 8.), (-12., 12.), (0., 12.), (-16., 0.))


def _sat_cases():
    cases = []
    for o in range(5):
        for e in range(5-o):
            for s in range(5-o-e):
                r = o+e+s
                exponents = [(0, 0, -r, -e-s, -s, 0, s)]
                coefficients = [0.0]
                for k in range(o+1):
                    for b in range(5-r):
                        if s+k+b == 0:
                            continue
                        exponents.append((4, b, 4-r, o+b, k+b, s+k+b, 0))
                        coefficients.append(log2(comb(o, k)*comb(4-r, b)))
                cases.append((o, e, s, np.asarray(exponents, dtype=float), np.asarray(coefficients, dtype=float)))
    return tuple(cases)


_SAT_CASES = _sat_cases()


def _sat_mode_factor(width):
    if width % 4:
        raise ValueError('Balanced four-label mode factor requires width divisible4')
    quarter = width//4
    multinomial = comb(width, quarter)*comb(width-quarter, quarter)*comb(width-2*quarter, quarter)
    return Fraction(4**width, multinomial)


_SAT_MODE = _sat_mode_factor(SAT_WIDTH)
SAT_LOG_MODE = log2(_SAT_MODE.numerator)-log2(_SAT_MODE.denominator)


def _sat_bernstein():
    lookup = {(o, e, s): i for i, (o, e, s, _, _) in enumerate(_SAT_CASES)}
    result = []
    for kO in range(5):
        for kE in range(5-kO):
            for kS in range(5-kO-kE):
                kZ = 4-kO-kE-kS
                indices, weights = [], []
                for o in range(kO+1):
                    for e in range(kE+1):
                        s = kS+kO+kE-o-e
                        indices.append(lookup[o, e, s])
                        weights.append(Fraction(comb(kO, o)*comb(kE, e)*3**(o+e), 4**(kO+kE)))
                result.append(((kZ, kS, kO, kE), np.asarray(indices), tuple(weights),
                               np.asarray([log2(float(weight)) for weight in weights])))
    return tuple(result)


_SAT_BERNSTEIN = _sat_bernstein()


def saturated_actions(logs, averaged_cone=False):
    rows = [_ns_logsum(exponents @ logs+coefficients, exponents)
            for _, _, _, exponents, coefficients in _SAT_CASES]
    if averaged_cone:
        values = np.asarray([row[0] for row in rows])
        gradients = np.asarray([row[1] for row in rows])
        hessians = np.asarray([row[2] for row in rows])
        rows = [_ns_smooth_max(values[indices]+logweights, gradients[indices], hessians[indices], 1.0)
                for _, indices, _, logweights in _SAT_BERNSTEIN]
    first_exp = np.asarray(((0, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 1, 0)), dtype=float)
    endpoint = [_ns_logsum(first_exp @ logs, first_exp)]
    for exponent in ((0, 0, -1, 0, 0, 1, 0), (0, 0, -1, -1, -1, 1, 0)):
        vector = np.asarray(exponent, dtype=float)
        endpoint.append((float(vector @ logs), vector, np.zeros((7, 7))))
    return rows, endpoint


def saturated_model(logs, averaged_cone=False):
    local, endpoint = saturated_actions(logs, averaged_cone)
    # The existing generic barrier expects four endpoint constraints. A
    # duplicate changes neither the maximum nor the positive dual argument.
    rows = local+endpoint+[endpoint[0]]
    return (np.asarray([row[0] for row in rows]), np.asarray([row[1] for row in rows]),
            np.asarray([row[2] for row in rows]))


def _sat_logchoose(n, k):
    if not 0 <= k <= n:
        return float('-inf')
    if k in (0, n):
        return 0.0
    return (lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/log(2)


def saturated_profile_data(profile=SAT_PROFILE, volume_lower=None, death_upper=None):
    h, rank, births, footprint = map(index, profile)
    if not (1 <= h <= SAT_M and 128*h <= rank <= 1 << 20 and 1 <= births <= footprint <= SAT_N):
        raise ValueError('Invalid fixed(h,R,f,F) profile')
    G, CH = footprint-births, footprint-rank+128*h
    if G % 4 or G > SAT_N-SAT_WIDTH:
        raise ValueError('Invalid canonical live volume')
    L0 = rank+1 if volume_lower is None else index(volume_lower)
    D0 = SAT_W-rank-1 if death_upper is None else index(death_upper)
    if L0 < 0:
        raise ValueError('Volume threshold must be nonnegative')
    empty = (L0 > footprint or D0 < max(0, births-SAT_WIDTH)
             or not (129*h <= CH <= 256*h and 0 <= CH-births <= G
                     and 0 <= 256*h-CH <= SAT_N-footprint))
    route = float('-inf') if empty else min(0.0, _sat_logchoose(G, CH-births)
        +_sat_logchoose(SAT_N-footprint, 256*h-CH)-255*_sat_logchoose(SAT_M, h))
    linear = np.asarray((-G, -births, 0., 0., 0., -L0, -D0))/SAT_POWER
    return dict(profile=(h, rank, births, footprint), G=G, C_H=CH, L0=L0, D0=D0,
                empty=empty, linear=linear, route_log2=route,
                default_dimension_class=(volume_lower is None and death_upper is None))


def saturated_cover_bound(profile=SAT_PROFILE, volume_lower=None, death_upper=None, logs=SAT_INITIAL,
                           averaged_cone=False):
    logs = np.asarray(logs, dtype=float)
    if logs.shape != (7,) or not np.all(np.isfinite(logs)) or logs[5] < 0 or logs[6] > 0:
        raise ValueError('Require seven finite logs with log_v>=0 and log_d<=0')
    data = saturated_profile_data(profile, volume_lower, death_upper)
    if data['empty']:
        return dict(profile=data['profile'], status='empty_threshold_class', full_code_certificate=False)
    local, endpoint = saturated_actions(logs, averaged_cone)
    values = [entry[0] for entry in local]
    growth, last = max(values), max(entry[0] for entry in endpoint)
    mode_correction = (SAT_UPDATES-1)*SAT_LOG_MODE if averaged_cone else 0.0
    correction = data['route_log2']+SAT_POWER*float(data['linear'] @ logs)+mode_correction
    order = sorted(range(35), key=values.__getitem__, reverse=True)
    if averaged_cone:
        tight = [dict(vertex_counts=dict(zip(('Z', 'S', 'O_mix', 'E_mix'), _SAT_BERNSTEIN[i][0])),
                      log2_ratio=values[i]) for i in order[:6]]
    else:
        tight = [dict(o=_SAT_CASES[i][0], e=_SAT_CASES[i][1], s=_SAT_CASES[i][2],
                      log2_ratio=values[i]) for i in order[:6]]
    return dict(profile=data['profile'], volume_lower=data['L0'], death_upper=data['D0'],
                status='saturated_support_global_count', logs=dict(zip(SAT_NAMES, map(float, logs))),
                log2_local_growth=growth, log2_endpoint_per_coordinate=last,
                external_log2_correction=correction,
                raw_log2_probability_bound=SAT_POWER*growth+SAT_WIDTH*last+correction,
                tightest_compositions=tight, averaged_reachable_cone=bool(averaged_cone),
                log2_mode_factor_per_update=SAT_LOG_MODE if averaged_cone else 0.0,
                total_log2_mode_correction=mode_correction,
                default_dimension_class=data['default_dimension_class'],
                covers_original_supports_by_saturation=True, saturated_volume_not_fixed_exactly=True,
                all_open_runs_filled=True, unused_active_endpoint_rejected=True,
                birth_set_and_normalization_unchanged=True, H_union_already_in_route=True,
                exact_actual_hits_retained=True, global_cone_multiplier_used=False,
                field_factor_used=False, floating_diagnostic_only=True,
                optimizer_optimality_proved=False, full_code_certificate=False)


def saturated_cover_search(profile=SAT_PROFILE, volume_lower=None, death_upper=None,
                            initial=SAT_INITIAL, bounds=SAT_BOUNDS, max_evaluations=1024, mu_initial=.01,
                            averaged_cone=False):
    """Bounded seven-log epigraph search, with a finite-box FLOAT dual diagnostic."""
    if (FB_LOCAL_COUNT, FB_ENDPOINT_COUNT) != (35, 4) or abs(FB_RATIO-SAT_RATIO) > 1e-15:
        raise ValueError('Incompatible read-only generic barrier primitives')
    initial, bounds = np.asarray(initial, dtype=float), np.asarray(bounds, dtype=float)
    max_evaluations = index(max_evaluations)
    if (initial.shape != (7,) or bounds.shape != (7, 2) or not np.all(np.isfinite(initial))
            or not np.all(np.isfinite(bounds)) or np.any(bounds[:, 0] >= bounds[:, 1])
            or bounds[5, 0] < 0 or bounds[6, 1] > 0 or not 64 <= max_evaluations <= 2048
            or not isfinite(mu_initial) or not 1e-6 <= mu_initial <= .1):
        raise ValueError('Invalid v/d marker box or bounded evaluation schedule')
    data = saturated_profile_data(profile, volume_lower, death_upper)
    if data['empty']:
        return dict(profile=data['profile'], status='empty_threshold_class', full_code_certificate=False)
    if averaged_cone:
        # The generic barrier's constant slot is named route_log2. Pass an
        # explicit copy including the mode correction; no globals are changed.
        data = dict(data)
        data['route_log2'] += (SAT_UPDATES-1)*SAT_LOG_MODE
    baseline = np.clip(initial, bounds[:, 0], bounds[:, 1])
    baseline_model = saturated_model(baseline, averaged_cone)
    best_logs, best_value = baseline.copy(), _fb_raw(baseline, baseline_model, data)
    best_lower, best_dual = float('-inf'), None
    margin = .001*(bounds[:, 1]-bounds[:, 0])
    logs = np.clip(baseline, bounds[:, 0]+margin, bounds[:, 1]-margin)
    model = saturated_model(logs, averaged_cone)
    point = np.concatenate((logs, (model[0][:35].max()+35*mu_initial,
                                   model[0][35:].max()+4*mu_initial/SAT_RATIO)))
    evaluations = iterations = stages = backtracks = 0
    termination = 'schedule_complete'

    def evaluate(candidate, mu):
        nonlocal evaluations, best_logs, best_value, best_lower, best_dual
        if evaluations >= max_evaluations:
            return None
        evaluations += 1
        if (not np.all(np.isfinite(candidate)) or np.any(candidate[:7] <= bounds[:, 0])
                or np.any(candidate[:7] >= bounds[:, 1])):
            return None
        result = _fb_barrier(candidate, mu, bounds, data, saturated_model(candidate[:7], averaged_cone))
        if result is not None:
            if result['raw'] < best_value:
                best_logs, best_value = candidate[:7].copy(), result['raw']
            lower, diagnostic = _fb_dual(candidate, bounds, data, result, SAT_NAMES)
            if lower > best_lower:
                best_lower, best_dual = lower, diagnostic
        return result

    for stage in range(14):
        mu, centered = mu_initial*(.2**stage), False
        for _ in range(64):
            current = evaluate(point, mu)
            if current is None:
                if evaluations >= max_evaluations:
                    termination = 'evaluation_budget'; break
                raise ArithmeticError('Accepted saturated-cover epigraph point lost feasibility')
            direction, _ = _fb_newton_direction(current['gradient'], current['hessian'])
            descent = float(current['gradient'] @ direction)
            if -descent <= max(1e-16, mu*1e-5):
                centered = True; break
            if not isfinite(descent) or descent >= 0:
                termination = 'numerical_direction_stall'; break
            step = 1.
            negative, positive = direction[:7] < 0, direction[:7] > 0
            if np.any(negative):
                step = min(step, float(np.min(-.99*current['lower'][negative]/direction[:7][negative])))
            if np.any(positive):
                step = min(step, float(np.min(.99*current['upper'][positive]/direction[:7][positive])))
            change = current['jacobian'] @ direction
            toward = change > 0
            if np.any(toward):
                step = min(step, float(np.min(.99*current['slack'][toward]/change[toward])))
            accepted = False
            for _ in range(40):
                trial_point = point+step*direction
                trial = evaluate(trial_point, mu)
                if trial is not None and trial['value'] <= current['value']+1e-4*step*descent:
                    point, accepted = trial_point, True; break
                if evaluations >= max_evaluations:
                    termination = 'evaluation_budget'; break
                step *= .5
                backtracks += 1
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
    if best_dual is None:
        raise ArithmeticError('No finite saturated-cover dual diagnostic')
    result = saturated_cover_bound(profile, volume_lower, death_upper, best_logs, averaged_cone)
    result.update(optimizer_evaluations=evaluations, optimizer_iterations=iterations,
                  optimizer_centered_stages=stages, optimizer_backtracks=backtracks,
                  optimizer_termination=termination, search_box=bounds.tolist(),
                  finite_box_dual_lower_diagnostic=SAT_POWER*best_lower,
                  finite_box_primal_dual_gap_diagnostic=result['raw_log2_probability_bound']-SAT_POWER*best_lower,
                  best_dual_diagnostic=best_dual, lower_diagnostic_is_exact_certificate=False,
                  scalar_family_impossibility_proved=False, unbounded_parameter_domain_certified=False)
    return result


def _sat_subsets(mask):
    subset = mask
    while True:
        yield subset
        if subset == 0:
            break
        subset = (subset-1)&mask


def _sat_canonical(support, partitions, width):
    live, births, previous = [], [], 0
    for state, partition in zip(support, partitions):
        union = sum(block for block in partition if block&previous)
        live.append(union)
        births.append(state&~union)
        previous = state
    allowed = [(live[t]&live[t+1])|births[t] for t in range(len(support)-1)]
    allowed.append(live[-1]|births[-1])
    return tuple(live), tuple(births), tuple(allowed)


def saturate_support_history(support, partitions, width):
    """Independent geometric map: fill occupied runs and every open allowed run."""
    support = tuple(support)
    live, births, allowed = _sat_canonical(support, partitions, width)
    saturated = [0]*len(support)
    for coordinate in range(width):
        bit, start = 1 << coordinate, None
        for t in range(len(support)+1):
            present = t < len(support) and bool(allowed[t]&bit)
            if present and start is None:
                start = t
            if not present and start is not None:
                if t == len(support) or any(support[j]&bit for j in range(start, t)):
                    for j in range(start, t):
                        saturated[j] |= bit
                start = None
    return tuple(saturated), live, births, allowed


def _sat_deaths(support):
    return sum((support[t]&~support[t+1]).bit_count() for t in range(len(support)-1))


def _sat_transitions(old, active, used, partition, z, w, v, d):
    for choice in range(1 << len(partition)):
        blocks = [block for i, block in enumerate(partition) if choice >> i & 1]
        new = sum(blocks)
        for births in _sat_subsets(new&~old):
            for starting_used in _sat_subsets(new&old&~active):
                new_used = (used&new)|starting_used|births
                if any(not block&new_used for block in blocks):
                    continue
                new_active = (old&new)|births
                weight = (z**new.bit_count()*w**births.bit_count()*v**new_used.bit_count()
                          *d**(used&~new).bit_count())
                yield new, new_active, new_used, weight


def _sat_multinomial(counts):
    remaining, value = sum(counts), 1
    for count in counts:
        value *= comb(remaining, count)
        remaining -= count
    return value


def _sat_quartic_value(probabilities, actions):
    total = Fraction(0)
    for i, (o, e, s, _, _) in enumerate(_SAT_CASES):
        counts = (4-o-e-s, o, e, s)
        weight = Fraction(_sat_multinomial(counts))
        for probability, count in zip(probabilities, counts):
            weight *= probability**count
        total += weight*actions[i]
    return total


def _sat_averaged_checks(physical_actions, logs):
    weight_checks = identities = physical_rows = mode_checks = derivative_checks = 0
    for _, _, weights, _ in _SAT_BERNSTEIN:
        assert sum(weights) == 1
        weight_checks += 1
    points = ((Fraction(1, 4),)*4,
              (Fraction(1, 2), Fraction(0), Fraction(1, 4), Fraction(1, 4)),
              (Fraction(0), Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)))
    points += tuple(tuple(Fraction(int(i == j)) for j in range(4)) for i in range(4))
    for actions in (physical_actions, [Fraction(i+1, i%3+1) for i in range(35)]):
        mixed = [sum(weight*actions[int(index)] for index, weight in zip(indices, weights))
                 for _, indices, weights, _ in _SAT_BERNSTEIN]
        for qZ, qS, qO, qE in points:
            probabilities = (qZ, 3*qO/4, 3*qE/4, qS+(qO+qE)/4)
            direct = _sat_quartic_value(probabilities, actions)
            bernstein = Fraction(0)
            for (counts, _, _, _), coefficient in zip(_SAT_BERNSTEIN, mixed):
                weight = Fraction(_sat_multinomial(counts))
                for probability, count in zip((qZ, qS, qO, qE), counts):
                    weight *= probability**count
                bernstein += weight*coefficient
            assert direct == bernstein <= max(mixed)
            identities += 1
    physical_mixed = [sum(weight*physical_actions[int(index)] for index, weight in zip(indices, weights))
                      for _, indices, weights, _ in _SAT_BERNSTEIN]
    lookup = {(o, e, s): i for i, (o, e, s, _, _) in enumerate(_SAT_CASES)}
    local_polynomial = {(o, e, s): Fraction(_sat_multinomial((4-o-e-s, o, e, s)))*physical_actions[i]
                        for i, (o, e, s, _, _) in enumerate(_SAT_CASES)}
    for width in (4, 8):
        mode_factor = _sat_mode_factor(width)
        least_probability = Fraction(1)
        for o in range(width+1):
            for e in range(width-o+1):
                for s in range(width-o-e+1):
                    counts = (width-o-e-s, o, e, s)
                    probability = Fraction(_sat_multinomial(counts))
                    for count in counts:
                        probability *= Fraction(count, width)**count
                    assert probability >= 1/mode_factor
                    least_probability = min(least_probability, probability)
                    mode_checks += 1
        assert least_probability == 1/mode_factor
        coefficient = {(0, 0, 0): Fraction(1)}
        for _ in range(width//4):
            following = {}
            for left, value in coefficient.items():
                for right, term in local_polynomial.items():
                    key = tuple(left[i]+right[i] for i in range(3))
                    following[key] = following.get(key, Fraction(0))+value*term
            coefficient = following
        partitions = [(15,)] if width == 4 else [
            (block, 255^block) for rest in combinations(range(1, 8), 3)
            for block in (1+sum(1 << i for i in rest),)]
        for u in range(0, width+1, 4):
            for o in range(u+1):
                for e in range(u-o+1):
                    s = u-o-e
                    if 4*s < u:
                        continue
                    old, active, used = (1 << u)-1, ((1 << (e+s))-1) << o, ((1 << s)-1) << (o+e)
                    physical = Fraction(0)
                    for partition in partitions:
                        action = Fraction(1)
                        for block in partition:
                            triple = ((block&old&~active).bit_count(), (block&active&~used).bit_count(),
                                      (block&used).bit_count())
                            action *= physical_actions[lookup[triple]]
                        physical += action
                    physical /= len(partitions)
                    counts = (width-u, o, e, s)
                    multinomial = _sat_multinomial(counts)
                    assert physical == coefficient[o, e, s]/multinomial
                    probabilities = tuple(Fraction(count, width) for count in counts)
                    mode_probability = Fraction(multinomial)
                    for probability, count in zip(probabilities, counts):
                        mode_probability *= probability**count
                    quartic = _sat_quartic_value(probabilities, physical_actions)
                    assert physical <= quartic**(width//4)/mode_probability
                    assert physical <= mode_factor*max(physical_mixed)**(width//4)
                    physical_rows += 1
    model = saturated_model(logs, True)
    for coordinate in range(7):
        delta = np.zeros(7); delta[coordinate] = 1e-5
        plus, minus = saturated_model(logs+delta, True), saturated_model(logs-delta, True)
        assert np.max(np.abs((plus[0]-minus[0])/(2e-5)-model[1][:, coordinate])) < 1e-6
        assert np.max(np.abs((plus[1]-minus[1])/(2e-5)-model[2][:, :, coordinate])) < 1e-5
        derivative_checks += 1
    data, bounds = saturated_profile_data(), np.asarray(SAT_BOUNDS)
    data = dict(data)
    data['route_log2'] += (SAT_UPDATES-1)*SAT_LOG_MODE
    point = np.concatenate((logs, (model[0][:35].max()+1, model[0][35:].max()+1)))
    evaluation = _fb_barrier(point, .001, bounds, data, model)
    lower, _ = _fb_dual(point, bounds, data, evaluation, SAT_NAMES)
    for candidate in (logs, np.asarray(SAT_INITIAL), bounds[:, 0], bounds[:, 1]):
        assert lower <= _fb_raw(candidate, saturated_model(candidate, True), data)+1e-10
    return dict(exact_bernstein_weight_sums=weight_checks, exact_quartic_cone_identities=identities,
                exact_physical_cone_rows=physical_rows, exact_balanced_mode_probabilities=mode_checks,
                averaged_gradient_hessian_coordinates=derivative_checks, averaged_dual_checks=4,
                mode_factor_applied_once_per_update=True)


def saturated_cover_checks():
    histories = strict_volume = strict_deaths = local = endpoint = moments = derivatives = duals = 0
    z, w, v, d = Fraction(1, 2), Fraction(1, 3), Fraction(3, 2), Fraction(1, 2)
    alpha, beta, gamma = Fraction(2), Fraction(3, 4), Fraction(5, 4)
    exact_saturated_moments = {}

    def check_history(support, partitions, width):
        nonlocal histories, strict_volume, strict_deaths
        saturated, U, B, A = saturate_support_history(support, partitions, width)
        new_U, new_B, new_A = _sat_canonical(saturated, partitions, width)
        assert U == new_U and B == new_B and A == new_A
        assert all(not support[t]&~saturated[t] for t in range(len(support)))
        old_L, new_L = sum(map(int.bit_count, support)), sum(map(int.bit_count, saturated))
        old_D, new_D = _sat_deaths(support), _sat_deaths(saturated)
        assert new_L >= old_L and new_D <= old_D
        assert saturate_support_history(saturated, partitions, width)[0] == saturated
        assert tuple(U[t]|B[t] for t in range(len(support))) == tuple(new_U[t]|new_B[t] for t in range(len(support)))
        histories += 1
        strict_volume += new_L > old_L
        strict_deaths += new_D < old_D
        return saturated, U, B

    for updates in range(1, 4):
        partitions, total = ((15,),)*updates, Fraction(0)
        for support in product(range(16), repeat=updates):
            saturated, U, B = check_history(support, partitions, 4)
            if saturated == support:
                total += z**sum(map(int.bit_count, U))*w**sum(map(int.bit_count, B))*v**sum(map(int.bit_count, support))*d**_sat_deaths(support)
        exact_saturated_moments[updates] = total
    P0, P1, P2 = (15, 240), (23, 232), (113, 142)
    for partitions in ((P0, P1, P2, P0), (P2, P1, P0, P2), (P1, P0, P2, P1)):
        for support in product((0, 1, 3, 15, 113, 255), repeat=4):
            check_history(support, partitions, 8)
    strict_support = (1, 3, 1, 3)
    strict_saturated, strict_U, strict_B = check_history(strict_support, ((15,),)*4, 4)
    assert strict_saturated == (1, 15, 15, 15)
    assert sum(map(int.bit_count, strict_support)) == 6
    assert sum(map(int.bit_count, strict_saturated)) == 13
    assert _sat_deaths(strict_support) == 1 and _sat_deaths(strict_saturated) == 0
    assert strict_U == (0, 15, 15, 15) and strict_B == (1, 0, 0, 0)
    assert strict_volume > 0 and strict_deaths > 0
    ratios = []
    for o, e, s, exponents, _ in _SAT_CASES:
        r = o+e+s
        old, active, used = (1 << r)-1, ((1 << (e+s))-1) << o, ((1 << s)-1) << (o+e)
        physical = sum(weight*alpha**new.bit_count()*beta**new_active.bit_count()*gamma**new_used.bit_count()
                       for new, new_active, new_used, weight in _sat_transitions(old, active, used, (15,), z, w, v, d))
        selected = sum(Fraction(comb(o, k)*comb(4-r, b))*(gamma*v)**(s+k)*(w*beta*gamma*v)**b
                       for k in range(o+1) for b in range(5-r) if s+k+b >= 1)
        assert physical == d**s+z**4*alpha**4*beta**r*selected
        multiplicities = [1]+[comb(o, k)*comb(4-r, b) for k in range(o+1) for b in range(5-r) if s+k+b >= 1]
        model = Fraction(0)
        for coefficient, powers in zip(multiplicities, exponents):
            term = Fraction(coefficient)
            for marker, power in zip((z, w, alpha, beta, gamma, v, d), powers):
                term *= marker**int(power)
            model += term
        assert model == physical/(alpha**r*beta**(e+s)*gamma**s)
        ratios.append(model)
        local += 1
        terminal = Fraction(0) if e else v**r*(1+w*v)**(4-r)
        endpoint_bound = max(1+w*v, v/alpha, v/(alpha*beta*gamma))**4*alpha**r*beta**(e+s)*gamma**s
        assert terminal <= endpoint_bound
        endpoint += 1
    states = {(0, 0, 0): Fraction(1)}
    for updates in range(1, 4):
        count = sum(weight*v**old.bit_count()*(1+w*v)**(4-old.bit_count())
                    for (old, active, used), weight in states.items() if active == used)
        assert count == exact_saturated_moments[updates]
        assert count <= max(ratios)**(updates-1)*max(1+w*v, v/alpha, v/(alpha*beta*gamma))**4
        moments += 1
        following = {}
        for (old, active, used), weight in states.items():
            for new, new_active, new_used, step in _sat_transitions(old, active, used, (15,), z, w, v, d):
                key = (new, new_active, new_used)
                following[key] = following.get(key, Fraction(0))+weight*step
        states = following
    logs = np.asarray(tuple(log2(float(value)) for value in (z, w, alpha, beta, gamma, v, d)))
    model = saturated_model(logs)
    for coordinate in range(7):
        delta = np.zeros(7); delta[coordinate] = 1e-5
        plus, minus = saturated_model(logs+delta), saturated_model(logs-delta)
        assert np.max(np.abs((plus[0]-minus[0])/(2e-5)-model[1][:, coordinate])) < 1e-6
        assert np.max(np.abs((plus[1]-minus[1])/(2e-5)-model[2][:, :, coordinate])) < 1e-5
        derivatives += 1
    data, bounds = saturated_profile_data(), np.asarray(SAT_BOUNDS)
    point = np.concatenate((logs, (model[0][:35].max()+1, model[0][35:].max()+1)))
    evaluation = _fb_barrier(point, .001, bounds, data, model)
    lower, _ = _fb_dual(point, bounds, data, evaluation, SAT_NAMES)
    for candidate in (logs, np.asarray(SAT_INITIAL), bounds[:, 0], bounds[:, 1]):
        assert lower <= _fb_raw(candidate, saturated_model(candidate), data)+1e-10
        duals += 1
    print(dict(independent_saturation_histories=histories, strict_volume_increases=strict_volume,
               strict_death_decreases=strict_deaths, exact_local_compositions=local,
               exact_free_endpoint_cases=endpoint, exact_saturated_history_moments=moments,
               gradient_hessian_coordinates=derivatives, finite_box_dual_checks=duals,
               footprint_and_births_preserved=True, saturation_idempotent=True,
               full_code_certificate=False, **_sat_averaged_checks(ratios, logs)), flush=True)


def saturated_cover_main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--checks', action='store_true')
    modes.add_argument('--probe', action='store_true')
    modes.add_argument('--search', action='store_true')
    parser.add_argument('--profile', type=int, nargs=4, default=SAT_PROFILE)
    parser.add_argument('--volume-lower', type=int)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--averaged-cone', action='store_true')
    parser.add_argument('--logs', type=float, nargs=7, default=SAT_INITIAL)
    parser.add_argument('--evaluations', type=int, default=1024)
    parser.add_argument('--mu', type=float, default=.01)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124)); timer.daemon = True; timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            saturated_cover_checks()
        elif args.search:
            print(saturated_cover_search(args.profile, args.volume_lower, args.death_upper, args.logs,
                                         max_evaluations=args.evaluations, mu_initial=args.mu,
                                         averaged_cone=args.averaged_cone), flush=True)
        else:
            print(saturated_cover_bound(args.profile, args.volume_lower, args.death_upper, args.logs,
                                        averaged_cone=args.averaged_cone), flush=True)
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    saturated_cover_main()
