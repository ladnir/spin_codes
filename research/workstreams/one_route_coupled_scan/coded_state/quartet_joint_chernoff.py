"""Full covered-hole bins bounded by positive MDS exponential actions.

Inject quartet_joint_mds.py and quartet_joint_subset.py definitions, with
their main blocks removed.  No existing helper is modified or imported.
The lift is phi(u,e)=2**(logpotential*u+log_hole_potential*e) for u>0,
and phi(0,e)=1, on disjoint fixed bins.  Thus the general zero operator
intertwines exactly.  Zero-u targets use only exact zero-u axis actions.
Destination Chernoff tilts are minimized per physical source BEFORE taking
source-bin maxima.  Quota minima occur only after regional composition.
All reported numerical probes are FLOAT diagnostics, not certificates.
"""

import argparse
from fractions import Fraction
from itertools import product
from math import comb, lgamma
import os
import shutil
import threading
import time

import numpy as np


JC_LD = np.longdouble
JC_W = 943718


def _joint_chernoff_dependencies():
    required = ('joint_mds_local', 'joint_subset_region', 'joint_subset_tail')
    if any(name not in globals() for name in required):
        raise RuntimeError('Inject quartet_joint_mds.py and quartet_joint_subset.py before this helper')


def _jc_cast(value):
    if isinstance(value, Fraction):
        return JC_LD(value.numerator) / JC_LD(value.denominator)
    return JC_LD(value)


def joint_chernoff_prepare(width, u_step=32, e_step=32, logpotential=0,
                           tilt_logs=(-1/16, 0, 1/16), log_hole_potential=0):
    if width < 4 or width % 4 or min(u_step, e_step) < 1:
        raise ValueError('Require a quartet width and positive bin widths')
    if not np.isfinite(logpotential) or not np.isfinite(log_hole_potential):
        raise ValueError('The fixed log2 potentials must be finite')
    states = np.array([(u, e) for u in range(width + 1) for e in range(width - u + 1)], dtype=int)
    raw_keys = [(0 if u == 0 else 1 + (u - 1) // u_step,
                 0 if e == 0 else 1 + (e - 1) // e_step) for u, e in states]
    keys = tuple(sorted(set(raw_keys)))
    index = {key: i for i, key in enumerate(keys)}
    bins = np.array([index[key] for key in raw_keys], dtype=int)
    order = np.argsort(bins, kind='stable')
    starts = np.r_[0, np.flatnonzero(np.diff(bins[order])) + 1]
    tilts = tuple(product(tuple(tilt_logs), repeat=2))
    if not tilts or any(not np.isfinite(a + b) for a, b in tilts):
        raise ValueError('Require a finite nonempty Chernoff tilt family')
    # Axes are exact coefficient slices, used only by matching destination
    # bins.  They do not replace or change the finite interior tilt family.
    candidates = (tuple(('interior', lt, ls) for lt, ls in tilts)
                  + tuple(('e_zero', lt, 0) for lt in dict.fromkeys(lt for lt, _ in tilts))
                  + tuple(('u_zero', 0, ls) for ls in dict.fromkeys(ls for _, ls in tilts)))
    minima = np.empty((len(keys), len(candidates)), dtype=JC_LD)
    allowed = np.empty((len(keys), len(candidates)), dtype=bool)
    for j, (kind, lt, ls) in enumerate(candidates):
        costs = JC_LD(lt) * states[:, 0] + JC_LD(ls) * states[:, 1]
        minima[:, j] = np.minimum.reduceat(costs[order], starts)
        # On u=0 the lift is one, so the separable base s0**e does not
        # represent that target potential.  Use only its exact axis family.
        allowed[:, j] = [(kind == 'u_zero') if key[0] == 0 else
                         (kind == 'interior' or (kind == 'e_zero' and key[1] == 0))
                         for key in keys]
    max_sum = np.maximum.reduceat(states.sum(axis=1)[order], starts)
    logphi = (JC_LD(logpotential) * states[:, 0]
              + JC_LD(log_hole_potential) * states[:, 1] * (states[:, 0] > 0))
    zero_vector = np.array([JC_LD(key[0] == 0) for key in keys])
    return dict(width=width, u_step=u_step, e_step=e_step, logpotential=JC_LD(logpotential),
                log_hole_potential=JC_LD(log_hole_potential),
                states=states, keys=keys, bins=bins, order=order, starts=starts,
                tilts=tilts, candidates=candidates, candidate_allowed=allowed,
                minima=minima, max_sum=max_sum, logphi=logphi,
                zero_vector=zero_vector)


def joint_chernoff_positive_action(width, m, b, eligible, outside, alpha, t, s,
                                   coverage=1, output_discount=1, exact=False):
    """Exact positive P(phi), P=K-Z, with no subtraction of the zero path.

    The k>=1 hit-quartet terms have old u>0.  For old u=0, a positive
    recurrence computes B^groups-B0^groups, B0=(1+z)^4.  Quota-uncovered
    labels have already been summed to one; formal conditioning is applied
    to the resulting positive regional operator, as in the parent helpers.
    """
    if width < 4 or width % 4:
        raise ValueError('Require a positive quartet width')
    cast = Fraction if exact else _jc_cast
    dtype = object if exact else JC_LD
    m, b, eligible, outside, alpha, t, s, coverage, output_discount = map(
        cast, (m, b, eligible, outside, alpha, t, s, coverage, output_discount))
    parameters = (m, b, eligible, outside, alpha, t, s, coverage, output_discount)
    if min(parameters) < 0:
        raise ValueError('Local weights and potential arguments must be nonnegative')
    if not exact and not np.all(np.isfinite(np.array(parameters, dtype=JC_LD))):
        raise ValueError('Local weights and potential arguments must be finite')
    local = joint_mds_local(m, b, eligible, outside, alpha, coverage,
                            exact=exact, output_discount=output_discount)
    states = tuple((u, e) for u in range(5) for e in range(5 - u))
    index = {state: j for j, state in enumerate(states)}
    hit = np.zeros(5, dtype=dtype)
    empty_positive = np.zeros(5, dtype=dtype)
    for u in range(1, 5):
        hit[u] = comb(4, u) * sum(local[index[u, 0], index[v, 4 - v]]
                                  * t ** v * s ** (4 - v) for v in range(5))
    for e in range(5):
        empty_positive[e] = comb(4, e) * sum(
            local[index[0, e], index[v, 0]] * t ** v for v in range(1, 5))
    empty_zero = np.array([cast(comb(4, e)) for e in range(5)], dtype=dtype)
    empty = empty_zero + empty_positive
    groups = width // 4
    left = np.zeros((width + 1, groups + 1), dtype=dtype)
    right = np.zeros_like(left)
    hp = np.array([cast(1)], dtype=dtype)
    bp = np.array([cast(1)], dtype=dtype)
    b0p = np.array([cast(1)], dtype=dtype)
    dp = np.array([cast(0)], dtype=dtype)
    for k in range(groups + 1):
        if k:
            left[:len(hp), k] = hp * comb(groups, k)
        right[:len(bp), groups - k] = bp
        if k < groups:
            hp = np.convolve(hp, hit)
            dp = np.convolve(dp, empty) + np.convolve(b0p, empty_positive)
            b0p = np.convolve(b0p, empty_zero)
            bp = np.convolve(bp, empty)
    coefficients = left @ right.T
    coefficients[0, :] += dp
    result = np.zeros_like(coefficients)
    for u in range(width + 1):
        denominator = np.array([cast(comb(width, u) * comb(width - u, e))
                                for e in range(width - u + 1)], dtype=dtype)
        result[u, :width - u + 1] = coefficients[u, :width - u + 1] / denominator
    if not exact and (not np.all(np.isfinite(result)) or np.any(result < 0)):
        raise ArithmeticError('Invalid positive MDS action')
    if not exact and min(m, b, eligible, coverage, output_discount, t) > 0:
        # Whenever A/E can occupy disjoint quartets, an all-H route and
        # all-new-A choices give a strictly positive nonempty transition.
        # This also applies to the terminal evaluation s=0.  Do not let
        # physical underflow become a false -inf Chernoff upper bound.
        for u in range(width + 1):
            e = np.arange(width - u + 1, dtype=int)
            feasible = (u + 3) // 4 + (e + 3) // 4 <= groups
            if np.any(result[u, :width - u + 1][feasible] <= 0):
                raise ArithmeticError('Positive physical action underflow; narrow the parameter range')
    if not exact and t == 0 and min(s, alpha, eligible, coverage) > 0:
        # The exact u_new=0 slice vanishes when old u=0.  For positive old
        # availability, all-new-E on an all-H route is a positive witness.
        for u in range(1, width + 1):
            e = np.arange(width - u + 1, dtype=int)
            feasible = (u + 3) // 4 + (e + 3) // 4 <= groups
            if np.any(result[u, :width - u + 1][feasible] <= 0):
                raise ArithmeticError('Positive u_new=0 slice underflow')
    return result


def joint_chernoff_comparison(prepared, m, b, eligible, outside, alpha,
                              coverage=1, output_discount=1, exact=False):
    """Return B=Z+Bplus and positive last-step terminal comparison."""
    _joint_chernoff_dependencies()
    width, states = prepared['width'], prepared['states']
    keys, bins = prepared['keys'], prepared['bins']
    dimension = len(keys)
    old_u, old_e = states[:, 0], states[:, 1]
    zero_vector = prepared['zero_vector']
    if exact:
        if any(int(x) != x for pair in prepared['tilts'] for x in pair):
            raise ValueError('Exact checks require integer log2 Chernoff tilts')
        if int(prepared['logpotential']) != prepared['logpotential']:
            raise ValueError('Exact checks require an integer base log2 potential')
        if int(prepared['log_hole_potential']) != prepared['log_hole_potential']:
            raise ValueError('Exact checks require an integer base hole log2 potential')
        F = Fraction
        phi = [F(2) ** (int(prepared['logpotential']) * int(u)
                       + int(prepared['log_hole_potential']) * int(e)) if u else F(1)
               for u, e in states]
        actions = []
        tilt_values = []
        for kind, lt, ls in prepared['candidates']:
            tau, sigma = F(2) ** int(lt), F(2) ** int(ls)
            action = joint_chernoff_positive_action(width, m, b, eligible, outside, alpha,
                F(0) if kind == 'u_zero' else F(2) ** int(prepared['logpotential']) * tau,
                F(0) if kind == 'e_zero' else
                (sigma if kind == 'u_zero' else F(2) ** int(prepared['log_hole_potential']) * sigma),
                coverage, output_discount, exact=True)
            actions.append([action[u, e] for u, e in states])
            tilt_values.append((tau, sigma))
        positive = np.full((dimension, dimension), F(0), dtype=object)
        for target in range(1, dimension):
            target_states = states[bins == target]
            constants = [min(tau ** int(v) * sigma ** int(f) for v, f in target_states)
                         for tau, sigma in tilt_values]
            candidates = np.flatnonzero(prepared['candidate_allowed'][target])
            for source, (u, e) in enumerate(states):
                if (u == 0 and keys[target][1] != 0) or (4 * ((u + 3) // 4) > prepared['max_sum'][target]):
                    continue
                upper = min(actions[j][source] / constants[j] for j in candidates) / phi[source]
                positive[bins[source], target] = max(positive[bins[source], target], upper)
        terminal_action = joint_chernoff_positive_action(width, m, b, eligible, outside,
            alpha, F(1), F(0), coverage, output_discount, exact=True)
        terminal = np.full(dimension, F(0), dtype=object)
        for source, (u, e) in enumerate(states):
            terminal[bins[source]] = max(terminal[bins[source]], terminal_action[u, e] / phi[source])
        comparison = positive.copy()
        comparison[:, 0] = [F(int(value)) for value in zero_vector]
        return comparison, terminal, positive
    action_logs = np.empty((len(states), len(prepared['candidates'])), dtype=JC_LD)
    terminal_cached = None
    for j, (kind, lt, ls) in enumerate(prepared['candidates']):
        new_t = JC_LD(0) if kind == 'u_zero' else np.exp2(prepared['logpotential'] + JC_LD(lt))
        new_s = JC_LD(0) if kind == 'e_zero' else np.exp2(
            JC_LD(ls) + (JC_LD(0) if kind == 'u_zero' else prepared['log_hole_potential']))
        action = joint_chernoff_positive_action(width, m, b, eligible, outside, alpha,
            new_t, new_s, coverage, output_discount)
        if new_t == 1 and new_s == 0:
            terminal_cached = action
        values = action[old_u, old_e]
        action_logs[:, j] = JC_LD('-inf')
        mask = values > 0
        action_logs[mask, j] = np.log2(values[mask])
    positive_logs = np.full((dimension, dimension), JC_LD('-inf'), dtype=JC_LD)
    for target in range(1, dimension):
        candidates = prepared['candidate_allowed'][target]
        upper = np.min(action_logs[:, candidates] - prepared['minima'][target, candidates][None, :], axis=1) - prepared['logphi']
        impossible = 4 * ((old_u + 3) // 4) > prepared['max_sum'][target]
        if keys[target][1] != 0:
            impossible |= old_u == 0
        upper[impossible] = JC_LD('-inf')
        positive_logs[:, target] = np.maximum.reduceat(upper[prepared['order']], prepared['starts'])
    positive = np.exp2(positive_logs)
    if not np.all(np.isfinite(positive)) or np.any(np.isfinite(positive_logs) & (positive == 0)):
        raise ArithmeticError('Bin bound overflow/positive underflow')
    terminal_action = terminal_cached
    if terminal_action is None:
        terminal_action = joint_chernoff_positive_action(width, m, b, eligible, outside,
            alpha, JC_LD(1), JC_LD(0), coverage, output_discount)
    phi = np.exp2(prepared['logphi'])
    if not np.all(np.isfinite(phi)) or np.any(phi <= 0):
        raise ArithmeticError('Invalid fixed bin potential')
    terminal = np.maximum.reduceat((terminal_action[old_u, old_e] / phi)[prepared['order']], prepared['starts'])
    comparison = positive.copy()
    comparison[:, 0] = zero_vector
    return comparison, terminal, positive


def joint_chernoff_checks():
    _joint_chernoff_dependencies()
    F = Fraction
    weights = (F(3, 4), F(5, 4), F(1, 3), F(2, 9), F(2, 3))
    entries = inequalities = terminal_checks = intertwining = action_checks = two_quartet_checks = 0
    for output_discount, log_hole_potential in product((F(1), F(3, 4)), (0, -1)):
        prepared = joint_chernoff_prepare(4, 2, 2, -1, (-1, 0, 1), log_hole_potential)
        states = prepared['states']
        local = joint_mds_local(*weights, F(4, 5), exact=True, output_discount=output_discount)
        v = np.array([F(int(u == 0)) for u, e in states], dtype=object)
        zero = np.zeros_like(local)
        zero[:, 0] = v
        positive = local - zero
        assert np.all(positive >= 0)
        for t, s in ((F(1, 2), F(2)), (F(1), F(0)), (F(0), F(2))):
            actual = joint_chernoff_positive_action(4, *weights, t, s, F(4, 5), output_discount, exact=True)
            expected = positive @ np.array([t ** int(u) * s ** int(e) for u, e in states], dtype=object)
            for source, (u, e) in enumerate(states):
                assert actual[u, e] == expected[source]
                action_checks += 1
            # Independent two-quartet count convolution.  Split the old
            # labels between two fixed quartets hypergeometrically, then
            # combine the local positive actions and exact zero actions.
            # This tests the general coefficient normalization at width 8,
            # without reusing the 129-column outer-product construction.
            actual_eight = joint_chernoff_positive_action(8, *weights, t, s,
                F(4, 5), output_discount, exact=True)
            lookup = {tuple(state): j for j, state in enumerate(states)}
            multiplicity = [comb(4, int(u)) * comb(4 - int(u), int(e)) for u, e in states]
            for u in range(9):
                for e in range(9 - u):
                    denominator = comb(8, u) * comb(8 - u, e)
                    truth = F(0)
                    for first, (u1, e1) in enumerate(states):
                        second = lookup.get((u - u1, e - e1))
                        if second is None:
                            continue
                        probability = F(multiplicity[first] * multiplicity[second], denominator)
                        truth += probability * (expected[first] * expected[second]
                            + expected[first] * v[second] + v[first] * expected[second])
                    assert actual_eight[u, e] == truth
                    two_quartet_checks += 1
        exact, end, _ = joint_chernoff_comparison(prepared, *weights, F(4, 5), output_discount, exact=True)
        floating, float_end, _ = joint_chernoff_comparison(prepared, *weights, F(4, 5), output_discount)
        lift = np.zeros((len(states), len(prepared['keys'])), dtype=object)
        for source, (u, e) in enumerate(states):
            lift[source, prepared['bins'][source]] = F(2) ** (-int(u) + log_hole_potential * int(e)) if u else F(1)
        zbin = np.zeros_like(exact)
        zbin[:, 0] = [F(int(value)) for value in prepared['zero_vector']]
        assert np.array_equal(zero @ lift, lift @ zbin)
        assert np.array_equal(zbin @ zbin, zbin)
        intertwining += 2
        target = local @ lift
        upper = lift @ exact
        for i in range(target.shape[0]):
            for j in range(target.shape[1]):
                assert target[i, j] <= upper[i, j]
                inequalities += 1
        # Axes sharpen positive-u targets.  Zero-u targets must retain
        # their axis even in the baseline: the separable interior lift is
        # not their piecewise lift when the hole potential differs from 1.
        interior_only = dict(prepared)
        interior_only['candidate_allowed'] = prepared['candidate_allowed'].copy()
        for target, key in enumerate(prepared['keys']):
            if key[0] > 0:
                interior_only['candidate_allowed'][target, len(prepared['tilts']):] = False
        old_comparison, _, _ = joint_chernoff_comparison(interior_only, *weights,
            F(4, 5), output_discount, exact=True)
        assert np.all(exact <= old_comparison)
        terminal_truth = positive @ np.array([F(int(e == 0)) for u, e in states], dtype=object)
        assert all(a <= b for a, b in zip(terminal_truth, lift @ end))
        terminal_checks += len(states)
        for i in range(len(exact)):
            assert abs(float_end[i] - _jc_cast(end[i])) < JC_LD('2e-13') * max(1, abs(_jc_cast(end[i])))
            for j in range(len(exact)):
                value = _jc_cast(exact[i, j])
                assert abs(floating[i, j] - value) < JC_LD('2e-13') * max(1, abs(value))
                entries += 1
        # Independent exact finite-region comparison, including the free endpoint.
        exact_region = np.linalg.matrix_power(exact, 16) - zbin
        exact_end = np.linalg.matrix_power(exact, 15) @ (np.array(zbin[:, 0]) + end) - zbin[:, 0]
        region, scale, last, last_scale = joint_subset_region(floating, float_end, prepared['zero_vector'])
        for i in range(len(exact)):
            for j in range(len(exact)):
                value = _jc_cast(exact_region[i, j])
                assert abs(region[i, j] * np.exp2(scale) - value) < JC_LD('2e-11') * max(1, abs(value))
            value = _jc_cast(exact_end[i])
            assert abs(last[i] * np.exp2(last_scale) - value) < JC_LD('2e-11') * max(1, abs(value))
    print(dict(exact_chernoff_entries=entries, exact_physical_bin_inequalities=inequalities,
               positive_action_checks=action_checks, exact_two_quartet_checks=two_quartet_checks,
               terminal_checks=terminal_checks,
               zero_intertwining_checks=intertwining, full_code_certificate=False), flush=True)


def joint_chernoff_value(h=2048, theta=1/16, logz=-.4375, u_step=32, e_step=32,
                         logpotential=0, tilt_logs=(-1/16, 0, 1/16), tilts=3,
                         coverage=0, outside=0, outputgap=0, progress=False,
                         log_hole_potential=0):
    _joint_chernoff_dependencies()
    if not 0 < h < 8192 or not 0 < theta <= 1 or logz > 0 or tilts < 1:
        raise ValueError('Invalid field/row/tilt parameters')
    if coverage > 0 or outside > 0 or outputgap < 0:
        raise ValueError('Require coverage/outside<=0 and outputgap>=0')
    if not all(np.isfinite(value) for value in (theta, logz, coverage, outside, outputgap)):
        raise ValueError('Field weights and event markers must be finite')
    theta, logz = map(JC_LD, (theta, logz))
    q = JC_LD(2) ** 32
    x = np.exp2(logz / theta)
    c = (q + 1) * x
    lam = 1 + c
    d = max(min(JC_LD(1), lam / (q + 1)), (1 + q * x * x) / lam)
    m, b, alpha = c ** theta, lam ** theta, d ** theta
    outside_factor = np.exp2(JC_LD(outside) - 32 * theta)
    actual = JC_LD(h) / 8192
    field = 1 / outside_factor
    low = actual / (field - (field - 1) * actual)
    choose = JC_LD((lgamma(8193) - lgamma(h + 1) - lgamma(8193 - h)) / float(np.log(JC_LD(2))))
    prepared = joint_chernoff_prepare(512, u_step, e_step, logpotential, tilt_logs, log_hole_potential)
    logs = endlogs = None
    for tilt, fraction in enumerate(np.linspace(0, 1, tilts)):
        started = time.monotonic()
        odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
        p = odds / (1 + odds)
        logpi = choose + h * np.log2(p) + (8192 - h) * np.log2(1 - p)
        comparison, terminal, _ = joint_chernoff_comparison(prepared, m, b, p,
            (1 - p) * outside_factor, alpha, np.exp2(JC_LD(coverage)), np.exp2(-JC_LD(outputgap)))
        region, scale, last, last_scale = joint_subset_region(comparison, terminal, prepared['zero_vector'])
        candidate = np.full_like(region, JC_LD('-inf'))
        mask = region > 0
        candidate[mask] = np.log2(region[mask]) + scale - logpi
        candidate_end = np.full_like(last, JC_LD('-inf'))
        mask = last > 0
        candidate_end[mask] = np.log2(last[mask]) + last_scale - logpi
        logs = candidate if logs is None else np.minimum(logs, candidate)
        endlogs = candidate_end if endlogs is None else np.minimum(endlogs, candidate_end)
        if progress:
            print(dict(chernoff_quota_tilt=tilt + 1, dimension=len(comparison),
                       potential_actions=len(prepared['candidates']) + 1,
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    scale, end_scale = logs.max(), endlogs.max()
    region, last = np.exp2(logs - scale), np.exp2(endlogs - end_scale)
    if np.any(np.isfinite(logs) & (region == 0)) or np.any(np.isfinite(endlogs) & (last == 0)):
        raise ArithmeticError('Regional normalization discarded a positive term')
    tail = joint_subset_tail(region, scale, last, end_scale, prepared['zero_vector'])
    bound = (theta - 32 * theta - JC_W * logz + choose - 32 * theta * 128 * h + tail
             - (1009255 - 128 * h) * outside - 1736703 * coverage + JC_W * outputgap)
    return float(bound), dict(dimension=len(prepared['keys']), full_physical_states=len(prepared['states']),
        u_step=u_step, e_step=e_step, logpotential=logpotential, log_hole_potential=log_hole_potential,
        chernoff_tilts=len(prepared['tilts']),
        axis_actions=len(prepared['candidates']) - len(prepared['tilts']),
        quota_tilts=tilts, logcoverage=coverage, logoutside=outside, outputgap=outputgap,
        exact_MDS_filter=True, occupied_129_retained=True, certificate=False)


def joint_chernoff_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--logz', type=float, default=-.4375)
    parser.add_argument('--u-step', type=int, default=32)
    parser.add_argument('--e-step', type=int, default=32)
    parser.add_argument('--logpotential', type=float, default=0.)
    parser.add_argument('--log-hole-potential', type=float, default=0.)
    parser.add_argument('--chernoff-logs', type=float, nargs='+', default=(-1/16, 0., 1/16))
    parser.add_argument('--tilts', type=int, default=3)
    parser.add_argument('--coverage', type=float, default=0.)
    parser.add_argument('--outside', type=float, default=0.)
    parser.add_argument('--outputgap', type=float, default=0.)
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
            joint_chernoff_checks()
        else:
            started = time.monotonic()
            value, detail = joint_chernoff_value(args.h, args.theta, args.logz, args.u_step,
                args.e_step, args.logpotential, args.chernoff_logs, args.tilts,
                args.coverage, args.outside, args.outputgap, progress=True,
                log_hole_potential=args.log_hole_potential)
            print(dict(h=args.h, theta=args.theta, logz=args.logz, log2_bound=value,
                       seconds=round(time.monotonic() - started, 3), **detail), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    joint_chernoff_main()
