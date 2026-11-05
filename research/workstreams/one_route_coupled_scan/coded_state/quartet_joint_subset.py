"""Positive global joint availability / covered-hole-subset transfer.

The state is (u,e): all u available positions and a uniformly retained
e-subset of the unavailable covered positions, with e in fixed levels J.
Every old retained hole must lie in a nonfull new covered quartet.  New
actual holes are Gamma(Aprev) minus Anew; retain floor_J(their number)
positions uniformly.  The auxiliary sampling law has total weight one.

The terminal step requires the ACTUAL new hole set to be empty.  Merely
requiring the retained subset to be empty would relax the free endpoint.
The zero operator is Z=v*e00^T, v_(u,e)=1[u=0], including e>0 states.

All coefficient arithmetic is positive.  A Horner contraction retains
exact new-hole degrees below max(J), plus one saturated tail.  Its cost
is O(n^3 max(J) |J|); the default J=(0,1) has cubic assembly cost.
States with ceil(u/4)+ceil(e/4)>n/4 are pruned: they have no separating
partition, hence identically zero outgoing kernel and terminal weight.

Fix candidate mask histories before averaging partitions.  Apply the
canonical field/count inequality before dropping unretained-hole
constraints.  All quota tilts use this same physical count basis.
No local helper import is required.  --probe is FLOAT, not a certificate.
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


JS_LD = np.longdouble
JS_WIDTH = 512
JS_W = 943718


def joint_subset_levels(width, levels):
    levels = tuple(int(value) for value in levels)
    if width < 4 or width % 4 or not levels or levels[0] != 0:
        raise ValueError('Require a positive quartet width and levels starting at zero')
    if levels != tuple(sorted(set(levels))) or levels[-1] > width:
        raise ValueError('Levels must be distinct, increasing, and at most the width')
    return levels


def joint_subset_states(width, levels):
    levels = joint_subset_levels(width, levels)
    return tuple((u, e) for u in range(width + 1) for e in levels
                 if (u + 3) // 4 + (e + 3) // 4 <= width // 4)


def _joint_integer_convolve(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for j, value in enumerate(right):
        if value:
            for i, coefficient in enumerate(left):
                result[i + j] += coefficient * value
    return result


def joint_subset_prepare(width, levels=(0, 1)):
    """Combinatorial data independent of every field and quota parameter."""
    levels = joint_subset_levels(width, levels)
    groups = width // 4
    states = joint_subset_states(width, levels)
    choose = tuple(np.array([JS_LD(comb(n, j)) for j in range(n + 1)])
                   for n in range(width + 1))
    occupancy = [[1]]
    for _ in range(groups):
        occupancy.append(_joint_integer_convolve(occupancy[-1], [0, 4, 6, 4, 1]))
    probabilities = np.zeros((width + 1, groups + 1), dtype=JS_LD)
    for k in range(groups + 1):
        for u, coefficient in enumerate(occupancy[k]):
            if coefficient:
                probabilities[u, k] = JS_LD(comb(groups, k) * coefficient) / JS_LD(comb(width, u))
    partial = [np.ones(1, dtype=JS_LD)]
    for _ in range(min(levels[-1], groups)):
        partial.append(np.convolve(partial[-1], np.array([1, 4, 6, 4], dtype=JS_LD)))
    # G_(s,e)(t)=sum_j C(s,j)c_(j,e) R(t)^j(1+t)^(4(s-j)).
    geometry = {}
    for e in levels:
        rows = []
        for s in range(groups + 1):
            polynomial = np.zeros(4 * s + 1, dtype=JS_LD)
            for j in range(min(e, s) + 1):
                if e >= len(occupancy[j]) or not occupancy[j][e]:
                    continue
                term = np.convolve(partial[j], choose[4 * (s - j)])
                polynomial[:len(term)] += JS_LD(comb(s, j) * occupancy[j][e]) * term
            rows.append(polynomial)
        geometry[e] = rows
    by_level = {}
    for e in levels:
        selected = np.array([i for i, (_, tag) in enumerate(states) if tag == e], dtype=int)
        sizes = np.array([states[i][0] for i in selected], dtype=int)
        by_level[e] = selected, sizes
    zero_vector = np.array([JS_LD(u == 0) for u, _ in states])
    return dict(width=width, levels=levels, states=states, choose=choose,
                probabilities=probabilities, geometry=geometry, by_level=by_level,
                zero_vector=zero_vector)


def _joint_horner_background(width, u, e, maximum, coefficients, scaled_geometry, b):
    """Evaluate sum_H c_H G_H(y)(w+b*y)^H with a saturated w tail."""
    groups = width // 4
    last_k = min(u, groups - (e + 3) // 4)
    first_k = (u + 3) // 4
    if last_k < first_k:
        return np.zeros((maximum + 1, width - u + 1), dtype=JS_LD)
    highest = 4 * last_k - u
    degree = width - 4 * last_k
    current = np.zeros((maximum + 1, width - u + 1), dtype=JS_LD)
    scratch = np.zeros_like(current)
    current[0, :degree + 1] = coefficients[last_k] * scaled_geometry[groups - last_k]
    for power in range(highest - 1, -1, -1):
        # Multiply by w+b*y.  The last row means w-degree at least maximum.
        scratch[:, 0] = 0
        np.multiply(current[:, :degree + 1], b, out=scratch[:, 1:degree + 2])
        if maximum:
            scratch[1:maximum, :degree + 1] += current[:maximum - 1, :degree + 1]
            scratch[maximum, :degree + 1] += current[maximum - 1, :degree + 1]
        scratch[maximum, :degree + 1] += current[maximum, :degree + 1]
        degree += 1
        current, scratch = scratch, current
        if (power + u) % 4 == 0:
            k = (power + u) // 4
            if first_k <= k <= last_k and coefficients[k]:
                polynomial = scaled_geometry[groups - k]
                current[0, :len(polynomial)] += coefficients[k] * polynomial
    return current


def joint_subset_matrix(width, m, b, p, r, alpha, levels=(0, 1), coverage=1.,
                        prepared=None, progress=False):
    """Return K, terminal-positive t-v, and the common zero vector v."""
    if prepared is None:
        prepared = joint_subset_prepare(width, levels)
    elif prepared['width'] != width or prepared['levels'] != tuple(levels):
        raise ValueError('Prepared geometry does not match width/levels')
    m, b, p, r, alpha, coverage = map(JS_LD, (m, b, p, r, alpha, coverage))
    if min(m, b, p, r, alpha, coverage) < 0:
        raise ValueError('All weights must be nonnegative; alpha may exceed one')
    a = m * p * coverage
    groups = width // 4
    maximum = levels[-1]
    states = prepared['states']
    indices = np.arange(width + 1, dtype=int)
    apowers, bpowers, alphapowers = a ** indices, b ** indices, alpha ** indices
    route_powers = (r * coverage) ** (4 * np.arange(groups + 1, dtype=int))
    scaled = {e: [poly * apowers[:len(poly)] for poly in prepared['geometry'][e]] for e in levels}
    scaled_sums = {e: np.array([poly.sum() for poly in rows]) for e, rows in scaled.items()}
    old_polynomials = [prepared['choose'][u] * bpowers[:u + 1] * alphapowers[u::-1]
                       for u in range(width + 1)]
    result = np.zeros((len(states), len(states)), dtype=JS_LD)
    terminal_positive = np.zeros(len(states), dtype=JS_LD)
    started = time.monotonic()
    for row, (u, e) in enumerate(states):
        coefficients = prepared['probabilities'][u] * route_powers / JS_LD(comb(width - u, e))
        background = _joint_horner_background(width, u, e, maximum, coefficients, scaled[e], b)
        if maximum == 0:
            outgoing = np.convolve(background[0], old_polynomials[u])[None, :]
        else:
            outgoing = np.zeros((maximum + 1, width + 1), dtype=JS_LD)
            for new_holes in range(maximum):
                for old_holes in range(min(u, new_holes) + 1):
                    shift = u - old_holes
                    factor = JS_LD(comb(u, old_holes)) * alphapowers[old_holes] * bpowers[shift]
                    outgoing[new_holes, shift:shift + width - u + 1] += factor * background[new_holes - old_holes]
            outgoing[maximum] = np.convolve(background[maximum], old_polynomials[u])
            for background_holes in range(maximum):
                needed = maximum - background_holes
                if needed <= u:
                    # Old polynomial terms with old-hole count >=needed are a prefix in y-degree.
                    addition = np.convolve(background[background_holes], old_polynomials[u][:u - needed + 1])
                    outgoing[maximum, :len(addition)] += addition
        for j, kept in enumerate(levels):
            columns, sizes = prepared['by_level'][kept]
            if kept == maximum:
                polynomial = outgoing[maximum]
            else:
                polynomial = outgoing[kept:levels[j + 1]].sum(axis=0)
            result[row, columns] = polynomial[sizes]
        if u == 0:
            result[row, 0] = 1  # Exact empty-coverage transition, including old retained holes.
            # Only k=0 is possible.  Remove the empty endpoint positively.
            terminal_positive[row] = scaled[e][groups][1:].sum() / JS_LD(comb(width, e))
        else:
            terminal_positive[row] = sum(coefficients[k] * bpowers[4 * k] * scaled_sums[e][groups - k]
                                         for k in range((u + 3) // 4, min(u, groups - (e + 3) // 4) + 1))
        if progress and (row + 1) % 128 == 0:
            print(dict(joint_subset_stage='matrix_rows', completed=row + 1, dimension=len(states),
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    if not np.all(np.isfinite(result)) or not np.all(np.isfinite(terminal_positive)):
        raise ArithmeticError('Nonfinite joint-subset transfer; diagnostic invalid')
    return result, terminal_positive, prepared['zero_vector'].copy()


def joint_subset_region(C, terminal_positive, zero_vector):
    """K^16-Z and K^15*t-v, using Z=v*e00^T and positive arithmetic."""
    P = C.copy()
    P[:, 0] -= zero_vector
    if np.any(P < 0):
        raise ArithmeticError('Invalid zero-coverage operator')
    Q = P.copy()
    scale = JS_LD(0)
    for _ in range(4):
        norm = Q.sum(axis=1).max()
        Q /= norm
        scale += np.log2(norm)
        common = max(scale, 2 * scale)
        new = (Q @ Q) * np.exp2(2 * scale - common)
        linear = np.exp2(scale - common)
        new += linear * zero_vector[:, None] * Q[0, :][None, :]
        new[:, 0] += linear * (Q @ zero_vector)
        Q, scale = new, common
    norm = Q.sum(axis=1).max()
    Q /= norm
    scale += np.log2(norm)
    end = terminal_positive.copy()
    norm = end.max()
    end /= norm
    end_scale = np.log2(norm)
    injection = P @ zero_vector
    for _ in range(15):
        common = max(JS_LD(0), end_scale)
        end = (C @ end) * np.exp2(end_scale - common) + injection * np.exp2(-common)
        norm = end.max()
        end /= norm
        end_scale = common + np.log2(norm)
    return Q, scale, end, end_scale


def joint_subset_tail(Q, scale, end, end_scale, zero_vector, regions=256, minimum=129):
    prefix = regions - 1
    row = np.zeros(len(Q), dtype=JS_LD)
    row[0] = 1
    row_scale = JS_LD(0)
    aa = np.empty(prefix + 1, dtype=JS_LD)
    bb = np.empty(prefix + 1, dtype=JS_LD)
    ln2 = np.log(JS_LD(2))
    for k in range(prefix + 1):
        left, right = row @ zero_vector, row @ end
        aa[k] = np.log2(left) + row_scale if left else -np.inf
        bb[k] = np.log2(right) + row_scale if right else -np.inf
        if k < prefix:
            row = row @ Q
            norm = row.sum()
            row /= norm
            row_scale += np.log2(norm)
    logs = []
    for last_positive in (False, True):
        poly = bb.copy() if last_positive else aa.copy()
        lower = minimum - int(last_positive)
        for zero in range(prefix - lower + 1):
            occupied = prefix - zero
            logs.append(poly[occupied] + scale * occupied + (end_scale if last_positive else 0))
            after = np.full(prefix + 1, JS_LD('-inf'))
            for j in range(prefix + 1):
                after[j:] = np.logaddexp(after[j:] * ln2,
                                         (poly[:prefix + 1 - j] + aa[j]) * ln2) / ln2
            poly = after
    largest = max(logs)
    return largest + np.log2(np.exp2(np.array(logs) - largest).sum())


def _joint_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        return [(first, 255 ^ first) for other in combinations(range(1, 8), 3)
                for first in [1 + sum(1 << j for j in other)]]
    raise ValueError('Independent mask checks use widths four and eight')


def _joint_exact_row(width, old, retained, levels, m, b, p, r, alpha, coverage, expand_tags=False):
    """Independent physical partition/mask enumeration, with no polynomial formula."""
    parts = _joint_partitions(width)
    census = {}
    terminal_census = {}
    for partition in parts:
        hit = 0
        for quad in partition:
            if old & quad:
                hit |= quad
        if retained & hit:
            continue
        for new in range(1 << width):
            covered = hit | new
            if any(retained & quad and covered & quad == quad for quad in partition):
                continue
            actual_holes = hit & ~new
            f = actual_holes.bit_count()
            g = max(level for level in levels if level <= f)
            born = (new & ~hit).bit_count()
            closed = (old & ~new).bit_count()
            key = (new.bit_count(), g, hit.bit_count(), born, closed)
            if expand_tags:
                positions = [j for j in range(width) if actual_holes >> j & 1]
                for chosen in combinations(positions, g):
                    mask = sum(1 << j for j in chosen)
                    full_key = (new, mask, *key[2:])
                    census[full_key] = census.get(full_key, Fraction(0)) + Fraction(1, comb(f, g))
            else:
                census[key] = census.get(key, 0) + 1
            if f == 0:
                terminal_census[key] = terminal_census.get(key, 0) + 1
    values = {}
    for key, count in census.items():
        new_count = key[0].bit_count() if expand_tags else key[0]
        hit_count, born, closed = key[2:]
        weight = ((r * coverage) ** hit_count * (m * p * coverage) ** born
                  * b ** (new_count - born) * alpha ** closed / len(parts))
        destination = key[:2]
        values[destination] = values.get(destination, Fraction(0)) + count * weight
    terminal = Fraction(0)
    for key, count in terminal_census.items():
        new_count, _, hit_count, born, closed = key
        terminal += count * (r * coverage) ** hit_count * (m * p * coverage) ** born * b ** (new_count - born) * alpha ** closed / len(parts)
    return values, terminal


def joint_subset_checks():
    F = Fraction
    m, b, p, r, alpha, coverage = F(3, 4), F(5, 4), F(1, 3), F(2, 5), F(2, 3), F(4, 5)
    entries = terminals = omitted = paths = zero_checks = occupied_checks = subset_checks = 0
    for width in (4, 8):
        for levels in ((0,), (0, 1), (0, 2, 4)):
            prepared = joint_subset_prepare(width, levels)
            states = prepared['states']
            index = {state: i for i, state in enumerate(states)}
            exact = [[F(0)] * len(states) for _ in states]
            exact_terminal = [F(0)] * len(states)
            for u in range(width + 1):
                for e in levels:
                    if u + e > width:
                        continue
                    old = (1 << u) - 1
                    retained = ((1 << e) - 1) << u
                    row, terminal = _joint_exact_row(width, old, retained, levels, m, b, p, r, alpha, coverage)
                    if (u, e) not in index:
                        assert not row and terminal == 0
                        omitted += 1
                        continue
                    source = index[u, e]
                    exact_terminal[source] = terminal
                    for destination, weight in row.items():
                        if destination in index:
                            exact[source][index[destination]] += weight
            C, terminal_positive, zero_vector = joint_subset_matrix(width, float(m), float(b), float(p), float(r),
                                                                    float(alpha), levels, float(coverage), prepared)
            for i in range(len(states)):
                truth = JS_LD(exact_terminal[i].numerator) / JS_LD(exact_terminal[i].denominator)
                assert abs(terminal_positive[i] + zero_vector[i] - truth) <= JS_LD('3e-13') * max(1, abs(truth))
                terminals += 1
                for j in range(len(states)):
                    truth = JS_LD(exact[i][j].numerator) / JS_LD(exact[i][j].denominator)
                    assert abs(C[i, j] - truth) <= JS_LD('3e-13') * max(1, abs(truth))
                    entries += 1
            zero, endpoint_zero, zvec = joint_subset_matrix(width, float(m), float(b), float(p), float(r),
                                                           float(alpha), levels, 0., prepared)
            expected_zero = np.zeros_like(zero)
            expected_zero[:, 0] = zvec
            assert np.array_equal(zero, expected_zero) and np.all(endpoint_zero == 0)
            if len(levels) > 1:
                assert any(zvec[i] == 1 and states[i][1] > 0 for i in range(len(states)))
            zero_checks += 1
            Q, scale, end, end_scale = joint_subset_region(C, terminal_positive, zero_vector)
            direct = np.linalg.matrix_power(C, 16) - expected_zero
            assert np.max(abs(Q * np.exp2(scale) - direct)) <= JS_LD('3e-12') * max(1, direct.max())
            full_terminal = terminal_positive + zero_vector
            direct_end = np.linalg.matrix_power(C, 15) @ full_terminal - zero_vector
            assert np.max(abs(end * np.exp2(end_scale) - direct_end)) <= JS_LD('3e-12') * max(1, direct_end.max())
            positive = Q * np.exp2(scale)
            final = end * np.exp2(end_scale)
            enumerated = {(0, 0): JS_LD(1)}
            for _ in range(2):
                after = {}
                for (occupied, i), weight in enumerated.items():
                    if zero_vector[i]:
                        after[occupied, 0] = after.get((occupied, 0), JS_LD(0)) + weight
                    for j in range(len(states)):
                        key = (occupied + 1, j)
                        after[key] = after.get(key, JS_LD(0)) + weight * positive[i, j]
                enumerated = after
            for minimum in (1, 2, 3):
                truth = sum(weight * ((zero_vector[i] if occupied >= minimum else 0)
                                      + (final[i] if occupied + 1 >= minimum else 0))
                            for (occupied, i), weight in enumerated.items())
                value = joint_subset_tail(Q, scale, end, end_scale, zero_vector, regions=3, minimum=minimum)
                assert abs(value - np.log2(truth)) <= JS_LD('3e-12')
                occupied_checks += 1
            if width == 4:
                physical = [(a_mask, e_mask) for a_mask in range(1 << width)
                            for e_mask in range(1 << width) if not (a_mask & e_mask)
                            and (a_mask.bit_count(), e_mask.bit_count()) in index]
                physical_index = {state: i for i, state in enumerate(physical)}
                matrix = [[F(0)] * len(physical) for _ in physical]
                endpoint = [F(0)] * len(physical)
                for i, (old, retained) in enumerate(physical):
                    row, endpoint[i] = _joint_exact_row(width, old, retained, levels, m, b, p, r,
                                                       alpha, coverage, expand_tags=True)
                    collapsed = [F(0)] * len(states)
                    for destination, weight in row.items():
                        if destination in physical_index:
                            matrix[i][physical_index[destination]] += weight
                            collapsed[index[destination[0].bit_count(), destination[1].bit_count()]] += weight
                    assert collapsed == exact[index[old.bit_count(), retained.bit_count()]]
                    subset_checks += 1
                row = [F(int(state == (0, 0))) for state in physical]
                reduced = [F(int(i == 0)) for i in range(len(states))]
                for _ in range(3):
                    assert sum(row[i] * endpoint[i] for i in range(len(physical))) == sum(reduced[i] * exact_terminal[i] for i in range(len(states)))
                    row = [sum(row[i] * matrix[i][j] for i in range(len(physical))) for j in range(len(physical))]
                    reduced = [sum(reduced[i] * exact[i][j] for i in range(len(states))) for j in range(len(states))]
                    paths += 1
    print(dict(joint_subset_matrix_entries=entries, exact_terminal_checks=terminals,
               dead_state_checks=omitted, physical_uniform_subset_rows=subset_checks,
               physical_composed_free_terminal_checks=paths, general_zero_operator_checks=zero_checks,
               occupied_region_checks=occupied_checks, full_code_certificate=False), flush=True)


def joint_subset_field_weights(theta, logz):
    theta, logz = JS_LD(theta), JS_LD(logz)
    if not 0 < theta <= 1 or logz > 0:
        raise ValueError('Require 0<theta<=1 and log2(z)<=0')
    q = JS_LD(2) ** 32
    x = np.exp2(logz / theta)
    c = (q + 1) * x
    lam = 1 + c
    d = max(min(JS_LD(1), lam / (q + 1)), (1 + q * x * x) / lam)
    return c ** theta, lam ** theta, d ** theta, d


def joint_subset_value(h=2048, theta=1/16, logz=-.4375, tilts=3, levels=(0, 1), progress=False, profile=False):
    if tilts < 1 or not 0 < h < 8192:
        raise ValueError('Require positive tilt count and 0<h<8192')
    m, b, alpha, d = joint_subset_field_weights(theta, logz)
    field = np.exp2(JS_LD(32 * theta))
    prepared = joint_subset_prepare(JS_WIDTH, levels)
    rows = 8192
    actual = JS_LD(h) / rows
    low = actual / (field - (field - 1) * actual)
    choose = (lgamma(rows + 1) - lgamma(h + 1) - lgamma(rows - h + 1)) / float(np.log(JS_LD(2)))
    logs = endlogs = None
    for tilt_index, fraction in enumerate(np.linspace(0., 1., tilts)):
        started = time.monotonic()
        odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
        p = odds / (1 + odds)
        logpi = choose + h * np.log2(p) + (rows - h) * np.log2(1 - p)
        r = p + (1 - p) / field
        C, terminal_positive, zero_vector = joint_subset_matrix(JS_WIDTH, m, b, p, r, alpha,
                                                                levels, prepared=prepared, progress=progress)
        if progress:
            print(dict(joint_subset_stage='matrix', tilt=tilt_index + 1, dimension=len(C),
                       seconds=round(time.monotonic() - started, 3)), flush=True)
        Q, scale, end, end_scale = joint_subset_region(C, terminal_positive, zero_vector)
        if not np.all(Q > 0) or not np.all(end > 0):
            raise ArithmeticError('Unexpected regional zero/underflow on pruned physical basis')
        candidate = np.log2(Q) + scale - logpi
        candidate_end = np.log2(end) + end_scale - logpi
        logs = candidate if logs is None else np.minimum(logs, candidate)
        endlogs = candidate_end if endlogs is None else np.minimum(endlogs, candidate_end)
        if progress:
            print(dict(joint_subset_stage='region', tilt=tilt_index + 1,
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    scale = logs.max()
    Q = np.exp2(logs - scale)
    norm = Q.sum(axis=1).max()
    Q /= norm
    scale += np.log2(norm)
    end_scale = endlogs.max()
    end = np.exp2(endlogs - end_scale)
    value = joint_subset_tail(Q, scale, end, end_scale, prepared['zero_vector'])
    bound = float(theta - 32 * theta - JS_W * logz + choose - 32 * theta * 128 * h + value)
    detail=dict(dimension=len(prepared['states']), levels=tuple(levels), tilts=tilts,
                alpha=float(alpha), closure_d=float(d))
    if profile:
        # All-positive-region term only; not statistics of actual code histories.
        forward=np.zeros(len(Q),dtype=JS_LD);forward[0]=1
        backward=end.copy();backward/=backward.max()
        for _ in range(128):forward=forward@Q;forward/=forward.sum()
        for _ in range(127):backward=Q@backward;backward/=backward.max()
        marginal=forward*backward;marginal/=marginal.sum()
        top=np.argsort(marginal)[-10:][::-1]
        detail['relaxed_midpoint_top_u_e']=[(prepared['states'][i],float(marginal[i])) for i in top]
        detail['relaxed_mean_available']=float(sum(u*marginal[i] for i,(u,e) in enumerate(prepared['states'])))
        detail['retained_hole_mass']=[(e,float(sum(marginal[i] for i,(_,ee) in enumerate(prepared['states']) if ee==e))) for e in levels]
    return bound,detail


def joint_subset_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--tilts', type=int, default=3)
    parser.add_argument('--levels', default='0,1')
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--logz', type=float, default=-.4375)
    parser.add_argument('--profile', action='store_true')
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
            joint_subset_checks()
            return
        started = time.monotonic()
        levels = tuple(int(value) for value in args.levels.split(','))
        bound, detail = joint_subset_value(args.h, args.theta, args.logz, args.tilts, levels, progress=True,profile=args.profile)
        print(dict(joint_subset_rows=args.h, theta=args.theta, logz=args.logz,
                   log2_diagnostic=bound, seconds=round(time.monotonic() - started, 3),
                   certificate=False, **detail), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    joint_subset_main()
