"""Disjoint-bin comparison for the adaptive Hall population operator (FLOAT).

The physical rows(a,w) have a divisible4 and0<=w<=a<=n. Rectangle bins
isolate the two corners(0,0),(n,n). Their positive lift
 H[s,I]=alpha0^a beta0^w 1[s in I]
is fixed for every destination Chernoff tilt and regional quota x.
For each destination J take the tilt minimum at EACH physical source,
remove exact intersection-infeasible destinations, then maximize over
the source bin. Exponential actions come from the frozen positive engine.

Compose16 updates before the quota bound x^-h/C(n*16,h) and entrywise
x minimum. The endpoint uses15 updates and q_x=(y+x)^w(1+x)^(n-w).
The whole bound is C*y^-threshold*phi0(n,n)*e_full Mbar^255 bbar.
There is no extra y^n factor. Bins and lift never change within a bound.

The helper reports destination-tilt overcount and source-bin variation
separately. Perron masses describe the comparison operator, not a
physical-history distribution. Arithmetic is positive longdouble with
structural-zero and range guards, but is not directed certification.
No jobs, files, or large allocations occur on import. Parent owns runs.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _ab_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_ab_thread_var] = '1'

import numpy as np

AB_LD = np.longdouble
AB_COMPARISON_WORKSPACE_LIMIT = 1 << 30


def _ab_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('adaptive_cover'):
        import quartet_adaptive_cover_checks as module
    else:
        import quartet_adaptive_population as module
    return getattr(module, name)


def _ab_normalize(array):
    if not np.all(np.isfinite(array)) or np.any(array < 0):
        raise ArithmeticError('Invalid nonnegative comparison array')
    positive = array[array > 0]
    if not len(positive):
        raise ArithmeticError('Zero comparison action')
    maximum = positive.max()
    result = array/maximum
    if np.any((array > 0) & (result == 0)):
        raise ArithmeticError('Positive comparison normalization underflow')
    return result, np.log2(maximum)


def _ab_from_logs(logs):
    finite = np.isfinite(logs)
    if np.any(np.isnan(logs)) or np.any(np.isposinf(logs)) or not np.any(finite):
        raise ArithmeticError('Invalid comparison logarithms')
    scale = logs[finite].max()
    result = np.zeros(logs.shape, dtype=AB_LD)
    result[finite] = np.exp2(logs[finite]-scale)
    if np.any(result[finite] <= 0) or not np.all(np.isfinite(result)):
        raise ArithmeticError('Positive comparison log alignment underflow/overflow')
    return result, scale


def _ab_logs(array, scale=0):
    result = np.full(array.shape, AB_LD('-inf'))
    positive = array > 0
    result[positive] = np.log2(array[positive])+scale
    return result


def _ab_product(left, left_scale, right, right_scale):
    # A positive operand pair must not disappear, even if another term masks it.
    if left[left > 0].min()*right[right > 0].min() <= 0:
        raise ArithmeticError('Positive matrix/vector product would underflow')
    support = (left > 0) @ (right > 0)
    result = left @ right
    if np.any(support & (result <= 0)) or np.any((~support) & (result != 0)):
        raise ArithmeticError('Comparison product changed structural support')
    result, extra = _ab_normalize(result)
    return result, left_scale+right_scale+extra


def _ab_power(matrix, scale, count):
    count = index(count)
    if count < 0:
        raise ValueError('Use a nonnegative power')
    result, result_scale = np.eye(matrix.shape[0], dtype=AB_LD), AB_LD(0)
    power, power_scale = matrix, scale
    while count:
        if count & 1:
            result, result_scale = _ab_product(result, result_scale, power, power_scale)
        count >>= 1
        if count:
            power, power_scale = _ab_product(power, power_scale, power, power_scale)
    return result, result_scale


def _ab_apply_power(matrix, scale, vector, vector_scale, count):
    count = index(count)
    if count < 0:
        raise ValueError('Use a nonnegative power')
    power, power_scale = matrix, scale
    while count:
        if count & 1:
            vector, vector_scale = _ab_product(power, power_scale, vector, vector_scale)
        count >>= 1
        if count:
            power, power_scale = _ab_product(power, power_scale, power, power_scale)
    return vector, vector_scale


def _ab_logsum_rows(logs):
    maximum = logs.max(axis=1)
    terms = np.exp2(logs-maximum[:, None])
    if np.any(np.isfinite(logs) & (terms <= 0)):
        raise ArithmeticError('Positive diagnostic sum underflow')
    return maximum+np.log2(terms.sum(axis=1, dtype=AB_LD))


def adaptive_bins_prepare(width=512, bin_width=64, alpha0=1, beta0=2,
                          tilt_offsets=(-1/2, 0, 1/2)):
    width, bin_width = index(width), index(bin_width)
    alpha0, beta0 = AB_LD(alpha0), AB_LD(beta0)
    if width < 4 or width > 512 or width % 4 or bin_width < 1:
        raise ValueError('Use quartet width4..512 and positive bin width')
    if min(alpha0, beta0) <= 0 or not np.isfinite(alpha0+beta0):
        raise ValueError('Use a finite positive common lift')
    offsets = tuple(AB_LD(value) for value in tilt_offsets)
    if not offsets or not all(np.isfinite(value) for value in offsets) or 0 not in offsets:
        raise ValueError('The bounded finite tilt grid must include zero')
    a = np.asarray([a for a in range(0, width+1, 4) for _ in range(a+1)], dtype=np.int16)
    w = np.asarray([w for a in range(0, width+1, 4) for w in range(a+1)], dtype=np.int16)
    labels, mapping, members, keys = [], {}, [], []
    for i, (aa, ww) in enumerate(zip(a, w)):
        key = ('empty',) if aa == ww == 0 else ('full',) if aa == ww == width else ('box', int(aa)//bin_width, int(ww)//bin_width)
        if key not in mapping:
            mapping[key] = len(members)
            members.append([])
            keys.append(key)
        label = mapping[key]
        labels.append(label)
        members[label].append(i)
    labels = np.asarray(labels, dtype=np.int32)
    members = tuple(np.asarray(row, dtype=np.int32) for row in members)
    # Reserve room for simultaneous row-bound/log/diagnostic temporaries and
    # scaled matrix products, plus the frozen action's full-square buffers.
    # This is an implementation workspace limit, not a mathematical bin limit.
    item_bytes = np.dtype(AB_LD).itemsize
    workspace = (item_bytes*(6*len(a)*len(members)+12*len(members)**2
                             +8*(width+1)**2+len(offsets)**2*len(members))
                 +(64 << 20))
    if workspace > AB_COMPARISON_WORKSPACE_LIMIT:
        raise ValueError('Estimated comparison workspace '+str(workspace)
                         +' bytes exceeds the1GiB workspace budget; use coarser bins')
    tilts = tuple(product(offsets, repeat=2))
    minima = np.empty((len(tilts), len(members)), dtype=AB_LD)
    old_a = np.arange(0, width+1, 4, dtype=np.int16)
    feasible = np.zeros((len(old_a), len(members)), dtype=bool)
    for j, rows in enumerate(members):
        target_a, target_w = a[rows], w[rows]
        feasible[:, j] = np.any((target_w[None, :] <= old_a[:, None])
                                & (old_a[:, None] <= width-target_a[None, :]+target_w[None, :]), axis=1)
        for t, (xi_a, xi_w) in enumerate(tilts):
            minima[t, j] = np.min(xi_a*target_a+xi_w*target_w)
    log_alpha, log_beta = np.log2(alpha0), np.log2(beta0)
    return dict(width=width, bin_width=bin_width, alpha0=alpha0, beta0=beta0,
                log2_alpha0=log_alpha, log2_beta0=log_beta, a=a, w=w,
                log2_lift=a*log_alpha+w*log_beta, bins=members, bin_ids=labels,
                bin_keys=tuple(keys), tilts=tilts, destination_minima=minima,
                feasibility_by_a=feasible, full_bin=int(labels[-1]),
                empty_bin=int(labels[0]), isolated_corners_only=True,
                estimated_workspace_bytes=workspace)


def adaptive_bins_evaluation(prepared, y=2, x=1, progress=False, retain_rows=False):
    started = time.monotonic()
    a, w, bins = prepared['a'], prepared['w'], prepared['bins']
    row_bounds = np.full((len(a), len(bins)), AB_LD('inf'))
    baseline = None
    for t, (xi_a, xi_w) in enumerate(prepared['tilts']):
        alpha = prepared['alpha0']*np.exp2(xi_a)
        beta = prepared['beta0']*np.exp2(xi_w)
        action = _ab_dependency('adaptive_population_action')(y, x, alpha, beta, prepared['width'])
        if not np.array_equal(a, action['a']) or not np.array_equal(w, action['w']):
            raise ValueError('Exponential action changed the physical row order')
        lifted = action['log2_action']-prepared['log2_lift']
        np.minimum(row_bounds, lifted[:, None]-prepared['destination_minima'][t], out=row_bounds)
        if xi_a == xi_w == 0:
            baseline = lifted.copy()
        if progress:
            print(dict(stage='adaptive_bin_action', x=float(x), tilt=(float(xi_a), float(xi_w)),
                       completed=t+1, total=len(prepared['tilts']),
                       action_seconds=round(action['seconds'], 3)), flush=True)
    row_bounds[~prepared['feasibility_by_a'][a//4]] = AB_LD('-inf')
    matrix_logs = np.asarray([np.max(row_bounds[rows], axis=0) for rows in bins], dtype=AB_LD)
    matrix, scale = _ab_from_logs(matrix_logs)
    physical_mass = _ab_logsum_rows(row_bounds)
    bin_mass = _ab_logsum_rows(matrix_logs)
    zero_tilt_bin_action = np.asarray([baseline[rows].max() for rows in bins], dtype=AB_LD)
    row_mass_inflation = bin_mass-zero_tilt_bin_action
    inflated_bin = int(np.argmax(row_mass_inflation))
    tilt_excess = physical_mass-baseline
    bin_excess = bin_mass[prepared['bin_ids']]-physical_mass
    endpoint_logs = (w*np.log2(AB_LD(y)+AB_LD(x))
                     +(prepared['width']-w)*np.log2(1+AB_LD(x))-prepared['log2_lift'])
    endpoint_logs = np.asarray([endpoint_logs[rows].max() for rows in bins], dtype=AB_LD)
    endpoint, endpoint_scale = _ab_from_logs(endpoint_logs)
    diagnostics = dict(x=float(x), bin_count=len(bins), physical_rows=len(a),
                       destination_tilt_cover_log2_excess_max=float(tilt_excess.max()),
                       destination_tilt_cover_log2_excess_mean=float(tilt_excess.mean()),
                       source_bin_max_log2_excess_max=float(bin_excess.max()),
                       source_bin_max_log2_excess_mean=float(bin_excess.mean()),
                       row_mass_inflation_log2_min=float(row_mass_inflation.min()),
                       row_mass_inflation_log2_max=float(row_mass_inflation.max()),
                       row_mass_inflation_max_bin=prepared['bin_keys'][inflated_bin],
                       exact_feasibility_zero_entries=int(np.count_nonzero(~np.isfinite(row_bounds))),
                       seconds=time.monotonic()-started,
                       losses_describe_comparison_not_construction=True)
    result = dict(matrix=matrix, scale=scale, endpoint=endpoint, endpoint_scale=endpoint_scale,
                  matrix_logs=matrix_logs, diagnostics=diagnostics)
    if retain_rows:
        result['physical_row_log_bounds'] = row_bounds
        result['physical_row_log_mass'] = baseline
    return result


def adaptive_bins_regional(prepared, y=2, x_values=(1, 2, 4), region_updates=16, quota=None, progress=False):
    region_updates = index(region_updates)
    slots = prepared['width']*region_updates
    quota = slots//2 if quota is None else index(quota)
    if region_updates < 2 or not 0 <= quota <= slots or not x_values:
        raise ValueError('Invalid regional length, quota, or evaluation family')
    log_choose = np.log2(AB_LD(comb(slots, quota)))
    ordinary_logs = terminal_logs = selectors = end_selectors = None
    summaries = []
    for index_x, x in enumerate(x_values):
        if not np.isfinite(x) or x <= 0:
            raise ValueError('Use positive finite quota evaluations')
        evaluation = adaptive_bins_evaluation(prepared, y, x, progress)
        matrix, scale = evaluation['matrix'], evaluation['scale']
        regional, regional_scale = _ab_power(matrix, scale, region_updates)
        endpoint, endpoint_scale = _ab_apply_power(matrix, scale, evaluation['endpoint'],
                                                  evaluation['endpoint_scale'], region_updates-1)
        correction = -quota*np.log2(AB_LD(x))-log_choose
        local_logs = _ab_logs(regional, regional_scale+correction)
        local_end = _ab_logs(endpoint, endpoint_scale+correction)
        if ordinary_logs is None:
            ordinary_logs, terminal_logs = local_logs, local_end
            selectors = np.full(local_logs.shape, index_x, dtype=np.int16)
            end_selectors = np.full(local_end.shape, index_x, dtype=np.int16)
        else:
            selectors[local_logs < ordinary_logs] = index_x
            end_selectors[local_end < terminal_logs] = index_x
            np.minimum(ordinary_logs, local_logs, out=ordinary_logs)
            np.minimum(terminal_logs, local_end, out=terminal_logs)
        summaries.append(evaluation['diagnostics'])
        if progress:
            print(dict(stage='adaptive_bin_region', x=float(x),
                       maximum_conditioned_entry_log2=float(local_logs.max()),
                       maximum_conditioned_endpoint_log2=float(local_end.max())), flush=True)
    ordinary, ordinary_scale = _ab_from_logs(ordinary_logs)
    endpoint, endpoint_scale = _ab_from_logs(terminal_logs)
    return dict(ordinary=ordinary, ordinary_scale=ordinary_scale, endpoint=endpoint,
                endpoint_scale=endpoint_scale, quota_log2_binomial=log_choose,
                x_values=tuple(x_values), summaries=summaries,
                quota_selector_counts=tuple(int(np.count_nonzero(selectors == i)) for i in range(len(x_values))),
                endpoint_selector_counts=tuple(int(np.count_nonzero(end_selectors == i)) for i in range(len(x_values))),
                common_disjoint_lift=True, compose_before_quota_minimum=True)


def _ab_perron_summary(prepared, regional):
    matrix = regional['ordinary']
    right = left = np.ones(len(prepared['bins']), dtype=AB_LD)
    for _ in range(64):
        right, _ = _ab_normalize(matrix @ right)
        left, _ = _ab_normalize(left @ matrix)
    mass = left*right
    if np.any(mass <= 0) or not np.all(np.isfinite(mass)):
        raise ArithmeticError('Positive Perron diagnostic product underflow/overflow')
    mass /= mass.sum(dtype=AB_LD)
    winner = int(np.argmax(mass))
    rows = prepared['bins'][winner]
    return dict(dominant_comparison_bin=prepared['bin_keys'][winner],
                dominant_comparison_bin_mass=float(mass[winner]),
                dominant_bin_a_range=(int(prepared['a'][rows].min()), int(prepared['a'][rows].max())),
                dominant_bin_w_range=(int(prepared['w'][rows].min()), int(prepared['w'][rows].max())),
                comparison_only_not_physical_history_distribution=True)


def adaptive_bins_bound(y=2, threshold=60518, alpha0=1, beta0=2, bin_width=64,
                        tilt_offsets=(-1/2, 0, 1/2), x_values=(1, 2, 4), progress=False):
    started = time.monotonic()
    prepared = adaptive_bins_prepare(512, bin_width, alpha0, beta0, tilt_offsets)
    regional = adaptive_bins_regional(prepared, y, x_values, 16, 4096, progress)
    endpoint, scale = _ab_apply_power(regional['ordinary'], regional['ordinary_scale'],
                                     regional['endpoint'], regional['endpoint_scale'], 255)
    full = prepared['full_bin']
    log_probability = (regional['quota_log2_binomial']-index(threshold)*np.log2(AB_LD(y))
                       +prepared['log2_lift'][-1]+scale+np.log2(endpoint[full]))
    return dict(y=float(y), threshold=index(threshold), alpha0=float(alpha0), beta0=float(beta0),
                bin_width=bin_width, bin_count=len(prepared['bins']), physical_rows=len(prepared['a']),
                estimated_workspace_bytes=prepared['estimated_workspace_bytes'],
                action_count=len(prepared['tilts'])*len(x_values),
                x_values=tuple(map(float, x_values)), tilt_offsets=tuple(map(float, tilt_offsets)),
                log2_probability_upper=float(log_probability),
                quota_selector_counts=regional['quota_selector_counts'],
                endpoint_selector_counts=regional['endpoint_selector_counts'],
                loss_diagnostics=regional['summaries'],
                comparison_diagnostics=_ab_perron_summary(prepared, regional),
                total_updates=4096, H_union_count=1, extra_y_to_width_factor=False,
                seconds=time.monotonic()-started, floating_diagnostic_only=True,
                full_code_certificate=False)


@lru_cache(maxsize=None)
def _ab_exact_kernel(width, y):
    """Independent labeled quartet/coin enumeration, averaged over physical partitions."""
    states = [(a, w) for a in range(0, width+1, 4) for w in range(a+1)]
    lookup = {state: i for i, state in enumerate(states)}
    partitions = _ab_dependency('_ap_partitions')(width)
    weights = _ab_dependency('adaptive_cover_weights')
    local_cache = {}
    kernel = [[[Fraction(0)]*(width+1) for _ in states] for _ in states]
    for source, (a, w) in enumerate(states):
        A, W = (1 << a)-1, (1 << w)-1
        for partition in partitions:
            running = {(0, 0, 0): Fraction(1)}
            for quartet in partition:
                old_A = sum(((A >> coordinate) & 1) << slot for slot, coordinate in enumerate(quartet))
                old_W = sum(((W >> coordinate) & 1) << slot for slot, coordinate in enumerate(quartet))
                key = (old_A, old_W)
                if key not in local_cache:
                    local = {}
                    u, v = old_A.bit_count(), old_W.bit_count()
                    for J in range(16):
                        j, r = J.bit_count(), (J & old_W).bit_count()
                        g, p, _, _, _ = weights(y, u, v, j, r)
                        for V, probability in ((0, p), (15, g)):
                            target = (V.bit_count(), (V & old_A).bit_count(), j)
                            local[target] = local.get(target, Fraction(0))+probability*y**((15 ^ J) & old_W).bit_count()
                    local_cache[key] = local
                after = {}
                for (aa, ww, jj), value in running.items():
                    for (da, dw, dj), local_value in local_cache[key].items():
                        target = (aa+da, ww+dw, jj+dj)
                        after[target] = after.get(target, Fraction(0))+value*local_value
                running = after
            for (aa, ww, j), value in running.items():
                kernel[source][lookup[(aa, ww)]][j] += value/len(partitions)
    return states, kernel


def _ab_mm(left, right):
    return [[sum(a*b for a, b in zip(row, column)) for column in zip(*right)] for row in left]


def _ab_mv(matrix, vector):
    return [sum(a*b for a, b in zip(row, vector)) for row in matrix]


def _ab_exact_bin_matrix(prepared, matrix, alpha0, beta0):
    states = list(zip(map(int, prepared['a']), map(int, prepared['w'])))
    lift = [alpha0**a*beta0**w for a, w in states]
    row_bounds = [[None]*len(prepared['bins']) for _ in states]
    for xi_a, xi_w in prepared['tilts']:
        fa, fw = Fraction(2)**int(xi_a), Fraction(2)**int(xi_w)
        potential = [value*fa**a*fw**w for value, (a, w) in zip(lift, states)]
        action = _ab_mv(matrix, potential)
        for j, rows in enumerate(prepared['bins']):
            denominator = min(fa**states[int(t)][0]*fw**states[int(t)][1] for t in rows)
            for s, (a, _) in enumerate(states):
                candidate = action[s]/(lift[s]*denominator) if prepared['feasibility_by_a'][a//4, j] else Fraction(0)
                old = row_bounds[s][j]
                row_bounds[s][j] = candidate if old is None else min(old, candidate)
    comparison = [[max(row_bounds[int(s)][j] for s in rows) for j in range(len(prepared['bins']))]
                  for rows in prepared['bins']]
    return comparison, row_bounds, lift


def adaptive_bins_checks():
    counters = dict(indicator=0, physical_bin=0, source_max=0, composition=0, quota=0, endpoint=0, full=0, floating=0)
    try:
        adaptive_bins_prepare(512, 1)
    except ValueError as error:
        assert 'workspace' in str(error)
    else:
        raise AssertionError('Unsupported fine-bin workspace was not rejected')
    for width in (4, 8):
        y, alpha0, beta0 = Fraction(2), Fraction(1), Fraction(2)
        prepared = adaptive_bins_prepare(width, 4, alpha0, beta0, (-1, 0, 1))
        assert prepared['bin_ids'].dtype == np.dtype(np.int32)
        states, kernel = _ab_exact_kernel(width, y)
        size, bins = len(states), prepared['bins']
        quota, C = width, Fraction(comb(2*width, width))
        B = [[sum(kernel[i][k][d]*kernel[k][j][quota-d]
                  for k in range(size) for d in range(width+1))/C for j in range(size)] for i in range(size)]
        terminal_polys = []
        for a, w in states:
            terminal_polys.append([sum(Fraction(comb(w, r)*comb(width-w, j-r))*y**(w-r)
                                       for r in range(max(0, j-(width-w)), min(w, j)+1)) for j in range(width+1)])
        terminal = [sum(kernel[i][k][d]*terminal_polys[k][quota-d]
                        for k in range(size) for d in range(width+1))/C for i in range(size)]
        region_min = endpoint_min = None
        for x in (Fraction(1), Fraction(2)):
            matrix = [[sum(coefficient*x**j for j, coefficient in enumerate(poly)) for poly in row] for row in kernel]
            comparison, row_bounds, lift = _ab_exact_bin_matrix(prepared, matrix, alpha0, beta0)
            H = [[lift[s] if int(prepared['bin_ids'][s]) == j else Fraction(0) for j in range(len(bins))] for s in range(size)]
            KH = _ab_mm(matrix, H)
            for s in range(size):
                for j in range(len(bins)):
                    assert KH[s][j] <= lift[s]*row_bounds[s][j]
                    assert row_bounds[s][j] <= comparison[int(prepared['bin_ids'][s])][j]
                    counters['physical_bin'] += 1
                    counters['source_max'] += 1
            for xi_a, xi_w in prepared['tilts']:
                fa, fw = Fraction(2)**int(xi_a), Fraction(2)**int(xi_w)
                for rows in bins:
                    denominator = min(fa**states[int(t)][0]*fw**states[int(t)][1] for t in rows)
                    for t in rows:
                        assert denominator <= fa**states[int(t)][0]*fw**states[int(t)][1]
                        counters['indicator'] += 1
            matrix2, comp2 = _ab_mm(matrix, matrix), _ab_mm(comparison, comparison)
            left, right = _ab_mm(matrix2, H), _ab_mm(H, comp2)
            assert all(a <= b for lr, rr in zip(left, right) for a, b in zip(lr, rr))
            counters['composition'] += size*len(bins)
            qx = [(y+x)**w*(1+x)**(width-w) for a, w in states]
            endpoint_bin = [max(qx[int(s)]/lift[int(s)] for s in rows) for rows in bins]
            local_region = [[value/(x**quota*C) for value in row] for row in comp2]
            local_end = [value/(x**quota*C) for value in _ab_mv(comparison, endpoint_bin)]
            if region_min is None:
                region_min, endpoint_min = local_region, local_end
            else:
                region_min = [[min(a, b) for a, b in zip(first, second)] for first, second in zip(region_min, local_region)]
                endpoint_min = [min(a, b) for a, b in zip(endpoint_min, local_end)]
            floating = adaptive_bins_evaluation(prepared, float(y), float(x))
            actual = floating['matrix']*np.exp2(floating['scale'])
            for i in range(len(bins)):
                for j in range(len(bins)):
                    exact = comparison[i][j]
                    if exact:
                        reference = AB_LD(exact.numerator)/exact.denominator
                        assert abs(actual[i, j]/reference-1) < AB_LD('1e-10')
                    else:
                        assert actual[i, j] == 0
                    counters['floating'] += 1
        BH, Hregion = _ab_mm(B, H), _ab_mm(H, region_min)
        assert all(a <= b for left, right in zip(BH, Hregion) for a, b in zip(left, right))
        counters['quota'] += size*len(bins)
        assert all(value <= lift[s]*endpoint_min[int(prepared['bin_ids'][s])] for s, value in enumerate(terminal))
        counters['endpoint'] += size
        exact_full = _ab_mv(B, terminal)[-1]
        upper_full = lift[-1]*_ab_mv(region_min, endpoint_min)[prepared['full_bin']]
        assert exact_full <= upper_full
        counters['full'] += 1
        floating_region = adaptive_bins_regional(prepared, float(y), (1, 2), 2, quota)
        out, out_scale = _ab_apply_power(floating_region['ordinary'], floating_region['ordinary_scale'],
                                        floating_region['endpoint'], floating_region['endpoint_scale'], 1)
        found = out[prepared['full_bin']]*np.exp2(out_scale)*AB_LD(lift[-1].numerator)/lift[-1].denominator
        reference = AB_LD(upper_full.numerator)/upper_full.denominator
        assert abs(found/reference-1) < AB_LD('1e-9')
        counters['floating'] += 1
    print(dict(**counters, workspace_guard_checks=1, exact_width8_partitions=35,
               all_operator_checks_tiny=True, full_code_certificate=False), flush=True)


def adaptive_bins_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--y', type=float, default=2.)
    parser.add_argument('--threshold', type=int, default=60518)
    parser.add_argument('--alpha0', type=float, default=1.)
    parser.add_argument('--beta0', type=float, default=2.)
    parser.add_argument('--bin-width', type=int, default=64)
    parser.add_argument('--tilt-step', type=float, default=1/2)
    parser.add_argument('--tilt-offsets', nargs='+', type=float)
    parser.add_argument('--x-values', nargs='+', type=float, default=(1., 2., 4.))
    args = parser.parse_args()
    if not np.isfinite(args.tilt_step) or args.tilt_step <= 0:
        raise ValueError('Use a finite positive tilt step')
    offsets = tuple(args.tilt_offsets) if args.tilt_offsets is not None else (-args.tilt_step, 0, args.tilt_step)
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            adaptive_bins_checks()
        else:
            print(adaptive_bins_bound(args.y, args.threshold, args.alpha0, args.beta0,
                                      args.bin_width, offsets, args.x_values, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_bins_main()
