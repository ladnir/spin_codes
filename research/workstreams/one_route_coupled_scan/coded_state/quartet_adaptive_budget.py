"""Nonlinear adaptive comparison preserving each physical row's total mass.

Use the frozen common lift and physical destination bounds m_sJ. The
zero-tilt action gives t_s=K(phi0)(s)/phi0(s). Set c_sJ=min(1,m_sJ/t_s),
preserving structural zeros and never renormalizing the caps. For v>=0,
sort its bin values once, then fill at most one unit of mass in that
order. This computes F_s(v)=max{p.v:0<=p<=c_s,sum p<=1}; exact caps
from the true mass have sum>=1, so this also attains the equality budget.
N_x(v)_I=max_(s in I)t_s F_s(v) is monotone and homogeneous.

Compose16 N_x maps before the regional quota correction and x minimum.
Intersect this vector bound with old Mbar v. Likewise intersect the
15-step nonlinear endpoint with the old endpoint. A positive witness
v with T(v)<=rho*v and endpoint<=mu*v gives the whole bound
 C*y^-threshold*phi0(full)*v_full*rho^255*mu.
The old linear full-power bound is retained as a control. The default
evaluates two nonlinear witnesses, beginning at64 old-matrix iterations.

All production arithmetic is FLOAT. Positive-range guards do not make
it a directed certificate. The explicit max(0,1-prefix) in the budget
algorithm is the mathematical remaining budget, not negative clipping
of a probability or a cap renormalization. Tiny numerical cap deficits
are reported; substantial deficits are rejected, never repaired.
No jobs or large allocations occur on import. Parent owns execution.
"""

import argparse
from fractions import Fraction
from itertools import product
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _bu_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_bu_thread_var] = '1'

import numpy as np

BU_LD = np.longdouble
BU_BATCH = 2048
BU_WORKSPACE_LIMIT = 3*(1 << 29)


def _bu_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_adaptive_bins as module
    return getattr(module, name)


def _bu_workspace(prepared, evaluations):
    rows, bins = len(prepared['a']), len(prepared['bins'])
    estimate = (prepared['estimated_workspace_bytes']
                +np.dtype(BU_LD).itemsize*(evaluations*(rows*bins+2*rows)
                                          +3*min(BU_BATCH, rows)*bins+8*rows))
    if estimate > BU_WORKSPACE_LIMIT:
        raise ValueError('Estimated nonlinear cache workspace exceeds the1.5GiB budget')
    return estimate


def _bu_map_data(prepared, evaluation, x):
    """Consume retained row logs as caps, keeping one common row-mass scale."""
    caps = evaluation.pop('physical_row_log_bounds')
    log_mass = evaluation.pop('physical_row_log_mass')
    supported = np.isfinite(caps)
    np.subtract(caps, log_mass[:, None], out=caps)
    np.minimum(caps, BU_LD(0), out=caps)  # Exactly min(1,m/t), in logarithms.
    np.exp2(caps, out=caps)
    if (np.any(caps[supported] <= 0) or np.any(caps[~supported] != 0)
            or np.any(caps > 1) or not np.all(np.isfinite(caps))):
        raise ArithmeticError('Physical budget caps changed support or range')
    cap_sums = caps.sum(axis=1, dtype=BU_LD)
    if cap_sums.min() < 1-BU_LD('1e-10'):
        raise ArithmeticError('Computed caps do not contain a unit-mass row')
    mass_scale = log_mass.max()
    mass = np.exp2(log_mass-mass_scale)
    if np.any(mass <= 0) or not np.all(np.isfinite(mass)):
        raise ArithmeticError('Physical row mass normalization underflow/overflow')
    return dict(x=BU_LD(x), caps=caps, mass=mass, mass_scale=mass_scale,
                endpoint=evaluation['endpoint'], endpoint_scale=evaluation['endpoint_scale'],
                cap_sum_min=float(cap_sums.min()), cap_sum_max=float(cap_sums.max()),
                cap_sum_shortage=float(max(BU_LD(0), 1-cap_sums.min())),
                structural_zero_caps=int(np.count_nonzero(~supported)),
                linear_diagnostics=evaluation['diagnostics'])


def _bu_buffers(prepared):
    rows, bins = len(prepared['a']), len(prepared['bins'])
    shape = (min(BU_BATCH, rows), bins)
    return (np.empty(shape, dtype=BU_LD), np.empty(shape, dtype=BU_LD),
            np.empty(shape, dtype=BU_LD), np.empty(rows, dtype=BU_LD))


def adaptive_budget_cache(prepared, y=2, x_values=(1, 2, 4), region_updates=16,
                          quota=None, progress=False):
    region_updates = index(region_updates)
    quota = prepared['width']*region_updates//2 if quota is None else index(quota)
    if region_updates < 2 or not 0 <= quota <= prepared['width']*region_updates or not x_values:
        raise ValueError('Invalid regional length, quota, or evaluations')
    workspace = _bu_workspace(prepared, len(x_values))
    log_C = np.log2(BU_LD(comb(prepared['width']*region_updates, quota)))
    old_logs = old_endpoint_logs = None
    maps = []
    for x in x_values:
        evaluation = _bu_dependency('adaptive_bins_evaluation')(prepared, y, x, progress, True)
        matrix, scale = evaluation['matrix'], evaluation['scale']
        region, region_scale = _bu_dependency('_ab_power')(matrix, scale, region_updates)
        endpoint, endpoint_scale = _bu_dependency('_ab_apply_power')(
            matrix, scale, evaluation['endpoint'], evaluation['endpoint_scale'], region_updates-1)
        correction = -quota*np.log2(BU_LD(x))-log_C
        local = _bu_dependency('_ab_logs')(region, region_scale+correction)
        local_end = _bu_dependency('_ab_logs')(endpoint, endpoint_scale+correction)
        if old_logs is None:
            old_logs, old_endpoint_logs = local, local_end
        else:
            np.minimum(old_logs, local, out=old_logs)
            np.minimum(old_endpoint_logs, local_end, out=old_endpoint_logs)
        maps.append(_bu_map_data(prepared, evaluation, x))
        if progress:
            print(dict(stage='adaptive_budget_cache', x=float(x),
                       cap_sum_min=maps[-1]['cap_sum_min'], cap_sum_max=maps[-1]['cap_sum_max']), flush=True)
    ordinary, ordinary_scale = _bu_dependency('_ab_from_logs')(old_logs)
    endpoint, endpoint_scale = _bu_dependency('_ab_from_logs')(old_endpoint_logs)
    return dict(prepared=prepared, maps=tuple(maps), y=BU_LD(y), quota=quota,
                region_updates=region_updates, log_C=log_C,
                ordinary=ordinary, ordinary_scale=ordinary_scale,
                endpoint=endpoint, endpoint_scale=endpoint_scale,
                workspace=_bu_buffers(prepared),
                estimated_workspace_bytes=workspace)


def adaptive_budget_apply(prepared, row_map, vector, scale=0, workspace=None):
    """One homogeneous N_x action, sorted once and evaluated in contiguous batches."""
    vector = np.asarray(vector, dtype=BU_LD)
    if vector.shape != (len(prepared['bins']),) or np.any(vector <= 0) or not np.all(np.isfinite(vector)):
        raise ValueError('Use a finite strictly positive bin vector')
    order = np.argsort(-vector, kind='stable')
    values = vector[order]
    rows, bins = row_map['caps'].shape
    capacity = min(BU_BATCH, rows)
    sorted_caps, prefix, allocation, physical = _bu_buffers(prepared) if workspace is None else workspace
    for start in range(0, rows, capacity):
        stop = min(rows, start+capacity)
        count = stop-start
        caps, cumulative, assigned = sorted_caps[:count], prefix[:count], allocation[:count]
        np.take(row_map['caps'][start:stop], order, axis=1, out=caps)
        np.cumsum(caps, axis=1, out=cumulative)
        assigned[:, 0] = 1
        np.subtract(1, cumulative[:, :-1], out=assigned[:, 1:])
        np.maximum(assigned, 0, out=assigned)  # The exact remaining unit budget.
        np.minimum(caps, assigned, out=assigned)
        positive = assigned[assigned > 0]
        if not len(positive) or positive.min()*values.min() <= 0:
            raise ArithmeticError('A positive fractional-knapsack product would underflow')
        if np.any(assigned < 0) or np.any(assigned > caps):
            raise ArithmeticError('Budget allocation left its cap box')
        if np.any(assigned.sum(axis=1, dtype=BU_LD) > 1+BU_LD('1e-12')):
            raise ArithmeticError('Floating budget allocation exceeded one')
        np.matmul(assigned, values, out=physical[start:stop])
    np.multiply(physical, row_map['mass'], out=physical)
    if np.any(physical <= 0) or not np.all(np.isfinite(physical)):
        raise ArithmeticError('Positive physical budget action underflow/overflow')
    result = np.asarray([physical[rows].max() for rows in prepared['bins']], dtype=BU_LD)
    result, extra = _bu_dependency('_ab_normalize')(result)
    return result, BU_LD(scale)+row_map['mass_scale']+extra


def adaptive_budget_region(cache, vector, scale=0, progress=False):
    """Valid regional vector bound, intersected with old Mbar v."""
    linear, linear_scale = _bu_dependency('_ab_product')(
        cache['ordinary'], cache['ordinary_scale'], vector, scale)
    old_logs = _bu_dependency('_ab_logs')(linear, linear_scale)
    best = old_logs.copy()
    selectors = np.full(len(vector), -1, dtype=np.int32)
    for j, row_map in enumerate(cache['maps']):
        current, current_scale = vector, BU_LD(scale)
        for _ in range(cache['region_updates']):
            current, current_scale = adaptive_budget_apply(cache['prepared'], row_map, current, current_scale, cache.get('workspace'))
        correction = -cache['quota']*np.log2(row_map['x'])-cache['log_C']
        logs = _bu_dependency('_ab_logs')(current, current_scale+correction)
        selectors[logs < best] = j
        np.minimum(best, logs, out=best)
        if progress:
            print(dict(stage='adaptive_budget_region', x=float(row_map['x']),
                       maximum_output_log2=float(logs.max())), flush=True)
    result, result_scale = _bu_dependency('_ab_from_logs')(best)
    return dict(vector=result, scale=result_scale,
                selector_counts=tuple(int(np.count_nonzero(selectors == j)) for j in range(-1, len(cache['maps']))),
                maximum_log2_gain_over_old_linear=float((old_logs-best).max()))


def adaptive_budget_endpoint(cache, progress=False):
    best = _bu_dependency('_ab_logs')(cache['endpoint'], cache['endpoint_scale'])
    old = best.copy()
    for row_map in cache['maps']:
        current, scale = row_map['endpoint'], row_map['endpoint_scale']
        for _ in range(cache['region_updates']-1):
            current, scale = adaptive_budget_apply(cache['prepared'], row_map, current, scale, cache.get('workspace'))
        correction = -cache['quota']*np.log2(row_map['x'])-cache['log_C']
        logs = _bu_dependency('_ab_logs')(current, scale+correction)
        np.minimum(best, logs, out=best)
        if progress:
            print(dict(stage='adaptive_budget_endpoint', x=float(row_map['x']),
                       maximum_output_log2=float(logs.max())), flush=True)
    result, scale = _bu_dependency('_ab_from_logs')(best)
    return dict(vector=result, scale=scale, maximum_log2_gain_over_old_linear=float((old-best).max()))


def adaptive_budget_bound(iterations=2, progress=False):
    """One fixed beta0=2,46-bin control; no marker or lift search."""
    iterations = index(iterations)
    if not 1 <= iterations <= 2:
        raise ValueError('This bounded control permits one or two witness evaluations')
    started = time.monotonic()
    prepared = _bu_dependency('adaptive_bins_prepare')(512, 64, 1, 2, (-1/2, 0, 1/2))
    cache = adaptive_budget_cache(prepared, 2, (1, 2, 4), 16, 4096, progress)
    prefix = cache['log_C']-60518+prepared['log2_lift'][-1]
    old, old_scale = _bu_dependency('_ab_apply_power')(
        cache['ordinary'], cache['ordinary_scale'], cache['endpoint'], cache['endpoint_scale'], 255)
    baseline = prefix+old_scale+np.log2(old[prepared['full_bin']])
    endpoint = adaptive_budget_endpoint(cache, progress)
    endpoint_logs = _bu_dependency('_ab_logs')(endpoint['vector'], endpoint['scale'])
    vector = np.ones(len(prepared['bins']), dtype=BU_LD)
    for _ in range(64):
        vector, _ = _bu_dependency('_ab_normalize')(cache['ordinary'] @ vector)
    best, source, results = baseline, 'old_linear_full_power', []
    for iteration in range(iterations):
        step_started = time.monotonic()
        regional = adaptive_budget_region(cache, vector, 0, progress)
        logs = _bu_dependency('_ab_logs')(regional['vector'], regional['scale'])
        log_vector = np.log2(vector)
        rho, mu = (logs-log_vector).max(), (endpoint_logs-log_vector).max()
        candidate = prefix+log_vector[prepared['full_bin']]+255*rho+mu
        if candidate < best:
            best, source = candidate, 'nonlinear_witness_'+str(iteration)
        results.append(dict(iteration=iteration, log2_probability_upper=float(candidate),
                            log2_rho=float(rho), log2_mu=float(mu),
                            regional_selector_counts=regional['selector_counts'],
                            maximum_region_log2_gain=regional['maximum_log2_gain_over_old_linear'],
                            seconds=time.monotonic()-step_started))
        vector = regional['vector']  # A new positive witness, with arbitrary common scale removed.
        if progress:
            print(dict(stage='adaptive_budget_witness', **results[-1]), flush=True)
    return dict(y=2, threshold=60518, alpha0=1, beta0=2, bin_width=64,
                bin_count=len(prepared['bins']), physical_rows=len(prepared['a']),
                x_values=(1, 2, 4), tilt_offsets=(-1/2, 0, 1/2),
                old_linear_full_power_log2_upper=float(baseline),
                best_log2_probability_upper=float(best), best_bound_source=source,
                witnesses=results, endpoint_maximum_log2_gain=endpoint['maximum_log2_gain_over_old_linear'],
                cap_diagnostics=[{key: row[key] for key in ('cap_sum_min', 'cap_sum_max', 'cap_sum_shortage', 'structural_zero_caps')} for row in cache['maps']],
                estimated_workspace_bytes=cache['estimated_workspace_bytes'], seconds=time.monotonic()-started,
                caps_renormalized=False, compose_before_quota_minimum=True, H_union_count=1,
                extra_y_to_width_factor=False, nonlinear_comparison_not_physical_distribution=True,
                floating_diagnostic_only=True, full_code_certificate=False)


def _bu_fraction_fill(caps, vector):
    remaining, value = Fraction(1), Fraction(0)
    for j in sorted(range(len(vector)), key=lambda j: vector[j], reverse=True):
        take = min(caps[j], remaining)
        value += take*vector[j]
        remaining -= take
    return value


def _bu_vertices(caps, vector):
    """Independent vertex enumeration of the capped at-most-unit simplex."""
    best, size = Fraction(0), len(caps)
    for bits in product((0, 1), repeat=size):
        point = [caps[j]*bits[j] for j in range(size)]
        if sum(point) <= 1:
            best = max(best, sum(p*v for p, v in zip(point, vector)))
    for free in range(size):
        others = [j for j in range(size) if j != free]
        for bits in product((0, 1), repeat=size-1):
            point = [Fraction(0)]*size
            for j, bit in zip(others, bits):
                point[j] = caps[j]*bit
            point[free] = 1-sum(point)
            if 0 <= point[free] <= caps[free]:
                best = max(best, sum(p*v for p, v in zip(point, vector)))
    return best


def _bu_exact_apply(prepared, caps, mass, vector):
    physical = [t*_bu_fraction_fill(c, vector) for c, t in zip(caps, mass)]
    return [max(physical[int(s)] for s in rows) for rows in prepared['bins']]


def _bu_float_fraction_array(values):
    return np.asarray([BU_LD(value.numerator)/value.denominator for value in values], dtype=BU_LD)


def adaptive_budget_checks():
    counters = dict(knapsack_vertices=0, physical_rows=0, linear_dominance=0,
                    iterated_composition=0, quota=0, endpoints=0, collatz=0, floating=0)
    for caps in ((Fraction(1, 2),)*3, (Fraction(0), Fraction(1), Fraction(1, 3)),
                 (Fraction(1, 4),)*3, (Fraction(2, 3), Fraction(1, 2), Fraction(1, 5)),
                 (Fraction(1), Fraction(1), Fraction(1))):
        for values in product((Fraction(0), Fraction(1), Fraction(3)), repeat=3):
            assert _bu_fraction_fill(caps, values) == _bu_vertices(caps, values)
            assert _bu_fraction_fill(caps, [2*v for v in values]) == 2*_bu_fraction_fill(caps, values)
            counters['knapsack_vertices'] += 1
    mm, mv = _bu_dependency('_ab_mm'), _bu_dependency('_ab_mv')
    for width in (4, 8):
        prepared = _bu_dependency('adaptive_bins_prepare')(width, 4, 1, 2, (-1, 0, 1))
        states, kernel = _bu_dependency('_ab_exact_kernel')(width, Fraction(2))
        size, bins, quota = len(states), prepared['bins'], width
        C = Fraction(comb(2*width, width))
        B = [[sum(kernel[i][k][d]*kernel[k][j][quota-d] for k in range(size) for d in range(width+1))/C
              for j in range(size)] for i in range(size)]
        qpoly = [[sum(Fraction(comb(w, r)*comb(width-w, j-r))*2**(w-r)
                      for r in range(max(0, j-(width-w)), min(w, j)+1)) for j in range(width+1)] for a, w in states]
        q = [sum(kernel[i][k][d]*qpoly[k][quota-d] for k in range(size) for d in range(width+1))/C for i in range(size)]
        vector = [Fraction(j+2, j+1) for j in range(len(bins))]
        regional = end = old_region = old_end = H = lift = None
        float_maps = []
        for x in (Fraction(1), Fraction(2)):
            matrix = [[sum(value*x**j for j, value in enumerate(poly)) for poly in row] for row in kernel]
            comparison, row_bounds, lift = _bu_dependency('_ab_exact_bin_matrix')(prepared, matrix, Fraction(1), Fraction(2))
            H = [[lift[s] if int(prepared['bin_ids'][s]) == j else Fraction(0) for j in range(len(bins))] for s in range(size)]
            KH = mm(matrix, H)
            mass = [value/lift[s] for s, value in enumerate(mv(matrix, lift))]
            caps = [[min(Fraction(1), bound/mass[s]) for bound in row] for s, row in enumerate(row_bounds)]
            for s in range(size):
                assert sum(caps[s]) >= 1
                probability = [value/(lift[s]*mass[s]) for value in KH[s]]
                assert sum(probability) == 1 and all(p <= c for p, c in zip(probability, caps[s]))
                counters['physical_rows'] += 1
            once = _bu_exact_apply(prepared, caps, mass, vector)
            twice = _bu_exact_apply(prepared, caps, mass, once)
            assert all(a <= b for a, b in zip(mv(matrix, mv(H, vector)), mv(H, once)))
            assert all(a <= b for a, b in zip(once, mv(comparison, vector)))
            assert all(a <= b for a, b in zip(mv(mm(matrix, matrix), mv(H, vector)), mv(H, twice)))
            counters['linear_dominance'] += len(bins)
            counters['iterated_composition'] += size
            qx = [(2+x)**w*(1+x)**(width-w) for a, w in states]
            ex = [max(qx[int(s)]/lift[int(s)] for s in rows) for rows in bins]
            local = [value/(x**quota*C) for value in twice]
            local_end = [value/(x**quota*C) for value in _bu_exact_apply(prepared, caps, mass, ex)]
            old_local = [[value/(x**quota*C) for value in row] for row in mm(comparison, comparison)]
            old_local_end = [value/(x**quota*C) for value in mv(comparison, ex)]
            if regional is None:
                regional, end, old_region, old_end = local, local_end, old_local, old_local_end
            else:
                regional = [min(a, b) for a, b in zip(regional, local)]
                end = [min(a, b) for a, b in zip(end, local_end)]
                old_region = [[min(a, b) for a, b in zip(first, second)] for first, second in zip(old_region, old_local)]
                old_end = [min(a, b) for a, b in zip(old_end, old_local_end)]
            evaluation = _bu_dependency('adaptive_bins_evaluation')(prepared, 2, float(x), False, True)
            float_map = _bu_map_data(prepared, evaluation, x)
            exact_caps = np.asarray([_bu_float_fraction_array(row) for row in caps])
            assert np.array_equal(float_map['caps'] == 0, exact_caps == 0)
            positive = exact_caps > 0
            assert np.max(np.abs(float_map['caps'][positive]/exact_caps[positive]-1)) < BU_LD('1e-10')
            assert np.max(np.abs(float_map['mass']*np.exp2(float_map['mass_scale'])/_bu_float_fraction_array(mass)-1)) < BU_LD('1e-10')
            float_maps.append(float_map)
            found, scale = adaptive_budget_apply(prepared, float_map, _bu_float_fraction_array(vector))
            assert np.max(np.abs(found*np.exp2(scale)/_bu_float_fraction_array(once)-1)) < BU_LD('1e-11')
            counters['floating'] += len(bins)
        regional = [min(a, b) for a, b in zip(regional, mv(old_region, vector))]
        end = [min(a, b) for a, b in zip(end, old_end)]
        assert all(a <= b for a, b in zip(mv(B, mv(H, vector)), mv(H, regional)))
        assert all(a <= b for a, b in zip(q, mv(H, end)))
        counters['quota'] += size
        counters['endpoints'] += size
        rho, mu = max(a/b for a, b in zip(regional, vector)), max(a/b for a, b in zip(end, vector))
        assert mv(mm(B, B), q)[-1] <= lift[-1]*vector[prepared['full_bin']]*rho**2*mu
        counters['collatz'] += 1
        cache = dict(prepared=prepared, maps=float_maps, region_updates=2, quota=quota,
                     log_C=np.log2(BU_LD(C.numerator)),
                     ordinary=np.asarray([_bu_float_fraction_array(row) for row in old_region]), ordinary_scale=BU_LD(0),
                     endpoint=_bu_float_fraction_array(old_end), endpoint_scale=BU_LD(0))
        found = adaptive_budget_region(cache, _bu_float_fraction_array(vector))
        found_end = adaptive_budget_endpoint(cache)
        assert np.max(np.abs(found['vector']*np.exp2(found['scale'])/_bu_float_fraction_array(regional)-1)) < BU_LD('1e-10')
        assert np.max(np.abs(found_end['vector']*np.exp2(found_end['scale'])/_bu_float_fraction_array(end)-1)) < BU_LD('1e-10')
        counters['floating'] += 2*len(bins)
    try:
        _bu_workspace(prepared, 10**9)
    except ValueError as error:
        assert 'workspace' in str(error)
    else:
        raise AssertionError('Oversized nonlinear cache was not rejected')
    print(dict(**counters, workspace_guard_checks=1, caps_never_renormalized=True, exact_width8_partitions=35,
               full_code_certificate=False), flush=True)


def adaptive_budget_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--iterations', type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            adaptive_budget_checks()
        else:
            print(adaptive_budget_bound(args.iterations, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_budget_main()
