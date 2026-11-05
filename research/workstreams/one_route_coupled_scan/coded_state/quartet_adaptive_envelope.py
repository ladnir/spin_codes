"""A bounded minimum-of-exponentials comparison on every physical population.

Fix common positive modes phi_l(a,w)=alpha_l^a beta_l^w for every quota
evaluation x. For V_c=min_l c_l phi_l, positivity gives
 K_x V_c <= F_x(c)=min_l c_l K_x phi_l.
Set d_l=max_s F_x(c)(s)/phi_l(s). Then F_x(c)<=V_d, and d is the least
coefficient vector whose fixed-slope minimum majorizes F_x(c). These
are physical function minima, not minima in an overlapping linear lift.

Compose16 such maps at fixed x, then set dbar_l=min_x x^-h d_x,l/C.
The min over x commutes with the physical min over modes. The input-one
endpoint is first enveloped, then propagated15 times at the same x.
Physical-row ratios rho=max V_dbar/V_c and mu=max V_endpoint/V_c give
 C*y^-threshold*V_c(full)*rho^255*mu.

All33,153 safe rows are retained; there is no source-bin aggregation. The
two possible losses are K(min)<=min(K), and fixed-slope re-enveloping.
Only the latter is measured by the projection diagnostics. Log envelopes
are concave piecewise affine; this family need not contain a good witness.
The old15-type coefficient-max obstruction does not apply to these exact
partition-averaged actions. Failed controls are not general no-go results.

Default: the same nine modes around(alpha0,beta0)=(1,2), offsets+-0.5,
quota x=1,2,4, marker y=2, threshold60518, and two finite witnesses.
Production is noncertifying FLOAT. No jobs or large arrays on import.
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

for _ae_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_ae_thread_var] = '1'

import numpy as np

AE_LD = np.longdouble
AE_WORKSPACE_LIMIT = 1 << 29


def _ae_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('_ab_'):
        import quartet_adaptive_bins as module
    else:
        import quartet_adaptive_population as module
    return getattr(module, name)


def adaptive_envelope_prepare(width=512, alpha0=1, beta0=2, offsets=(-1/2, 0, 1/2)):
    width = index(width)
    alpha0, beta0 = AE_LD(alpha0), AE_LD(beta0)
    if width < 4 or width > 512 or width % 4 or min(alpha0, beta0) <= 0 or not np.isfinite(alpha0+beta0):
        raise ValueError('Invalid quartet width or base potential')
    offsets = tuple(AE_LD(value) for value in offsets)
    if not offsets or 0 not in offsets or not all(np.isfinite(value) for value in offsets):
        raise ValueError('Use a finite mode family containing the central mode')
    a = np.asarray([a for a in range(0, width+1, 4) for _ in range(a+1)], dtype=np.int16)
    w = np.asarray([w for a in range(0, width+1, 4) for w in range(a+1)], dtype=np.int16)
    shifts = tuple(product(offsets, repeat=2))
    modes = np.asarray([(np.log2(alpha0)+da, np.log2(beta0)+dw) for da, dw in shifts], dtype=AE_LD)
    if len(a)*len(modes)*np.dtype(AE_LD).itemsize > AE_WORKSPACE_LIMIT//4:
        raise ValueError('Mode-basis workspace exceeds the bounded control budget')
    basis = a[:, None]*modes[:, 0]+w[:, None]*modes[:, 1]
    central = shifts.index((AE_LD(0), AE_LD(0)))
    return dict(width=width, a=a, w=w, modes=modes, log_basis=basis,
                central_mode=central, alpha0=alpha0, beta0=beta0, offsets=offsets)


def adaptive_envelope_cache(prepared, y=2, x_values=(1, 2, 4), region_updates=16,
                            quota=None, progress=False, action_callback=None, endpoint_callback=None):
    steps = index(region_updates)
    width, rows, modes = prepared['width'], len(prepared['a']), len(prepared['modes'])
    quota = width*steps//2 if quota is None else index(quota)
    if steps < 2 or not 0 <= quota <= width*steps or y <= 1 or not x_values:
        raise ValueError('Invalid marker, quota, or regional length')
    estimate = np.dtype(AE_LD).itemsize*((len(x_values)+8)*rows*modes+16*(width+1)**2)
    if estimate > AE_WORKSPACE_LIMIT:
        raise ValueError('Action-cache workspace exceeds the512MiB control budget')
    if action_callback is not None and not callable(action_callback):
        raise ValueError('The optional action callback must be callable')
    if endpoint_callback is not None and not callable(endpoint_callback):
        raise ValueError('The optional endpoint callback must be callable')
    action_function = _ae_dependency('adaptive_population_action') if action_callback is None else action_callback
    actions = np.empty((len(x_values), rows, modes), dtype=AE_LD)
    for xi, x in enumerate(x_values):
        if x <= 0 or not np.isfinite(x):
            raise ValueError('Use positive finite quota evaluations')
        for mode, (la, lb) in enumerate(prepared['modes']):
            action = action_function(y, x, np.exp2(la), np.exp2(lb), width)
            if not np.array_equal(prepared['a'], action['a']) or not np.array_equal(prepared['w'], action['w']):
                raise ValueError('An exponential action changed physical row order')
            actions[xi, :, mode] = action['log2_action']
            if progress:
                print(dict(stage='adaptive_envelope_action', x=float(x), mode=mode,
                           log2_alpha=float(la), log2_beta=float(lb),
                           seconds=round(action['seconds'], 3)), flush=True)
    if not np.all(np.isfinite(actions)):
        raise ArithmeticError('Invalid cached positive action logarithms')
    return dict(prepared=prepared, y=AE_LD(y), x_values=tuple(AE_LD(x) for x in x_values),
                actions=actions, region_updates=steps, quota=quota,
                log_C=np.log2(AE_LD(comb(width*steps, quota))),
                endpoint_callback=endpoint_callback,
                estimated_workspace_bytes=estimate)


def _ae_endpoint_logs(cache, x):
    prepared = cache['prepared']
    callback = cache.get('endpoint_callback')
    if callback is None:
        return (prepared['w']*np.log2(cache['y']+x)
                +(prepared['width']-prepared['w'])*np.log2(1+x))
    values = np.asarray(callback(prepared, cache['y'], x), dtype=AE_LD)
    if values.shape != prepared['a'].shape or not np.all(np.isfinite(values)):
        raise ValueError('Endpoint callback changed the physical row order/shape or returned invalid logs')
    return values


def _ae_normalize_coefficients(coefficients):
    if not np.all(np.isfinite(coefficients)):
        raise ArithmeticError('Invalid positive-envelope coefficient logarithms')
    scale = coefficients.max()
    return coefficients-scale, scale


def _ae_values(prepared, coefficients, scale=0):
    return np.min(prepared['log_basis']+coefficients, axis=1)+scale


def _ae_project(prepared, physical_logs):
    coefficients = np.max(physical_logs[:, None]-prepared['log_basis'], axis=0)
    envelope = _ae_values(prepared, coefficients)
    excess = envelope-physical_logs
    if excess.min() < -AE_LD('1e-8'):
        raise ArithmeticError('Envelope projection did not majorize its input')
    location = int(np.argmax(excess))
    normalized, scale = _ae_normalize_coefficients(coefficients)
    return normalized, scale, dict(maximum_projection_log2_excess=float(excess[location]),
                                   projection_max_row=(int(prepared['a'][location]), int(prepared['w'][location])))


def adaptive_envelope_step(cache, x_index, coefficients, scale=0):
    physical = np.min(cache['actions'][x_index]+coefficients, axis=1)
    result, extra, diagnostics = _ae_project(cache['prepared'], physical)
    return result, AE_LD(scale)+extra, diagnostics


def _ae_iterate(cache, x_index, coefficients, scale, steps):
    worst = dict(maximum_projection_log2_excess=0., projection_max_row=(0, 0))
    summed = 0.
    for _ in range(steps):
        coefficients, scale, diagnostic = adaptive_envelope_step(cache, x_index, coefficients, scale)
        summed += diagnostic['maximum_projection_log2_excess']
        if diagnostic['maximum_projection_log2_excess'] > worst['maximum_projection_log2_excess']:
            worst = diagnostic
    return coefficients, scale, dict(**worst, sum_step_max_projection_log2_excess=summed)


def adaptive_envelope_region(cache, coefficients, scale=0, progress=False):
    best, selectors, diagnostics = None, None, []
    for xi, x in enumerate(cache['x_values']):
        current, current_scale, diagnostic = _ae_iterate(cache, xi, coefficients.copy(), scale, cache['region_updates'])
        candidate = current+current_scale-cache['quota']*np.log2(x)-cache['log_C']
        if best is None:
            best = candidate
            selectors = np.full(len(candidate), xi, dtype=np.int32)
        else:
            selectors[candidate < best] = xi
            np.minimum(best, candidate, out=best)
        diagnostics.append(dict(x=float(x), **diagnostic))
        if progress:
            print(dict(stage='adaptive_envelope_region', **diagnostics[-1]), flush=True)
    result, result_scale = _ae_normalize_coefficients(best)
    return dict(coefficients=result, scale=result_scale, quota_x_choices_per_mode=tuple(map(int, selectors)),
                projection_diagnostics=diagnostics)


def adaptive_envelope_endpoint(cache, progress=False):
    prepared, best, selectors, diagnostics = cache['prepared'], None, None, []
    for xi, x in enumerate(cache['x_values']):
        terminal = _ae_endpoint_logs(cache, x)
        coefficients, scale, initial_projection = _ae_project(prepared, terminal)
        coefficients, scale, diagnostic = _ae_iterate(cache, xi, coefficients, scale, cache['region_updates']-1)
        candidate = coefficients+scale-cache['quota']*np.log2(x)-cache['log_C']
        if best is None:
            best = candidate
            selectors = np.full(len(candidate), xi, dtype=np.int32)
        else:
            selectors[candidate < best] = xi
            np.minimum(best, candidate, out=best)
        diagnostics.append(dict(x=float(x), initial_endpoint_projection=initial_projection, **diagnostic))
        if progress:
            print(dict(stage='adaptive_envelope_endpoint', **diagnostics[-1]), flush=True)
    coefficients, scale = _ae_normalize_coefficients(best)
    return dict(coefficients=coefficients, scale=scale, quota_x_choices_per_mode=tuple(map(int, selectors)),
                projection_diagnostics=diagnostics)


def _ae_mode_activity(prepared, coefficients):
    active = np.argmin(prepared['log_basis']+coefficients, axis=1)
    return dict(active_physical_mode_counts=tuple(int(np.count_nonzero(active == j)) for j in range(len(coefficients))),
                active_mode_empty=int(active[0]), active_mode_full=int(active[-1]))


def _ae_single_mode_controls(cache, threshold, regions):
    prepared = cache['prepared']
    one_step = np.max(cache['actions']-prepared['log_basis'][None, :, :], axis=1)
    region_logs, endpoint_logs = [], []
    for xi, x in enumerate(cache['x_values']):
        correction = -cache['quota']*np.log2(x)-cache['log_C']
        region_logs.append(cache['region_updates']*one_step[xi]+correction)
        terminal = _ae_endpoint_logs(cache, x)
        mu = np.max(terminal[:, None]-prepared['log_basis'], axis=0)
        endpoint_logs.append(mu+(cache['region_updates']-1)*one_step[xi]+correction)
    rho = np.min(np.asarray(region_logs), axis=0)
    endpoint = np.min(np.asarray(endpoint_logs), axis=0)
    bounds = (cache['log_C']-threshold*np.log2(cache['y'])+prepared['log_basis'][-1]
              +(regions-1)*rho+endpoint)
    return dict(one_update_exact_action_collatz_logs=tuple(tuple(map(float, row)) for row in one_step),
                whole_bound_by_mode=tuple(map(float, bounds)), best_mode=int(np.argmin(bounds)),
                best_bound=float(bounds.min()), _best_value=bounds.min(), not_a_family_impossibility_claim=True)


def adaptive_envelope_cached_bound(cache, threshold=60518, regions=256, iterations=2, progress=False):
    """Compare finite positive witnesses using one already fixed operator cache."""
    iterations = index(iterations)
    threshold, regions = index(threshold), index(regions)
    if not 1 <= iterations <= 2:
        raise ValueError('Use one or two bounded positive witnesses')
    if threshold < 0 or regions < 1:
        raise ValueError('Use a nonnegative threshold and positive region count')
    started = time.monotonic()
    prepared = cache['prepared']
    endpoint = adaptive_envelope_endpoint(cache, progress)
    endpoint_values = _ae_values(prepared, endpoint['coefficients'], endpoint['scale'])
    # The central mode is represented exactly by a minimum of dominating modes.
    central = prepared['log_basis'][:, prepared['central_mode']]
    initial = np.max(central[:, None]-prepared['log_basis'], axis=0)
    coefficients, _ = _ae_normalize_coefficients(initial)
    single = _ae_single_mode_controls(cache, threshold, regions)
    best, best_source = single.pop('_best_value'), 'single_mode_'+str(single['best_mode'])
    prefix = cache['log_C']-threshold*np.log2(cache['y'])
    results = []
    for iteration in range(iterations):
        current_values = _ae_values(prepared, coefficients)
        region = adaptive_envelope_region(cache, coefficients, progress=progress)
        region_values = _ae_values(prepared, region['coefficients'], region['scale'])
        ratios = region_values-current_values
        endpoint_ratios = endpoint_values-current_values
        worst = int(np.argmax(ratios))
        rho, mu = ratios[worst], endpoint_ratios.max()
        candidate = prefix+current_values[-1]+(regions-1)*rho+mu
        if candidate < best:
            best, best_source = candidate, 'envelope_witness_'+str(iteration)
        results.append(dict(iteration=iteration, log2_probability_upper=float(candidate),
                            physical_log2_rho=float(rho), physical_log2_mu=float(mu),
                            dominant_ratio_physical_row=(int(prepared['a'][worst]), int(prepared['w'][worst])),
                            quota_x_choices_per_mode=region['quota_x_choices_per_mode'],
                            projection_diagnostics=region['projection_diagnostics'],
                            **_ae_mode_activity(prepared, coefficients)))
        if progress:
            print(dict(stage='adaptive_envelope_witness', **results[-1]), flush=True)
        coefficients = region['coefficients']
    return dict(y=float(cache['y']), threshold=threshold, regions=regions,
                modes=tuple(tuple(map(float, row)) for row in prepared['modes']),
                x_values=tuple(map(float, cache['x_values'])), physical_rows=len(prepared['a']), bins_used=False,
                best_log2_probability_upper=float(best), best_bound_source=best_source,
                witnesses=results, single_mode_controls=single,
                endpoint_projection_diagnostics=endpoint['projection_diagnostics'],
                endpoint_quota_x_choices_per_mode=endpoint['quota_x_choices_per_mode'],
                endpoint_mode_activity=_ae_mode_activity(prepared, endpoint['coefficients']),
                estimated_workspace_bytes=cache['estimated_workspace_bytes'], seconds=time.monotonic()-started,
                H_union_count=1, extra_y_to_width_factor=False, physical_collatz_ratios=True,
                projection_loss_is_not_total_relaxation_loss=True,
                floating_diagnostic_only=True, full_code_certificate=False)


def adaptive_envelope_bound(iterations=2, progress=False):
    """Preserved baseline control with the original operator and endpoint."""
    if not 1 <= index(iterations) <= 2:
        raise ValueError('Use one or two bounded positive witnesses')
    started = time.monotonic()
    prepared = adaptive_envelope_prepare()
    cache = adaptive_envelope_cache(prepared, progress=progress)
    result = adaptive_envelope_cached_bound(cache, 60518, 256, iterations, progress)
    result['seconds'] = time.monotonic()-started
    return result


def _ae_exact_values(basis, coefficients):
    return [min(c*p for c, p in zip(coefficients, row)) for row in basis]


def _ae_exact_step(matrix, basis, coefficients):
    action = _ae_dependency('_ab_mm')(matrix, basis)
    physical = [min(c*p for c, p in zip(coefficients, row)) for row in action]
    projected = [max(physical[s]/basis[s][j] for s in range(len(basis))) for j in range(len(coefficients))]
    return projected, physical


def adaptive_envelope_checks():
    counters = dict(physical_minimum=0, projection=0, regional_quota=0, endpoint=0,
                    physical_collatz=0, coefficient_ratio_control=0, floating=0)
    mm, mv = _ae_dependency('_ab_mm'), _ae_dependency('_ab_mv')
    for width in (4, 8):
        prepared = adaptive_envelope_prepare(width, 1, 2, (-1, 0, 1))
        states, kernel = _ae_dependency('_ab_exact_kernel')(width, Fraction(2))
        basis = [[Fraction(2)**(int(la)*a+int(lb)*w) for la, lb in prepared['modes']] for a, w in states]
        # Two modes cross inside the physical state set; the remaining modes
        # are safely inactive. This exercises K(min)<min(K), not just one mode.
        coefficients = [Fraction(2)**(4*width) for _ in prepared['modes']]
        coefficients[0] = Fraction(2)**(width//2)
        coefficients[prepared['central_mode']] = Fraction(1)
        values = _ae_exact_values(basis, coefficients)
        size, quota, C = len(states), width, Fraction(comb(2*width, width))
        B = [[sum(kernel[i][k][d]*kernel[k][j][quota-d] for k in range(size) for d in range(width+1))/C
              for j in range(size)] for i in range(size)]
        qpoly = [[sum(Fraction(comb(w, r)*comb(width-w, j-r))*2**(w-r)
                      for r in range(max(0, j-(width-w)), min(w, j)+1)) for j in range(width+1)] for a, w in states]
        q = [sum(kernel[i][k][d]*qpoly[k][quota-d] for k in range(size) for d in range(width+1))/C for i in range(size)]
        regional = endpoint = None
        actions = []
        for x in (Fraction(1), Fraction(2)):
            matrix = [[sum(value*x**j for j, value in enumerate(poly)) for poly in row] for row in kernel]
            once, intermediate = _ae_exact_step(matrix, basis, coefficients)
            once_values = _ae_exact_values(basis, once)
            actual = mv(matrix, values)
            assert all(a <= b <= c for a, b, c in zip(actual, intermediate, once_values))
            assert any(a < b for a, b in zip(actual, intermediate))
            counters['physical_minimum'] += size
            counters['projection'] += size
            twice, _ = _ae_exact_step(matrix, basis, once)
            local = [value/(x**quota*C) for value in twice]
            qx = [(2+x)**w*(1+x)**(width-w) for a, w in states]
            end_coeff = [max(qx[s]/basis[s][j] for s in range(size)) for j in range(len(coefficients))]
            end_coeff, _ = _ae_exact_step(matrix, basis, end_coeff)
            local_end = [value/(x**quota*C) for value in end_coeff]
            if regional is None:
                regional, endpoint = local, local_end
            else:
                regional = [min(a, b) for a, b in zip(regional, local)]
                endpoint = [min(a, b) for a, b in zip(endpoint, local_end)]
            action = mm(matrix, basis)
            actions.append(np.asarray([[np.log2(AE_LD(v.numerator)/v.denominator) for v in row] for row in action]))
        regional_values, endpoint_values = _ae_exact_values(basis, regional), _ae_exact_values(basis, endpoint)
        assert all(a <= b for a, b in zip(mv(B, values), regional_values))
        assert all(a <= b for a, b in zip(q, endpoint_values))
        counters['regional_quota'] += size
        counters['endpoint'] += size
        rho = max(a/b for a, b in zip(regional_values, values))
        mu = max(a/b for a, b in zip(endpoint_values, values))
        assert mv(mm(B, B), q)[-1] <= values[-1]*rho**2*mu
        assert rho <= max(a/b for a, b in zip(regional, coefficients))
        counters['physical_collatz'] += 1
        counters['coefficient_ratio_control'] += 1
        cache = dict(prepared=prepared, actions=np.asarray(actions), x_values=(AE_LD(1), AE_LD(2)),
                     y=AE_LD(2), quota=quota, region_updates=2, log_C=np.log2(AE_LD(C.numerator)))
        logs = np.asarray([np.log2(AE_LD(c.numerator)/c.denominator) for c in coefficients])
        found = adaptive_envelope_region(cache, logs)
        found_end = adaptive_envelope_endpoint(cache)
        expected = np.asarray([np.log2(AE_LD(v.numerator)/v.denominator) for v in regional_values])
        expected_end = np.asarray([np.log2(AE_LD(v.numerator)/v.denominator) for v in endpoint_values])
        assert np.max(np.abs(_ae_values(prepared, found['coefficients'], found['scale'])-expected)) < AE_LD('1e-10')
        assert np.max(np.abs(_ae_values(prepared, found_end['coefficients'], found_end['scale'])-expected_end)) < AE_LD('1e-10')
        counters['floating'] += 2*size
    print(dict(**counters, exact_width8_partitions=35, physical_row_ratio_not_coefficient_ratio=True,
               fixed_modes_across_quota_evaluations=True, full_code_certificate=False), flush=True)


def adaptive_envelope_main():
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
            adaptive_envelope_checks()
        else:
            print(adaptive_envelope_bound(args.iterations, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_envelope_main()
