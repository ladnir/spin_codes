"""One cached full-population regional witness for a general proof network.

The default is the ORIGINAL network(1,1,2), marker y=3/2, threshold60518,
and quota evaluations7/5,33/20,9/5. The encoder and setup are unchanged.
The positive potential has log2 phi(a,w)=(236a+232w)/512.

For each fixed x, compute K_x^16 phi and save K_x^15 phi. Only after
composition take physical, common-state minima of x^-h K_x^16 phi/C0.
The general free endpoint is
 q_x(a,w)=(1+x)^(n-a)(y^(source-temporal)+x)^(a-w)(y^source+x)^w.
With mu_x=max(q_x/phi), the endpoint is bounded by the corresponding
minimum of x^-h mu_x K_x^15 phi/C0. Thus the geometry-tail comparison is
 C0*y^-threshold*phi(full)*rho^(regions-1)*nu.
There is one row-set union and NO additional y^n endpoint factor.

This caller requires input capacity=2*source and half regional quotas;
these hypotheses remove the balanced-label offset. It supports one
witness only. It imports the unchanged Horner core and T-power cache;
no bins, amplitude envelope, quota-profile projection, or monkeypatching.
Every physical population remains present. Arithmetic is FLOAT-only,
even though the polynomial action identity is exact. No full certificate
or failure probability lower bound follows from an unsuccessful control.
Parent owns all numerical execution; import starts no job.
"""

import argparse
from fractions import Fraction as F
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _gr_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_gr_thread_var] = '1'

import numpy as np

GR_LD = np.longdouble
GR_XS = (F(7, 5), F(33, 20), F(9, 5))


def _gr_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('_hr_'):
        import quartet_horner_region as module
    elif name == 'adaptive_vector_cache_prepare':
        import quartet_adaptive_vector_cache as module
    elif name.startswith('general_flow') or name == '_pt_capacities':
        import quartet_partial_temporal as module
    elif name == '_ap_partitions':
        import quartet_adaptive_population as module
    else:
        import quartet_adaptive_vector as module
    return getattr(module, name)


def _gr_ld(value):
    value = F(value)
    return GR_LD(value.numerator)/GR_LD(value.denominator)


def _gr_cached_action(core, log_values):
    prepared, powers = core
    return _gr_dependency('adaptive_vector_action')(
        prepared, log_values, log_scale=0, progress=False, power_cache=powers)


def general_horner_region_prepare(y=F(3, 2), capacities=(1, 1, 2), width=512,
                                  region_updates=16, x_values=GR_XS,
                                  deadline=None, progress=False):
    """Prepare common rows and three independent, input-free T caches."""
    y, width, steps = F(y), index(width), index(region_updates)
    capacities = _gr_dependency('_pt_capacities')(capacities)
    xs = tuple(F(x) for x in x_values)
    if (y <= 1 or capacities[2] != 2*capacities[0] or width < 4 or width > 512
            or width % 4 or steps < 2 or steps > 16 or not 1 <= len(xs) <= 3
            or any(x <= 0 for x in xs)):
        raise ValueError('Require balanced c=2a, y>1, quartet width4..512, steps2..16, and1..3 positive x')
    prepare = _gr_dependency('adaptive_vector_prepare')
    cache_powers = _gr_dependency('adaptive_vector_cache_prepare')
    deadline_check = _gr_dependency('_hr_deadline')
    a = np.asarray([a for a in range(0, width+1, 4) for _ in range(a+1)], dtype=np.int16)
    w = np.asarray([w for a in range(0, width+1, 4) for w in range(a+1)], dtype=np.int16)
    cores, retained = [], 0
    for x in xs:
        deadline_check(deadline)
        started = time.monotonic()
        core = prepare(_gr_ld(y), _gr_ld(x), width, capacities)
        if not np.array_equal(core['a'], a) or not np.array_equal(core['w'], w):
            raise ValueError('General core changed physical row ordering')
        powers = cache_powers(core, retained_cache_bytes=retained,
                              workspace_limit_bytes=1 << 30, progress=False)
        retained += powers['cache_bytes']
        cores.append((core, powers))
        if progress:
            print(dict(stage='general_horner_prepare', x=str(x),
                       retained_cache_MiB=round(retained/(1 << 20), 3),
                       seconds=round(time.monotonic()-started, 3)), flush=True)
    quota = width*steps//2
    log_C = np.log2(GR_LD(comb(width*steps, quota)))
    if not np.isfinite(log_C):
        raise ArithmeticError('Quota binomial logarithm overflow')
    return dict(width=width, steps=steps, quota=quota, log_C=log_C,
                capacities=capacities, y=_gr_ld(y), rational_y=y,
                x_values=tuple(map(_gr_ld, xs)), rational_x_values=xs,
                prepared=cores, action=_gr_cached_action, a=a, w=w,
                full=len(a)-1, retained_cache_bytes=retained)


def general_horner_region_witness(cache, log_phi, regions=256, threshold=60518,
                                  deadline=None, progress=False):
    """One fixed physical potential, with all quota minima after composition."""
    started = time.monotonic()
    regions, threshold = index(regions), index(threshold)
    if regions < 1 or threshold < 0:
        raise ValueError('Require positive regions and a nonnegative integer threshold')
    finite = _gr_dependency('_hr_finite_logs')
    phi = finite(log_phi, cache['a'].shape, 'general potential')
    F_upper = np.full(phi.shape, GR_LD('inf'))
    endpoint = np.full(phi.shape, GR_LD('inf'))
    choices = np.zeros(phi.shape, dtype=np.int16)
    summaries = []
    for xi, (x, core) in enumerate(zip(cache['x_values'], cache['prepared'])):
        x_started = time.monotonic()
        saved, final = _gr_dependency('_hr_chain')(cache, core, phi, deadline)
        correction = cache['quota']*np.log2(x)+cache['log_C']
        candidate = final-correction
        choices[candidate < F_upper] = xi
        np.minimum(F_upper, candidate, out=F_upper)
        q = _gr_dependency('general_flow_endpoint_logs')(
            cache, cache['y'], x, cache['capacities'])
        q = finite(q, phi.shape, 'general endpoint')
        mu = (q-phi).max()
        local_end = saved+mu-correction
        np.minimum(endpoint, local_end, out=endpoint)
        row = int(np.argmax(candidate-phi))
        summary = dict(x=str(cache['rational_x_values'][xi]),
                       seconds=round(time.monotonic()-x_started, 3),
                       log2_bulk_ratio=float((candidate-phi)[row]),
                       bulk_row=(int(cache['a'][row]), int(cache['w'][row])),
                       log2_endpoint_multiplier=float(mu),
                       log2_endpoint_ratio=float((local_end-phi).max()))
        summaries.append(summary)
        if progress:
            print(dict(stage='general_horner_x', **summary), flush=True)
    F_upper = finite(F_upper, phi.shape, 'general regional upper')
    endpoint = finite(endpoint, phi.shape, 'general endpoint upper')
    ratios, ends = F_upper-phi, endpoint-phi
    row, low, endrow = int(np.argmax(ratios)), int(np.argmin(ratios)), int(np.argmax(ends))
    rho, nu = ratios[row], ends[endrow]
    tail = (cache['log_C']-threshold*np.log2(cache['y'])+phi[cache['full']]
            +(regions-1)*rho+nu)
    summary = dict(seconds=round(time.monotonic()-started, 3),
                   capacities=cache['capacities'], y=str(cache['rational_y']),
                   threshold=threshold, rows=len(phi), regions=regions,
                   updates_per_region=cache['steps'], quota=cache['quota'],
                   log2_rho=float(rho), log2_nu=float(nu),
                   log2_phi_full=float(phi[cache['full']]), log2_tail_upper=float(tail),
                   rho_row=(int(cache['a'][row]), int(cache['w'][row])),
                   endpoint_row=(int(cache['a'][endrow]), int(cache['w'][endrow])),
                   log2_comparison_lower_floor=float(ratios[low]),
                   comparison_lower_row=(int(cache['a'][low]), int(cache['w'][low])),
                   lower_floor_is_for_this_nonlinear_comparison_only=True,
                   quota_choice_counts=[int(np.count_nonzero(choices == i)) for i in range(len(summaries))],
                   cached_T_powers=True, exact_physical_action_algebra=True,
                   floating_diagnostic_only=True, full_code_certificate=False)
    return dict(summary=summary, per_x=summaries, log2_F=F_upper,
                log2_endpoint_upper=endpoint, log2_tail_upper=tail)


def general_horner_region_control(y=F(3, 2), capacities=(1, 1, 2), threshold=60518,
                                  x_values=GR_XS, log_p=236, log_h=468,
                                  deadline=None, progress=True):
    started = time.monotonic()
    cache = general_horner_region_prepare(y, capacities, x_values=x_values,
                                          deadline=deadline, progress=progress)
    log_p, log_h = GR_LD(log_p), GR_LD(log_h)
    if not np.isfinite(log_p+log_h):
        raise ValueError('Potential pure-state logs must be finite')
    phi = (log_p*cache['a'].astype(GR_LD)+(log_h-log_p)*cache['w'])/cache['width']
    result = general_horner_region_witness(cache, phi, threshold=threshold,
                                          deadline=deadline, progress=progress)
    summary = dict(result['summary'])
    summary.update(stage='general_horner_complete', seconds=round(time.monotonic()-started, 3),
                   pure_state_log_potential=(0, float(log_p), float(log_h)),
                   retained_cache_MiB=round(cache['retained_cache_bytes']/(1 << 20), 3),
                   per_x=result['per_x'])
    return summary


def _gr_exact_kernel(width, y, capacities):
    """Tiny physical label/coin/partition enumeration with a quota counter."""
    states = [(a, w) for a in range(0, width+1, 4) for w in range(a+1)]
    lookup = {state: i for i, state in enumerate(states)}
    partitions = _gr_dependency('_ap_partitions')(width)
    weights = _gr_dependency('general_flow_weights')
    source_capacity, temporal, _ = capacities
    kernel = [[[F(0)]*(width+1) for _ in states] for _ in states]
    local_cache = {}
    for source, (a, w) in enumerate(states):
        A, W = (1 << a)-1, (1 << w)-1
        for partition in partitions:
            running = {(0, 0, 0): F(1)}
            for quartet in partition:
                aa = sum(((A >> v) & 1) << i for i, v in enumerate(quartet))
                ww = sum(((W >> v) & 1) << i for i, v in enumerate(quartet))
                key = aa, ww
                if key not in local_cache:
                    local = {}
                    for J in range(16):
                        g, p, _, _, _ = weights(y, capacities, aa.bit_count(), ww.bit_count(),
                                                J.bit_count(), (J & aa).bit_count(), (J & ww).bit_count())
                        cost = ((source_capacity-temporal)*((15 ^ J) & aa).bit_count()
                                +temporal*((15 ^ J) & ww).bit_count())
                        for hole, probability in ((False, p), (True, g)):
                            target = (4 if hole else 0, aa.bit_count() if hole else 0, J.bit_count())
                            local[target] = local.get(target, F(0))+probability*y**cost
                    local_cache[key] = local
                after = {}
                for (a0, w0, j0), value in running.items():
                    for (da, dw, dj), factor in local_cache[key].items():
                        target = a0+da, w0+dw, j0+dj
                        after[target] = after.get(target, F(0))+value*factor
                running = after
            for (aa, ww, j), value in running.items():
                kernel[source][lookup[aa, ww]][j] += value/len(partitions)
    return states, kernel


def general_horner_region_checks():
    counts = dict(physical_coefficients=0, quota_rows=0, endpoint_rows=0,
                  cached_chain_rows=0, unit_endpoint_rows=0, scale_rows=0, full_comparisons=0)
    mv = lambda matrix, vector: [sum(a*b for a, b in zip(row, vector)) for row in matrix]
    logs = _gr_dependency('_av_fraction_logs')
    xs = (F(7, 5), F(9, 5))
    for capacities in ((1, 1, 2), (2, 1, 4)):
        source, temporal, _ = capacities
        y = F(3, 2)
        for width in (4, 8):
            states, polys = _gr_exact_kernel(width, y, capacities)
            size, quota, C = len(states), width, F(comb(2*width, width))
            # K_x*1=q_x coefficientwise; no final reset or extra y^n.
            qpoly = [[sum(polys[i][j][d] for j in range(size)) for d in range(width+1)]
                     for i in range(size)]
            B = [[sum(polys[i][k][d]*polys[k][j][quota-d]
                      for k in range(size) for d in range(width+1))/C
                  for j in range(size)] for i in range(size)]
            b = [sum(polys[i][k][d]*qpoly[k][quota-d]
                     for k in range(size) for d in range(width+1))/C for i in range(size)]
            counts['physical_coefficients'] += size*size*(width+1)
            cache = general_horner_region_prepare(y, capacities, width, 2, xs)
            phi = [F(i*i+3, i+2) for i in range(size)]
            upper, endpoint = None, None
            for x in xs:
                matrix = [[sum(v*x**d for d, v in enumerate(poly)) for poly in row] for row in polys]
                q = [(1+x)**(width-a)*(y**(source-temporal)+x)**(a-w)*(y**source+x)**w
                     for a, w in states]
                assert mv(matrix, [F(1)]*size) == q
                counts['unit_endpoint_rows'] += size
                first, second = mv(matrix, phi), mv(matrix, mv(matrix, phi))
                mu = max(a/bb for a, bb in zip(q, phi))
                local, end = [v/(C*x**quota) for v in second], [mu*v/(C*x**quota) for v in first]
                upper = local if upper is None else list(map(min, upper, local))
                endpoint = end if endpoint is None else list(map(min, endpoint, end))
            assert all(a <= b0 for a, b0 in zip(mv(B, phi), upper))
            assert all(a <= b0 for a, b0 in zip(b, endpoint))
            counts['quota_rows'] += size
            counts['endpoint_rows'] += size
            rho, nu = max(v/p for v, p in zip(upper, phi)), max(v/p for v, p in zip(endpoint, phi))
            assert mv(B, b)[-1] <= phi[-1]*rho*nu
            counts['full_comparisons'] += 1
            result = general_horner_region_witness(cache, logs(phi), regions=2, threshold=7)
            for key, exact in (('log2_F', upper), ('log2_endpoint_upper', endpoint)):
                assert np.max(np.abs(result[key]-logs(exact))) < GR_LD('1e-9')
                counts['cached_chain_rows'] += size
            expected = np.log2(_gr_ld(C*y**(-7)*phi[-1]*rho*nu))
            assert abs(result['log2_tail_upper']-expected) < GR_LD('1e-9')
            shifted = general_horner_region_witness(cache, logs(phi)+23, regions=2, threshold=7)
            assert abs(shifted['log2_tail_upper']-result['log2_tail_upper']) < GR_LD('1e-9')
            counts['scale_rows'] += size
    print(dict(**counts, general_endpoint_checked=True, threshold_log_y_checked=True,
               two_networks=True, no_large_check_arrays=True, full_code_certificate=False), flush=True)


def general_horner_region_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--source', type=int, default=1)
    parser.add_argument('--temporal', type=int, default=1)
    parser.add_argument('--input', type=int, default=2)
    parser.add_argument('--y', default='3/2')
    parser.add_argument('--xs', default='7/5,33/20,9/5')
    parser.add_argument('--threshold', type=int, default=60518)
    parser.add_argument('--log-p', type=float, default=236.)
    parser.add_argument('--log-h', type=float, default=468.)
    parser.add_argument('--seconds', type=int, default=240)
    args = parser.parse_args()
    if not 30 <= args.seconds <= 240:
        raise ValueError('Use a bounded runtime between30 and240 seconds')
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    started = time.monotonic()
    timer = threading.Timer(args.seconds-5, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            general_horner_region_checks()
        else:
            print(general_horner_region_control(
                F(args.y), (args.source, args.temporal, args.input), args.threshold,
                tuple(F(x) for x in args.xs.split(',')), args.log_p, args.log_h,
                started+args.seconds-10, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    general_horner_region_main()
