"""One bounded destination-conditioned quota witness on all physical rows.

Use disjoint destination groups J0={w<n/2}, J1={w>=n/2}, and one
positive common potential phi. For epsilon=2^-8192 set
 f_J=phi*(1_J+epsilon*1_Jc).
Their sum is kappa*phi, kappa=1+epsilon. Positivity and coefficient
Cauchy give kappa*Bphi <=sum_J min_x x^-h K_x^16 f_J/C0.
This caller drops division by kappa, giving a weaker valid bound.

For the exact endpoint b=[x^h]K_x^15 q_x/C0, split the destination
of those15 actions by the same groups. With mu_x=max(q_x/phi),
 kappa*b <=sum_J min_x mu_x*x^-h K_x^15 f_J/C0.
The multiplier is q_x/phi, NEVER q_x/f_J. Each group may use its
own fixed quota grid. All16 or15 actions precede the quota minimum;
then the two positive group bounds are summed. There is no physical
source maximum until the final scalar Collatz ratios.

Production: y2,threshold60518,phi=2^((315a+372w)/512), low grid
(1.6,1.86,2.04), high grid(2.04,2.45,2.76), one witness only.
One audited T-power cache is retained at a time and reused for its
16-step chain. The original core and regional caller are not changed.
Arithmetic is noncertifying FLOAT; the algebraic action is exact.
No arrays/jobs/files on import. Parent owns every execution, serially.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _dq_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_dq_thread_var] = '1'

import numpy as np

DQ_LD = np.longdouble
DQ_GRIDS = (('1.6', '1.86', '2.04'), ('2.04', '2.45', '2.76'))


def _dq_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('_hr_'):
        import quartet_horner_region as module
    elif name.startswith('_ab_'):
        import quartet_adaptive_bins as module
    elif name == 'adaptive_vector_cache_prepare':
        import quartet_adaptive_vector_cache as module
    else:
        import quartet_adaptive_vector as module
    return getattr(module, name)


def _dq_core_prepare(y, x, width, capacities):
    prepared = _dq_dependency('adaptive_vector_prepare')(y, x, width, capacities)
    power_cache = _dq_dependency('adaptive_vector_cache_prepare')(
        prepared, retained_cache_bytes=0, workspace_limit_bytes=1 << 30, progress=False)
    return dict(a=prepared['a'], w=prepared['w'], prepared=prepared,
                power_cache=power_cache, cache_bytes=power_cache['cache_bytes'])


def _dq_core_action(core, logs):
    return _dq_dependency('adaptive_vector_action')(
        core['prepared'], logs, log_scale=0, progress=False, power_cache=core['power_cache'])


def _dq_logadd(left, right):
    """Positive two-term physical sum; no vanishing summand is accepted."""
    peak = np.maximum(left, right)
    first, second = np.exp2(left-peak), np.exp2(right-peak)
    if np.any(first <= 0) or np.any(second <= 0) or not np.all(np.isfinite(first+second)):
        raise ArithmeticError('A positive destination-group term would underflow/overflow')
    result = peak+np.log2(first+second)
    if not np.all(np.isfinite(result)):
        raise ArithmeticError('Destination-group log sum failed')
    return result


def destination_quota_bound(width=512, region_updates=16, regions=256, threshold=60518,
                            epsilon_bits=8192, group_x_values=DQ_GRIDS, log_phi=None,
                            prepare_callback=None, action_callback=None,
                            deadline=None, progress=False):
    """Return one witness; internal vectors are never printed by the CLI.

    Callback contract matches the old regional caller: prepare accepts
    y,x,width,capacities and returns a,w; action(core,logs) returns
    absolute log2_action with zero implicit external scale. Optional
    small dimensions/callbacks exist for independent physical tests.
    """
    started = time.monotonic()
    width, steps, regions = index(width), index(region_updates), index(regions)
    threshold, epsilon_bits = index(threshold), index(epsilon_bits)
    if width < 4 or width > 512 or width % 4 or steps < 2 or regions < 1:
        raise ValueError('Use quartet width4..512, at least two updates and positive region count')
    if threshold < 0 or epsilon_bits < 1 or len(group_x_values) != 2:
        raise ValueError('Use nonnegative threshold, positive regularization bits and two grids')
    grids = tuple(tuple(DQ_LD(x) for x in values) for values in group_x_values)
    if any(not values or any(x <= 0 or not np.isfinite(x) for x in values) for values in grids):
        raise ValueError('Each destination group needs a finite positive quota grid')
    prepare = _dq_core_prepare if prepare_callback is None else prepare_callback
    action = _dq_core_action if action_callback is None else action_callback
    if not callable(prepare) or not callable(action):
        raise ValueError('Core adapters must be callable')
    a = np.asarray([a for a in range(0, width+1, 4) for _ in range(a+1)], dtype=np.int16)
    w = np.asarray([w for a in range(0, width+1, 4) for w in range(a+1)], dtype=np.int16)
    finite_logs, check_deadline = _dq_dependency('_hr_finite_logs'), _dq_dependency('_hr_deadline')
    phi = ((315*a.astype(DQ_LD)+372*w.astype(DQ_LD))/512 if log_phi is None
           else finite_logs(log_phi, a.shape, 'destination potential'))
    masks = (w < width//2, w >= width//2)
    assert np.all(masks[0] ^ masks[1]) and all(np.any(mask) for mask in masks)
    quota = width*steps//2
    log_C = np.log2(DQ_LD(comb(width*steps, quota)))
    if not np.isfinite(log_C):
        raise ArithmeticError('Regional quota normalizer overflow')
    chain = _dq_dependency('_hr_chain')
    context = dict(width=width, steps=steps, a=a, w=w, action=action)
    aggregate, endpoint = None, None
    group_summaries = []
    max_cache_bytes = 0
    for group, (inside, xs) in enumerate(zip(masks, grids)):
        f_logs = phi.copy()
        f_logs[~inside] -= epsilon_bits
        group_bound = np.full(phi.shape, DQ_LD('inf'))
        group_end = np.full(phi.shape, DQ_LD('inf'))
        choices = np.zeros(phi.shape, dtype=np.int16)
        end_choices = np.zeros(phi.shape, dtype=np.int16)
        for xi, x in enumerate(xs):
            check_deadline(deadline)
            x_started = time.monotonic()
            # No previous core/cache remains live here. The cache is reused
            # for this chain only, so retained_cache_bytes=0 is accurate.
            core = prepare(2, x, width, (1, 1, 2))
            if not np.array_equal(core['a'], a) or not np.array_equal(core['w'], w):
                raise ValueError('Destination action changed common physical row ordering')
            max_cache_bytes = max(max_cache_bytes, int(core.get('cache_bytes', 0)))
            saved, final = chain(context, core, f_logs, deadline)
            correction = quota*np.log2(x)+log_C
            candidate = final-correction
            take = candidate < group_bound
            choices[take] = xi
            np.minimum(group_bound, candidate, out=group_bound)
            q_logs = w*np.log2(2+x)+(width-w)*np.log2(1+x)
            mu = (q_logs-phi).max()
            candidate_end = saved+mu-correction
            take_end = candidate_end < group_end
            end_choices[take_end] = xi
            np.minimum(group_end, candidate_end, out=group_end)
            if progress:
                print(dict(stage='destination_quota_x', group=group, x=float(x),
                           seconds=round(time.monotonic()-x_started, 3),
                           log2_group_bulk_ratio=float((candidate-phi).max()),
                           log2_endpoint_multiplier=float(mu),
                           log2_group_endpoint_ratio=float((candidate_end-phi).max()),
                           cache_MiB=round(int(core.get('cache_bytes', 0))/(1 << 20), 3)), flush=True)
            # In particular, the next preparation does not retain three
            # caches accidentally through a container or closure.
            del core
        group_bound = finite_logs(group_bound, phi.shape, 'complete destination bulk')
        group_end = finite_logs(group_end, phi.shape, 'complete destination endpoint')
        aggregate = group_bound if aggregate is None else _dq_logadd(aggregate, group_bound)
        endpoint = group_end if endpoint is None else _dq_logadd(endpoint, group_end)
        summary = dict(group=group, destination_states=int(np.count_nonzero(inside)),
                       x_values=[float(x) for x in xs],
                       quota_choice_counts=[int(np.count_nonzero(choices == i)) for i in range(len(xs))],
                       endpoint_choice_counts=[int(np.count_nonzero(end_choices == i)) for i in range(len(xs))])
        group_summaries.append(summary)
        if progress:
            print(dict(stage='destination_quota_group', **summary), flush=True)
    ratios, end_ratios = aggregate-phi, endpoint-phi
    row, low_row, end_row = int(np.argmax(ratios)), int(np.argmin(ratios)), int(np.argmax(end_ratios))
    rho, nu = ratios[row], end_ratios[end_row]
    tail = log_C-threshold+phi[-1]+(regions-1)*rho+nu
    summary = dict(stage='destination_quota_complete', seconds=round(time.monotonic()-started, 3),
                   rows=len(a), y=2, threshold=threshold, quota=quota, regions=regions,
                   updates_per_region=steps, epsilon_bits=epsilon_bits,
                   destination_split_w=width//2, kappa_denominator_dropped=True,
                   log2_phi_full=float(phi[-1]), log2_rho=float(rho), log2_nu=float(nu),
                   log2_comparison_lower_floor=float(ratios[low_row]),
                   lower_floor_is_for_this_comparison_only=True,
                   rho_row=(int(a[row]), int(w[row])),
                   minimum_row=(int(a[low_row]), int(w[low_row])),
                   endpoint_row=(int(a[end_row]), int(w[end_row])),
                   log2_tail_upper=float(tail), completed_witnesses=1,
                   max_live_T_caches=1, max_cache_MiB=max_cache_bytes/(1 << 20),
                   exact_physical_action_algebra=True, floating_diagnostic_only=True,
                   full_code_certificate=False)
    return dict(summary=summary, groups=group_summaries,
                log2_bulk_upper=aggregate, log2_endpoint_upper=endpoint, log2_tail_upper=tail)


def _dq_fraction_logs(values):
    return np.asarray([np.log2(DQ_LD(v.numerator)/DQ_LD(v.denominator)) for v in values], dtype=DQ_LD)


def destination_quota_checks():
    counts = dict(physical_coefficients=0, grouped_quota=0, direct_endpoint=0,
                  shortcut_endpoint=0, regularization_identity=0, kappa_omission=0,
                  composition=0, floating_rows=0, scale_rows=0, quota_order=0)
    mv = _dq_dependency('_ab_mv')
    for width in (4, 8):
        states, polys = _dq_dependency('_ab_exact_kernel')(width, Fraction(2))
        size, quota = len(states), width
        C, epsilon = Fraction(comb(2*width, width)), Fraction(1, 8)
        kappa = 1+epsilon
        B = [[sum(polys[i][k][d]*polys[k][j][quota-d]
                  for k in range(size) for d in range(width+1))/C
              for j in range(size)] for i in range(size)]
        counts['physical_coefficients'] += size*size
        q_polys = [[sum(Fraction(comb(w, r)*comb(width-w, j-r))*2**(w-r)
                        for r in range(max(0, j-(width-w)), min(w, j)+1))
                    for j in range(width+1)] for a, w in states]
        b = [sum(polys[i][k][d]*q_polys[k][quota-d]
                 for k in range(size) for d in range(width+1))/C for i in range(size)]
        xs_by_group = ((Fraction(1), Fraction(2)), (Fraction(2), Fraction(3)))
        matrices = {x: [[sum(value*x**j for j, value in enumerate(poly)) for poly in row]
                        for row in polys] for x in (Fraction(1), Fraction(2), Fraction(3))}
        masks = ([w < width//2 for a, w in states], [w >= width//2 for a, w in states])
        trace = []

        def toy_prepare(y, x, n, capacities):
            assert y == 2 and n == width and capacities == (1, 1, 2)
            return dict(a=np.asarray([a for a, w in states]),
                        w=np.asarray([w for a, w in states]), x=int(x),
                        matrix=np.asarray([[DQ_LD(v.numerator)/v.denominator for v in row]
                                           for row in matrices[Fraction(int(x))]], dtype=DQ_LD))

        def toy_action(core, logs):
            trace.append(core['x'])
            peak = logs.max()
            values = core['matrix'] @ np.exp2(logs-peak)
            assert np.all(values > 0)
            return dict(log2_action=np.log2(values)+peak, a=core['a'], w=core['w'])

        for phi in ([Fraction(2)**a*Fraction(3, 2)**w for a, w in states],
                    [Fraction(i*i+3, i+2) for i in range(size)]):
            total, total_end = [Fraction(0)]*size, [Fraction(0)]*size
            group_inputs, regular_endpoints = [], []
            for inside, xs in zip(masks, xs_by_group):
                mask = [Fraction(1) if yes else epsilon for yes in inside]
                f = [p*r for p, r in zip(phi, mask)]
                group_inputs.append(f)
                true_end = [sum(polys[i][k][d]*q_polys[k][quota-d]*mask[k]
                                for k in range(size) for d in range(width+1))/C for i in range(size)]
                regular_endpoints.append(true_end)
                bound, end_bound = None, None
                for x in xs:
                    matrix = matrices[x]
                    first = mv(matrix, f)
                    candidate = [v/(x**quota*C) for v in mv(matrix, first)]
                    qx = [(2+x)**w*(1+x)**(width-w) for a, w in states]
                    mu = max(q/p for q, p in zip(qx, phi))
                    direct = [v/(x**quota*C) for v in mv(matrix, [q*r for q, r in zip(qx, mask)])]
                    shortcut = [mu*v/(x**quota*C) for v in first]
                    assert all(v <= upper for v, upper in zip(mv(B, f), candidate))
                    assert all(v <= upper for v, upper in zip(true_end, direct))
                    assert all(v <= upper for v, upper in zip(direct, shortcut))
                    counts['grouped_quota'] += size
                    counts['direct_endpoint'] += size
                    counts['shortcut_endpoint'] += size
                    bound = candidate if bound is None else list(map(min, bound, candidate))
                    end_bound = shortcut if end_bound is None else list(map(min, end_bound, shortcut))
                total = [a+b for a, b in zip(total, bound)]
                total_end = [a+b for a, b in zip(total_end, end_bound)]
            for i in range(size):
                assert group_inputs[0][i]+group_inputs[1][i] == kappa*phi[i]
                assert regular_endpoints[0][i]+regular_endpoints[1][i] == kappa*b[i]
                assert mv(B, phi)[i] <= total[i]/kappa <= total[i]
                assert b[i] <= total_end[i]/kappa <= total_end[i]
                counts['regularization_identity'] += 2
                counts['kappa_omission'] += 2
            rho = max(v/p for v, p in zip(total, phi))
            nu = max(v/p for v, p in zip(total_end, phi))
            assert mv(B, b)[-1] <= phi[-1]*rho*nu
            counts['composition'] += 1
            trace.clear()
            result = destination_quota_bound(width, 2, 2, 7, 3, xs_by_group,
                                             _dq_fraction_logs(phi), toy_prepare, toy_action)
            assert trace == [1, 1, 2, 2, 2, 2, 3, 3]
            counts['quota_order'] += 1
            for key, reference in (('log2_bulk_upper', total), ('log2_endpoint_upper', total_end)):
                assert np.max(np.abs(result[key]-_dq_fraction_logs(reference))) < DQ_LD('1e-10')
                counts['floating_rows'] += size
            scale = Fraction(9, 4)
            scaled = destination_quota_bound(width, 2, 2, 7, 3, xs_by_group,
                                             _dq_fraction_logs([scale*p for p in phi]), toy_prepare, toy_action)
            assert abs(scaled['log2_tail_upper']-result['log2_tail_upper']) < DQ_LD('1e-10')
            assert np.max(np.abs(scaled['log2_bulk_upper']-result['log2_bulk_upper']-
                                 np.log2(DQ_LD(9)/4))) < DQ_LD('1e-10')
            assert np.max(np.abs(scaled['log2_endpoint_upper']-result['log2_endpoint_upper'])) < DQ_LD('1e-10')
            counts['scale_rows'] += size
    print(dict(**counts, exact_width8_partitions=35, all_checks_tiny=True,
               full_code_certificate=False), flush=True)


def destination_quota_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--seconds', type=int, default=300)
    args = parser.parse_args()
    if not 30 <= args.seconds <= 300:
        raise ValueError('Use a bounded runtime between30 and300 seconds')
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
            destination_quota_checks()
        else:
            result = destination_quota_bound(deadline=started+args.seconds-10, progress=True)
            print(result['summary'], flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    destination_quota_main()
