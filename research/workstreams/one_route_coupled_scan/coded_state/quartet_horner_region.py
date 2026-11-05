"""Bounded regional caller for the exact-algebra physical Horner action.

At each fixed x=1,2,4, apply K_x sixteen times to the SAME positive
physical potential. Only afterward take the common-population minimum
of x^-4096 K_x^16 phi/C0. Save step15 for the endpoint majorant
 x^-4096 mu_x K_x^15 phi/C0, mu_x=max(q_x/phi).
Then C0*2^-60518*phi(full)*rho^255*nu is a geometry-tail comparison.
There is no bin/envelope projection or extra terminal2^512 factor.

Initial phi(a,w)=2^((a+w)/2); its endpoint multiplier is exactly
(1+x)^512. Any later physical potential recomputes mu_x over all rows.
One witness is the default. A second must be explicitly requested and
fit the measured remaining runtime. It uses fixed geometric half damping
phi_next=sqrt(phi*F(phi)). All caches remain in this process.

The physical minimum F(phi)/phi is a lower floor for scalar supersolutions
of THIS nonlinear quota comparison F, not a lower bound on physical B.
If F(phi)>=lambda*phi and F(psi)<=rho*psi, monotonicity, homogeneity,
iteration and positive finite-state comparability imply lambda<=rho.
Half damping cannot increase max F(phi)/phi: with R=F(phi)/phi<=rho,
F(sqrt(phi*F(phi)))<=sqrt(rho)*F(phi), so the new ratio is at most
sqrt(rho*R)<=rho. Endpoint/full-tail improvement is not guaranteed.

This is noncertifying FLOAT arithmetic despite the exact polynomial
identity. Positive terms/scale guards belong to the action core; this
caller additionally checks every physical logarithm and row ordering.
No jobs, arrays, files, or core imports occur merely from importing this
module. Parent alone runs checks and production controls, serially.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _hr_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_hr_thread_var] = '1'

import numpy as np

HR_LD = np.longdouble
HR_MEMORY_LIMIT = 2 << 30


def _hr_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('_ab_'):
        import quartet_adaptive_bins as module
    else:
        import quartet_adaptive_vector as module
    return getattr(module, name)


def _hr_core_prepare(y, x, width, capacities):
    """Isolated stable core adapter; no global state or implicit scale."""
    return _hr_dependency('adaptive_vector_prepare')(y, x, width, capacities)


def _hr_core_action(cache, log_values):
    """Return dict(log2_action=absolute logs of K_x*2**log_values)."""
    return _hr_dependency('adaptive_vector_action')(cache, log_values, log_scale=0, progress=False)


def _hr_finite_logs(values, shape, where):
    result = np.asarray(values, dtype=HR_LD)
    if result.shape != shape or not np.all(np.isfinite(result)):
        raise ArithmeticError('Invalid strictly positive physical log-vector at '+where)
    return result


def _hr_deadline(deadline):
    if deadline is not None and time.monotonic() >= deadline:
        raise TimeoutError('Regional witness incomplete at the bounded deadline')


def horner_region_prepare(width=512, region_updates=16, x_values=(1, 2, 4),
                          prepare_callback=None, action_callback=None):
    """Prepare common physical rows and one reusable core cache per x.

    The callback contract is prepare(y,x,width,capacities)->cache with
    a,w arrays, and action(cache,log_values)->dict(log2_action=...).
    Input/output action logs contain no implicit external common scale.
    Small widths/regions are supported only for the independent checks.
    """
    width, steps = index(width), index(region_updates)
    xs = tuple(HR_LD(x) for x in x_values)
    if width < 4 or width > 512 or width % 4 or steps < 2:
        raise ValueError('Use quartet width4..512 and at least two regional updates')
    if not xs or any(x <= 0 or not np.isfinite(x) for x in xs):
        raise ValueError('Quota evaluations must be finite and positive')
    prepare = _hr_core_prepare if prepare_callback is None else prepare_callback
    action = _hr_core_action if action_callback is None else action_callback
    if not callable(prepare) or not callable(action):
        raise ValueError('Core adapters must be callable')
    a = np.asarray([a for a in range(0, width+1, 4) for _ in range(a+1)], dtype=np.int16)
    w = np.asarray([w for a in range(0, width+1, 4) for w in range(a+1)], dtype=np.int16)
    prepared = []
    for x in xs:
        cache = prepare(2, x, width, (1, 1, 2))
        if not np.array_equal(cache['a'], a) or not np.array_equal(cache['w'], w):
            raise ValueError('Core preparation changed physical population order')
        prepared.append(cache)
    quota = width*steps//2
    # Production C(8192,4096) fits extended longdouble. The core rejects
    # insufficient exponent range for the width512 production action.
    log_C = np.log2(HR_LD(comb(width*steps, quota)))
    if not np.isfinite(log_C):
        raise ArithmeticError('Quota binomial logarithm overflow')
    return dict(width=width, steps=steps, x_values=xs, a=a, w=w,
                prepared=prepared, action=action, quota=quota, log_C=log_C,
                full=len(a)-1, y=HR_LD(2))


def _hr_chain(cache, core, log_phi, deadline):
    """Normalize each action and carry its one external scale explicitly."""
    scale = log_phi.max()
    working = log_phi-scale
    step15 = None
    for step in range(1, cache['steps']+1):
        _hr_deadline(deadline)
        result = cache['action'](core, working)
        if 'a' in result and not np.array_equal(result['a'], cache['a']):
            raise ValueError('Core action changed old a ordering')
        if 'w' in result and not np.array_equal(result['w'], cache['w']):
            raise ValueError('Core action changed old w ordering')
        action_logs = _hr_finite_logs(result['log2_action'], log_phi.shape, 'action')
        peak = action_logs.max()
        working = action_logs-peak
        scale += peak
        if not np.isfinite(scale):
            raise ArithmeticError('Physical chain common scale overflow')
        if step == cache['steps']-1:
            step15 = working.copy()+scale
    _hr_deadline(deadline)
    return step15, working+scale


def horner_region_witness(cache, log_phi, regions=256, threshold=60518,
                          endpoint_upper=None, deadline=None, progress=False):
    """One complete witness; arrays remain internal unless caller requests them.

    Optional endpoint_upper is any already valid absolute log upper on
    the SAME exact first-region endpoint. Its physical minimum remains
    valid and may improve a later witness. No old potential is reused
    inside an x chain.
    """
    started = time.monotonic()
    regions, threshold = index(regions), index(threshold)
    if regions < 1 or threshold < 0:
        raise ValueError('Use positive region count and a nonnegative deficit threshold')
    phi = _hr_finite_logs(log_phi, cache['a'].shape, 'potential')
    F = np.full(phi.shape, HR_LD('inf'))
    endpoint = np.full(phi.shape, HR_LD('inf'))
    choices = np.zeros(phi.shape, dtype=np.int16)
    summaries = []
    for xi, (x, core) in enumerate(zip(cache['x_values'], cache['prepared'])):
        x_started = time.monotonic()
        saved, final = _hr_chain(cache, core, phi, deadline)
        correction = cache['quota']*np.log2(x)+cache['log_C']
        candidate = final-correction
        use = candidate < F
        choices[use] = xi
        np.minimum(F, candidate, out=F)
        q_logs = cache['w']*np.log2(2+x)+(cache['width']-cache['w'])*np.log2(1+x)
        mu = (q_logs-phi).max()
        candidate_end = saved+mu-correction
        np.minimum(endpoint, candidate_end, out=endpoint)
        worst = int(np.argmax(candidate-phi))
        summary = dict(x=float(x), seconds=round(time.monotonic()-x_started, 3),
                       log2_bulk_ratio=float((candidate-phi)[worst]),
                       bulk_row=(int(cache['a'][worst]), int(cache['w'][worst])),
                       log2_endpoint_multiplier=float(mu),
                       log2_endpoint_ratio=float((candidate_end-phi).max()))
        summaries.append(summary)
        if progress:
            print(dict(stage='horner_region_x', **summary), flush=True)
    if endpoint_upper is not None:
        previous = _hr_finite_logs(endpoint_upper, phi.shape, 'previous endpoint upper')
        np.minimum(endpoint, previous, out=endpoint)
    F = _hr_finite_logs(F, phi.shape, 'complete regional upper')
    endpoint = _hr_finite_logs(endpoint, phi.shape, 'complete endpoint upper')
    row = int(np.argmax(F-phi))
    minimum_row = int(np.argmin(F-phi))
    end_row = int(np.argmax(endpoint-phi))
    rho, nu = (F-phi)[row], (endpoint-phi)[end_row]
    tail = cache['log_C']-threshold+phi[cache['full']]+(regions-1)*rho+nu
    summary = dict(seconds=round(time.monotonic()-started, 3),
                   rows=len(phi), updates_per_region=cache['steps'], regions=regions,
                   threshold=threshold, y=2, quota=cache['quota'],
                   log2_rho=float(rho), log2_nu=float(nu),
                   log2_comparison_lower_floor=float((F-phi)[minimum_row]),
                   comparison_lower_row=(int(cache['a'][minimum_row]), int(cache['w'][minimum_row])),
                   lower_floor_is_for_nonlinear_comparison_only=True,
                   log2_phi_full=float(phi[cache['full']]), log2_tail_upper=float(tail),
                   rho_row=(int(cache['a'][row]), int(cache['w'][row])),
                   endpoint_row=(int(cache['a'][end_row]), int(cache['w'][end_row])),
                   quota_choice_counts=[int(np.count_nonzero(choices == i)) for i in range(len(summaries))],
                   exact_physical_action_algebra=True, floating_diagnostic_only=True,
                   full_code_certificate=False)
    return dict(summary=summary, per_x=summaries, log2_F=F,
                log2_endpoint_upper=endpoint, log2_tail_upper=tail)


def horner_region_control(iterations=1, deadline=None, progress=True):
    iterations = index(iterations)
    if iterations not in (1, 2):
        raise ValueError('This bounded control permits exactly one or two requested witnesses')
    started = time.monotonic()
    cache = horner_region_prepare()
    phi = (cache['a'].astype(HR_LD)+cache['w'])/2
    endpoint, best = None, None
    completed = []
    for iteration in range(iterations):
        _hr_deadline(deadline)
        if iteration and deadline is not None:
            # The first completed witness is the only runtime estimate.
            # Retain it rather than starting a chain unlikely to finish.
            estimate = completed[-1]['seconds']*HR_LD('1.15')+5
            if time.monotonic()+estimate >= deadline:
                if progress:
                    print(dict(stage='horner_region_stop', reason='second_witness_exceeds_remaining_budget',
                               completed_witnesses=len(completed)), flush=True)
                break
        result = horner_region_witness(cache, phi, endpoint_upper=endpoint,
                                      deadline=deadline, progress=progress)
        summary = dict(iteration=iteration+1, **result['summary'])
        completed.append(summary)
        if progress:
            print(dict(stage='horner_region_witness', **summary), flush=True)
        if best is None or result['log2_tail_upper'] < best['log2_tail_upper']:
            best = dict(summary=summary, log2_tail_upper=result['log2_tail_upper'])
        endpoint = result['log2_endpoint_upper']
        # Fixed geometric half damping contracts the Collatz interval
        # for a positive monotone homogeneous F. This does not imply
        # that the endpoint or complete tail bound improves.
        phi = (phi+result['log2_F'])/2
        phi -= phi.max()
    return dict(stage='horner_region_complete', seconds=round(time.monotonic()-started, 3),
                requested_witnesses=iterations, completed_witnesses=len(completed),
                best=best['summary'], floating_diagnostic_only=True, full_code_certificate=False)


def _hr_fraction_logs(values):
    return np.asarray([np.log2(HR_LD(v.numerator)/HR_LD(v.denominator)) for v in values], dtype=HR_LD)


def horner_region_checks():
    """Independent width4/8 coefficient quotas and physical endpoints."""
    counts = dict(two_step_coefficients=0, quota_rows=0, endpoint_rows=0,
                  initial_endpoint=0, scale_rows=0, floating_rows=0,
                  full_comparisons=0, fixed_x_order=0, damping_rows=0,
                  lower_floor_rows=0)
    mv = _hr_dependency('_ab_mv')
    for width in (4, 8):
        states, polys = _hr_dependency('_ab_exact_kernel')(width, Fraction(2))
        size, quota = len(states), width
        C = Fraction(comb(2*width, width))
        B = [[sum(polys[i][k][d]*polys[k][j][quota-d]
                  for k in range(size) for d in range(width+1))/C
              for j in range(size)] for i in range(size)]
        counts['two_step_coefficients'] += size*size
        q_polys = [[sum(Fraction(comb(w, r)*comb(width-w, j-r))*2**(w-r)
                        for r in range(max(0, j-(width-w)), min(w, j)+1))
                    for j in range(width+1)] for a, w in states]
        b = [sum(polys[i][k][d]*q_polys[k][quota-d]
                 for k in range(size) for d in range(width+1))/C for i in range(size)]
        matrices = {x: [[sum(v*x**j for j, v in enumerate(poly)) for poly in row]
                        for row in polys] for x in (Fraction(1), Fraction(2), Fraction(4))}
        trace = []

        def toy_prepare(y, x, n, capacities):
            assert y == 2 and n == width and capacities == (1, 1, 2)
            matrix = matrices[Fraction(int(x))]
            return dict(a=np.asarray([a for a, w in states]),
                        w=np.asarray([w for a, w in states]), x=int(x),
                        matrix=np.asarray([[HR_LD(v.numerator)/v.denominator for v in row]
                                           for row in matrix], dtype=HR_LD))

        def toy_action(core, log_values):
            trace.append(core['x'])
            peak = log_values.max()
            values = core['matrix'] @ np.exp2(log_values-peak)
            assert np.all(values > 0)
            return dict(log2_action=np.log2(values)+peak, a=core['a'], w=core['w'])

        cache = horner_region_prepare(width, 2, (1, 2, 4), toy_prepare, toy_action)
        for phi in ([Fraction(2)**a*Fraction(3, 2)**w for a, w in states],
                    [Fraction(i*i+2, i+3) for i in range(size)]):
            F, endpoint = None, None
            for x, matrix in matrices.items():
                first, second = mv(matrix, phi), mv(matrix, mv(matrix, phi))
                qx = [(2+x)**w*(1+x)**(width-w) for a, w in states]
                mu = max(v/p for v, p in zip(qx, phi))
                local = [v/(x**quota*C) for v in second]
                end = [mu*v/(x**quota*C) for v in first]
                F = local if F is None else list(map(min, F, local))
                endpoint = end if endpoint is None else list(map(min, endpoint, end))
                for (a, w), q in zip(states, qx):
                    assert q*q <= (1+x)**(2*width)*2**(a+w)
                    counts['initial_endpoint'] += 1
                assert qx[0] == (1+x)**width
            assert all(v <= bound for v, bound in zip(mv(B, phi), F))
            assert all(v <= bound for v, bound in zip(b, endpoint))
            counts['quota_rows'] += size
            counts['endpoint_rows'] += size
            rho = max(v/p for v, p in zip(F, phi))
            nu = max(v/p for v, p in zip(endpoint, phi))
            exact = mv(B, b)[-1]
            assert exact <= phi[-1]*rho*nu
            counts['full_comparisons'] += 1
            trace.clear()
            result = horner_region_witness(cache, _hr_fraction_logs(phi), regions=2, threshold=7)
            assert trace == [1, 1, 2, 2, 4, 4]
            counts['fixed_x_order'] += 1
            for key, reference in (('log2_F', F), ('log2_endpoint_upper', endpoint)):
                assert np.max(np.abs(result[key]-_hr_fraction_logs(reference))) < HR_LD('1e-10')
                counts['floating_rows'] += size
            ratio_logs = result['log2_F']-_hr_fraction_logs(phi)
            assert abs(result['summary']['log2_comparison_lower_floor']-
                       float(ratio_logs.min())) < 1e-10
            counts['lower_floor_rows'] += size
            damped_phi = (_hr_fraction_logs(phi)+result['log2_F'])/2
            damped = horner_region_witness(cache, damped_phi, regions=2, threshold=7)
            damped_ratios = damped['log2_F']-damped_phi
            assert damped_ratios.max() <= ratio_logs.max()+HR_LD('1e-10')
            assert damped_ratios.min() >= ratio_logs.min()-HR_LD('1e-10')
            counts['damping_rows'] += size
            scale = Fraction(9, 4)
            scaled = horner_region_witness(cache, _hr_fraction_logs([scale*p for p in phi]),
                                           regions=2, threshold=7,
                                           endpoint_upper=result['log2_endpoint_upper'])
            assert abs(scaled['log2_tail_upper']-result['log2_tail_upper']) < HR_LD('1e-10')
            assert np.max(np.abs(scaled['log2_F']-result['log2_F']-
                                 np.log2(HR_LD(9)/4))) < HR_LD('1e-10')
            assert np.max(np.abs(scaled['log2_endpoint_upper']-result['log2_endpoint_upper'])) < HR_LD('1e-10')
            counts['scale_rows'] += size
    try:
        _hr_deadline(time.monotonic()-1)
    except TimeoutError:
        pass
    else:
        raise AssertionError('Expired deadline was not rejected')
    print(dict(**counts, deadline_guard=1, exact_width8_partitions=35,
               checks_are_tiny=True, full_code_certificate=False), flush=True)


def horner_region_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--iterations', type=int, choices=(1, 2), default=1)
    parser.add_argument('--seconds', type=int, default=300)
    args = parser.parse_args()
    if not 30 <= args.seconds <= 300:
        raise ValueError('Use a bounded runtime between30 and300 seconds')
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (HR_MEMORY_LIMIT, HR_MEMORY_LIMIT))
    started = time.monotonic()
    timer = threading.Timer(args.seconds-5, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            horner_region_checks()
        else:
            print(horner_region_control(args.iterations, started+args.seconds-10, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    horner_region_main()
