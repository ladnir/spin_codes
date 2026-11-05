"""One guarded positive polynomial action for the exact regional cover quota.

For K(x)=H diag(w_b(x)), compute [x^h]K(x)^R v for a strictly positive
potential v. Each step uses129 direct positive polynomial convolutions,
then one longdouble H @ coefficient-matrix product. No FFT or subtraction
is used. All states/degrees share one positive scalar normalization.

A common tilt x0 stores coefficient_j*x0^j. After step j retain only
 max(0,h-n*(R-j)) <= degree <= min(h,n*j).
Discarded degrees cannot contribute to the final quota. The tracked
offset indexes the physical degree; no additional offset gauge is removed.
At extraction subtract h*log2(x0) once. Intermediate structural zeros
in w_b products remain zero; every represented positive coefficient and
positive product is checked for underflow/overflow.

The default control computes ONE Bv for y2, where B=[x^4096]K(x)^16.
It obtains v from a bounded positive iteration of the existing Cauchy
regional comparison, then uses lambda=max(Bv/v). The already valid Cauchy
endpoint q gives mu=max(q/v), hence e_128 B^255 q<=v_128 lambda^255 mu.
The free-endpoint y512 is charged once, and the H union once. An exact
endpoint action is exposed separately, not silently added to the first job.

This is FLOAT longdouble arithmetic, not directed certification. Runtime
is unmeasured until the parent performs the first serial bounded control.
No jobs or large arrays are created on import. The standalone guard is
175 seconds,2GiB address space,and20GiB free /tmp space.
"""

import argparse
from math import comb, lgamma, log
from operator import index
import os
import shutil
import threading
import time

for _cq_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_cq_thread_var] = '1'

import numpy as np


CQ_LD = np.longdouble


def _cq_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('_sh_'):
        import quartet_selector_hit_moment as module
    else:
        import quartet_cover_holes as module
    return getattr(module, name)


def _cq_positive(array, stage):
    if np.any(array <= 0) or not np.all(np.isfinite(array)):
        raise ArithmeticError('Positive coefficient underflow/overflow at '+stage)


def _cq_positive_powers(poly, count):
    """Small positive powers with their scalar scales, before shared alignment."""
    normal = _cq_dependency('_sh_normalize')
    local, local_scale = normal(poly)
    powers = [(np.ones(1, dtype=CQ_LD), CQ_LD(0))]
    for _ in range(count):
        previous, scale = powers[-1]
        if previous.min()*local.min() <= 0:
            raise ArithmeticError('Positive quartet-power product would underflow')
        product = np.convolve(previous, local)
        _cq_positive(product, 'quartet powers')
        product, extra = normal(product)
        powers.append((product, scale+local_scale+extra))
    return powers


def cover_hole_quota_weights(prepared, tilt=2):
    """Build all w_b(x0*x) positively, with one shared coefficient scale."""
    width, tilt = prepared['width'], CQ_LD(tilt)
    if not np.isfinite(tilt) or tilt <= 0:
        raise ValueError('The common polynomial tilt must be finite and positive')
    groups = width//4
    tilted = np.asarray([tilt**j for j in range(5)], dtype=CQ_LD)
    G = prepared['G_coefficients']*tilted
    # P(x0*x)=x*Pbar(x); the x0 factor stays in Pbar's constant.
    Pbar = prepared['P_coefficients'][1:]*tilted[1:]
    _cq_positive(G, 'tilted G coefficients')
    _cq_positive(Pbar, 'tilted Pbar coefficients')
    Gpowers, Ppowers = _cq_positive_powers(G, groups), _cq_positive_powers(Pbar, groups)
    logs = np.full((groups+1, width+1), CQ_LD('-inf'))
    starts = []
    for b in range(groups+1):
        left, left_scale = Gpowers[b]
        right, right_scale = Ppowers[groups-b]
        if left.min()*right.min() <= 0:
            raise ArithmeticError('Positive population polynomial product would underflow')
        coefficients = np.convolve(left, right)
        _cq_positive(coefficients, 'population polynomial')
        start = groups-b
        starts.append(start)
        logs[b, start:] = (np.log2(coefficients)+left_scale+right_scale
                           +np.log2(CQ_LD(comb(groups, b))))
    supported = np.isfinite(logs)
    scale = logs[supported].max()
    weights = np.zeros_like(logs)
    weights[supported] = np.exp2(logs[supported]-scale)
    _cq_positive(weights[supported], 'common population-polynomial scale')
    return dict(coefficients=weights, scale=scale, starts=tuple(starts), tilt=tilt,
                width=width, y=prepared['y'])


def _cq_band(width, steps, quota, completed):
    return max(0, quota-width*(steps-completed)), min(quota, width*completed)


def cover_hole_quota_action(prepared, vector, region_updates=16, quota=4096,
                            tilt=2, endpoint=False, progress=False, weights=None):
    """Return scaled exact-quota Bv, or the separate endpoint coefficient action.

    endpoint=False: [x^h] K(x)^R v.
    endpoint=True:  [x^h] (1+x)^n K(x)^(R-1) v; take v=1 for q.
    All coefficients are evaluated in FLOAT, despite the exact algebraic quota.
    """
    width, steps, quota = prepared['width'], index(region_updates), index(quota)
    if steps < 1 or not 0 <= quota <= width*steps:
        raise ValueError('Invalid regional step count or quota')
    vector = np.asarray(vector, dtype=CQ_LD)
    states = width//4+1
    if vector.shape != (states,):
        raise ValueError('Potential has the wrong population dimension')
    _cq_positive(vector, 'initial potential')
    if width == 512 and np.finfo(CQ_LD).maxexp < 16000:
        raise RuntimeError('Production quota action requires genuinely extended longdouble')
    weights = cover_hole_quota_weights(prepared, tilt) if weights is None else weights
    if weights['width'] != width or weights['y'] != prepared['y'] or weights['tilt'] != CQ_LD(tilt):
        raise ValueError('Prepared polynomial weights do not match this action')
    overlap, overlap_scale = prepared['overlap'], CQ_LD(prepared['overlap_scale'])
    _cq_positive(overlap, 'overlap matrix')
    overlap_min = overlap.min()
    capacity = states*(min(quota, width*steps)+1)
    current_storage = np.empty(capacity, dtype=CQ_LD)
    weighted_storage = np.empty(capacity, dtype=CQ_LD)
    after_storage = np.empty(capacity, dtype=CQ_LD)
    current = current_storage[:states].reshape(states, 1)
    maximum = vector.max()
    current[:, 0] = vector/maximum
    scale, low, high = np.log2(maximum), 0, 0
    _cq_positive(current, 'normalized initial potential')
    elapsed_convolution = elapsed_mix = 0.0
    started = time.monotonic()
    labels = np.asarray([CQ_LD(comb(width, j))*weights['tilt']**j for j in range(width+1)], dtype=CQ_LD) if endpoint else None
    if endpoint:
        _cq_positive(labels, 'tilted first-layer labels')
        labels, label_scale = _cq_dependency('_sh_normalize')(labels)
    for completed in range(1, steps+1):
        next_low, next_high = _cq_band(width, steps, quota, completed)
        length = next_high-next_low+1
        weighted = weighted_storage[:states*length].reshape(states, length)
        weighted.fill(0)
        label_step = bool(endpoint and completed == steps)
        local_scale = label_scale if label_step else weights['scale']
        weighted_min = None
        phase = time.monotonic()
        for b in range(states):
            shift = 0 if label_step else weights['starts'][b]
            local = labels if label_step else weights['coefficients'][b, shift:]
            convolution_low = low+shift
            left = max(next_low, convolution_low)
            right = min(next_high, high+width)
            if left > right:
                continue  # Exact support, not numerical clipping.
            begin, end = left-convolution_low, right-convolution_low+1
            source = current[b]
            # Only indices participating in a retained coefficient are checked.
            source_lo, source_hi = max(0, begin-(len(local)-1)), min(len(source)-1, end-1)
            local_lo, local_hi = max(0, begin-(len(source)-1)), min(len(local)-1, end-1)
            if source[source_lo:source_hi+1].min()*local[local_lo:local_hi+1].min() <= 0:
                raise ArithmeticError('A retained positive convolution product would underflow')
            convolution = np.convolve(source, local)
            retained = convolution[begin:end]
            _cq_positive(retained, 'retained population convolution')
            local_minimum = retained.min()
            weighted_min = local_minimum if weighted_min is None else min(weighted_min, local_minimum)
            weighted[b, left-next_low:right-next_low+1] = retained
        elapsed_convolution += time.monotonic()-phase
        after = after_storage[:states*length].reshape(states, length)
        phase = time.monotonic()
        if label_step:
            after[:] = weighted
            action_scale = CQ_LD(0)
        else:
            if weighted_min is None or weighted_min*overlap_min <= 0:
                raise ArithmeticError('A positive overlap product would underflow')
            np.matmul(overlap, weighted, out=after)
            action_scale = overlap_scale
        elapsed_mix += time.monotonic()-phase
        _cq_positive(after, 'quota state/degree action')
        maximum = after.max()
        np.divide(after, maximum, out=after)
        _cq_positive(after, 'shared quota normalization')
        scale += local_scale+action_scale+np.log2(maximum)
        current_storage, after_storage = after_storage, current_storage
        current = current_storage[:states*length].reshape(states, length)
        low, high = next_low, next_high
        if progress:
            print(dict(stage='positive_quota_action', completed=completed, total_steps=steps,
                       first_layer_labels=label_step, degree_low=low, degree_high=high,
                       seconds=round(time.monotonic()-started, 3),
                       convolution_seconds=round(elapsed_convolution, 3),
                       mixing_seconds=round(elapsed_mix, 3)), flush=True)
    if low != quota or high != quota:
        raise ArithmeticError('Quota band did not end at the requested coefficient')
    return dict(vector=current[:, 0].copy(), scale=scale-quota*np.log2(weights['tilt']),
                width=width, region_updates=steps, quota=quota, endpoint=bool(endpoint),
                tilt=weights['tilt'], common_scale_all_states_and_degrees=True,
                positive_direct_convolution=True, exact_degree_band_pruning=True,
                convolution_seconds=elapsed_convolution, mixing_seconds=elapsed_mix,
                seconds=time.monotonic()-started, floating_diagnostic_only=True)


def cover_hole_quota_control(y=2, x_values=(1, 2, 4), tilt=2, witness_iterations=64, progress=False):
    """One exact-quota Bv action; retain the existing Cauchy endpoint bound."""
    iterations = index(witness_iterations)
    if not 1 <= iterations <= 128:
        raise ValueError('Use one to128 bounded potential iterations')
    prepared = _cq_dependency('cover_hole_prepare')(y, 512)
    regional = _cq_dependency('cover_hole_regional')(prepared, x_values, 16, 4096, progress)
    vector = np.ones(129, dtype=CQ_LD)
    for _ in range(iterations):
        vector, _ = _cq_dependency('_sh_normalize')(regional['ordinary'] @ vector)
    cauchy_action = regional['ordinary'] @ vector
    _cq_positive(cauchy_action, 'Cauchy potential action')
    cauchy_growth = regional['ordinary_scale']+np.log2(cauchy_action/vector).max()
    action = cover_hole_quota_action(prepared, vector, 16, 4096, tilt, False, progress)
    ratios = action['scale']+np.log2(action['vector'])-np.log2(vector)
    growth = ratios.max()
    endpoint_ratio = regional['endpoint_scale']+(np.log2(regional['endpoint'])-np.log2(vector)).max()
    raw = 255*growth+endpoint_ratio+np.log2(vector[-1])
    log_y = prepared['log_y']
    log_C = CQ_LD((lgamma(8193)-2*lgamma(4097))/log(2))
    conditional = raw+512*log_y-256*log_C
    tail = conditional+log_C-104860*log_y
    return dict(y=float(prepared['y']), x_values=tuple(map(float, regional['x_values'])), tilt=float(tilt),
                witness_iterations=iterations, potential_log2_min=float(np.log2(vector).min()),
                exact_quota_bulk_log2_collatz=float(growth),
                cauchy_bulk_log2_collatz=float(cauchy_growth),
                bulk_log2_collatz_gain=float(cauchy_growth-growth),
                maximum_ratio_population=4*int(np.argmax(ratios)),
                minimum_ratio_log2=float(ratios.min()),
                endpoint_log2_ratio=float(endpoint_ratio),
                endpoint_method='existing_Cauchy_upper_vector',
                log2_raw_coefficient_moment_upper=float(raw),
                log2_conditional_cover_moment_upper=float(conditional),
                log2_probability_upper=float(tail), total_updates=4096, bulk_regions=255,
                global_endpoint_charged_once=True, H_union_count=1,
                quota_action_seconds=action['seconds'],
                convolution_seconds=action['convolution_seconds'], mixing_seconds=action['mixing_seconds'],
                exact_quota_algebra_with_FLOAT_arithmetic=True,
                no_FFT=True, no_float64_action=True, floating_diagnostic_only=True,
                full_code_certificate=False)


def _cq_reference(prepared, rational_y, vector, steps, quota, endpoint):
    """Unpruned Fraction polynomial-matrix action, independent of FLOAT scales."""
    from fractions import Fraction
    width = prepared['width']
    kernel = _cq_dependency('_ch_kernel_polynomial')(width, rational_y)
    multiply = _cq_dependency('_ch_poly_product')
    current = [[Fraction(value)] for value in vector]
    kernel_steps = steps-1 if endpoint else steps
    for _ in range(kernel_steps):
        size = len(current[0])+width
        after = [[Fraction(0)]*size for _ in current]
        for a in range(len(current)):
            for b in range(len(current)):
                product_poly = multiply(kernel[a][b], current[b])
                for degree, coefficient in enumerate(product_poly):
                    after[a][degree] += coefficient
        current = after
    if endpoint:
        labels = [Fraction(comb(width, j)) for j in range(width+1)]
        current = [multiply(poly, labels) for poly in current]
    return [poly[quota] for poly in current]


def cover_hole_quota_checks():
    from fractions import Fraction
    coefficient_checks = action_checks = band_checks = collatz_checks = 0
    for width in (4, 8):
        for rational_y in (Fraction(2), Fraction(4)):
            prepared = _cq_dependency('cover_hole_prepare')(float(rational_y), width)
            kernel = _cq_dependency('_ch_kernel_polynomial')(width, rational_y)
            vector = [Fraction(i+2, i+1) for i in range(width//4+1)]
            float_vector = np.asarray([CQ_LD(v.numerator)/v.denominator for v in vector])
            for tilt in (1, 2, 3):
                weights = cover_hole_quota_weights(prepared, tilt)
                # H_0b=1, so the source-population-zero row is exactly w_b.
                for b, poly in enumerate(kernel[0]):
                    for degree, exact in enumerate(poly):
                        found = weights['coefficients'][b, degree]*np.exp2(weights['scale'])
                        reference = CQ_LD(exact.numerator)/exact.denominator*CQ_LD(tilt)**degree
                        if exact == 0:
                            assert found == 0
                        else:
                            assert abs(found/reference-1) < CQ_LD('1e-12')
                        coefficient_checks += 1
                for steps, quota in ((2, width), (3, width+1), (2, 0), (2, 2*width)):
                    for endpoint in (False, True):
                        reference = _cq_reference(prepared, rational_y, vector, steps, quota, endpoint)
                        action = cover_hole_quota_action(prepared, float_vector, steps, quota, tilt, endpoint, weights=weights)
                        found = action['vector']*np.exp2(action['scale'])
                        expected = np.asarray([CQ_LD(v.numerator)/v.denominator for v in reference])
                        assert np.max(np.abs(found/expected-1)) < CQ_LD('1e-11')
                        action_checks += len(vector)
                        for completed in range(steps+1):
                            low, high = _cq_band(width, steps, quota, completed)
                            for degree in range(width*completed+1):
                                reachable = degree <= quota <= degree+width*(steps-completed)
                                assert reachable == (low <= degree <= high)
                                band_checks += 1
                # An exact positive potential gives an entrywise two-region bound.
                exact_action = _cq_reference(prepared, rational_y, vector, 2, width, False)
                lam = max(value/potential for value, potential in zip(exact_action, vector))
                twice = _cq_reference(prepared, rational_y, exact_action, 2, width, False)
                assert all(value <= lam*lam*potential for value, potential in zip(twice, vector))
                collatz_checks += len(vector)
    print(dict(exact_weight_coefficient_comparisons=coefficient_checks,
               positive_quota_action_comparisons=action_checks,
               exact_degree_band_checks=band_checks,
               exact_collatz_composition_checks=collatz_checks,
               rational_markers=(2, 4), common_tilts=(1, 2, 3),
               exact_endpoint_callable_checked=True, full_code_certificate=False), flush=True)


def cover_hole_quota_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--y', type=float, default=2.)
    parser.add_argument('--x-values', nargs='+', type=float, default=(1., 2., 4.))
    parser.add_argument('--tilt', type=float, default=2.)
    parser.add_argument('--witness-iterations', type=int, default=64)
    args = parser.parse_args()
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
            cover_hole_quota_checks()
        else:
            result = cover_hole_quota_control(args.y, args.x_values, args.tilt, args.witness_iterations, True)
            print(dict(seconds=round(time.monotonic()-started, 3), **result), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    cover_hole_quota_main()
