"""Positive Horner action for an arbitrary strictly positive population vector.

Let T,S be the15-term homogeneous old-color quartics for a hole and a
selected quartet, including old-type multinomial coefficients and all
current label weights. For m=n/4 and positive f(4k,r), the exact numerator
of K_x f is
 sum_k C(m,k) [T^k coefficientwise*f(4k,Wdegree+Edegree)] S^(m-k).
The marked degree is the new intersection with OLD A, not just old W.

Use P0=1,H0=f(0,0). At stepk:
 Pk=T*Pprev; Bk=C(m,k)*Pk*f(4k,Wdegree+Edegree); Hk=S*Hprev+Bk.
Both products use15 positive shifts. P and H have separate common scales;
the positive B/H addition aligns those scales without discarding terms.
Finally divide once by the old population multinomial. There is no dense
33,153-square matrix, FFT, bin, or min-exponential envelope.

Input is a finite log2 vector in fixed order a=0,4,...,n; w=0,...,a.
Every input state is strictly positive. Zeros outside each homogeneous
coefficient triangle are structural and remain zero. Positive products,
coefficients, and scale alignments are range-checked. Output includes
absolute action logs and a normalized physical vector plus common scale.
This is noncertifying FLOAT arithmetic, not a regional/global bound.
No jobs or large arrays occur on import; parent owns execution.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _av_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_av_thread_var] = '1'

import numpy as np

AV_LD = np.longdouble
AV_TERMS = tuple((w, e) for w in range(5) for e in range(5-w))


def _av_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('_pt_') or name.startswith('general_flow'):
        import quartet_partial_temporal as module
    else:
        import quartet_adaptive_population as module
    return getattr(module, name)


def adaptive_vector_prepare(y=2, x=2, width=512, capacities=(1, 1, 2)):
    width = index(width)
    source, temporal, sink = _av_dependency('_pt_capacities')(capacities)
    y, x = AV_LD(y), AV_LD(x)
    if width < 4 or width > 512 or width % 4 or y <= 1 or x <= 0 or not np.isfinite(y+x):
        raise ValueError('Use quartet width4..512, finite y>1 and x>0')
    if width == 512 and np.finfo(AV_LD).maxexp < 16000:
        raise RuntimeError('Production vector action requires extended longdouble')
    log_y, log_x = np.log2(y), np.log2(x)
    hole_logs, selected_logs = [], []
    probability_logs = _av_dependency('_pt_probability_logs')
    positive_logsum = _av_dependency('_ap_logadd')
    for w, e in AV_TERMS:
        u = w+e
        hole, selected = [], []
        for r in range(w+1):
            for s in range(e+1):
                for t in range(5-u):
                    j = r+s+t
                    gain = sink*j-(source-temporal)*(r+s)-temporal*r
                    loss = max(0, 4*source+temporal*u-gain)
                    cost = (source-temporal)*(u-r-s)+temporal*(w-r)
                    log_g, log_p = probability_logs(log_y, loss, gain)
                    common = (np.log2(AV_LD(comb(w, r)*comb(e, s)*comb(4-u, t)))
                              +j*log_x+cost*log_y)
                    if np.isfinite(log_g):
                        hole.append(common+log_g)
                    if np.isfinite(log_p):
                        selected.append(common+log_p)
        if not hole or not selected:
            raise ValueError('This Horner engine requires strictly positive T and S old-type coefficients')
        multiplicity = np.log2(AV_LD(comb(4, u)*comb(u, w)))
        hole_logs.append(multiplicity+positive_logsum(hole))
        selected_logs.append(multiplicity+positive_logsum(selected))
    hole_logs, selected_logs = np.asarray(hole_logs, dtype=AV_LD), np.asarray(selected_logs, dtype=AV_LD)
    hole_scale, selected_scale = hole_logs.max(), selected_logs.max()
    hole, selected = np.exp2(hole_logs-hole_scale), np.exp2(selected_logs-selected_scale)
    if min(hole.min(), selected.min()) <= 0 or not np.all(np.isfinite(hole+selected)):
        raise ArithmeticError('Positive split-quartic normalization failed')
    a_rows = np.asarray([a for a in range(0, width+1, 4) for _ in range(a+1)], dtype=np.int16)
    w_rows = np.asarray([w for a in range(0, width+1, 4) for w in range(a+1)], dtype=np.int16)
    logs = np.zeros(width+1, dtype=AV_LD)
    logs[1:] = np.cumsum(np.log2(np.arange(1, width+1, dtype=AV_LD)), dtype=AV_LD)
    denominator = logs[width]-logs[w_rows]-logs[a_rows-w_rows]-logs[width-a_rows]
    axes = np.arange(width+1, dtype=np.int16)
    return dict(y=y, x=x, width=width, groups=width//4, capacities=(source, temporal, sink),
                T=hole, T_scale=hole_scale, S=selected, S_scale=selected_scale,
                T_log_coefficients=hole_logs, S_log_coefficients=selected_logs,
                a=a_rows, w=w_rows, old_multinomial_logs=denominator,
                degree_sums=axes[:, None]+axes[None, :],
                binomial_logs=np.asarray([np.log2(AV_LD(comb(width//4, k))) for k in range(width//4+1)]))


def _av_triangle_range(array, valid, stage):
    positive = array[valid]
    if np.any(positive <= 0) or not np.all(np.isfinite(positive)):
        raise ArithmeticError('Positive Horner coefficient underflow/overflow at '+stage)
    if np.any(array[~valid] != 0):
        raise ArithmeticError('Horner structural zero changed at '+stage)
    return positive.min(), positive.max()


def _av_shift_sum(current, coefficients, after_storage, product_storage):
    """One15-shift multiplication; caller carries scales and normalization."""
    size, next_size = current.shape[0], current.shape[0]+4
    after = after_storage[:next_size*next_size].reshape(next_size, next_size)
    scratch = product_storage[:size*size].reshape(size, size)
    after.fill(0)
    for (w, e), coefficient in zip(AV_TERMS, coefficients):
        np.multiply(current, coefficient, out=scratch)
        destination = after[w:w+size, e:e+size]
        np.add(destination, scratch, out=destination)
    return after


def adaptive_vector_action(prepared, log_vector, log_scale=0, progress=False, power_cache=None):
    """Return exact-algebra K_x f for arbitrary strictly positive f, in FLOAT."""
    started = time.monotonic()
    width, groups = prepared['width'], prepared['groups']
    log_vector = np.asarray(log_vector, dtype=AV_LD)
    if log_vector.shape != prepared['a'].shape or not np.all(np.isfinite(log_vector)) or not np.isfinite(log_scale):
        raise ValueError('Use one finite log2 value per positive physical input state')
    input_peak = log_vector.max()
    working = log_vector-input_peak
    external_scale = AV_LD(log_scale)+input_peak
    capacity = (width+1)**2
    if power_cache is None:
        p_storage, p_after = np.empty(capacity, dtype=AV_LD), np.empty(capacity, dtype=AV_LD)
        p_storage[0] = 1
        P = p_storage[:1].reshape(1, 1)
    elif (power_cache.get('prepared') is not prepared
          or len(power_cache.get('powers', ())) != groups+1
          or power_cache.get('T_scale') != prepared['T_scale']
          or not np.array_equal(power_cache.get('T'), prepared['T'])):
        raise ValueError('Use a validated T-power cache for this exact prepared object')
    h_storage, h_after = np.empty(capacity, dtype=AV_LD), np.empty(capacity, dtype=AV_LD)
    product_storage, b_storage = np.empty(capacity, dtype=AV_LD), np.empty(capacity, dtype=AV_LD)
    h_storage[0] = 1
    H = h_storage[:1].reshape(1, 1)
    p_scale, h_scale = AV_LD(0), working[0]
    p_min = h_min = AV_LD(1)
    T, S = prepared['T'], prepared['S']
    for k in range(1, groups+1):
        degree, size = 4*k, 4*k+1
        sums = prepared['degree_sums'][:size, :size]
        valid = sums <= degree
        if h_min*S.min() <= 0:
            raise ArithmeticError('A positive Horner shift product would underflow')
        if power_cache is None:
            if p_min*T.min() <= 0:
                raise ArithmeticError('A positive T-power shift product would underflow')
            next_P = _av_shift_sum(P, T, p_after, product_storage)
            p_min_raw, p_peak = _av_triangle_range(next_P, valid, 'T power')
            p_min = p_min_raw/p_peak
            if p_min <= 0:
                raise ArithmeticError('T-power normalization would underflow')
            np.divide(next_P, p_peak, out=next_P)
            p_scale += prepared['T_scale']+np.log2(p_peak)
            p_storage, p_after = p_after, p_storage
            P = p_storage[:size*size].reshape(size, size)
        else:
            P = power_cache['powers'][k]
            p_scale, p_min = power_cache['scales'][k], power_cache['minima'][k]
            if (P.shape != (size, size) or P.dtype != np.dtype(AV_LD)
                    or P.flags.writeable or p_min <= 0
                    or not np.isfinite(p_scale) or not np.isfinite(p_min)):
                raise ValueError('Invalid cached normalized T power')
        next_H = _av_shift_sum(H, S, h_after, product_storage)
        h_raw_min, h_peak = _av_triangle_range(next_H, valid, 'selected Horner term')
        h_term_scale = h_scale+prepared['S_scale']
        offset = 2*k*k-k  # First physical state with next population a'=4k.
        local_logs = working[offset:offset+size]
        local_scale = local_logs.max()
        weights = np.exp2(local_logs-local_scale)
        if np.any(weights <= 0) or not np.all(np.isfinite(weights)) or p_min*weights.min() <= 0:
            raise ArithmeticError('Positive next-intersection weighting would underflow')
        B = b_storage[:size*size].reshape(size, size)
        # Valid indices are exactly Wdegree+Edegree<=4k. Clipped indices are
        # outside that triangle, where P is identically zero; no coefficient
        # or probability is clipped.
        np.take(weights, sums, out=B, mode='clip')
        np.multiply(P, B, out=B)
        b_min, b_peak = _av_triangle_range(B, valid, 'new-hole Horner term')
        b_scale = p_scale+local_scale+prepared['binomial_logs'][k]
        common = max(h_term_scale+np.log2(h_peak), b_scale+np.log2(b_peak))
        h_factor, b_factor = np.exp2(h_term_scale-common), np.exp2(b_scale-common)
        if (not np.isfinite(h_factor) or not np.isfinite(b_factor)
                or h_factor <= 0 or b_factor <= 0
                or h_raw_min*h_factor <= 0 or b_min*b_factor <= 0):
            raise ArithmeticError('A positive Horner summand would disappear in scale alignment')
        np.multiply(next_H, h_factor, out=next_H)
        np.multiply(B, b_factor, out=B)
        np.add(next_H, B, out=next_H)
        h_min_raw, h_peak = _av_triangle_range(next_H, valid, 'aligned Horner sum')
        h_min = h_min_raw/h_peak
        if h_min <= 0:
            raise ArithmeticError('Horner normalization would underflow')
        np.divide(next_H, h_peak, out=next_H)
        h_scale = common+np.log2(h_peak)
        h_storage, h_after = h_after, h_storage
        H = h_storage[:size*size].reshape(size, size)
        if progress and (k % 32 == 0 or k == groups):
            print(dict(stage='adaptive_vector_horner', quartets=k, degree=degree,
                       seconds=round(time.monotonic()-started, 3)), flush=True)
    a, w = prepared['a'], prepared['w']
    unscaled_action = h_scale+np.log2(H[w, a-w])-prepared['old_multinomial_logs']
    if not np.all(np.isfinite(unscaled_action)):
        raise ArithmeticError('Invalid physical action logarithms')
    peak = unscaled_action.max()
    output_logs = unscaled_action-peak
    vector = np.exp2(output_logs)
    if np.any(vector <= 0) or not np.all(np.isfinite(vector)):
        raise ArithmeticError('Positive physical vector normalization underflow/overflow')
    return dict(a=a, w=w, vector=vector, scale=external_scale+peak,
                log2_scale=external_scale+peak,
                log2_vector=output_logs, log2_action=unscaled_action+external_scale,
                log2_ratio=unscaled_action-working, width=width,
                y=prepared['y'], x=prepared['x'], capacities=prepared['capacities'],
                arbitrary_strictly_positive_vector=True, next_intersection_degree='W+E',
                old_multinomial_divisions=1, independent_T_and_H_scales=True,
                cached_T_powers=power_cache is not None,
                structural_triangle_zeros_preserved=True, floating_diagnostic_only=True,
                regional_or_global_bound=False, seconds=time.monotonic()-started)


@lru_cache(maxsize=None)
def _av_physical_local(y, x, capacities, A, W):
    source, temporal, _ = capacities
    weights = _av_dependency('general_flow_weights')
    hole = selected = Fraction(0)
    for J in range(16):
        g, p, _, _, _ = weights(y, capacities, A.bit_count(), W.bit_count(), J.bit_count(),
                                (J & A).bit_count(), (J & W).bit_count())
        cost = (source-temporal)*((15 ^ J) & A).bit_count()+temporal*((15 ^ J) & W).bit_count()
        common = x**J.bit_count()*y**cost
        hole += common*g
        selected += common*p
    return hole, selected


def _av_physical_kernel(width, y, x, capacities):
    states = [(a, w) for a in range(0, width+1, 4) for w in range(a+1)]
    positions = {state: i for i, state in enumerate(states)}
    partitions = _av_dependency('_ap_partitions')(width)
    matrix = [[Fraction(0)]*len(states) for _ in states]
    for source, (a, w) in enumerate(states):
        A, W = (1 << a)-1, (1 << w)-1
        for partition in partitions:
            distribution = {(0, 0): Fraction(1)}
            for quartet in partition:
                local_A = sum(((A >> coordinate) & 1) << slot for slot, coordinate in enumerate(quartet))
                local_W = sum(((W >> coordinate) & 1) << slot for slot, coordinate in enumerate(quartet))
                hole, selected = _av_physical_local(y, x, capacities, local_A, local_W)
                after = {}
                for (aa, ww), weight in distribution.items():
                    for target, factor in (((aa, ww), selected), ((aa+4, ww+local_A.bit_count()), hole)):
                        after[target] = after.get(target, Fraction(0))+weight*factor
                distribution = after
            for target, value in distribution.items():
                matrix[source][positions[target]] += value/len(partitions)
    return states, matrix


def _av_fraction_logs(values):
    return np.asarray([np.log2(AV_LD(v.numerator)/v.denominator) for v in values], dtype=AV_LD)


def adaptive_vector_checks():
    counts = dict(split_coefficients=0, physical_actions=0, nonlinear_input_rows=0,
                  linear_combination_rows=0, unit_rows=0, exponential_rows=0, external_scale_rows=0)
    for width in (4, 8):
        for capacities, y, xs in (((1, 1, 2), Fraction(2), (Fraction(1), Fraction(2))),
                                  ((2, 1, 4), Fraction(3, 2), (Fraction(1), Fraction(9, 4)))):
            for x in xs:
                prepared = adaptive_vector_prepare(float(y), float(x), width, capacities)
                states, matrix = _av_physical_kernel(width, y, x, capacities)
                for j, (ww, ee) in enumerate(AV_TERMS):
                    u = ww+ee
                    hole, selected = _av_physical_local(y, x, capacities, (1 << u)-1, (1 << ww)-1)
                    multiplicity = comb(4, u)*comb(u, ww)
                    for key, exact in (('T', hole*multiplicity), ('S', selected*multiplicity)):
                        found = prepared[key][j]*np.exp2(prepared[key+'_scale'])
                        assert abs(found/(AV_LD(exact.numerator)/exact.denominator)-1) < AV_LD('1e-12')
                        counts['split_coefficients'] += 1
                first = [Fraction(i+2, i+1) for i in range(len(states))]
                second = [Fraction(3+a//4+w*w, 2+w) for a, w in states]
                combination = [a+2*b for a, b in zip(first, second)]
                exponential = [Fraction(3, 2)**a*Fraction(2, 3)**w for a, w in states]
                inputs = ([Fraction(1)]*len(states), first, second, combination, exponential)
                results = []
                for values in inputs:
                    exact = [sum(coefficient*value for coefficient, value in zip(row, values)) for row in matrix]
                    result = adaptive_vector_action(prepared, _av_fraction_logs(values))
                    expected = _av_fraction_logs(exact)
                    assert np.max(np.abs(result['log2_action']-expected)) < AV_LD('1e-10')
                    assert np.max(np.abs(result['log2_vector']+result['scale']-expected)) < AV_LD('1e-10')
                    counts['physical_actions'] += len(states)
                    results.append(result)
                counts['nonlinear_input_rows'] += 2*len(states)
                linear = np.exp2(results[1]['log2_action'])+2*np.exp2(results[2]['log2_action'])
                assert np.max(np.abs(np.exp2(results[3]['log2_action'])/linear-1)) < AV_LD('1e-10')
                counts['linear_combination_rows'] += len(states)
                endpoint = _av_dependency('general_flow_endpoint_logs')(prepared, float(y), float(x), capacities)
                assert np.max(np.abs(results[0]['log2_action']-endpoint)) < AV_LD('1e-10')
                counts['unit_rows'] += len(states)
                reference = _av_dependency('general_flow_population_action')(
                    float(y), float(x), 1.5, 2/3, width, capacities=capacities)
                assert np.max(np.abs(results[4]['log2_action']-reference['log2_action'])) < AV_LD('1e-10')
                counts['exponential_rows'] += len(states)
                shifted = adaptive_vector_action(prepared, _av_fraction_logs(first)-17, log_scale=23)
                assert np.max(np.abs(shifted['log2_action']-results[1]['log2_action']-6)) < AV_LD('1e-10')
                counts['external_scale_rows'] += len(states)
    print(dict(**counts, exact_width8_partitions=35, strictly_positive_nonexponential_inputs=True,
               no_dense_population_matrix=True, full_code_certificate=False), flush=True)


def adaptive_vector_control(y=2, x=2, unit=False, progress=False):
    started = time.monotonic()
    prepared = adaptive_vector_prepare(y, x)
    logs = np.zeros(len(prepared['a']), dtype=AV_LD) if unit else AV_LD('0.5')*(prepared['a']+prepared['w'])
    result = adaptive_vector_action(prepared, logs, progress=progress)
    action_seconds = result['seconds']
    if unit:
        expected = _av_dependency('general_flow_endpoint_logs')(prepared, y, x, (1, 1, 2))
        reference_seconds = 0.
    else:
        reference = _av_dependency('adaptive_population_action')(y, x, np.sqrt(AV_LD(2)), np.sqrt(AV_LD(2)), 512)
        expected, reference_seconds = reference['log2_action'], reference['seconds']
    error = np.max(np.abs(result['log2_action']-expected))
    if error > AV_LD('1e-8'):
        raise ArithmeticError('Production Horner/reference action identity failed')
    return dict(y=float(y), x=float(x), physical_rows=len(prepared['a']),
                input_mode='unit' if unit else 'alpha=beta=sqrt(2)',
                maximum_all_row_log2_error=float(error),
                minimum_log2_action=float(result['log2_action'].min()),
                maximum_log2_action=float(result['log2_action'].max()),
                horner_seconds=action_seconds, reference_seconds=reference_seconds,
                seconds=time.monotonic()-started, arbitrary_vector_engine=True,
                one_update_only=True, floating_diagnostic_only=True, full_code_certificate=False)


def adaptive_vector_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    mode.add_argument('--unit-control', action='store_true')
    parser.add_argument('--y', type=float, default=2.)
    parser.add_argument('--x', type=float, default=2.)
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
            adaptive_vector_checks()
        else:
            print(adaptive_vector_control(args.y, args.x, args.unit_control, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_vector_main()
