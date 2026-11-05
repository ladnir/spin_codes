"""Optional normalized T-power cache for the positive Horner action.

The default adaptive_vector_action remains uncached. This helper prepares
the same T^k powers once and passes them through its power_cache argument.
Each read-only power retains its common log2 scale and positive minimum.
No input-vector weight enters the cache. H, the binomial factor, and all
scale alignment are still evaluated by the original action implementation.

At width512, storing all minimal square arrays costs about174MiB. The
workspace guard includes already-retained caches and a conservative ten
full-array allowance for preparation/action temporaries. Callers retaining
several x caches must pass their total bytes as retained_cache_bytes.
This is a FLOAT acceleration, not directed or certified arithmetic.
Import performs no numerical job; parent owns all execution.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _avc_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_avc_thread_var] = '1'

import numpy as np

AVC_LD = np.longdouble


def _avc_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_adaptive_vector as module
    return getattr(module, name)


def adaptive_vector_cache_prepare(prepared, retained_cache_bytes=0,
                                  workspace_limit_bytes=1 << 30, progress=False):
    """Cache validated T powers; the prepared object must be reused unchanged."""
    started = time.monotonic()
    width, groups = prepared['width'], prepared['groups']
    retained_cache_bytes, workspace_limit_bytes = index(retained_cache_bytes), index(workspace_limit_bytes)
    if (retained_cache_bytes < 0 or workspace_limit_bytes <= 0
            or workspace_limit_bytes > 1536*(1 << 20)):
        raise ValueError('Use nonnegative retained bytes and a workspace limit at most1536MiB')
    itemsize = np.dtype(AVC_LD).itemsize
    power_bytes = itemsize*sum((4*k+1)**2 for k in range(groups+1))
    metadata_bytes = itemsize*(2*(groups+1)+len(prepared['T']))
    cache_bytes = power_bytes+metadata_bytes
    temporary_bytes = 10*itemsize*(width+1)**2
    estimated_workspace_bytes = retained_cache_bytes+cache_bytes+temporary_bytes
    if cache_bytes > 256*(1 << 20) or estimated_workspace_bytes > workspace_limit_bytes:
        raise MemoryError('T-power cache and reserved action workspace exceed the explicit memory budget')
    shift_sum, triangle_range = _avc_dependency('_av_shift_sum'), _avc_dependency('_av_triangle_range')
    after_storage = np.empty((width+1)**2, dtype=AVC_LD)
    product_storage = np.empty((width+1)**2, dtype=AVC_LD)
    first = np.ones((1, 1), dtype=AVC_LD)
    first.flags.writeable = False
    powers = [first]
    scales, minima = np.zeros(groups+1, dtype=AVC_LD), np.ones(groups+1, dtype=AVC_LD)
    T = prepared['T']
    for k in range(1, groups+1):
        if minima[k-1]*T.min() <= 0:
            raise ArithmeticError('A positive cached T-power shift product would underflow')
        after = shift_sum(powers[-1], T, after_storage, product_storage)
        size = 4*k+1
        valid = prepared['degree_sums'][:size, :size] <= 4*k
        low, peak = triangle_range(after, valid, 'cached T power')
        minimum = low/peak
        if minimum <= 0 or not np.isfinite(minimum):
            raise ArithmeticError('Cached T-power normalization would underflow')
        np.divide(after, peak, out=after)
        # Copy only the current square: retaining a view would pin the entire
        # scratch buffer, and reusing that buffer would corrupt earlier powers.
        power = after.copy()
        power.flags.writeable = False
        powers.append(power)
        scales[k] = scales[k-1]+prepared['T_scale']+np.log2(peak)
        minima[k] = minimum
        if not np.isfinite(scales[k]):
            raise ArithmeticError('Invalid cached T-power scale')
        if progress and (k % 32 == 0 or k == groups):
            print(dict(stage='adaptive_vector_cache', quartets=k,
                       retained_MiB=round(sum(p.nbytes for p in powers)/(1 << 20), 3),
                       seconds=round(time.monotonic()-started, 3)), flush=True)
    T_snapshot = T.copy()
    for array in (T_snapshot, scales, minima):
        array.flags.writeable = False
    return dict(prepared=prepared, powers=tuple(powers), scales=scales, minima=minima,
                T=T_snapshot, T_scale=prepared['T_scale'], cache_bytes=cache_bytes,
                estimated_workspace_bytes=estimated_workspace_bytes,
                workspace_limit_bytes=workspace_limit_bytes,
                input_independent=True, structural_triangle_zeros_preserved=True,
                floating_diagnostic_only=True, seconds=time.monotonic()-started)


def _avc_fraction_powers(prepared, y, x, capacities):
    """Independent small exact coefficient convolution, including old Mult4."""
    physical = _avc_dependency('_av_physical_local')
    coefficients = {}
    for w in range(5):
        for e in range(5-w):
            u = w+e
            hole, _ = physical(y, x, capacities, (1 << u)-1, (1 << w)-1)
            coefficients[w, e] = comb(4, u)*comb(u, w)*hole
    powers = [{(0, 0): Fraction(1)}]
    for _ in range(prepared['groups']):
        after = {}
        for (w, e), value in powers[-1].items():
            for (ww, ee), coefficient in coefficients.items():
                target = w+ww, e+ee
                after[target] = after.get(target, Fraction(0))+value*coefficient
        powers.append(after)
    return powers


def adaptive_vector_cache_checks():
    prepare, action = _avc_dependency('adaptive_vector_prepare'), _avc_dependency('adaptive_vector_action')
    fraction_logs, physical_kernel = _avc_dependency('_av_fraction_logs'), _avc_dependency('_av_physical_kernel')
    counts = dict(exact_power_coefficients=0, cached_uncached_rows=0, physical_rows=0,
                  scale_rows=0, immutable_arrays=0, mismatch_guards=0, memory_guards=0)
    for width in (4, 8):
        for capacities, y, x in (((1, 1, 2), Fraction(2), Fraction(2)),
                                  ((2, 1, 4), Fraction(3, 2), Fraction(9, 4))):
            prepared = prepare(float(y), float(x), width, capacities)
            cache = adaptive_vector_cache_prepare(prepared)
            exact_powers = _avc_fraction_powers(prepared, y, x, capacities)
            for k, exact in enumerate(exact_powers):
                actual = cache['powers'][k]*np.exp2(cache['scales'][k])
                for (w, e), value in exact.items():
                    expected = AVC_LD(value.numerator)/value.denominator
                    assert abs(actual[w, e]/expected-1) < AVC_LD('1e-11')
                    counts['exact_power_coefficients'] += 1
                assert not cache['powers'][k].flags.writeable
                counts['immutable_arrays'] += 1
            states, matrix = physical_kernel(width, y, x, capacities)
            first = [Fraction(3+a//4+w*w, 2+w) for a, w in states]
            second = [Fraction(i+2, i+1) for i in range(len(states))]
            inputs = ([Fraction(1)]*len(states), first, second,
                      [a+2*b for a, b in zip(first, second)])
            for values in inputs:
                logs = fraction_logs(values)
                uncached = action(prepared, logs)
                cached = action(prepared, logs, power_cache=cache)
                assert np.max(np.abs(cached['log2_action']-uncached['log2_action'])) < AVC_LD('1e-11')
                counts['cached_uncached_rows'] += len(states)
                exact = [sum(coefficient*value for coefficient, value in zip(row, values)) for row in matrix]
                assert np.max(np.abs(cached['log2_action']-fraction_logs(exact))) < AVC_LD('1e-10')
                counts['physical_rows'] += len(states)
            shifted = action(prepared, fraction_logs(first)-17, log_scale=23, power_cache=cache)
            unshifted = action(prepared, fraction_logs(first), power_cache=cache)
            assert np.max(np.abs(shifted['log2_action']-unshifted['log2_action']-6)) < AVC_LD('1e-11')
            counts['scale_rows'] += len(states)
            other = prepare(float(y), float(x), width, capacities)
            try:
                action(other, np.zeros(len(states), dtype=AVC_LD), power_cache=cache)
            except ValueError:
                counts['mismatch_guards'] += 1
            else:
                raise AssertionError('A different prepared object accepted the cache')
            try:
                adaptive_vector_cache_prepare(prepared, workspace_limit_bytes=1)
            except MemoryError:
                counts['memory_guards'] += 1
            else:
                raise AssertionError('Explicit cache memory guard failed')
    print(dict(**counts, uncached_default_preserved=True, no_large_check_arrays=True,
               full_code_certificate=False), flush=True)


def adaptive_vector_cache_control(y=2, x=2, progress=False):
    started = time.monotonic()
    prepare, action = _avc_dependency('adaptive_vector_prepare'), _avc_dependency('adaptive_vector_action')
    prepared = prepare(y, x)
    cache = adaptive_vector_cache_prepare(prepared, progress=progress)
    logs = AVC_LD('0.5')*(prepared['a']+prepared['w'])
    cached = action(prepared, logs, power_cache=cache, progress=progress)
    uncached = action(prepared, logs, progress=progress)
    error = np.max(np.abs(cached['log2_action']-uncached['log2_action']))
    if error > AVC_LD('1e-8'):
        raise ArithmeticError('Cached/uncached production action mismatch')
    return dict(y=float(y), x=float(x), physical_rows=len(logs),
                cache_MiB=cache['cache_bytes']/(1 << 20),
                estimated_workspace_MiB=cache['estimated_workspace_bytes']/(1 << 20),
                cache_prepare_seconds=cache['seconds'], cached_seconds=cached['seconds'],
                uncached_seconds=uncached['seconds'], maximum_all_row_log2_error=float(error),
                seconds=time.monotonic()-started, one_update_only=True,
                floating_diagnostic_only=True, full_code_certificate=False)


def adaptive_vector_cache_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
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
            adaptive_vector_cache_checks()
        else:
            print(adaptive_vector_cache_control(args.y, args.x, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_vector_cache_main()
