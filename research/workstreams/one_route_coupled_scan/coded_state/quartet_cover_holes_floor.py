"""Constant-hole-population lower controls for the exact-quota cover comparison.

Fix one hole-quartet count a for every real layer T,...,2. Write
 c_a=[x^h]w_a(x)^R,
 e_a=[x^h](1+x)^n w_a(x)^(R-1),
where R is the regional length and T=R*regions. The contribution of this
single population path to e_full B^(regions-1)q is exactly
 y^(4a) H_aa^(T-2) c_a^(regions-1) e_a.
The first overlap H_full,a=y^(4a) is separate. The input-1 label factor
has no cover coin. The remaining free-endpoint factor y^n is external.

All coefficient operations below are positive direct scalar convolutions.
A common polynomial tilt is removed once, at the requested coefficient;
safe degree bands discard only terms that cannot reach that coefficient.
Production performs no matrix powers or Perron iterations.

This is a FLOAT diagnostic of a mathematically exact path contribution.
A positive comparison-tail floor would obstruct this pair-cover bound at
the chosen y, even with exact quotas and an optimal potential. It would
NOT lower-bound the actual failure probability, obstruct the exact-own
cover bound, or rule out other markers. Floating arithmetic is not a
directed certificate. A separate integer/rational histogram witness,
included in --checks, certifies a y2 pair-comparison floor >=2^42132.
No jobs or large arrays are created on import.
"""

import argparse
from math import comb, lgamma, log
from operator import index
import os
import shutil
import threading
import time

for _cf_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_cf_thread_var] = '1'

import numpy as np


CF_LD = np.longdouble


def _cf_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('cover_hole_quota'):
        import quartet_cover_holes_quota as module
    else:
        import quartet_cover_holes as module
    return getattr(module, name)


def _cf_positive(array, stage):
    if np.any(array <= 0) or not np.all(np.isfinite(array)):
        raise ArithmeticError('Positive coefficient underflow/overflow at '+stage)


def _cf_coefficient_log2(factors, quota, tilt):
    """Coefficient of a product of (positive contiguous array, scale, offset).

    Each array represents coefficients already tilted by the SAME tilt.
    Arrays contain no structural zeros; their offsets encode lower support.
    A return value -inf means the requested coefficient is structurally zero.
    """
    quota, tilt = index(quota), CF_LD(tilt)
    if quota < 0 or tilt <= 0 or not np.isfinite(tilt):
        raise ValueError('Use a nonnegative quota and finite positive tilt')
    normalized = []
    minimum = maximum = 0
    for array, scale, offset in factors:
        array = np.asarray(array, dtype=CF_LD)
        offset = index(offset)
        if array.ndim != 1 or not len(array) or offset < 0 or not np.isfinite(scale):
            raise ValueError('Invalid positive polynomial factor')
        _cf_positive(array, 'scalar factor')
        peak = array.max()
        reduced = array/peak
        _cf_positive(reduced, 'scalar factor normalization')
        normalized.append((reduced, CF_LD(scale)+np.log2(peak), offset))
        minimum += offset
        maximum += offset+len(array)-1
    if not minimum <= quota <= maximum:
        return CF_LD('-inf')
    current, scale, low, high = np.ones(1, dtype=CF_LD), CF_LD(0), 0, 0
    remaining_min, remaining_max = minimum, maximum
    for local, local_scale, offset in normalized:
        local_high = offset+len(local)-1
        remaining_min -= offset
        remaining_max -= local_high
        convolution_low = low+offset
        next_low = max(convolution_low, quota-remaining_max)
        next_high = min(high+local_high, quota-remaining_min)
        if next_low > next_high:
            return CF_LD('-inf')
        begin, end = next_low-convolution_low, next_high-convolution_low+1
        source_lo = max(0, begin-len(local)+1)
        source_hi = min(len(current)-1, end-1)
        local_lo = max(0, begin-len(current)+1)
        local_hi = min(len(local)-1, end-1)
        if current[source_lo:source_hi+1].min()*local[local_lo:local_hi+1].min() <= 0:
            raise ArithmeticError('A retained positive scalar product would underflow')
        product = np.convolve(current, local)[begin:end]
        _cf_positive(product, 'retained scalar convolution')
        peak = product.max()
        current = product/peak
        _cf_positive(current, 'scalar convolution normalization')
        scale += local_scale+np.log2(peak)
        low, high = next_low, next_high
    if low != quota or high != quota or len(current) != 1:
        raise ArithmeticError('Scalar degree band did not reach the coefficient')
    return scale+np.log2(current[0])-quota*np.log2(tilt)


def cover_hole_constant_path(prepared, population, region_updates=16, quota=4096,
                             regions=256, threshold=104860, tilt=2, weights=None):
    """Evaluate the constant-a path, including conditional/global normalizers.

    population counts hole QUARTETS, so the hole-slot count is4*population.
    The returned floors are algebraic path lower bounds evaluated in FLOAT.
    """
    width = prepared['width']
    a, steps, quota, regions = map(index, (population, region_updates, quota, regions))
    threshold = index(threshold)
    if not 0 <= a <= width//4 or steps < 2 or regions < 1 or not 0 <= quota <= width*steps:
        raise ValueError('Invalid population, regional length, quota, or region count')
    if threshold < 0:
        raise ValueError('The deficit threshold must be nonnegative')
    if width == 512 and np.finfo(CF_LD).maxexp < 16000:
        raise RuntimeError('Production path control requires extended longdouble')
    weights = _cf_dependency('cover_hole_quota_weights')(prepared, tilt) if weights is None else weights
    if weights['width'] != width or weights['y'] != prepared['y'] or weights['tilt'] != CF_LD(tilt):
        raise ValueError('Prepared polynomial weights do not match the path')
    shift = weights['starts'][a]
    factor = (weights['coefficients'][a, shift:], weights['scale'], shift)
    labels = np.asarray([CF_LD(comb(width, j))*CF_LD(tilt)**j for j in range(width+1)])
    _cf_positive(labels, 'first-layer labels')
    log_bulk = _cf_coefficient_log2([factor]*steps, quota, tilt)
    log_end = _cf_coefficient_log2([factor]*(steps-1)+[(labels, CF_LD(0), 0)], quota, tilt)
    log_y = CF_LD(prepared['log_y'])
    log_overlap = CF_LD(prepared['overlap_scale'])+np.log2(prepared['overlap'][a, a])
    total_updates = steps*regions
    # This expression also handles regions=1; there is still a first overlap.
    raw = 4*a*log_y+(total_updates-2)*log_overlap+log_end
    if regions > 1:
        raw += (regions-1)*log_bulk
    slots = width*steps
    log_choose = CF_LD((lgamma(slots+1)-lgamma(quota+1)-lgamma(slots-quota+1))/log(2))
    conditional = raw+width*log_y-regions*log_choose
    tail = conditional+log_choose-threshold*log_y
    return dict(hole_quartets=a, hole_slots=4*a, width=width,
                region_updates=steps, regions=regions, quota=quota, threshold=threshold,
                y=float(prepared['y']), tilt=float(tilt),
                log2_overlap_aa=float(log_overlap),
                log2_bulk_scalar_coefficient=float(log_bulk),
                log2_endpoint_scalar_coefficient=float(log_end),
                log2_first_bulk_path=float(4*a*log_y+(steps-1)*log_overlap+log_bulk),
                log2_repeated_bulk_path=float(steps*log_overlap+log_bulk),
                log2_endpoint_path=float((steps-1)*log_overlap+log_end),
                log2_raw_coefficient_path_floor=float(raw),
                log2_conditional_comparison_path_floor=float(conditional),
                log2_H_unioned_comparison_tail_floor=float(tail),
                overlap_count=total_updates-2, bulk_coefficient_count=regions-1,
                exact_quota_algebra_with_FLOAT_arithmetic=True,
                actual_failure_probability_lower_bound=False,
                directed_lower_certificate=False, floating_diagnostic_only=True)


def cover_hole_constant_path_probe(y=2, populations=(96, 104, 112, 120), tilt=2):
    prepared = _cf_dependency('cover_hole_prepare')(y, 512)
    weights = _cf_dependency('cover_hole_quota_weights')(prepared, tilt)
    rows = []
    for population in populations:
        started = time.monotonic()
        result = cover_hole_constant_path(prepared, population, tilt=tilt, weights=weights)
        result['seconds'] = time.monotonic()-started
        print(result, flush=True)
        rows.append(result)
    return rows


def _cf_fraction_log2_floor(value):
    """Exact integer logarithm bound, independent of all FLOAT path routines."""
    numerator, denominator = value.numerator, value.denominator
    if numerator <= 0:
        raise ValueError('Use a positive rational value')
    exponent = numerator.bit_length()-denominator.bit_length()
    if (numerator < denominator << exponent if exponent >= 0
            else numerator << (-exponent) < denominator):
        exponent -= 1
    return exponent


def _cf_histogram_term(coefficients, histogram):
    """One positive multinomial term, evaluated with exact rational arithmetic."""
    from fractions import Fraction
    if len(coefficients) != len(histogram) or any(count < 0 for count in histogram):
        raise ValueError('Invalid polynomial histogram')
    remaining, value = sum(histogram), Fraction(1)
    for coefficient, count in zip(coefficients, histogram):
        value *= comb(remaining, count)*coefficient**count
        remaining -= count
    return value


def cover_hole_exact_pair_floor_witness():
    """Certify a y2 pair-cover comparison floor >=2^42132, with no floats.

    Choose96 hole quartets at every real layer and one overlap term I=300.
    For each exact regional coefficient retain one label histogram only.
    This proves a limitation of the pair-cover comparison at y2 for every
    potential/precision. It does not lower-bound the actual failure event.
    """
    from fractions import Fraction
    g, p = _cf_dependency('cover_hole_probabilities_exact')(Fraction(2))
    assert g == (Fraction(1), Fraction(127, 511), Fraction(21, 341),
                 Fraction(31, 2047), Fraction(1, 273))
    G = tuple(comb(4, j)*g[j] for j in range(5))
    P = tuple(comb(4, j)*p[j] for j in range(5))
    bulk_holes, bulk_selected = (178, 508, 544, 260, 46), (0, 28, 106, 216, 162)
    end_holes, end_selected = (167, 476, 510, 244, 43), (0, 26, 100, 202, 152)
    for steps, holes, selected, label_degree in (
            (16, bulk_holes, bulk_selected, 0), (15, end_holes, end_selected, 256)):
        assert sum(holes) == 96*steps and sum(selected) == 32*steps
        assert sum(j*(holes[j]+selected[j]) for j in range(5))+label_degree == 4096
    overlap = Fraction(comb(384, 300)*comb(128, 84)*2**300, comb(512, 384))
    bulk = (comb(128, 96)**16*_cf_histogram_term(G, bulk_holes)
            *_cf_histogram_term(P, bulk_selected))
    endpoint = (comb(128, 96)**15*_cf_histogram_term(G, end_holes)
                *_cf_histogram_term(P, end_selected)*comb(512, 256))
    overlap_bits = _cf_fraction_log2_floor(overlap)
    bulk_bits = _cf_fraction_log2_floor(bulk)
    endpoint_bits = _cf_fraction_log2_floor(endpoint)
    quota_bits_upper = comb(8192, 4096).bit_length()
    assert overlap_bits >= 291 and bulk_bits >= 4070 and endpoint_bits >= 4322
    assert quota_bits_upper <= 8186
    exponent = 896-104860+4094*overlap_bits+255*bulk_bits+endpoint_bits-255*quota_bits_upper
    assert exponent >= 42132
    return dict(y=2, hole_quartets=96, retained_overlap=300,
                overlap_log2_lower=overlap_bits, bulk_coefficient_log2_lower=bulk_bits,
                endpoint_coefficient_log2_lower=endpoint_bits,
                quota_binomial_log2_upper=quota_bits_upper,
                H_unioned_pair_comparison_tail_log2_lower=exponent,
                certified_lower_exponent=42132, exact_integer_rational_certificate=True,
                independent_of_potential_and_precision=True,
                only_y2_pair_cover_comparison=True,
                exact_own_cover_obstruction=False,
                actual_failure_probability_lower_bound=False)


def cover_hole_constant_path_checks():
    from fractions import Fraction
    product = _cf_dependency('_ch_poly_product')
    power = _cf_dependency('_ch_poly_power')
    overlap = _cf_dependency('_ch_overlap_exact')
    coefficient_checks = entry_checks = path_checks = tilt_checks = 0
    for width in (4, 8):
        groups = width//4
        for y in (Fraction(2), Fraction(4)):
            prepared = _cf_dependency('cover_hole_prepare')(float(y), width)
            kernel = _cf_dependency('_ch_kernel_polynomial')(width, y)
            labels = [Fraction(comb(width, degree)) for degree in range(width+1)]
            for steps in (2, 3):
                quota = width
                # Independent exact regional matrix and free endpoint, retaining all paths.
                size = groups+1
                matrix_poly = [[[Fraction(int(a == b))] for b in range(size)] for a in range(size)]
                endpoint_poly = [[Fraction(1)] for _ in range(size)]
                for completed in range(steps):
                    after = [[[Fraction(0)]*(width*(completed+1)+1) for _ in range(size)] for _ in range(size)]
                    for a in range(size):
                        for b in range(size):
                            for c in range(size):
                                term = product(kernel[a][c], matrix_poly[c][b])
                                for degree, value in enumerate(term):
                                    after[a][b][degree] += value
                    matrix_poly = after
                    if completed < steps-1:
                        endpoint_after = [[Fraction(0)]*(width*(completed+1)+1) for _ in range(size)]
                        for a in range(size):
                            for b in range(size):
                                for degree, value in enumerate(product(kernel[a][b], endpoint_poly[b])):
                                    endpoint_after[a][degree] += value
                        endpoint_poly = endpoint_after
                B = [[poly[quota] for poly in row] for row in matrix_poly]
                q = [product(poly, labels)[quota] for poly in endpoint_poly]
                for a in range(size):
                    w = kernel[0][a]  # H_0a=1, independently from exact HG overlaps.
                    c = power(w, steps)[quota]
                    e = product(power(w, steps-1), labels)[quota]
                    haa = overlap(width, a, a, y)
                    first = y**(4*a)*haa**(steps-1)*c
                    middle = haa**steps*c
                    last = haa**(steps-1)*e
                    assert first <= B[groups][a] and middle <= B[a][a] and last <= q[a]
                    entry_checks += 3
                    for regions in (1, 2, 3):
                        path = y**(4*a)*haa**(steps*regions-2)*c**(regions-1)*e
                        vector = q
                        for _ in range(regions-1):
                            vector = [sum(B[i][j]*vector[j] for j in range(size)) for i in range(size)]
                        assert path <= vector[groups]
                        if regions >= 2:
                            assert path == first*middle**(regions-2)*last
                        path_checks += 1
                    for tilt in (1, 2, 3):
                        weights = _cf_dependency('cover_hole_quota_weights')(prepared, tilt)
                        shift = weights['starts'][a]
                        local = (weights['coefficients'][a, shift:], weights['scale'], shift)
                        found_c = _cf_coefficient_log2([local]*steps, quota, tilt)
                        label_float = np.asarray([CF_LD(v.numerator)/v.denominator*CF_LD(tilt)**j for j, v in enumerate(labels)])
                        found_e = _cf_coefficient_log2([local]*(steps-1)+[(label_float, CF_LD(0), 0)], quota, tilt)
                        for found, exact in ((found_c, c), (found_e, e)):
                            if exact:
                                expected = CF_LD(exact.numerator)/exact.denominator
                                assert abs(np.exp2(found)/expected-1) < CF_LD('1e-11')
                            else:
                                assert np.isneginf(found)
                            coefficient_checks += 1
                        result = cover_hole_constant_path(prepared, a, steps, quota, 3, 1, tilt, weights)
                        exact_path = y**(4*a)*haa**(3*steps-2)*c**2*e
                        if exact_path:
                            expected = np.log2(CF_LD(exact_path.numerator)/exact_path.denominator)
                            assert abs(CF_LD(result['log2_raw_coefficient_path_floor'])-expected) < CF_LD('1e-10')
                        else:
                            assert np.isneginf(result['log2_raw_coefficient_path_floor'])
                        tilt_checks += 1
    # Structural impossibility stays zero; it is not a numerical floor.
    assert np.isneginf(_cf_coefficient_log2([(np.ones(2), 0, 3)], 2, 1))
    print(dict(exact_scalar_coefficient_comparisons=coefficient_checks,
               exact_regional_path_entry_checks=entry_checks,
               exact_full_moment_path_checks=path_checks,
               common_tilt_and_formula_checks=tilt_checks,
               structural_zero_checks=1,
               exact_pair_floor_witness=cover_hole_exact_pair_floor_witness(),
               path_diagnostic_FLOAT_but_witness_exact=True), flush=True)


def cover_hole_constant_path_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--y', type=float, default=2.)
    parser.add_argument('--populations', nargs='+', type=int, default=(96, 104, 112, 120))
    parser.add_argument('--tilt', type=float, default=2.)
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
            cover_hole_constant_path_checks()
        else:
            cover_hole_constant_path_probe(args.y, args.populations, args.tilt)
    finally:
        timer.cancel()


if __name__ == '__main__':
    cover_hole_constant_path_main()
