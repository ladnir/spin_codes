r"""Positive129-state global control for the normalized dual-Hall random cover.

At a quartet with j non-H inputs, select it with the asymmetric cover
probability p_j=(1-y^(-2j))/(1-y^(-(8+j))) for j>0. Its hole probability
is g_j=y^(-2j)*(1-y^(-(8-j)))/(1-y^(-(8+j))); g0=1,p0=0.
Use the positive polynomials G(x)=sum C4j*g_j*x^j and
P(x)=sum_(j>=1) C4j*p_j*x^j. No subtractive evaluation of P is used.

For V_t the union of hole quartets, the number of uncovered sources is
at most n+|V_T|+sum_(u=2)^(T-1)|V_u intersect V_(u+1)|. Conditional on
hole count b, its union is a uniform4b-subset. In the common physical
hole-population basis a,b=0..n/4, put
 H_ab=E[y^|A intersect B|], |A|=4a,|B|=4b,
 K_x[a,b]=H_ab*C(n/4,b)*G(x)^b*P(x)^(n/4-b).

Reverse time starts at a=n/4. Its first overlap charges |V_T|. The
external y^n charges the final source layer ONCE in the full4096-update
calculation. Input layer1 has no cover coin; it supplies (1+x)^n.

For M=n*16, quota h=M/2, C=C(M,h), define
 B=[x^h]K_x^16, q=[x^h](1+x)^n*K_x^15*1.
The exact quota-conditioned cover moment is bounded by
 y^n*C^-256*e_(n/4)*B^255*q. Replace B and q by entrywise minima of
positive Cauchy evaluations x^-h AFTER each full composition. One H-set
union and the threshold delta>=104860 give the tail upper bound
 y^(512-104860)*C^-255*e_128*B_upper^255*q_upper.

This is a FLOAT diagnostic with exact small rational checks. It does
not assemble the larger dual-Hall operator or claim a global certificate.
The production probe always uses4096 updates, not repeated64-windows.
Parent owns every numerical job. No jobs or sampling run on import.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, lgamma, log
from operator import index
import os
import shutil
import threading
import time

for _ch_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_ch_thread_var] = '1'

import numpy as np


CH_LD = np.longdouble
CH_WIDTH, CH_UPDATES, CH_REGION, CH_QUOTA, CH_THRESHOLD = 512, 4096, 16, 4096, 104860


def _ch_dependency(name):
    if name in globals():
        return globals()[name]
    if name == 'temporal_selector_overlap':
        import quartet_selector_temporal as module
    elif name.startswith('_rr_'):
        import quartet_regional_rank_cap as module
    else:
        import quartet_selector_hit_moment as module
    return getattr(module, name)


def cover_hole_probabilities_exact(y):
    y = Fraction(y)
    if y <= 1:
        raise ValueError('Require y>1')
    holes, selected = [Fraction(1)], [Fraction(0)]
    for j in range(1, 5):
        A, B = y**(-(8-j)), y**(-2*j)
        holes.append(B*(1-A)/(1-A*B))
        selected.append((1-B)/(1-A*B))
    return tuple(holes), tuple(selected)


def cover_hole_prepare(y=2, width=512):
    """Reuse the positive overlap Gram construction once for all x tilts."""
    width, y = index(width), CH_LD(y)
    if width < 4 or width > 512 or width % 4 or not np.isfinite(y) or y <= 1:
        raise ValueError('Require quartet width4..512 and finite y>1')
    log_y, ln_y = np.log2(y), np.log(y)
    holes, selected = np.ones(5, dtype=CH_LD), np.zeros(5, dtype=CH_LD)
    for j in range(1, 5):
        denominator = -np.expm1(-(8+j)*ln_y)
        holes[j] = np.exp(-2*j*ln_y)*(-np.expm1(-(8-j)*ln_y))/denominator
        selected[j] = -np.expm1(-2*j*ln_y)/denominator
    if (np.any(holes <= 0) or np.any(selected[1:] <= 0)
            or not np.all(np.isfinite(holes)) or not np.all(np.isfinite(selected))):
        raise ArithmeticError('Cover probabilities underflowed or overflowed')
    overlap, scale = _ch_dependency('temporal_selector_overlap')(log_y, width)
    coefficients = np.asarray([comb(4, j) for j in range(5)], dtype=CH_LD)
    return dict(width=width, y=y, log_y=log_y, overlap=overlap, overlap_scale=scale,
                G_coefficients=coefficients*holes, P_coefficients=coefficients*selected)


def cover_hole_kernel(prepared, x):
    """Positive scalar polynomial evaluations and column population weights."""
    x = CH_LD(x)
    if not np.isfinite(x) or x <= 0:
        raise ValueError('Cauchy tilt x must be finite and positive')
    G, P = CH_LD(0), CH_LD(0)
    gc, pc = prepared['G_coefficients'], prepared['P_coefficients']
    for j in (4, 3, 2, 1, 0):
        G, P = G*x+gc[j], P*x+pc[j]
    if min(G, P) <= 0 or not np.isfinite(G+P):
        raise ArithmeticError('Positive cover polynomial evaluation became invalid')
    groups = prepared['width']//4
    logs = np.asarray([np.log2(CH_LD(comb(groups, b)))+b*np.log2(G)+(groups-b)*np.log2(P)
                       for b in range(groups+1)], dtype=CH_LD)
    weight_scale = logs.max()
    weights = np.exp2(logs-weight_scale)
    if np.any(weights <= 0) or not np.all(np.isfinite(weights)):
        raise ArithmeticError('Positive hole-population weight disappeared')
    matrix, extra = _ch_dependency('_sh_normalize')(prepared['overlap']*weights[None, :])
    return matrix, prepared['overlap_scale']+weight_scale+extra


def cover_hole_regional(prepared, x_values=None, region_updates=16, quota=4096, progress=False):
    """Minimize Cauchy upper operators only after complete regional composition."""
    width = prepared['width']
    region_updates, quota = index(region_updates), index(quota)
    if region_updates < 2 or not 0 <= quota <= width*region_updates:
        raise ValueError('Invalid region length or quota')
    if x_values is None:
        y = prepared['y']
        x_values = (CH_LD(1), y, y*y)
    x_values = tuple(dict.fromkeys(map(CH_LD, x_values)))
    if not 1 <= len(x_values) <= 9 or any(not np.isfinite(x) or x <= 0 for x in x_values):
        raise ValueError('Use one to nine positive finite Cauchy tilts')
    ordinary = endpoint = None
    ordinary_scale = endpoint_scale = CH_LD(0)
    summaries = []
    for x in x_values:
        matrix, scale = cover_hole_kernel(prepared, x)
        data = dict(kernel=matrix, kernel_scale=scale,
                    terminal=np.ones(width//4+1, dtype=CH_LD),
                    terminal_scale=width*np.log2(1+x))
        bulk, bs, last, ls = _ch_dependency('_rr_region_blocks')(data, region_updates)
        correction = -quota*np.log2(x)
        bs, ls = bs+correction, ls+correction
        if ordinary is None:
            ordinary, ordinary_scale, endpoint, endpoint_scale = bulk, bs, last, ls
        else:
            ordinary, ordinary_scale = _ch_dependency('_sh_minimum')(ordinary, ordinary_scale, bulk, bs)
            endpoint, endpoint_scale = _ch_dependency('_sh_minimum')(endpoint, endpoint_scale, last, ls)
        summaries.append(dict(x=float(x), bulk_scale=float(bs), endpoint_scale=float(ls)))
        if progress:
            print(dict(stage='cover_hole_regional_tilt', **summaries[-1]), flush=True)
    return dict(ordinary=ordinary, ordinary_scale=ordinary_scale,
                endpoint=endpoint, endpoint_scale=endpoint_scale,
                x_values=x_values, tilt_summaries=summaries,
                width=width, region_updates=region_updates, quota=quota)


def cover_hole_raw_moment(regional, regions):
    """Raw coefficient upper moment from the FULL-hole start population."""
    regions = index(regions)
    if regions < 1:
        raise ValueError('At least one region is required')
    # Existing positive endpoint power starts at index0. One fixed basis
    # permutation makes index0 mean full holes; it changes no physical min.
    matrix = np.ascontiguousarray(regional['ordinary'][::-1, ::-1])
    endpoint = np.ascontiguousarray(regional['endpoint'][::-1])
    return _ch_dependency('_rr_power_endpoint')(matrix, regional['ordinary_scale'],
        endpoint, regional['endpoint_scale'], regions-1)


def cover_hole_bound(y=2, x_values=None, progress=False):
    """Full4096-update half-row class; no window reset or repeated endpoint loss."""
    prepared = cover_hole_prepare(y, CH_WIDTH)
    regional = cover_hole_regional(prepared, x_values, CH_REGION, CH_QUOTA, progress)
    raw = CH_LD(cover_hole_raw_moment(regional, CH_UPDATES//CH_REGION))
    M, regions = CH_WIDTH*CH_REGION, CH_UPDATES//CH_REGION
    log_C = CH_LD((lgamma(M+1)-lgamma(CH_QUOTA+1)-lgamma(M-CH_QUOTA+1))/log(2))
    endpoint = CH_WIDTH*prepared['log_y']
    conditional = raw+endpoint-regions*log_C
    tail = conditional+log_C-CH_THRESHOLD*prepared['log_y']
    return dict(y=float(prepared['y']), x_values=tuple(map(float, regional['x_values'])),
                width=CH_WIDTH, updates=CH_UPDATES, regions=regions, updates_per_region=CH_REGION,
                quota_nonH_per_region=CH_QUOTA, integer_deficit_threshold=CH_THRESHOLD,
                log2_raw_coefficient_moment=float(raw), global_endpoint_log2=float(endpoint),
                log2_quota_binomial=float(log_C), log2_H_union=float(log_C),
                log2_conditional_cover_moment_upper=float(conditional),
                log2_probability_upper=float(tail), bulk_regions=regions-1,
                endpoint_cover_updates=CH_REGION-1, first_input_labels_counted=True,
                full_hole_initial_population=CH_WIDTH//4,
                endpoint_charged_once_globally=True, common_physical_basis=True,
                minimum_after_regional_composition=True,
                floating_diagnostic_only=True, marker_optimality_claim=False,
                full_code_certificate=False)


def _ch_choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def _ch_overlap_exact(width, a, b, y):
    return sum((Fraction(_ch_choose(4*a, ell)*_ch_choose(width-4*a, 4*b-ell), comb(width, 4*b))*y**ell
                for ell in range(max(0, 4*b-(width-4*a)), min(4*a, 4*b)+1)), Fraction(0))


def _ch_poly_product(left, right):
    result = [Fraction(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    result[i+j] += a*b
    return result


def _ch_poly_power(poly, count):
    result = [Fraction(1)]
    for _ in range(count):
        result = _ch_poly_product(result, poly)
    return result


def _ch_kernel_polynomial(width, y):
    holes, selected = cover_hole_probabilities_exact(y)
    G = [comb(4, j)*holes[j] for j in range(5)]
    P = [comb(4, j)*selected[j] for j in range(5)]
    groups = width//4
    weights = [_ch_poly_product(_ch_poly_power(G, b), _ch_poly_power(P, groups-b)) for b in range(groups+1)]
    return [[[coefficient*comb(groups, b)*_ch_overlap_exact(width, a, b, y) for coefficient in weights[b]]
             for b in range(groups+1)] for a in range(groups+1)]


def _ch_evaluate(poly, x):
    result = Fraction(0)
    for value in reversed(poly):
        result = result*x+value
    return result


def _ch_matmul(left, right):
    size = len(left)
    return [[sum((left[i][k]*right[k][j] for k in range(size)), Fraction(0)) for j in range(size)]
            for i in range(size)]


def _ch_partitions(width):
    if width == 4:
        return ((15,),)
    return tuple((mask, 255 ^ mask) for triple in combinations(range(1, 8), 3)
                 for mask in (1+sum(1 << coordinate for coordinate in triple),))


def _ch_physical_kernel_polynomial(width, y):
    holes, selected = cover_hole_probabilities_exact(y)
    groups = width//4
    result = [[[Fraction(0)]*(width+1) for _ in range(groups+1)] for _ in range(groups+1)]
    partitions = _ch_partitions(width)
    for partition in partitions:
        for J in range(1 << width):
            labels = [(J & quartet).bit_count() for quartet in partition]
            for choice in range(1 << groups):
                V, probability = 0, Fraction(1)
                for group, quartet in enumerate(partition):
                    is_hole = (choice >> group) & 1
                    probability *= holes[labels[group]] if is_hole else selected[labels[group]]
                    if is_hole:
                        V |= quartet
                if not probability:
                    continue
                b, degree = V.bit_count()//4, J.bit_count()
                for a in range(groups+1):
                    old = (1 << (4*a))-1
                    result[a][b][degree] += probability*y**(old & V).bit_count()
    count = len(partitions)
    return [[[value/count for value in poly] for poly in row] for row in result]


def _ch_small_region_exact(y):
    """Width4, two-update region, quota4; exact polynomial coefficients."""
    K = _ch_kernel_polynomial(4, y)
    bulk = [[Fraction(0)]*2 for _ in range(2)]
    end = [Fraction(0)]*2
    labels = [Fraction(comb(4, j)) for j in range(5)]
    for a in range(2):
        for b in range(2):
            for c in range(2):
                bulk[a][b] += _ch_poly_product(K[a][c], K[c][b])[4]
            end[a] += _ch_poly_product(labels, K[a][b])[4]
    return K, bulk, end


def _ch_cover_sum_width4(J, y):
    """Explicit cover coins, exact missing sources, and overlap upper moment."""
    holes, selected = cover_hole_probabilities_exact(y)
    updates = len(J)
    exact = upper = total_probability = Fraction(0)
    missing_checks = 0
    for choice in range(1 << (updates-1)):
        C, V, probability = [J[0]]+[0]*(updates-1), [0]*updates, Fraction(1)
        for t in range(1, updates):
            j = J[t].bit_count()
            is_hole = (choice >> (t-1)) & 1
            probability *= holes[j] if is_hole else selected[j]
            V[t] = 15 if is_hole else 0
            C[t] = 0 if is_hole else J[t]
        if not probability:
            continue
        predecessors = list(C)
        for t in range(1, updates):
            if C[t]:
                predecessors[t-1] |= 15
                if t >= 2:
                    predecessors[t-2] |= 15
        missing = 4*updates-sum(mask.bit_count() for mask in predecessors)
        bound = 4
        old = 15
        for t in range(updates-1, 0, -1):
            bound += (old & V[t]).bit_count()
            old = V[t]
        assert missing <= bound
        exact += probability*y**missing
        upper += probability*y**bound
        total_probability += probability
        missing_checks += 1
    assert total_probability == 1 and exact <= upper
    return upper, missing_checks


def _ch_fraction_array(matrix):
    return np.asarray([[CH_LD(value.numerator)/value.denominator for value in row] for row in matrix], dtype=CH_LD)


def cover_hole_checks():
    probability_checks = overlap_checks = coefficient_checks = missing_checks = 0
    quota_checks = cauchy_checks = float_checks = 0
    for y in (Fraction(2), Fraction(4)):
        holes, selected = cover_hole_probabilities_exact(y)
        for j in range(5):
            assert holes[j]+selected[j] == 1
            assert holes[j] > 0 and (selected[j] > 0 if j else selected[j] == 0)
            probability_checks += 1
        for width in (4, 8):
            groups = width//4
            for a in range(groups+1):
                old = (1 << (4*a))-1
                for b in range(groups+1):
                    total = Fraction(0)
                    for positions in combinations(range(width), 4*b):
                        new = sum(1 << i for i in positions)
                        total += y**(old & new).bit_count()
                    assert total/comb(width, 4*b) == _ch_overlap_exact(width, a, b, y)
                    overlap_checks += 1
            formula = _ch_kernel_polynomial(width, y)
            physical = _ch_physical_kernel_polynomial(width, y)
            assert formula == physical
            coefficient_checks += (groups+1)**2*(width+1)
            prepared = cover_hole_prepare(float(y), width)
            for x in (Fraction(1), y, y*y):
                exact = [[_ch_evaluate(poly, x) for poly in row] for row in formula]
                kernel, scale = cover_hole_kernel(prepared, float(x))
                reference = _ch_fraction_array(exact)
                assert np.max(np.abs(kernel*np.exp2(scale)-reference)/reference) < CH_LD('1e-12')
                float_checks += (groups+1)**2
        K, bulk, end = _ch_small_region_exact(y)
        minimum_bulk = minimum_end = None
        for x in (Fraction(1), y, y*y):
            evaluated = [[_ch_evaluate(poly, x) for poly in row] for row in K]
            square = _ch_matmul(evaluated, evaluated)
            B = [[value/x**4 for value in row] for row in square]
            q = [(1+x)**4*sum(row)/x**4 for row in evaluated]
            for a in range(2):
                assert end[a] <= q[a]
                for b in range(2):
                    assert bulk[a][b] <= B[a][b]
                    cauchy_checks += 1
                cauchy_checks += 1
            minimum_bulk = B if minimum_bulk is None else [[min(minimum_bulk[a][b], B[a][b]) for b in range(2)] for a in range(2)]
            minimum_end = q if minimum_end is None else [min(minimum_end[a], q[a]) for a in range(2)]
        exact_raw = sum(bulk[1][b]*end[b] for b in range(2))
        comparison_raw = sum(minimum_bulk[1][b]*minimum_end[b] for b in range(2))
        assert exact_raw <= comparison_raw
        # Independent physical four-layer enumeration with two exact half quotas.
        region_labels = [(mask & 15, mask >> 4) for positions in combinations(range(8), 4)
                         for mask in (sum(1 << position for position in positions),)]
        physical_raw = Fraction(0)
        for first, second in product(region_labels, repeat=2):
            value, checked = _ch_cover_sum_width4(first+second, y)
            physical_raw += value
            missing_checks += checked
        assert physical_raw == y**4*exact_raw
        quota_checks += 1
        prepared = cover_hole_prepare(float(y), 4)
        regional = cover_hole_regional(prepared, (1, float(y), float(y*y)), 2, 4)
        reference = _ch_fraction_array(minimum_bulk)
        assert np.max(np.abs(regional['ordinary']*np.exp2(regional['ordinary_scale'])-reference)/reference) < CH_LD('1e-12')
        reference_end = np.asarray([CH_LD(value.numerator)/value.denominator for value in minimum_end])
        assert np.max(np.abs(regional['endpoint']*np.exp2(regional['endpoint_scale'])-reference_end)/reference_end) < CH_LD('1e-12')
        value = cover_hole_raw_moment(regional, 2)
        expected = np.log2(CH_LD(comparison_raw.numerator)/comparison_raw.denominator)
        assert abs(value-expected) < CH_LD('1e-10')
        float_checks += 7
        # First-update labels are free; no first-update cover coin exists.
        for J in range(16):
            upper, checked = _ch_cover_sum_width4((J,), y)
            assert upper == y**4
            missing_checks += checked
    print(dict(exact_cover_probability_checks=probability_checks,
               exact_overlap_checks=overlap_checks,
               exact_physical_kernel_coefficients=coefficient_checks,
               exact_missing_source_checks=missing_checks,
               exact_regional_quota_compositions=quota_checks,
               exact_regional_cauchy_checks=cauchy_checks,
               floating_toy_checks=float_checks, rational_markers=(2, 4),
               global_endpoint_charged_once=True, full_code_certificate=False), flush=True)


def cover_hole_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--y', type=float, choices=(2., 4.), default=2.)
    parser.add_argument('--x-values', nargs='+', type=float)
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
            cover_hole_checks()
        else:
            result = cover_hole_bound(args.y, args.x_values, progress=True)
            print(dict(seconds=round(time.monotonic()-started, 3), **result), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    cover_hole_main()
