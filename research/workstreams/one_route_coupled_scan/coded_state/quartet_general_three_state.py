"""Bounded general-network gate on the three pure population states.

Retain O=(0,0), P=(n,0), H=(n,n) at EVERY update. This gives the exact
positive restriction L_x of the full adaptive population operator K_x.
Local coefficients come directly from the general loss/gain/cost rule.
The default proof network is(2,1,4), y=3/2, integer threshold184116.
Neither the encoder nor its routing distribution changes.

The rowwise comparison is inf_x x^-h L_x^16 v/C0. A subeigenvector for
this restricted map extends by zero to the full physical map, so its
growth lower bound obstructs FULL sourcewise Cauchy supersolutions.
Numerical derivative brackets and convex tangent bounds are FLOAT-only;
they nominate an obstruction for rational replay, not a certificate.
The corresponding restricted UPPER bound controls only paths confined
to O/P/H. It is not an upper bound for the full operator.

The restricted free endpoint is[x^h]L_x^15 q_x/C0, with
 q_x=((1+x)^n,(y^(a-b)+x)^n,(y^a+x)^n).
The full endpoint is at least1 because its cover cost is nonnegative.
There is one C0 row-set union and no additional y^n endpoint factor.
Default grid(1,y^a,y^(2a)) and at most12 continuous half-updates.
No numerical jobs or arrays on import; parent owns all execution.
"""

import argparse
from fractions import Fraction as F
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _gn_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_gn_thread_var] = '1'

import numpy as np

GN_LD = np.longdouble
GN_EDGES = ((1, 0, 0), (0, 0, 1), (3, 1, 0), (2, 1, 2), (5, 2, 0), (4, 2, 2))


def _gn_logfraction(value):
    value = F(value)
    return np.log2(GN_LD(value.numerator)/GN_LD(value.denominator))


def _gn_local(y, capacities, u, w):
    source, temporal, sink = capacities
    hole, selected = [F(0)]*5, [F(0)]*5
    for r in range(w+1):
        for s in range(u-w+1):
            for t in range(5-u):
                j = r+s+t
                gain = sink*j-(source-temporal)*(r+s)-temporal*r
                loss = max(0, 4*source+temporal*u-gain)
                if gain == 0:
                    g, p = F(1), F(0)
                elif loss == 0:
                    g, p = F(0), F(1)
                else:
                    denominator = y**(loss+gain)-1
                    g, p = (y**loss-1)/denominator, (y**(loss+gain)-y**loss)/denominator
                cost = (source-temporal)*(u-r-s)+temporal*(w-r)
                weight = comb(w, r)*comb(u-w, s)*comb(4-u, t)*y**cost
                hole[j] += weight*g
                selected[j] += weight*p
    return tuple(hole), tuple(selected)


def general_three_state_prepare(capacities=(2, 1, 4), y=F(3, 2), width=512, updates=16):
    capacities = tuple(map(index, capacities))
    source, temporal, sink = capacities
    y, width, updates = F(y), index(width), index(updates)
    if (source < 1 or not 0 <= temporal <= source or sink < source+temporal or y <= 1
            or width < 4 or width > 512 or width % 4 or updates < 2 or (width*updates) % 2):
        raise ValueError('Invalid general network, marker, quartet width, or region')
    coefficients = sum((_gn_local(y, capacities, *state) for state in ((0, 0), (4, 0), (4, 4))), ())
    logs = np.asarray([[_gn_logfraction(value) if value else GN_LD('-inf') for value in row]
                       for row in coefficients], dtype=GN_LD)
    quota = width*updates//2
    normalizer = comb(width*updates, quota)
    r0, short = y**source, y**(source-temporal)
    alpha, beta = (short+r0)/(1+r0), 2*r0/(short+r0)
    initial = np.asarray((0, width*_gn_logfraction(alpha), width*_gn_logfraction(alpha*beta)), dtype=GN_LD)
    return dict(capacities=capacities, y=y, log_y=_gn_logfraction(y), width=width,
                groups=width//4, updates=updates, quota=quota, C=normalizer,
                log_C=np.log2(GN_LD(normalizer)), coefficients=coefficients,
                coefficient_logs=logs, q_bases=(F(1), short, r0),
                initial=initial, gauge=(alpha, beta), grid=(F(1), r0, r0*r0))


def _gn_matrix(cache, log_x):
    logs = cache['coefficient_logs']+np.arange(5, dtype=GN_LD)*log_x
    peaks = logs.max(axis=1)
    values = np.exp2(logs-peaks[:, None])
    if np.any(np.isfinite(logs) & (values <= 0)):
        raise ArithmeticError('A positive local polynomial term underflowed')
    totals = values.sum(axis=1, dtype=GN_LD)
    edge_logs = cache['groups']*(peaks+np.log2(totals))
    edge_slopes = cache['groups']*(values @ np.arange(5, dtype=GN_LD))/totals
    scale = edge_logs.max()
    matrix, derivative = np.zeros((3, 3), dtype=GN_LD), np.zeros((3, 3), dtype=GN_LD)
    for polynomial, row, column in GN_EDGES:
        value = np.exp2(edge_logs[polynomial]-scale)
        if value <= 0 or not np.isfinite(value):
            raise ArithmeticError('A positive three-state edge underflowed')
        matrix[row, column] = value
        derivative[row, column] = value*edge_slopes[polynomial]
    return matrix, derivative, scale


def _gn_product(left, right):
    a, da, sa = left
    b, db, sb = right
    if a[a > 0].min()*b[b > 0].min() <= 0:
        raise ArithmeticError('A positive matrix product term would underflow')
    value, derivative = a @ b, da @ b+a @ db
    reachable = np.zeros((3, 3), dtype=bool)
    for k in range(3):
        reachable |= (a[:, k] > 0)[:, None] & (b[k] > 0)[None, :]
    if (np.any(reachable & (value <= 0)) or not np.all(np.isfinite(value+derivative))
            or np.any(derivative < 0)):
        raise ArithmeticError('Invalid positive matrix/derivative product')
    peak = value.max()
    if value[reachable].min()/peak <= 0:
        raise ArithmeticError('Positive matrix normalization would underflow')
    return value/peak, derivative/peak, sa+sb+np.log2(peak)


def _gn_power(base, exponent):
    result = (np.eye(3, dtype=GN_LD), np.zeros((3, 3), dtype=GN_LD), GN_LD(0))
    while exponent:
        if exponent & 1:
            result = _gn_product(result, base)
        exponent >>= 1
        if exponent:
            base = _gn_product(base, base)
    return result


def _gn_evaluate(cache, log_x, log_phi, endpoint=False):
    """Return all corrected log values and x*d/dx log-polynomial slopes."""
    matrix = _gn_power(_gn_matrix(cache, log_x), cache['updates']-int(endpoint))
    p, dp, scale = matrix
    if endpoint:
        x = np.exp2(GN_LD(log_x))
        bases = np.asarray([GN_LD(v.numerator)/v.denominator for v in cache['q_bases']])
        vector_logs = cache['width']*np.log2(bases+x)
        vector_slopes = cache['width']*x/(bases+x)
    else:
        vector_logs, vector_slopes = np.asarray(log_phi, dtype=GN_LD), np.zeros(3, dtype=GN_LD)
    peak = vector_logs.max()
    vector = np.exp2(vector_logs-peak)
    if np.any(vector <= 0) or p[p > 0].min()*vector.min() <= 0:
        raise ArithmeticError('Positive vector product would underflow')
    value = p @ vector
    derivative = dp @ vector+p @ (vector*vector_slopes)
    logs = np.log2(value)+scale+peak-cache['quota']*log_x-cache['log_C']
    if not endpoint:
        logs -= log_phi
    slopes = derivative/value-cache['quota']
    if not np.all(np.isfinite(logs+slopes)):
        raise ArithmeticError('Invalid positive convex row evaluation')
    return logs, slopes


def _gn_minimize(evaluate, iterations=40):
    """FLOAT derivative bracket plus convex tangent lower and sampled upper."""
    lo, hi = GN_LD(-2), GN_LD(3)
    fl, sl = evaluate(lo)
    fh, sh = evaluate(hi)
    for _ in range(5):
        if sl < 0 < sh:
            break
        if sl >= 0:
            lo -= 2
            fl, sl = evaluate(lo)
        if sh <= 0:
            hi += 2
            fh, sh = evaluate(hi)
    if not sl < 0 < sh:
        raise ArithmeticError('Could not bracket the all-positive-x convex minimum')
    best, best_t = (fl, lo) if fl <= fh else (fh, hi)
    for _ in range(iterations):
        middle = (lo+hi)/2
        value, slope = evaluate(middle)
        if value < best:
            best, best_t = value, middle
        if slope < 0:
            lo, fl, sl = middle, value, slope
        else:
            hi, fh, sh = middle, value, slope
    lower = max(fl+sl*(hi-lo), fh-sh*(hi-lo))
    if lower > best+GN_LD('1e-7'):
        raise ArithmeticError('Convex minimum enclosure is inconsistent')
    return dict(lower=lower, upper=best, log_x=best_t, log_bracket=(lo, hi),
                slopes=(sl, sh), x=float(np.exp2(best_t)), floating_bracket_only=True)


def _gn_rows(cache, log_phi, endpoint=False):
    rows = []
    for row in range(3):
        def evaluate(log_x):
            values, slopes = _gn_evaluate(cache, log_x, log_phi, endpoint)
            return values[row], slopes[row]
        rows.append(_gn_minimize(evaluate))
    return rows


def general_three_state_control(capacities=(2, 1, 4), y=F(3, 2), threshold=184116,
                                half_updates=6, log_p=None, log_h=None, progress=False):
    started = time.monotonic()
    threshold, half_updates = index(threshold), index(half_updates)
    if threshold < 0 or not 0 <= half_updates <= 12 or ((log_p is None) != (log_h is None)):
        raise ValueError('Use nonnegative threshold, at most12 half-updates, and both optional gauge logs')
    cache = general_three_state_prepare(capacities, y)
    phi = cache['initial'].copy() if log_p is None else np.asarray((0, log_p, log_h), dtype=GN_LD)
    if not np.all(np.isfinite(phi)):
        raise ValueError('Gauge logs must be finite')
    endpoint_rows = _gn_rows(cache, phi, endpoint=True)
    endpoint_upper = np.asarray([row['upper'] for row in endpoint_rows])
    summaries = []
    for iteration in range(half_updates+1):
        grid_values = np.asarray([_gn_evaluate(cache, _gn_logfraction(x), phi)[0] for x in cache['grid']])
        grid_ratio = grid_values.min(axis=0)
        rows = _gn_rows(cache, phi)
        lower = np.asarray([row['lower'] for row in rows])
        upper = np.asarray([row['upper'] for row in rows])
        floor = lower.min()
        rho, nu = upper.max(), (endpoint_upper-phi).max()
        common = cache['log_C']-threshold*cache['log_y']
        summary = dict(iteration=iteration, log_phi=[float(v) for v in phi],
                       grid=[str(x) for x in cache['grid']],
                       grid_row_ratio_logs=[float(v) for v in grid_ratio],
                       grid_comparison_floor_display=float(grid_ratio.min()),
                       continuous_row_lower_displays=[float(v) for v in lower],
                       continuous_row_upper_displays=[float(v) for v in upper],
                       continuous_minimizers=[row['x'] for row in rows],
                       continuous_comparison_floor_display=float(floor),
                       restricted_upper_growth_display=float(rho),
                       inherited_full_sourcewise_tail_floor_display=float(common+255*floor),
                       restricted_tail_upper_display=float(common+phi[2]+255*rho+nu),
                       restricted_endpoint_upper_logs=[float(v) for v in endpoint_upper],
                       endpoint_multiplier_log=float(nu),
                       candidate_for_exact_obstruction_replay=bool(common+255*floor > 0),
                       directed_or_exact_numerical_certificate=False,
                       full_physical_upper_bound=False)
        summaries.append(summary)
        if progress:
            print(dict(stage='general_three_state', **summary), flush=True)
        phi += upper/2
        phi -= phi[0]
    return dict(capacities=cache['capacities'], y=str(cache['y']), threshold=threshold,
                initial_gauge=tuple(map(str, cache['gauge'])), completed_witnesses=len(summaries),
                necessary_bulk_budget_before_endpoint=float((threshold*cache['log_y']-cache['log_C'])/255),
                final=summaries[-1], endpoint_lower_for_full_operator=1,
                row_set_union='binom(8192,4096), once', full_start='H=(512,512)',
                extra_terminal_factor=False, fixed_encoder=True,
                seconds=time.monotonic()-started, no_go_proved=False, full_code_certificate=False)


def _gn_poly_product(left, right):
    after = [F(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            after[i+j] += a*b
    return after


def _gn_poly_power(poly, power):
    result = [F(1)]
    for _ in range(power):
        result = _gn_poly_product(result, poly)
    return result


def general_three_state_checks():
    counts = dict(local_coefficients=0, physical_entries=0, power_derivatives=0,
                  quota_upper_rows=0, endpoint_rows=0, convex_minima=0, gauge_endpoints=0)
    if '_av_physical_kernel' in globals():
        physical_kernel = globals()['_av_physical_kernel']
    else:
        from quartet_adaptive_vector import _av_physical_kernel as physical_kernel
    for capacities, y in (((1, 1, 2), F(2)), ((2, 1, 4), F(3, 2))):
        a, b, _ = capacities
        for width in (4, 8):
            cache = general_three_state_prepare(capacities, y, width, 2)
            for u, w in ((0, 0), (4, 0), (4, 4)):
                hole, selected = _gn_local(y, capacities, u, w)
                for j in range(5):
                    expected = comb(4, j)*y**(((a-b)*u+b*w)*(4-j)//4)
                    assert hole[j]+selected[j] == expected
                    counts['local_coefficients'] += 1
            for x in (F(1), y**a, y**(2*a)):
                states, full = physical_kernel(width, y, x, capacities)
                positions = [states.index(state) for state in ((0, 0), (width, 0), (width, width))]
                normalized, derivative, scale = _gn_matrix(cache, _gn_logfraction(x))
                actual = normalized*np.exp2(scale)
                for i in range(3):
                    for j in range(3):
                        exact = full[positions[i]][positions[j]]
                        expected = GN_LD(exact.numerator)/exact.denominator
                        assert abs(actual[i, j]-expected) <= GN_LD('1e-10')*max(1, expected)
                        counts['physical_entries'] += 1
                phi = np.asarray((0, 1, 3), dtype=GN_LD)
                values, slopes = _gn_evaluate(cache, _gn_logfraction(x), phi)
                # Independent rational degree-polynomial two-step expansion.
                matrix = [[[F(0)]]*3 for _ in range(3)]
                for poly, i, j in GN_EDGES:
                    matrix[i][j] = _gn_poly_power(cache['coefficients'][poly], width//4)
                for row in range(3):
                    polynomial = [F(0)]*(2*width+1)
                    for middle in range(3):
                        for column in range(3):
                            product = _gn_poly_product(matrix[row][middle], matrix[middle][column])
                            for degree, coefficient in enumerate(product):
                                polynomial[degree] += coefficient*(1 << (0, 1, 3)[column])
                    exact = sum(coefficient*x**degree for degree, coefficient in enumerate(polynomial))
                    exact_derivative = sum(degree*coefficient*x**degree for degree, coefficient in enumerate(polynomial))
                    expected = _gn_logfraction(exact/(x**cache['quota']*cache['C']*(1 << (0, 1, 3)[row])))
                    assert abs(values[row]-expected) < GN_LD('1e-9')
                    assert abs(slopes[row]-GN_LD(exact_derivative.numerator)/exact_derivative.denominator/(GN_LD(exact.numerator)/exact.denominator)+cache['quota']) < GN_LD('1e-9')
                    counts['power_derivatives'] += 1
                    coefficient = polynomial[cache['quota']]/cache['C']/(1 << (0, 1, 3)[row])
                    assert coefficient == 0 or _gn_logfraction(coefficient) <= values[row]+GN_LD('1e-9')
                    counts['quota_upper_rows'] += 1
                endpoint, _ = _gn_evaluate(cache, _gn_logfraction(x), phi, endpoint=True)
                q = [(base+x)**width for base in cache['q_bases']]
                restricted = [[full[positions[i]][positions[j]] for j in range(3)] for i in range(3)]
                for row in range(3):
                    exact = sum(restricted[row][j]*q[j] for j in range(3))/(x**cache['quota']*cache['C'])
                    assert abs(endpoint[row]-_gn_logfraction(exact)) < GN_LD('1e-9')
                    counts['endpoint_rows'] += 1
            r0 = y**a
            alpha, beta = cache['gauge']
            assert (1+r0)**width*alpha**width == (y**(a-b)+r0)**width
            assert (1+r0)**width*(alpha*beta)**width == (2*r0)**width
            counts['gauge_endpoints'] += 2
    for a, c in ((1, 1), (4, 9), (9, 4)):
        def evaluate(t):
            x = np.exp2(t)
            value = a+3*x+c*x*x
            return np.log2(value)-t, (3*x+2*c*x*x)/value-1
        result = _gn_minimize(evaluate)
        expected = np.log2(3+2*np.sqrt(GN_LD(a*c)))
        assert result['lower'] <= expected+GN_LD('1e-9')
        assert result['upper'] >= expected-GN_LD('1e-9')
        counts['convex_minima'] += 1
    print(dict(**counts, exact_fraction_physical_reference=True,
               continuous_enclosures_FLOAT_only=True, full_code_certificate=False), flush=True)


def general_three_state_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--source', type=int, default=2)
    parser.add_argument('--temporal', type=int, default=1)
    parser.add_argument('--input', type=int, default=4)
    parser.add_argument('--y', type=F, default=F(3, 2))
    parser.add_argument('--threshold', type=int, default=184116)
    parser.add_argument('--half-updates', type=int, default=6)
    parser.add_argument('--log-p', type=float)
    parser.add_argument('--log-h', type=float)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(90, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            general_three_state_checks()
        else:
            print(general_three_state_control((args.source, args.temporal, args.input), args.y,
                                             args.threshold, args.half_updates, args.log_p,
                                             args.log_h, True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    general_three_state_main()
