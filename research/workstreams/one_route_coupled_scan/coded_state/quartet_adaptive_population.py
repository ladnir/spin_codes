"""Positive adaptive-cover action on one exponential population potential.

Reverse states are nested hole sets W subset A, with populations(a,w).
For phi(A,W)=alpha^a beta^w, form the homogeneous positive quartic
 P(W,E,Z)=sum_(w+e+z=4) Mult(4;w,e,z) h_(w+e,w)(x) W^w E^e Z^z,
 h_uw(x)=sum_j,r C(w,r)C(4-w,j-r)x^j y^(w-r)
                         *[p(u,w,j,r)+g(u,w,j,r)alpha^4 beta^u].
Then (K_x phi)(a,w) is the coefficient of W^w E^(a-w) Z^(n-a)
in P^(n/4), divided ONCE by Mult(n;w,a-w,n-a).

The engine uses15 positive shifts per quartet power, with three reusable
longdouble buffers and one common normalization across the valid triangle.
It does not use FFTs, coefficientwise maxima, or source-state aggregation.
Structural zeros outside the triangle stay zero. Positive terms and all
valid coefficients are guarded against underflow/overflow. Multinomial
division is performed in log space after the polynomial power.

Returned rows have a=0,4,...,n and0<=w<=a, in fixed lexicographic order.
Both log2(K_x phi) and log2(K_x phi/phi) are returned, with their a,w
arrays and potential logs. This is ONLY one exponential-mode action;
it is not an arbitrary-vector action, a regional quota operator, or an
iterated/global bound. All production arithmetic is noncertifying FLOAT.
No jobs or large arrays are created on import. Parent owns execution.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _ap_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_ap_thread_var] = '1'

import numpy as np


AP_LD = np.longdouble
AP_TERMS = tuple((w, e) for w in range(5) for e in range(5-w))


def _ap_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_adaptive_cover_checks as module
    return getattr(module, name)


def _ap_positive(array, stage):
    if np.any(array <= 0) or not np.all(np.isfinite(array)):
        raise ArithmeticError('Positive adaptive coefficient underflow/overflow at '+stage)


def _ap_logadd(terms):
    values = np.asarray(terms, dtype=AP_LD)
    if not len(values) or not np.all(np.isfinite(values)):
        raise ArithmeticError('Invalid positive local terms')
    maximum = values.max()
    weights = np.exp2(values-maximum)
    _ap_positive(weights, 'local positive sum')
    return maximum+np.log2(weights.sum(dtype=AP_LD))


def _ap_probability_logs(log_y, u, w, j, r):
    """Stable log2(g),log2(p), with -inf only for forced structural choices."""
    if j == 0:
        return AP_LD(0), AP_LD('-inf')
    b = 2*j-r
    a = max(0, 4+u-b)
    if a == 0:
        return AP_LD('-inf'), AP_LD(0)
    ln_y = log_y*np.log(AP_LD(2))
    denominator = -np.expm1(-(a+b)*ln_y)
    absent_numerator = -np.expm1(-a*ln_y)
    present_numerator = -np.expm1(-b*ln_y)
    _ap_positive(np.asarray((denominator, absent_numerator, present_numerator)), 'adaptive probabilities')
    return (-b*log_y+np.log2(absent_numerator)-np.log2(denominator),
            np.log2(present_numerator)-np.log2(denominator))


def adaptive_population_quartic(y=2, x=1, alpha=1, beta=1):
    """Return15 normalized positive coefficients and their common log2 scale."""
    y, x, alpha, beta = map(AP_LD, (y, x, alpha, beta))
    if not all(np.isfinite(value) for value in (y, x, alpha, beta)) or y <= 1 or min(x, alpha, beta) <= 0:
        raise ValueError('Require finite y>1 and positive x,alpha,beta')
    log_y, log_x, log_alpha, log_beta = map(np.log2, (y, x, alpha, beta))
    coefficients = []
    for w, e in AP_TERMS:
        u = w+e
        terms = []
        for j in range(5):
            for r in range(max(0, j-(4-w)), min(w, j)+1):
                log_g, log_p = _ap_probability_logs(log_y, u, w, j, r)
                common = (np.log2(AP_LD(comb(w, r)*comb(4-w, j-r)))
                          +j*log_x+(w-r)*log_y)
                if np.isfinite(log_p):
                    terms.append(common+log_p)
                if np.isfinite(log_g):
                    terms.append(common+log_g+4*log_alpha+u*log_beta)
        coefficients.append(np.log2(AP_LD(comb(4, u)*comb(u, w)))+_ap_logadd(terms))
    logs = np.asarray(coefficients, dtype=AP_LD)
    scale = logs.max()
    normalized = np.exp2(logs-scale)
    _ap_positive(normalized, 'common quartic scale')
    return dict(coefficients=normalized, scale=scale, terms=AP_TERMS,
                log2_coefficients=logs, y=y, x=x, alpha=alpha, beta=beta,
                log2_y=log_y, log2_x=log_x,
                log2_alpha=log_alpha, log2_beta=log_beta)


def adaptive_population_action(y=2, x=1, alpha=1, beta=1, width=512, progress=False, quartic=None):
    """Compute one exponential action, optionally from a validated positive quartic.

    The optional table changes coefficients only; the15-shift engine and old
    multinomial normalization are shared. A supplied table is responsible
    for its operator semantics and any model-specific unit-action identity.
    """
    width = index(width)
    if width < 4 or width > 512 or width % 4:
        raise ValueError('Use quartet width4..512')
    if width == 512 and np.finfo(AP_LD).maxexp < 16000:
        raise RuntimeError('Production action requires genuinely extended longdouble')
    started = time.monotonic()
    supplied_quartic = quartic is not None
    if quartic is None:
        quartic = adaptive_population_quartic(y, x, alpha, beta)
    else:
        quartic = dict(quartic)
        for key, value in zip(('y', 'x', 'alpha', 'beta'), map(AP_LD, (y, x, alpha, beta))):
            if not np.isfinite(value) or value <= (1 if key == 'y' else 0) or quartic[key] != value:
                raise ValueError('Supplied quartic marker metadata does not match the action')
            if not np.isfinite(quartic['log2_'+key]) or abs(quartic['log2_'+key]-np.log2(value)) > AP_LD('1e-12'):
                raise ValueError('Supplied quartic logarithmic metadata is inconsistent')
        quartic['coefficients'] = np.asarray(quartic['coefficients'], dtype=AP_LD)
        quartic['log2_coefficients'] = np.asarray(quartic['log2_coefficients'], dtype=AP_LD)
        if (tuple(quartic['terms']) != AP_TERMS or quartic['coefficients'].shape != (15,)
                or quartic['log2_coefficients'].shape != (15,) or not np.isfinite(quartic['scale'])):
            raise ValueError('Supplied quartic must have the fixed15 positive monomials')
        _ap_positive(quartic['coefficients'], 'supplied quartic')
        if (not np.all(np.isfinite(quartic['log2_coefficients']))
                or np.max(np.abs(np.log2(quartic['coefficients'])+quartic['scale']
                                 -quartic['log2_coefficients'])) > AP_LD('1e-10')):
            raise ValueError('Supplied quartic coefficient scale is inconsistent')
    coefficients, local_scale = quartic['coefficients'], quartic['scale']
    coefficient_minimum = coefficients.min()
    capacity = (width+1)**2
    current_storage = np.empty(capacity, dtype=AP_LD)
    after_storage = np.empty(capacity, dtype=AP_LD)
    product_storage = np.empty(capacity, dtype=AP_LD)
    current_storage[0] = 1
    current = current_storage[:1].reshape(1, 1)
    degree, scale, current_minimum = 0, AP_LD(0), AP_LD(1)
    axes = np.arange(width+1, dtype=np.int16)
    degree_sums = axes[:, None]+axes[None, :]
    for completed in range(1, width//4+1):
        size, next_size = degree+1, degree+5
        after = after_storage[:next_size*next_size].reshape(next_size, next_size)
        after.fill(0)
        product_buffer = product_storage[:size*size].reshape(size, size)
        if current_minimum*coefficient_minimum <= 0:
            raise ArithmeticError('A positive quartet-power product would underflow')
        # Each source square is contiguous; invalid-triangle entries are exact zeros.
        # Reuse one multiplication buffer instead of allocating15 temporaries per power.
        for (w, e), coefficient in zip(AP_TERMS, coefficients):
            np.multiply(current, coefficient, out=product_buffer)
            destination = after[w:w+size, e:e+size]
            np.add(destination, product_buffer, out=destination)
        degree += 4
        valid = degree_sums[:next_size, :next_size] <= degree
        positive = after[valid]
        _ap_positive(positive, 'complete valid coefficient triangle')
        if np.any(after[~valid] != 0):
            raise ArithmeticError('A structural coefficient zero became positive')
        maximum = positive.max()
        current_minimum = positive.min()/maximum
        if current_minimum <= 0 or not np.isfinite(current_minimum):
            raise ArithmeticError('Common polynomial normalization would underflow')
        np.divide(after, maximum, out=after)
        scale += local_scale+np.log2(maximum)
        current_storage, after_storage = after_storage, current_storage
        current = current_storage[:next_size*next_size].reshape(next_size, next_size)
        if progress and (completed % 32 == 0 or completed == width//4):
            print(dict(stage='adaptive_single_mode_power', quartets=completed,
                       degree=degree, seconds=round(time.monotonic()-started, 3)), flush=True)
    # Fixed physical order: a increasing in multiples4, then w=0,...,a.
    a_rows = np.repeat(np.arange(0, width+1, 4, dtype=np.int16), np.arange(0, width+1, 4)+1)
    w_rows = np.concatenate([np.arange(a+1, dtype=np.int16) for a in range(0, width+1, 4)])
    e_rows = a_rows-w_rows
    factorial_logs = np.zeros(width+1, dtype=AP_LD)
    factorial_logs[1:] = np.cumsum(np.log2(np.arange(1, width+1, dtype=AP_LD)), dtype=AP_LD)
    log_multinomial = (factorial_logs[width]-factorial_logs[w_rows]
                       -factorial_logs[e_rows]-factorial_logs[width-a_rows])
    log_action = scale+np.log2(current[w_rows, e_rows])-log_multinomial
    log_potential = a_rows*quartic['log2_alpha']+w_rows*quartic['log2_beta']
    if not np.all(np.isfinite(log_action)) or not np.all(np.isfinite(log_potential)):
        raise ArithmeticError('Invalid log-space population action')
    unit_error = None
    if not supplied_quartic and quartic['alpha'] == 1 and quartic['beta'] == 1:
        expected = (w_rows*np.log2(quartic['y']+quartic['x'])
                    +(width-w_rows)*np.log2(1+quartic['x']))
        unit_error = np.max(np.abs(log_action-expected))
        if unit_error > AP_LD('1e-8'):
            raise ArithmeticError('The all-row unit-potential identity failed')
    return dict(a=a_rows, w=w_rows, log2_action=log_action,
                log2_ratio=log_action-log_potential, log2_potential=log_potential,
                log2_alpha=quartic['log2_alpha'], log2_beta=quartic['log2_beta'],
                log2_x=quartic['log2_x'], log2_y=quartic['log2_y'],
                y=quartic['y'], x=quartic['x'], alpha=quartic['alpha'], beta=quartic['beta'],
                width=width, row_order='a=0,4,...,n; w=0,...,a',
                polynomial_log2_scale=scale, quartic_log2_coefficients=quartic['log2_coefficients'],
                common_positive_scale=True, all_physical_rows_retained=True,
                single_exponential_mode_only=True, arbitrary_vector_action=False,
                regional_or_global_bound=False, floating_diagnostic_only=True,
                unit_potential_all_row_log2_error=unit_error,
                supplied_positive_quartic=supplied_quartic,
                seconds=time.monotonic()-started)


def _ap_exact_quartic(y, x, alpha, beta):
    local_coefficients = _ap_dependency('adaptive_cover_local_coefficients')
    result = {}
    for w, e in AP_TERMS:
        u = w+e
        normalized = local_coefficients(y, alpha, beta, u, w)
        value = sum(coefficient*x**j for j, coefficient in enumerate(normalized))*alpha**u*beta**w
        result[(w, e)] = comb(4, u)*comb(u, w)*value
    return result


def _ap_exact_power(quartic, groups):
    result = {(0, 0): Fraction(1)}
    for _ in range(groups):
        after = {}
        for (w, e), coefficient in result.items():
            for (dw, de), local in quartic.items():
                key = (w+dw, e+de)
                after[key] = after.get(key, Fraction(0))+coefficient*local
        result = after
    return result


def _ap_partitions(width):
    if width == 4:
        return (((0, 1, 2, 3),),)
    assert width == 8
    result = []
    for rest in combinations(range(1, 8), 3):
        first = (0,)+rest
        second = tuple(i for i in range(8) if i not in first)
        result.append((first, second))
    assert len(result) == 35
    return tuple(result)


@lru_cache(maxsize=None)
def _ap_physical_quartet(y, x, alpha, beta, A, W):
    """Explicit labeled J masks and both adaptive choices, without coefficients."""
    assert not (W & ~A)
    weights = _ap_dependency('adaptive_cover_weights')
    u, w, value = A.bit_count(), W.bit_count(), Fraction(0)
    for J in range(16):
        j, r = J.bit_count(), (J & W).bit_count()
        hole, present, _, _, _ = weights(y, u, w, j, r)
        for V, probability in ((0, present), (15, hole)):
            next_A, next_W = V, V & A
            value += (probability*x**j*y**((15 ^ J) & W).bit_count()
                      *alpha**next_A.bit_count()*beta**next_W.bit_count())
    return value


def _ap_physical_action(width, y, x, alpha, beta, a, w):
    A, W = (1 << a)-1, (1 << w)-1
    value = Fraction(0)
    partitions = _ap_partitions(width)
    for partition in partitions:
        term = Fraction(1)
        for quartet in partition:
            local_A = sum(((A >> coordinate) & 1) << slot for slot, coordinate in enumerate(quartet))
            local_W = sum(((W >> coordinate) & 1) << slot for slot, coordinate in enumerate(quartet))
            term *= _ap_physical_quartet(y, x, alpha, beta, local_A, local_W)
        value += term
    return value/len(partitions)


def adaptive_population_checks():
    local_checks = coefficient_checks = physical_checks = action_checks = unit_checks = 0
    for width in (4, 8):
        for y in (Fraction(2), Fraction(4)):
            for x, alpha, beta in ((Fraction(1), Fraction(1), Fraction(1)),
                                   (Fraction(2), Fraction(3, 2), Fraction(2, 3)),
                                   (Fraction(1, 2), Fraction(2, 3), Fraction(3, 2))):
                quartic = _ap_exact_quartic(y, x, alpha, beta)
                power = _ap_exact_power(quartic, width//4)
                floating = adaptive_population_action(float(y), float(x), float(alpha), float(beta), width)
                for term, found in zip(AP_TERMS, floating['quartic_log2_coefficients']):
                    exact = quartic[term]
                    expected = AP_LD(exact.numerator)/exact.denominator
                    assert abs(np.exp2(found)/expected-1) < AP_LD('1e-12')
                    local_checks += 1
                assert len(power) == (width+1)*(width+2)//2 and all(value > 0 for value in power.values())
                coefficient_checks += len(power)
                for position, (a_value, w_value) in enumerate(zip(floating['a'], floating['w'])):
                    a, w = int(a_value), int(w_value)
                    exact = power[(w, a-w)]/(comb(width, a)*comb(a, w))
                    physical = _ap_physical_action(width, y, x, alpha, beta, a, w)
                    assert exact == physical
                    physical_checks += 1
                    expected = AP_LD(exact.numerator)/exact.denominator
                    assert abs(np.exp2(floating['log2_action'][position])/expected-1) < AP_LD('1e-11')
                    potential = alpha**a*beta**w
                    ratio = exact/potential
                    expected_ratio = AP_LD(ratio.numerator)/ratio.denominator
                    assert abs(np.exp2(floating['log2_ratio'][position])/expected_ratio-1) < AP_LD('1e-11')
                    action_checks += 2
                    if alpha == beta == 1:
                        assert exact == (y+x)**w*(1+x)**(width-w)
                        unit_checks += 1
                assert np.array_equal(floating['a'], np.asarray([a for a in range(0, width+1, 4) for _ in range(a+1)]))
                assert np.array_equal(floating['w'], np.asarray([w for a in range(0, width+1, 4) for w in range(a+1)]))
    print(dict(exact_local_quartic_comparisons=local_checks,
               positive_exact_power_coefficients=coefficient_checks,
               physical_partition_action_equalities=physical_checks,
               FLOAT_action_and_ratio_comparisons=action_checks,
               exact_unit_potential_identities=unit_checks,
               exact_width8_partition_count=35, single_exponential_mode_only=True,
               full_code_certificate=False), flush=True)


def adaptive_population_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--y', type=float, default=2.)
    parser.add_argument('--x', type=float, default=1.)
    parser.add_argument('--alpha', type=float, default=1.)
    parser.add_argument('--beta', type=float, default=1.)
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
            adaptive_population_checks()
        else:
            result = adaptive_population_action(args.y, args.x, args.alpha, args.beta, progress=True)
            maximum = int(np.argmax(result['log2_ratio']))
            print(dict(y=args.y, x=args.x, alpha=args.alpha, beta=args.beta,
                       physical_rows=len(result['a']),
                       minimum_log2_action=float(result['log2_action'].min()),
                       maximum_log2_action=float(result['log2_action'].max()),
                       maximum_log2_ratio=float(result['log2_ratio'][maximum]),
                       maximum_ratio_row=(int(result['a'][maximum]), int(result['w'][maximum])),
                       unit_potential_all_row_log2_error=(None if result['unit_potential_all_row_log2_error'] is None
                                                        else float(result['unit_potential_all_row_log2_error'])),
                       seconds=result['seconds'],
                       single_exponential_mode_only=True, regional_or_global_bound=False,
                       floating_diagnostic_only=True), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_population_main()
