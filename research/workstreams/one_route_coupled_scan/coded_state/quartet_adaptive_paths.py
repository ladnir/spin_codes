r"""Separated-defect diagonal paths of the actual adaptive population operator.

Let n=4m and fix k with9k<=m. At old population(a,w)=(n-4k,n-8k),
the E=A\W and Z=A^c defects each have4k positions. A fresh partition
separates all8k defects with probability
 Psep=C(m,8k)*4^(8k)/C(n,8k).
There are4k quartets of type(u,w)=(4,3),4k of type(3,3), andm-8k
pure(4,4) quartets. Keep every mixed quartet a hole and select exactlyk
pure quartets. The next population equals the old population.

For T_uw(x)=sum_J x^|J| y^|W\J| g and S_uw with p instead of g,
the positive one-step diagonal contribution is
 J_k=Psep*C(m-8k,k)*T44^(m-9k)*T43^(4k)*T33^(4k)*S44^k.
Thus[x^h]J_k^R/C(nR,h) lower-bounds the exact regional diagonal.
The probe evaluates this coefficient by positive scalar convolutions,
factoring J_k^R into quartic factors. It does not take a maximum over
old types, iterate a potential, or assemble a population matrix.

Only a diagonal contribution is reported. No initial/terminal path is
included, so this is NOT a whole-moment or failure-probability floor.
Production arithmetic is FLOAT and not a directed lower certificate.
Exact small checks enumerate partitions, labels, and selected quartets.
No jobs or large allocations occur on import; parent owns execution.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _adp_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_adp_thread_var] = '1'

import numpy as np

ADP_LD = np.longdouble


def _adp_dependency(name):
    if name in globals():
        return globals()[name]
    if name == 'adaptive_cover_weights':
        import quartet_adaptive_cover_checks as module
    elif name == '_cf_coefficient_log2':
        import quartet_cover_holes_floor as module
    else:
        import quartet_adaptive_population as module
    return getattr(module, name)


def adaptive_path_local(y, u, w, selected=False):
    """Exact local polynomial for one fixed old type; no multinomial factor."""
    y = Fraction(y)
    if y <= 1 or not 0 <= w <= u <= 4:
        raise ValueError('Invalid marker or old quartet type')
    weights = _adp_dependency('adaptive_cover_weights')
    coefficients = []
    for j in range(5):
        value = Fraction(0)
        for r in range(max(0, j-(4-w)), min(w, j)+1):
            g, p, _, _, _ = weights(y, u, w, j, r)
            value += comb(w, r)*comb(4-w, j-r)*y**(w-r)*(p if selected else g)
        coefficients.append(value)
    return tuple(coefficients)


def _adp_multiply(left, right):
    result = [Fraction(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                result[i+j] += a*b
    return result


def _adp_power(poly, exponent):
    result = [Fraction(1)]
    for _ in range(exponent):
        result = _adp_multiply(result, poly)
    return result


def adaptive_path_general_polynomial(width, e, z, k, y=2):
    """Exact general separated-defect contribution, also for arbitrary toy colors.

    Old(a,w)=(n-z,n-z-e); target(a',w')=(n-4k,n-z-4k).
    Toy old a need not be a reachable full-history multiple of four.
    Production uses e=z=4k, where both populations coincide.
    """
    width, e, z, k = map(index, (width, e, z, k))
    groups, defects = width//4, e+z
    if width < 4 or width % 4 or min(e, z, k) < 0 or defects > groups or k > groups-defects:
        raise ValueError('Invalid separated-defect geometry')
    separation = Fraction(comb(groups, defects)*4**defects, comb(width, defects))
    polynomial = [separation*comb(groups-defects, k)]
    for old_type, selected, count in (((4, 4), False, groups-defects-k),
                                      ((4, 3), False, e), ((3, 3), False, z),
                                      ((4, 4), True, k)):
        polynomial = _adp_multiply(polynomial, _adp_power(adaptive_path_local(y, *old_type, selected), count))
    return tuple(polynomial)


def _adp_factor(poly, tilt):
    shift = next((i for i, value in enumerate(poly) if value), None)
    if shift is None or any(value <= 0 for value in poly[shift:]):
        raise ValueError('Production factors must have positive contiguous support')
    values = np.asarray([ADP_LD(value.numerator)/value.denominator*tilt**j
                         for j, value in enumerate(poly) if j >= shift], dtype=ADP_LD)
    if np.any(values <= 0) or not np.all(np.isfinite(values)):
        raise ArithmeticError('Tilted adaptive path factor underflow/overflow')
    return values, ADP_LD(0), shift


def adaptive_path_diagonal(k, width=512, region_updates=16, quota=4096, y=2, tilt=2):
    """One exact-algebra regional diagonal lower contribution, evaluated in FLOAT."""
    k, width, steps, quota = map(index, (k, width, region_updates, quota))
    groups = width//4
    y, tilt = Fraction(y), ADP_LD(tilt)
    if width < 4 or width > 512 or width % 4 or k < 1 or 9*k > groups:
        raise ValueError('Require1<=k<=floor((width/4)/9)')
    if steps < 1 or not 0 <= quota <= width*steps or y <= 1 or tilt <= 0 or not np.isfinite(tilt):
        raise ValueError('Invalid marker, common tilt, regional length, or quota')
    if width == 512 and np.finfo(ADP_LD).maxexp < 16000:
        raise RuntimeError('Production path control requires extended longdouble')
    counts = (steps*(groups-9*k), steps*4*k, steps*4*k, steps*k)
    polys = (adaptive_path_local(y, 4, 4), adaptive_path_local(y, 4, 3),
             adaptive_path_local(y, 3, 3), adaptive_path_local(y, 4, 4, True))
    factors = []
    for poly, count in zip(polys, counts):
        factors.extend([_adp_factor(poly, tilt)]*count)
    started = time.monotonic()
    coefficient_log = _adp_dependency('_cf_coefficient_log2')(factors, quota, tilt)
    separation = Fraction(comb(groups, 8*k)*4**(8*k), comb(width, 8*k))
    log_separation = np.log2(ADP_LD(separation.numerator))-np.log2(ADP_LD(separation.denominator))
    log_choices = np.log2(ADP_LD(comb(groups-8*k, k)))
    coefficient_log += steps*(log_separation+log_choices)
    log_quota = np.log2(ADP_LD(comb(width*steps, quota)))
    return dict(k=k, width=width, region_updates=steps, quota=quota, y=float(y), tilt=float(tilt),
                old_and_new_population=(width-4*k, width-8*k),
                log2_separation_probability=float(log_separation),
                log2_pure_selection_multiplicity=float(log_choices),
                log2_regional_path_coefficient=float(coefficient_log),
                log2_quota_binomial=float(log_quota),
                log2_conditioned_diagonal_lower_control=float(coefficient_log-log_quota),
                quartic_factor_counts=counts, positive_scalar_convolutions=True,
                actual_adaptive_operator_diagonal=True, initial_and_terminal_paths_included=False,
                whole_moment_lower_bound=False, actual_failure_probability_lower_bound=False,
                directed_lower_certificate=False, floating_diagnostic_only=True,
                seconds=time.monotonic()-started)


def _adp_physical_polynomial(width, e, z, k, y):
    """Enumerate physical partitions, all J masks, and retained pure selections."""
    full = (1 << width)-1
    Z = (1 << z)-1
    E = ((1 << e)-1) << z
    A, W = full ^ Z, full ^ (Z | E)
    partitions = _adp_dependency('_ap_partitions')(width)
    weights = _adp_dependency('adaptive_cover_weights')
    result = [Fraction(0)]*(width+1)
    separated = branches = 0
    for partition in partitions:
        masks = tuple(sum(1 << coordinate for coordinate in quartet) for quartet in partition)
        if any(((E | Z) & mask).bit_count() > 1 for mask in masks):
            continue
        separated += 1
        pure = [mask for mask in masks if not ((E | Z) & mask)]
        for selected in combinations(pure, k):
            U = 0
            for mask in selected:
                U |= mask
            V = full ^ U
            assert (V.bit_count(), (V & A).bit_count()) == (width-4*k, width-z-4*k)
            branches += 1
            for J in range(1 << width):
                weight = Fraction(1)
                for mask in masks:
                    u, w = (A & mask).bit_count(), (W & mask).bit_count()
                    j, r = (J & mask).bit_count(), (J & W & mask).bit_count()
                    g, p, _, _, _ = weights(y, u, w, j, r)
                    weight *= y**(w-r)*(p if U & mask else g)
                result[J.bit_count()] += weight/len(partitions)
    assert Fraction(separated, len(partitions)) == Fraction(comb(width//4, e+z)*4**(e+z), comb(width, e+z))
    assert branches == separated*comb(width//4-e-z, k)
    return tuple(result), branches


def adaptive_path_checks():
    local_checks = physical_checks = branch_checks = scalar_checks = diagonal_checks = 0
    weights = _adp_dependency('adaptive_cover_weights')
    for y in (Fraction(2), Fraction(4)):
        for u in range(5):
            for w in range(u+1):
                A, W = (1 << u)-1, (1 << w)-1
                for selected in (False, True):
                    physical = [Fraction(0)]*5
                    for J in range(16):
                        g, p, _, _, _ = weights(y, u, w, J.bit_count(), (J & W).bit_count())
                        physical[J.bit_count()] += y**((15 ^ J) & W).bit_count()*(p if selected else g)
                    assert tuple(physical) == adaptive_path_local(y, u, w, selected)
                    local_checks += 5
        for width in (4, 8):
            groups = width//4
            for e in range(groups+1):
                for z in range(groups-e+1):
                    for k in range(groups-e-z+1):
                        exact, branches = _adp_physical_polynomial(width, e, z, k, y)
                        formula = adaptive_path_general_polynomial(width, e, z, k, y)
                        assert exact == formula
                        physical_checks += len(formula)
                        branch_checks += branches
    # A genuine diagonal needsm>=9k. These bounded cases test the quota
    # factorization independently of the general physical one-step toys above.
    for width in (36, 40):
        k, steps, quota, y = 1, 2, width, Fraction(2)
        polynomial = adaptive_path_general_polynomial(width, 4, 4, k, y)
        exact = _adp_power(polynomial, steps)[quota]/comb(width*steps, quota)
        for tilt in (1, 2, 3):
            result = adaptive_path_diagonal(k, width, steps, quota, y, tilt)
            reference = ADP_LD(exact.numerator)/exact.denominator
            assert abs(np.exp2(ADP_LD(result['log2_conditioned_diagonal_lower_control']))/reference-1) < ADP_LD('1e-10')
            scalar_checks += 1
        # Construct a separated colored configuration and enumerate its pure choices.
        groups = width//4
        E = sum(1 << (4*i) for i in range(4))
        Z = sum(1 << (4*i) for i in range(4, 8))
        full, A = (1 << width)-1, ((1 << width)-1) ^ Z
        W = A ^ E
        masks = [15 << (4*i) for i in range(groups)]
        types = [((A & q).bit_count(), (W & q).bit_count()) for q in masks]
        assert types.count((4, 3)) == 4 and types.count((3, 3)) == 4 and types.count((4, 4)) == groups-8
        for selected in masks[8:]:
            V = full ^ selected
            assert (V.bit_count(), (V & A).bit_count()) == (A.bit_count(), W.bit_count())
            diagonal_checks += 1
    print(dict(exact_local_coefficients=local_checks, physical_partition_polynomial_coefficients=physical_checks,
               retained_physical_selector_branches=branch_checks, exact_quota_factorization_controls=scalar_checks,
               explicit_separated_diagonal_returns=diagonal_checks, exact_width8_partition_count=35,
               toy_old_colors_need_not_be_reachable_history_states=True,
               whole_moment_lower_bound=False), flush=True)


def adaptive_path_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--k-values', nargs='+', type=int, default=(4, 8, 10, 12))
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
            adaptive_path_checks()
        else:
            for k in args.k_values:
                print(adaptive_path_diagonal(k, tilt=args.tilt), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_path_main()
