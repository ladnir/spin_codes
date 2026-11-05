"""Bounded field correction for an independently chosen ACTUAL death marker.

The supplied moment must use kernel marker d^D, where log_death_marker
is log2(d), independently of x and theta.  It must already include the
regional H-quota conditioning, but not the global C(8192,h) union or
fixed-profile marker/field corrections.  It may be an exact-birth
coefficient moment; then no nu^-f correction is applied.

For t_x=qx^2/(1+qx), B_x=(gamma/q)q^-R x^-W(1+qx)^Ls, the scalar
field correction is sup_D [theta*log2(B_x*t_x^D)-D*log2(d)].
With the optional probability cap it becomes
sup_D [theta*min(0,log2(B_x)+D*log2(t_x))-D*log2(d)].
No D union is incurred inside one interval.  Bounds for disjoint D
intervals are SUMMED, never minimized.  Actual census complements are
used only with their separate failure budgets; no availability cuts enter.

This helper evaluates an already supplied positive-moment bound.  It
does not build a second expensive kernel or optimize marker grids.  All
production evaluations are FLOAT diagnostics, not distance certificates.
Parent runners may inject the actual-state helper; normal local import
is a fallback.  No extra dependencies beyond the existing NumPy stack.
"""

import argparse
from fractions import Fraction
from itertools import product
from math import comb, lgamma, log
from operator import index
import os
import shutil
import threading

import numpy as np


AD_LD = np.longdouble
AD_Q = 1 << 32
AD_W = 943718
AD_PROFILE = (4096, 786432, 786432, 4, 1048576)
AD_INTERVALS = ((0, 157286), (157287, 190053), (190054, 262144))


def _ad_ld(value):
    return AD_LD(value.numerator) / value.denominator if isinstance(value, Fraction) else AD_LD(value)


def _ad_census():
    if 'actual_state_census_interval' in globals():
        return globals()['actual_state_census_interval']
    from quartet_actual_state_profiles import actual_state_census_interval
    return actual_state_census_interval


def actual_death_scalar_correction(log_B, log_true_t, theta, log_death_marker,
                                    interval, cap_probability=True):
    """Uniform scalar field/marker correction; only endpoints and one kink matter."""
    lo, hi = map(index, interval)
    log_B, log_true_t, theta, log_death_marker = map(_ad_ld,
        (log_B, log_true_t, theta, log_death_marker))
    if lo < 0 or lo > hi or not 0 <= theta <= 1 or not all(np.isfinite(v)
            for v in (log_B, log_true_t, theta, log_death_marker)):
        raise ValueError('Invalid field marker or death interval')
    candidates = {lo, hi}
    if cap_probability and log_true_t != 0:
        crossing = -log_B / log_true_t
        if lo - 1 <= crossing <= hi + 1:
            left = int(np.floor(crossing))
            candidates.update(v for v in (left, left + 1) if lo <= v <= hi)
    best = None
    maximizing = None
    for D in sorted(candidates):
        field = log_B + D * log_true_t
        if cap_probability:
            field = min(AD_LD(0), field)
        value = theta * field - D * log_death_marker
        if best is None or value > best:
            best, maximizing = value, D
    return dict(log2_field_marker_correction=float(best), maximizing_D=maximizing,
                residual_coefficient=float(theta * log_true_t - log_death_marker),
                theta_log_B=float(theta * log_B), probability_cap_used=cap_probability,
                candidate_death_counts=sorted(candidates), no_death_count_union=True)


def actual_death_profile_bound(log_moment, profile=AD_PROFILE,
                               intervals=AD_INTERVALS, log_death_marker=-1,
                               exponent=16, theta=1/16, log_s=1, log_z=0,
                               log_nu=-18, log_eta=-2, birth_extracted=False,
                               cap_probability=True):
    """One bounded fixed-profile evaluation from a quota-conditioned moment.

    Set birth_extracted=True for a supplied exact-f coefficient, and set
    log_nu=0.  Otherwise the moment must already contain nu^factual.
    The default three intervals cover the audited hard-profile D range.
    A custom partial list reports covers_refined_interval=False.
    """
    h, R, Ls, factual, F = map(index, profile)
    census = _ad_census()(h, R, Ls, factual, F)
    result = dict(profile=(h, R, Ls, factual, F), statistics='actual_state_support',
                  scope='profile_remainder_outside_actual_prior_classes',
                  census_status=census['status'], floating_diagnostic_only=True,
                  full_code_certificate=False)
    if census['status'] != 'remaining':
        return dict(**result, slices=[], log2_sum_bound=None, covers_refined_interval=True)
    lower, upper = census['refined_death_interval']
    requested = sorted((index(lo), index(hi)) for lo, hi in intervals)
    if any(lo < 0 or lo > hi for lo, hi in requested) or any(
            requested[j][0] <= requested[j - 1][1] for j in range(1, len(requested))):
        raise ValueError('Death slices must be nonnegative, nonempty, and disjoint')
    retained = [(max(lo, lower), min(hi, upper)) for lo, hi in requested if hi >= lower and lo <= upper]
    values = (log_moment, log_death_marker, exponent, theta, log_s, log_z, log_nu, log_eta)
    log_moment, log_death_marker, exponent, theta, log_s, log_z, log_nu, log_eta = map(AD_LD, values)
    if not all(np.isfinite(v) for v in (log_moment, log_death_marker, exponent, theta, log_s, log_z, log_nu, log_eta)):
        raise ValueError('Moment and markers must be finite')
    if exponent < 0 or not 0 <= theta <= 1 or (birth_extracted and log_nu != 0):
        raise ValueError('Require exponent>=0, theta in [0,1], and log_nu=0 after exact birth extraction')
    log_lambda = np.logaddexp2(AD_LD(0), 32 - exponent)
    log_true_t = 32 - 2 * exponent - log_lambda
    log_B = -31 - 32 * R + AD_W * exponent + Ls * log_lambda
    eligible_union = AD_LD((lgamma(8193) - lgamma(h + 1) - lgamma(8193 - h)) / log(2.0))
    marker = (-log_s * Ls - log_z * (F - factual) - log_eta * (R - 128 * h)
              - (AD_LD(0) if birth_extracted else log_nu * factual))
    common = log_moment + eligible_union + marker
    slices = []
    for interval in retained:
        field = actual_death_scalar_correction(log_B, log_true_t, theta,
                                               log_death_marker, interval, cap_probability)
        slices.append(dict(interval=interval, log2_bound=float(common + field['log2_field_marker_correction']), **field))
    if slices:
        entries = np.array([item['log2_bound'] for item in slices], dtype=AD_LD)
        maximum = entries.max()
        total = float(maximum + np.log2(np.exp2(entries - maximum).sum()))
    else:
        total = float('-inf')
    covers = bool(retained and retained[0][0] == lower and retained[-1][1] == upper
                  and all(retained[j][0] == retained[j - 1][1] + 1 for j in range(1, len(retained))))
    return dict(**result, refined_death_interval=(lower, upper),
                census_improves_physical_interval=(lower, upper) != (census['dmin'], census['dmax']),
                log_moment=float(log_moment), log_death_marker=float(log_death_marker),
                log_true_t=float(log_true_t), scalar_log_B=float(log_B), theta=float(theta),
                eligible_union=float(eligible_union), marker_correction=float(marker),
                birth_coefficient_extracted=birth_extracted, separate_routing_factor_used=False,
                slices=slices, log2_sum_bound=total, covers_refined_interval=covers,
                combination='sum_of_disjoint_slice_bounds', global_profile_union_not_yet_applied=True)


def actual_death_checks():
    F = Fraction
    scalar_checks = coefficient_checks = weighted_checks = 0
    for B, t, d, theta in product((-4, 0, 5), (-3, -1, 0), (-4, -1, 1), (F(0), F(1, 2), F(1))):
        for cap in (False, True):
            actual = actual_death_scalar_correction(B, t, theta, d, (1, 12), cap)
            reference = max(theta * (min(0, B + t * D) if cap else B + t * D) - d * D for D in range(1, 13))
            assert abs(AD_LD(actual['log2_field_marker_correction']) - AD_LD(reference.numerator) / reference.denominator) < AD_LD('1e-12')
            scalar_checks += 1
    # Exact projective coefficient rewrite and terminal supremum, for toy q.
    # This verifies the alternative identity; no broad coefficient optimizer runs.
    for q in (4, 16):
        for Ls in range(1, 13):
            for D in range(Ls + 1):
                for terminal in range(min(4, Ls - D) + 1):
                    births, survivors = D + terminal, Ls - D - terminal
                    for W in range(2 * D + terminal, 2 * Ls + 2):
                        limit = min(survivors, W - 2 * D - terminal)
                        delta = max(0, Ls + D - W)
                        polynomial = sum(comb(survivors, k) * q ** (births + k) for k in range(limit + 1))
                        leading = q ** (Ls - delta) * comb(survivors, delta)
                        assert q ** (births + limit) * comb(survivors, limit) == leading
                        if survivors < q:
                            assert polynomial <= leading / (1 - F(survivors, q))
                        assert comb(survivors, delta) <= comb(Ls - D, delta)
                        coefficient_checks += 1
    # Nonuniform exact toy weights: uniform correction times the WHOLE
    # free-marker moment bounds each slice; slice bounds must be summed.
    counts = {D: 1 + D * D for D in range(5)}
    for d in (F(1, 8), F(1, 2), F(2)):
        moment = sum(weight * d ** D for D, weight in counts.items())
        field = {D: min(F(1), F(3, 2) * F(1, 4) ** D) for D in counts}
        summed = F(0)
        for lo, hi in ((0, 1), (2, 2), (3, 4)):
            correction = max(field[D] / d ** D for D in range(lo, hi + 1))
            actual = sum(counts[D] * field[D] for D in range(lo, hi + 1))
            assert actual <= moment * correction
            summed += moment * correction
            weighted_checks += 1
        assert sum(counts[D] * field[D] for D in counts) <= summed
    census = _ad_census()(*AD_PROFILE)
    assert census['status'] == 'remaining' and census['refined_death_interval'] == (0, 262144)
    assert AD_W - AD_PROFILE[2] == 157286
    assert AD_INTERVALS[0][0] == 0 and AD_INTERVALS[-1][1] == 262144
    assert all(AD_INTERVALS[j][0] == AD_INTERVALS[j - 1][1] + 1 for j in (1, 2))
    print(dict(exact_scalar_supremum_checks=scalar_checks, exact_projective_coefficient_checks=coefficient_checks,
               exact_weighted_slice_checks=weighted_checks, hard_profile_census_interval=(0, 262144),
               census_improvement_claimed=False, kernel_runs=0, full_code_certificate=False), flush=True)


def actual_death_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--bound', action='store_true')
    parser.add_argument('--moment-log2', type=float)
    parser.add_argument('--interval', nargs=2, action='append', type=int)
    parser.add_argument('--log-death-marker', type=float, default=-1)
    parser.add_argument('--exponent', type=float, default=16)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--log-s', type=float, default=1)
    parser.add_argument('--log-z', type=float, default=0)
    parser.add_argument('--log-nu', type=float, default=-18)
    parser.add_argument('--log-eta', type=float, default=-2)
    parser.add_argument('--birth-extracted', action='store_true')
    parser.add_argument('--uncapped-field', action='store_true')
    args = parser.parse_args()
    if args.bound and args.moment_log2 is None:
        parser.error('--bound requires a supplied --moment-log2; this helper runs no kernel')
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            actual_death_checks()
        else:
            print(actual_death_profile_bound(args.moment_log2, intervals=args.interval or AD_INTERVALS,
                log_death_marker=args.log_death_marker, exponent=args.exponent, theta=args.theta,
                log_s=args.log_s, log_z=args.log_z, log_nu=args.log_nu, log_eta=args.log_eta,
                birth_extracted=args.birth_extracted, cap_probability=not args.uncapped_field), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    actual_death_main()
