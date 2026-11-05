"""Geometric profile-box adapter for the exact MDS/full-E comparison.

Inject the MDS, joint-subset, joint-Chernoff and rank-box definitions with
their main blocks removed.  This adapter adds no field/output prefactor:
the moment is E sum_V a^-R b^L u^f chi^F, with alpha=1.  The existing
box_diagnostic supplies the independent optimized conditional field bound.
Canonical feasibility is screened before any numerical moment is formed.
Passing that screen is not a proof that a footprint or codeword exists.
All large results are FLOAT diagnostics.  Run numerical jobs serially.
"""

import argparse
from fractions import Fraction
from math import comb, lgamma
import os
import shutil
import threading
import time

import numpy as np


JBOX_LD = np.longdouble
JBOX_N = 1 << 21
JBOX_K = 1 << 20
JBOX_W = 943718
JBOX_OBSERVED = (903775, 891035, 71426, 1132763)
JBOX_CANDIDATE = (900000, 900000, 20000, 1150000)
JBOX_LEGACY = ((900000, 1000000, 40000, 1200000),
               (900000, 1100000, 20000, 1200000),
               (900000, 900000, 20000, 1200000))


def _joint_boxes_dependencies():
    required = ('joint_mds_local', 'joint_chernoff_prepare', 'joint_chernoff_comparison',
                'joint_subset_region', 'joint_subset_tail', 'box_diagnostic')
    if any(name not in globals() for name in required):
        raise RuntimeError('Inject MDS, subset, Chernoff and rank-box helper definitions first')


def joint_boxes_canonical_screen(h, box, bad_weight=JBOX_W):
    """Integer interval propagation of necessary canonical inequalities.

    A nonempty result means only 'not excluded by these aggregate tests'.
    The returned bounds can safely replace the supplied box in a bad-word
    count.  No feasibility assumption is introduced into the kernel.
    """
    if not isinstance(h, int) or not 1 <= h <= 8192 or len(box) != 4:
        raise ValueError('Require an integer eligible-row count and four profile intervals')
    bounds = [list(map(int, interval)) for interval in box]
    if any(tuple(interval) != tuple(converted) for interval, converted in zip(box, bounds)):
        raise ValueError('Profile bounds must be integers')
    constraints = [
        ('R >= 128h', (-1, 0, 0, 0), -128 * h),
        ('R <= K', (1, 0, 0, 0), JBOX_K),
        ('F-R >= h', (1, 0, 0, -1), -h),
        ('F-R <= 128h', (-1, 0, 0, 1), 128 * h),
        ('f <= L', (0, -1, 1, 0), 0),
        ('L <= F', (0, 1, 0, -1), 0),
        ('F <= 4L+f', (0, -4, -1, 1), 0),
        ('f <= F-L+512', (0, 1, 1, -1), 512),
        ('2F-L-f <= N-512', (0, -1, -1, 2), JBOX_N - 512),
        ('F <= 2R', (-2, 0, 0, 1), 0),
        ('f >= 1', (0, 0, -1, 0), -1),
        ('F <= N', (0, 0, 0, 1), JBOX_N),
    ]
    if bad_weight is not None:
        constraints.append(('2f-512 <= W', (0, 0, 2, 0), int(bad_weight) + 512))
    changed = True
    while changed:
        changed = False
        for name, coefficients, rhs in constraints:
            if any(lo > hi for lo, hi in bounds):
                return dict(status='empty', reason=name, bounds=tuple(map(tuple, bounds)))
            minimum = sum(coefficient * bounds[j][0 if coefficient >= 0 else 1]
                          for j, coefficient in enumerate(coefficients))
            if minimum > rhs:
                return dict(status='empty', reason=name, bounds=tuple(map(tuple, bounds)))
            for j, coefficient in enumerate(coefficients):
                if not coefficient:
                    continue
                own_minimum = coefficient * bounds[j][0 if coefficient >= 0 else 1]
                available = rhs - (minimum - own_minimum)
                if coefficient > 0:
                    upper = available // coefficient
                    if upper < bounds[j][1]:
                        bounds[j][1] = upper
                        changed = True
                else:
                    lower = -((-available) // coefficient)  # ceil(available/coefficient)
                    if lower > bounds[j][0]:
                        bounds[j][0] = lower
                        changed = True
    return dict(status='not_excluded', reason='necessary aggregate conditions only',
                bounds=tuple(map(tuple, bounds)))


def joint_boxes_minimum(markers, box):
    loga, logb, logu, logchi = markers
    return sum(coefficient * interval[0 if coefficient >= 0 else 1]
               for coefficient, interval in zip((-loga, logb, logu, logchi), box))


def joint_boxes_geometry_logmoment(h, loga, logb, logu, logchi=0, tilts=3,
                                   u_step=32, e_step=32, logpotential=0,
                                   log_hole_potential=0, tilt_logs=(-1/16, 0, 1/16),
                                   progress=False):
    """Only the geometric moment; arbitrary finite real profile markers."""
    _joint_boxes_dependencies()
    if not 0 < h < 8192 or tilts < 1:
        raise ValueError('Require 0<h<8192 and a positive quota-tilt count')
    markers = tuple(map(JBOX_LD, (loga, logb, logu, logchi)))
    if not all(np.isfinite(value) for value in markers):
        raise ValueError('Geometric log2 markers must be finite')
    loga, logb, logu, logchi = markers
    a, b, birth_ratio, chi = map(np.exp2, markers)
    m = b * birth_ratio
    if not all(np.isfinite(value) and value > 0 for value in (a, b, m, chi)):
        raise ValueError('Geometric weights overflowed or underflowed')
    actual = JBOX_LD(h) / 8192
    low = actual / (actual + a * (1 - actual))
    choose = JBOX_LD((lgamma(8193) - lgamma(h + 1) - lgamma(8193 - h))
                    / float(np.log(JBOX_LD(2))))
    prepared = joint_chernoff_prepare(512, u_step, e_step, logpotential,
                                      tilt_logs, log_hole_potential)
    logs = endlogs = None
    for tilt, fraction in enumerate(np.linspace(0, 1, tilts)):
        started = time.monotonic()
        odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
        p = odds / (1 + odds)
        logpi = choose + h * np.log2(p) + (8192 - h) * np.log2(1 - p)
        # alpha=1 and output_discount=1: no field/output count is hidden here.
        comparison, terminal, _ = joint_chernoff_comparison(prepared, m, b, p,
            (1 - p) / a, JBOX_LD(1), chi, JBOX_LD(1))
        region, scale, last, last_scale = joint_subset_region(comparison, terminal, prepared['zero_vector'])
        candidate = np.full_like(region, JBOX_LD('-inf'))
        mask = region > 0
        candidate[mask] = np.log2(region[mask]) + scale - logpi
        candidate_end = np.full_like(last, JBOX_LD('-inf'))
        mask = last > 0
        candidate_end[mask] = np.log2(last[mask]) + last_scale - logpi
        logs = candidate if logs is None else np.minimum(logs, candidate)
        endlogs = candidate_end if endlogs is None else np.minimum(endlogs, candidate_end)
        if progress:
            print(dict(joint_box_quota_tilt=tilt + 1, dimension=len(comparison),
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    scale, end_scale = logs.max(), endlogs.max()
    region, last = np.exp2(logs - scale), np.exp2(endlogs - end_scale)
    if np.any(np.isfinite(logs) & (region == 0)) or np.any(np.isfinite(endlogs) & (last == 0)):
        raise ArithmeticError('Geometric regional normalization discarded positive terms')
    tail = joint_subset_tail(region, scale, last, end_scale, prepared['zero_vector'])
    return choose - 128 * h * loga + tail


def joint_boxes_checks():
    _joint_boxes_dependencies()
    F = Fraction
    a, bulk, birth_ratio, chi, p = F(4), F(2), F(1, 2), F(3, 4), F(1, 3)
    local = joint_mds_local(bulk * birth_ratio, bulk, p, (1 - p) / a,
                            F(1), chi, exact=True)
    states = tuple((u, e) for u in range(5) for e in range(5 - u))
    index = {state: j for j, state in enumerate(states)}
    expected = np.full_like(local, F(0))
    # Independently enumerate masks and raw routing labels.  Profile weights
    # are a^-outside * b^available * u^forced * chi^covered, with no theta.
    for row, (u, e) in enumerate(states):
        old_a, old_e = (1 << u) - 1, ((1 << e) - 1) << u
        if old_a and old_e:
            continue
        for new_a in range(16):
            if not old_a and old_e and new_a == 15:
                continue
            covered = 15 if old_a else new_a
            new_e = covered & ~new_a
            available, forced, footprint = new_a.bit_count(), (0 if old_a else new_a.bit_count()), covered.bit_count()
            for eligible_mask in range(16):
                outside = covered & ~eligible_mask
                if old_a:
                    if (new_e & ~eligible_mask).bit_count() >= u:
                        continue
                elif outside:
                    continue
                route = p ** eligible_mask.bit_count() * (1 - p) ** (4 - eligible_mask.bit_count())
                expected[row, index[available, new_e.bit_count()]] += (
                    route * a ** -outside.bit_count() * bulk ** available
                    * birth_ratio ** forced * chi ** footprint)
    assert np.array_equal(local, expected)
    half = 1024
    box = lambda center: tuple((value - half, value + half) for value in center)
    for center in JBOX_LEGACY:
        assert joint_boxes_canonical_screen(2048, box(center))['status'] == 'empty'
    assert joint_boxes_canonical_screen(2048, box(JBOX_OBSERVED))['status'] == 'not_excluded'
    assert joint_boxes_canonical_screen(2048, box(JBOX_CANDIDATE))['status'] == 'not_excluded'
    old_hard = box(JBOX_LEGACY[-1])
    markers = (F(2), F(25, 16), F(-7, 16), F(0))
    assert joint_boxes_minimum(markers, old_hard) == -406596
    # Check interval propagation's rounding in both coefficient directions.
    one = tuple((value, value) for value in JBOX_CANDIDATE)
    assert joint_boxes_canonical_screen(2048, one)['bounds'] == one
    print(dict(exact_geometric_local_entries=local.size, empty_legacy_boxes=3,
               aggregate_not_excluded_boxes=2, exact_box_minimum=True,
               field_prefactors_in_geometry=False, full_code_certificate=False), flush=True)


def joint_boxes_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--center', type=int, nargs=4, metavar=('R', 'L', 'f', 'F'))
    parser.add_argument('--halfwidth', type=int, default=1024)
    parser.add_argument('--loga', type=float, default=2.)
    parser.add_argument('--logb', type=float)
    parser.add_argument('--logu', type=float, default=-.4375)
    parser.add_argument('--logchi', type=float, default=0.)
    parser.add_argument('--u-step', type=int, default=32)
    parser.add_argument('--e-step', type=int, default=32)
    parser.add_argument('--logpotential', type=float, default=0.)
    parser.add_argument('--log-hole-potential', type=float, default=0.)
    parser.add_argument('--chernoff-logs', type=float, nargs='+', default=(-1/16, 0., 1/16))
    parser.add_argument('--tilts', type=int, default=3)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            joint_boxes_checks()
            return
        _joint_boxes_dependencies()
        if args.halfwidth < 0:
            raise ValueError('Box half-width must be nonnegative')
        centers = [tuple(args.center)] if args.center else [JBOX_OBSERVED, JBOX_CANDIDATE, *JBOX_LEGACY]
        candidates = []
        for center in centers:
            box = tuple((value - args.halfwidth, value + args.halfwidth) for value in center)
            screen = joint_boxes_canonical_screen(args.h, box)
            print(dict(h=args.h, center_R_L_f_F=center, canonical_screen=screen), flush=True)
            if screen['status'] != 'empty':
                candidates.append((center, screen['bounds']))
        if not candidates:
            print(dict(skipped_geometric_moment=True, reason='Every requested box is excluded'), flush=True)
            return
        # This is the frozen old geometric marker, not a field-cost factor.
        logb = (JBOX_LD(1) / 16 * np.log2(1 + (JBOX_LD(2) ** 32 + 1) * JBOX_LD(2) ** -7)
                if args.logb is None else JBOX_LD(args.logb))
        markers = (JBOX_LD(args.loga), logb, JBOX_LD(args.logu), JBOX_LD(args.logchi))
        started = time.monotonic()
        moment = joint_boxes_geometry_logmoment(args.h, *markers, tilts=args.tilts,
            u_step=args.u_step, e_step=args.e_step, logpotential=args.logpotential,
            log_hole_potential=args.log_hole_potential, tilt_logs=args.chernoff_logs, progress=True)
        for center, box in candidates:
            print(dict(h=args.h, center_R_L_f_F=center, effective_box=box,
                       geometric_log2_moment=float(moment), markers=tuple(map(float, markers)),
                       **box_diagnostic(moment, markers, box)), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3),
                   field_prefactors_in_geometry=False, full_code_certificate=False), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    joint_boxes_main()
