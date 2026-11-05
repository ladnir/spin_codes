"""Bounded FLOAT search for profile-specific run/occupancy weights.

Inject quartet_run_occupancy.py and quartet_global_profile_probe.py.
The positive-vector residual includes the free-terminal factor. Floating
arithmetic is exploratory only: it does not certify any profile or code.
Both availability and actual-state variants are explicitly labelled.
"""

import argparse
import os
import shutil
import threading
import time
import numpy as np
from fractions import Fraction


def moment_search_nonempty_region(matrix, scale):
    """Positive normalized K^16-Z; no nearly equal scalar subtraction."""
    positive = matrix.copy()
    positive[0, 0] = 0  # K00=1 is exactly the single empty update.
    for _ in range(4):
        square = positive @ positive
        edges = np.zeros_like(positive)
        edges[0, :] += positive[0, :]
        edges[:, 0] += positive[:, 0]
        common = max(2 * scale, scale)
        positive = square * np.exp2(2 * scale - common) + edges * np.exp2(scale - common)
        norm = positive.sum(axis=-1).max()
        positive /= norm
        scale = common + np.log2(norm)
    return positive, scale


def moment_search_upper(prepared, log_s, log_t, log_z, iterations=96, log_birth=None,
                        minimum_regions=None, log_region=0.):
    if log_birth is None:
        matrix, scale = run_occupancy_matrix(512, log_s, log_t, log_z, prepared)
    else:
        if log_birth > 0:
            raise ValueError('Positive homogenization requires birth marker <=1')
        s, t, z, birth = (np.longdouble(2) ** a for a in (log_s, log_t, log_z, log_birth))
        matrix, _ = direct_available_matrix(512, s * birth, s, 1, z, t, prepared=prepared)
        norm = matrix.sum(axis=-1).max()
        matrix /= norm
        scale = np.log2(norm)
    updates, region_correction = 4096, 0.
    if minimum_regions is not None:
        positive, positive_scale = moment_search_nonempty_region(matrix, scale)
        updates = 256
        if minimum_regions == 256:
            matrix, scale = positive, positive_scale
        else:
            scale = max(np.longdouble(0), positive_scale + log_region)
            matrix = positive * np.exp2(positive_scale + log_region - scale)
            matrix[0, 0] += np.exp2(-scale)
            region_correction = minimum_regions * log_region
    vector = np.ones(513, dtype=np.longdouble)
    for _ in range(iterations):
        vector = matrix @ vector
        vector /= vector.max()
    if np.any(vector <= 0) or not np.all(np.isfinite(vector)):
        raise ArithmeticError('Invalid FLOAT witness')
    residual = (matrix @ vector) / vector
    growth = scale + np.log2(residual.max())
    boundary = np.log2(vector[0] / vector.min())
    return float(updates * growth + boundary - region_correction), float(growth), float(boundary)


def moment_search_profile(profile, exponent=16, theta=1., actual=False, rounds=4, birth=False, regions=False):
    h, R, L, f, F = profile
    C = F - R + 128 * h
    route = min(0., gp_logchoose(F - f, C - f)
                + gp_logchoose(GP_N - F, 256 * h - C) - 255 * gp_logchoose(8192, h))
    x = np.longdouble(2) ** -exponent
    q = np.longdouble(1 << 32)
    if actual:
        lam = 1 + q * x
        d = q * x * x / lam
        field = -31 - 32 * R + exponent * GP_W + L * np.log2(lam)
    else:
        c = (q + 1) * x
        lam = 1 + c
        d = max(min(np.longdouble(1), lam / (q + 1)), (1 + q * x * x) / lam)
        field = -31 - 32 * R + exponent * GP_W + f * np.log2(c) + (L - f) * np.log2(lam)
    log_t = float(theta * np.log2(d))
    prepared = direct_available_prepare(512) if birth else run_occupancy_prepare(512)
    cache = {}
    G = F - f
    rho = max(129, (C + h - 1) // h,
              (R - 128 * h + 8192 - h - 1) // (8192 - h) if h < 8192 else 0) if regions else None
    def value(s, z, b=0., k=0.):
        key = s, z, b, k
        if key not in cache:
            moment, growth, boundary = moment_search_upper(prepared, s, log_t, z, log_birth=b if birth else None,
                                                          minimum_regions=rho, log_region=k)
            cache[key] = dict(log_s=s, log_z=z, log_birth=b, log_region=k, moment=moment, growth=growth, boundary=boundary,
                              total=float(moment - s * L - z * G - b * f + route + theta * field))
        return cache[key]
    # A small deterministic pattern search; the returned bound, not stationarity,
    # is the diagnostic. Retain the measured residual after every proposal.
    starts = ((.25, -.25), (1., -1.), (2., -2.), (4., -4.),
              (-2., 2.), (-4., 4.), (-8., 8.), (-10., 10.))
    best = min((value(s, z, b, 8. if regions and rho < 256 else 0.) for s, z in starts
                for b in ((0., -4., -8.) if birth else (0.,))),
               key=lambda a: a['total'])
    step = .5
    for _ in range(rounds):
        previous = best
        s, z, b, k = best['log_s'], best['log_z'], best['log_birth'], best['log_region']
        candidates = [best] + [value(s + ds * step, z + dz * step, b, k)
                             for ds, dz in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1))]
        if birth:
            candidates += [value(s, z, b + db, k) for db in (-2., 2.) if b + db <= 0]
        if regions and rho < 256:
            candidates += [value(s, z, b, k + dk) for dk in (-4., 4.) if k + dk >= 0]
        best = min(candidates, key=lambda a: a['total'])
        if best is previous:
            step *= .5
    return dict(profile=profile, actual_state_support=actual, exponent=exponent, theta=theta,
                log_t=log_t, route=route, field_prefactor=float(theta * field),
                evaluations=len(cache), retaining_actual_forced_births=birth, minimum_occupied_regions=rho,
                best=best, float_diagnostic_only=True,
                full_code_certificate=False)


def moment_search_checks():
    entries = regions = 0
    for width in (4, 8):
        for s, t, z in ((Fraction(3, 4), Fraction(1, 2), Fraction(5, 4)),
                        (Fraction(5, 4), Fraction(1, 16), Fraction(3, 4))):
            exact, _ = _direct_exact_matrix(width, s, s, Fraction(1), z, t, Fraction(1))
            exact = np.array(exact, dtype=object)
            direct, _ = direct_available_matrix(width, s, s, 1, z, t)
            old, scale = run_occupancy_matrix(width, np.log2(np.longdouble(float(s))),
                np.log2(np.longdouble(float(t))), np.log2(np.longdouble(float(z))))
            truth = np.array([[np.longdouble(a.numerator) / np.longdouble(a.denominator)
                               for a in row] for row in exact])
            assert np.all(abs(direct - truth) <= 1e-13 * np.maximum(1, truth))
            assert np.all(abs(old * np.exp2(scale) - truth) <= 1e-13 * np.maximum(1, truth))
            entries += (width + 1) ** 2
            exact16 = exact.copy()
            for _ in range(4):
                exact16 = exact16 @ exact16
            exact16[0, 0] -= 1
            positive, pscale = moment_search_nonempty_region(old, scale)
            truth16 = np.array([[np.longdouble(a.numerator) / np.longdouble(a.denominator)
                                 for a in row] for row in exact16])
            assert np.all(abs(positive * np.exp2(pscale) - truth16) <= 1e-12 * np.maximum(1, truth16))
            regions += (width + 1) ** 2
    print(dict(exact_reference_birth_neutral_entries=entries,
               exact_reference_positive_region_entries=regions,
               floating_comparisons_with_stated_tolerances=True, full_code_certificate=False), flush=True)


def moment_search_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', nargs=5, type=int, default=(7168, 917504, 1310720, 33280, 1720320))
    parser.add_argument('--exponent', type=int, default=16)
    parser.add_argument('--theta', type=float, default=1.)
    parser.add_argument('--actual', action='store_true')
    parser.add_argument('--birth', action='store_true')
    parser.add_argument('--regions', action='store_true')
    parser.add_argument('--rounds', type=int, default=4)
    parser.add_argument('--checks', action='store_true')
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            moment_search_checks()
            return
        print(dict(seconds=round(time.monotonic() - started, 3),
                   **moment_search_profile(tuple(args.profile), args.exponent, args.theta, args.actual, args.rounds, args.birth, args.regions)),
              flush=True)
        print(dict(elapsed_seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    moment_search_main()
