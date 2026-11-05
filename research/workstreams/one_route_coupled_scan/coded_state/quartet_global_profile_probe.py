"""Bounded FLOAT map of the complement of the global rank certificates.

No code certificate is claimed.  A fixed exact run/occupancy moment is
combined with the canonical closure interval, all prior linear census
classes, forced-birth routing, and fractional conditional field bounds.
The finite grid identifies candidate gaps; it cannot establish coverage.
Run only on the designated remote numerical host, serially and file-free.
"""

import argparse
from collections import Counter
from math import ceil, floor, lgamma, log2
import os
import shutil
import threading
import time

GP_N = 1 << 21
GP_K = 1 << 20
GP_W = 943718
GP_Q = 1 << 32
GP_S = log2(75 / 64)
GP_T = log2((GP_Q + 65536) / (256 * (GP_Q + 257)))
GP_M = 22323  # independently integer-verified; this probe remains FLOAT


def gp_census(j):
    if j == 0:
        return GP_N
    if j == 3:
        return 361984
    return 512 * j + (3 * GP_N + (1 << (j + 1)) - 1) // (1 << (j + 1))


def gp_interval(h, R, L, f, F):
    """Exact necessary screens; intersect complements of prior certificates."""
    C = F - R + 128 * h
    O = R - 128 * h
    if not (1 <= h <= 8192 and 128 * h <= R <= GP_K
            and 1 <= f <= L <= F <= GP_N and h <= F - R <= 128 * h
            and f <= C and 2 * f - 512 <= GP_W and (F - f) % 4 == 0
            and F - f <= GP_N - 512 and L - 512 <= F - f <= 4 * L
            and 2 * F - L - f <= GP_N - 512
            and F - f + 4 * max(0, O - L + f) <= 4 * L):
        return 'empty', None
    if F >= 1736704 or R >= 1009256 or R - L >= 65537:
        return 'prior_class', None
    lo, hi = max(0, f - 512), min(F - L, L, GP_N - L)
    # Endpoint capacity: f-D <= L-G/4-max(0,O-L+f).
    terminal_max = min(512, L - (F - f) // 4 - max(0, O - L + f))
    if terminal_max < 0:
        return 'empty', None
    lo = max(lo, f - terminal_max)
    if 2 * R + max(0, f - 512) >= 2 * 1009257:
        return 'prior_class', None
    for j in range(17):
        c = gp_census(j)
        # Complement: 32(R-L)-2jD <= c+44.
        if j:
            lo = max(lo, -((c + 44 - 32 * (R - L)) // (2 * j)))
        elif 32 * (R - L) > c + 44:
            return 'prior_class', None
        # Complement: (15-2j)D <=16W+c+91-32R+16L.
        a = 15 - 2 * j
        b = 16 * GP_W + c + 91 - 32 * R + 16 * L
        if a > 0:
            hi = min(hi, b // a)
        else:
            lo = max(lo, -((-b) // a))
    if lo > hi:
        return 'prior_class', None
    return 'remaining', (lo, hi)


def gp_logchoose(n, k):
    if not 0 <= k <= n:
        return float('-inf')
    if k == 0 or k == n:
        return 0.
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / 0.6931471805599453


def gp_field_markers():
    result = []
    for exponent in range(33):
        x = 2. ** -exponent
        c = (GP_Q + 1) * x
        lam = 1 + c
        d = max(min(1., lam / (GP_Q + 1)), (1 + GP_Q * x * x) / lam)
        ld = log2(d)
        theta = GP_T / ld if ld < GP_T else 1.
        result.append((exponent, log2(c), log2(lam), ld, theta))
    return result


GP_MARKERS = gp_field_markers()


def gp_value(h, R, L, f, F):
    status, interval = gp_interval(h, R, L, f, F)
    if status != 'remaining':
        return dict(status=status)
    C = F - R + 128 * h
    route = min(0., gp_logchoose(F - f, C - f)
                + gp_logchoose(GP_N - F, 256 * h - C) - 255 * gp_logchoose(8192, h))
    geometry = GP_M - GP_S * L + (F - f) / 4 + route
    lo, hi = interval
    # theta=0 uses only geometry and the closure upper endpoint.
    best, choice = -GP_T * hi, (0, 0.)
    for e, lc, ll, ld, middle in GP_MARKERS:
        field = -31 - 32 * R + e * GP_W + f * lc + (L - f) * ll
        for theta in (1., middle):
            residual = theta * ld - GP_T
            value = theta * field + residual * (hi if residual >= 0 else lo)
            if value < best:
                best, choice = value, (e, theta)
    total = geometry + best
    occupancy_only = None
    if 'occupancy_geometry_logbound' in globals():
        occupancy_only = occupancy_geometry_logbound(L, F - f) + route
        total = min(total, occupancy_only)
    return dict(status='remaining', h=h, R=R, L=L, f=f, F=F, D_interval=interval,
                route=route, geometry=geometry, field_and_interval=best,
                occupancy_only=occupancy_only,
                value=total, choice=choice, float_diagnostic_only=True)


def gp_checks():
    assert gp_interval(2048, 900000, 900000, 20000, 1200000)[0] == 'empty'
    status, bounds = gp_interval(2048, 900000, 900000, 20000, 1150000)
    assert status == 'remaining' and bounds == (19488, 117951)
    # Check signed integer complement arithmetic directly.
    tests = 0
    for a in tuple(range(-17, 0)) + tuple(range(1, 18)):
        for b in range(-51, 52):
            for d in range(31):
                bound = b // a if a > 0 else -((-b) // a)
                assert (a * d <= b) == (d <= bound if a > 0 else d >= bound)
                tests += 1
    print(dict(exact_complement_checks=tests, full_code_certificate=False), flush=True)


def gp_probe(top=12):
    counts = Counter()
    worst = []
    hvalues = (1, 4, 16, 64, 128, 256, 512, 768, 1024, 1536, 2048,
               2560, 3072, 3584, 4096, 5120, 6144, 7168)
    for h in hvalues:
        local = []
        for R in sorted({128 * h, 1009255, *range(max(128 * h, 32768), 1009256, 32768)}):
            if R > 1009255:
                continue
            for extra_num in (1, 16, 32, 64, 96, 112, 120, 124, 128):
                F = R + h * extra_num
                if F >= 1736704:
                    continue
                cross = R + (F - GP_W) // 4
                for L in sorted({R, GP_W, F, cross, cross + 1,
                                 *range(max(1, R - 65536), F + 1, 32768)}):
                    if not 1 <= L <= F:
                        continue
                    fmin = max(1, 2 * F - L - GP_N + 512, F - 4 * L)
                    fmax = min(L, F - L + 512, (GP_W + 512) // 2, F - R + 128 * h)
                    for rawf in sorted({fmin, fmax, 512, 16384, 65536, 131072, 262144}):
                        f = rawf + (F - rawf) % 4
                        if not fmin <= f <= fmax:
                            continue
                        value = gp_value(h, R, L, f, F)
                        counts[value['status']] += 1
                        if value['status'] == 'remaining':
                            counts['float_negative' if value['value'] < -172 else 'float_open'] += 1
                            local.append(value)
        local.sort(key=lambda x: x['value'], reverse=True)
        worst.extend(local[:top])
        print(dict(h=h, candidates=len(local), worst=local[0] if local else None), flush=True)
    worst.sort(key=lambda x: x['value'], reverse=True)
    print(dict(counts=dict(counts), finite_grid_only=True, full_code_certificate=False), flush=True)
    for entry in worst[:top]:
        print(entry, flush=True)


def gp_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', action='store_true')
    parser.add_argument('--probe', action='store_true')
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
            gp_checks()
        elif args.probe:
            gp_probe()
        else:
            parser.error('Choose --checks or --probe')
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    gp_main()
