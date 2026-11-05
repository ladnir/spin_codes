"""Small tensor-norm diagnostic for the exact quartet availability interface.

The local physical operator has 81 configurations.  Its weighted Euclidean
norm equals that of the 15-state orthonormal count matrix.  Orthogonal
symmetrization gives a width-512 norm bound by the 128th power of the local
norm.  Regional quota conditioning costs one inverse binomial probability.
This ALL-REGION norm control does not retain the >=129 occupied-region
restriction.  It is a FLOAT proposal, never a directed certificate.

Standalone; no helper imports or subprocesses.  Run numerical jobs serially.
"""

import argparse
from fractions import Fraction
from itertools import product
from math import comb, lgamma
import os
import shutil
import threading
import time

import numpy as np


SP_LD = np.longdouble
SP_N = 1 << 21
SP_W = 943718
SP_FCAP = 1736703
SP_RCAP = 1009255
SP_STATES = tuple((u, e) for u in range(5) for e in range(5 - u))
SP_INDEX = {state: index for index, state in enumerate(SP_STATES)}
SP_MULTIPLICITY = tuple(comb(4, u) * comb(4 - u, e) for u, e in SP_STATES)


def _sp_ld(value):
    if isinstance(value, Fraction):
        return SP_LD(value.numerator) / SP_LD(value.denominator)
    return SP_LD(value)


def joint_spectral_local(m, b, p, r, alpha, coverage=1, exact=False):
    """Positive local count kernel, including all zero rows and endpoints."""
    cast = Fraction if exact else _sp_ld
    dtype = object if exact else SP_LD
    m, b, p, r, alpha, coverage = map(cast, (m, b, p, r, alpha, coverage))
    if min(m, b, p, r, alpha, coverage) < 0:
        raise ValueError('Local weights must be nonnegative')
    matrix = np.zeros((15, 15), dtype=dtype)
    birth, route4 = m * p * coverage, (r * coverage) ** 4
    for row, (u, e) in enumerate(SP_STATES):
        if u:
            if e:
                continue  # Old A and E cannot share this quartet.
            for v in range(5):
                coefficient = sum(comb(u, k) * comb(4 - u, v - k) * alpha ** (u - k)
                                  for k in range(max(0, v - (4 - u)), min(u, v) + 1))
                matrix[row, SP_INDEX[v, 4 - v]] = route4 * b ** v * coefficient
        else:
            # A full fresh birth quartet is forbidden when old E is present.
            for v in range(5 if e == 0 else 4):
                matrix[row, SP_INDEX[v, 0]] = comb(4, v) * birth ** v
    if not exact and (np.any(matrix < 0) or not np.all(np.isfinite(matrix))):
        raise ArithmeticError('Invalid positive local matrix')
    return matrix


def joint_spectral_log_norm(local, logt=0, logs=0):
    """log2 ||S D^-1 K D S^-1||_2, S_ii=sqrt(multinomial_i).

    Singular values use a scaled float64 matrix because numpy.linalg does
    not support longdouble SVD.  Results remain explicitly diagnostic.
    """
    logt, logs = map(SP_LD, (logt, logs))
    logs_matrix = np.full((15, 15), SP_LD('-inf'), dtype=SP_LD)
    positive = local > 0
    logs_matrix[positive] = np.log2(local[positive])
    counts = np.array(SP_STATES, dtype=SP_LD)
    diagonal = counts[:, 0] * logt + counts[:, 1] * logs
    log_mult = np.log2(np.array(SP_MULTIPLICITY, dtype=SP_LD)) / 2
    logs_matrix += ((log_mult - diagonal)[:, None]
                    - (log_mult - diagonal)[None, :])
    scale = logs_matrix.max()
    if not np.isfinite(scale):
        return SP_LD('-inf')
    shifted = logs_matrix[positive] - scale
    # Do not silently discard a positive coefficient in a floating proposal.
    if np.any(shifted < -1000):
        raise ArithmeticError('Scaled SVD would lose positive entries; narrow the potential range')
    scaled = np.zeros((15, 15), dtype=np.float64)
    scaled[positive] = np.asarray(np.exp2(shifted), dtype=np.float64)
    singular = np.linalg.svd(scaled, compute_uv=False)[0]
    if not np.isfinite(singular) or singular <= 0:
        raise ArithmeticError('Invalid local singular value')
    return scale + np.log2(SP_LD(singular))


def _sp_physical_local(m, b, p, r, alpha, coverage):
    """Independent rational configuration operator on 3^4 local states."""
    states = tuple(product(range(3), repeat=4))  # 0=uncovered, 1=A, 2=E
    matrix = np.full((81, 81), Fraction(0), dtype=object)
    masks = [(sum(1 << j for j, value in enumerate(state) if value == 1),
              sum(1 << j for j, value in enumerate(state) if value == 2)) for state in states]
    for i, (old_a, old_e) in enumerate(masks):
        for j, (new_a, new_e) in enumerate(masks):
            if old_a:
                if old_e or new_e != (15 ^ new_a):
                    continue
                matrix[i, j] = ((r * coverage) ** 4 * b ** new_a.bit_count()
                                * alpha ** (old_a & ~new_a).bit_count())
            elif new_e == 0 and not (old_e and new_a == 15):
                matrix[i, j] = (m * p * coverage) ** new_a.bit_count()
    return states, masks, matrix


def joint_spectral_checks():
    F = Fraction
    entries = norm_checks = terminal_checks = zero_checks = 0
    for alpha, coverage in ((F(2, 3), F(4, 5)), (F(0), F(1)),
                            (F(2), F(1)), (F(2, 3), F(0))):
        m, b, p, r = F(3, 4), F(5, 4), F(1, 3), F(2, 5)
        exact = joint_spectral_local(m, b, p, r, alpha, coverage, exact=True)
        _, masks, physical = _sp_physical_local(m, b, p, r, alpha, coverage)
        for row, (old_a, old_e) in enumerate(masks):
            expected = [F(0)] * 15
            for col, (new_a, new_e) in enumerate(masks):
                expected[SP_INDEX[new_a.bit_count(), new_e.bit_count()]] += physical[row, col]
            assert expected == list(exact[SP_INDEX[old_a.bit_count(), old_e.bit_count()]])
            entries += 15
        floating = joint_spectral_local(m, b, p, r, alpha, coverage)
        for lt, ls in ((0, 0), (-1, 1), (1, -1)):
            # Compare the full 81-configuration norm against the 15-count norm.
            values = np.array([[_sp_ld(value) for value in row] for row in physical], dtype=np.float64)
            diagonal = np.array([2. ** (lt * a.bit_count() + ls * e.bit_count()) for a, e in masks])
            transformed = values * diagonal[None, :] / diagonal[:, None]
            truth = np.log2(np.linalg.svd(transformed, compute_uv=False)[0])
            actual = joint_spectral_log_norm(floating, lt, ls)
            assert abs(actual - truth) < SP_LD('2e-12')
            norm_checks += 1
        if coverage == 0:
            # The transformed zero operator has norm (1+s^-2)^(4/2).
            for ls in (-1, 0, 1):
                actual = joint_spectral_log_norm(floating, 0, ls)
                expected = 2 * np.log2(SP_LD(1) + SP_LD(2) ** (-2 * ls))
                assert abs(actual - expected) < SP_LD('2e-12')
                zero_checks += 1
    for width in (4, 8):
        for t in (F(1, 2), F(1), F(2)):
            squared = sum(F(comb(width, u)) * t ** (-2 * u) for u in range(width + 1))
            assert squared == (1 + t ** -2) ** width
            terminal_checks += 1
    print(dict(exact_physical_count_entries=entries, norm_normalization_checks=norm_checks,
               general_zero_norm_checks=zero_checks, free_terminal_checks=terminal_checks,
               all_regions_only=True, occupied_129_retained=False,
               full_code_certificate=False), flush=True)


def joint_spectral_value(h=2048, theta=1/16, logz=-.4375, logt=0, logs=0,
                         coverage=0, outside=0, p=None, tilts=9):
    """All-region norm bound; coverage/outside arguments are log2 markers.

    The known remainder has F<=1736703 and C_out<=1009255-128h.
    Nonpositive markers therefore contribute -Fcap*coverage-Ccap*outside
    to the external logarithm.  Potentials t,s are common to every tilt.
    """
    if not 0 < h < 8192 or not 0 < theta <= 1 or logz > 0 or tilts < 1:
        raise ValueError('Invalid field/row/tilt parameters')
    if coverage > 0 or outside > 0:
        raise ValueError('Coverage and outside log2 markers must be nonpositive')
    if p is not None and not 0 < p < 1:
        raise ValueError('Explicit p must lie strictly between zero and one')
    theta, logz, logt, logs, coverage, outside = map(SP_LD,
        (theta, logz, logt, logs, coverage, outside))
    q = SP_LD(2) ** 32
    x = np.exp2(logz / theta)
    c = (q + 1) * x
    lam = 1 + c
    d = max(min(SP_LD(1), lam / (q + 1)), (1 + q * x * x) / lam)
    m, b, alpha = c ** theta, lam ** theta, d ** theta
    chi = np.exp2(coverage)
    outside_weight = np.exp2(outside - 32 * theta)
    actual = SP_LD(h) / 8192
    effective_field = 1 / outside_weight
    low = actual / (effective_field - (effective_field - 1) * actual)
    choose = SP_LD((lgamma(8193) - lgamma(h + 1) - lgamma(8193 - h)) / float(np.log(SP_LD(2))))
    candidates = [SP_LD(p)] if p is not None else []
    if p is None:
        for fraction in np.linspace(0, 1, tilts):
            odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
            candidates.append(odds / (1 + odds))
    best = None
    for probability in candidates:
        logpi = choose + h * np.log2(probability) + (8192 - h) * np.log2(1 - probability)
        r = probability + (1 - probability) * outside_weight
        local = joint_spectral_local(m, b, probability, r, alpha, chi)
        lognorm = joint_spectral_log_norm(local, logt, logs)
        # 16 updates * 128 quartets, then one region's quota conditioning.
        logregion = 2048 * lognorm - logpi
        if best is None or logregion < best[0]:
            best = (logregion, probability, lognorm, logpi)
    terminal = 256 * np.logaddexp2(SP_LD(0), -2 * logt)
    outside_cap = SP_RCAP - 128 * h
    external = (theta - 32 * theta - SP_W * logz + choose - 32 * theta * 128 * h
                - outside_cap * outside - SP_FCAP * coverage)
    bound = external + 256 * best[0] + terminal
    # An invariant diagonal gives a p-independent floor on THIS diagnostic.
    diagonal_floor = (external + SP_N * (np.log2(b) + coverage)
                      + (SP_N - 256 * h) * (outside - 32 * theta) + terminal)
    return float(bound), dict(h=h, theta=float(theta), logz=float(logz),
        logt=float(logt), logs=float(logs), logcoverage=float(coverage), logoutside=float(outside),
        p=float(best[1]), log2_local_norm=float(best[2]), log2_quota_probability=float(best[3]),
        log2_terminal=float(terminal), log2_diagonal_floor=float(diagonal_floor),
        tilts=len(candidates), local_dimension=15, all_regions_only=True,
        occupied_129_retained=False, certificate=False)


def joint_spectral_search(h=2048, theta=1/16, logz=-.4375, tilts=9):
    """One explicitly bounded 144-point marked proposal grid (1296 tiny SVDs)."""
    best = []
    for logt, logs, coverage, outside in product((-1., 0., 1.), (-1., 0., 1.),
            (0., -1/16, -1/8, -1/4), (0., -1/32, -1/16, -1/8)):
        value, detail = joint_spectral_value(h, theta, logz, logt, logs, coverage, outside, tilts=tilts)
        best.append((value, detail))
    best.sort(key=lambda item: item[0])
    return best[:5]


def joint_spectral_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    mode.add_argument('--search', action='store_true')
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--logz', type=float, default=-.4375)
    parser.add_argument('--logt', type=float, default=0.)
    parser.add_argument('--logs', type=float, default=0.)
    parser.add_argument('--coverage', type=float, default=0.)
    parser.add_argument('--outside', type=float, default=0.)
    parser.add_argument('--p', type=float)
    parser.add_argument('--tilts', type=int, default=9)
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
            joint_spectral_checks()
        elif args.search:
            for bound, detail in joint_spectral_search(args.h, args.theta, args.logz, args.tilts):
                print(dict(log2_bound=bound, **detail), flush=True)
            print(dict(search_points=144, seconds=round(time.monotonic() - started, 3)), flush=True)
        else:
            bound, detail = joint_spectral_value(args.h, args.theta, args.logz, args.logt,
                args.logs, args.coverage, args.outside, args.p, args.tilts)
            print(dict(log2_bound=bound, seconds=round(time.monotonic() - started, 3), **detail), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    joint_spectral_main()
