"""Tiny quota-coefficient envelope on the restricted O,P,H kernel only.

At time t, f_t(s,r)<=min_l d_l(s)*x_l^-r on the feasible count interval
 I_t=[max(0,h-n*(T-t)),min(h,n*t)], and f_t is zero outside it.
Apply c_l=L(x_l)d_l, then set
 d'_l(s)=max_{r in I_t} min_j c_j(s)*(x_l/x_j)^r.
The continuous maximum is at an endpoint or a pairwise line crossing;
it conservatively bounds the integer maximum. On the interval the
represented minimum is unchanged. The possible losses occur when
moving a minimum through L and extending an envelope outside I_t.

This gives an UPPER comparison for the restricted three-state quota
matrix, not a lower bound or an upper bound for the full physical
operator. The optional displayed tail expression concerns only
e_H B_restricted^256 1. For the full physical operator b=B1 is exact,
but its first-layer endpoint is not the restricted L_x*1.

Fixed seven x values, y2, one initial witness and at most12 geometric
half-updates. No full physical kernel, large arrays, files, or import
jobs. FLOAT controls are not certificates. Parent runs everything.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _qe_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_qe_thread_var] = '1'

import numpy as np

QE_LD = np.longdouble
QE_X = ('1', '1.6', '1.86', '2.04', '2.45', '2.76', '4')
QE_RAW = ((255, 204, 0, 0, 0), (0, 816, 1530, 1020, 255),
          tuple(comb(4, j)*(256//4**j-1) for j in range(5)),
          tuple(comb(4, j)*(256-256//4**j) for j in range(5)),
          tuple(comb(4, j)*(4096//4**j-16//2**j) for j in range(5)),
          tuple(comb(4, j)*(4096//2**j-4096//4**j) for j in range(5)))
QE_EDGES = ((1, 0, 0), (0, 0, 1), (3, 1, 0), (2, 1, 2), (5, 2, 0), (4, 2, 2))


def _qe_logs(values):
    return np.asarray([np.log2(QE_LD(v.numerator)/QE_LD(v.denominator)) for v in values], dtype=QE_LD)


def _qe_logchoose(n, k):
    value = comb(n, k)
    shift = max(0, value.bit_length()-64)
    # FLOAT display/comparison only; no directed rounding claim.
    return shift+np.log2(QE_LD(value >> shift))


def _qe_matvec(log_matrices, log_vectors):
    """Matrices indexed(mode,source,destination), vectors(state,mode)."""
    terms = log_matrices+log_vectors.T[:, None, :]
    peaks = terms.max(axis=2)
    if not np.all(np.isfinite(peaks)):
        raise ArithmeticError('Every restricted kernel row must have positive mass')
    positive = np.isfinite(terms)
    weights = np.exp2(terms-peaks[:, :, None])
    if np.any(positive & (weights <= 0)):
        raise ArithmeticError('Positive three-state action term underflow')
    result = (peaks+np.log2(weights.sum(axis=2, dtype=QE_LD))).T
    if not np.all(np.isfinite(result)):
        raise ArithmeticError('Invalid positive action logs')
    return result


def _qe_project(log_c, log_x, lower, upper):
    """Continuous quota maximum: endpoints plus all clipped crossings."""
    if lower > upper:
        raise ValueError('Empty feasible coefficient interval')
    projected = np.empty_like(log_c)
    for state, row in enumerate(log_c):
        candidates = [QE_LD(lower), QE_LD(upper)]
        for j in range(len(log_x)):
            for k in range(j):
                crossing = (row[j]-row[k])/(log_x[j]-log_x[k])
                candidates.append(max(QE_LD(lower), min(QE_LD(upper), crossing)))
        r = np.asarray(candidates, dtype=QE_LD)
        envelope = (row[None, :]-r[:, None]*log_x).min(axis=1)
        projected[state] = (envelope[:, None]+r[:, None]*log_x).max(axis=0)
    if not np.all(np.isfinite(projected)):
        raise ArithmeticError('Quota-envelope projection produced invalid logs')
    return projected


def quota_envelope_action(log_matrices, log_x, log_phi, degree, updates, quota, log_C):
    log_phi = np.asarray(log_phi, dtype=QE_LD)
    log_x = np.asarray(log_x, dtype=QE_LD)
    if (not np.all(np.isfinite(log_phi)) or not np.all(np.isfinite(log_x))
            or len(np.unique(log_x)) != len(log_x)):
        raise ValueError('Use positive finite input and distinct positive quota evaluations')
    d = np.repeat(log_phi[:, None], len(log_x), axis=1)
    for t in range(1, updates+1):
        c = _qe_matvec(log_matrices, d)
        lower, upper = max(0, quota-degree*(updates-t)), min(quota, degree*t)
        d = _qe_project(c, log_x, lower, upper)
    return (d-quota*log_x).min(axis=1)-log_C


def _qe_unsplit(log_matrices, log_x, log_phi, updates, quota, log_C):
    d = np.repeat(np.asarray(log_phi, dtype=QE_LD)[:, None], len(log_x), axis=1)
    for _ in range(updates):
        d = _qe_matvec(log_matrices, d)
    return (d-quota*log_x).min(axis=1)-log_C


def _qe_production_matrices():
    xs = np.asarray(QE_X, dtype=QE_LD)
    matrices = np.full((len(xs), 3, 3), QE_LD('-inf'))
    for xi, x in enumerate(xs):
        values = [sum(QE_LD(coefficient)*x**j for j, coefficient in enumerate(raw))/255
                  for raw in QE_RAW]
        for poly, source, destination in QE_EDGES:
            matrices[xi, source, destination] = 128*np.log2(values[poly])
    return matrices, np.log2(xs)


def three_state_quota_control(half_updates=12, progress=True):
    half_updates = index(half_updates)
    if not 0 <= half_updates <= 12:
        raise ValueError('Use zero through12 explicitly bounded half-updates')
    started = time.monotonic()
    matrices, log_x = _qe_production_matrices()
    phi = np.asarray((0, 315, 687), dtype=QE_LD)
    log_C = _qe_logchoose(8192, 4096)
    best, summaries = None, []
    for iteration in range(half_updates+1):
        upper = quota_envelope_action(matrices, log_x, phi, 512, 16, 4096, log_C)
        unsplit = _qe_unsplit(matrices, log_x, phi, 16, 4096, log_C)
        ratios, old_ratios = upper-phi, unsplit-phi
        row, lower_row = int(np.argmax(ratios)), int(np.argmin(ratios))
        rho = ratios[row]
        restricted_tail = log_C-60518+phi[2]-phi.min()+256*rho
        summary = dict(iteration=iteration, log2_quota_envelope_growth_upper=float(rho),
                       maximum_row=('O', 'P', 'H')[row],
                       log2_envelope_comparison_floor=float(ratios[lower_row]),
                       minimum_row=('O', 'P', 'H')[lower_row],
                       log2_unsplit_growth_upper=float(old_ratios.max()),
                       log2_unsplit_comparison_floor=float(old_ratios.min()),
                       log2_phi_relative_H=[float(v-phi[2]) for v in phi],
                       log2_restricted_moment_tail_upper=float(restricted_tail),
                       restricted_operator_upper_only=True,
                       full_physical_tail_upper=False, floating_diagnostic_only=True)
        summaries.append(summary)
        if best is None or rho < best[0]:
            best = (rho, summary)
        if progress:
            print(dict(stage='three_state_quota_envelope', **summary), flush=True)
        phi = (phi+upper)/2
        phi -= phi.max()
    return dict(stage='three_state_quota_complete', seconds=round(time.monotonic()-started, 3),
                witnesses=len(summaries), quota_grid=list(QE_X), best=best[1],
                no_full_physical_operator_bound=True, full_code_certificate=False)


def _qe_convolve(a, b):
    out = [Fraction(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def _qe_restricted_polynomials(width):
    matrix = [[[Fraction(0)]*(width+1) for _ in range(3)] for _ in range(3)]
    for poly, source, destination in QE_EDGES:
        power = [Fraction(1)]
        for _ in range(width//4):
            power = _qe_convolve(power, [Fraction(c, 255) for c in QE_RAW[poly]])
        matrix[source][destination] = power
    return matrix


def _qe_exact_step(kernel, values):
    size, degree = len(kernel), len(kernel[0][0])-1
    after = [[Fraction(0)]*(len(values[0])+degree) for _ in range(size)]
    for i in range(size):
        for j in range(size):
            for r, v in enumerate(values[j]):
                for k, coefficient in enumerate(kernel[i][j]):
                    after[i][r+k] += coefficient*v
    return after


def three_state_quota_checks():
    counts = dict(coefficients=0, pruned_coefficients=0, exact_integer_envelope=0,
                  continuous_envelope=0, projection_identity=0, scaling=0, damping=0)
    generic = [[[Fraction(1+(i+1)*(j+2)+k*k, 2+i+j+k) for k in range(3)]
                for j in range(2)] for i in range(2)]
    cases = ((generic, 3, 3), (_qe_restricted_polynomials(4), 4, 8),
             (_qe_restricted_polynomials(8), 4, 16))
    xs = (Fraction(1), Fraction(2), Fraction(4))
    log_x = np.asarray((0, 1, 2), dtype=QE_LD)
    for kernel, updates, quota in cases:
        size, degree = len(kernel), len(kernel[0][0])-1
        evaluated = [[[sum(c*x**k for k, c in enumerate(poly)) for poly in row]
                      for row in kernel] for x in xs]
        logs = np.full((len(xs), size, size), QE_LD('-inf'))
        for l, matrix in enumerate(evaluated):
            for i in range(size):
                for j in range(size):
                    if matrix[i][j]:
                        logs[l, i, j] = _qe_logs([matrix[i][j]])[0]
        for phi in ([Fraction(i+2, i+1) for i in range(size)],
                    [Fraction(2)**(3*i) for i in range(size)]):
            full = [[p] for p in phi]
            pruned = [[p] for p in phi]
            d = [[p]*len(xs) for p in phi]
            for t in range(1, updates+1):
                full, pruned = _qe_exact_step(kernel, full), _qe_exact_step(kernel, pruned)
                lower, upper = max(0, quota-degree*(updates-t)), min(quota, degree*t)
                for i in range(size):
                    for r in range(len(pruned[i])):
                        if lower <= r <= upper:
                            assert pruned[i][r] == full[i][r]
                            counts['pruned_coefficients'] += 1
                        else:
                            pruned[i][r] = Fraction(0)
                c = [[sum(evaluated[l][i][j]*d[j][l] for j in range(size))
                      for l in range(len(xs))] for i in range(size)]
                old_logs = np.asarray([_qe_logs(row) for row in c])
                continuous = _qe_project(old_logs, log_x, lower, upper)
                d = [[max(min(c[i][j]*(xs[l]/xs[j])**r for j in range(len(xs)))
                          for r in range(lower, upper+1)) for l in range(len(xs))]
                     for i in range(size)]
                for i in range(size):
                    for r in range(lower, upper+1):
                        exact_bound = min(d[i][l]*xs[l]**(-r) for l in range(len(xs)))
                        assert full[i][r] <= exact_bound
                        counts['exact_integer_envelope'] += 1
                        before = min(c[i][l]*xs[l]**(-r) for l in range(len(xs)))
                        after_log = (continuous[i]-r*log_x).min()
                        assert abs(after_log-_qe_logs([before])[0]) < QE_LD('1e-10')
                        counts['projection_identity'] += 1
            C = Fraction(comb(degree*updates, quota))
            exact = [row[quota]/C for row in full]
            final_integer = [min(d[i][l]*xs[l]**(-quota) for l in range(len(xs)))/C
                             for i in range(size)]
            found = quota_envelope_action(logs, log_x, _qe_logs(phi), degree, updates, quota,
                                           _qe_logs([C])[0])
            assert np.all(found >= _qe_logs(exact)-QE_LD('1e-10'))
            assert np.all(found >= _qe_logs(final_integer)-QE_LD('1e-10'))
            counts['coefficients'] += size
            counts['continuous_envelope'] += size
            scaled = quota_envelope_action(logs, log_x, _qe_logs([p*Fraction(7, 3) for p in phi]),
                                            degree, updates, quota, _qe_logs([C])[0])
            assert np.max(np.abs(scaled-found-_qe_logs([Fraction(7, 3)])[0])) < QE_LD('1e-10')
            counts['scaling'] += size
            damped_phi = (_qe_logs(phi)+found)/2
            damped = quota_envelope_action(logs, log_x, damped_phi, degree, updates, quota,
                                            _qe_logs([C])[0])
            assert (damped-damped_phi).max() <= (found-_qe_logs(phi)).max()+QE_LD('1e-10')
            counts['damping'] += size
    print(dict(**counts, generic_two_state=True, restricted_widths=(4, 8),
               no_production_operator_test=True, full_code_certificate=False), flush=True)


def three_state_quota_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--half-updates', type=int, default=12)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            three_state_quota_checks()
        else:
            print(three_state_quota_control(args.half_updates), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    three_state_quota_main()
