"""Tiny binomial-reference quota bands for the restricted three-state kernel.

Represent coefficient f_t(s,r)<=C(nt,r)*v_(t,I)(s) on disjoint bands.
For each old band J and target r, intersect J with[r-n,r]. At rational
x=A/B the maximum C(n(t-1),u)*x^u is at the clipped integer mode
floor((n(t-1)+1)*A/(A+B)). Multiply that scalar by L_x*v_J, minimize
over x BEFORE summing old bands, divide by C(nt,r), then maximize r
within each new band. At the final singleton r=h, normalization is
already complete; do NOT divide by C(nT,h) again.

Optional anti-concentration multiplies each old-band/x contribution by
min(1,(hi-lo+1)*eta_x) BEFORE the x minimum. eta_x bounds every atom
of each evaluated transition's label-degree law. It is an entrywise
evaluated-operator inequality, not a polynomial coefficientwise factor.
The default has no such factor and preserves the original comparison.

Production remains only O,P,H: width512,16 updates, seven fixed x,
bandwidth8, initial potential(1,2^315,2^687), one witness by default.
At most3 explicitly requested half-updates. Every target-r operation
is vectorized within one old band's reachable slice; no all-band cube.
This is a FLOAT UPPER comparison for the restricted quota operator,
not a lower bound or upper bound for the full physical operator.
The displayed B_restricted^256*1 expression has only that scope.
No full population action, files, import jobs, or construction change.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time

for _qb_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_qb_thread_var] = '1'

import numpy as np

QB_LD = np.longdouble
QB_X = (Fraction(1), Fraction(8, 5), Fraction(93, 50), Fraction(51, 25),
        Fraction(49, 20), Fraction(69, 25), Fraction(4))


def _qb_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('quota_anticoncentration'):
        import quartet_quota_anticoncentration as module
    else:
        import quartet_three_state_quota_envelope as module
    return getattr(module, name)


def _qb_interval(degree, updates, quota, t):
    return max(0, quota-degree*(updates-t)), min(quota, degree*t)


def _qb_bands(lower, upper, width):
    return [(max(lower, k*width), min(upper, (k+1)*width-1))
            for k in range(lower//width, upper//width+1)]


def _qb_mode(total, numerator, denominator, lower, upper):
    if not 0 <= lower <= upper <= total or min(numerator, denominator) <= 0:
        raise ValueError('Invalid binomial-mode interval or positive rational evaluation')
    return max(lower, min(upper, ((total+1)*numerator)//(numerator+denominator)))


def _qb_logchoose_rows(degree, updates):
    rows = []
    for t in range(updates+1):
        total = degree*t
        row = np.zeros(total+1, dtype=QB_LD)
        half = total//2
        r = np.arange(1, half+1, dtype=QB_LD)
        row[1:half+1] = np.cumsum(np.log2((total-r+1)/r), dtype=QB_LD)
        # Enforce the exact symmetry rather than summing negative increments
        # after the middle and accumulating cancellation near the endpoint.
        row[total-np.arange(half+1)] = row[:half+1]
        rows.append(row)
    return rows


def _qb_actions(log_matrices, old_values):
    """One positive action per(oldband,x), all three sources retained."""
    terms = log_matrices[None, :, :, :]+old_values[:, None, None, :]
    peaks = terms.max(axis=3)
    if not np.all(np.isfinite(peaks)):
        raise ArithmeticError('An old-band kernel row has no positive mass')
    weights = np.exp2(terms-peaks[:, :, :, None])
    if np.any(np.isfinite(terms) & (weights <= 0)):
        raise ArithmeticError('A positive band action term underflowed')
    return peaks+np.log2(weights.sum(axis=3, dtype=QB_LD))


def _qb_accumulate(destination, incoming):
    peak = np.maximum(destination, incoming)
    small = np.exp2(np.minimum(destination, incoming)-peak)
    if np.any(np.isfinite(destination) & (small <= 0)):
        raise ArithmeticError('A positive old-band contribution underflowed')
    destination[:] = peak+np.log2(1+small)


def quota_band_action(log_matrices, xs, log_phi, degree, updates, quota, band_width,
                      logchoose=None, deadline=None, progress=False, eta=None):
    degree, updates, quota, band_width = map(index, (degree, updates, quota, band_width))
    if degree < 1 or updates < 1 or not 0 <= quota <= degree*updates or band_width < 1:
        raise ValueError('Invalid degree, update count, quota, or band width')
    xs = tuple(Fraction(x) for x in xs)
    if not xs or any(x <= 0 for x in xs) or len(set(xs)) != len(xs):
        raise ValueError('Use distinct positive rational evaluations')
    log_x = np.asarray([np.log2(QB_LD(x.numerator)/x.denominator) for x in xs], dtype=QB_LD)
    if eta is None:
        log_eta = None
    else:
        eta = tuple(Fraction(v) for v in eta)
        if len(eta) != len(xs) or any(not 0 < v <= 1 for v in eta):
            raise ValueError('Use one positive anti-concentration bound at most1 per evaluation')
        log_eta = np.asarray([np.log2(QB_LD(v.numerator)/v.denominator) for v in eta], dtype=QB_LD)
    phi = np.asarray(log_phi, dtype=QB_LD)
    if not np.all(np.isfinite(phi)) or log_matrices.shape != (len(xs), len(phi), len(phi)):
        raise ValueError('Invalid positive physical vector or evaluated kernel shape')
    logchoose = _qb_logchoose_rows(degree, updates) if logchoose is None else logchoose
    old_bands, old_values = [(0, 0)], phi[None, :].copy()
    started, max_bands = time.monotonic(), 1
    for t in range(1, updates+1):
        if deadline is not None and time.monotonic() >= deadline:
            raise TimeoutError('Quota-band witness incomplete at its bounded deadline')
        lower, upper = _qb_interval(degree, updates, quota, t)
        target_r = np.arange(lower, upper+1, dtype=np.int64)
        accum = np.full((len(target_r), len(phi)), QB_LD('-inf'))
        actions = _qb_actions(log_matrices, old_values)
        total = degree*(t-1)
        modes = [((total+1)*x.numerator)//(x.numerator+x.denominator) for x in xs]
        for bi, (old_lower, old_upper) in enumerate(old_bands):
            if bi % 32 == 0 and deadline is not None and time.monotonic() >= deadline:
                raise TimeoutError('Quota-band witness incomplete at its bounded deadline')
            rlo, rhi = max(lower, old_lower), min(upper, old_upper+degree)
            if rlo > rhi:
                continue
            sl = slice(rlo-lower, rhi-lower+1)
            r = target_r[sl]
            lo, hi = np.maximum(old_lower, r-degree), np.minimum(old_upper, r)
            length_logs = np.log2((hi-lo+1).astype(QB_LD)) if log_eta is not None else None
            candidate = np.full((len(r), len(phi)), QB_LD('inf'))
            for l, mode in enumerate(modes):
                u = np.minimum(np.maximum(mode, lo), hi)
                scalar = logchoose[t-1][u]+(u-r)*log_x[l]
                if log_eta is not None:
                    scalar += np.minimum(QB_LD(0), length_logs+log_eta[l])
                np.minimum(candidate, scalar[:, None]+actions[bi, l], out=candidate)
            _qb_accumulate(accum[sl], candidate)
        if not np.all(np.isfinite(accum)):
            raise ArithmeticError('A feasible target coefficient has no finite positive upper')
        accum -= logchoose[t][target_r, None]
        new_bands = _qb_bands(lower, upper, band_width)
        values = np.empty((len(new_bands), len(phi)), dtype=QB_LD)
        for bi, (rlo, rhi) in enumerate(new_bands):
            values[bi] = accum[rlo-lower:rhi-lower+1].max(axis=0)
        old_bands, old_values = new_bands, values
        max_bands = max(max_bands, len(new_bands))
        if progress and (t % 4 == 0 or t == updates):
            print(dict(stage='three_state_quota_bands_step', step=t, bands=len(new_bands),
                       feasible_degrees=len(target_r), seconds=round(time.monotonic()-started, 3)), flush=True)
    assert old_bands == [(quota, quota)]
    return dict(log2_upper=old_values[0].copy(), max_bands=max_bands,
                seconds=time.monotonic()-started, final_binomial_division_already_included=True)


def three_state_quota_band_control(band_width=8, half_updates=0, deadline=None,
                                   anticoncentration=False):
    band_width, half_updates = index(band_width), index(half_updates)
    if not 1 <= band_width <= 512 or not 0 <= half_updates <= 3:
        raise ValueError('Use bandwidth1..512 and at most3 requested half-updates')
    matrices, _ = _qb_dependency('_qe_production_matrices')()
    logchoose = _qb_logchoose_rows(512, 16)
    phi = np.asarray((0, 315, 687), dtype=QB_LD)
    log_C = logchoose[16][4096]
    eta = None
    if anticoncentration:
        eta = tuple(_qb_dependency('quota_anticoncentration_eta')(x, quartets=128) for x in QB_X)
        print(dict(stage='three_state_quota_band_eta', x_values=[str(x) for x in QB_X],
                   eta=[float(v) for v in eta], exact_rational_eta=True), flush=True)
    best, started = None, time.monotonic()
    for iteration in range(half_updates+1):
        result = quota_band_action(matrices, QB_X, phi, 512, 16, 4096, band_width,
                                   logchoose, deadline, True, eta)
        upper = result['log2_upper']
        ratios = upper-phi
        hi, lo = int(np.argmax(ratios)), int(np.argmin(ratios))
        rho = ratios[hi]
        restricted_tail = log_C-60518+phi[2]-phi.min()+256*rho
        summary = dict(iteration=iteration, seconds=round(result['seconds'], 3),
                       band_width=band_width, max_bands=result['max_bands'],
                       anticoncentration=bool(anticoncentration),
                       log2_growth_upper=float(rho), maximum_row=('O', 'P', 'H')[hi],
                       log2_band_comparison_floor=float(ratios[lo]), minimum_row=('O', 'P', 'H')[lo],
                       log2_phi_relative_H=[float(v-phi[2]) for v in phi],
                       log2_restricted_moment_tail_upper=float(restricted_tail),
                       restricted_operator_upper_only=True, full_physical_tail_upper=False,
                       floating_diagnostic_only=True)
        print(dict(stage='three_state_quota_bands_witness', **summary), flush=True)
        if best is None or rho < best[0]:
            best = (rho, summary)
        phi = (phi+upper)/2
        phi -= phi.max()
    return dict(stage='three_state_quota_bands_complete', seconds=round(time.monotonic()-started, 3),
                completed_witnesses=half_updates+1, best=best[1],
                full_physical_operator_not_bounded=True, full_code_certificate=False)


def _qb_exact_reference(kernel, xs, phi, updates, quota, band_width, counters, eta=None):
    size, degree = len(phi), len(kernel[0][0])-1
    matrices = [[[sum(c*x**k for k, c in enumerate(poly)) for poly in row]
                 for row in kernel] for x in xs]
    old_bands, old_values = [(0, 0)], [phi]
    full = [[p] for p in phi]
    for t in range(1, updates+1):
        full = _qb_dependency('_qe_exact_step')(kernel, full)
        lower, upper = _qb_interval(degree, updates, quota, t)
        total = degree*(t-1)
        values_by_r = []
        for r in range(lower, upper+1):
            bound = [Fraction(0)]*size
            for (old_lo, old_hi), vector in zip(old_bands, old_values):
                lo, hi = max(old_lo, r-degree), min(old_hi, r)
                if lo > hi:
                    continue
                candidates = []
                for l, (x, matrix) in enumerate(zip(xs, matrices)):
                    # Deliberate independent scan, not the clipped-mode formula.
                    multiplier = max(Fraction(comb(total, u))*x**u for u in range(lo, hi+1))*x**(-r)
                    if eta is not None:
                        multiplier *= min(Fraction(1), (hi-lo+1)*eta[l])
                    candidates.append([multiplier*sum(matrix[i][j]*vector[j] for j in range(size))
                                       for i in range(size)])
                for i in range(size):
                    bound[i] += min(row[i] for row in candidates)
            bound = [v/comb(degree*t, r) for v in bound]
            for i in range(size):
                assert full[i][r]/comb(degree*t, r) <= bound[i]
                counters['coefficient_band_inequality'] += 1
            values_by_r.append(bound)
        old_bands = _qb_bands(lower, upper, band_width)
        old_values = [[max(values_by_r[r-lower][i] for r in range(lo, hi+1)) for i in range(size)]
                      for lo, hi in old_bands]
    exact = [row[quota]/comb(degree*updates, quota) for row in full]
    assert all(v <= u for v, u in zip(exact, old_values[0]))
    counters['final_normalization'] += size
    return old_values[0], matrices


def three_state_quota_band_checks():
    counts = dict(mode_clipping=0, coefficient_band_inequality=0,
                  final_normalization=0, floating_rows=0, scale_rows=0, logchoose=0,
                  baseline_eta_one=0, refined_atoms=0, refined_float_rows=0)
    modes_x = (Fraction(1, 4), Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2), Fraction(4))
    for total in range(9):
        for lo in range(total+1):
            for hi in range(lo, total+1):
                for x in modes_x:
                    mode = _qb_mode(total, x.numerator, x.denominator, lo, hi)
                    assert comb(total, mode)*x**mode == max(comb(total, u)*x**u for u in range(lo, hi+1))
                    counts['mode_clipping'] += 1
    logrows = _qb_logchoose_rows(8, 3)
    logs = _qb_dependency('_qe_logs')
    for t, row in enumerate(logrows):
        for r, value in enumerate(row):
            assert abs(value-logs([Fraction(comb(8*t, r))])[0]) < QB_LD('1e-12')
            counts['logchoose'] += 1
    xs = (Fraction(1), Fraction(2), Fraction(4))
    for degree, updates, band_width in ((2, 3, 1), (4, 3, 2), (8, 2, 3)):
        quota = degree*updates//2
        kernel = [[[Fraction(1+(i+1)*(j+2)+k*k, 2+i+j+k) for k in range(degree+1)]
                   for j in range(2)] for i in range(2)]
        for phi in ([Fraction(1), Fraction(7, 3)], [Fraction(5, 2), Fraction(2)]):
            exact, matrices = _qb_exact_reference(kernel, xs, phi, updates, quota, band_width, counts)
            log_matrices = np.asarray([[logs(row) for row in matrix] for matrix in matrices])
            found = quota_band_action(log_matrices, xs, logs(phi), degree, updates, quota, band_width)
            assert np.max(np.abs(found['log2_upper']-logs(exact))) < QB_LD('1e-10')
            counts['floating_rows'] += len(phi)
            neutral = quota_band_action(log_matrices, xs, logs(phi), degree, updates, quota, band_width,
                                        eta=(Fraction(1),)*len(xs))
            assert np.array_equal(neutral['log2_upper'], found['log2_upper'])
            counts['baseline_eta_one'] += len(phi)
            scaled = quota_band_action(log_matrices, xs, logs([p*Fraction(9, 4) for p in phi]),
                                       degree, updates, quota, band_width)
            assert np.max(np.abs(scaled['log2_upper']-found['log2_upper']-logs([Fraction(9, 4)])[0])) < QB_LD('1e-10')
            counts['scale_rows'] += len(phi)
    # Actual restricted quartet kernels, not generic polynomials, are needed
    # for the refined eta assumption. The exact reference independently scans
    # every coefficient and old-band interval instead of using the mode rule.
    for degree, band_width in ((8, 1), (16, 2)):
        kernel = _qb_dependency('_qe_restricted_polynomials')(degree)
        eta = tuple(_qb_dependency('quota_anticoncentration_eta')(x, quartets=degree//4) for x in xs)
        phi = [Fraction(1), Fraction(8), Fraction(64)]
        exact, matrices = _qb_exact_reference(kernel, xs, phi, 2, degree, band_width, counts, eta)
        log_matrices = np.full((len(xs), 3, 3), QB_LD('-inf'))
        for l, (x, matrix) in enumerate(zip(xs, matrices)):
            for i in range(3):
                for j in range(3):
                    if matrix[i][j]:
                        log_matrices[l, i, j] = logs([matrix[i][j]])[0]
                    for k, coefficient in enumerate(kernel[i][j]):
                        assert coefficient*x**k <= eta[l]*matrix[i][j]
                        counts['refined_atoms'] += 1
        found = quota_band_action(log_matrices, xs, logs(phi), degree, 2, degree, band_width, eta=eta)
        assert np.max(np.abs(found['log2_upper']-logs(exact))) < QB_LD('1e-10')
        counts['refined_float_rows'] += 3
    print(dict(**counts, generic_positive_degrees=(2, 4, 8),
               no_full_physical_operator_test=True, full_code_certificate=False), flush=True)


def three_state_quota_band_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--band-width', type=int, default=8)
    parser.add_argument('--half-updates', type=int, default=0)
    parser.add_argument('--anticoncentration', action='store_true')
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    started = time.monotonic()
    timer = threading.Timer(85, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            three_state_quota_band_checks()
        else:
            print(three_state_quota_band_control(args.band_width, args.half_updates, started+80,
                                                 args.anticoncentration), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    three_state_quota_band_main()
