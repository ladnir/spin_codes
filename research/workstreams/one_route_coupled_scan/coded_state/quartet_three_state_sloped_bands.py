"""Tiny sloped quota bands with an exact capped-mass support function.

Old bands bound f(u,s)<=C(M,u)*min_z d_Bz(s)*z^u. For fixed target r,
old interval J, slope z and evaluation x, use weights
 w_u=C(M,u)*(x*z)^u, u in J intersect[r-n,r].
Their support under0<=p_u<=eta_x,sum p<=1 is eta times the largest
floor(1/eta) weights plus the residual times the next. If fewer weights
exist, cap every available weight. The mode uses x*z; the atom bound
is eta_x, NOT eta_(xz). Multiply by x^-r K_x d_Bz, minimize over(z,x)
before summing bands, then normalize by C(nt,r). New intercepts are
max_r W(r,s)*z^-r; their reconstructed envelope may inflate W.

Only the restricted three-state kernel is used. Defaultwidth256,
z=1/2,1,2,seven fixedx,one witness,gauge logphi=(0,142,645).
No full physical operator claim. Ordinary FLOAT is not certification.
One optional dictionary diagnostic reports Hilbert-distance distortion
without replacing any action. No jobs or large arrays on import.
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

for _qs_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_qs_thread_var] = '1'

import numpy as np

QS_LD = np.longdouble
QS_X = (Fraction(1), Fraction(8, 5), Fraction(93, 50), Fraction(51, 25),
        Fraction(49, 20), Fraction(69, 25), Fraction(4))
QS_Z = (Fraction(1, 2), Fraction(1), Fraction(2))


def _qs_dependency(name):
    if name in globals():
        return globals()[name]
    if name.startswith('_qb_'):
        import quartet_three_state_quota_bands as module
    elif name.startswith('quota_anticoncentration'):
        import quartet_quota_anticoncentration as module
    else:
        import quartet_three_state_quota_envelope as module
    return getattr(module, name)


def _qs_logfraction(value):
    return np.log2(QS_LD(value.numerator)/value.denominator)


def _qs_psi_logs(logchoose, total, lower, upper, tilt, eta):
    """Vectorized top-weight search around the clipped exact binomial mode."""
    tilt, eta = Fraction(tilt), Fraction(eta)
    lower, upper = np.asarray(lower, dtype=np.int64), np.asarray(upper, dtype=np.int64)
    if (tilt <= 0 or not 0 < eta <= 1 or lower.shape != upper.shape
            or np.any(lower < 0) or np.any(upper > total) or np.any(lower > upper)):
        raise ValueError('Invalid capped binomial-weight intervals')
    log_tilt, log_eta = _qs_logfraction(tilt), _qs_logfraction(eta)
    q = eta.denominator//eta.numerator
    residual = 1-q*eta
    mode = ((total+1)*tilt.numerator)//(tilt.numerator+tilt.denominator)
    mode = np.minimum(np.maximum(mode, lower), upper)
    result = logchoose[mode]+mode*log_tilt+log_eta
    left, right = mode-1, mode+1
    limit = min(int((upper-lower+1).max()), q+int(residual > 0))
    for rank in range(2, limit+1):
        li, ri = np.maximum(left, lower), np.minimum(right, upper)
        lv = np.where(left >= lower, logchoose[li]+li*log_tilt, QS_LD('-inf'))
        rv = np.where(right <= upper, logchoose[ri]+ri*log_tilt, QS_LD('-inf'))
        take_left = (lv >= rv) & np.isfinite(lv)
        take_right = (~take_left) & np.isfinite(rv)
        value = np.maximum(lv, rv)
        valid = np.isfinite(value)
        left -= take_left.astype(np.int64)
        right += take_right.astype(np.int64)
        factor = log_eta if rank <= q else _qs_logfraction(residual)
        incoming = value[valid]+factor
        old = result[valid]
        peak = np.maximum(old, incoming)
        small = np.exp2(np.minimum(old, incoming)-peak)
        if np.any(small <= 0):
            raise ArithmeticError('A positive capped-mass term underflowed')
        result[valid] = peak+np.log2(1+small)
    return result


def _qs_dictionary(log_inputs, limit):
    """Diagnostic only: one representative per input, no K-through-min."""
    rows = log_inputs.reshape(-1, log_inputs.shape[-1])
    chosen = [0]
    differences = rows-rows[0]
    nearest = differences.max(axis=1)-differences.min(axis=1)
    while len(chosen) < min(limit, len(rows)):
        candidate = int(np.argmax(nearest))
        if nearest[candidate] == 0:
            break
        chosen.append(candidate)
        differences = rows-rows[candidate]
        distances = differences.max(axis=1)-differences.min(axis=1)
        np.minimum(nearest, distances, out=nearest)
    return len(chosen), nearest.max()


def sloped_band_action(log_matrices, xs, zs, eta, log_phi, degree, updates, quota,
                       band_width, logchoose=None, deadline=None, progress=False,
                       dictionary_size=0):
    degree, updates, quota, band_width = map(index, (degree, updates, quota, band_width))
    xs, zs, eta = tuple(map(Fraction, xs)), tuple(map(Fraction, zs)), tuple(map(Fraction, eta))
    phi = np.asarray(log_phi, dtype=QS_LD)
    if (degree < 1 or updates < 1 or not 0 <= quota <= degree*updates or band_width < 1
            or not xs or not zs or len(eta) != len(xs) or min(xs+zs) <= 0
            or any(not 0 < e <= 1 for e in eta) or not np.all(np.isfinite(phi))
            or log_matrices.shape != (len(xs), len(phi), len(phi))):
        raise ValueError('Invalid sloped-band dimensions, rational markers, or input')
    lx, lz = np.asarray([_qs_logfraction(x) for x in xs]), np.asarray([_qs_logfraction(z) for z in zs])
    logchoose = _qs_dependency('_qb_logchoose_rows')(degree, updates) if logchoose is None else logchoose
    bands = [(0, 0)]
    intercepts = np.repeat(phi[None, None, :], len(zs), axis=1)
    action_count = 0
    projection_sum = projection_max = QS_LD(0)
    dictionary_actions = 0
    dictionary_distortion_sum = QS_LD(0)
    projection_location = None
    started = time.monotonic()
    for t in range(1, updates+1):
        if deadline is not None and time.monotonic() >= deadline:
            raise TimeoutError('Sloped-band witness incomplete at deadline')
        lo_target, hi_target = _qs_dependency('_qb_interval')(degree, updates, quota, t)
        targets = np.arange(lo_target, hi_target+1, dtype=np.int64)
        accum = np.full((len(targets), len(phi)), QS_LD('-inf'))
        terms = log_matrices[None, None, :, :, :]+intercepts[:, :, None, None, :]
        peaks = terms.max(axis=4)
        weights = np.exp2(terms-peaks[:, :, :, :, None])
        if np.any(np.isfinite(terms) & (weights <= 0)) or not np.all(np.isfinite(peaks)):
            raise ArithmeticError('Positive sloped-band action underflow or empty row')
        actions = peaks+np.log2(weights.sum(axis=4, dtype=QS_LD))
        action_count += len(bands)*len(zs)*len(xs)
        if dictionary_size:
            representatives, distortion = _qs_dictionary(intercepts, dictionary_size)
            dictionary_actions += representatives*len(xs)
            dictionary_distortion_sum += distortion
        total = degree*(t-1)
        for bi, (old_lo, old_hi) in enumerate(bands):
            if bi % 8 == 0 and deadline is not None and time.monotonic() >= deadline:
                raise TimeoutError('Sloped-band witness incomplete at deadline')
            rlo, rhi = max(lo_target, old_lo), min(hi_target, old_hi+degree)
            if rlo > rhi:
                continue
            selection = slice(rlo-lo_target, rhi-lo_target+1)
            r = targets[selection]
            lower, upper = np.maximum(old_lo, r-degree), np.minimum(old_hi, r)
            candidate = np.full((len(r), len(phi)), QS_LD('inf'))
            for zi, z in enumerate(zs):
                for xi, x in enumerate(xs):
                    psi = _qs_psi_logs(logchoose[t-1], total, lower, upper, x*z, eta[xi])
                    scalar = psi-r*lx[xi]
                    np.minimum(candidate, scalar[:, None]+actions[bi, zi, xi], out=candidate)
            _qs_dependency('_qb_accumulate')(accum[selection], candidate)
        if not np.all(np.isfinite(accum)):
            raise ArithmeticError('A feasible coefficient has no finite positive bound')
        accum -= logchoose[t][targets, None]
        new_bands = _qs_dependency('_qb_bands')(lo_target, hi_target, band_width)
        new_values = np.empty((len(new_bands), len(zs), len(phi)), dtype=QS_LD)
        step_projection = QS_LD(0)
        for bi, (rlo, rhi) in enumerate(new_bands):
            selection = slice(rlo-lo_target, rhi-lo_target+1)
            r = targets[selection]
            degreewise = accum[selection]
            new_values[bi] = (degreewise[:, None, :]-r[:, None, None]*lz[None, :, None]).max(axis=0)
            envelope = (new_values[bi][None, :, :]+r[:, None, None]*lz[None, :, None]).min(axis=1)
            inflation = envelope-degreewise
            if inflation.min() < QS_LD('-1e-9'):
                raise ArithmeticError('New intercepts failed to majorize their degreewise input')
            value = inflation.max()
            step_projection = max(step_projection, value)
            if value > projection_max:
                rr, state = np.unravel_index(np.argmax(inflation), inflation.shape)
                projection_max, projection_location = value, (t, int(r[rr]), int(state))
        projection_sum += step_projection
        bands, intercepts = new_bands, new_values
        if progress:
            print(dict(stage='sloped_band_step', step=t, bands=len(bands),
                       projection_max_bits=float(step_projection), cumulative_actions=action_count,
                       seconds=round(time.monotonic()-started, 3)), flush=True)
    assert bands == [(quota, quota)]
    upper = (intercepts[0]+quota*lz[:, None]).min(axis=0)
    return dict(log2_upper=upper, action_count=action_count,
                projection_max_bits=projection_max, projection_sum_max_bits=projection_sum,
                projection_location=projection_location,
                dictionary_actions=dictionary_actions if dictionary_size else None,
                dictionary_sum_distortion_bits=dictionary_distortion_sum if dictionary_size else None,
                seconds=time.monotonic()-started)


def three_state_sloped_control(band_width=256, gauge='passing', dictionary_size=0, deadline=None):
    band_width, dictionary_size = index(band_width), index(dictionary_size)
    if not 1 <= band_width <= 512 or not 0 <= dictionary_size <= 8:
        raise ValueError('Use bandwidth1..512 and optional dictionary size0..8')
    matrices, _ = _qs_dependency('_qe_production_matrices')()
    phi = np.asarray((0, 142, 645) if gauge == 'passing' else (0, 315, 687), dtype=QS_LD)
    eta = tuple(_qs_dependency('quota_anticoncentration_eta')(x, quartets=128) for x in QS_X)
    result = sloped_band_action(matrices, QS_X, QS_Z, eta, phi, 512, 16, 4096, band_width,
                                deadline=deadline, progress=True, dictionary_size=dictionary_size)
    ratios = result['log2_upper']-phi
    rho = ratios.max()
    log_C = _qs_dependency('_qe_logchoose')(8192, 4096)
    return dict(stage='three_state_sloped_complete', seconds=round(result['seconds'], 3),
                band_width=band_width, quota_slopes=[str(z) for z in QS_Z], gauge=gauge,
                log2_growth_upper=float(rho), log2_comparison_floor=float(ratios.min()),
                maximum_row=('O', 'P', 'H')[int(np.argmax(ratios))],
                minimum_row=('O', 'P', 'H')[int(np.argmin(ratios))],
                log2_restricted_moment_tail_upper=float(log_C-60518+phi[2]-phi.min()+256*rho),
                action_count=result['action_count'], projection_max_bits=float(result['projection_max_bits']),
                projection_sum_max_bits=float(result['projection_sum_max_bits']),
                projection_location=result['projection_location'],
                dictionary_actions=result['dictionary_actions'],
                dictionary_sum_distortion_bits=(float(result['dictionary_sum_distortion_bits'])
                                               if dictionary_size else None),
                dictionary_diagnostic_only=True, restricted_operator_upper_only=True,
                full_physical_tail_upper=False, floating_diagnostic_only=True, full_code_certificate=False)


def _qs_exact_psi(weights, eta):
    ordered = sorted(weights, reverse=True)
    q = eta.denominator//eta.numerator
    value = eta*sum(ordered[:q])
    if q < len(ordered):
        value += (1-q*eta)*ordered[q]
    return value


def _qs_vertices(weights, eta):
    q = eta.denominator//eta.numerator
    if q >= len(weights):
        return eta*sum(weights)
    best = Fraction(0)
    for full in combinations(range(len(weights)), q):
        base = eta*sum(weights[i] for i in full)
        for other in range(len(weights)):
            if other not in full:
                best = max(best, base+(1-q*eta)*weights[other])
    return best


def _qs_reference(kernel, xs, zs, eta, phi, updates, quota, band_width, counts):
    states, degree = len(phi), len(kernel[0][0])-1
    evaluated = [[[sum(c*x**k for k, c in enumerate(poly)) for poly in row] for row in kernel] for x in xs]
    bands, ds = [(0, 0)], [[phi[:] for _ in zs]]
    full = [[p] for p in phi]
    for t in range(1, updates+1):
        full = _qs_dependency('_qe_exact_step')(kernel, full)
        lower, upper = _qs_dependency('_qb_interval')(degree, updates, quota, t)
        total = degree*(t-1)
        W = []
        for r in range(lower, upper+1):
            bound = [Fraction(0)]*states
            for bi, (old_lo, old_hi) in enumerate(bands):
                lo, hi = max(old_lo, r-degree), min(old_hi, r)
                if lo > hi:
                    continue
                candidates = []
                for zi, z in enumerate(zs):
                    for xi, x in enumerate(xs):
                        psi = _qs_exact_psi([comb(total, u)*(x*z)**u for u in range(lo, hi+1)], eta[xi])
                        candidates.append([x**(-r)*psi*sum(evaluated[xi][s][j]*ds[bi][zi][j]
                                                          for j in range(states)) for s in range(states)])
                bound = [v+min(row[s] for row in candidates) for s, v in enumerate(bound)]
            bound = [v/comb(degree*t, r) for v in bound]
            assert all(full[s][r]/comb(degree*t, r) <= bound[s] for s in range(states))
            counts['coefficient_bounds'] += states
            W.append(bound)
        bands = _qs_dependency('_qb_bands')(lower, upper, band_width)
        ds = [[[max(W[r-lower][s]*z**(-r) for r in range(lo, hi+1)) for s in range(states)]
               for z in zs] for lo, hi in bands]
        for bi, (lo, hi) in enumerate(bands):
            for r in range(lo, hi+1):
                assert all(W[r-lower][s] <= min(ds[bi][zi][s]*z**r for zi, z in enumerate(zs))
                           for s in range(states))
                counts['projection_bounds'] += states
    final = [min(ds[0][zi][s]*z**quota for zi, z in enumerate(zs)) for s in range(states)]
    assert all(full[s][quota]/comb(degree*updates, quota) <= final[s] for s in range(states))
    counts['normalization'] += states
    return final, evaluated


def three_state_sloped_checks():
    counts = dict(psi_extreme_vertices=0, top_weight_logs=0, coefficient_bounds=0,
                  projection_bounds=0, normalization=0, floating_rows=0, scale_rows=0)
    logs = _qs_dependency('_qe_logs')
    for total in range(7):
        intervals = [(lo, hi) for lo in range(total+1) for hi in range(lo, min(total, lo+4)+1)]
        lo = np.asarray([a for a, b in intervals]); hi = np.asarray([b for a, b in intervals])
        logchoose = logs([Fraction(comb(total, u)) for u in range(total+1)])
        for tilt in (Fraction(1, 2), Fraction(1), Fraction(2)):
            for eta in (Fraction(1), Fraction(1, 2), Fraction(2, 5), Fraction(1, 4)):
                exact = []
                for left, right in intervals:
                    weights = [Fraction(comb(total, u))*tilt**u for u in range(left, right+1)]
                    value = _qs_exact_psi(weights, eta)
                    assert value == _qs_vertices(weights, eta)
                    exact.append(value)
                    counts['psi_extreme_vertices'] += 1
                actual = _qs_psi_logs(logchoose, total, lo, hi, tilt, eta)
                assert np.max(np.abs(actual-logs(exact))) < QS_LD('1e-10')
                counts['top_weight_logs'] += len(intervals)
    xs, zs = (Fraction(1), Fraction(2), Fraction(4)), QS_Z
    generic = [[[Fraction(1+(i+1)*(j+2)+k*k, 2+i+j+k) for k in range(3)] for j in range(2)] for i in range(2)]
    cases = ((generic, (Fraction(1),)*3, 3, 3, 2),)
    for width in (8, 16):
        kernel = _qs_dependency('_qe_restricted_polynomials')(width)
        eta = tuple(_qs_dependency('quota_anticoncentration_eta')(x, quartets=width//4) for x in xs)
        cases += ((kernel, eta, 2, width, 3),)
    for kernel, eta, updates, quota, band_width in cases:
        size, degree = len(kernel), len(kernel[0][0])-1
        phi = [Fraction(2)**(2*i) for i in range(size)]
        exact, evaluated = _qs_reference(kernel, xs, zs, eta, phi, updates, quota, band_width, counts)
        matrices = np.full((len(xs), size, size), QS_LD('-inf'))
        for xi, matrix in enumerate(evaluated):
            for i in range(size):
                for j in range(size):
                    if matrix[i][j]:
                        matrices[xi, i, j] = logs([matrix[i][j]])[0]
        found = sloped_band_action(matrices, xs, zs, eta, logs(phi), degree, updates, quota, band_width)
        assert np.max(np.abs(found['log2_upper']-logs(exact))) < QS_LD('1e-10')
        counts['floating_rows'] += size
        scaled = sloped_band_action(matrices, xs, zs, eta, logs([v*Fraction(7, 3) for v in phi]),
                                    degree, updates, quota, band_width)
        assert np.max(np.abs(scaled['log2_upper']-found['log2_upper']-logs([Fraction(7, 3)])[0])) < QS_LD('1e-10')
        counts['scale_rows'] += size
    print(dict(**counts, actual_restricted_widths=(8, 16), no_full_physical_operator_test=True,
               full_code_certificate=False), flush=True)


def three_state_sloped_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--control', action='store_true')
    parser.add_argument('--band-width', type=int, default=256)
    parser.add_argument('--gauge', choices=('passing', 'original'), default='passing')
    parser.add_argument('--dictionary-size', type=int, default=0)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    started = time.monotonic()
    timer = threading.Timer(85, lambda: os._exit(124)); timer.daemon = True; timer.start()
    try:
        if args.checks:
            three_state_sloped_checks()
        else:
            print(three_state_sloped_control(args.band_width, args.gauge, args.dictionary_size, started+80), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    three_state_sloped_main()
