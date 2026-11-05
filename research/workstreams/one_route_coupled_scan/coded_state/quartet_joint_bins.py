"""Positive availability-bin comparison for the joint covered-hole transfer.

Load after quartet_joint_subset.py with its main block removed.  This file
streams its physical Horner rows and never allocates the physical matrix.
For fixed retained-hole count, positive available counts are grouped in
fixed intervals.  Available count zero remains a singleton in every group.

Let H_(s,I)=phi_s*1[s in I], phi_s=t**u_s.  With P=K-Z, define
Bplus_(I,J)=max_(s in I) (P H)_(s,J)/phi_s and
cp_I=max_(s in I) (terminal_positive)_s/phi_s.  Then K H<=H B,
where B=Zbin+Bplus, and terminal<=H(vbin+cp).  The singleton zero bins
give Z H=H Zbin and v=H vbin exactly.  Consequently the positive regional
and free-terminal comparisons compose.  Take minima across quota tilts
only AFTER regional composition and quota conditioning.  Bins and t are
fixed across all quota tilts.  All reported probes are FLOAT diagnostics.
"""

import argparse
from fractions import Fraction
from math import comb, lgamma
import os
import shutil
import threading
import time

import numpy as np


JB_LD = np.longdouble
JB_WIDTH = 512
JB_W = 943718


def _joint_bins_dependencies():
    needed = ('joint_subset_prepare', '_joint_horner_background', 'joint_subset_field_weights',
              'joint_subset_region', 'joint_subset_tail')
    if any(name not in globals() for name in needed):
        raise RuntimeError('Inject quartet_joint_subset.py definitions before using the binned helper')


def joint_bins_prepare(width, levels=(0, 1, 16), u_step=4, logpotential=0.):
    _joint_bins_dependencies()
    if not isinstance(u_step, int) or u_step < 1:
        raise ValueError('Availability bin width must be a positive integer')
    physical = joint_subset_prepare(width, levels)
    keys = sorted({(0 if u == 0 else 1 + (u - 1) // u_step, e) for u, e in physical['states']})
    key_index = {key: j for j, key in enumerate(keys)}
    source_bins = np.array([key_index[0 if u == 0 else 1 + (u - 1) // u_step, e]
                            for u, e in physical['states']], dtype=int)
    phi = np.exp2(JB_LD(logpotential) * np.arange(width + 1, dtype=int))
    if not np.all(np.isfinite(phi)) or not np.all(phi > 0):
        raise ValueError('Potential must be finite and strictly positive throughout the width')
    target_groups = {}
    descriptors = [None] * len(keys)
    for e in physical['levels']:
        _, sizes = physical['by_level'][e]
        last = int(sizes[-1])
        starts = np.array([0] + list(range(1, last + 1, u_step)), dtype=int)
        bins = np.array([key_index[0 if start == 0 else 1 + (int(start) - 1) // u_step, e]
                         for start in starts], dtype=int)
        for start, target in zip(starts, bins):
            finish = 0 if start == 0 else min(last, int(start) + u_step - 1)
            descriptors[target] = (int(start), finish, e)
        target_groups[e] = sizes, starts, bins
    zero_vector = np.array([JB_LD(key[0] == 0) for key in keys])
    return dict(width=width, levels=tuple(levels), u_step=u_step, logpotential=JB_LD(logpotential),
                physical=physical, keys=tuple(keys), descriptors=tuple(descriptors),
                source_bins=source_bins, phi=phi, target_groups=target_groups,
                zero_vector=zero_vector)


def _joint_bin_physical_rows(prepared, m, b, p, r, alpha, coverage):
    """Yield exact low-degree / saturated-tail row polynomials, without a matrix."""
    physical = prepared['physical']
    width = prepared['width']
    levels = prepared['levels']
    groups = width // 4
    maximum = levels[-1]
    m, b, p, r, alpha, coverage = map(JB_LD, (m, b, p, r, alpha, coverage))
    if min(m, b, p, r, alpha, coverage) < 0:
        raise ValueError('All physical weights must be nonnegative; alpha may exceed one')
    a = m * p * coverage
    exponents = np.arange(width + 1, dtype=int)
    apowers, bpowers, alphapowers = a ** exponents, b ** exponents, alpha ** exponents
    route_powers = (r * coverage) ** (4 * np.arange(groups + 1, dtype=int))
    scaled = {e: [poly * apowers[:len(poly)] for poly in physical['geometry'][e]] for e in levels}
    scaled_sums = {e: np.array([poly.sum() for poly in rows]) for e, rows in scaled.items()}
    old_polynomials = [physical['choose'][u] * bpowers[:u + 1] * alphapowers[u::-1]
                       for u in range(width + 1)]
    for source, (u, e) in enumerate(physical['states']):
        coefficients = physical['probabilities'][u] * route_powers / JB_LD(comb(width - u, e))
        background = _joint_horner_background(width, u, e, maximum, coefficients, scaled[e], b)
        if maximum == 0:
            outgoing = np.convolve(background[0], old_polynomials[u])[None, :]
        else:
            outgoing = np.zeros((maximum + 1, width + 1), dtype=JB_LD)
            for new_holes in range(maximum):
                for old_holes in range(min(u, new_holes) + 1):
                    shift = u - old_holes
                    factor = JB_LD(comb(u, old_holes)) * alphapowers[old_holes] * bpowers[shift]
                    outgoing[new_holes, shift:shift + width - u + 1] += factor * background[new_holes - old_holes]
            outgoing[maximum] = np.convolve(background[maximum], old_polynomials[u])
            for background_holes in range(maximum):
                needed = maximum - background_holes
                if needed <= u:
                    addition = np.convolve(background[background_holes], old_polynomials[u][:u - needed + 1])
                    outgoing[maximum, :len(addition)] += addition
        if u == 0:
            # The complete zero-coverage term is exactly one at destination (0,0).
            outgoing[0, 0] = 1
            terminal_positive = scaled[e][groups][1:].sum() / JB_LD(comb(width, e))
        else:
            terminal_positive = sum(coefficients[k] * bpowers[4 * k] * scaled_sums[e][groups - k]
                                    for k in range((u + 3) // 4, min(u, groups - (e + 3) // 4) + 1))
        yield source, u, e, outgoing, terminal_positive


def joint_bins_matrix(width, m, b, p, r, alpha, levels=(0, 1, 16), u_step=4,
                       logpotential=0., coverage=1., prepared=None, progress=False):
    """Stream max(P H)/phi and max(terminal_positive/phi), then add exact Zbin."""
    if prepared is None:
        prepared = joint_bins_prepare(width, levels, u_step, logpotential)
    elif (prepared['width'] != width or prepared['levels'] != tuple(levels)
          or prepared['u_step'] != u_step or prepared['logpotential'] != JB_LD(logpotential)):
        raise ValueError('Prepared bins do not match width, levels, bin width, or potential')
    dimension = len(prepared['keys'])
    positive = np.zeros((dimension, dimension), dtype=JB_LD)
    terminal_positive = np.zeros(dimension, dtype=JB_LD)
    physical_count = len(prepared['physical']['states'])
    started = time.monotonic()
    for source, u, e, outgoing, terminal in _joint_bin_physical_rows(prepared, m, b, p, r, alpha, coverage):
        aggregate = np.zeros(dimension, dtype=JB_LD)
        for j, kept in enumerate(prepared['levels']):
            sizes, starts, bins = prepared['target_groups'][kept]
            polynomial = outgoing[kept] if kept == levels[-1] else outgoing[kept:levels[j + 1]].sum(axis=0)
            weighted = polynomial[sizes] * prepared['phi'][sizes]
            aggregate[bins] = np.add.reduceat(weighted, starts)
        aggregate /= prepared['phi'][u]
        if u == 0:
            aggregate[0] = 0  # Remove exact Z before taking any source maxima.
        compressed = prepared['source_bins'][source]
        np.maximum(positive[compressed], aggregate, out=positive[compressed])
        terminal_positive[compressed] = max(terminal_positive[compressed], terminal / prepared['phi'][u])
        if progress and (source + 1) % 128 == 0:
            print(dict(joint_bins_stage='physical_rows', completed=source + 1, physical_dimension=physical_count,
                       compressed_dimension=dimension, seconds=round(time.monotonic() - started, 3)), flush=True)
    positive[:, 0] += prepared['zero_vector']
    if not np.all(np.isfinite(positive)) or not np.all(np.isfinite(terminal_positive)):
        raise ArithmeticError('Nonfinite compressed comparison; diagnostic invalid')
    return positive, terminal_positive, prepared['zero_vector'].copy()


def joint_bins_checks():
    _joint_bins_dependencies()
    F = Fraction
    m, b, p, r, alpha, coverage = F(3, 4), F(5, 4), F(1, 3), F(2, 5), F(2, 3), F(4, 5)
    comparisons = singleton_checks = zero_checks = endpoint_checks = region_checks = exact_comparisons = 0
    for width in (4, 8):
        for levels in ((0, 1), (0, 2, 4)):
            physical_prepared = joint_subset_prepare(width, levels)
            physical_states = physical_prepared['states']
            physical_index = {state: i for i, state in enumerate(physical_states)}
            exact = [[F(0)] * len(physical_states) for _ in physical_states]
            exact_terminal = [F(0)] * len(physical_states)
            for source, (u, e) in enumerate(physical_states):
                row, exact_terminal[source] = _joint_exact_row(width, (1 << u) - 1, ((1 << e) - 1) << u,
                                                              levels, m, b, p, r, alpha, coverage)
                for destination, weight in row.items():
                    if destination in physical_index:
                        exact[source][physical_index[destination]] += weight
            K, tp, v = joint_subset_matrix(width, float(m), float(b), float(p), float(r), float(alpha),
                                           levels, float(coverage), physical_prepared)
            P = K.copy()
            P[:, 0] -= v
            physical_Q, physical_scale, physical_end, physical_end_scale = joint_subset_region(K, tp, v)
            for step in (1, 2, 4):
                for logpotential in (0., 1.):
                    prepared = joint_bins_prepare(width, levels, step, logpotential)
                    B, cp, vb = joint_bins_matrix(width, float(m), float(b), float(p), float(r), float(alpha),
                                                  levels, step, logpotential, float(coverage), prepared)
                    H = np.zeros((len(physical_states), len(B)), dtype=JB_LD)
                    for s, (u, _) in enumerate(physical_states):
                        H[s, prepared['source_bins'][s]] = prepared['phi'][u]
                    Bplus = B.copy()
                    Bplus[:, 0] -= vb
                    left, right = P @ H, H @ Bplus
                    assert np.max(left - right) <= JB_LD('4e-13') * max(1, abs(right).max())
                    assert np.max(tp - H @ cp) <= JB_LD('4e-13') * max(1, abs(H @ cp).max())
                    assert np.array_equal(v, H @ vb)
                    for s, (u, _) in enumerate(physical_states):
                        phi_s = F(2) ** (u if logpotential else 0)
                        source_bin = prepared['source_bins'][s]
                        expected_endpoint = (exact_terminal[s] - int(u == 0)) / phi_s
                        actual_endpoint = JB_LD(expected_endpoint.numerator) / JB_LD(expected_endpoint.denominator)
                        assert actual_endpoint <= cp[source_bin] + JB_LD('4e-13') * max(1, abs(actual_endpoint))
                        endpoint_checks += 1
                        for target_bin in range(len(B)):
                            total = F(0)
                            for t, (new_u, _) in enumerate(physical_states):
                                if prepared['source_bins'][t] != target_bin:
                                    continue
                                probability = exact[s][t] - int(u == 0 and t == 0)
                                total += probability * F(2) ** (new_u if logpotential else 0) / phi_s
                            truth = JB_LD(total.numerator) / JB_LD(total.denominator)
                            assert truth <= Bplus[source_bin, target_bin] + JB_LD('4e-13') * max(1, abs(truth))
                            exact_comparisons += 1
                    if step == 1:
                        diagonal = np.array([prepared['phi'][u] for u, _ in physical_states])
                        transformed = K * diagonal[None, :] / diagonal[:, None]
                        assert np.max(abs(B - transformed)) <= JB_LD('4e-13') * max(1, abs(transformed).max())
                        assert np.max(abs(cp - tp / diagonal)) <= JB_LD('4e-13') * max(1, abs(tp / diagonal).max())
                        singleton_checks += 1
                    zero, endpoint_zero, zv = joint_bins_matrix(width, float(m), float(b), float(p), float(r),
                                                               float(alpha), levels, step, logpotential, 0., prepared)
                    expected_zero = np.zeros_like(zero)
                    expected_zero[:, 0] = zv
                    assert np.array_equal(zero, expected_zero) and np.all(endpoint_zero == 0)
                    physical_zero = np.zeros_like(K)
                    physical_zero[:, 0] = v
                    assert np.array_equal(physical_zero @ H, H @ expected_zero)
                    zero_checks += 1
                    Q, scale, end, end_scale = joint_subset_region(B, cp, vb)
                    left = (physical_Q * np.exp2(physical_scale)) @ H
                    right = H @ (Q * np.exp2(scale))
                    assert np.max(left - right) <= JB_LD('4e-12') * max(1, abs(right).max())
                    left_end = physical_end * np.exp2(physical_end_scale)
                    right_end = H @ (end * np.exp2(end_scale))
                    assert np.max(left_end - right_end) <= JB_LD('4e-12') * max(1, abs(right_end).max())
                    comparisons += len(physical_states) * len(B)
                    region_checks += 1
    print(dict(joint_bin_entry_comparisons=comparisons, exact_physical_row_inequalities=exact_comparisons,
               exact_endpoint_inequalities=endpoint_checks, singleton_bin_identities=singleton_checks,
               exact_zero_intertwining_checks=zero_checks, regional_comparison_checks=region_checks,
               full_code_certificate=False), flush=True)


def joint_bins_value(h=2048, theta=1/16, logz=-.4375, tilts=3, levels=(0, 1, 16),
                     u_step=4, logpotential=0., progress=False,logoutside=0.,logcoverage=0.,
                     mdsgap=0.):
    _joint_bins_dependencies()
    if tilts < 1 or not 0 < h < 8192:
        raise ValueError('Require positive tilt count and 0<h<8192')
    if logoutside>0 or logcoverage>0 or (logoutside and h>7884):
        raise ValueError('Upper-cap markers require nonpositive logs and the remaining rank class')
    if mdsgap < 0:
        raise ValueError('The aggregate MDS marker must have nonnegative log2 weight')
    m, b, alpha, d = joint_subset_field_weights(theta, logz)
    field = np.exp2(JB_LD(32 * theta))
    prepared = joint_bins_prepare(JB_WIDTH, levels, u_step, logpotential)
    rows = 8192
    actual = JB_LD(h) / rows
    low = actual / (field - (field - 1) * actual)
    choose = (lgamma(rows + 1) - lgamma(h + 1) - lgamma(rows - h + 1)) / float(np.log(JB_LD(2)))
    logs = endlogs = None
    for tilt_index, fraction in enumerate(np.linspace(0., 1., tilts)):
        started = time.monotonic()
        odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
        p = odds / (1 + odds)
        logpi = choose + h * np.log2(p) + (rows - h) * np.log2(1 - p)
        outside = (1 - p)*np.exp2(JB_LD(logoutside)) / field
        route_available = p + outside
        route_hole = p + outside * np.exp2(JB_LD(-4 * mdsgap))
        # One canonical marker: eta^(4L-F+f-4E_out) >= 1.
        # The ratio is only an evaluation shortcut. In the labeled polynomial,
        # available and hole coordinates have separate positive route factors;
        # the same quota-conditioned operator is used at every p.
        marked_m = m * np.exp2(JB_LD(5 * mdsgap))
        marked_b = b * np.exp2(JB_LD(4 * mdsgap)) * route_available / route_hole
        B, cp, vb = joint_bins_matrix(JB_WIDTH, marked_m, marked_b, p, route_hole,
                                     alpha, levels, u_step, logpotential,
                                     coverage=np.exp2(JB_LD(logcoverage - mdsgap)),
                                     prepared=prepared, progress=progress)
        if progress:
            print(dict(joint_bins_stage='matrix', tilt=tilt_index + 1, dimension=len(B),
                       seconds=round(time.monotonic() - started, 3)), flush=True)
        Q, scale, end, end_scale = joint_subset_region(B, cp, vb)
        if not np.all(Q > 0) or not np.all(end > 0):
            raise ArithmeticError('Unexpected regional zero/underflow in compressed bound')
        candidate = np.log2(Q) + scale - logpi
        candidate_end = np.log2(end) + end_scale - logpi
        # Common disjoint lift H permits these minima only after region conditioning.
        logs = candidate if logs is None else np.minimum(logs, candidate)
        endlogs = candidate_end if endlogs is None else np.minimum(endlogs, candidate_end)
        if progress:
            print(dict(joint_bins_stage='region', tilt=tilt_index + 1,
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    scale = logs.max()
    Q = np.exp2(logs - scale)
    norm = Q.sum(axis=1).max()
    Q /= norm
    scale += np.log2(norm)
    end_scale = endlogs.max()
    end = np.exp2(endlogs - end_scale)
    value = joint_subset_tail(Q, scale, end, end_scale, prepared['zero_vector'])
    coverage_cap = min(1736703, 1009255 + 128 * h)
    bound = float(theta - 32 * theta - JB_W * logz + choose - 32 * theta * 128 * h + value
                  -(1009255-128*h)*logoutside-coverage_cap*logcoverage)
    return bound, dict(dimension=len(prepared['keys']), physical_dimension=len(prepared['physical']['states']),
                       levels=tuple(levels), u_step=u_step, logpotential=logpotential, tilts=tilts,
                       alpha=float(alpha), closure_d=float(d),logoutside=logoutside,
                       logcoverage=logcoverage, coverage_cap=coverage_cap, mdsgap=mdsgap)


def joint_bins_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--tilts', type=int, default=3)
    parser.add_argument('--levels', default='0,1,16')
    parser.add_argument('--u-step', type=int, default=4)
    parser.add_argument('--logpotential', type=float, default=0.)
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--logz', type=float, default=-.4375)
    parser.add_argument('--outside', type=float, default=0.)
    parser.add_argument('--coverage', type=float, default=0.)
    parser.add_argument('--mdsgap', type=float, default=0.)
    parser.add_argument('--seconds', type=int, default=175)
    args = parser.parse_args()
    if not 1 <= args.seconds <= 350:
        parser.error('--seconds must be in 1..350')
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(args.seconds, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            joint_bins_checks()
            return
        started = time.monotonic()
        levels = tuple(int(value) for value in args.levels.split(','))
        bound, detail = joint_bins_value(args.h, args.theta, args.logz, args.tilts, levels,
                                          args.u_step, args.logpotential, progress=True,
                                          logoutside=args.outside,logcoverage=args.coverage,
                                          mdsgap=args.mdsgap)
        print(dict(joint_bins_rows=args.h, theta=args.theta, logz=args.logz,
                   log2_diagnostic=bound, seconds=round(time.monotonic() - started, 3),
                   certificate=False, **detail), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    joint_bins_main()
