"""Positive actual-S / eligible-label / accepted-MDS global moment.

The state is the actual support population u, not availability.  Markers
s,z,nu,eta count Ls, G=F-factual, factual, and covered-outside mass
Cout=R-128h.  The closure marker is t=(qx^2/(1+qx))^theta.

For a quartet with r old and v new state positions, let
  D_rv(t)=[Y^v](t+Y)^r(1+Y)^(4-r).
At eligible-label counter xi, put e=p*xi, o=(1-p)*eta,
U=1-p+p*xi.  The local polynomial has empty-old row
  (U+e*nu*s*Y)^4
and live row r>=1
  C(4,r) X^r z^4 sum_v D_rv(t)(sY)^v(e+o)^v
    sum_{j<=min(4-v,r-1)} C(4-v,j)e^(4-v-j)o^j.
The cap counts outside-H NEW ZEROS, not all outside positions.  It follows
from the accepted superregular A map and zero outer inputs outside H.

K_uv=[X^u Y^v]P^(width/4)/C(width,u).  At xi=1 every coefficient is
strictly positive for 0<p<1 and positive markers.  Assembly uses positive
25-shift convolutions with reusable buffers and rescaling; no FFT or
subtraction.  The final vector is all ones, with no terminal reset.

The profile prefactor includes C(8192,h)/pi_h(p)^256, all four marker
corrections, and the fractional conditional field factor.  THERE IS NO
SEPARATE ROUTING Q FACTOR.  All probes and Collatz bounds are FLOAT
diagnostics, not directed/exact certificates.  Small checks use independent
physical old/new/H masks and exact rational arithmetic.  No local imports,
subprocesses or data files.  Run numerical jobs serially.

Regional mode instead forms Q_p=K_p^16-Z and minimizes Q_p/pi_h(p)
entrywise in the SAME actual-popcount basis.  At formal xi the exact
empty-region term is (1-p+p*xi)^8192 Z.  Extracting xi^h and dividing
by pi_h(p) leaves Z exactly; the positive conditioned remainder is
independent of p and bounded by every Q_p/pi_h(p).  Thus R=Z+Qmin is
a common regional upper kernel, and its external prefactor omits pi^-256.
Qmin[0,0] counts nonzero regional excursions and must never be erased.

CLI --floor-majorant explicitly rounds normalized positive entries UP at
2^-4000, including each convolution normalization, and uses logarithmic
regional/mixed-term scaling.  These positive enlargements are mathematically
monotone, but the remaining floating arithmetic is not directed or certified.
The resulting bounds and --diagnostics describe the enlarged comparison
operator, not a physical-history distribution.  Default arithmetic is unchanged.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, lgamma, log
from operator import index
import os
import shutil
import threading
import time

import numpy as np


AM_LD = np.longdouble
AM_N = 1 << 21
AM_W = 943718
AM_Q = 1 << 32


def _am_floor_bits(bits):
    if bits is None:
        return None
    bits = index(bits)
    if not 1 <= bits <= 4000:
        raise ValueError('Positive-floor precision must be in 1..4000 bits')
    return bits


def _am_floor_positive(matrix, bits):
    """Add positive mass; preserve structural zeros.  FLOAT, not directed arithmetic."""
    if bits is not None:
        cutoff = np.exp2(-AM_LD(bits))
        selected = (matrix > 0) & (matrix < cutoff)
        matrix[selected] = cutoff


def _am_logs_majorant(log_entries, bits):
    """Represent positive log entries with a common scale and upward floor."""
    if bits is None:
        raise ValueError('Log-majorant assembly requires an explicit positive floor')
    if np.any(np.isnan(log_entries)) or np.any(np.isposinf(log_entries)):
        raise ArithmeticError('Invalid log entries in positive majorant')
    selected = np.isfinite(log_entries)
    if not np.any(selected):
        raise ArithmeticError('The positive majorant is empty')
    scale = log_entries[selected].max()
    matrix = np.zeros_like(log_entries)
    matrix[selected] = np.exp2(np.maximum(log_entries[selected] - scale, -AM_LD(bits)))
    return matrix, scale


def _am_choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def actual_mds_local_exact(s, t, z, nu, eta, p, xi=1):
    """Exact local coefficient table, including the formal H-label counter."""
    s, t, z, nu, eta, p, xi = map(Fraction, (s, t, z, nu, eta, p, xi))
    if min(s, t, z, nu, eta, xi) < 0 or not 0 <= p <= 1:
        raise ValueError('Exact local markers must be nonnegative and p in [0,1]')
    e, o, uncovered = p * xi, (1 - p) * eta, 1 - p + p * xi
    result = [[Fraction(0) for _ in range(5)] for _ in range(5)]
    for v in range(5):
        result[0][v] = comb(4, v) * (e * nu * s) ** v * uncovered ** (4 - v)
    for r in range(1, 5):
        for v in range(5):
            death = sum(Fraction(comb(r, kept) * _am_choose(4 - r, v - kept)) * t ** (r - kept)
                        for kept in range(r + 1))
            labels = sum(Fraction(comb(4 - v, outside)) * e ** (4 - v - outside) * o ** outside
                         for outside in range(min(4 - v, r - 1) + 1))
            result[r][v] = comb(4, r) * death * s ** v * z ** 4 * (e + o) ** v * labels
    return result


def _am_logsum(values):
    values = np.asarray(values, dtype=AM_LD)
    top = values.max()
    terms = np.exp2(values - top)
    if np.any(terms <= 0) or not np.all(np.isfinite(terms)):
        raise ArithmeticError('Local log sum lost a positive term')
    return top + np.log2(terms.sum())


def actual_mds_local_scaled(log_s, log_t, log_z, log_nu, log_eta, p, floor_bits=None):
    """Positive xi=1 local table normalized by one scalar power of two."""
    log_s, log_t, log_z, log_nu, log_eta, p = map(AM_LD, (log_s, log_t, log_z, log_nu, log_eta, p))
    if not all(np.isfinite(v) for v in (log_s, log_t, log_z, log_nu, log_eta, p)):
        raise ValueError('All logarithmic markers and p must be finite')
    if not 0 < p < 1:
        raise ValueError('The strictly positive floating kernel requires 0<p<1')
    lp, lo = np.log2(p), np.log2(1 - p) + log_eta
    all_labels = np.logaddexp2(lp, lo)
    logs = np.empty((5, 5), dtype=AM_LD)
    for v in range(5):
        logs[0, v] = np.log2(AM_LD(comb(4, v))) + v * (lp + log_nu + log_s)
    for r in range(1, 5):
        for v in range(5):
            deaths = [np.log2(AM_LD(comb(r, kept) * comb(4 - r, v - kept))) + (r - kept) * log_t
                      for kept in range(max(0, v - (4 - r)), min(r, v) + 1)]
            labels = [np.log2(AM_LD(comb(4 - v, outside))) + (4 - v - outside) * lp + outside * lo
                      for outside in range(min(4 - v, r - 1) + 1)]
            logs[r, v] = (np.log2(AM_LD(comb(4, r))) + _am_logsum(deaths)
                          + v * log_s + 4 * log_z + v * all_labels + _am_logsum(labels))
    scale = logs.max()
    floor_bits = _am_floor_bits(floor_bits)
    table = np.exp2(logs - scale if floor_bits is None else np.maximum(logs - scale, -AM_LD(floor_bits)))
    if not np.all(np.isfinite(table)) or np.any(table <= 0):
        raise ArithmeticError('Local scaling underflowed or overflowed a positive coefficient')
    return table, scale


def _am_positive(matrix, stage, allow_zero=False):
    if not np.all(np.isfinite(matrix)) or (np.any(matrix < 0) if allow_zero else np.any(matrix <= 0)):
        invalid = ~np.isfinite(matrix) | (matrix < 0 if allow_zero else matrix <= 0)
        coordinates = np.argwhere(invalid)[:4].tolist()
        raise ArithmeticError('Invalid or underflowed positive data at ' + stage
                              + '; first coordinates=' + str(coordinates))


def _am_operator_diagnostics(matrix, scale, right, right_action, iterations=32,
                              top=8, selector_data=None):
    """Bounded approximate Perron diagnostics of the COMPARISON operator.

    Entrywise p minima do not define one physical Bernoulli law.  The
    population and flow masses below are diagnostics of the upper matrix,
    not a distribution of actual histories or a certificate of dominance.
    """
    iterations, top = index(iterations), index(top)
    if not 1 <= iterations <= 96 or not 1 <= top <= 32:
        raise ValueError('Diagnostics require 1..96 iterations and top in 1..32')
    scale = AM_LD(scale)
    right_logs = np.log2(right)
    right_residual = scale + np.log2(right_action / right)
    rows = np.argsort(right_residual)[-top:][::-1]
    result = dict(interpretation='approximate_Perron_diagnostic_of_comparison_operator',
                  physical_history_distribution=False, additional_left_iterations=iterations,
                  right_residual_log2_min=float(right_residual.min()),
                  right_residual_log2_max=float(right_residual.max()),
                  largest_right_residual_rows=[dict(population=int(u),
                      log2_ratio=float(right_residual[u]), log2_right_relative=float(right_logs[u])) for u in rows])
    if selector_data is not None:
        result['regional_minimum_selectors'] = selector_data['summary']
    left = np.ones(len(matrix), dtype=AM_LD)
    for _ in range(iterations):
        left = matrix.T @ left
        _am_positive(left, 'diagnostic left iteration')
        left /= left.max()
    left_action = matrix.T @ left
    _am_positive(left_action, 'diagnostic left residual')
    left_residual = scale + np.log2(left_action / left)
    left_logs = np.log2(left)
    population_log = left_logs + right_logs
    population = np.exp2(population_log - population_log.max())
    population /= population.sum()
    populations = np.arange(len(matrix), dtype=AM_LD)
    mean = (populations * population).sum()
    top_populations = np.argsort(population)[-top:][::-1]
    result.update(left_residual_log2_min=float(left_residual.min()),
                  left_residual_log2_max=float(left_residual.max()),
                  approximate_mean_population=float(mean),
                  approximate_population_stddev=float(np.sqrt(((populations - mean) ** 2 * population).sum())),
                  approximate_zero_population_mass=float(population[0]),
                  approximate_full_population_mass=float(population[-1]),
                  top_populations=[dict(population=int(u), mass=float(population[u]),
                      log2_left_relative=float(left_logs[u]), log2_right_relative=float(right_logs[u])) for u in top_populations])
    flow_logs = left_logs[:, None] + np.log2(matrix) + right_logs[None, :]
    flow_max = flow_logs.max()
    flow = np.exp2(flow_logs - flow_max)
    flow /= flow.sum()
    top_edges = np.argsort(flow, axis=None)[-top:][::-1]
    edges = []
    for position in top_edges:
        u, v = np.unravel_index(position, flow.shape)
        entry = dict(source_population=int(u), target_population=int(v),
                     comparison_flow_mass=float(flow[u, v]))
        if selector_data is not None:
            mask = int(selector_data['masks'][u, v])
            entry['positive_part_minimizing_p_values'] = [float(p) for i, p in enumerate(selector_data['p_values']) if mask >> i & 1]
        edges.append(entry)
    result['top_comparison_transitions'] = edges
    result['diagnostic_underflowed_flow_entries'] = int(np.count_nonzero(flow == 0))
    if selector_data is not None:
        # Only entry (0,0) also contains the exact empty-region operator Z.
        # Split its flow using the saved positive entry, without subtraction.
        positive_flow = flow.copy()
        positive_flow[0, 0] *= selector_data['positive00_scaled'] / matrix[0, 0]
        zero_flow = flow[0, 0] * selector_data['empty00_scaled'] / matrix[0, 0]
        masks = selector_data['masks']
        result['empty_region_flow_mass'] = float(zero_flow)
        result['empty_split_approximate_under_floor'] = selector_data.get('empty_split_approximate_under_floor', False)
        result['positive_flow_mass_by_unique_selector'] = [dict(p=float(p),
            mass=float(positive_flow[masks == (1 << i)].sum())) for i, p in enumerate(selector_data['p_values'])]
        result['positive_flow_mass_with_tied_selectors'] = float(sum(
            positive_flow[masks == pattern['mask']].sum()
            for pattern in selector_data['summary']['patterns'] if len(pattern['p_values']) > 1))
        result['selector_counts_at_top_populations'] = [dict(population=int(u),
            entries_attaining_min_by_p=[int(np.count_nonzero(masks[u] & (1 << i)))
                                       for i in range(len(selector_data['p_values']))]) for u in top_populations]
        result['zero_boundary_self_loop_log2'] = float(scale + np.log2(matrix[0, 0]))
    return result


def actual_mds_matrix(width=512, log_s=1, log_t=-1, log_z=0, log_nu=0,
                      log_eta=0, p=.5, progress=False, floor_bits=None):
    """Return normalized K and its log2 scale; O(width^3) positive assembly.

    Three fixed buffers avoid per-shift matrix allocation in the hot loop.
    Every stage is rescaled, and all mathematically positive entries must
    remain finite and strictly positive.  Unsupported marker ranges fail.
    """
    width = index(width)
    if width < 4 or width % 4:
        raise ValueError('Width must be a positive integer multiple of four')
    floor_bits = _am_floor_bits(floor_bits)
    local, local_scale = actual_mds_local_scaled(log_s, log_t, log_z, log_nu, log_eta, p, floor_bits)
    shape = width + 1, width + 1
    current = np.zeros(shape, dtype=AM_LD)
    after = np.empty(shape, dtype=AM_LD)
    scratch = np.empty(shape, dtype=AM_LD)
    current[0, 0] = 1
    scale = AM_LD(0)
    started = time.monotonic()
    for group in range(width // 4):
        size = 4 * group + 1
        destination_size = size + 4
        target = after[:destination_size, :destination_size]
        target.fill(0)
        source = current[:size, :size]
        work = scratch[:size, :size]
        if source.min() * local.min() <= 0:
            raise ArithmeticError('A positive convolution product underflowed')
        for r in range(5):
            for v in range(5):
                np.multiply(source, local[r, v], out=work)
                np.add(after[r:r + size, v:v + size], work,
                       out=after[r:r + size, v:v + size])
        _am_positive(target, 'quartet polynomial convolution')
        norm = target.max()
        np.divide(target, norm, out=target)
        _am_floor_positive(target, floor_bits)
        _am_positive(target, 'normalized quartet convolution')
        scale += local_scale + np.log2(norm)
        current, after = after, current
        if progress and ((group + 1) % 32 == 0 or group + 1 == width // 4):
            print(dict(stage='actual_mds_matrix', quartets=group + 1,
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    inverse_choose = np.array([AM_LD(1) / AM_LD(comb(width, u)) for u in range(width + 1)])
    np.multiply(current, inverse_choose[:, None], out=current)
    _am_positive(current, 'old-mask normalization')
    norm = current.sum(axis=1).max()
    np.divide(current, norm, out=current)
    _am_floor_positive(current, floor_bits)
    scale += np.log2(norm)
    _am_positive(current, 'final kernel normalization')
    return current, scale


def actual_mds_collatz(matrix, scale, updates=4096, iterations=96,
                       diagnostic_top=0, diagnostic_iterations=32, selector_data=None):
    """FLOAT Collatz upper diagnostic, including v0/min(v) for the free end."""
    updates, iterations = index(updates), index(iterations)
    if updates < 0 or iterations < 1:
        raise ValueError('Require nonnegative updates and positive iterations')
    _am_positive(matrix, 'Collatz input')
    vector = np.ones(len(matrix), dtype=AM_LD)
    for _ in range(iterations):
        vector = matrix @ vector
        _am_positive(vector, 'Collatz vector')
        vector /= vector.max()
    action = matrix @ vector
    _am_positive(action, 'Collatz final action')
    growth = AM_LD(scale) + np.log2((action / vector).max())
    boundary = np.log2(vector[0]) - np.log2(vector.min())
    result = dict(log2_moment=float(updates * growth + boundary),
                log2_growth=float(growth), log2_free_endpoint=float(boundary),
                method='floating_collatz', includes_empty_history=True)
    if diagnostic_top:
        try:
            result['operator_diagnostics'] = _am_operator_diagnostics(matrix, scale, vector, action,
                diagnostic_iterations, diagnostic_top, selector_data)
        except ArithmeticError as error:
            result['operator_diagnostics'] = dict(status='diagnostic_range_failure', message=str(error),
                                                   bound_calculation_unchanged=True)
    return result


def actual_mds_positive_part_power(positive, scale, updates, floor_bits=None):
    """Power the explicitly supplied Q=R-Z; never erase its (0,0) entry.

    Q_2m=Q_m^2+Q_m Z+Z Q_m.  The input can have Q00=0 for a single
    update, or Q00>0 for a region containing internal nonzero excursions.
    """
    updates = index(updates)
    if updates < 1 or updates & (updates - 1):
        raise ValueError('Positive-power mode requires a positive power-of-two update count')
    positive = positive.copy()
    floor_bits = _am_floor_bits(floor_bits)
    _am_floor_positive(positive, floor_bits)
    scale = AM_LD(scale)
    _am_positive(positive, 'explicit positive part', allow_zero=True)
    for _ in range(updates.bit_length() - 1):
        square = positive @ positive
        _am_positive(square, 'matrix square')
        edges = np.zeros_like(positive)
        edges[0, :] += positive[0, :]
        edges[:, 0] += positive[:, 0]
        if floor_bits is not None:
            log_square = np.log2(square) + 2 * scale
            log_edges = np.full_like(edges, AM_LD('-inf'))
            np.log2(edges, out=log_edges, where=edges > 0)
            log_edges += scale
            positive, scale = _am_logs_majorant(np.logaddexp2(log_square, log_edges), floor_bits)
            continue
        common = max(2 * scale, scale)
        a, b = np.exp2(2 * scale - common), np.exp2(scale - common)
        if min(a, b) <= 0 or not np.isfinite(a + b):
            raise ArithmeticError('Nonzero power scaling lost a positive mixed term; use Collatz')
        scaled_square, scaled_edges = square * a, edges * b
        _am_positive(scaled_square, 'scaled positive square')
        if not np.all(np.isfinite(scaled_edges)) or np.any((edges > 0) & (scaled_edges <= 0)):
            raise ArithmeticError('Nonzero power scaling lost a positive edge contribution')
        positive = scaled_square + scaled_edges
        scale = common
        norm = positive.sum(axis=1).max()
        positive /= norm
        scale += np.log2(norm)
        _am_positive(positive, 'power normalization')
    return positive, scale


def actual_mds_power_moment(matrix, scale, updates=4096, nonzero=False, floor_bits=None):
    """One-update K power diagnostic; K00=1 permits removing that entry here.

    For an already composed region use actual_mds_positive_part_power
    with its explicit positive part instead; its Q00 must be preserved.
    """
    updates = index(updates)
    if updates < 1 or updates & (updates - 1):
        raise ValueError('Power mode requires a positive power-of-two update count')
    positive = matrix.copy()
    floor_bits = _am_floor_bits(floor_bits)
    _am_floor_positive(positive, floor_bits)
    scale = AM_LD(scale)
    if nonzero:
        positive[0, 0] = 0
        positive, scale = actual_mds_positive_part_power(positive, scale, updates, floor_bits)
    else:
        for _ in range(updates.bit_length() - 1):
            positive = positive @ positive
            _am_positive(positive, 'full matrix square')
            norm = positive.sum(axis=1).max()
            positive /= norm
            _am_floor_positive(positive, floor_bits)
            scale = 2 * scale + np.log2(norm)
            _am_positive(positive, 'full power normalization')
    total = positive[0].sum()
    if total <= 0 or not np.isfinite(total):
        raise ArithmeticError('Invalid free-endpoint power sum')
    return dict(log2_moment=float(scale + np.log2(total)),
                method='positive_nonzero_powers' if nonzero else 'positive_full_powers',
                includes_empty_history=not nonzero)


def _am_logchoose(n, k):
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / log(2.0)


def actual_mds_full_state_floor(profile=(4096, 786432, 786432, 4, 1048576),
                                exponent=16, theta=1/16, p=.25, log_s=1,
                                log_z=0, log_nu=0, log_eta=None):
    """Cheap analytic path floor for the RELAXED moment, evaluated in FLOAT.

    S_1=...=S_T=[512] survives every MDS cutoff.  Its contribution is
      (p*nu*s)^512 [s*z*(p+(1-p)*eta)]^(512*(T-1)).
    This need not have the requested profile; that is precisely why it
    diagnoses loss from a marker bound.  It is not a codeword witness.
    """
    h, R, Ls, factual, F = map(index, profile)
    exponent, theta, p, log_s, log_z, log_nu = map(AM_LD, (exponent, theta, p, log_s, log_z, log_nu))
    log_eta = -np.log2(AM_LD(9)) if log_eta is None else AM_LD(log_eta)
    if not (1 <= h <= 8192 and 0 < p < 1 and exponent >= 0 and 0 <= theta <= 1
            and all(np.isfinite(v) for v in (exponent, theta, p, log_s, log_z, log_nu, log_eta))):
        raise ValueError('Invalid full-state floor parameters')
    outside, G = R - 128 * h, F - factual
    eligible_union = AM_LD(_am_logchoose(8192, h))
    quota = eligible_union + h * np.log2(p) + (8192 - h) * np.log2(1 - p)
    log_lambda = np.logaddexp2(AM_LD(0), 32 - exponent)
    prefactor = (eligible_union - 256 * quota - log_s * Ls - log_z * G
                 - log_nu * factual - log_eta * outside
                 + theta * (-31 - 32 * R + AM_W * exponent + Ls * log_lambda))
    labels = np.logaddexp2(np.log2(p), np.log2(1 - p) + log_eta)
    floor = 512 * (np.log2(p) + log_nu + log_s) + (AM_N - 512) * (log_s + log_z + labels)
    return dict(full_state_path_log2_moment_floor=float(floor),
                full_state_path_log2_profile_floor=float(prefactor + floor),
                floor_rules_out_target=bool(prefactor + floor > -172),
                target_exponent=-172, floating_diagnostic_only=True,
                path_is_not_a_codeword_or_profile_witness=True)


def actual_mds_profile(profile=(4096, 786432, 786432, 4, 1048576),
                       exponent=16, theta=1/16, p=.25, log_s=1, log_z=0,
                       log_nu=0, log_eta=None, iterations=96,
                       method='collatz', progress=False, diagnostic_top=0, diagnostic_iterations=32,
                       floor_bits=None):
    """Single fixed-profile FLOAT diagnostic; eligible routing is already inside."""
    h, R, Ls, factual, F = map(index, profile)
    if not (1 <= h <= 8192 and 128 * h <= R <= 1 << 20
            and 1 <= factual <= Ls <= F <= AM_N):
        raise ValueError('Invalid actual-state profile ranges')
    C_H, outside, G = F - R + 128 * h, R - 128 * h, F - factual
    if not (129 * h <= C_H <= 256 * h and factual <= C_H
            and outside <= 128 * (8192 - h) and G % 4 == 0):
        raise ValueError('Profile violates covered-row or forced-birth support')
    exponent, theta, p = AM_LD(exponent), AM_LD(theta), AM_LD(p)
    if not (np.isfinite(exponent) and exponent >= 0 and 0 <= theta <= 1 and 0 < p < 1):
        raise ValueError('Require exponent>=0, theta in [0,1], and 0<p<1')
    log_eta = -np.log2(AM_LD(9)) if log_eta is None else AM_LD(log_eta)
    log_s, log_z, log_nu = map(AM_LD, (log_s, log_z, log_nu))
    log_x = -exponent
    log_lambda = np.logaddexp2(AM_LD(0), 32 + log_x)
    log_t = theta * (32 + 2 * log_x - log_lambda)
    if diagnostic_top and method != 'collatz':
        raise ValueError('The bounded population diagnostic currently requires Collatz mode')
    matrix, scale = actual_mds_matrix(512, log_s, log_t, log_z, log_nu, log_eta, p, progress, floor_bits)
    if method == 'collatz':
        moment = actual_mds_collatz(matrix, scale, 4096, iterations, diagnostic_top, diagnostic_iterations)
    elif method in ('powers', 'nonzero'):
        moment = actual_mds_power_moment(matrix, scale, 4096, method == 'nonzero', floor_bits)
    else:
        raise ValueError('Unknown moment method')
    eligible_union = AM_LD(_am_logchoose(8192, h))
    quota = eligible_union + h * np.log2(p) + (8192 - h) * np.log2(1 - p)
    field = theta * (-31 - 32 * R - AM_W * log_x + Ls * log_lambda)
    marker_correction = -log_s * Ls - log_z * G - log_nu * factual - log_eta * outside
    prefactor = eligible_union - 256 * quota + marker_correction + field
    value = AM_LD(moment['log2_moment']) + prefactor
    floor = actual_mds_full_state_floor(profile, exponent, theta, p, log_s, log_z, log_nu, log_eta)
    return dict(statistics='actual_state_support', profile=(h, R, Ls, factual, F),
                Cout=outside, G=G, exponent=float(exponent), theta=float(theta), p=float(p),
                log_s=float(log_s), log_t=float(log_t), log_z=float(log_z),
                log_nu=float(log_nu), log_eta=float(log_eta),
                eligible_union=float(eligible_union), log2_regional_quota=float(quota),
                marker_correction=float(marker_correction), field_prefactor=float(field),
                combined_prefactor=float(prefactor), value=float(value),
                separate_routing_factor_used=False, floating_diagnostic_only=True,
                positive_floor_bits=floor_bits, conservative_floor_diagnostic=floor_bits is not None,
                full_code_certificate=False,
                full_state_path_log2_moment_floor=floor['full_state_path_log2_moment_floor'],
                full_state_path_log2_profile_floor=floor['full_state_path_log2_profile_floor'],
                floor_rules_out_target=floor['floor_rules_out_target'], **moment)


def actual_mds_regional_kernel(width=512, region_updates=16, eligible=4096,
                               p_values=(1/5, 1/3, 1/2), log_s=1, log_t=-1,
                               log_z=0, log_nu=0, log_eta=-2, progress=False,
                               selector_data=None, floor_bits=None):
    """Compose, quota-adjust, THEN minimize in the common actual-popcount basis.

    Each Q_p=K_p^region_updates-Z bounds the nonempty conditional region
    after division by pi_eligible(p).  Empty regions have their exact Z,
    not Z/pi.  Regional Q00>0 counts internal excursions and is retained.
    The returned positive part is entrywise min_p Q_p/pi, with one scale.
    """
    width, region_updates, eligible = map(index, (width, region_updates, eligible))
    floor_bits = _am_floor_bits(floor_bits)
    slots = width * region_updates
    if region_updates < 1 or region_updates & (region_updates - 1) or not 0 <= eligible <= slots:
        raise ValueError('Require a power-of-two region length and a supported quota')
    values = tuple(AM_LD(p) for p in p_values)
    if not values or any(not 0 < p < 1 for p in values):
        raise ValueError('Regional tilts must lie strictly between zero and one')
    if selector_data is not None and len(values) > 63:
        raise ValueError('Selector diagnostics support at most 63 tilts')
    choose = AM_LD(_am_logchoose(slots, eligible))
    common_positive = None
    common_scale = None
    details = []
    selectors = None
    for tilt_index, p in enumerate(values):
        started = time.monotonic()
        if progress:
            print(dict(stage='actual_mds_region_tilt', p=float(p)), flush=True)
        kernel, scale = actual_mds_matrix(width, log_s, log_t, log_z, log_nu, log_eta, p, progress, floor_bits)
        one_positive = kernel.copy()
        one_positive[0, 0] = 0  # Valid ONLY at the one-update K00=1 initialization.
        positive, qscale = actual_mds_positive_part_power(one_positive, scale, region_updates, floor_bits)
        quota = choose + eligible * np.log2(p) + (slots - eligible) * np.log2(1 - p)
        qscale -= quota
        if common_positive is None:
            common_positive, common_scale = positive, qscale
            if selector_data is not None:
                selectors = np.ones(positive.shape, dtype=np.uint64)
        else:
            if floor_bits is None:
                shared = max(common_scale, qscale)
                left = common_positive * np.exp2(common_scale - shared)
                right = positive * np.exp2(qscale - shared)
                _am_positive(left, 'regional old bound rescaling', allow_zero=region_updates == 1)
                _am_positive(right, 'regional new bound rescaling', allow_zero=region_updates == 1)
            else:
                left = np.full_like(common_positive, AM_LD('-inf'))
                right = np.full_like(positive, AM_LD('-inf'))
                np.log2(common_positive, out=left, where=common_positive > 0)
                np.log2(positive, out=right, where=positive > 0)
                left += common_scale
                right += qscale
            if selector_data is not None:
                smaller, tied = right < left, right == left
                selectors[smaller] = 1 << tilt_index
                selectors[tied] |= np.uint64(1 << tilt_index)
            if floor_bits is None:
                common_positive = np.minimum(left, right)
                norm = common_positive.sum(axis=1).max()
                common_positive /= norm
                common_scale = shared + np.log2(norm)
            else:
                common_positive, common_scale = _am_logs_majorant(np.minimum(left, right), floor_bits)
            _am_positive(common_positive, 'regional minimum normalization', allow_zero=region_updates == 1)
        details.append(dict(p=float(p), log2_quota=float(quota),
                            positive_scale=float(qscale), seconds=round(time.monotonic() - started, 3)))
    if region_updates > 1 and common_positive[0, 0] <= 0:
        raise ArithmeticError('Regional nonzero excursions were lost at (0,0)')
    if selector_data is not None:
        patterns, counts = np.unique(selectors, return_counts=True)
        summary_patterns = [dict(mask=int(mask), entry_count=int(count),
            p_values=[float(p) for i, p in enumerate(values) if int(mask) >> i & 1])
            for mask, count in zip(patterns, counts)]
        selector_data.update(masks=selectors, p_values=values,
            summary=dict(comparisons='FLOAT entries in the common scaled basis',
                positive_floor_bits=floor_bits,
                integer_counts_exact_for_recorded_comparisons=True,
                analytic_argmin_certificate=False, total_entries=int(selectors.size),
                patterns=summary_patterns,
                tied_entries=sum(item['entry_count'] for item in summary_patterns if len(item['p_values']) > 1),
                entries_attaining_min_by_p=[dict(p=float(p),
                    entries=int(np.count_nonzero(selectors & (1 << i)))) for i, p in enumerate(values)]))
    return common_positive, common_scale, details


def actual_mds_regional_profile(profile=(4096, 786432, 786432, 4, 1048576),
                                exponent=16, theta=1/16, p_values=(1/5, 1/3, 1/2),
                                log_s=1, log_z=0, log_nu=0, log_eta=-2,
                                iterations=96, method='collatz', log_region=0,
                                minimum_regions=0, progress=False,
                                diagnostic_top=0, diagnostic_iterations=32, floor_bits=None):
    """Three-tilt regional control; external prefactor has NO quota divisor.

    For a supplied lower bound rho on occupied regions, kappa>=1 weights
    Qmin by kappa and contributes kappa^-rho externally.  By default rho=0.
    The profile guarantees max(129,ceil(C_H/h),ceil(Cout/(8192-h))) occupied
    regions; omit the last ratio when h=8192.  Larger supplied values need
    a separate proof and are not accepted by this standalone control.
    """
    h, R, Ls, factual, F = map(index, profile)
    minimum_regions = index(minimum_regions)
    if not (1 <= h <= 8192 and 128 * h <= R <= 1 << 20
            and 1 <= factual <= Ls <= F <= AM_N and (F - factual) % 4 == 0):
        raise ValueError('Invalid actual-state profile ranges')
    C_H, outside, G = F - R + 128 * h, R - 128 * h, F - factual
    if not (129 * h <= C_H <= 256 * h and factual <= C_H and outside <= 128 * (8192 - h)):
        raise ValueError('Invalid covered-row or forced-birth support')
    guaranteed_regions = max(129, (C_H + h - 1) // h)
    if h < 8192:
        guaranteed_regions = max(guaranteed_regions, (outside + 8192 - h - 1) // (8192 - h))
    exponent, theta, log_s, log_z, log_nu, log_eta, log_region = map(
        AM_LD, (exponent, theta, log_s, log_z, log_nu, log_eta, log_region))
    if not (exponent >= 0 and 0 <= theta <= 1 and log_region >= 0
            and 0 <= minimum_regions <= guaranteed_regions
            and all(np.isfinite(v) for v in (exponent, theta, log_s, log_z, log_nu, log_eta, log_region))):
        raise ValueError('Invalid markers or unsupported occupied-region lower bound')
    log_lambda = np.logaddexp2(AM_LD(0), 32 - exponent)
    log_t = theta * (32 - 2 * exponent - log_lambda)
    if diagnostic_top and method != 'collatz':
        raise ValueError('The bounded population diagnostic currently requires Collatz mode')
    selector_data = {} if diagnostic_top else None
    positive, scale, tilts = actual_mds_regional_kernel(512, 16, h, p_values,
        log_s, log_t, log_z, log_nu, log_eta, progress, selector_data, floor_bits)
    scale += log_region
    # This Qmin already has a positive (0,0) entry.  NEVER overwrite it.
    excursion_log = scale + np.log2(positive[0, 0])
    if method == 'nonzero':
        powered, power_scale = actual_mds_positive_part_power(positive, scale, 256, floor_bits)
        moment = dict(log2_moment=float(power_scale + np.log2(powered[0].sum())),
                      method='regional_positive_nonzero_powers', includes_empty_history=False)
    else:
        if floor_bits is None:
            shared = max(AM_LD(0), scale)
            regional = positive * np.exp2(scale - shared)
            _am_positive(regional, 'weighted regional positive part')
            empty_weight = np.exp2(-shared)
            if empty_weight <= 0:
                raise ArithmeticError('Regional assembly underflowed the exact empty operator')
            if selector_data is not None:
                selector_data.update(positive00_scaled=regional[0, 0], empty00_scaled=empty_weight)
            regional[0, 0] += empty_weight
        else:
            logs = np.log2(positive) + scale
            live00_log = logs[0, 0]
            logs[0, 0] = np.logaddexp2(live00_log, AM_LD(0))
            combined00_log = logs[0, 0]
            regional, shared = _am_logs_majorant(logs, floor_bits)
            if selector_data is not None:
                # If this entry was floored, its added comparison mass is
                # split using the unfloored local ratio solely for diagnostics.
                selector_data.update(positive00_scaled=regional[0, 0] * np.exp2(live00_log - combined00_log),
                    empty00_scaled=regional[0, 0] * np.exp2(-combined00_log),
                    empty_split_approximate_under_floor=True)
        if method == 'collatz':
            moment = actual_mds_collatz(regional, shared, 256, iterations,
                diagnostic_top, diagnostic_iterations, selector_data)
            moment['method'] = 'regional_floating_collatz'
        elif method == 'powers':
            moment = actual_mds_power_moment(regional, shared, 256, nonzero=False, floor_bits=floor_bits)
            moment['method'] = 'regional_positive_full_powers'
        else:
            raise ValueError('Unknown regional moment method')
    eligible_union = AM_LD(_am_logchoose(8192, h))
    field = theta * (-31 - 32 * R + AM_W * exponent + Ls * log_lambda)
    marker = (-log_s * Ls - log_z * G - log_nu * factual
              - log_eta * outside - log_region * minimum_regions)
    prefactor = eligible_union + marker + field  # Quotas were charged inside each Q_p.
    if diagnostic_top and 'zero_boundary_self_loop_log2' in moment.get('operator_diagnostics', {}):
        self_loop = moment['operator_diagnostics']['zero_boundary_self_loop_log2']
        moment['operator_diagnostics'].update(zero_boundary_path_log2_comparison_floor=256 * self_loop,
            zero_boundary_path_log2_profile_comparison_floor=float(prefactor + 256 * self_loop))
    return dict(statistics='actual_state_support', profile=tuple(profile),
                regional_tilts=tilts, regional_minimum_after_composition=True,
                regional_positive_00_log2=float(excursion_log),
                p_values=[float(p) for p in p_values], log_eta=float(log_eta),
                log_s=float(log_s), log_z=float(log_z), log_nu=float(log_nu),
                log_t=float(log_t), exponent=float(exponent), theta=float(theta),
                minimum_occupied_regions=minimum_regions, log_region=float(log_region),
                eligible_union=float(eligible_union), field_prefactor=float(field),
                marker_correction=float(marker), combined_prefactor=float(prefactor),
                value=float(AM_LD(moment['log2_moment']) + prefactor),
                external_quota_divisor_used=False, separate_routing_factor_used=False,
                positive_floor_bits=floor_bits, conservative_floor_diagnostic=floor_bits is not None,
                floating_diagnostic_only=True, full_code_certificate=False, **moment)


def _am_physical_local(s, t, z, nu, eta, p, xi=1):
    """Independent labeled old/new/H enumeration; return mask and count tables."""
    s, t, z, nu, eta, p, xi = map(Fraction, (s, t, z, nu, eta, p, xi))
    masks = [[Fraction(0) for _ in range(16)] for _ in range(16)]
    count = [[[Fraction(0) for _ in range(5)] for _ in range(5)] for _ in range(5)]
    for old in range(16):
        r = old.bit_count()
        for new in range(16):
            v = new.bit_count()
            covered = 15 if old else new
            base = s ** v * t ** (old & ~new).bit_count() * z ** (4 if old else 0)
            if not old:
                base *= nu ** v
            for eligible in range(16):
                if (not old and new & ~eligible) or (old and (15 & ~eligible & ~new).bit_count() > r - 1):
                    continue
                h = eligible.bit_count()
                weight = base * eta ** (covered & ~eligible).bit_count() * p ** h * (1 - p) ** (4 - h)
                count[r][v][h] += weight
                masks[old][new] += weight * xi ** h
    return masks, count


def _am_exact_power(local, groups):
    current = [[Fraction(1)]]
    for _ in range(groups):
        size = len(current)
        after = [[Fraction(0) for _ in range(size + 4)] for _ in range(size + 4)]
        for u in range(size):
            for v in range(size):
                for r in range(5):
                    for w in range(5):
                        after[u + r][v + w] += current[u][v] * local[r][w]
        current = after
    width = 4 * groups
    return [[current[u][v] / comb(width, u) for v in range(width + 1)] for u in range(width + 1)]


def _am_partitions8():
    return [(first, 255 ^ first) for rest in combinations(range(1, 8), 3)
            for first in [1 + sum(1 << j for j in rest)]]


def _am_width8_reference(mask_weights):
    result = [[Fraction(0) for _ in range(9)] for _ in range(9)]
    partitions = _am_partitions8()
    for partition in partitions:
        maps = []
        for quartet in partition:
            positions = [j for j in range(8) if quartet >> j & 1]
            maps.append([sum(((mask >> j) & 1) << k for k, j in enumerate(positions)) for mask in range(256)])
        for u in range(9):
            old = (1 << u) - 1
            a, b = maps[0][old], maps[1][old]
            for new in range(256):
                result[u][new.bit_count()] += mask_weights[a][maps[0][new]] * mask_weights[b][maps[1][new]]
    return [[value / len(partitions) for value in row] for row in result]


def _am_matmul(left, right):
    size = len(left)
    return [[sum(left[i][k] * right[k][j] for k in range(size)) for j in range(size)] for i in range(size)]


def actual_mds_checks():
    local_checks = endpoint_checks = matrix_checks = quota_checks = nonzero_checks = regional_checks = 0
    F = Fraction
    parameters = [(F(3, 4), F(2, 3), F(5, 4), F(1, 2), F(7, 8), F(1, 3), F(3, 2)),
                  (F(3, 4), F(2, 3), F(5, 4), F(1, 2), F(7, 8), F(1, 3), F(1)),
                  (F(1), F(0), F(1), F(1), F(1), F(1, 2), F(1)),
                  (F(1), F(1), F(1), F(1), F(1), F(0), F(1)),
                  (F(1), F(1), F(1), F(1), F(1), F(1), F(1)),
                  (F(1), F(1), F(1), F(1), F(1), F(1, 2), F(0))]
    for params in parameters:
        local = actual_mds_local_exact(*params)
        masks, counter = _am_physical_local(*params)
        xi = params[-1]
        for r in range(5):
            for v in range(5):
                assert local[r][v] == sum(counter[r][v][h] * xi ** h for h in range(5))
                local_checks += 1
        if params[-2] == 0:
            for r in range(1, 5):
                for v in range(5):
                    assert bool(local[r][v]) == (v >= 5 - r)
                    endpoint_checks += 1
    params = parameters[1]
    local = actual_mds_local_exact(*params)
    masks, _ = _am_physical_local(*params)
    expected8 = _am_width8_reference(masks)
    actual8 = _am_exact_power(local, 2)
    assert actual8 == expected8
    matrix_checks += 81
    actual4 = _am_exact_power(local, 1)
    s, t, z, nu, eta, p, _ = params
    assert actual4[0][4] == (p * nu * s) ** 4
    assert actual4[4][4] == (s * z * (p + (1 - p) * eta)) ** 4
    endpoint_checks += 2
    for u in range(5):
        old = (1 << u) - 1
        for v in range(5):
            assert actual4[u][v] == sum(masks[old][new] for new in range(16) if new.bit_count() == v)
            matrix_checks += 1
    # Free-terminal sum: two actual-mask updates, no extra expansion or reset.
    direct = sum(masks[0][first] * masks[first][last] for first in range(16) for last in range(16))
    square = _am_matmul(actual4, actual4)
    assert sum(square[0]) == direct
    endpoint_checks += 1
    # A two-update toy region has eight H labels.  Extract its total H count
    # in the SAME u basis, then condition by its binomial probability.
    conditional = None
    regional_upper = []
    for p in (F(1, 3), F(2, 3)):
        _, coefficient = _am_physical_local(F(3, 4), F(2, 3), F(5, 4), F(1, 2), F(7, 8), p)
        one = [[[coefficient[u][v][h] / comb(4, u) for h in range(5)] for v in range(5)] for u in range(5)]
        region = [[[Fraction(0) for _ in range(9)] for _ in range(5)] for _ in range(5)]
        for u, v, middle in product(range(5), repeat=3):
            for h1, h2 in product(range(5), repeat=2):
                region[u][v][h1 + h2] += one[u][middle][h1] * one[middle][v][h2]
        normalized = [[[region[u][v][h] / (comb(8, h) * p ** h * (1 - p) ** (8 - h))
                        for h in range(9)] for v in range(5)] for u in range(5)]
        if conditional is not None:
            assert normalized == conditional
            quota_checks += 5 * 5 * 9
        conditional = normalized
        one_kernel = _am_exact_power(actual_mds_local_exact(F(3, 4), F(2, 3), F(5, 4),
                                      F(1, 2), F(7, 8), p), 1)
        unconditioned_region = _am_matmul(one_kernel, one_kernel)
        quota5 = comb(8, 5) * p ** 5 * (1 - p) ** 3
        upper = [[(unconditioned_region[u][v] - int(u == 0 and v == 0)) / quota5
                  for v in range(5)] for u in range(5)]
        for u in range(5):
            for v in range(5):
                conditioned_positive = normalized[u][v][5] - int(u == 0 and v == 0)
                assert 0 <= conditioned_positive <= upper[u][v]
                regional_checks += 1
        regional_upper.append(upper)
    minimum = [[min(bound[u][v] for bound in regional_upper) for v in range(5)] for u in range(5)]
    assert minimum[0][0] > 0  # Genuine regional excursions must survive.
    conditioned_region = [[conditional[u][v][5] for v in range(5)] for u in range(5)]
    comparison_region = [[minimum[u][v] + int(u == 0 and v == 0) for v in range(5)] for u in range(5)]
    composed_conditional = _am_matmul(conditioned_region, conditioned_region)
    composed_comparison = _am_matmul(comparison_region, comparison_region)
    positive_composition = _am_matmul(minimum, minimum)
    for u in range(5):
        for v in range(5):
            positive_composition[u][v] += (minimum[u][0] if v == 0 else 0) + (minimum[0][v] if u == 0 else 0)
            assert composed_conditional[u][v] <= composed_comparison[u][v]
            assert positive_composition[u][v] == composed_comparison[u][v] - int(u == 0 and v == 0)
            regional_checks += 2
    # Exact nonzero recurrence control; subtraction appears only in this reference test.
    positive = [row[:] for row in actual4]
    positive[0][0] = 0
    full = [row[:] for row in actual4]
    for _ in range(2):
        full = _am_matmul(full, full)
        after = _am_matmul(positive, positive)
        for u in range(5):
            for v in range(5):
                after[u][v] += (positive[u][0] if v == 0 else 0) + (positive[0][v] if u == 0 else 0)
                assert after[u][v] == full[u][v] - int(u == 0 and v == 0)
                nonzero_checks += 1
        positive = after
    # The floating assembly is only compared with the exact toy matrix here.
    s, t, z, nu, eta, p, _ = params
    logs = [np.log2(AM_LD(value.numerator) / value.denominator) for value in (s, t, z, nu, eta)]
    floating, scale = actual_mds_matrix(8, *logs, AM_LD(p.numerator) / p.denominator)
    reference = np.array([[AM_LD(value.numerator) / value.denominator for value in row] for row in actual8])
    restored = floating * np.exp2(scale)
    assert np.max(np.abs(restored - reference) / reference) < AM_LD('1e-12')
    regional_float, regional_scale, _ = actual_mds_regional_kernel(4, 2, 5,
        (AM_LD(1) / 3, AM_LD(2) / 3), *logs)
    regional_reference = np.array([[AM_LD(value.numerator) / value.denominator for value in row] for row in minimum])
    assert np.max(np.abs(regional_float * np.exp2(regional_scale) - regional_reference) / regional_reference) < AM_LD('1e-12')
    explicit_power, explicit_scale = actual_mds_positive_part_power(regional_float, regional_scale, 2)
    positive_reference = np.array([[AM_LD(value.numerator) / value.denominator for value in row] for row in positive_composition])
    assert np.max(np.abs(explicit_power * np.exp2(explicit_scale) - positive_reference) / positive_reference) < AM_LD('1e-12')
    selectors = {}
    diagnosed_region, diagnosed_scale, _ = actual_mds_regional_kernel(4, 2, 5,
        (AM_LD(1) / 3, AM_LD(2) / 3), *logs, selector_data=selectors)
    assert sum(item['entry_count'] for item in selectors['summary']['patterns']) == 25
    shared = max(AM_LD(0), diagnosed_scale)
    diagnosed_matrix = diagnosed_region * np.exp2(diagnosed_scale - shared)
    selectors.update(positive00_scaled=diagnosed_matrix[0, 0], empty00_scaled=np.exp2(-shared))
    diagnosed_matrix[0, 0] += selectors['empty00_scaled']
    diagnosed = actual_mds_collatz(diagnosed_matrix, shared, updates=4, iterations=8,
        diagnostic_top=3, diagnostic_iterations=8, selector_data=selectors)
    assert 'operator_diagnostics' in diagnosed
    assert diagnosed['operator_diagnostics']['regional_minimum_selectors']['total_entries'] == 25
    # Deliberately coarse floors make the majorant effect visible on a toy.
    # This is a FLOAT monotonicity check, separate from all exact comparisons.
    floored, floor_scale = actual_mds_matrix(8, *logs, AM_LD(p.numerator) / p.denominator, floor_bits=2)
    assert np.all(floored * np.exp2(floor_scale) >= reference * (1 - AM_LD('1e-12')))
    region_floor, region_floor_scale, _ = actual_mds_regional_kernel(4, 2, 5,
        (AM_LD(1) / 3, AM_LD(2) / 3), *logs, floor_bits=2)
    assert np.all(region_floor * np.exp2(region_floor_scale) >= regional_reference * (1 - AM_LD('1e-12')))
    floor_power, floor_power_scale = actual_mds_positive_part_power(region_floor, region_floor_scale, 2, floor_bits=2)
    assert np.all(floor_power * np.exp2(floor_power_scale) >= positive_reference * (1 - AM_LD('1e-12')))
    print(dict(exact_local_coefficients=local_checks, exact_matrix_entries=matrix_checks,
               exact_endpoint_checks=endpoint_checks, exact_quota_conditioning_entries=quota_checks,
               exact_nonzero_power_entries=nonzero_checks, floating_matrix_entries=81,
               exact_regional_comparison_checks=regional_checks,
               floating_regional_entries=25, floating_explicit_positive_power_entries=25,
               floating_selector_diagnostic_checks=2, floating_floor_majorant_entries=131,
               full_code_certificate=False), flush=True)


def actual_mds_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    mode.add_argument('--floor-only', action='store_true')
    parser.add_argument('--profile', nargs=5, type=int, default=(4096, 786432, 786432, 4, 1048576),
                        metavar=('H', 'R', 'LS', 'FACTUAL', 'F'))
    parser.add_argument('--exponent', type=float, default=16)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--p', type=float, default=.25)
    parser.add_argument('--log-s', type=float, default=1)
    parser.add_argument('--log-z', type=float, default=0)
    parser.add_argument('--log-nu', type=float, default=0)
    eta = parser.add_mutually_exclusive_group()
    eta.add_argument('--log-eta', type=float)
    eta.add_argument('--eta', type=float)
    parser.add_argument('--iterations', type=int, default=96)
    parser.add_argument('--regional', action='store_true',
                        help='Compose and quota-condition each tilt before the common-basis minimum')
    parser.add_argument('--p-values', nargs=3, type=float, default=(1/5, 1/3, 1/2))
    parser.add_argument('--log-region', type=float, default=0)
    parser.add_argument('--minimum-regions', type=int, default=0)
    parser.add_argument('--diagnostics', action='store_true',
                        help='Bounded left/right comparison-population and entry-selector summaries')
    parser.add_argument('--diagnostic-top', type=int, default=8)
    parser.add_argument('--diagnostic-iterations', type=int, default=32)
    parser.add_argument('--floor-majorant', action='store_true',
                        help='Explicit FLOAT majorant: floor normalized positive entries upward at 2^-4000')
    method = parser.add_mutually_exclusive_group()
    method.add_argument('--full-powers', action='store_true')
    method.add_argument('--nonzero', action='store_true')
    args = parser.parse_args()
    if args.regional and args.floor_only:
        parser.error('--floor-only is the unconditioned single-p control, not a regional bound')
    if args.diagnostics and (args.nonzero or args.full_powers):
        parser.error('--diagnostics currently requires the default Collatz mode')
    if args.diagnostics and not (1 <= args.diagnostic_top <= 32 and 1 <= args.diagnostic_iterations <= 96):
        parser.error('Diagnostics require top in 1..32 and iterations in 1..96')
    if args.eta is not None and args.eta <= 0:
        parser.error('--eta must be positive')
    log_eta = np.log2(AM_LD(args.eta)) if args.eta is not None else args.log_eta
    if args.regional and log_eta is None:
        log_eta = -2
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
            actual_mds_checks()
            return
        selected = 'nonzero' if args.nonzero else 'powers' if args.full_powers else 'collatz'
        diagnostic_top = args.diagnostic_top if args.diagnostics else 0
        floor_bits = 4000 if args.floor_majorant else None
        if args.regional:
            result = actual_mds_regional_profile(tuple(args.profile), args.exponent, args.theta,
                tuple(args.p_values), args.log_s, args.log_z, args.log_nu, log_eta,
                args.iterations, selected, args.log_region, args.minimum_regions, progress=True,
                diagnostic_top=diagnostic_top, diagnostic_iterations=args.diagnostic_iterations,
                floor_bits=floor_bits)
            print(dict(seconds=round(time.monotonic() - started, 3), **result), flush=True)
            return
        floor = actual_mds_full_state_floor(tuple(args.profile), args.exponent, args.theta, args.p,
                                            args.log_s, args.log_z, args.log_nu, log_eta)
        print(dict(stage='analytic_full_state_path_floor', **floor), flush=True)
        if args.floor_only:
            return
        result = actual_mds_profile(tuple(args.profile), args.exponent, args.theta, args.p,
                                    args.log_s, args.log_z, args.log_nu, log_eta,
                                    args.iterations, selected, progress=True,
                                    diagnostic_top=diagnostic_top, diagnostic_iterations=args.diagnostic_iterations,
                                    floor_bits=floor_bits)
        print(dict(seconds=round(time.monotonic() - started, 3), **result), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    actual_mds_main()
