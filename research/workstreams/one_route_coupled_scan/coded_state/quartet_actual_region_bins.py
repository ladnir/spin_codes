"""Regional coverage bins before the actual-state eligible-quota minimum.

Fix an actual-support history before averaging its quartet partitions and
eligible-slot labels.  In one region F=G+f lies in 0..width*region_updates.
Replacing z by z*beta and nu by nu*beta marks that regional F by beta^F.
For a nonempty coverage bin I=[a,b], its conditional operator is bounded by

    beta^(-a) Q_p(beta)/pi_h(p), beta >= 1;
    beta^(-b) Q_p(beta)/pi_h(p), beta <= 1,

where Q_p(beta)=K_p(z*beta,nu*beta)^region_updates-Z and Z=e0 e0^T.
All candidates act in the SAME physical-popcount basis.  Minimize entries
over (p,beta) separately for each bin, then sum the disjoint bins.  Finally
take the entrywise minimum with the unsplit beta=1 comparison.  Add Z only
after these operations.  Nonempty regional (0,0) excursions are retained.
The endpoint vector is all ones; there is no terminal reset or expansion.

This ordering retains regional coverage information before the quota
minimum.  It does not extract the exact global F coefficient.  The base
global z,nu markers and fixed-profile correction are still needed.  No
extra bin union factor appears.  This is a FLOAT diagnostic, not a code
certificate.  Optional positive floors have the existing helper's FLOAT
majorant semantics, without directed-rounding guarantees.

Only this new file is owned by this control.  The existing MDS helper can
be supplied as a module/dictionary under actual_mds_helper, or imported.
All numerical jobs belong to the parent runner.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time
from types import SimpleNamespace

import numpy as np


RB_LD = np.longdouble
RB_PROFILE = (4096, 786432, 786432, 4, 1048576)


def _rb_helper(helper=None):
    if helper is None:
        helper = globals().get('actual_mds_helper')
    if helper is None:
        import quartet_actual_mds_moment as helper
    return SimpleNamespace(**helper) if isinstance(helper, dict) else helper


def _rb_bins(slots, bins=None, bin_width=256):
    slots, bin_width = index(slots), index(bin_width)
    if slots < 1 or bin_width < 1:
        raise ValueError('Require positive slots and bin width')
    if bins is None:
        bins = tuple((lo, min(slots, lo + bin_width - 1))
                     for lo in range(1, slots + 1, bin_width))
    else:
        bins = tuple((index(lo), index(hi)) for lo, hi in bins)
    if not bins or bins[0][0] != 1 or bins[-1][1] != slots or any(
            lo > hi for lo, hi in bins) or any(
            bins[j][0] != bins[j - 1][1] + 1 for j in range(1, len(bins))):
        raise ValueError('Coverage bins must partition every integer from 1 through slots')
    return bins


def _rb_logs(matrix, scale):
    matrix, scale = np.asarray(matrix, dtype=RB_LD), RB_LD(scale)
    if not np.isfinite(scale) or not np.all(np.isfinite(matrix)) or np.any(matrix < 0):
        raise ArithmeticError('Invalid scaled positive comparison')
    logs = np.full_like(matrix, RB_LD('-inf'))
    np.log2(matrix, out=logs, where=matrix > 0)
    logs += scale
    return logs


def _rb_materialize(logs, helper, floor_bits):
    if floor_bits is not None:
        return helper._am_logs_majorant(logs, floor_bits)
    selected = np.isfinite(logs)
    if not np.any(selected) or np.any(np.isnan(logs)) or np.any(np.isposinf(logs)):
        raise ArithmeticError('Empty or invalid regional comparison')
    scale = logs[selected].max()
    matrix = np.zeros_like(logs)
    matrix[selected] = np.exp2(logs[selected] - scale)
    if np.any(matrix[selected] <= 0) or not np.all(np.isfinite(matrix)):
        raise ArithmeticError('Regional comparison rescaling lost a positive entry; use explicit floors')
    return matrix, scale


class ActualRegionBinComparison:
    """Reusable log-domain bin update for quota-adjusted, scaled Q matrices.

    add() receives a nonempty regional upper matrix already divided by its
    quota probability and already carrying beta^F.  This also permits a
    future exact-birth coefficient caller: include beta^j in degree j's
    scale when only z was tilted.  This class does not build kernels,
    change birth degrees, add Z, or charge global profile corrections.
    """

    def __init__(self, size, slots, bins=None, bin_width=256):
        self.size = index(size)
        if self.size < 1:
            raise ValueError('Require a positive state count')
        self.bins = _rb_bins(slots, bins, bin_width)
        shape = self.size, self.size
        self.bounds = np.full((len(self.bins), *shape), RB_LD('inf'), dtype=RB_LD)
        self.unsplit = np.full(shape, RB_LD('inf'), dtype=RB_LD)
        self.scratch = np.empty(shape, dtype=RB_LD)
        self.candidates = self.unsplit_candidates = 0

    def add(self, positive, scale, log_beta):
        log_beta = RB_LD(log_beta)
        if not np.isfinite(log_beta) or np.shape(positive) != (self.size, self.size):
            raise ValueError('Invalid coverage tilt or comparison shape')
        logs = _rb_logs(positive, scale)
        for bound, (lo, hi) in zip(self.bounds, self.bins):
            endpoint = lo if log_beta >= 0 else hi
            np.subtract(logs, endpoint * log_beta, out=self.scratch)
            np.minimum(bound, self.scratch, out=bound)
        if log_beta == 0:
            np.minimum(self.unsplit, logs, out=self.unsplit)
            self.unsplit_candidates += 1
        self.candidates += 1

    def finish(self, helper=None, floor_bits=None, comparison_data=None):
        helper = _rb_helper(helper)
        floor_bits = helper._am_floor_bits(floor_bits)
        if not self.candidates or not self.unsplit_candidates:
            raise ValueError('At least one beta=1 candidate is required for the unsplit control')
        summed = self.bounds[0].copy()
        for bound in self.bounds[1:]:
            np.logaddexp2(summed, bound, out=summed)
        np.minimum(summed, self.unsplit, out=summed)
        positive, scale = _rb_materialize(summed, helper, floor_bits)
        baseline, baseline_scale = _rb_materialize(self.unsplit, helper, floor_bits)
        # Floors are monotone here: x<=y implies max(x,2^-b max(x)) <=
        # max(y,2^-b max(y)).  Thus final materialization preserves the cap.
        represented = _rb_logs(positive, scale)
        baseline_logs = _rb_logs(baseline, baseline_scale)
        supported = np.isfinite(baseline_logs)
        if not np.all(np.isfinite(represented[supported])):
            raise ArithmeticError('A supported regional comparison entry disappeared')
        gains = baseline_logs[supported] - represented[supported]
        tolerance = RB_LD('1e-10')
        if np.any(gains < -tolerance):
            raise ArithmeticError('The binned comparison exceeded its unsplit cap')
        details = dict(coverage_bins=self.bins, candidates=self.candidates,
            unsplit_candidates=self.unsplit_candidates,
            bins_partition_nonempty_coverage=True, empty_operator_kept_separate=True,
            entrywise_unsplit_cap=True, entrywise_gain_min_bits=float(gains.min()),
            entrywise_gain_max_bits=float(gains.max()),
            entrywise_gain_mean_bits=float(gains.mean()),
            improved_entries=int(np.count_nonzero(gains > tolerance)),
            total_entries=int(gains.size), positive_floor_bits=floor_bits,
            floating_diagnostic_only=True)
        if np.isfinite(baseline_logs[0, 0]):
            details['positive_00_gain_bits'] = float(baseline_logs[0, 0] - represented[0, 0])
        if comparison_data is not None:
            comparison_data.update(unsplit_positive=baseline, unsplit_scale=baseline_scale)
        return positive, scale, details


def actual_region_binned_kernel(width=512, region_updates=16, eligible=4096,
        p_values=(1/5, 1/3, 1/2), log_s=1, log_t=-1, log_z=0, log_nu=-18,
        log_eta=-2, coverage_bins=None, bin_width=256, log_betas=(-1, 0, 1),
        progress=False, floor_bits=None, comparison_data=None, helper=None):
    """Nine positive regional builds by default; bin updates reuse those builds."""
    helper = _rb_helper(helper)
    width, region_updates, eligible = map(index, (width, region_updates, eligible))
    slots = width * region_updates
    if width < 4 or width % 4 or region_updates < 1 or region_updates & (region_updates - 1):
        raise ValueError('Require quartet width and a positive power-of-two region length')
    if not 0 <= eligible <= slots:
        raise ValueError('Unsupported eligible quota')
    p_values = tuple(RB_LD(p) for p in p_values)
    log_betas = tuple(RB_LD(b) for b in log_betas)
    if not p_values or any(not 0 < p < 1 for p in p_values):
        raise ValueError('Every p must lie strictly between zero and one')
    if not log_betas or not all(np.isfinite(b) for b in log_betas) or 0 not in log_betas:
        raise ValueError('Finite coverage tilts must include log_beta=0')
    floor_bits = helper._am_floor_bits(floor_bits)
    comparison = ActualRegionBinComparison(width + 1, slots, coverage_bins, bin_width)
    choose = RB_LD(helper._am_logchoose(slots, eligible))
    builds = []
    for log_beta in log_betas:
        for p in p_values:
            started = time.monotonic()
            if progress:
                print(dict(stage='actual_region_coverage_tilt', p=float(p),
                           log_beta=float(log_beta)), flush=True)
            kernel, scale = helper.actual_mds_matrix(width, log_s, log_t,
                RB_LD(log_z) + log_beta, RB_LD(log_nu) + log_beta,
                log_eta, p, progress, floor_bits)
            kernel[0, 0] = 0  # One-update K00=1; NEVER repeat this after composition.
            positive, scale = helper.actual_mds_positive_part_power(
                kernel, scale, region_updates, floor_bits)
            quota = choose + eligible * np.log2(p) + (slots - eligible) * np.log2(1 - p)
            scale -= quota
            comparison.add(positive, scale, log_beta)
            builds.append(dict(p=float(p), log_beta=float(log_beta),
                log2_quota=float(quota), positive_scale=float(scale),
                seconds=round(time.monotonic() - started, 3)))
    positive, scale, details = comparison.finish(helper, floor_bits, comparison_data)
    if region_updates > 1 and positive[0, 0] <= 0:
        raise ArithmeticError('Regional nonempty return-to-zero histories were lost')
    details['regional_builds'] = builds
    return positive, scale, details


def _rb_moment(positive, scale, helper, method, iterations, floor_bits):
    if method == 'nonzero':
        powered, power_scale = helper.actual_mds_positive_part_power(positive, scale, 256, floor_bits)
        return dict(log2_moment=float(power_scale + np.log2(powered[0].sum())),
                    method='regional_positive_nonzero_powers', includes_empty_history=False)
    logs = _rb_logs(positive, scale)
    logs[0, 0] = np.logaddexp2(logs[0, 0], RB_LD(0))
    regional, regional_scale = _rb_materialize(logs, helper, floor_bits)
    if method == 'collatz':
        return helper.actual_mds_collatz(regional, regional_scale, 256, iterations)
    if method == 'powers':
        return helper.actual_mds_power_moment(regional, regional_scale, 256, False, floor_bits)
    raise ValueError('Unknown regional evaluation method')


def actual_region_binned_profile(profile=RB_PROFILE, exponent=16, theta=1/16,
        p_values=(1/5, 1/3, 1/2), log_s=1, log_z=0, log_nu=-18, log_eta=-2,
        bin_width=256, log_betas=(-1, 0, 1), iterations=96, method='collatz',
        log_region=0, minimum_regions=0, progress=False, floor_bits=None, helper=None):
    """Fixed-profile FLOAT control, including the unsplit evaluation for comparison."""
    helper = _rb_helper(helper)
    h, R, Ls, factual, F = map(index, profile)
    minimum_regions = index(minimum_regions)
    if not (1 <= h <= 8192 and 128 * h <= R <= 1 << 20
            and 1 <= factual <= Ls <= F <= 1 << 21 and (F - factual) % 4 == 0):
        raise ValueError('Invalid actual-support profile')
    C_H, outside, G = F - R + 128 * h, R - 128 * h, F - factual
    if not (129 * h <= C_H <= 256 * h and factual <= C_H and outside <= 128 * (8192 - h)):
        raise ValueError('Invalid covered-row support or forced births')
    rho = max(129, (C_H + h - 1) // h)
    if h < 8192:
        rho = max(rho, (outside + 8192 - h - 1) // (8192 - h))
    exponent, theta, log_s, log_z, log_nu, log_eta, log_region = map(RB_LD,
        (exponent, theta, log_s, log_z, log_nu, log_eta, log_region))
    if not (exponent >= 0 and 0 <= theta <= 1 and log_region >= 0
            and 0 <= minimum_regions <= rho and all(np.isfinite(v) for v in
            (exponent, theta, log_s, log_z, log_nu, log_eta, log_region))):
        raise ValueError('Invalid markers or unsupported occupied-region threshold')
    log_lambda = np.logaddexp2(RB_LD(0), 32 - exponent)
    log_t = theta * (32 - 2 * exponent - log_lambda)
    control = {}
    positive, scale, comparison = actual_region_binned_kernel(512, 16, h, p_values,
        log_s, log_t, log_z, log_nu, log_eta, bin_width=bin_width, log_betas=log_betas,
        progress=progress, floor_bits=floor_bits, comparison_data=control, helper=helper)
    scale += log_region
    binned = _rb_moment(positive, scale, helper, method, iterations, floor_bits)
    baseline = _rb_moment(control['unsplit_positive'], control['unsplit_scale'] + log_region,
                          helper, method, iterations, floor_bits)
    # Finite Collatz iteration can choose different witnesses.  Both are
    # valid upper diagnostics; retain their minimum as the final control.
    moment = min(binned['log2_moment'], baseline['log2_moment'])
    eligible_union = RB_LD(helper._am_logchoose(8192, h))
    field = theta * (-31 - 32 * R + 943718 * exponent + Ls * log_lambda)
    marker = -log_s * Ls - log_z * G - log_nu * factual - log_eta * outside - log_region * minimum_regions
    prefactor = eligible_union + field + marker
    return dict(statistics='actual_state_support', profile=tuple(profile),
        comparison=comparison, binned_moment=binned, unsplit_moment=baseline,
        log2_moment=float(moment), binned_value=float(prefactor + binned['log2_moment']),
        unsplit_value=float(prefactor + baseline['log2_moment']), value=float(prefactor + moment),
        overall_gain_bits=float(baseline['log2_moment'] - moment),
        combined_prefactor=float(prefactor), eligible_union=float(eligible_union),
        field_prefactor=float(field), marker_correction=float(marker),
        regional_positive_00_log2=float(scale + np.log2(positive[0, 0])),
        log_s=float(log_s), log_t=float(log_t), log_z=float(log_z),
        log_nu=float(log_nu), log_eta=float(log_eta), exponent=float(exponent), theta=float(theta),
        log_region=float(log_region), minimum_occupied_regions=minimum_regions,
        external_quota_divisor_used=False, separate_routing_factor_used=False,
        exact_global_coverage_extracted=False, extra_bin_union_factor_used=False,
        positive_floor_bits=floor_bits, floating_diagnostic_only=True, full_code_certificate=False)


def _rb_toy_histogram(s, t, z, nu, eta):
    """Independent width-four physical masks, with exact H and F counters."""
    masks = [[{} for _ in range(16)] for _ in range(16)]
    for old in range(16):
        r = old.bit_count()
        for new in range(16):
            covered = 15 if old else new
            base = s ** new.bit_count() * t ** (old & ~new).bit_count()
            base *= z ** (4 if old else 0) * nu ** (0 if old else new.bit_count())
            entries = masks[old][new]
            for H in range(16):
                if (not old and new & ~H) or (old and (15 & ~H & ~new).bit_count() > r - 1):
                    continue
                key = covered.bit_count(), H.bit_count()
                entries[key] = entries.get(key, Fraction(0)) + base * eta ** (covered & ~H).bit_count()
    region = {}
    for u in range(5):
        old = (1 << u) - 1
        for middle in range(16):
            for new in range(16):
                for (F1, h1), a in masks[old][middle].items():
                    for (F2, h2), b in masks[middle][new].items():
                        key = u, new.bit_count(), F1 + F2, h1 + h2
                        region[key] = region.get(key, Fraction(0)) + a * b
    return region


def actual_region_bins_checks(helper=None):
    helper = _rb_helper(helper)
    Q = Fraction
    s, t, z, nu, eta = Q(3, 4), Q(2, 3), Q(5, 4), Q(1, 2), Q(7, 8)
    bins = ((1, 2), (3, 4), (5, 6), (7, 8))
    histogram = _rb_toy_histogram(s, t, z, nu, eta)
    quota = 5
    true_bins = [[[sum(histogram.get((u, v, F, quota), Q(0)) for F in range(lo, hi + 1)) / comb(8, quota)
                   for v in range(5)] for u in range(5)] for lo, hi in bins]
    candidates, baseline = [[] for _ in bins], []
    checks = 0
    for beta in (Q(1, 2), Q(1), Q(2)):
        for p in (Q(1, 3), Q(2, 3)):
            kernel = helper._am_exact_power(helper.actual_mds_local_exact(s, t, z * beta, nu * beta, eta, p), 1)
            square = helper._am_matmul(kernel, kernel)
            probability = comb(8, quota) * p ** quota * (1 - p) ** (8 - quota)
            upper = [[(square[u][v] - int(u == 0 and v == 0)) / probability
                      for v in range(5)] for u in range(5)]
            if beta == 1:
                baseline.append(upper)
            for j, (lo, hi) in enumerate(bins):
                factor = beta ** (-(lo if beta >= 1 else hi))
                candidate = [[factor * value for value in row] for row in upper]
                candidates[j].append(candidate)
                for u in range(5):
                    for v in range(5):
                        assert true_bins[j][u][v] <= candidate[u][v]
                        checks += 1
    unsplit = [[min(a[u][v] for a in baseline) for v in range(5)] for u in range(5)]
    result = [[min(unsplit[u][v], sum(min(a[u][v] for a in group) for group in candidates))
               for v in range(5)] for u in range(5)]
    true_positive = [[sum(a[u][v] for a in true_bins) for v in range(5)] for u in range(5)]
    for u in range(5):
        for v in range(5):
            assert true_positive[u][v] <= result[u][v] <= unsplit[u][v]
            checks += 1
    assert histogram[(0, 0, 0, quota)] == comb(8, quota)
    assert true_positive[0][0] > 0 and result[0][0] > 0
    physical = [[true_positive[u][v] + int(u == 0 and v == 0) for v in range(5)] for u in range(5)]
    upper = [[result[u][v] + int(u == 0 and v == 0) for v in range(5)] for u in range(5)]
    physical_power, upper_power = helper._am_matmul(physical, physical), helper._am_matmul(upper, upper)
    assert sum(physical_power[0]) <= sum(upper_power[0])  # Free terminal, two regions.
    logs = [np.log2(RB_LD(a.numerator) / a.denominator) for a in (s, t, z, nu, eta)]
    data = {}
    floating, scale, info = actual_region_binned_kernel(4, 2, 5,
        (RB_LD(1) / 3, RB_LD(2) / 3), *logs, coverage_bins=bins,
        log_betas=(-1, 0, 1), comparison_data=data, helper=helper)
    restored = floating * np.exp2(scale)
    reference = np.array([[RB_LD(a.numerator) / a.denominator for a in row] for row in result])
    assert np.max(np.abs(restored - reference) / reference) < RB_LD('1e-11')
    assert info['candidates'] == 6 and info['unsplit_candidates'] == 2
    floored, fscale, _ = actual_region_binned_kernel(4, 2, 5,
        (RB_LD(1) / 3, RB_LD(2) / 3), *logs, coverage_bins=bins,
        log_betas=(-1, 0, 1), comparison_data=data, helper=helper, floor_bits=2)
    true_float = np.array([[RB_LD(a.numerator) / a.denominator for a in row] for row in true_positive])
    assert np.all(floored * np.exp2(fscale) >= true_float * (1 - RB_LD('1e-11')))
    assert np.all(floored * np.exp2(fscale) <= data['unsplit_positive'] * np.exp2(data['unsplit_scale']) * (1 + RB_LD('1e-11')))
    print(dict(exact_bin_comparison_entries=checks, exact_empty_and_excursion_checks=2,
        exact_free_endpoint_checks=1, floating_exact_comparison_entries=25,
        floating_floor_comparison_entries=50, full_code_certificate=False), flush=True)


def actual_region_bins_main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--checks', action='store_true')
    modes.add_argument('--probe', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=RB_PROFILE)
    parser.add_argument('--exponent', type=float, default=16)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--log-s', type=float, default=1)
    parser.add_argument('--log-z', type=float, default=0)
    parser.add_argument('--log-nu', type=float, default=-18)
    parser.add_argument('--log-eta', type=float, default=-2)
    parser.add_argument('--p', action='append', type=float)
    parser.add_argument('--log-beta', action='append', type=float)
    parser.add_argument('--bin-width', type=int, default=256)
    parser.add_argument('--iterations', type=int, default=96)
    parser.add_argument('--method', choices=('collatz', 'powers', 'nonzero'), default='collatz')
    parser.add_argument('--log-region', type=float, default=0)
    parser.add_argument('--minimum-regions', type=int, default=0)
    parser.add_argument('--floor-majorant', action='store_true')
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
            actual_region_bins_checks()
        else:
            print(actual_region_binned_profile(tuple(args.profile), args.exponent, args.theta,
                args.p or (1/5, 1/3, 1/2), args.log_s, args.log_z, args.log_nu, args.log_eta,
                args.bin_width, args.log_beta or (-1, 0, 1), args.iterations, args.method,
                args.log_region, args.minimum_regions, progress=True,
                floor_bits=4000 if args.floor_majorant else None), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    actual_region_bins_main()
