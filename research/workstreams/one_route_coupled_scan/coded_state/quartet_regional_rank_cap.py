"""Positive129-state FLOAT control for regional full-footprint capacity.

Charge the CURRENT update: v_t=|U_t|+|B_t|. The integrated-birth kernel
has z^(4i) on source population4i, and positive selected-quartet polynomial
 (1+w+beta X)^4-1.
Its free endpoint is t_i=(z*beta)^(4i)*(1+w)^(n-4i), including B_T.

For a regional threshold c and eta>=0, the exact regional factor is
2^[-eta*(v_region-c)+]. Positivity bounds its operator by
 min(K(z,w)^ell, 2^(eta*c)*K(z*2^-eta,w*2^-eta)^ell).
Use255 ordinary16-transition regions and a final15-transition/endpoint
vector. These minima are in the common physical population basis.

For an actual profile, C_H<=sum_regions min(h,v_region), hence
 sum_regions(v_region-h)+ <= F-C_H=R-128h.
The fixed threshold c=h therefore covers the whole profile, without a
threshold union. fixed_h_capacity minimizes only over eta choices for
this same event, and includes eta0. An arbitrary c call bounds only the
subclass with sum(v_region-c)+<=the specified excess budget; it is not a
whole-profile statement or a minimum-over-c cover.

The profile correction is the existing capped route factor, minus
G logz+f logw+(L-f)logbeta, plus eta*(R-128h). Births are integrated;
there is no birth binomial, hit correction, or second H union. No C_A
marker or death restriction is used. Arithmetic is FLOAT, not directed.
Load separately from quartet_selector_hit_certificate.py: the helper
reuses the FLOAT selector-hit positive normalization/minimum primitives.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import comb, lgamma, log, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np

try:
    _sh_normalize
except NameError:
    from quartet_selector_hit_moment import _sh_normalize, _sh_minimum


RR_LD = np.longdouble
RR_WIDTH, RR_UPDATES, RR_REGION, RR_M = 512, 4096, 16, 8192
RR_PROFILE = (4096, 819200, 838861, 65536, 1081344)
RR_LOG_BETA = log2(3)
RR_LOG_Z, RR_LOG_W = -1.01-log2(1.5), -1.0


def regional_rank_current_kernel(log_beta=RR_LOG_BETA, log_z=RR_LOG_Z, log_w=RR_LOG_W, width=512):
    """Positive current-volume kernel, supporting any beta,z,w>0."""
    width = index(width)
    log_beta, log_z, log_w = map(RR_LD, (log_beta, log_z, log_w))
    if width < 4 or width % 4 or width > 512 or not np.all(np.isfinite((log_beta, log_z, log_w))):
        raise ValueError('Require quartet width4..512 and finite marker logs')
    groups = width//4
    log_one = np.logaddexp2(RR_LD(0), log_w)
    constant_terms = np.asarray([np.log2(RR_LD(comb(4, ell)))+ell*log_w for ell in range(1, 5)], dtype=RR_LD)
    maximum = constant_terms.max()
    log_constant = maximum+np.log2(np.exp2(constant_terms-maximum).sum())
    local_logs = np.asarray((log_constant, 2+log_beta+3*log_one,
                             np.log2(RR_LD(6))+2*log_beta+2*log_one,
                             2+3*log_beta+log_one, 4*log_beta), dtype=RR_LD)
    local_scale = local_logs.max()
    local = np.exp2(local_logs-local_scale)
    if np.any(local <= 0) or not np.all(np.isfinite(local)):
        raise ArithmeticError('Positive birth polynomial underflow or overflow')
    old = np.arange(0, width+1, 4)
    denominator = np.asarray([np.log2(RR_LD(comb(width, int(u)))) for u in old], dtype=RR_LD)
    matrix_logs = np.empty((groups+1, groups+1), dtype=RR_LD)
    selected, selected_scale = np.ones(1, dtype=RR_LD), RR_LD(0)
    for j in range(groups+1):
        remaining = width-4*j
        unselected = np.asarray([RR_LD(comb(remaining, k)) for k in range(remaining+1)], dtype=RR_LD)
        unselected, unselected_scale = _sh_normalize(unselected)
        coefficients = np.convolve(selected, unselected)[old]
        if np.any(coefficients <= 0) or not np.all(np.isfinite(coefficients)):
            raise ArithmeticError('Invalid positive integrated-birth coefficient')
        matrix_logs[:, j] = (np.log2(coefficients)+selected_scale+unselected_scale-denominator
                              +np.log2(RR_LD(comb(groups, j)))+old*log_z)
        if j < groups:
            selected, extra = _sh_normalize(np.convolve(selected, local))
            selected_scale += local_scale+extra
    scale = matrix_logs.max()
    matrix = np.exp2(matrix_logs-scale)
    if np.any(matrix <= 0) or not np.all(np.isfinite(matrix)):
        raise ArithmeticError('Current integrated-birth matrix underflow or overflow')
    terminal_logs = old*(log_z+log_beta)+(width-old)*log_one
    terminal_scale = terminal_logs.max()
    terminal = np.exp2(terminal_logs-terminal_scale)
    if np.any(terminal <= 0) or not np.all(np.isfinite(terminal)):
        raise ArithmeticError('Current-footprint free endpoint became invalid')
    return dict(kernel=matrix, kernel_scale=scale, terminal=terminal, terminal_scale=terminal_scale,
                width=width, log_beta=log_beta, log_z=log_z, log_w=log_w)


def _rr_region_blocks(kernel_data, length):
    """Share positive square powers between K^length and K^(length-1)*terminal."""
    length = index(length)
    if length < 2:
        raise ValueError('Regions must contain at least two updates')
    powers = [(kernel_data['kernel'], kernel_data['kernel_scale'])]
    while 1 << len(powers) <= length:
        matrix, scale = powers[-1]
        square, extra = _sh_normalize(matrix @ matrix)
        powers.append((square, 2*scale+extra))
    ordinary, ordinary_scale = None, RR_LD(0)
    count, bit = length, 0
    while count:
        if count & 1:
            matrix, scale = powers[bit]
            if ordinary is None:
                ordinary, ordinary_scale = matrix.copy(), scale
            else:
                ordinary, extra = _sh_normalize(ordinary @ matrix)
                ordinary_scale += scale+extra
        count >>= 1
        bit += 1
    terminal, terminal_scale = kernel_data['terminal'].copy(), kernel_data['terminal_scale']
    count, bit = length-1, 0
    while count:
        if count & 1:
            matrix, scale = powers[bit]
            terminal, extra = _sh_normalize(matrix @ terminal)
            terminal_scale += scale+extra
        count >>= 1
        bit += 1
    return ordinary, ordinary_scale, terminal, terminal_scale


def _rr_power_endpoint(matrix, scale, terminal, terminal_scale, count):
    count = index(count)
    if count < 0:
        raise ValueError('Negative region count')
    vector, vector_scale = terminal.copy(), RR_LD(terminal_scale)
    power, power_scale = matrix.copy(), RR_LD(scale)
    while count:
        if count & 1:
            vector, extra = _sh_normalize(power @ vector)
            vector_scale += power_scale+extra
        count >>= 1
        if count:
            power, extra = _sh_normalize(power @ power)
            power_scale = 2*power_scale+extra
    value = vector_scale+np.log2(vector[0])
    if not np.isfinite(value):
        raise ArithmeticError('Non-finite regional footprint moment')
    return value


def regional_rank_cap_prepare(log_beta=RR_LOG_BETA, log_z=RR_LOG_Z, log_w=RR_LOG_W,
                              etas=(0, 1/16), width=512, updates=4096, region_updates=16):
    """Cache ordinary/final block comparisons; marker tilts are independent of c."""
    width, updates, region_updates = map(index, (width, updates, region_updates))
    if region_updates < 2 or updates < region_updates or updates % region_updates:
        raise ValueError('Require a whole number of regions, each with at least two updates')
    etas = tuple(sorted(set(map(float, etas)) | {0.0}))
    if len(etas) > 16 or not all(np.isfinite(eta) and 0 <= eta <= 16 for eta in etas):
        raise ValueError('Require at most16 finite eta values in[0,16]')
    blocks = {}
    for eta in etas:
        kernel = regional_rank_current_kernel(log_beta, log_z-eta, log_w-eta, width)
        blocks[eta] = _rr_region_blocks(kernel, region_updates)
    return dict(blocks=blocks, etas=etas, width=width, updates=updates,
                region_updates=region_updates, regions=updates//region_updates,
                log_beta=float(log_beta), log_z=float(log_z), log_w=float(log_w))


def regional_rank_cap_moment(prepared, threshold, eta):
    """Upper moment with the per-region factor2^[-eta*(v-threshold)+]."""
    threshold, eta = index(threshold), float(eta)
    if not 0 <= threshold <= prepared['width']*prepared['region_updates'] or eta not in prepared['blocks']:
        raise ValueError('Invalid threshold or eta not present in prepared blocks')
    plain, ps, final, fs = prepared['blocks'][0.0]
    if eta:
        marked, ms, marked_final, mfs = prepared['blocks'][eta]
        ordinary, ordinary_scale = _sh_minimum(plain, ps, marked, ms+RR_LD(eta)*threshold)
        endpoint, endpoint_scale = _sh_minimum(final, fs, marked_final, mfs+RR_LD(eta)*threshold)
    else:
        ordinary, ordinary_scale, endpoint, endpoint_scale = plain, ps, final, fs
    value = _rr_power_endpoint(ordinary, ordinary_scale, endpoint, endpoint_scale, prepared['regions']-1)
    return dict(threshold=threshold, eta=eta, log2_moment_upper=float(value),
                ordinary_regions=prepared['regions']-1,
                ordinary_transitions=prepared['region_updates'],
                final_transitions=prepared['region_updates']-1,
                final_births_included=True, current_footprint_charged=True,
                common_physical_population_basis=True, floating_diagnostic_only=True)


def regional_rank_threshold_class(prepared, threshold, excess_budget):
    """Arbitrary-c subclass only; minimize eta bounds for ONE specified event."""
    excess_budget = index(excess_budget)
    if excess_budget < 0:
        return dict(status='empty_threshold_class', full_profile_coverage=False)
    candidates = []
    for eta in prepared['etas']:
        result = regional_rank_cap_moment(prepared, threshold, eta)
        result['log2_threshold_class_moment_upper'] = result['log2_moment_upper']+eta*excess_budget
        candidates.append(result)
    best = min(candidates, key=lambda result: result['log2_threshold_class_moment_upper'])
    return dict(status='fixed_threshold_subclass', threshold=int(threshold), excess_budget=excess_budget,
                candidates=candidates, best=best, full_profile_coverage=False,
                minimum_over_different_thresholds_is_not_a_cover=True)


def _rr_logchoose(n, k):
    if not 0 <= k <= n:
        return float('-inf')
    if k in (0, n):
        return 0.0
    return (lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/log(2)


def fixed_h_capacity(profile=RR_PROFILE, log_beta=RR_LOG_BETA, log_z=RR_LOG_Z, log_w=RR_LOG_W,
                     etas=(0, 1/16), prepared=None):
    """Whole-profile geometry control at the justified fixed threshold c=h."""
    h, rank, volume, births, footprint = map(index, profile)
    N = RR_WIDTH*RR_UPDATES
    if not (1 <= h <= RR_M and 128*h <= rank <= 1 << 20 and 1 <= births <= volume <= footprint <= N):
        raise ValueError('Invalid actual-state profile')
    G, m, CH = footprint-births, volume-births, footprint-rank+128*h
    if G % 4 or G > N-RR_WIDTH or not np.isfinite(log_beta) or log_beta < 0:
        raise ValueError('Invalid live-selector volume or overlap marker beta<1')
    if not (129*h <= CH <= 256*h and 0 <= CH-births <= G and 0 <= 256*h-CH <= N-footprint):
        return dict(profile=tuple(profile), status='empty_profile', full_code_certificate=False)
    if prepared is None:
        prepared = regional_rank_cap_prepare(log_beta, log_z, log_w, etas)
    if any(prepared[key] != value for key, value in (
            ('width', RR_WIDTH), ('updates', RR_UPDATES), ('region_updates', RR_REGION),
            ('log_beta', float(log_beta)), ('log_z', float(log_z)), ('log_w', float(log_w)))):
        raise ValueError('Prepared parameters do not match the production profile markers')
    route = min(0.0, _rr_logchoose(G, CH-births)+_rr_logchoose(N-footprint, 256*h-CH)-255*_rr_logchoose(RR_M, h))
    common = route-G*log_z-births*log_w-m*log_beta
    result = regional_rank_threshold_class(prepared, h, rank-128*h)
    for candidate in result['candidates']:
        candidate['raw_log2_probability_bound'] = candidate['log2_threshold_class_moment_upper']+common
    result.update(profile=tuple(map(int, profile)), status='whole_profile_fixed_h_capacity',
                  full_profile_coverage=True, covers_all_deaths=True, allowed_run_closure_retained=False,
                  route_log2=route, common_log2_profile_correction=common,
                  raw_log2_probability_bound=result['best']['raw_log2_probability_bound'],
                  threshold_union_required=False, H_union_already_in_route=True,
                  births_integrated=True, full_code_certificate=False)
    return result


def _rr_partitions(width):
    if width == 4:
        return [(15,)]
    return [(block, 255^block) for rest in combinations(range(1, 8), 3)
            for block in (1+sum(1 << j for j in rest),)]


def _rr_subsets(mask):
    subset = mask
    while True:
        yield subset
        if subset == 0:
            break
        subset = (subset-1)&mask


def _rr_exact_kernel(width, beta, z, w):
    size, partitions = width//4+1, _rr_partitions(width)
    kernel = [[{} for _ in range(size)] for _ in range(size)]
    terminal = [{} for _ in range(size)]
    for i in range(size):
        u, old = 4*i, (1 << (4*i))-1
        for partition in partitions:
            for choice in range(1 << len(partition)):
                blocks = [block for j, block in enumerate(partition) if choice >> j & 1]
                new = sum(blocks)
                for births in _rr_subsets(new&~old):
                    if any(not (block&(old|births)) for block in blocks):
                        continue
                    v = u+births.bit_count()
                    weight = z**u*beta**(old&new).bit_count()*w**births.bit_count()/len(partitions)
                    polynomial = kernel[i][len(blocks)]
                    polynomial[v] = polynomial.get(v, Fraction(0))+weight
        for count in range(width-u+1):
            terminal[i][u+count] = z**u*beta**u*comb(width-u, count)*w**count
    return kernel, terminal


def _rr_poly_add_product(target, first, second):
    for i, left in first.items():
        for j, right in second.items():
            target[i+j] = target.get(i+j, Fraction(0))+left*right


def _rr_exact_region(kernel, terminal, length):
    size = len(kernel)
    power = [[{0: Fraction(int(i == j))} if i == j else {} for j in range(size)] for i in range(size)]
    final = None
    for step in range(length+1):
        if step == length-1:
            final = [{} for _ in range(size)]
            for i in range(size):
                for j in range(size):
                    _rr_poly_add_product(final[i], power[i][j], terminal[j])
        if step == length:
            break
        following = [[{} for _ in range(size)] for _ in range(size)]
        for i in range(size):
            for j in range(size):
                for k in range(size):
                    _rr_poly_add_product(following[i][j], power[i][k], kernel[k][j])
        power = following
    return power, final


def regional_rank_cap_checks():
    kernels = endpoints = regional = moments = capacity = 0
    beta, z, w = Fraction(3, 2), Fraction(2, 3), Fraction(1, 4)
    for width in (4, 8):
        exact, terminal = _rr_exact_kernel(width, beta, z, w)
        numeric = regional_rank_current_kernel(log2(float(beta)), log2(float(z)), log2(float(w)), width)
        for i in range(width//4+1):
            for j in range(width//4+1):
                value = sum(exact[i][j].values())
                observed = np.log2(numeric['kernel'][i, j])+numeric['kernel_scale']
                assert abs(float(observed)-log2(float(value))) < 1e-10
                kernels += 1
            value = sum(terminal[i].values())
            observed = np.log2(numeric['terminal'][i])+numeric['terminal_scale']
            assert abs(float(observed)-log2(float(value))) < 1e-10
            endpoints += 1
        blocks, final = _rr_exact_region(exact, terminal, 2)
        prepared = regional_rank_cap_prepare(log2(float(beta)), log2(float(z)), log2(float(w)),
                                             (0, 1), width=width, updates=4, region_updates=2)
        for threshold in (0, width//2, width, 2*width):
            filtered = [[sum(value*Fraction(1, 2)**max(v-threshold, 0) for v, value in entry.items())
                         for entry in row] for row in blocks]
            ending = [sum(value*Fraction(1, 2)**max(v-threshold, 0) for v, value in entry.items()) for entry in final]
            upper = [[min(sum(entry.values()), sum(value*Fraction(1, 2)**(v-threshold) for v, value in entry.items()))
                      for entry in row] for row in blocks]
            end_upper = [min(sum(entry.values()), sum(value*Fraction(1, 2)**(v-threshold) for v, value in entry.items()))
                         for entry in final]
            for i in range(len(blocks)):
                for j in range(len(blocks)):
                    assert filtered[i][j] <= upper[i][j]
                    regional += 1
                assert ending[i] <= end_upper[i]
                regional += 1
            exact_moment = sum(filtered[0][j]*ending[j] for j in range(len(blocks)))
            comparison_moment = sum(upper[0][j]*end_upper[j] for j in range(len(blocks)))
            numeric_moment = regional_rank_cap_moment(prepared, threshold, 1)['log2_moment_upper']
            assert exact_moment <= comparison_moment
            assert abs(numeric_moment-log2(float(comparison_moment))) < 1e-9
            moments += 1
            for budget in (0, width, 2*width):
                restricted = Fraction(0)
                for j in range(len(blocks)):
                    for v1, first in blocks[0][j].items():
                        for v2, last in final[j].items():
                            if max(v1-threshold, 0)+max(v2-threshold, 0) <= budget:
                                restricted += first*last
                assert restricted <= 2**budget*exact_moment
                moments += 1
        neutral = regional_rank_cap_moment(prepared, width, 0)['log2_moment_upper']
        unmarked = sum(sum(blocks[0][j].values())*sum(final[j].values()) for j in range(len(blocks)))
        assert abs(neutral-log2(float(unmarked))) < 1e-9
        moments += 1
    for h in range(9):
        for first in range(9):
            for second in range(9):
                maximum_H = min(h, first)+min(h, second)
                for covered_H in range(maximum_H+1):
                    assert max(first-h, 0)+max(second-h, 0) <= first+second-covered_H
                    capacity += 1
    print(dict(exact_current_kernel_entries=kernels, exact_final_birth_endpoints=endpoints,
               exact_regional_capacity_comparisons=regional, exact_two_region_moments=moments,
               exact_fixed_h_capacity_identities=capacity, no_threshold_union=True,
               full_code_certificate=False), flush=True)


def regional_rank_cap_main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--checks', action='store_true')
    modes.add_argument('--probe', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=RR_PROFILE)
    parser.add_argument('--log-beta', type=float, default=RR_LOG_BETA)
    parser.add_argument('--log-z', type=float, default=RR_LOG_Z)
    parser.add_argument('--log-w', type=float, default=RR_LOG_W)
    parser.add_argument('--etas', type=float, nargs='+', default=(0, 1/16))
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            regional_rank_cap_checks()
        else:
            print(fixed_h_capacity(args.profile, args.log_beta, args.log_z, args.log_w, args.etas), flush=True)
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    regional_rank_cap_main()
