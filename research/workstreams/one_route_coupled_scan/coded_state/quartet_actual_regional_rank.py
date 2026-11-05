r"""Actual-support population moment with regional full-footprint capacity.

For old/new support populations a,b, the positive local polynomial is

 P(X,Y)=(1+w*v*Y)^4
       +z^4 sum_(r=1)^4 C(4,r) X^r (d+v*Y)^r (1+v*Y)^(4-r),
 K_ab=[X^a Y^b]P^(n/4)/C(n,a).

Every update charges CURRENT G=|Gamma(Sprev)|, f=|Snew\Gamma(Sprev)|,
L=|Snew|, and D=|Sprev\Snew|. All old/new supports and exact live hits
are retained. This counts unsaturated supports too; it is a separate
upper control, not a claim to improve every saturated-support bound.

For a fixed profile(h,R,f,F), the route probability is at most Q_route
uniformly over its fixed support/partition histories. That event also
requires sum_regions(v_region-h)+<=Cout=R-128h, where v_region=G_r+f_r.
For eta>=0 the exact penalized regional population operator is bounded by

 B_eta=min(K(z,w)^16, 2^(eta*h) K(z*2^-eta,w*2^-eta)^16)

entrywise AFTER composition in the same physical population basis.
The moment is e0 B_eta^256 1: all4096 updates are charged, without a
terminal reset or an extra endpoint update. Its profile correction is

 log2 Q_route-G log2 z-f log2 w-L0 log2 v-D0 log2 d+eta*Cout,

with v>=1 and0<d<=1. Defaults L0=R+1,D0=W-R-1 cover the entire stated
dimension-obstruction class after the separately proved birth pruning.
Q_route already includes the H union. No additional quota divisor,
birth-choice factor, amplitude count, or field factor is used.

All numerical bounds are FLOAT diagnostics. Reuse the existing positive
power and validation primitives; the only new hot loop assembles this
different25-coefficient polynomial. No optimizer or numerical jobs on
import. Parent runs all tests/probes serially, with the standalone guard.
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

for _ar_thread_var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_ar_thread_var] = '1'

import numpy as np


AR_LD = np.longdouble
AR_WIDTH, AR_UPDATES, AR_REGION, AR_M, AR_W = 512, 4096, 16, 8192, 943718
AR_N = AR_WIDTH*AR_UPDATES
AR_PROFILE = (4096, 819200, 65536, 1081344)
AR_LOGS = (-.8163, -2.35998, .83827, -1.25818)  # z,w,v,d; inherited scalar control only


def _ar_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_actual_mds_moment as module
    return getattr(module, name)


def _ar_validate_logs(log_z, log_w, log_v, log_d):
    values = tuple(map(AR_LD, (log_z, log_w, log_v, log_d)))
    if not np.all(np.isfinite(values)):
        raise ValueError('All four marker logs must be finite')
    return values


def actual_regional_rank_local_exact(z, w, v, d):
    """Independent rational coefficient formula; no eligible labels here."""
    z, w, v, d = map(Fraction, (z, w, v, d))
    if min(z, w, v, d) < 0:
        raise ValueError('Exact toy markers must be nonnegative')
    table = [[Fraction(0) for _ in range(5)] for _ in range(5)]
    for b in range(5):
        table[0][b] = comb(4, b)*(w*v)**b
    for r in range(1, 5):
        for b in range(5):
            table[r][b] = comb(4, r)*z**4*v**b*sum(
                Fraction(comb(r, kept)*comb(4-r, b-kept))*d**(r-kept)
                for kept in range(max(0, b-(4-r)), min(r, b)+1))
    return table


def _ar_local_scaled(log_z, log_w, log_v, log_d):
    log_z, log_w, log_v, log_d = _ar_validate_logs(log_z, log_w, log_v, log_d)
    logsum = _ar_dependency('_am_logsum')
    logs = np.empty((5, 5), dtype=AR_LD)
    for b in range(5):
        logs[0, b] = np.log2(AR_LD(comb(4, b)))+b*(log_w+log_v)
    for r in range(1, 5):
        for b in range(5):
            terms = [np.log2(AR_LD(comb(r, k)*comb(4-r, b-k)))+(r-k)*log_d
                     for k in range(max(0, b-(4-r)), min(r, b)+1)]
            logs[r, b] = np.log2(AR_LD(comb(4, r)))+4*log_z+b*log_v+logsum(terms)
    scale = logs.max()
    table = np.exp2(logs-scale)
    _ar_dependency('_am_positive')(table, 'actual regional local polynomial')
    return table, scale


def actual_regional_rank_matrix(width=512, log_z=AR_LOGS[0], log_w=AR_LOGS[1],
                                log_v=AR_LOGS[2], log_d=AR_LOGS[3], progress=False):
    """Return normalized K and its log2 scale using three fixed buffers."""
    width = index(width)
    if width < 4 or width > 512 or width % 4:
        raise ValueError('Width must be a multiple of four in4..512')
    positive = _ar_dependency('_am_positive')
    local, local_scale = _ar_local_scaled(log_z, log_w, log_v, log_d)
    current = np.zeros((width+1, width+1), dtype=AR_LD)
    after, scratch = np.empty_like(current), np.empty_like(current)
    current[0, 0] = 1
    scale, local_min = AR_LD(0), local.min()
    started = time.monotonic()
    for group in range(width//4):
        size, target_size = 4*group+1, 4*group+5
        source = current[:size, :size]
        target, work = after[:target_size, :target_size], scratch[:size, :size]
        target.fill(0)
        if source.min()*local_min <= 0:
            raise ArithmeticError('A positive polynomial product would underflow')
        for r in range(5):
            for b in range(5):
                np.multiply(source, local[r, b], out=work)
                np.add(after[r:r+size, b:b+size], work, out=after[r:r+size, b:b+size])
        positive(target, 'actual regional polynomial convolution')
        norm = target.max()
        np.divide(target, norm, out=target)
        positive(target, 'actual regional convolution normalization')
        scale += local_scale+np.log2(norm)
        current, after = after, current
        if progress and ((group+1) % 32 == 0 or group+1 == width//4):
            print(dict(stage='actual_regional_matrix', quartets=group+1,
                       seconds=round(time.monotonic()-started, 3)), flush=True)
    inverse_choose = np.asarray([AR_LD(1)/AR_LD(comb(width, u)) for u in range(width+1)])
    np.multiply(current, inverse_choose[:, None], out=current)
    positive(current, 'actual regional old-support normalization')
    norm = current.sum(axis=1).max()
    np.divide(current, norm, out=current)
    positive(current, 'actual regional final normalization')
    return current, scale+np.log2(norm)


def _ar_region_power(matrix, scale, count):
    """Reuse Q-power arithmetic; restore the exact empty region Z once."""
    count = index(count)
    if count < 2 or count & (count-1):
        raise ValueError('Region length must be a power of two at least2')
    positive = matrix.copy()
    positive[0, 0] = 0  # Only ONE UPDATE has K00=1; internal excursions remain.
    positive, scale = _ar_dependency('actual_mds_positive_part_power')(positive, scale, count)
    logs = np.log2(positive)+scale
    logs[0, 0] = np.logaddexp2(logs[0, 0], AR_LD(0))
    scale = logs.max()
    region = np.exp2(logs-scale)
    _ar_dependency('_am_positive')(region, 'restored full actual region')
    return region, scale


def actual_regional_rank_region(width=512, region_updates=16, threshold=4096, eta=0,
                                log_z=AR_LOGS[0], log_w=AR_LOGS[1],
                                log_v=AR_LOGS[2], log_d=AR_LOGS[3], progress=False):
    """One upper regional operator; each branch uses the same population basis."""
    threshold, region_updates = index(threshold), index(region_updates)
    eta = AR_LD(eta)
    if not 0 <= threshold <= width*region_updates or not np.isfinite(eta) or eta < 0:
        raise ValueError('Require a regional threshold in range and finite eta>=0')
    plain, ps = actual_regional_rank_matrix(width, log_z, log_w, log_v, log_d, progress)
    plain, ps = _ar_region_power(plain, ps, region_updates)
    if eta == 0:
        return plain, ps
    tilted, ts = actual_regional_rank_matrix(width, AR_LD(log_z)-eta, AR_LD(log_w)-eta,
                                            log_v, log_d, progress)
    tilted, ts = _ar_region_power(tilted, ts, region_updates)
    entries = np.minimum(np.log2(plain)+ps, np.log2(tilted)+ts+eta*threshold)
    scale = entries.max()
    comparison = np.exp2(entries-scale)
    _ar_dependency('_am_positive')(comparison, 'common-population regional minimum')
    return comparison, scale


def _ar_logchoose(n, k):
    if not 0 <= k <= n:
        return float('-inf')
    if k in (0, n):
        return 0.0
    return (lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/log(2)


def actual_regional_rank_profile(profile=AR_PROFILE, volume_lower=None, death_upper=None,
                                 log_z=AR_LOGS[0], log_w=AR_LOGS[1],
                                 log_v=AR_LOGS[2], log_d=AR_LOGS[3], eta=0, progress=False):
    """One fixed(h,R,f,F) whole-threshold-class FLOAT upper diagnostic."""
    h, rank, births, footprint = map(index, profile)
    log_z, log_w, log_v, log_d = _ar_validate_logs(log_z, log_w, log_v, log_d)
    eta = AR_LD(eta)
    if log_v < 0 or log_d > 0 or not np.isfinite(eta) or eta < 0:
        raise ValueError('Threshold mode requires v>=1,0<d<=1,eta>=0')
    if not (1 <= h <= AR_M and 128*h <= rank <= 1 << 20 and 1 <= births <= footprint <= AR_N):
        raise ValueError('Invalid fixed(h,R,f,F) profile')
    G, CH, Cout = footprint-births, footprint-rank+128*h, rank-128*h
    L0 = rank+1 if volume_lower is None else index(volume_lower)
    D0 = AR_W-rank-1 if death_upper is None else index(death_upper)
    if G % 4 or G > AR_N-AR_WIDTH or L0 < 0:
        raise ValueError('Invalid canonical volume or lower threshold')
    empty = (L0 > footprint or D0 < max(0, births-AR_WIDTH)
             or not (129*h <= CH <= 256*h and 0 <= CH-births <= G
                     and 0 <= 256*h-CH <= AR_N-footprint))
    if empty:
        return dict(profile=tuple(profile), status='empty_threshold_class', full_code_certificate=False)
    route = min(0.0, _ar_logchoose(G, CH-births)+_ar_logchoose(AR_N-footprint, 256*h-CH)
                -255*_ar_logchoose(AR_M, h))
    region, scale = actual_regional_rank_region(AR_WIDTH, AR_REGION, h, eta,
                                               log_z, log_w, log_v, log_d, progress)
    # Full powers retain regional (0,0) excursions. Never reset this entry.
    moment = _ar_dependency('actual_mds_power_moment')(region, scale, AR_UPDATES//AR_REGION)
    correction = AR_LD(route)-G*log_z-births*log_w-L0*log_v-D0*log_d+eta*Cout
    return dict(profile=tuple(profile), status='actual_support_regional_capacity',
                volume_lower=L0, death_upper=D0, G=G, C_H=CH, Cout=Cout,
                logs=dict(z=float(log_z), w=float(log_w), v=float(log_v), d=float(log_d)),
                eta=float(eta), route_log2=route, external_log2_correction=float(correction),
                log2_moment=moment['log2_moment'],
                raw_log2_probability_bound=float(correction+AR_LD(moment['log2_moment'])),
                regions=256, updates_per_region=16, total_updates=4096,
                terminal='all_ones_without_reset', exact_actual_hits_retained=True,
                H_union_already_in_route=True, separate_quota_divisor=False,
                saturated_histories_only=False, field_factor_used=False,
                default_dimension_class=volume_lower is None and death_upper is None,
                full_threshold_class_covered=True, floating_diagnostic_only=True,
                full_code_certificate=False)


def _ar_partitions(width):
    if width == 4:
        return ((15,),)
    if width == 8:
        return tuple((1+sum(1 << i for i in triple), 255 ^ (1+sum(1 << i for i in triple)))
                     for triple in combinations(range(1, 8), 3))
    raise ValueError('Independent physical checks use width4 or8')


def _ar_weight(old, new, partition, markers):
    z, w, v, d = markers
    live = 0
    for quartet in partition:
        if old & quartet:
            live |= quartet
    G, births = live.bit_count(), (new & ~live).bit_count()
    return z**G*w**births*v**new.bit_count()*d**(old & ~new).bit_count(), G+births


def _ar_exact_matrix(width, markers):
    # The existing exact helper already divides its source rows by C(n,a).
    return _ar_dependency('_am_exact_power')(actual_regional_rank_local_exact(*markers), width//4)


def _ar_fraction_array(matrix):
    return np.asarray([[AR_LD(x.numerator)/x.denominator for x in row] for row in matrix], dtype=AR_LD)


def actual_regional_rank_checks():
    """Exact independent mask/partition tests, followed by FLOAT toy comparisons."""
    markers = tuple(map(Fraction, ('3/4', '2/3', '5/4', '1/2')))
    local_checks = matrix_checks = endpoint_checks = region_checks = cap_checks = 0
    table = actual_regional_rank_local_exact(*markers)
    physical = [[Fraction(0) for _ in range(5)] for _ in range(5)]
    for old in range(16):
        for new in range(16):
            physical[old.bit_count()][new.bit_count()] += _ar_weight(old, new, (15,), markers)[0]
    assert physical == table
    local_checks += 25
    exact_matrices = {}
    for width in (4, 8):
        matrix = _ar_exact_matrix(width, markers)
        exact_matrices[width] = matrix
        partitions = _ar_partitions(width)
        for a in range(width+1):
            old = (1 << a)-1
            row = [Fraction(0) for _ in range(width+1)]
            for partition in partitions:
                for new in range(1 << width):
                    row[new.bit_count()] += _ar_weight(old, new, partition, markers)[0]
            row = [entry/len(partitions) for entry in row]
            assert row == matrix[a]
            matrix_checks += width+1
        assert matrix[0][0] == 1
        assert matrix[width][width] == (markers[0]*markers[2])**width
        square = _ar_dependency('_am_matmul')(matrix, matrix)
        z, w, v, d = markers
        assert square[0][0] == (1+z**4*((1+w*v*d)**4-1))**(width//4)
        endpoint_checks += 3
    matrix4 = exact_matrices[4]
    power = [[Fraction(int(i == j)) for j in range(5)] for i in range(5)]
    for updates in (1, 2, 3):
        power = _ar_dependency('_am_matmul')(power, matrix4)
        total = Fraction(0)
        for history in product(range(16), repeat=updates):
            old, weight = 0, Fraction(1)
            for new in history:
                weight *= _ar_weight(old, new, (15,), markers)[0]
                old = new
            total += weight  # No terminal death or live-volume charge.
        assert total == sum(power[0])
        endpoint_checks += 1
    # Exact two-update region with threshold3 and penalty2^{-(F_r-3)+}.
    threshold, base = 3, Fraction(2)
    cap_physical = [[Fraction(0) for _ in range(16)] for _ in range(16)]
    for old in range(16):
        for mid in range(16):
            first, f1 = _ar_weight(old, mid, (15,), markers)
            for new in range(16):
                second, f2 = _ar_weight(mid, new, (15,), markers)
                cap_physical[old][new] += first*second*base**(-max(0, f1+f2-threshold))
    capped = [[sum(cap_physical[(1 << a)-1][b] for b in range(16) if b.bit_count() == j)
               for j in range(5)] for a in range(5)]
    for old in range(16):
        for j in range(5):
            assert sum(cap_physical[old][b] for b in range(16) if b.bit_count() == j) == capped[old.bit_count()][j]
            region_checks += 1
    plain = _ar_dependency('_am_matmul')(matrix4, matrix4)
    tilted_markers = (markers[0]/base, markers[1]/base, markers[2], markers[3])
    tilted = _ar_exact_matrix(4, tilted_markers)
    tilted = _ar_dependency('_am_matmul')(tilted, tilted)
    comparison = [[min(plain[a][b], base**threshold*tilted[a][b]) for b in range(5)] for a in range(5)]
    for a in range(5):
        for b in range(5):
            assert capped[a][b] <= comparison[a][b]
            cap_checks += 1
    # Regional population composition preserves the full physical free endpoint.
    physical_twice = _ar_dependency('_am_matmul')(cap_physical, cap_physical)
    capped_twice = _ar_dependency('_am_matmul')(capped, capped)
    comparison_twice = _ar_dependency('_am_matmul')(comparison, comparison)
    assert sum(physical_twice[0]) == sum(capped_twice[0]) <= sum(comparison_twice[0])
    endpoint_checks += 1
    # A full-H capacity consequence uses ALL footprint, including empty-old births.
    for h in range(5):
        eligible = [mask for mask in range(16) if mask.bit_count() == h]
        for footprint in range(16):
            for H in eligible:
                covered = (footprint & H).bit_count()
                assert covered <= min(h, footprint.bit_count())
                cap_checks += 1
    logs = tuple(np.log2(AR_LD(x.numerator)/x.denominator) for x in markers)
    floating_checks = 0
    for width in (4, 8):
        matrix, scale = actual_regional_rank_matrix(width, *logs)
        reference = _ar_fraction_array(exact_matrices[width])
        assert np.max(np.abs(matrix*np.exp2(scale)-reference)/reference) < AR_LD('1e-12')
        floating_checks += (width+1)**2
    region, scale = actual_regional_rank_region(4, 2, threshold, 1, *logs)
    reference = _ar_fraction_array(comparison)
    assert np.max(np.abs(region*np.exp2(scale)-reference)/reference) < AR_LD('1e-12')
    result = _ar_dependency('actual_mds_power_moment')(region, scale, 2)
    expected = sum(comparison_twice[0])
    expected_log = np.log2(AR_LD(expected.numerator)/expected.denominator)
    assert abs(AR_LD(result['log2_moment'])-expected_log) < AR_LD('1e-10')
    floating_checks += 26
    print(dict(exact_local_coefficients=local_checks, exact_physical_matrix_entries=matrix_checks,
               exact_free_endpoint_checks=endpoint_checks, exact_regional_lumpability=region_checks,
               exact_capacity_checks=cap_checks, floating_toy_checks=floating_checks,
               full_code_certificate=False), flush=True)


def actual_regional_rank_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--profile', nargs=4, type=int, default=AR_PROFILE, metavar=('H', 'R', 'F_BIRTH', 'F'))
    parser.add_argument('--volume-lower', type=int)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--log-z', type=float, default=AR_LOGS[0])
    parser.add_argument('--log-w', type=float, default=AR_LOGS[1])
    parser.add_argument('--log-v', type=float, default=AR_LOGS[2])
    parser.add_argument('--log-d', type=float, default=AR_LOGS[3])
    parser.add_argument('--eta', type=float, default=0)
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
            actual_regional_rank_checks()
        else:
            result = actual_regional_rank_profile(tuple(args.profile), args.volume_lower, args.death_upper,
                args.log_z, args.log_w, args.log_v, args.log_d, args.eta, progress=True)
            print(dict(seconds=round(time.monotonic()-started, 3), **result), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    actual_regional_rank_main()
