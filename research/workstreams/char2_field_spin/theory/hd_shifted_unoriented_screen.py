"""FLOAT two-variable bounds covering every orientation at fixed (r,a,b).

For singleton count s and full-pair count f, divide each oriented local
coefficient by binom(s,n10), then maximize entrywise over n10. The result
is a valid unsigned local envelope. This removes E/O from the regional
operator without sampling or independently resampling singleton labels.
All reported numerical bounds are nondirected screens.
"""

from functools import lru_cache
from itertools import product
from math import comb, log2

import numpy as np

from hd_grouped_packet_screen import LIMIT, log_multinomial
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_general_screen import general_coefficients
from hd_shifted_good_pair_screen import (
    good_pair_coefficients, good_step_coefficients, paired_frequency,
)
from hd_shifted_pairs_count import one_pairing_count, profile_count_upper


DEFAULT_PROFILE = (992, 450, 420)
DEFAULT_MARK = .49367570922234283
DEFAULT_ALPHA = .09802545440976157
DEFAULT_ODDS = ((-4.838869172554973, -2.2088270710712097),
                (-4.102509760155703, -2.296457984764048))


def unsigned_histograms(pairs, weight, even_full, odd_full):
    """Return both (empty,singleton,full) histograms."""
    return tuple((pairs-weight+full, weight-2*full, full)
                 for full in (even_full, odd_full))


@lru_cache(4)
def unsigned_coefficients(z, alpha):
    """B[:,s,f] = max_a C[:,a,s-a,f] / binom(s,a), entrywise."""
    oriented = general_coefficients(z, alpha, (0., 0., 0.))
    result = np.zeros((4, 33, 33))
    for singleton in range(33):
        for first in range(singleton+1):
            np.maximum(result[:, singleton, :],
                       oriented[:, first, singleton-first, :]/comb(singleton, first),
                       out=result[:, singleton, :])
    return result


def low_degree_frequency(coefficients, odds, grid, target):
    """Evaluate a small unsigned polynomial, folding degrees modulo grid."""
    result = np.zeros((4, grid, grid))
    for singleton in range(min(len(coefficients)-1, target[0])+1):
        for full in range(min(coefficients.shape[1]-1, target[1])+1):
            result[:, singleton % grid, full % grid] += (
                coefficients[singleton, full].ravel()
                * np.exp2(singleton*odds[0]+full*odds[1]))
    return np.fft.fftn(result, axes=(1, 2))


def positive_pair_matrix(z, alpha, odds):
    """Evaluate the unsigned positive two-step polynomial at its markers."""
    dirty = unsigned_coefficients(z, alpha).copy()
    dirty[:, :2, :] = 0
    marker = np.exp2(np.arange(33)[:, None]*odds[0]+np.arange(33)[None, :]*odds[1])
    dirty = np.sum(dirty*marker[None, :, :], axis=(1, 2))
    def evaluate(coefficients):
        weight = np.exp2(np.arange(len(coefficients))[:, None]*odds[0]
                         + np.arange(coefficients.shape[1])[None, :]*odds[1])
        return np.sum(coefficients*weight[:, :, None, None], axis=(0, 1)).ravel()
    good = evaluate(good_step_coefficients(z, alpha))
    composed = evaluate(good_pair_coefficients(z, alpha))
    return paired_frequency(dirty, good, composed).reshape((2, 2))


def regional(z, alpha, histogram, odds, *, grid=64, method='alias', cushion=1e-12):
    if len(histogram) != 3 or min(histogram) < 0 or sum(histogram) != 4096:
        raise ValueError('Expected an unsigned histogram of 4096 pairs')
    if not 0 < z < 1 or not 0 < alpha <= 1 or len(odds) != 2:
        raise ValueError('Invalid mark, alpha, or markers')
    target = tuple(histogram[1:])
    normalizer = log_multinomial(histogram)+np.dot(target, odds)
    if method == 'global':
        positive = positive_pair_matrix(z, alpha, odds)
        logs = np.log2(positive)
        return np.array(power(tuple(map(tuple, logs)), 64))-normalizer, dict(method=method)
    if method != 'alias' or grid < 64 or grid & (grid-1):
        raise ValueError('Expected alias/global and a power-of-two grid >=64')
    dirty = unsigned_coefficients(z, alpha).copy()
    dirty[:, :2, :] = 0
    dirty[:, target[0]+1:, :] = 0
    dirty[:, :, target[1]+1:] = 0
    dirty *= np.exp2(np.arange(33)[None, :, None]*odds[0]
                     + np.arange(33)[None, None, :]*odds[1])
    dirty = np.fft.fftn(dirty, s=(grid, grid), axes=(1, 2))
    good = low_degree_frequency(good_step_coefficients(z, alpha), odds, grid, target)
    composed = low_degree_frequency(good_pair_coefficients(z, alpha), odds, grid, target)
    values = paired_frequency(dirty, good, composed)
    a, b, c, d = values[:, 0, 0].real
    eigenvalue = (a+d+np.sqrt((a-d)**2+4*b*c))/2
    vector = np.array([b, eigenvalue-a] if eigenvalue-a >= eigenvalue-d else [eigenvalue-d, c])
    if not eigenvalue > 0 or not np.all(vector > 0):
        raise ArithmeticError('Degenerate positive similarity scaling')
    vector /= max(vector)
    values *= (vector[None, :]/vector[:, None]/eigenvalue).ravel()[:, None, None]
    for _ in range(6):
        a, b, c, d = values
        trace, bc = a+d, b*c
        values = np.stack((a*a+bc, b*trace, c*trace, d*d+bc))
    selected = np.fft.ifftn(values, axes=(1, 2))[:, target[0] % grid, target[1] % grid]
    if max(abs(selected.imag)) > 1e-10 or min(selected.real) < 100*max(1e-12, cushion):
        raise ArithmeticError('FLOAT alias cancellation diagnostic failed')
    logs = np.log2(selected.real+cushion).reshape((2, 2))
    logs += 64*log2(eigenvalue)+np.log2(vector[:, None]/vector[None, :])-normalizer
    return logs, dict(method=method, grid=grid,
                      minimum_normalized_alias=float(min(selected.real)),
                      maximum_imaginary_alias=float(max(abs(selected.imag))), directed=False)


def bound(profile=DEFAULT_PROFILE, z=DEFAULT_MARK, alpha=DEFAULT_ALPHA, odds=DEFAULT_ODDS,
          *, grid=64, method='alias', matching_count=False):
    r, a, b = profile
    count = (one_pairing_count(4096, r, a) if matching_count
             else profile_count_upper(4096, r, a, b))
    if count <= 0 or r <= 0:
        raise ValueError('Empty or zero-message support class')
    histograms = unsigned_histograms(4096, *profile)
    first, d0 = regional(z, alpha, histograms[0], odds[0], grid=grid, method=method)
    second, d1 = regional(z, alpha, histograms[1], odds[1], grid=grid, method=method)
    run = power(multiply(tuple(map(tuple, first)), tuple(map(tuple, second))), 128)
    count_bits = log2(count)
    prefix = count_bits+alpha*(r*(2048-128*log2(255))-log2(255)-LIMIT*log2(z))
    return dict(profile=list(profile), profile_coordinates=['r', 'a', 'b'],
                mark=z, alpha=alpha, marker_log_odds=[list(row) for row in odds],
                support_count_method='single_matching_upper' if matching_count else 'joint_triple_upper',
                log2_support_count=count_bits, log2_upper_screen=float(prefix+total(run[0])),
                method=method, grid=grid, regional_diagnostics=[d0, d1],
                all_even_odd_splits=True, directed=False, whole_code_certificate=False)


def global_box_bound(ranges, z, alpha, odds):
    """FLOAT convex bound on a rectangular all-valid (r,a,b) domain.

    The single-matching support count makes the profile exponent convex:
    -127 log multinomial(H0)-128 log multinomial(H1), plus affine terms.
    Its maximum over this box is attained at a corner. No sampled points
    stand in for the interior. Every box corner must have valid histograms;
    this helper does not replace an invalid corner by a feasible point.
    """
    if len(ranges) != 3 or any(lo > hi for lo, hi in ranges):
        raise ValueError('Three ordered profile intervals required')
    maximum, worst = -float('inf'), None
    for corner in product(*(tuple(dict.fromkeys(interval)) for interval in ranges)):
        if corner[0] <= 0 or any(min(hist) < 0 for hist in unsigned_histograms(4096, *corner)):
            raise ValueError('All box corners must have valid histograms')
        value = bound(corner, z, alpha, odds, method='global', matching_count=True)['log2_upper_screen']
        if value > maximum:
            maximum, worst = value, corner
    classes = int(np.prod([hi-lo+1 for lo, hi in ranges]))
    return dict(ranges=[list(row) for row in ranges], worst_corner=list(worst),
                class_log2_upper_screen=maximum, class_count=classes,
                union_log2_upper_screen=maximum+log2(classes),
                all_even_odd_splits=True, method='convex-global-corners', directed=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid', type=int, default=64)
    parser.add_argument('--global-bound', action='store_true')
    args = parser.parse_args()
    print(json.dumps(bound(grid=args.grid, method='global' if args.global_bound else 'alias')), flush=True)
