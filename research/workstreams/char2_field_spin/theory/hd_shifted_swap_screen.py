"""FLOAT unsigned regional bounds for independently swapped routed pairs.

Swaps are independent fair setup bits for every pair in every region.
The one-step coefficients average already-fractional packet caps. The
orientation-invariant good-good polynomial is unchanged. No numerical
result in this module is a directed certificate or a whole-code claim.
"""
from functools import lru_cache
from math import log2

import numpy as np

from hd_grouped_packet_screen import LIMIT, local_logs, log_multinomial
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_bulk_search import LOG_DENSITY, PositiveWitness
from hd_shifted_good_pair_screen import (good_pair_coefficients,
    good_step_coefficients, paired_frequency)
from hd_shifted_pairs_count import profile_count_upper
from hd_shifted_pairs_screen import half_geometry, three_type_geometry
from hd_shifted_unoriented_screen import low_degree_frequency, unsigned_histograms


CONSTRUCTION = 'toeplitz-gf256-h4-alternating-swapped-pairs-k20-rate-half'


@lru_cache(4)
def swap_coefficients(z, alpha):
    """Return B[:,s,f] = 2**(-s) sum_u C[:,u,s-u,f], unmarked."""
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError('Invalid mark or fractional exponent')
    counts, weights, degrees, indices = three_type_geometry()
    halves = weights-counts[:, 0]-counts[:, 1]
    weights = (halves[:, None]+halves[None, :]).ravel()
    degrees = (degrees//33**2+(degrees//33) % 33)*33+degrees % 33
    logs = local_logs(z).reshape((-1, 4))
    result = np.empty((4, 33**2))
    for entry in range(4):
        result[entry] = np.bincount(degrees, weights=np.exp2(weights+alpha*logs[indices, entry]),
                                    minlength=33**2)
    return result.reshape((4, 33, 33))


@lru_cache(1)
def _half_terms():
    halves, terms, _ = half_geometry()
    rows = [(3*i+min(2, s), s, f, 2**(n-s)*(1 if a == b else 2))
            for i, ((a, b), row) in enumerate(zip(halves, terms)) for s, f, n in row]
    return tuple(np.array([row[j] for row in rows], dtype=int if j < 3 else float)
                 for j in range(4))


@lru_cache(4)
def _positive_caps(z, alpha):
    return np.exp2(alpha*local_logs(z).reshape((-1, 4))[half_geometry()[2]])


def _evaluate_small(coefficients, odds):
    marker = np.exp2(np.arange(len(coefficients))[:, None]*odds[0]
                     +np.arange(coefficients.shape[1])[None, :]*odds[1])
    return np.sum(coefficients*marker[:, :, None, None], axis=(0, 1)).ravel()


def positive_pair_matrix(z, alpha, odds):
    """Evaluate the two-step operator without rebuilding all coefficients.

    The 153 sorted half-occupancies are split into zero, one, or at least
    two singleton pairs. This computes dirty contributions directly and
    avoids subtracting nearly equal positive matrices on sparse profiles.
    """
    if not 0 < z < 1 or not 0 < alpha <= 1 or len(odds) != 2:
        raise ValueError('Invalid mark, alpha, or marker pair')
    indices, singleton, full, counts = _half_terms()
    weights = counts*np.exp2(singleton*odds[0]+full*odds[1])
    halves = np.bincount(indices, weights=weights, minlength=153*3).reshape((153, 3))
    clean, one, dirty = halves.T
    combined = (dirty[:, None]*halves.sum(axis=1)[None, :]
                +(clean+one)[:, None]*dirty[None, :]+one[:, None]*one[None, :])
    dirty_matrix = combined.ravel() @ _positive_caps(z, alpha)
    good = _evaluate_small(good_step_coefficients(z, alpha), odds)
    composed = _evaluate_small(good_pair_coefficients(z, alpha), odds)
    return paired_frequency(dirty_matrix, good, composed).reshape((2, 2))


class SwapPositiveWitness(PositiveWitness):
    @classmethod
    def create(cls, z, alpha, odds, centre=(0, 0, 0)):
        odds = tuple(tuple(map(float, row)) for row in odds)
        matrices, regions = {}, []
        for row in odds:
            if row not in matrices:
                positive = positive_pair_matrix(z, alpha, row)
                if not np.all(positive > 0) or not np.all(np.isfinite(positive)):
                    raise ArithmeticError('Invalid positive matrix evaluation')
                matrices[row] = power(tuple(map(tuple, np.log2(positive))), 64)
            regions.append(matrices[row])
        numerator = float(total(power(multiply(regions[0], regions[1]), 128)[0]))
        intercept = numerator-alpha*(log2(255)+LIMIT*log2(z))
        (s0, f0), (s1, f1) = odds
        slopes = (alpha*LOG_DENSITY+1-128*(s0+s1),
                  -2+256*s0-128*f0, 256*s1-128*f1)
        return cls(float(z), float(alpha), odds, tuple(centre), intercept, slopes, numerator)

    def serializable(self):
        return dict(super().serializable(), construction=CONSTRUCTION, pair_swaps=True,
                    orientation_aggregation='independent-fair-swaps')


def regional(z, alpha, histogram, odds, *, grid=64, method='alias', cushion=1e-12):
    if len(histogram) != 3 or min(histogram) < 0 or sum(histogram) != 4096:
        raise ValueError('Expected an unsigned histogram of 4096 pairs')
    if not 0 < z < 1 or not 0 < alpha <= 1 or len(odds) != 2:
        raise ValueError('Invalid mark, alpha, or markers')
    target = tuple(map(int, histogram[1:]))
    normalizer = log_multinomial(histogram)+np.dot(target, odds)
    if method == 'global':
        positive = positive_pair_matrix(z, alpha, odds)
        if not np.all(np.isfinite(positive)) or not np.all(positive > 0):
            raise ArithmeticError('Invalid positive matrix evaluation')
        return np.array(power(tuple(map(tuple, np.log2(positive))), 64))-normalizer, dict(method=method)
    if method != 'alias' or grid < 64 or grid & (grid-1):
        raise ValueError('Expected alias/global and a power-of-two grid >=64')
    dirty = swap_coefficients(z, alpha).copy()
    dirty[:, :2, :] = 0
    dirty[:, target[0]+1:, :] = 0
    dirty[:, :, target[1]+1:] = 0
    dirty *= np.exp2(np.arange(33)[None, :, None]*odds[0]
                     +np.arange(33)[None, None, :]*odds[1])
    dirty = np.fft.fftn(dirty, s=(grid, grid), axes=(1, 2))
    good = low_degree_frequency(good_step_coefficients(z, alpha), odds, grid, target)
    composed = low_degree_frequency(good_pair_coefficients(z, alpha), odds, grid, target)
    values = paired_frequency(dirty, good, composed)
    peak = float(np.max(values[:, 0, 0].real))
    if not 0 < peak < float('inf') or not np.all(np.isfinite(values)):
        raise ArithmeticError('Nonfinite unscaled alias matrix')
    # Large finite boundary markers can overflow the quadratic eigenvalue
    # formula even though the matrix itself is representable.
    values /= peak
    a, b, c, d = values[:, 0, 0].real
    eigenvalue = (a+d+np.sqrt((a-d)**2+4*b*c))/2
    vector = np.array([b, eigenvalue-a] if eigenvalue-a >= eigenvalue-d else [eigenvalue-d, c])
    if not 0 < eigenvalue < float('inf') or not np.all(vector > 0) or not np.all(np.isfinite(vector)):
        raise ArithmeticError('Degenerate positive similarity scaling')
    vector /= max(vector)
    values *= (vector[None, :]/vector[:, None]/eigenvalue).ravel()[:, None, None]
    for _ in range(6):
        a, b, c, d = values
        trace, bc = a+d, b*c
        values = np.stack((a*a+bc, b*trace, c*trace, d*d+bc))
    selected = np.fft.ifftn(values, axes=(1, 2))[:, target[0] % grid, target[1] % grid]
    if (not np.all(np.isfinite(selected)) or max(abs(selected.imag)) > 1e-10
            or min(selected.real) < 100*max(1e-12, cushion)):
        raise ArithmeticError('FLOAT alias cancellation diagnostic failed')
    logs = np.log2(selected.real+cushion).reshape((2, 2))
    logs += 64*(log2(peak)+log2(eigenvalue))+np.log2(vector[:, None]/vector[None, :])-normalizer
    return logs, dict(method=method, grid=grid, minimum_normalized_alias=float(min(selected.real)),
                      maximum_imaginary_alias=float(max(abs(selected.imag))), directed=False)


@lru_cache(256)
def support_bits(profile):
    count = profile_count_upper(4096, *profile)
    if count <= 0 or profile[0] <= 0:
        raise ValueError('Empty or zero-message support class')
    return log2(count)


def bound(profile, z, alpha, odds, *, grid=64, method='global'):
    profile = tuple(map(int, profile))
    count_bits = support_bits(profile)
    histograms = unsigned_histograms(4096, *profile)
    first, d0 = regional(z, alpha, histograms[0], odds[0], grid=grid, method=method)
    second, d1 = regional(z, alpha, histograms[1], odds[1], grid=grid, method=method)
    run = power(multiply(tuple(map(tuple, first)), tuple(map(tuple, second))), 128)
    prefix = count_bits+alpha*(profile[0]*LOG_DENSITY-log2(255)-LIMIT*log2(z))
    return dict(construction=CONSTRUCTION, pair_swaps=True, orientation_aggregation='independent-fair-swaps',
                profile=list(profile), mark=z, alpha=alpha, marker_log_odds=[list(row) for row in odds],
                log2_support_count=count_bits, support_count_method='joint_triple_upper',
                log2_upper_screen=float(prefix+total(run[0])), method=method, grid=grid,
                regional_diagnostics=[d0, d1], all_even_odd_splits=True,
                directed=False, whole_code_certificate=False)
