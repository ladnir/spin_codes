"""FLOAT mixed-profile screen composing low-impurity steps before alpha.

Good steps contain at most one singleton pair. Their two-step matrices are
composed before entrywise alpha. All remaining cases use the existing
single-step bound. Positive coefficients make the decomposition sound in
exact arithmetic; these FLOAT computations are not directed receipts.
"""

from functools import lru_cache
from math import comb, log2

import numpy as np

from hd_grouped_packet_screen import LIMIT, OCC, local_logs, log_multinomial
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_general_count import matching_histograms
from hd_shifted_general_screen import general_coefficients, support_bits


@lru_cache(None)
def good_geometry():
    """425 groups: 153 clean and 272 with one specified singleton type."""
    groups = []
    for a in range(17):
        for b in range(a, 17):
            groups.append((0, a+b, comb(16, a)*comb(16, b)*(1 if a == b else 2),
                           (a, a, b, b)))
    for a in range(16):
        for b in range(17):
            groups.append((1, a+b, 32*comb(15, a)*comb(16, b),
                           tuple(sorted((a, a+1, b, b)))))
    lookup = {tuple(occupancy): index for index, occupancy in enumerate(OCC)}
    impurity = np.array([row[0] for row in groups], dtype=np.int32)
    full = np.array([row[1] for row in groups], dtype=np.int32)
    counts = np.array([row[2] for row in groups], dtype=np.int64)
    indices = np.array([lookup[row[3]] for row in groups], dtype=np.int32)
    return impurity, full, counts, indices


@lru_cache(4)
def good_step_coefficients(z, alpha):
    """Unmarked G0(y)+X G1(y), where X=x10+x01, y=x11."""
    impurity, full, counts, indices = good_geometry()
    degree = impurity*33 + full
    logs = local_logs(z).reshape((-1, 4))[indices]
    weighted = np.log2(counts)[:, None] + alpha*logs
    result = np.empty((66, 4))
    for entry in range(4):
        result[:, entry] = np.bincount(degree, weights=np.exp2(weighted[:, entry]), minlength=66)
    return result.reshape((2, 33, 2, 2))


@lru_cache(4)
def good_pair_coefficients(z, alpha):
    """Unmarked ordered good-good composition before entrywise alpha.

    Axes are singleton count (0..2), full count (0..64), state, state.
    Expanding X**k retains both singleton orientations and their exact
    multiplicities. This evaluates 425**2 ordered geometry pairs.
    """
    impurity, full, counts, indices = good_geometry()
    logs = local_logs(z).reshape((-1, 2, 2))[indices]
    composed = np.logaddexp2(logs[:, None, :, 0, None] + logs[None, :, None, 0, :],
                            logs[:, None, :, 1, None] + logs[None, :, None, 1, :])
    weights = np.log2(counts)
    weighted = weights[:, None, None, None] + weights[None, :, None, None] + alpha*composed
    degree = ((impurity[:, None]+impurity[None, :])*65
              + full[:, None]+full[None, :]).ravel()
    result = np.empty((195, 2, 2))
    for row in range(2):
        for column in range(2):
            result[:, row, column] = np.bincount(
                degree, weights=np.exp2(weighted[:, :, row, column]).ravel(), minlength=195)
    return result.reshape((3, 65, 2, 2))


def square_good_coefficients(good):
    """Coefficient reference G**2, with alpha already applied to each step."""
    result = np.zeros((3, 65, 2, 2))
    for left in range(2):
        for right in range(2):
            for row in range(2):
                for column in range(2):
                    result[left+right, :, row, column] += (
                        np.convolve(good[left, :, row, 0], good[right, :, 0, column])
                        + np.convolve(good[left, :, row, 1], good[right, :, 1, column]))
    return result


def matrix_product(left, right):
    """Multiply arrays of 2x2 matrices stored as four leading entries."""
    a, b, c, d = left
    e, f, g, h = right
    return np.stack((a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h))


def paired_frequency(dirty, good, composed_good):
    """C_good2 + GD + DG + DD, retaining every matrix product's order."""
    result = composed_good.copy()
    result += matrix_product(good, dirty)
    result += matrix_product(dirty, good)
    result += matrix_product(dirty, dirty)
    return result


def low_impurity_frequency(coefficients, odds, grid, target):
    """Evaluate a low-degree X=x10+x01 polynomial on the 3D FFT grid.

    Degrees >=grid are folded, never truncated. In particular a two-step
    full degree of 64 contributes at index zero on a grid of size 64.
    Target truncation removes only monomials that cannot contribute.
    """
    if grid < 2 or grid & (grid-1):
        raise ValueError('Expected a power-of-two grid')
    roots = np.exp(-2j*np.pi*np.arange(grid)/grid)
    singleton10 = np.exp2(odds[0])*roots[:, None, None]
    singleton01 = np.exp2(odds[1])*roots[None, :, None]
    result = np.zeros((4, grid, grid, grid), dtype=np.complex128)
    for impurity, layer in enumerate(coefficients):
        orientation = np.zeros((grid, grid, 1), dtype=np.complex128)
        for first in range(impurity+1):
            second = impurity-first
            if first <= target[0] and second <= target[1]:
                orientation += comb(impurity, first)*singleton10**first*singleton01**second
        if not np.any(orientation):
            continue
        maximum = min(len(layer)-1, target[2])
        weighted = (layer[:maximum+1].reshape((-1, 4))
                    * np.exp2(np.arange(maximum+1)*odds[2])[:, None])
        folded = np.zeros((grid, 4))
        np.add.at(folded, np.arange(maximum+1) % grid, weighted)
        frequencies = np.fft.fft(folded, axis=0).T
        result += orientation[None, :, :, :]*frequencies[:, None, None, :]
    return result


def regional(z, alpha, histogram, odds, *, grid=64, compose_before_alpha=True, cushion=1e-12):
    """FLOAT regional coefficient bound from 64 improved two-step blocks."""
    if len(histogram) != 4 or min(histogram) < 0 or sum(histogram) != 4096:
        raise ValueError('Expected a nonnegative histogram of 4096 pairs')
    if not 0 < z < 1 or not 0 < alpha <= 1 or len(odds) != 3:
        raise ValueError('Invalid mark, alpha, or markers')
    if grid < 64 or grid & (grid-1):
        raise ValueError('Expected a power-of-two grid at least 64')
    target = tuple(int(n) for n in histogram[1:])
    dirty_coefficients = general_coefficients(z, alpha, odds)
    dirty_coefficients[:, 0, 0, :] = 0
    dirty_coefficients[:, 1, 0, :] = 0
    dirty_coefficients[:, 0, 1, :] = 0
    for axis, maximum in enumerate(target, 1):
        if maximum < 32:
            index = [slice(None)]*4
            index[axis] = slice(maximum+1, None)
            dirty_coefficients[tuple(index)] = 0
    dirty = np.fft.fftn(dirty_coefficients, s=(grid, grid, grid), axes=(1, 2, 3))
    del dirty_coefficients
    unmarked_good = good_step_coefficients(z, alpha)
    good = low_impurity_frequency(unmarked_good, odds, grid, target)
    unmarked_pair = (good_pair_coefficients(z, alpha) if compose_before_alpha
                     else square_good_coefficients(unmarked_good))
    composed = low_impurity_frequency(unmarked_pair, odds, grid, target)
    values = paired_frequency(dirty, good, composed)
    del dirty, good, composed
    a, b, c, d = values[:, 0, 0, 0].real
    eigenvalue = (a+d+np.sqrt((a-d)**2+4*b*c))/2
    vector = np.array([b, eigenvalue-a] if eigenvalue-a >= eigenvalue-d else [eigenvalue-d, c])
    if not eigenvalue > 0 or not np.all(vector > 0):
        raise ArithmeticError('Degenerate positive similarity scaling')
    vector /= np.max(vector)
    factors = (vector[None, :]/vector[:, None]/eigenvalue).ravel()
    values *= factors[:, None, None, None]
    for _ in range(6):
        a, b, c, d = values
        trace, bc = a+d, b*c
        values = np.stack((a*a+bc, b*trace, c*trace, d*d+bc))
    aliases = np.fft.ifftn(values, axes=(1, 2, 3))
    selected = aliases[(slice(None),)+tuple(n % grid for n in target)]
    if np.max(abs(selected.imag)) > 1e-10 or np.min(selected.real) < -1e-10:
        raise ArithmeticError('FLOAT alias cancellation diagnostic failed')
    if np.min(selected.real) < 100*max(1e-12, cushion):
        raise ArithmeticError('FLOAT alias coefficient is below its cancellation floor')
    logs = np.log2(selected.real+cushion).reshape((2, 2))
    logs += 64*log2(eigenvalue)+np.log2(vector[:, None]/vector[None, :])
    logs -= log_multinomial(histogram)+np.dot(target, odds)
    diagnostic = dict(grid=grid, composition_steps=2, good_maximum_singletons=1,
                      good_groups=425, ordered_good_geometry_pairs=425**2,
                      composed_before_alpha=compose_before_alpha,
                      minimum_normalized_alias=float(np.min(selected.real)),
                      maximum_imaginary_alias=float(np.max(abs(selected.imag))),
                      normalized_numerical_cushion=cushion, directed=False)
    return logs, diagnostic


def bound(profile, z, alpha, odds, *, grid=64, exact_count=True, compose_before_alpha=True):
    histograms = matching_histograms(4096, *profile)
    first, d0 = regional(z, alpha, histograms[0], odds[0], grid=grid,
                         compose_before_alpha=compose_before_alpha)
    second, d1 = regional(z, alpha, histograms[1], odds[1], grid=grid,
                          compose_before_alpha=compose_before_alpha)
    run = power(multiply(tuple(map(tuple, first)), tuple(map(tuple, second))), 128)
    count_bits = support_bits(tuple(profile), exact_count)
    active = sum(profile[:2])
    prefix = count_bits+alpha*(active*(2048-128*log2(255))-log2(255)-LIMIT*log2(z))
    return dict(profile=list(profile), active_groups=active, mark=z, alpha=alpha,
                marker_log_odds=[list(row) for row in odds], log2_support_count=count_bits,
                exact_support_count=exact_count, grid=grid,
                log2_upper_screen=float(prefix+total(run[0])), regional_diagnostics=[d0, d1],
                coefficient_method='good-singleton-paired-alias',
                directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid', type=int, default=64)
    parser.add_argument('--unpaired-reference', action='store_true')
    args = parser.parse_args()
    print(json.dumps(bound((512, 480, 450, 420), .49367570922234283, .09802545440976157,
        ((-5.408942950886019, -6.4540781593556, -2.2088270710712097),
         (-5.44207015077779, -4.8278390848501544, -2.296457984764048)), grid=args.grid,
        compose_before_alpha=not args.unpaired_reference)), flush=True)
