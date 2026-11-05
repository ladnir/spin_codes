"""Nonzero-amplitude caps for the alternating, independently swapped pair route.

Fix a regional support, before sampling its route. Let m be its size, f its
number of full pairs, and s=m-2*f its number of singleton pairs. The matching
is fixed: even/odd regions use their prescribed adjacent-pair matching.
The remaining randomness is a uniform permutation of the 4096 pairs,
independent fair swaps, and independent uniform nonzero input amplitudes.
The latter law is conditional on the XOR outer output mask, supplied by D2.

Let M_o(z) be the existing nonzero-amplitude two-state packet cap for four
RS[16,4] block occupancies o. If n(s,f,o) counts oriented arrangements in a
32-pair packet, define B[s,f]=sum_o n(s,f,o)*M_o(z)/2**s. Then

  R[m,f] = [X**s Y**f] B(X,Y)**128 / (C(4096,s)*C(4096-s,f)).

This is an upper transfer matrix, not an exact transition probability.
The coefficient extraction and route normalization themselves are exact.
The two states and refresh-density factors are unchanged from
hd_xor_multigroup; there is no fractional power or outer-density factor.
Truncating to s+2*f<=maximum preserves every requested coefficient.

Entrywise max_f R[m,f] dominates every fixed support of size m. Such maxima
may be composed in chronological regional order even when the source-pair
geometries and XOR output masks are correlated. This is a relaxation of
those fixed geometries, not independent resampling across regions/columns.
Callers must still use the actual correlated XOR mask law and source counts.
No previous scalar-route receipt, union bound, or whole-code claim transfers.
"""

from fractions import Fraction as F
from math import comb, log2

import numpy as np
from flint import arb

from hd_bucket_feedback import profiles
from hd_common_refresh_screen import principal_defects
from hd_grouped_quiet import quiet_cap
from hd_shifted_swap_sparse import (encode_arb, power_arb, power_screen,
                                    sparse_geometry, valid_degrees)
from hd_xor_multigroup import (nonzero_kernel_probability, nonzero_matrices,
                               profile_arrays)


MAXIMUM = 17
PAIR_COUNT = 4096
PACKETS = 128
ROUTE = 'alternating-adjacent-pairs-independent-fair-swaps'


def _guard(z, maximum):
    if isinstance(maximum, bool) or not isinstance(maximum, int) or not 0 <= maximum <= MAXIMUM:
        raise ValueError('Paired-route occupancy maximum must lie in 0..17')
    if not 0 < z <= 1:
        raise ValueError('Weight mark must lie in (0,1]')


def route_metadata(maximum, directed=False):
    _guard(1, maximum)
    return dict(route=ROUTE, maximum_regional_occupancy=maximum,
                field_order=256, regional_length=8192, regional_pairs=4096,
                packets_per_region=128, packet_width=64, local_rank=4,
                fixed_support_conditioning=True, full_pair_count_retained=True,
                pair_permutation_averaged=True, independent_fair_swaps=True,
                independent_nonzero_amplitudes=True, alpha=1,
                source_geometry_resampled=False, fixed_lane_assumption=False,
                scalar_route_receipts_transferred=False, directed=bool(directed),
                whole_code_certificate=False)


def _rational(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def _upper_min(*values):
    """Minimum of certified upper endpoints, not an overlapping-ball order."""
    return min(value.upper() for value in values)


def nonzero_packet_caps_arb(z, maximum):
    """Existing forced-nonzero packet caps, indexed by occupancy profile.

    This exposes the profile-wise directed formulas used by
    nonzero_arb_coefficients, whose public result has already summed over
    scalar-route supports. All entries returned here are upper endpoints.
    """
    z = F(z)
    _guard(z, maximum)
    q, h, length, state = 256, 4, 16, 16
    zero, one, zz = arb(0), arb(1), _rational(z)
    kernel = [_rational(nonzero_kernel_probability(j)) for j in range(17)]
    quiet = [zz**j*kernel[j] for j in range(17)]
    free = [(zz+(one-zz)/q)**j for j in range(17)]
    uniform_quiet = [sum((comb(j, c)*zz**(j-c)*(one-zz)**c/q**(c+min(h, j-c))
                         for c in range(j+1)), zero) for j in range(17)]
    condition = [(arb(q)/(q-1))**j for j in range(17)]
    ret = [zero]
    defects = principal_defects(h)
    for j in range(1, 17):
        rank_cap = zero
        for c in range(j+1):
            for e in range(min(h, length-j+1)):
                first = arb(q)**(-max(0, c-min(defects[c], max(0, h-j+c))))
                second = arb(q)**(max(0, j-c-h)-j)*(q**min(min(h, j), h-e)-1)
                selected = _upper_min(one-arb(q)**(-min(h, j)), first, second)
                rank_cap += comb(j, c)*comb(length-j, e)*zz**(length-c-e)*(one-zz)**(c+e)*selected
        ret.append(_upper_min(condition[j]*(free[j]-uniform_quiet[j]),
                              (one-kernel[j])*zz**max(1, length-h+1-j),
                              condition[j]*rank_cap))
    base = one/q**h
    excess = [sum((comb(n, c)*zz**(n-c)*(one-zz)**c*(one/q**c-base)
                   for c in range(min(h, n+1))), zero) for n in range(17)]
    dim = q**state-1
    refresh_density = arb(dim)/(q-1)**state
    inverse = _rational(quiet_cap(z, q=q, state=state, width=64, h=h))
    result = {}
    for profile, _ in profiles(16, 64, 4):
        degree = sum(profile)
        if degree > maximum:
            continue
        packet_ret, packet_total = zero, one
        live_nonzero, live_total, feedback_quiet, free_product = zero, one, one, one
        for occupancy in profile:
            packet_ret = packet_ret*quiet[occupancy]+packet_total*ret[occupancy]
            packet_total *= quiet[occupancy]+ret[occupancy]
            live_nonzero = live_nonzero*base+live_total*excess[length-occupancy]
            live_total *= base+excess[length-occupancy]
            feedback_quiet *= kernel[occupancy]
            free_product *= free[occupancy]
        live = free_product*live_nonzero/(one-one/q**state)*(arb(q)/(q-1))**degree
        packet_ret = _upper_min(packet_ret, dim*live)
        actual_live = refresh_density*live
        if degree == 0:
            actual_live = _upper_min(actual_live, inverse)
        entries = (zz**degree*feedback_quiet,
                   refresh_density*zz**degree*(one-feedback_quiet),
                   packet_ret/dim, actual_live)
        result[profile] = tuple(value.upper() for value in entries)
        if any(value < 0 or not value.is_finite() for value in result[profile]):
            raise ArithmeticError('Invalid nonzero-amplitude packet cap')
    return result


def pair_coefficients_screen(z, maximum):
    """FLOAT unmarked four-entry bivariate packet coefficients."""
    _guard(z, maximum)
    _, _, values = nonzero_matrices(float(z))
    occupancies, _ = profile_arrays()
    table = {tuple(occ): value.ravel() for occ, value in zip(occupancies, values)}
    result = np.zeros((4, maximum+1, maximum//2+1))
    for s, f, count, occupancy in sparse_geometry(maximum):
        result[:, s, f] += (count/2**s)*table[occupancy]
    return result


def arb_pair_coefficients(z, maximum):
    """Directed positive B[s,f]; caller controls the Arb precision."""
    table = nonzero_packet_caps_arb(z, maximum)
    result = {}
    for s, f, count, occupancy in sparse_geometry(maximum):
        values = result.setdefault((s, f), [arb(0) for _ in range(4)])
        weight = _rational(F(count, 2**s))
        for entry in range(4):
            values[entry] += weight*table[occupancy][entry]
    return {key: tuple(value.upper() for value in values) for key, values in result.items()}


def pair_regions_screen(z, maximum):
    """FLOAT log2 R[m,f], retaining the fixed full-pair histogram."""
    _guard(z, maximum)
    marker = max(1, maximum)/PAIR_COUNT
    base = pair_coefficients_screen(z, maximum)
    for s, f in valid_degrees(maximum):
        base[:, s, f] *= marker**(s+f)
    powered = power_screen(base, PACKETS, maximum)
    result = {}
    with np.errstate(divide='ignore'):
        for s, f in valid_degrees(maximum):
            normalizer = log2(comb(PAIR_COUNT, s)*comb(PAIR_COUNT-s, f))+(s+f)*log2(marker)
            result[s+2*f, f] = (np.log2(powered[:, s, f])-normalizer).reshape(2, 2)
    return result


def arb_pair_regions(z, maximum):
    """Directed R[m,f] upper endpoints at the caller's current precision."""
    z = F(z)
    _guard(z, maximum)
    marker = F(max(1, maximum), PAIR_COUNT)
    base = encode_arb(arb_pair_coefficients(z, maximum), maximum, marker)
    powered = power_arb(base, PACKETS, maximum)
    stride = 2*maximum+1
    result = {}
    for s, f in valid_degrees(maximum):
        normalizer = comb(PAIR_COUNT, s)*comb(PAIR_COUNT-s, f)*_rational(marker)**(s+f)
        result[s+2*f, f] = tuple((powered[i][s+stride*f]/normalizer).upper() for i in range(4))
    return result


def feasible_full_pairs(m, source_size=None, source_full_pairs=None):
    """Possible full-pair counts in a size-m subset of one fixed source set.

    Optional source_size=r and source_full_pairs=F preserve the information
    U subset G, |G|=r. No independent choice of source geometry is implied.
    """
    if isinstance(m, bool) or not isinstance(m, int) or not 0 <= m <= MAXIMUM:
        raise ValueError('Occupancy must lie in 0..17')
    if source_size is None and source_full_pairs is None:
        return range(m//2+1)
    if (not isinstance(source_size, int) or not isinstance(source_full_pairs, int)
            or not m <= source_size <= MAXIMUM or not 0 <= source_full_pairs <= source_size//2):
        raise ValueError('Invalid fixed source geometry')
    return range(max(0, m-source_size+source_full_pairs), min(source_full_pairs, m//2)+1)


def pair_region_operators(z, maximum, *, source_full_pairs=None):
    """FLOAT log2 matrix array (m,2,2), maximized over feasible full pairs.

    If supplied, source_full_pairs describes a fixed maximum-element source
    set under this region's matching; each regional support is its subset.
    """
    regions = pair_regions_screen(z, maximum)
    return np.array([np.max([regions[m, f] for f in feasible_full_pairs(
        m, maximum if source_full_pairs is not None else None, source_full_pairs)], axis=0)
                     for m in range(maximum+1)])


def arb_pair_region_matrices(z, maximum, *, source_full_pairs=None):
    """Directed list of flattened four-entry Arb upper matrices, indexed m."""
    regions = arb_pair_regions(z, maximum)
    return [tuple(max(regions[m, f][entry].upper() for f in feasible_full_pairs(
        m, maximum if source_full_pairs is not None else None, source_full_pairs))
                  for entry in range(4)) for m in range(maximum+1)]
