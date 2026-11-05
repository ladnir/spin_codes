"""Positive regional powers for the proposed coupled five-cone recurrence.

Rows and columns are Z,L,R,B,P as defined in HD_COUPLED_TWO_STREAM_CONES.
The packet caps are unpowered first moments. This module adds no local
inequality: it averages the paired geometry, forms the ordered polynomial
power B**128, and divides each requested coefficient by ONE global pair
multinomial. These FLOAT matrices are screens, not directed receipts.

Only nonnegative bidegrees occur. Keeping s+2*f<=maximum and f<=full_max
therefore preserves every requested coefficient throughout the power.
The coupled recurrence has no independent-half tensor factorization.
"""
from math import comb, log2

import numpy as np
from scipy.signal import convolve2d

from hd_shifted_swap_sparse import sparse_geometry


MAXIMUM = 64
STATE_ORDER = ('Z', 'L', 'R', 'B', 'P')


def checked_scope(maximum, full_max=None):
    if type(maximum) is not int or not 0 <= maximum <= MAXIMUM:
        raise ValueError('Explicit active-weight maximum in 0..64 required')
    full_max = maximum//2 if full_max is None else full_max
    if type(full_max) is not int or not 0 <= full_max <= maximum//2:
        raise ValueError('Full-pair cutoff must lie in 0..floor(maximum/2)')
    return full_max


def _checked_coefficients(coefficients, maximum, full_max):
    full_max = checked_scope(maximum, full_max)
    values = np.asarray(coefficients, dtype=float)
    if (values.ndim != 4 or values.shape[0] != values.shape[1]
            or not 1 <= values.shape[0] <= 16
            or values.shape[2:] != (maximum+1, full_max+1)
            or np.any(~np.isfinite(values)) or np.any(values < 0)):
        raise ValueError('Complete finite nonnegative square matrix coefficients required')
    forbidden = np.arange(maximum+1)[:, None]+2*np.arange(full_max+1) > maximum
    if np.any(values[:, :, forbidden] != 0):
        raise ValueError('Coefficients outside the declared down-set must vanish')
    return values


def product_coefficients(left, right, maximum, *, full_max=None):
    """Ordered matrix product with direct positive scalar convolutions.

    Exact-zero polynomial blocks are skipped. No symmetry, stationarity,
    or persistent sparsity of powers is assumed.
    """
    full_max = checked_scope(maximum, full_max)
    left = _checked_coefficients(left, maximum, full_max)
    right = _checked_coefficients(right, maximum, full_max)
    if left.shape != right.shape:
        raise ValueError('Matrix dimensions and declared coefficient scopes must agree')
    dimension = left.shape[0]
    result = np.zeros_like(left)
    left_nonzero = np.any(left != 0, axis=(2, 3))
    right_nonzero = np.any(right != 0, axis=(2, 3))
    for row in range(dimension):
        for middle in np.flatnonzero(left_nonzero[row]):
            for column in np.flatnonzero(right_nonzero[middle]):
                result[row, column] += convolve2d(left[row, middle], right[middle, column])[
                    :maximum+1, :full_max+1]
    forbidden = np.arange(maximum+1)[:, None]+2*np.arange(full_max+1) > maximum
    result[:, :, forbidden] = 0.
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise ArithmeticError('Invalid positive polynomial matrix product')
    return result


def power_coefficients(base, exponent, maximum, *, full_max=None):
    """Positive ordered power; exponent zero returns the polynomial identity."""
    full_max = checked_scope(maximum, full_max)
    base = _checked_coefficients(base, maximum, full_max)
    if type(exponent) is not int or exponent < 0:
        raise ValueError('Nonnegative integer packet exponent required')
    result = None
    while exponent:
        if exponent & 1:
            result = base.copy() if result is None else product_coefficients(
                result, base, maximum, full_max=full_max)
        exponent >>= 1
        if exponent:
            base = product_coefficients(base, base, maximum, full_max=full_max)
    if result is None:
        result = np.zeros_like(base)
        index = np.arange(base.shape[0])
        result[index, index, 0, 0] = 1.
    return result


def packet_from_table(table, maximum, *, marker=1., full_max=None, geometry=None):
    """Average oriented placements once, without powering caps or counts.

    Table values are raw 5x5 caps. Geometry rows (s,f,n,occupancy) count
    labelled oriented arrangements; their coefficient is n*2**(-s).
    Optional explicit rows support exact tiny independent tests.
    """
    full_max = checked_scope(maximum, full_max)
    if isinstance(marker, bool) or not np.isfinite(marker) or not 0 < marker <= 1:
        raise ValueError('Finite coefficient marker in (0,1] required')
    rows = sparse_geometry(maximum) if geometry is None else geometry
    result = np.zeros((5, 5, maximum+1, full_max+1))
    for single, full, count, occupancy in rows:
        if (any(type(x) is not int or x < 0 for x in (single, full, count))
                or single+2*full > maximum):
            raise ValueError('Nonnegative in-scope geometry counts and degrees required')
        if full > full_max:
            continue
        values = np.asarray(table[occupancy], dtype=float)
        if values.shape != (5, 5) or np.any(~np.isfinite(values)) or np.any(values < 0):
            raise ValueError('Finite nonnegative raw 5x5 packet cap required')
        result[:, :, single, full] += count/2**single*marker**(single+full)*values
    return _checked_coefficients(result, maximum, full_max)


def normalize_coefficients(coefficients, maximum, *, pairs, marker, full_max=None):
    """Return log matrices keyed by (active weight, full pairs).

    The denominator is C(pairs,s)*C(pairs-s,f)*marker**(s+f).
    Orientation averaging has already occurred inside each packet.
    """
    full_max = checked_scope(maximum, full_max)
    values = _checked_coefficients(coefficients, maximum, full_max)
    if (type(pairs) is not int or pairs < 1 or maximum > 2*pairs
            or isinstance(marker, bool) or not np.isfinite(marker) or not 0 < marker <= 1):
        raise ValueError('Valid global pair count and finite positive marker required')
    result = {}
    with np.errstate(divide='ignore'):
        for full in range(full_max+1):
            for single in range(maximum-2*full+1):
                if single+full > pairs:
                    continue
                denominator = (log2(comb(pairs, single)*comb(pairs-single, full))
                               +(single+full)*log2(marker))
                result[single+2*full, full] = np.log2(values[:, :, single, full])-denominator
    return result


def regions_from_packet(packet, maximum, *, packets=128, pairs_per_packet=32,
                        marker=1., full_max=None):
    """Explicit packet interface; all packet positions update the SAME pair state."""
    full_max = checked_scope(maximum, full_max)
    packet = _checked_coefficients(packet, maximum, full_max)
    if packet.shape[:2] != (5, 5):
        raise ValueError('The coupled provider requires the five declared cone states')
    if (type(packets) is not int or packets < 1 or type(pairs_per_packet) is not int
            or pairs_per_packet < 1):
        raise ValueError('Positive packet and pair-slot counts required')
    pairs = packets*pairs_per_packet
    # Fail before any polynomial powers if the requested normalization is invalid.
    normalize_coefficients(packet, maximum, pairs=pairs, marker=marker, full_max=full_max)
    powered = power_coefficients(packet, packets, maximum, full_max=full_max)
    regions = normalize_coefficients(powered, maximum, pairs=pairs, marker=marker,
                                    full_max=full_max)
    return dict(regions=regions, metadata=dict(state_order=list(STATE_ORDER),
        initial_state='Z', terminal_vector=[1]*5, total_packets=packets,
        pairs_per_packet=pairs_per_packet, total_pairs=pairs,
        coefficient_maximum=maximum, coefficient_full_max=full_max,
        coefficient_marker=marker, global_normalization_count=1,
        orientation_averaged_in_packet=True, source_geometry_resampled=False,
        states_persist_across_regions=True, independent_half_normalization=False,
        tensor_factorization=False, regional_alpha=1, directed=False,
        whole_code_certificate=False))


def _packet_table(z):
    # Lazy import keeps generic polynomial and toy tests independent of cap builds.
    from hd_coupled_state_packet import raw_packet_table
    return raw_packet_table(z)


def raw_region_bank(z, maximum=18, *, full_max=None):
    """Build one coupled packet table and its globally normalized regional power."""
    full_max = checked_scope(maximum, full_max)
    if maximum == 0 or isinstance(z, bool) or not np.isfinite(z) or not 0 < z < 1:
        raise ValueError('Positive occupancy maximum and finite mark in (0,1) required')
    marker = maximum/4096
    packet = packet_from_table(_packet_table(float(z)), maximum, marker=marker,
                               full_max=full_max)
    result = regions_from_packet(packet, maximum, marker=marker, full_max=full_max)
    result['metadata'].update(mark=float(z), packet_banks_built=1,
        amplitudes='independent-uniform-nonzero', packet_width=64,
        state_symbols_per_half=16, proposal='coupled-second-order-state',
        recurrence="Y=X+AU; (Uprime,Vprime)=(V,M(U+V+CX))")
    return result


def raw_regions(z, maximum=18, *, full_max=None):
    """Convenience view of raw_region_bank: (m,f) -> log 5x5 matrix."""
    return raw_region_bank(z, maximum, full_max=full_max)['regions']
