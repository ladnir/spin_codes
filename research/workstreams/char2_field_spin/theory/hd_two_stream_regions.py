"""RAW NZ two-stream regional comparison, isolated from pinned providers.

Even and odd packet positions update separate states. A region contains
128 packets (64 per stream); its parity is unchanged at the next region.
State order is00,01,10,11, with the even-packet state first. Both states
start at0 and persist across regions. These are FLOAT caps, not receipts.

For the SAME tilted raw packet polynomial B, let H=B^64. The two-stream
coefficient is [x^s y^f](H tensor H), using polynomial convolution, while
the one-stream control is [x^s y^f](H H). Both use ONE global4096-slot
multinomial denominator and ONE marker cancellation. Fair orientation
2^-s is already in each packet coefficient; halves are never normalized.
"""
from math import comb, log2

import numpy as np
from scipy.signal import convolve2d

from hd_shifted_swap_sparse import power_screen, product_screen, valid_degrees
from hd_two_sum_gf16_factorized import nonzero_matrices, profile_arrays, sparse_geometry


MAXIMUM = 64
STATE_ORDER = ('00', '01', '10', '11')
LUMPED_STATE_ORDER = ('zero', 'one-live-average', 'both-live')


def checked_maximum(maximum):
    if type(maximum) is not int or not 1 <= maximum <= MAXIMUM:
        raise ValueError('Explicit occupancy maximum1..64 required')
    return maximum


def packet_from_table(table, maximum, *, marker=1., geometry=None):
    """Build the SAME raw B used by forced_nonzero_regions, alpha1 only.

    The optional explicit geometry is for small independent tests. Rows
    are(s,f,oriented_count,occupancy); counts precede fair-swap averaging.
    No fractional power is applied to either matrices or counts.
    """
    checked_maximum(maximum)
    if not np.isfinite(marker) or not 0 < marker <= 1:
        raise ValueError('Positive finite coefficient marker at most1 required')
    rows = sparse_geometry(maximum) if geometry is None else geometry
    result = np.zeros((4, maximum+1, maximum//2+1))
    for single, full, count, occupancy in rows:
        if (any(type(x) is not int or x < 0 for x in (single, full, count))
                or single+2*full > maximum):
            raise ValueError('Nonnegative in-scope degrees and integer geometry counts required')
        values = np.asarray(table[occupancy], dtype=float)
        if values.size != 4 or np.any(~np.isfinite(values)) or np.any(values < 0):
            raise ValueError('Four finite nonnegative raw packet caps per geometry required')
        result[:, single, full] += count/2**single*values.reshape(4)*marker**(single+full)
    if np.any(~np.isfinite(result)):
        raise ArithmeticError('Raw packet coefficient overflow')
    return result


def raw_packet(z, maximum=18, *, marker=None):
    """Build raw forced-NZ packet coefficients once; no regional powers."""
    checked_maximum(maximum)
    if isinstance(z, bool) or not np.isfinite(z) or not 0 < z < 1:
        raise ValueError('Finite weight mark in(0,1) required')
    marker = maximum/4096 if marker is None else marker
    _, _, matrices = nonzero_matrices(float(z))
    occupancies, _ = profile_arrays()
    table = {tuple(occupancy): matrix.ravel()
             for occupancy, matrix in zip(occupancies, matrices)}
    return packet_from_table(table, maximum, marker=marker)


def _checked_coefficients(coefficients, maximum, entries=4):
    checked_maximum(maximum)
    values = np.asarray(coefficients, dtype=float)
    shape = entries, maximum+1, maximum//2+1
    if values.shape != shape or np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError('Complete finite nonnegative matrix coefficient array required')
    forbidden = np.arange(maximum+1)[:, None]+2*np.arange(maximum//2+1) > maximum
    if np.any(values[:, forbidden] != 0):
        raise ValueError('Coefficients outside the declared active-weight down-set must vanish')
    return values


def tensor_square_coefficients(half, maximum):
    """Polynomial H tensor H in state order00,01,10,11, truncated once.

    Ten distinct scalar convolutions suffice. Reusing conv(H_i,H_j) for
    conv(H_j,H_i) also preserves exact simultaneous stream-swap symmetry
    in the returned FLOAT array. No matrix entry products are omitted.
    """
    half = _checked_coefficients(half, maximum)
    result = np.zeros((16, maximum+1, maximum//2+1))
    forbidden = np.arange(maximum+1)[:, None]+2*np.arange(maximum//2+1) > maximum
    for left in range(4):
        for right in range(left, 4):
            coefficient = convolve2d(half[left], half[right])[:maximum+1, :maximum//2+1]
            coefficient[forbidden] = 0.
            row, column = 2*(left//2)+right//2, 2*(left%2)+right%2
            result[4*row+column] = coefficient
            row, column = 2*(right//2)+left//2, 2*(right%2)+left%2
            result[4*row+column] = coefficient
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise ArithmeticError('Invalid positive tensor coefficient')
    return result


def normalize_coefficients(coefficients, maximum, *, pairs, marker, full_max=None):
    """ONE global normalization; no independent2048-slot half averages."""
    if (type(pairs) is not int or pairs < 1 or maximum > 2*pairs
            or not np.isfinite(marker) or not 0 < marker <= 1
            or (full_max is not None and (type(full_max) is not int or not 0 <= full_max <= maximum//2))):
        raise ValueError('Valid total-pair geometry, marker, and optional full-count cutoff required')
    entries = np.shape(coefficients)[0] if np.ndim(coefficients) == 3 else 0
    if entries not in (4, 16):
        raise ValueError('Two-state or four-state coefficients required')
    values = _checked_coefficients(coefficients, maximum, entries)
    dimension = 2 if entries == 4 else 4
    result = {}
    with np.errstate(divide='ignore'):
        for single, full in valid_degrees(maximum):
            if single+full > pairs or (full_max is not None and full > full_max):
                continue
            normalizer = (log2(comb(pairs, single)*comb(pairs-single, full))
                          +(single+full)*log2(marker))
            result[single+2*full, full] = (np.log2(values[:, single, full])-normalizer).reshape(dimension, dimension)
    return result


def lump_stream_symmetry(log_matrix):
    """Exact orbit quotient in log form; one-live is the normalized average.

    The input must commute with01<->10. The quotient has three states:
    00, (01+10)/2,11. Initial00-to-terminalones products agree for every
    chronological sequence of such matrices, even when they do not commute.
    """
    matrix = np.asarray(log_matrix, dtype=float)
    if matrix.shape != (4, 4) or np.any(np.isnan(matrix)) or np.any(np.isposinf(matrix)):
        raise ValueError('Finite-or-zero4x4 log matrix required')
    swap = np.array([0, 2, 1, 3])
    if not np.array_equal(matrix, matrix[np.ix_(swap, swap)]):
        raise ValueError('Exact simultaneous stream-swap symmetry is required')
    return np.array([[matrix[0, 0], 1+matrix[0, 1], matrix[0, 3]],
                     [matrix[1, 0], np.logaddexp2(matrix[1, 1], matrix[1, 2]), matrix[1, 3]],
                     [matrix[3, 0], 1+matrix[3, 1], matrix[3, 3]]])


def controls_from_packet(packet, maximum, *, packets=128, pairs_per_packet=32,
                         marker=1., full_max=None):
    """Explicit-packet interface for tests and matched provider construction.

    `packet` already includes its marker and its singleton orientation
    average. The two streams receive alternating packet positions.
    One-stream128 control is the ordered product of the same two halves.
    """
    packet = _checked_coefficients(packet, maximum)
    if (type(packets) is not int or packets < 2 or packets % 2
            or type(pairs_per_packet) is not int or pairs_per_packet < 1):
        raise ValueError('An even positive packet count and positive slots per packet required')
    pairs = packets*pairs_per_packet
    # Validate normalization arguments before polynomial work.
    normalize_coefficients(packet, maximum, pairs=pairs, marker=marker, full_max=full_max)
    half = power_screen(packet, packets//2, maximum)
    two = tensor_square_coefficients(half, maximum)
    one = product_screen(half, half, maximum)
    one_regions = normalize_coefficients(one, maximum, pairs=pairs, marker=marker, full_max=full_max)
    two_regions = normalize_coefficients(two, maximum, pairs=pairs, marker=marker, full_max=full_max)
    return dict(one_stream=one_regions, two_stream=two_regions,
                two_stream_lumped={key: lump_stream_symmetry(matrix) for key, matrix in two_regions.items()},
                metadata=dict(total_packets=packets, packets_per_stream=packets//2,
                    pairs_per_packet=pairs_per_packet, total_pairs=pairs,
                    coefficient_maximum=maximum, returned_full_max=full_max,
                    coefficient_marker=marker, state_order=list(STATE_ORDER),
                    lumped_state_order=list(LUMPED_STATE_ORDER),
                    initial_state='00', terminal_vector=[1, 1, 1, 1],
                    streams_persist_across_regions=True, orientation_averaged_in_packet=True,
                    global_normalization_count=1, independent_half_normalization=False,
                    matched_one_stream_from_same_packet=True, regional_alpha=1,
                    directed=False, whole_code_certificate=False))


def raw_region_controls(z, maximum=18, *, full_max=None):
    """Matched one/two-stream RAW NZ regional log banks through maximum.

    Keys are(active_groups=m,full_pairs=f), singleton count=m-2f. The
    optional full_max filters outputs only; both halves retain the complete
    requested degree down-set. Existing provider defaults are untouched.
    """
    checked_maximum(maximum)
    if full_max is not None and (type(full_max) is not int or not 0 <= full_max <= maximum//2):
        raise ValueError('Optional full-pair cutoff must be in0..floor(maximum/2)')
    marker = maximum/4096
    packet = raw_packet(z, maximum, marker=marker)
    result = controls_from_packet(packet, maximum, marker=marker, full_max=full_max)
    result['metadata'].update(mark=float(z), amplitudes='independent-uniform-nonzero',
        packet_width=64, packet_state_symbols=16, stream_assignment='packet_index_mod2',
        one_stream_state_order=['0', '1'], packet_banks_built=1)
    return result


def two_stream_regions(z, maximum=18, *, full_max=None):
    """Convenience4-state view; raw_region_controls exposes matched controls."""
    return raw_region_controls(z, maximum, full_max=full_max)['two_stream']
