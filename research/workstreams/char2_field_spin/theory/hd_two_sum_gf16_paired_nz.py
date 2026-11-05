"""FLOAT NZ regional provider composing good-good packets before alpha.

Good packets have at most one singleton source pair. Their cap is invariant
under that singleton's fair orientation. The existing good_geometry counts
therefore require no extra2**(-s). For each ORDERED pair of good geometries,
compose the two raw NZ cap matrices, then power each entry. The other cases
retain G*D+D*G+D*D with the single-packet alpha caps.

The64 positive two-packet blocks preserve all128 packet positions and the
same unpowered route normalization. This new provider neither edits nor
relabels the pinned packetwise-alpha provider. Arithmetic remains FLOAT.
"""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from hd_grouped_packet_screen import OCC
from hd_shifted_good_pair_screen import good_geometry
from hd_shifted_swap_sparse import power_screen, product_screen, sparse_geometry
from hd_two_sum_gf16_fractional_state import (checked_alpha,
    fractional_packet_coefficients, normalize_regional_coefficients)
from hd_xor_multigroup import nonzero_matrices, profile_arrays


MAXIMUM = 64


def checked_maximum(maximum):
    if type(maximum) is not int or not 1 <= maximum <= MAXIMUM:
        raise ValueError('Target occupancy maximum1..64 required')
    return maximum


@lru_cache(8)
def relevant_good_geometry(maximum):
    """Unoriented good counts:425 full rows before target down-set truncation."""
    checked_maximum(maximum)
    singles, full, counts, indices = good_geometry()
    return tuple((int(s), int(f), int(n), tuple(map(int, OCC[index])))
                 for s, f, n, index in zip(singles, full, counts, indices)
                 if s+2*f <= maximum)


def checked_table_entries(table, occupancies):
    values = np.asarray([table[occupancy] for occupancy in occupancies], dtype=float)
    if values.shape != (len(occupancies), 4) or np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError('The explicit packet table must contain four finite nonnegative caps per profile')
    return values.reshape(-1, 2, 2)


def good_composition_from_table(table, alpha, maximum, marker):
    """C=sum_g,h count_g count_h (M_g M_h)**alpha, before regional power.

    Geometry multiplicities and coefficient markers are outside alpha.
    Products are kept ordered; exchanging packet positions is not assumed.
    """
    checked_maximum(maximum)
    alpha = checked_alpha(alpha)
    if not 0 < marker <= 1:
        raise ValueError('Positive coefficient marker at most1 required')
    rows = relevant_good_geometry(maximum)
    singles = np.array([row[0] for row in rows], dtype=np.int32)
    full = np.array([row[1] for row in rows], dtype=np.int32)
    counts = np.array([row[2] for row in rows], dtype=float)
    matrices = checked_table_entries(table, [row[3] for row in rows])
    composed = np.matmul(matrices[:, None], matrices[None, :])
    composed **= float(alpha)
    s = singles[:, None]+singles[None, :]
    f = full[:, None]+full[None, :]
    keep = s+2*f <= maximum
    weights = counts[:, None]*counts[None, :]*marker**(s+f)
    nf = maximum//2+1
    degrees = (s*nf+f)[keep]
    result = np.empty((4, maximum+1, nf))
    for entry in range(4):
        values = (weights*composed[:, :, entry//2, entry%2])[keep]
        result[entry] = np.bincount(degrees, weights=values,
                                   minlength=(maximum+1)*nf).reshape(maximum+1, nf)
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise ArithmeticError('Invalid positive good-good coefficient array')
    return result


def paired_components_from_table(table, alpha, maximum, *, marker=None):
    """Public explicit-table interface for independent coefficient tests."""
    checked_maximum(maximum)
    alpha = checked_alpha(alpha)
    marker = maximum/4096 if marker is None else marker
    single = fractional_packet_coefficients(table, sparse_geometry(maximum), alpha, maximum, marker)
    good, dirty = single.copy(), single.copy()
    good[:, 2:, :] = 0.
    dirty[:, :2, :] = 0.
    composed = good_composition_from_table(table, alpha, maximum, marker)
    base = (composed+product_screen(good, dirty, maximum)
            +product_screen(dirty, good, maximum)+product_screen(dirty, dirty, maximum))
    if np.any(~np.isfinite(base)) or np.any(base < 0):
        raise ArithmeticError('Invalid positive two-packet polynomial')
    return dict(single=single, good=good, dirty=dirty, composed=composed, paired=base)


def regions_from_table(table, alpha, maximum, *, marker=None):
    """All regional histograms from an explicit raw NZ cap table."""
    checked_maximum(maximum)
    marker = maximum/4096 if marker is None else marker
    pieces = paired_components_from_table(table, alpha, maximum, marker=marker)
    run = power_screen(pieces['paired'], 64, maximum)
    return normalize_regional_coefficients(run, maximum, marker)


def paired_nz_regions(z, alpha, maximum=32):
    """New paired-composition NZ provider, compatible with fractional_regions."""
    checked_maximum(maximum)
    checked_alpha(alpha)
    if not 0 < z < 1:
        raise ValueError('Weight mark in(0,1) required')
    _, _, matrices = nonzero_matrices(float(z))
    occupancies, _ = profile_arrays()
    table = {tuple(occ): matrix.ravel() for occ, matrix in zip(occupancies, matrices)}
    return regions_from_table(table, F(alpha), maximum)


def provider_metadata():
    return dict(amplitudes='independent-uniform-nonzero', local_packet_width=64,
        packets_per_region=128, composition_steps=2, paired_blocks=64,
        good_maximum_singletons=1, ordered_good_composition=True,
        good_pair_composed_before_alpha=True, other_cases='GD+DG+DD',
        extra_good_orientation_divisor=False, route_denominator_powered=False,
        regional_pairs=4096, pair_swaps=True, fixed_route_fractional_chronology=True,
        directed=False, proposal_requires_fresh_directed_replay=True, whole_code_certificate=False)
