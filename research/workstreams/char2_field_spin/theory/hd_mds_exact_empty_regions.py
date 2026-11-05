"""FLOAT regional adapter changing only the empty-packet11 cap.

The live reference remains uniform on the nonzero16-symbol states.
The exact inverse-refresh weight law combines with the actual emission
enumerator to bound output density. All nonempty packet entries, route
counts, pair swaps, chronology, and regional normalizers stay unchanged.
This is an isolated experiment; frozen providers and receipts are not edited.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import comb, log2

import numpy as np

from hd_grouped_quiet import weight_transition
from hd_mds_empty_weight import emission_by_weight
from hd_shifted_swap_sparse import sparse_geometry, power_screen, valid_degrees
from hd_xor_multigroup import nonzero_matrices, profile_arrays


def exact_empty_cap(mark, *, check_time=lambda: None):
    if not isinstance(mark, (int, F)) or not 0 < mark <= 1:
        raise ValueError('Exact rational weight mark required')
    means = emission_by_weight(F(mark), check_time=check_time)
    transition = weight_transition()
    # Inverse refresh has the same weight transition, since H^-1 is MDS.
    candidates = tuple(sum((p*mean for p, mean in zip(row, means[1:])), F(0))
                       for row in transition)
    return max(candidates), candidates, means


def forced_nonzero_regions(mark, maximum=18, *, check_time=lambda: None):
    mark = F(mark)
    if type(maximum) is not int or not 1 <= maximum <= 64 or not 0 < mark < 1:
        raise ValueError('Exact mark in(0,1) and maximum1..64 required')
    cap, candidates, means = exact_empty_cap(mark, check_time=check_time)
    check_time()
    degrees, _, matrices = nonzero_matrices(float(mark))
    occupancy, _ = profile_arrays()
    empty = np.flatnonzero(degrees == 0)
    if len(empty) != 1 or tuple(occupancy[empty[0]]) != (0, 0, 0, 0):
        raise AssertionError('Exactly one empty packet profile required')
    before = matrices.copy()
    old = float(matrices[empty[0], 1, 1])
    matrices[empty[0], 1, 1] = min(old, float(cap))
    changed = np.argwhere(matrices != before)
    if changed.tolist() not in ([], [[int(empty[0]), 1, 1]]):
        raise AssertionError('Only empty-packet11 may change')
    table = {tuple(occ): matrix.ravel() for occ, matrix in zip(occupancy, matrices)}
    base = np.zeros((4, maximum+1, maximum//2+1))
    marker = maximum/4096
    for single, full, count, occ in sparse_geometry(maximum):
        base[:, single, full] += count/2**single*table[occ]*marker**(single+full)
    check_time()
    powered = power_screen(base, 128, maximum)
    if np.any(~np.isfinite(powered)) or np.any(powered < 0):
        raise ArithmeticError('Invalid FLOAT positive regional coefficients')
    result = {}
    with np.errstate(divide='ignore'):
        for single, full in valid_degrees(maximum):
            normalizer = (log2(comb(4096, single)*comb(4096-single, full))
                          +(single+full)*log2(marker))
            result[single+2*full, full] = (np.log2(powered[:, single, full])-normalizer).reshape(2, 2)
    return result, dict(old_empty11_log2=log2(old), new_empty11_log2=log2(cap),
        changed_packet_entries=changed.tolist(), exact_cap=str(cap),
        exact_cap_sha256=sha256(str(cap).encode()).hexdigest(),
        worst_output_weight=1+candidates.index(cap),
        reference_cone_unchanged=True, nonempty_packet_caps_unchanged=True,
        packet_profile_count=len(matrices), mark=str(mark), maximum=maximum,
        directed=False, whole_code_certificate=False)
