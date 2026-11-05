"""RAW eight-cone regions for XOR after refresh; isolated proposal provider.

The packet recurrence is Y=X+AU, (U',V')=(V,M(U+CX)+V). Rows/columns
are Z,L,R,B,D,P,F,G. Packet caps are unpowered first moments, not a
probability matrix. In particular F retains conditional flat density,
whereas P and G bound arbitrary probability mass on their declared sets.

The unchanged paired geometry contributes count/2**s once. Ordered B**128
coefficients receive ONE4096-slot multinomial and marker cancellation.
Positive generic products are reused from the frozen five-cone provider;
its five-state packet wrapper and its construction identity are NOT used.
"""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_shifted_swap_sparse import sparse_geometry
from hd_coupled_state_regions import (checked_scope, normalize_coefficients,
    power_coefficients, product_coefficients)


STATE_ORDER = ('Z', 'L', 'R', 'B', 'D', 'P', 'F', 'G')
CONSTRUCTION = 'proposed-g16-post-refresh-coupled-inner'
FROZEN_MARK = Fraction(691, 700)
DIRECTORY = Path(__file__).resolve().parent
SOURCE_PINS = {
    'hd_post_refresh_packet.py': 'b71bc026e32fb673478958df2fc0ec243c704024f796988d773d30c4fb607fbd',
    'hd_coupled_state_packet.py': 'faab3a80e582f0c1864b79f63b256788d99a5034f15bfbd9e781940b43bc2a43',
    'hd_coupled_state_regions.py': 'd843965b9e74d1aa86b80d09c7ab4e3a0afa594dc2c91057a7fa7497c97fa01f',
}


def packet_from_table(table, maximum, *, marker=1., full_max=None, geometry=None):
    """Positive labelled packet average; no local bounds are changed here."""
    full_max = checked_scope(maximum, full_max)
    if isinstance(marker, bool) or not np.isfinite(marker) or not 0 < marker <= 1:
        raise ValueError('Finite coefficient marker in (0,1] required')
    rows = sparse_geometry(maximum) if geometry is None else geometry
    result = np.zeros((8, 8, maximum+1, full_max+1))
    for single, full, count, occupancy in rows:
        if (any(type(x) is not int or x < 0 for x in (single, full, count))
                or single+2*full > maximum):
            raise ValueError('Nonnegative in-scope geometry degrees and counts required')
        if full > full_max:
            continue
        matrix = np.asarray(table[occupancy], dtype=float)
        if matrix.shape != (8, 8) or np.any(~np.isfinite(matrix)) or np.any(matrix < 0):
            raise ValueError('Finite nonnegative 8x8 raw packet cap required')
        result[:, :, single, full] += count/2**single*marker**(single+full)*matrix
    if np.any(~np.isfinite(result)):
        raise ArithmeticError('Raw positive packet coefficients overflowed')
    return result


def regions_from_packet(packet, maximum, *, packets=128, pairs_per_packet=32,
                        marker=1., full_max=None):
    """Explicit eight-cone packet interface for independent tiny controls."""
    full_max = checked_scope(maximum, full_max)
    if np.shape(packet) != (8, 8, maximum+1, full_max+1):
        raise ValueError('Complete eight-state packet coefficients required')
    if (type(packets) is not int or packets < 1 or type(pairs_per_packet) is not int
            or pairs_per_packet < 1):
        raise ValueError('Positive packet and pair-slot counts required')
    pairs = packets*pairs_per_packet
    # Also validates finite nonnegative coefficients and the active-weight down-set.
    normalize_coefficients(packet, maximum, pairs=pairs, marker=marker, full_max=full_max)
    powered = power_coefficients(packet, packets, maximum, full_max=full_max)
    regions = normalize_coefficients(powered, maximum, pairs=pairs, marker=marker,
                                    full_max=full_max)
    return dict(regions=regions, metadata=dict(construction=CONSTRUCTION,
        state_order=list(STATE_ORDER), initial_state='Z', terminal_vector=[1]*8,
        total_packets=packets, pairs_per_packet=pairs_per_packet, total_pairs=pairs,
        coefficient_maximum=maximum, coefficient_full_max=full_max,
        coefficient_marker=marker, global_normalization_count=1,
        orientation_averaged_in_packet=True, source_geometry_resampled=False,
        states_persist_across_regions=True, independent_half_normalization=False,
        tensor_factorization=False, regional_alpha=1,
        cone_P_arbitrary_diagonal_probability=True,
        cone_G_arbitrary_first_nonzero_probability=True,
        cone_F_conditional_second_flat_density=True,
        probability_rows_assumed=False, directed=False, whole_code_certificate=False))


def _packet_table(z):
    # No production cap table is built by importing this module or its tests.
    from hd_post_refresh_packet import raw_packet_table
    return raw_packet_table(z)


def raw_region_bank(z, maximum=18, *, full_max=None):
    """Return all requested (m,f) log8x8 caps and explicit proposal metadata."""
    full_max = checked_scope(maximum, full_max)
    if maximum == 0 or isinstance(z, bool) or not np.isfinite(z) or not 0 < z < 1:
        raise ValueError('Positive occupancy maximum and finite mark in (0,1) required')
    marker = maximum/4096
    packet = packet_from_table(_packet_table(float(z)), maximum, marker=marker,
                               full_max=full_max)
    result = regions_from_packet(packet, maximum, marker=marker, full_max=full_max)
    result['metadata'].update(mark=float(z), packet_banks_built=1,
        amplitudes='independent-uniform-nonzero', packet_width=64,
        state_symbols_per_half=16, recurrence="Y=X+AU; (Uprime,Vprime)=(V,M(U+CX)+V)")
    return result


def raw_regions(z, maximum=18, *, full_max=None):
    return raw_region_bank(z, maximum, full_max=full_max)['regions']


def serialize_bank(bank, *, maximum=18, full_max=1):
    """Complete frozen-scope bank; null encodes exactly zero, never a floor."""
    full_max = checked_scope(maximum, full_max)
    expected = {(m, f) for m in range(maximum+1)
                for f in range(min(full_max, m//2)+1)}
    if set(bank) != expected:
        raise ValueError('Complete and exact requested histogram scope required')
    records = []
    for m, f in sorted(expected):
        matrix = np.asarray(bank[m, f], dtype=float)
        if matrix.shape != (8, 8) or np.any(np.isnan(matrix)) or np.any(np.isposinf(matrix)):
            raise ValueError('Finite-or-zero eight-state log matrix required')
        records.append(dict(active_weight=m, full_pairs=f,
            log2_matrix=[[None if x == -np.inf else float(x) for x in row] for row in matrix]))
    return records


def build_saved_bank(*, max_seconds=30):
    """One fixed z=691/700, maximum18, F1 bank; no outer bound or fitting.

    The elapsed-time guard is checked after the indivisible positive
    packet/power build. An over-budget bank is retained but labelled
    incomplete, so a saved-bank gate can fail closed rather than rebuild it.
    """
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive resource budget at most thirty seconds required')
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed dependency changed: '+name)
    started = perf_counter()
    provided = raw_region_bank(float(FROZEN_MARK), 18, full_max=1)
    metadata = provided['metadata']
    required = dict(construction=CONSTRUCTION, state_order=list(STATE_ORDER),
        global_normalization_count=1, regional_alpha=1, directed=False,
        whole_code_certificate=False)
    if any(metadata.get(key) != value for key, value in required.items()):
        raise ValueError('Provider metadata does not identify the reviewed post-refresh law')
    records = serialize_bank(provided['regions'])
    elapsed = perf_counter()-started
    result = dict(construction=CONSTRUCTION, mark=str(FROZEN_MARK),
        maximum_regional_occupancy=18, full_pair_maximum=1, regional_alpha=1,
        state_order=list(STATE_ORDER), initial_state='Z', terminal_vector=[1]*8,
        scope='raw regional caps through m18 and f<=1 only; no outer-class bound',
        regional_histogram_count=len(records), regional_log2_matrices=records,
        provider_metadata=metadata, zero_matrix_entries_encoded_as_none=True,
        raw_packet_banks_built=1, parameter_fits=0, reviewed_source_sha256=SOURCE_PINS,
        compute_budget_seconds=max_seconds, elapsed_seconds=elapsed,
        complete_requested_bank=elapsed < max_seconds, directed=False,
        whole_code_certificate=False, source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    if not result['complete_requested_bank']:
        result['incomplete_reason'] = 'The one indivisible raw bank build exceeded its resource budget'
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=30)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen raw bank')
    answer = build_saved_bank(max_seconds=args.max_seconds)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=answer['complete_requested_bank'],
        histograms=answer['regional_histogram_count'], elapsed_seconds=answer['elapsed_seconds'],
        source_sha256=answer['source_sha256']), allow_nan=False))
