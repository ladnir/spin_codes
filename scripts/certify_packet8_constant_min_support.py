#!/usr/bin/env python3
"""Exact 40-packet floor for final constant-packet outer words.

The exact-length layout replaces one band-zero coordinate in each of 128
distinct data blocks by a graph coordinate.  A nonzero data block therefore
has at least 21 surviving coordinates.  A nonzero graph word has at least 22
ones, all placed in physical band-zero holes.

If a final word has ``h <= 39`` active 0xff packets, the number of active
logical vertices (data blocks, plus the graph block when nonzero) is at most
14.  Every active packet is an eight-subset of those vertices, and two such
subsets intersect in at least two vertices.  Packets in different bands would
therefore put the same pair of data blocks in a common tile in two bands,
contradicting pair capacity one.  (A graph vertex occurs only in band zero, so
an intersection with another band still contains at least two data vertices.)

All active packets must consequently lie in one band.  A data block then has
at most one original one outside that band, namely its possible band-zero
puncture.  The exact fixed-band outside distances 5,4,5 rule this out.  Thus a
nonzero final constant-packet word has at least 40 active packets.
"""

from __future__ import annotations

from certificate_spectra import load_ebch128_spectrum
from analyze_systematic_group_kernel import SPECTRUM
from certify_bch_band_projections import (
    BANDS,
    minimum_dependency,
    parity_check_columns,
)
from certify_three_band_tile_map import blocks, certify_pair_capacity
from probe_bch_forward_xor_circuit import LENGTH, target_outputs


PACKET_BITS = 8
DATA_DISTANCE = 22
PUNCTURED_DATA_DISTANCE = DATA_DISTANCE - 1
GRAPH_DISTANCE = 22
CLAIMED_PACKET_FLOOR = 40


def outside_distances() -> tuple[int, int, int]:
    outputs = target_outputs()
    result = []
    for begin, end in BANDS:
        outside = list(range(0, begin)) + list(range(end, LENGTH))
        syndromes, rank = parity_check_columns(outside, outputs)
        if rank != 64:
            raise SystemExit("constant support: outside projection lost rank")
        result.append(len(minimum_dependency(syndromes)))
    return tuple(result)


def main() -> None:
    spectrum = [(weight, count) for weight, count in load_ebch128_spectrum(SPECTRUM) if count]
    minimum_weight = next(weight for weight, count in spectrum if weight and count)
    if minimum_weight != DATA_DISTANCE:
        raise SystemExit("constant support: EBCH distance changed")
    labels = blocks()
    certify_pair_capacity(labels)
    distances = outside_distances()
    if distances != (5, 4, 5):
        raise SystemExit("constant support: fixed-band outside distances changed")

    maximum_h = CLAIMED_PACKET_FLOOR - 1
    data_vertices_graph_zero = PACKET_BITS * maximum_h // PUNCTURED_DATA_DISTANCE
    data_vertices_graph_nonzero = (
        (PACKET_BITS * maximum_h - GRAPH_DISTANCE) // PUNCTURED_DATA_DISTANCE
    )
    total_vertices_graph_nonzero = data_vertices_graph_nonzero + 1
    maximum_vertices = max(data_vertices_graph_zero, total_vertices_graph_nonzero)
    minimum_packet_intersection = 2 * PACKET_BITS - maximum_vertices
    if maximum_vertices != 14 or minimum_packet_intersection < 2:
        raise SystemExit("constant support: vertex/intersection arithmetic changed")
    # A sole-band word may lose at most one band-zero coordinate.  Every exact
    # outside distance is at least four, so no nonzero data block can survive.
    if min(distances) <= 1:
        raise SystemExit("constant support: puncture could hide an outside support")

    print("packet-8 final constant-word minimum-support certificate")
    print(f"ebch_distance={minimum_weight} punctured_data_floor={PUNCTURED_DATA_DISTANCE}")
    print(f"graph_nonzero_floor={GRAPH_DISTANCE} graph_positions=PHYSICAL_BAND_ZERO")
    print(f"fixed_band_outside_distances={','.join(map(str, distances))}")
    print("three_band_pair_capacity=1 PASS")
    print(f"largest_excluded_active_packet_count={maximum_h}")
    print(f"maximum_active_logical_vertices={maximum_vertices}")
    print(f"minimum_cross_band_packet_intersection={minimum_packet_intersection}")
    print(f"nonzero_constant_packet_support_floor={CLAIMED_PACKET_FLOOR} PASS")
    print("status=EXACT_INTEGER_CONSTANT_PACKET_MINIMUM_SUPPORT_CERTIFICATE")


if __name__ == "__main__":
    main()
