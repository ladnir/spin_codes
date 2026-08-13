#!/usr/bin/env python3
"""Exact reduced-size rank probe for the {0x00,0xff} packet subspace.

Use the first two fixed bands, of sizes 42 and 43.  A constant packet has one
binary variable shared by its eight lanes.  For ``T`` tiles there are therefore

    8 * T * (42+43) = 680 T

packet variables.  The combined 85-coordinate EBCH projection has rank 64 and
21 parity checks.  Substitute each block's selected packet variables into
those checks after independently sampling its within-band coordinate
permutations and every physical group's balanced random lane-to-packet map.

The resulting binary matrix is reduced exactly with integer XOR elimination.
The all-ones EBCH word always supplies one solution, so nullity one is the
strongest possible outcome.  Runs with ``T<256`` are reduced-size diagnostics,
not a proof for the full construction.
"""

from __future__ import annotations

import argparse
import random

from certify_bch_band_projections import BANDS, parity_check_columns
from probe_bch_forward_xor_circuit import target_outputs


LANES = 64
PACKET_LANES = 8
BAND_SIZES = (42, 43)
CHECKS = sum(BAND_SIZES) - 64


def pair_syndromes() -> list[int]:
    coordinates = list(range(BANDS[0][0], BANDS[0][1])) + list(
        range(BANDS[1][0], BANDS[1][1])
    )
    syndromes, rank = parity_check_columns(coordinates, target_outputs())
    if rank != 64 or len(syndromes) != sum(BAND_SIZES):
        raise SystemExit("constant packet rank: combined projection changed")
    if max(syndromes).bit_length() != CHECKS:
        raise SystemExit("constant packet rank: parity-check dimension changed")
    return syndromes


def variable_index(band: int, tile: int, column: int, packet: int, tiles: int) -> int:
    band_offset = 0 if band == 0 else tiles * BAND_SIZES[0] * 8
    return band_offset + (tile * BAND_SIZES[band] + column) * 8 + packet


def exact_rank(rows, columns: int) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        value = row
        while value:
            pivot = value.bit_length() - 1
            prior = basis.get(pivot)
            if prior is None:
                basis[pivot] = value
                break
            value ^= prior
        if len(basis) == columns:
            return columns
    return len(basis)


def build_rows(tiles: int, seed: int):
    rng = random.Random(seed)
    syndromes = pair_syndromes()

    # packet_maps[band][tile][column][lane] is the packet containing that lane.
    packet_maps = []
    for band, columns in enumerate(BAND_SIZES):
        band_maps = []
        for tile in range(tiles):
            tile_maps = []
            for column in range(columns):
                order = list(range(LANES))
                rng.shuffle(order)
                lane_to_packet = [0] * LANES
                for position, lane in enumerate(order):
                    lane_to_packet[lane] = position // PACKET_LANES
                tile_maps.append(lane_to_packet)
            band_maps.append(tile_maps)
        packet_maps.append(band_maps)

    for tile0 in range(tiles):
        for lane in range(LANES):
            tile1 = (tile0 + 9 * lane) % tiles
            selected_variables = []
            for band, (tile, columns) in enumerate(
                ((tile0, BAND_SIZES[0]), (tile1, BAND_SIZES[1]))
            ):
                permutation = list(range(columns))
                rng.shuffle(permutation)
                for coordinate in range(columns):
                    physical_column = permutation[coordinate]
                    packet = packet_maps[band][tile][physical_column][lane]
                    selected_variables.append(
                        variable_index(
                            band, tile, physical_column, packet, tiles
                        )
                    )
            for check in range(CHECKS):
                row = 0
                for coordinate, syndrome in enumerate(syndromes):
                    if syndrome >> check & 1:
                        row ^= 1 << selected_variables[coordinate]
                if row:
                    yield row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tiles", default="1,2,4,8,16")
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--seed-base", type=int, default=20260810)
    args = parser.parse_args()
    tile_counts = [int(value) for value in args.tiles.split(",")]
    if any(value <= 0 or value > 256 or value & (value - 1) for value in tile_counts):
        raise SystemExit("constant packet rank: tiles must be powers of two through 256")
    if args.seeds <= 0:
        raise SystemExit("constant packet rank: seeds must be positive")

    print("constant-packet two-band exact rank probe")
    print("combined_projection=[85,64] parity_checks=21")
    for tiles in tile_counts:
        columns = tiles * sum(BAND_SIZES) * 8
        for trial in range(args.seeds):
            seed = args.seed_base + 1000 * tiles + trial
            rank = exact_rank(build_rows(tiles, seed), columns)
            nullity = columns - rank
            if nullity < 1:
                raise SystemExit("constant packet rank: lost all-ones solution")
            print(
                f"tiles={tiles} seed={seed} variables={columns} "
                f"rank={rank} nullity={nullity}"
            )
    print("status=EXACT_REDUCED_SIZE_RANDOM_INSTANCE_RANK_DIAGNOSTIC")


if __name__ == "__main__":
    main()
