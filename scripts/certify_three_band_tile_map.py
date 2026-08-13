#!/usr/bin/env python3
"""Exact certificate for the fixed three-band tile map.

There are 256 base tiles and 64 lanes.  Block ``(t,l)`` is sent to

    (t, t + 9 l, t + 20 l) mod 256.

The script exhaustively certifies:

* every band tile contains exactly 64 blocks;
* every pair of band labels identifies at most one block;
* there is no four-block Pasch/intercalate (two disjoint collisions in
  every band); and
* the exact number of K4-minus-one-edge motifs, split according to the
  band containing the single collision.

The four-block enumeration is over alternating rectangles in each pair of
bands.  Such a rectangle supplies two disjoint collision pairs in each of
those bands.  Its two diagonals are the only remaining block pairs.  Zero,
one, or two diagonal collisions in the third band respectively give an
ordinary rectangle, K4-minus-one-edge, or Pasch motif.  Every unordered
four-block set is enumerated exactly once for its two dense bands.
"""

from __future__ import annotations

from itertools import combinations


TILES = 256
LANES = 64
SLOPES = (0, 9, 20)
EXPECTED_K4_MINUS_BY_SPARSE_BAND = (79616, 122880, 79360)


def tile_labels(base: int, lane: int) -> tuple[int, int, int]:
    return tuple((base + slope * lane) % TILES for slope in SLOPES)


def blocks() -> list[tuple[int, int, int]]:
    return [
        tile_labels(base, lane)
        for base in range(TILES)
        for lane in range(LANES)
    ]


def certify_balance(labels: list[tuple[int, int, int]]) -> None:
    for band in range(3):
        counts = [0] * TILES
        for block in labels:
            counts[block[band]] += 1
        if any(count != LANES for count in counts):
            raise AssertionError(f"band {band} is not exactly balanced")


def certify_pair_capacity(labels: list[tuple[int, int, int]]) -> None:
    for first, second in combinations(range(3), 2):
        occupied: set[tuple[int, int]] = set()
        for block in labels:
            pair = (block[first], block[second])
            if pair in occupied:
                raise AssertionError(
                    f"bands {(first, second)} have pair capacity greater than one"
                )
            occupied.add(pair)


def projection_matrix(
    labels: list[tuple[int, int, int]], row_band: int, column_band: int
) -> tuple[list[int], list[int]]:
    """Return row-neighbor bitsets and a dense label-pair-to-block table."""

    row_bits = [0] * TILES
    table = [-1] * (TILES * TILES)
    for block_index, block in enumerate(labels):
        row = block[row_band]
        column = block[column_band]
        offset = row * TILES + column
        if table[offset] != -1:
            raise AssertionError("projection matrix is not binary")
        table[offset] = block_index
        row_bits[row] |= 1 << column
    return row_bits, table


def set_bits(value: int) -> list[int]:
    result: list[int] = []
    while value:
        low = value & -value
        result.append(low.bit_length() - 1)
        value ^= low
    return result


def enumerate_rectangles(
    labels: list[tuple[int, int, int]],
    row_band: int,
    column_band: int,
    sparse_band: int,
) -> tuple[int, int, int]:
    """Count rectangles by the number of colliding third-band diagonals."""

    row_bits, table = projection_matrix(labels, row_band, column_band)
    by_diagonal_collisions = [0, 0, 0]
    for first_row in range(TILES):
        for second_row in range(first_row + 1, TILES):
            common_columns = set_bits(row_bits[first_row] & row_bits[second_row])
            for first_column, second_column in combinations(common_columns, 2):
                top_left = table[first_row * TILES + first_column]
                top_right = table[first_row * TILES + second_column]
                bottom_left = table[second_row * TILES + first_column]
                bottom_right = table[second_row * TILES + second_column]
                if min(top_left, top_right, bottom_left, bottom_right) < 0:
                    raise AssertionError("incomplete rectangle")

                diagonal_collisions = int(
                    labels[top_left][sparse_band]
                    == labels[bottom_right][sparse_band]
                ) + int(
                    labels[top_right][sparse_band]
                    == labels[bottom_left][sparse_band]
                )
                by_diagonal_collisions[diagonal_collisions] += 1
    return tuple(by_diagonal_collisions)


def main() -> None:
    labels = blocks()
    if len(labels) != TILES * LANES or len(set(labels)) != len(labels):
        raise AssertionError("block triples are not distinct")

    certify_balance(labels)
    certify_pair_capacity(labels)

    rectangle_counts: dict[int, tuple[int, int, int]] = {}
    for sparse_band in range(3):
        dense_bands = [band for band in range(3) if band != sparse_band]
        rectangle_counts[sparse_band] = enumerate_rectangles(
            labels, dense_bands[0], dense_bands[1], sparse_band
        )

    k4_minus = tuple(rectangle_counts[band][1] for band in range(3))
    pasch = sum(rectangle_counts[band][2] for band in range(3))
    if pasch != 0:
        raise AssertionError(f"found {pasch} Pasch/intercalate motifs")
    if k4_minus != EXPECTED_K4_MINUS_BY_SPARSE_BAND:
        raise AssertionError(
            f"K4-minus counts changed: {k4_minus} != "
            f"{EXPECTED_K4_MINUS_BY_SPARSE_BAND}"
        )

    print("fixed three-band tile-map certificate")
    print(f"tiles={TILES} lanes={LANES} blocks={len(labels)} slopes={SLOPES}")
    print("balance_per_band=64,64,64 PASS")
    print("pair_capacity=1,1,1 PASS")
    for sparse_band in range(3):
        ordinary, one_diagonal, two_diagonals = rectangle_counts[sparse_band]
        print(
            f"sparse_band={sparse_band} rectangles_by_diagonal_collisions="
            f"{ordinary},{one_diagonal},{two_diagonals}"
        )
    print(f"pasch_count={pasch} PASS")
    print(f"k4_minus_by_sparse_band={','.join(map(str, k4_minus))}")
    print(f"k4_minus_total={sum(k4_minus)} PASS")


if __name__ == "__main__":
    main()
