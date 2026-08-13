#!/usr/bin/env python3
"""Exact rank and outside-distance checks for the 42/43/43 BCH bands.

The fused outer kernel uses fixed contiguous coordinate bands, so its proof
must not substitute the hypergeometric law of a fresh random coordinate
split.  This script checks the exact binary projection maps and finds the
minimum nonzero weight in the complement of every band by meet-in-the-middle
parity-check search.  Every reported dependency is lifted back to a BCH
message and verified against the committed generator.
"""

from __future__ import annotations

import itertools

from probe_bch_forward_xor_circuit import DIMENSION, LENGTH, target_outputs


BANDS = ((0, 42), (42, 85), (85, 128))


def row_reduce(rows: list[int], columns: int) -> tuple[list[int], list[int]]:
    rows = [row for row in rows if row]
    pivots: list[int] = []
    rank = 0
    for column in range(columns):
        source = next(
            (index for index in range(rank, len(rows)) if rows[index] >> column & 1),
            None,
        )
        if source is None:
            continue
        rows[rank], rows[source] = rows[source], rows[rank]
        for index in range(len(rows)):
            if index != rank and rows[index] >> column & 1:
                rows[index] ^= rows[rank]
        pivots.append(column)
        rank += 1
        if rank == len(rows):
            break
    return rows[:rank], pivots


def projection_rows(coordinates: list[int], outputs: list[int]) -> list[int]:
    rows = []
    for message_bit in range(DIMENSION):
        row = 0
        for local, coordinate in enumerate(coordinates):
            if outputs[coordinate] >> message_bit & 1:
                row |= 1 << local
        rows.append(row)
    return rows


def parity_check_columns(
    coordinates: list[int], outputs: list[int]
) -> tuple[list[int], int]:
    reduced, pivots = row_reduce(projection_rows(coordinates, outputs), len(coordinates))
    pivot_set = set(pivots)
    free = [column for column in range(len(coordinates)) if column not in pivot_set]
    null_rows = []
    for free_column in free:
        row = 1 << free_column
        for reduced_row, pivot in zip(reduced, pivots):
            if reduced_row >> free_column & 1:
                row |= 1 << pivot
        null_rows.append(row)
    syndromes = [
        sum(((row >> column) & 1) << bit for bit, row in enumerate(null_rows))
        for column in range(len(coordinates))
    ]
    return syndromes, len(pivots)


def xor_syndrome(indices: tuple[int, ...], syndromes: list[int]) -> int:
    result = 0
    for index in indices:
        result ^= syndromes[index]
    return result


def minimum_dependency(
    syndromes: list[int], maximum_weight: int = 8
) -> tuple[int, ...]:
    columns = len(syndromes)
    for weight in range(1, maximum_weight + 1):
        left_weight = weight // 2
        right_weight = weight - left_weight
        left_sums: dict[int, tuple[int, ...]] = {}
        for indices in itertools.combinations(range(columns), left_weight):
            left_sums.setdefault(xor_syndrome(indices, syndromes), indices)
        for indices in itertools.combinations(range(columns), right_weight):
            other = left_sums.get(xor_syndrome(indices, syndromes))
            if other is not None and not set(other).intersection(indices):
                return other + indices
    raise SystemExit("band projection: no dependency found through requested weight")


def solve_message(
    coordinates: list[int], support: tuple[int, ...], outputs: list[int]
) -> int:
    support_set = set(support)
    equations = [
        outputs[coordinate] | ((local in support_set) << DIMENSION)
        for local, coordinate in enumerate(coordinates)
    ]
    pivot_rows: dict[int, int] = {}
    for equation in equations:
        value = equation
        while value & ((1 << DIMENSION) - 1):
            pivot = (value & ((1 << DIMENSION) - 1)).bit_length() - 1
            if pivot in pivot_rows:
                value ^= pivot_rows[pivot]
            else:
                pivot_rows[pivot] = value
                break
        else:
            if value >> DIMENSION:
                raise SystemExit("band projection: inconsistent dependency lift")

    message = 0
    for pivot in sorted(pivot_rows):
        equation = pivot_rows[pivot]
        lower = equation & ((1 << pivot) - 1)
        rhs = equation >> DIMENSION & 1
        if (lower & message).bit_count() & 1:
            rhs ^= 1
        if rhs:
            message |= 1 << pivot
    for local, coordinate in enumerate(coordinates):
        actual = (outputs[coordinate] & message).bit_count() & 1
        if actual != (local in support_set):
            raise SystemExit("band projection: lifted message does not match support")
    return message


def encoded_word(message: int, outputs: list[int]) -> int:
    return sum(
        (((form & message).bit_count() & 1) << coordinate)
        for coordinate, form in enumerate(outputs)
    )


def main() -> None:
    outputs = target_outputs()
    print("fixed 42/43/43 EBCH band projections")
    for band, (begin, end) in enumerate(BANDS):
        inside = list(range(begin, end))
        outside = list(range(0, begin)) + list(range(end, LENGTH))
        _inside_syndromes, inside_rank = parity_check_columns(inside, outputs)
        outside_syndromes, outside_rank = parity_check_columns(outside, outputs)
        dependency = minimum_dependency(outside_syndromes)
        message = solve_message(outside, dependency, outputs)
        word = encoded_word(message, outputs)
        outside_weight = sum((word >> coordinate) & 1 for coordinate in outside)
        inside_weight = sum((word >> coordinate) & 1 for coordinate in inside)
        if outside_weight != len(dependency) or message == 0:
            raise SystemExit("band projection: invalid minimum-distance witness")
        print(
            f"band={band} range=[{begin},{end}) inside_rank={inside_rank} "
            f"outside_rank={outside_rank} outside_distance={outside_weight} "
            f"witness_inside_weight={inside_weight} "
            f"witness_outside_coordinates="
            f"{tuple(outside[index] for index in dependency)}"
        )
    print("status=EXACT_BINARY_CERTIFICATE")


if __name__ == "__main__":
    main()
