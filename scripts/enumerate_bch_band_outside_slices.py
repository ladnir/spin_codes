#!/usr/bin/env python3
"""Enumerate low outside-weight slices of the fixed BCH coordinate bands.

For each 42/43-coordinate band, the complementary 85/86-coordinate
projection is an injective [n,64] code.  A parity-check meet-in-the-middle
enumerates every projected codeword through a requested outside weight and
records the corresponding weight inside the omitted band.  This supplies the
fixed-band boundary data needed by the three-band outer proof.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
from collections import Counter
from pathlib import Path

from certify_bch_band_projections import BANDS, parity_check_columns
from probe_bch_forward_xor_circuit import DIMENSION, LENGTH, target_outputs


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "ebch128_fixed_band_outside_slices.csv"


def independent_equations(
    coordinates: list[int], outputs: list[int]
) -> tuple[list[int], list[int]]:
    basis: dict[int, tuple[int, int]] = {}
    selected: list[int] = []
    for local, coordinate in enumerate(coordinates):
        value = outputs[coordinate]
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot][0]
            else:
                basis[pivot] = (value, local)
                selected.append(local)
                break
        if len(selected) == DIMENSION:
            break
    if len(selected) != DIMENSION:
        raise SystemExit("outside slices: complement projection is not injective")

    masks = [outputs[coordinates[local]] for local in selected]
    transforms = [1 << index for index in range(DIMENSION)]
    row = 0
    for column in range(DIMENSION):
        source = next(
            (index for index in range(row, DIMENSION) if masks[index] >> column & 1),
            None,
        )
        if source is None:
            raise SystemExit("outside slices: selected equations are singular")
        masks[row], masks[source] = masks[source], masks[row]
        transforms[row], transforms[source] = transforms[source], transforms[row]
        for index in range(DIMENSION):
            if index != row and masks[index] >> column & 1:
                masks[index] ^= masks[row]
                transforms[index] ^= transforms[row]
        row += 1
    inverse_rows = [0] * DIMENSION
    for mask, transform in zip(masks, transforms):
        if mask.bit_count() != 1:
            raise SystemExit("outside slices: inverse reduction did not reach identity")
        inverse_rows[mask.bit_length() - 1] = transform
    return selected, inverse_rows


def subset_records(
    syndromes: list[int], weight: int
) -> list[tuple[int, int]]:
    records = []
    for subset in itertools.combinations(range(len(syndromes)), weight):
        syndrome = 0
        mask = 0
        for index in subset:
            syndrome ^= syndromes[index]
            mask |= 1 << index
        records.append((syndrome, mask))
    records.sort(key=lambda row: row[0])
    return records


def syndrome_groups(records: list[tuple[int, int]]):
    begin = 0
    while begin < len(records):
        end = begin + 1
        while end < len(records) and records[end][0] == records[begin][0]:
            end += 1
        yield records[begin][0], records[begin:end]
        begin = end


def decode_message(
    support: int, selected: list[int], inverse_rows: list[int]
) -> int:
    values = 0
    for equation, outside_local in enumerate(selected):
        if support >> outside_local & 1:
            values |= 1 << equation
    return sum(
        (((inverse & values).bit_count() & 1) << message_bit)
        for message_bit, inverse in enumerate(inverse_rows)
    )


def enumerate_band(
    band: int, maximum_weight: int, outputs: list[int]
) -> Counter[tuple[int, int]]:
    begin, end = BANDS[band]
    inside = list(range(begin, end))
    outside = list(range(0, begin)) + list(range(end, LENGTH))
    syndromes, rank = parity_check_columns(outside, outputs)
    if rank != DIMENSION:
        raise SystemExit("outside slices: complement projection lost rank")
    selected, inverse_rows = independent_equations(outside, outputs)

    histogram: Counter[tuple[int, int]] = Counter()
    records = {
        weight: list(syndrome_groups(subset_records(syndromes, weight)))
        for weight in range((maximum_weight + 1) // 2 + 1)
    }
    for outside_weight in range(1, maximum_weight + 1):
        left_weight = outside_weight // 2
        right_weight = outside_weight - left_weight
        left = records[left_weight]
        right = records[right_weight]
        supports: set[int] = set()
        left_index = 0
        right_index = 0
        while left_index < len(left) and right_index < len(right):
            left_syndrome, left_rows = left[left_index]
            right_syndrome, right_rows = right[right_index]
            if left_syndrome < right_syndrome:
                left_index += 1
                continue
            if right_syndrome < left_syndrome:
                right_index += 1
                continue
            for _syndrome, left_mask in left_rows:
                for _syndrome, right_mask in right_rows:
                    if left_mask & right_mask:
                        continue
                    supports.add(left_mask | right_mask)
            left_index += 1
            right_index += 1
        for support in supports:
            message = decode_message(support, selected, inverse_rows)
            word = sum(
                (((form & message).bit_count() & 1) << coordinate)
                for coordinate, form in enumerate(outputs)
            )
            actual_outside = sum((word >> coordinate) & 1 for coordinate in outside)
            if actual_outside != outside_weight:
                raise SystemExit("outside slices: decoded word failed verification")
            inside_weight = sum((word >> coordinate) & 1 for coordinate in inside)
            histogram[(outside_weight, inside_weight)] += 1
    return histogram


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-weight", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 1 <= args.maximum_weight <= 8:
        raise SystemExit("outside slices: maximum weight must lie in 1..8")

    outputs = target_outputs()
    rows = []
    for band in range(len(BANDS)):
        histogram = enumerate_band(band, args.maximum_weight, outputs)
        for (outside_weight, inside_weight), count in sorted(histogram.items()):
            rows.append((band, outside_weight, inside_weight, count))
        print(
            f"band={band} codewords_through_weight_{args.maximum_weight}="
            f"{sum(histogram.values())}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("band", "outside_weight", "inside_weight", "count"))
        writer.writerows(rows)
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    print(f"output={args.output}")
    print(f"sha256={digest}")
    print("status=EXACT_MEET_IN_THE_MIDDLE_ENUMERATION")


if __name__ == "__main__":
    main()
