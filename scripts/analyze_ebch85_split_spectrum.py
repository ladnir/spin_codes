#!/usr/bin/env python3
"""Exact 42/43 split spectrum of the first-two-band EBCH projection.

The projection of the committed EBCH [128,64,22] code onto coordinates 0..84
has rank 64 and dual dimension 21.  Enumerate all 2^21 dual words by Gray code
and apply the exact bivariate MacWilliams transform to obtain

    A[a,b] = #{ projected codewords with band weights (a,b) }.

All spectrum construction and validation use integers.  An optional CSV is a
reproducible exact artifact for the packet-constant configuration probe.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

from certify_bch_band_projections import BANDS, parity_check_columns
from probe_bch_forward_xor_circuit import target_outputs


N0 = 42
N1 = 43
DUAL_DIMENSION = N0 + N1 - 64


def parity_rows() -> list[int]:
    coordinates = list(range(BANDS[0][0], BANDS[0][1])) + list(
        range(BANDS[1][0], BANDS[1][1])
    )
    syndromes, rank = parity_check_columns(coordinates, target_outputs())
    if rank != 64 or max(syndromes).bit_length() != DUAL_DIMENSION:
        raise SystemExit("ebch85 split: projection rank changed")
    rows = []
    for check in range(DUAL_DIMENSION):
        row = 0
        for coordinate, syndrome in enumerate(syndromes):
            if syndrome >> check & 1:
                row |= 1 << coordinate
        rows.append(row)
    return rows


def dual_split_spectrum() -> list[list[int]]:
    spectrum = [[0] * (N1 + 1) for _ in range(N0 + 1)]
    rows = parity_rows()
    word = 0
    previous_gray = 0
    mask0 = (1 << N0) - 1
    for index in range(1 << DUAL_DIMENSION):
        gray = index ^ (index >> 1)
        if index:
            changed = gray ^ previous_gray
            word ^= rows[changed.bit_length() - 1]
        previous_gray = gray
        weight0 = (word & mask0).bit_count()
        weight1 = (word >> N0).bit_count()
        spectrum[weight0][weight1] += 1
    if sum(map(sum, spectrum)) != 1 << DUAL_DIMENSION:
        raise SystemExit("ebch85 split: dual spectrum mass mismatch")
    return spectrum


def krawtchouk(length: int) -> list[list[int]]:
    table = [[0] * (length + 1) for _ in range(length + 1)]
    for output_weight in range(length + 1):
        for input_weight in range(length + 1):
            table[output_weight][input_weight] = sum(
                (-1) ** intersection
                * math.comb(input_weight, intersection)
                * math.comb(length - input_weight, output_weight - intersection)
                for intersection in range(
                    max(0, output_weight - (length - input_weight)),
                    min(output_weight, input_weight) + 1,
                )
            )
    return table


def primal_split_spectrum(dual: list[list[int]]) -> list[list[int]]:
    k0 = krawtchouk(N0)
    k1 = krawtchouk(N1)
    nonzero_dual = [
        (weight0, weight1, count)
        for weight0, row in enumerate(dual)
        for weight1, count in enumerate(row)
        if count
    ]
    primal = [[0] * (N1 + 1) for _ in range(N0 + 1)]
    denominator = 1 << DUAL_DIMENSION
    for weight0 in range(N0 + 1):
        for weight1 in range(N1 + 1):
            numerator = sum(
                count * k0[weight0][dual0] * k1[weight1][dual1]
                for dual0, dual1, count in nonzero_dual
            )
            if numerator % denominator:
                raise SystemExit(
                    f"ebch85 split: nonintegral MacWilliams cell {weight0},{weight1}"
                )
            primal[weight0][weight1] = numerator // denominator
    return primal


def validate(primal: list[list[int]]) -> None:
    if sum(map(sum, primal)) != 1 << 64 or primal[0][0] != 1:
        raise SystemExit("ebch85 split: primal mass mismatch")
    for weight0 in range(N0 + 1):
        if sum(primal[weight0]) != math.comb(N0, weight0) << (64 - N0):
            raise SystemExit(f"ebch85 split: band-zero marginal {weight0} mismatch")
    for weight1 in range(N1 + 1):
        if sum(primal[weight0][weight1] for weight0 in range(N0 + 1)) != (
            math.comb(N1, weight1) << (64 - N1)
        ):
            raise SystemExit(f"ebch85 split: band-one marginal {weight1} mismatch")
    for weight0 in range(N0 + 1):
        for weight1 in range(N1 + 1):
            if primal[weight0][weight1] != primal[N0 - weight0][N1 - weight1]:
                raise SystemExit("ebch85 split: complement symmetry mismatch")


def digest(primal: list[list[int]]) -> str:
    value = hashlib.sha256()
    for weight0, row in enumerate(primal):
        for weight1, count in enumerate(row):
            if count:
                value.update(bytes((weight0, weight1)))
                value.update(count.to_bytes(16, "little"))
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-csv", type=Path)
    args = parser.parse_args()
    dual = dual_split_spectrum()
    primal = primal_split_spectrum(dual)
    validate(primal)
    if args.write_csv is not None:
        args.write_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.write_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(("band0_weight", "band1_weight", "count"))
            for weight0, row in enumerate(primal):
                for weight1, count in enumerate(row):
                    if count:
                        writer.writerow((weight0, weight1, count))
    positive = [
        (weight0 + weight1, weight0, weight1, count)
        for weight0, row in enumerate(primal)
        for weight1, count in enumerate(row)
        if count and (weight0 or weight1)
    ]
    print("exact EBCH first-two-band split spectrum")
    print(f"projection=[85,64] dual_dimension={DUAL_DIMENSION}")
    print(f"minimum_pair_weight={min(positive)[0]}")
    print(f"nonzero_cells={sum(count != 0 for row in primal for count in row)}")
    print(f"primal_split_sha256={digest(primal)}")
    print("mass_and_marginals=PASS")
    print("complement_symmetry=PASS")
    print("status=EXACT_INTEGER_EBCH85_SPLIT_SPECTRUM")


if __name__ == "__main__":
    main()
