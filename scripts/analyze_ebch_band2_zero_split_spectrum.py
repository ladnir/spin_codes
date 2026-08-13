#!/usr/bin/env python3
"""Exact 42/43 split spectrum of EBCH words zero on band two.

The committed EBCH [128,64,22] projection onto coordinates 85..127 has
rank 43, so its kernel has dimension 21.  Enumerate that message kernel by
Gray code, encode every word exactly, and count its weights in bands zero and
one.  The resulting CSV has the same schema as the full [85,64] projected
spectrum and can be consumed by the structured constant-packet bound.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

from certify_bch_band_projections import BANDS, encoded_word, row_reduce
from probe_bch_forward_xor_circuit import DIMENSION, target_outputs


KERNEL_DIMENSION = 21


def kernel_basis(rows: list[int]) -> list[int]:
    reduced, pivots = row_reduce(rows, DIMENSION)
    free = [column for column in range(DIMENSION) if column not in set(pivots)]
    basis = []
    for free_column in free:
        word = 1 << free_column
        for row, pivot in zip(reduced, pivots):
            if row >> free_column & 1:
                word |= 1 << pivot
        if any((row & word).bit_count() & 1 for row in rows):
            raise SystemExit("band2-zero spectrum: invalid kernel basis")
        basis.append(word)
    if len(basis) != KERNEL_DIMENSION:
        raise SystemExit("band2-zero spectrum: kernel dimension changed")
    return basis


def spectrum() -> list[list[int]]:
    outputs = target_outputs()
    begin2, end2 = BANDS[2]
    basis = kernel_basis(outputs[begin2:end2])
    n0 = BANDS[0][1] - BANDS[0][0]
    n1 = BANDS[1][1] - BANDS[1][0]
    result = [[0] * (n1 + 1) for _ in range(n0 + 1)]
    message = 0
    previous_gray = 0
    for index in range(1 << KERNEL_DIMENSION):
        gray = index ^ (index >> 1)
        if index:
            changed = gray ^ previous_gray
            message ^= basis[changed.bit_length() - 1]
        previous_gray = gray
        word = encoded_word(message, outputs)
        if word >> begin2:
            raise SystemExit("band2-zero spectrum: encoded word leaked into band two")
        weight0 = (word & ((1 << n0) - 1)).bit_count()
        weight1 = ((word >> BANDS[1][0]) & ((1 << n1) - 1)).bit_count()
        result[weight0][weight1] += 1
    if sum(map(sum, result)) != 1 << KERNEL_DIMENSION or result[0][0] != 1:
        raise SystemExit("band2-zero spectrum: mass mismatch")
    return result


def digest(table: list[list[int]]) -> str:
    value = hashlib.sha256()
    for weight0, row in enumerate(table):
        for weight1, count in enumerate(row):
            if count:
                value.update(bytes((weight0, weight1)))
                value.update(count.to_bytes(8, "little"))
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-csv", type=Path)
    args = parser.parse_args()
    table = spectrum()
    if args.write_csv is not None:
        args.write_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.write_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(("band0_weight", "band1_weight", "count"))
            for weight0, row in enumerate(table):
                for weight1, count in enumerate(row):
                    if count:
                        writer.writerow((weight0, weight1, count))
    positive = [
        (weight0 + weight1, weight0, weight1, count)
        for weight0, row in enumerate(table)
        for weight1, count in enumerate(row)
        if count and (weight0 or weight1)
    ]
    print("exact EBCH band-two-zero split spectrum")
    print(f"subcode_dimension={KERNEL_DIMENSION}")
    print(f"minimum_pair_weight={min(positive)[0]}")
    print(f"nonzero_cells={sum(count != 0 for row in table for count in row)}")
    print(f"split_sha256={digest(table)}")
    print("mass_and_zero_band=PASS")
    print("status=EXACT_INTEGER_EBCH_BAND2_ZERO_SPLIT_SPECTRUM")


if __name__ == "__main__":
    main()
