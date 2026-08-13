#!/usr/bin/env python3
"""Exact 43/43 split spectrum of the band-(1,2) EBCH projection."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

from analyze_ebch85_split_spectrum import krawtchouk
from certify_bch_band_projections import BANDS, parity_check_columns
from probe_bch_forward_xor_circuit import target_outputs


N0 = 43
N1 = 43
DIMENSION = 64
DUAL_DIMENSION = N0 + N1 - DIMENSION


def parity_rows() -> list[int]:
    coordinates = list(range(BANDS[1][0], BANDS[2][1]))
    syndromes, rank = parity_check_columns(coordinates, target_outputs())
    if rank != DIMENSION or max(syndromes).bit_length() != DUAL_DIMENSION:
        raise SystemExit("ebch86 band12: projection rank changed")
    return [
        sum(((syndrome >> check) & 1) << coordinate for coordinate, syndrome in enumerate(syndromes))
        for check in range(DUAL_DIMENSION)
    ]


def dual_spectrum() -> list[list[int]]:
    table = [[0] * (N1 + 1) for _ in range(N0 + 1)]
    rows = parity_rows()
    mask0 = (1 << N0) - 1
    word = 0
    previous_gray = 0
    for index in range(1 << DUAL_DIMENSION):
        gray = index ^ (index >> 1)
        if index:
            changed = gray ^ previous_gray
            word ^= rows[changed.bit_length() - 1]
        previous_gray = gray
        table[(word & mask0).bit_count()][(word >> N0).bit_count()] += 1
    if sum(map(sum, table)) != 1 << DUAL_DIMENSION:
        raise SystemExit("ebch86 band12: dual mass mismatch")
    return table


def primal_spectrum(dual: list[list[int]]) -> list[list[int]]:
    k0 = krawtchouk(N0)
    k1 = krawtchouk(N1)
    nonzero = [
        (a, b, count)
        for a, row in enumerate(dual)
        for b, count in enumerate(row)
        if count
    ]
    table = [[0] * (N1 + 1) for _ in range(N0 + 1)]
    denominator = 1 << DUAL_DIMENSION
    for a in range(N0 + 1):
        for b in range(N1 + 1):
            numerator = sum(
                count * k0[a][dual_a] * k1[b][dual_b]
                for dual_a, dual_b, count in nonzero
            )
            if numerator % denominator:
                raise SystemExit(f"ebch86 band12: nonintegral cell {a},{b}")
            table[a][b] = numerator // denominator
    return table


def validate(table: list[list[int]]) -> None:
    if sum(map(sum, table)) != 1 << DIMENSION or table[0][0] != 1:
        raise SystemExit("ebch86 band12: primal mass mismatch")
    marginal_factor = 1 << (DIMENSION - N0)
    for a in range(N0 + 1):
        if sum(table[a]) != math.comb(N0, a) * marginal_factor:
            raise SystemExit(f"ebch86 band12: first marginal {a} mismatch")
    for b in range(N1 + 1):
        if sum(table[a][b] for a in range(N0 + 1)) != (
            math.comb(N1, b) * marginal_factor
        ):
            raise SystemExit(f"ebch86 band12: second marginal {b} mismatch")


def digest(table: list[list[int]]) -> str:
    value = hashlib.sha256()
    for a, row in enumerate(table):
        for b, count in enumerate(row):
            if count:
                value.update(bytes((a, b)))
                value.update(count.to_bytes(16, "little"))
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-csv", type=Path)
    args = parser.parse_args()
    table = primal_spectrum(dual_spectrum())
    validate(table)
    if args.write_csv is not None:
        args.write_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.write_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(("band0_weight", "band1_weight", "count"))
            for a, row in enumerate(table):
                for b, count in enumerate(row):
                    if count:
                        writer.writerow((a, b, count))
    positive = [
        (a + b, a, b, count)
        for a, row in enumerate(table)
        for b, count in enumerate(row)
        if count and (a or b)
    ]
    print("exact EBCH band-(1,2) split spectrum")
    print(f"projection=[86,64] dual_dimension={DUAL_DIMENSION}")
    print(f"minimum_pair_weight={min(positive)[0]}")
    print(f"nonzero_cells={sum(count != 0 for row in table for count in row)}")
    print(f"primal_split_sha256={digest(table)}")
    print("mass_and_marginals=PASS")
    print("status=EXACT_INTEGER_EBCH86_BAND12_SPLIT_SPECTRUM")


if __name__ == "__main__":
    main()
