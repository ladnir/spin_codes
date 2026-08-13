#!/usr/bin/env python3
"""Exact averaged 41/43 split spectrum after one band-zero puncture.

Start from the exact 42/43 spectrum ``A[a,b]`` of the first-two-band EBCH
projection.  Choose one of the 42 band-zero coordinates uniformly and delete
it.  The number of (codeword, deleted-coordinate) pairs producing observed
weights ``a,b`` is

    P[a,b] = (42-a) A[a,b] + (a+1) A[a+1,b].

The membership probability of a fixed observed support after the random
within-band coordinate permutation and random puncture is therefore

    P[a,b] / (42 C(41,a) C(43,b)).

All construction and validation below use exact integers.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


N0_FULL = 42
N0 = 41
N1 = 43


def load_full(path: Path) -> list[list[int]]:
    table = [[0] * (N1 + 1) for _ in range(N0_FULL + 1)]
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            table[int(row["band0_weight"])][int(row["band1_weight"])] = int(
                row["count"]
            )
    if sum(map(sum, table)) != 1 << 64:
        raise SystemExit("punctured split: full spectrum mass mismatch")
    return table


def punctured_pairs(full: list[list[int]]) -> list[list[int]]:
    table = [[0] * (N1 + 1) for _ in range(N0 + 1)]
    for a in range(N0 + 1):
        for b in range(N1 + 1):
            table[a][b] = (N0_FULL - a) * full[a][b] + (a + 1) * full[a + 1][b]
    return table


def validate(table: list[list[int]]) -> None:
    if sum(map(sum, table)) != N0_FULL * (1 << 64):
        raise SystemExit("punctured split: pair mass mismatch")
    if table[0][0] != N0_FULL:
        raise SystemExit("punctured split: zero multiplicity mismatch")
    for a in range(N0 + 1):
        for b in range(N1 + 1):
            if table[a][b] != table[N0 - a][N1 - b]:
                raise SystemExit("punctured split: complement symmetry mismatch")


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
    root = Path(__file__).resolve().parent.parent
    parser.add_argument(
        "--full-spectrum",
        type=Path,
        default=root / "out" / "ebch85_band01_split_spectrum.csv",
    )
    parser.add_argument("--write-csv", type=Path)
    args = parser.parse_args()
    table = punctured_pairs(load_full(args.full_spectrum))
    validate(table)
    if args.write_csv is not None:
        args.write_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.write_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ("band0_weight", "band1_weight", "pair_count", "denominator")
            )
            for a, row in enumerate(table):
                for b, count in enumerate(row):
                    if count:
                        denominator = N0_FULL * math.comb(N0, a) * math.comb(N1, b)
                        writer.writerow((a, b, count, denominator))
    positive = [
        (a + b, a, b, count)
        for a, row in enumerate(table)
        for b, count in enumerate(row)
        if count and (a or b)
    ]
    print("exact averaged punctured EBCH band-(0,1) split spectrum")
    print("projection=[84,64] split=41/43 puncture_choices=42")
    print(f"minimum_pair_weight={min(positive)[0]}")
    print(f"nonzero_cells={sum(count != 0 for row in table for count in row)}")
    print(f"punctured_pair_sha256={digest(table)}")
    print("mass_and_complement_symmetry=PASS")
    print("status=EXACT_INTEGER_EBCH84_PUNCTURED_SPLIT_SPECTRUM")


if __name__ == "__main__":
    main()
