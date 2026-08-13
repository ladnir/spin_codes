#!/usr/bin/env python3
"""Enumerate the exact weight spectrum of the 24-bit EBCH graph block.

The graph check r=Rm is embedded as the first 24 message coordinates of the
committed cyclic [128,64,22] encoder, with the remaining 40 message
coordinates zero.  Gray-code enumeration visits all 2^24 graph words with one
generator-row XOR per word.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

from analyze_nosinger_bch_kernel import generator_rows


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "ebch128_graph24_spectrum.csv"
DIMENSION = 24


def enumerate_spectrum() -> list[int]:
    rows = generator_rows()[:DIMENSION]
    counts = [0] * 129
    word = 0
    counts[0] = 1
    for index in range(1, 1 << DIMENSION):
        word ^= rows[(index & -index).bit_length() - 1]
        counts[word.bit_count()] += 1
    if sum(counts) != 1 << DIMENSION:
        raise SystemExit("graph24 spectrum: cardinality mismatch")
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    rows = generator_rows()[:DIMENSION]
    digest = hashlib.sha256(
        b"".join(row.to_bytes(16, "little") for row in rows)
    ).hexdigest()
    counts = enumerate_spectrum()
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("weight", "count"))
        writer.writerows((weight, count) for weight, count in enumerate(counts) if count)
    minimum = next(weight for weight, count in enumerate(counts) if weight and count)
    print("graph24_spectrum_status=EXACT_GRAY_ENUMERATION")
    print(f"generator_rows_sha256={digest}")
    print(f"dimension={DIMENSION}")
    print(f"cardinality={sum(counts)}")
    print(f"minimum_distance={minimum}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
