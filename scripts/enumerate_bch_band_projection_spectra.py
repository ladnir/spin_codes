#!/usr/bin/env python3
"""Exact spectra of the three fixed-band complement projections.

Each 85/86-coordinate complement projection has dimension 64 and dual
dimension only 21/22.  Enumerate the dual by Gray code, then apply the binary
MacWilliams transform to obtain the complete primal spectrum exactly.
"""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path

from certify_bch_band_projections import BANDS, parity_check_columns
from probe_bch_forward_xor_circuit import DIMENSION, LENGTH, target_outputs


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "ebch128_fixed_band_projection_spectra.csv"


def parity_rows(syndromes: list[int]) -> list[int]:
    codimension = max(syndromes, default=0).bit_length()
    return [
        sum(((syndrome >> bit) & 1) << column for column, syndrome in enumerate(syndromes))
        for bit in range(codimension)
    ]


def dual_spectrum(rows: list[int], length: int) -> list[int]:
    spectrum = [0] * (length + 1)
    word = 0
    previous_gray = 0
    spectrum[0] = 1
    for index in range(1, 1 << len(rows)):
        gray = index ^ (index >> 1)
        changed = gray ^ previous_gray
        word ^= rows[changed.bit_length() - 1]
        spectrum[word.bit_count()] += 1
        previous_gray = gray
    return spectrum


def krawtchouk(length: int, output_weight: int, input_weight: int) -> int:
    return sum(
        (-1) ** selected
        * math.comb(input_weight, selected)
        * math.comb(length - input_weight, output_weight - selected)
        for selected in range(
            max(0, output_weight - (length - input_weight)),
            min(output_weight, input_weight) + 1,
        )
    )


def macwilliams(dual: list[int], dimension: int) -> list[int]:
    length = len(dual) - 1
    denominator = sum(dual)
    primal = []
    for output_weight in range(length + 1):
        numerator = sum(
            count * krawtchouk(length, output_weight, input_weight)
            for input_weight, count in enumerate(dual)
            if count
        )
        quotient, remainder = divmod(numerator, denominator)
        if remainder:
            raise SystemExit("band spectrum: nonintegral MacWilliams coefficient")
        primal.append(quotient)
    if sum(primal) != 1 << dimension or primal[0] != 1 or any(value < 0 for value in primal):
        raise SystemExit("band spectrum: invalid primal spectrum")
    return primal


def main() -> None:
    outputs = target_outputs()
    rows = []
    for band, (begin, end) in enumerate(BANDS):
        outside = list(range(0, begin)) + list(range(end, LENGTH))
        syndromes, rank = parity_check_columns(outside, outputs)
        if rank != DIMENSION:
            raise SystemExit("band spectrum: outside projection lost rank")
        checks = parity_rows(syndromes)
        if len(checks) != len(outside) - DIMENSION:
            raise SystemExit("band spectrum: parity-check dimension mismatch")
        dual = dual_spectrum(checks, len(outside))
        primal = macwilliams(dual, DIMENSION)
        minimum = next(weight for weight, count in enumerate(primal) if weight and count)
        for weight, count in enumerate(primal):
            if count:
                rows.append((band, "primal", weight, count))
        for weight, count in enumerate(dual):
            if count:
                rows.append((band, "dual", weight, count))
        print(
            f"band={band} outside_length={len(outside)} dual_dimension={len(checks)} "
            f"outside_distance={minimum}"
        )

    with DEFAULT_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("band", "code", "weight", "count"))
        writer.writerows(rows)
    digest = hashlib.sha256(DEFAULT_OUTPUT.read_bytes()).hexdigest()
    print(f"output={DEFAULT_OUTPUT}")
    print(f"sha256={digest}")
    print("status=EXACT_DUAL_ENUMERATION_AND_MACWILLIAMS")


if __name__ == "__main__":
    main()
