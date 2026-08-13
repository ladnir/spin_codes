#!/usr/bin/env python3
"""Exact fixed-band split rows for inside weights 0, 1 and complements.

For a band of size c, fixing its coordinate pattern leaves an affine coset of
dimension 64-c.  Enumerate the zero coset once, translate it by every unit
band pattern, and record the outside-weight distribution.  Complementing by
the all-ones BCH word also gives inside weights c and c-1.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np

from certify_bch_band_projections import BANDS, row_reduce, solve_message
from probe_bch_forward_xor_circuit import DIMENSION, LENGTH, target_outputs


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "ebch128_fixed_band_inside_boundary.csv"


def kernel_basis(equations: list[int]) -> list[int]:
    reduced, pivots = row_reduce(equations, DIMENSION)
    pivot_set = set(pivots)
    result = []
    for free in range(DIMENSION):
        if free in pivot_set:
            continue
        value = 1 << free
        for row, pivot in zip(reduced, pivots):
            if row >> free & 1:
                value |= 1 << pivot
        result.append(value)
    return result


def encode(message: int, outputs: list[int]) -> int:
    return sum(
        (((form & message).bit_count() & 1) << coordinate)
        for coordinate, form in enumerate(outputs)
    )


def kernel_words(basis: list[int], outputs: list[int]) -> tuple[np.ndarray, np.ndarray]:
    rows = [encode(message, outputs) for message in basis]
    count = 1 << len(rows)
    low = np.empty(count, dtype=np.uint64)
    high = np.empty(count, dtype=np.uint64)
    low[0] = 0
    high[0] = 0
    word = 0
    previous_gray = 0
    for index in range(1, count):
        gray = index ^ (index >> 1)
        changed = gray ^ previous_gray
        word ^= rows[changed.bit_length() - 1]
        low[index] = word & ((1 << 64) - 1)
        high[index] = word >> 64
        previous_gray = gray
    return low, high


def outside_histogram(
    low: np.ndarray,
    high: np.ndarray,
    representative: int,
    outside_mask: int,
    outside_length: int,
) -> np.ndarray:
    low_mask = np.uint64(outside_mask & ((1 << 64) - 1))
    high_mask = np.uint64(outside_mask >> 64)
    representative_low = np.uint64(representative & ((1 << 64) - 1))
    representative_high = np.uint64(representative >> 64)
    weights = (
        np.bitwise_count((low ^ representative_low) & low_mask)
        + np.bitwise_count((high ^ representative_high) & high_mask)
    )
    return np.bincount(weights.astype(np.int64), minlength=outside_length + 1)


def main() -> None:
    outputs = target_outputs()
    rows = []
    for band, (begin, end) in enumerate(BANDS):
        inside = list(range(begin, end))
        outside = list(range(0, begin)) + list(range(end, LENGTH))
        basis = kernel_basis(outputs[begin:end])
        if len(basis) != DIMENSION - len(inside):
            raise SystemExit("inside boundary: unexpected kernel dimension")
        low, high = kernel_words(basis, outputs)
        outside_mask = sum(1 << coordinate for coordinate in outside)

        histograms: dict[int, np.ndarray] = {}
        histograms[0] = outside_histogram(
            low, high, 0, outside_mask, len(outside)
        )
        unit_histogram = np.zeros(len(outside) + 1, dtype=np.int64)
        for inside_local in range(len(inside)):
            message = solve_message(inside, (inside_local,), outputs)
            representative = encode(message, outputs)
            unit_histogram += outside_histogram(
                low, high, representative, outside_mask, len(outside)
            )
        histograms[1] = unit_histogram

        for inside_weight, histogram in histograms.items():
            for outside_weight, count in enumerate(histogram):
                if count:
                    rows.append((band, inside_weight, outside_weight, int(count)))
                    rows.append(
                        (
                            band,
                            len(inside) - inside_weight,
                            len(outside) - outside_weight,
                            int(count),
                        )
                    )
        print(
            f"band={band} inside_size={len(inside)} kernel_dimension={len(basis)} "
            f"zero_row_count={int(np.sum(histograms[0]))} "
            f"one_row_count={int(np.sum(histograms[1]))}"
        )

    rows.sort()
    with DEFAULT_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("band", "inside_weight", "outside_weight", "count"))
        writer.writerows(rows)
    digest = hashlib.sha256(DEFAULT_OUTPUT.read_bytes()).hexdigest()
    print(f"output={DEFAULT_OUTPUT}")
    print(f"sha256={digest}")
    print("status=EXACT_AFFINE_COSET_ENUMERATION")


if __name__ == "__main__":
    main()
