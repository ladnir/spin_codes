#!/usr/bin/env python3
"""Monte Carlo variance probe for a small exact BA analogue.

The base code is the direct sum of four RM(1,3)=[8,4,4] constituents, so the
resulting [32,16] code has the same four-block architecture as EBCH128x4.
Every sampled code is enumerated exactly.  The experiment estimates shell
variance after successive uniform-interleaver accumulator stages and compares
it with the variance proxy Var(A_w)<=E[A_w] used by a uniform random generator.

This experiment is diagnostic.  It cannot prove or disprove the B=512 claim.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "small_ba32_spectrum_variance_probe.json"
B = 32
K = 16


def rm13_basis() -> list[int]:
    points = [(x >> 2 & 1, x >> 1 & 1, x & 1) for x in range(8)]
    rows = []
    for coordinate in range(4):
        word = 0
        for index, point in enumerate(points):
            value = 1 if coordinate == 0 else point[coordinate - 1]
            word |= value << index
        rows.append(word)
    return rows


def direct_sum_basis() -> np.ndarray:
    local = rm13_basis()
    return np.asarray(
        [word << (8 * block) for block in range(4) for word in local],
        dtype=np.uint32,
    )


def permute_word(word: int, permutation: np.ndarray) -> int:
    result = 0
    for output, source in enumerate(permutation):
        result |= ((word >> int(source)) & 1) << output
    return result


def accumulate_word(word: int) -> int:
    # Parallel prefix XOR in the low B bits.
    result = word
    shift = 1
    while shift < B:
        result ^= result << shift
        shift <<= 1
    return result & ((1 << B) - 1)


def transform_basis(basis: np.ndarray, permutation: np.ndarray) -> np.ndarray:
    return np.asarray(
        [accumulate_word(permute_word(int(word), permutation)) for word in basis],
        dtype=np.uint32,
    )


def half_table(basis: np.ndarray) -> np.ndarray:
    table = np.zeros(1 << 8, dtype=np.uint32)
    for value in range(1, 1 << 8):
        least = value & -value
        index = least.bit_length() - 1
        table[value] = table[value ^ least] ^ basis[index]
    return table


def spectrum(basis: np.ndarray, popcount16: np.ndarray) -> np.ndarray:
    left = half_table(basis[:8])
    right = half_table(basis[8:])
    words = left[:, None] ^ right[None, :]
    weights = popcount16[(words & 0xFFFF).astype(np.uint16)] + popcount16[
        (words >> 16).astype(np.uint16)
    ]
    return np.bincount(weights.ravel(), minlength=B + 1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2000)
    parser.add_argument("--maximum-accumulators", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0xBA512)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    popcount16 = np.asarray([value.bit_count() for value in range(1 << 16)], dtype=np.uint8)
    sums = np.zeros((args.maximum_accumulators + 1, B + 1), dtype=np.float64)
    sums2 = np.zeros_like(sums)
    maxima = np.zeros_like(sums)
    rng = np.random.default_rng(args.seed)
    base = direct_sum_basis()
    for sample in range(args.samples):
        basis = base.copy()
        for stage in range(args.maximum_accumulators + 1):
            counts = spectrum(basis, popcount16).astype(np.float64)
            sums[stage] += counts
            sums2[stage] += counts * counts
            maxima[stage] = np.maximum(maxima[stage], counts)
            if stage < args.maximum_accumulators:
                basis = transform_basis(basis, rng.permutation(B))
        if (sample + 1) % 100 == 0:
            print(f"samples,{sample + 1},{args.samples}", flush=True)

    rows = []
    for stage in range(args.maximum_accumulators + 1):
        mean = sums[stage] / args.samples
        variance = sums2[stage] / args.samples - mean * mean
        shell_rows = []
        for weight in range(1, B + 1):
            ratio = variance[weight] / mean[weight] if mean[weight] > 0 else math.nan
            shell_rows.append(
                {
                    "weight": weight,
                    "mean": float(mean[weight]),
                    "variance": float(max(0.0, variance[weight])),
                    "variance_over_mean": None if math.isnan(ratio) else float(ratio),
                    "maximum_observed": int(maxima[stage, weight]),
                }
            )
        finite = [row for row in shell_rows if row["variance_over_mean"] is not None]
        worst = max(finite, key=lambda row: float(row["variance_over_mean"]))
        central = [row for row in finite if 8 <= int(row["weight"]) <= 24]
        rows.append(
            {
                "accumulator_stages": stage,
                "worst_shell": worst,
                "worst_central_shell": max(
                    central, key=lambda row: float(row["variance_over_mean"])
                ),
                "shells": shell_rows,
            }
        )

    result = {
        "schema": "small-ba32-spectrum-variance-probe-v1",
        "status": "MONTE_CARLO_DIAGNOSTIC_WITH_EXACT_PER_SAMPLE_ENUMERATION",
        "parameters": {
            "base": "direct sum of four RM(1,3) [8,4,4] codes",
            "length": B,
            "dimension": K,
            "samples": args.samples,
            "seed": args.seed,
            "maximum_accumulator_stages": args.maximum_accumulators,
        },
        "stages": rows,
        "limitations": [
            "The sampler statistics are not confidence intervals.",
            "The B=32 analogue can expose a false universal variance lemma but cannot certify B=512 behavior.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            [
                {
                    "stage": row["accumulator_stages"],
                    "worst": row["worst_shell"],
                    "worst_central": row["worst_central_shell"],
                }
                for row in rows
            ],
            indent=2,
        )
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
