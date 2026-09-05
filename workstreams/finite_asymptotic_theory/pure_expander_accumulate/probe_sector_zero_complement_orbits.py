#!/usr/bin/env python3
"""Measure sector-zero cancellation after grouping message complements.

For odd sparse-row degree, complementing either message flips its complete
pre-accumulator output word.  Averaging the four corresponding pair kernels
does not change the total covariance.  This diagnostic measures whether that
exact four-term grouping makes a one-sided bound viable.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from probe_primal_schur_diagonal_target import centered_symmetric_diagonal
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "sector_zero_complement_orbit_probe.json"


def evaluate(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    levels: list[int],
    progress_every: int,
) -> dict[str, object]:
    if right_degree % 2 != 1:
        raise ValueError("complement grouping requires odd row degree")
    dtype = np.longdouble
    denominator = math.comb(message_bits, right_degree)
    biases = np.asarray(
        [
            dtype(krawtchouk(message_bits, right_degree, weight))
            / dtype(denominator)
            for weight in range(message_bits + 1)
        ],
        dtype=dtype,
    )
    scale = dtype(2) ** message_bits
    signed = {level: dtype(0) for level in levels}
    positive = {level: dtype(0) for level in levels}
    absolute = {level: dtype(0) for level in levels}
    types = 0
    started = time.perf_counter()
    for n11 in range(message_bits + 1):
        n01_rows: list[int] = []
        n10_rows: list[int] = []
        for n01 in range(message_bits - n11 + 1):
            for n10 in range(message_bits - n11 - n01 + 1):
                n01_rows.append(n01)
                n10_rows.append(n10)
        n01 = np.asarray(n01_rows, dtype=np.int16)
        n10 = np.asarray(n10_rows, dtype=np.int16)
        h1 = n10 + n11
        h2 = n01 + n11
        h3 = n01 + n10
        b1, b2, b3 = biases[h1], biases[h2], biases[h3]
        probability = (
            (1 + b1 + b2 + b3) / 4,
            (1 + b1 - b2 - b3) / 4,
            (1 - b1 + b2 - b3) / 4,
            (1 - b1 - b2 + b3) / 4,
        )
        determinant = (b3 - b1 * b2) / 4
        multiplicities = np.asarray(
            [
                math.comb(message_bits, n11)
                * math.comb(message_bits - n11, int(d01))
                * math.comb(message_bits - n11 - int(d01), int(d10))
                for d01, d10 in zip(n01, n10, strict=True)
            ],
            dtype=dtype,
        )
        type_scale = multiplicities * scale
        a, b, c, d = probability
        transforms = (
            (probability, determinant),
            ((c, d, a, b), -determinant),
            ((b, a, d, c), -determinant),
            ((d, c, b, a), determinant),
        )
        for level in levels:
            orbit = np.zeros_like(a)
            for transformed, transformed_determinant in transforms:
                centered, _absolute, _telescoping, _rational = (
                    centered_symmetric_diagonal(
                        transformed,
                        transformed_determinant,
                        output_bits,
                        level,
                    )
                )
                orbit += centered
            contribution = type_scale * orbit / 4
            signed[level] += np.sum(contribution, dtype=dtype)
            positive[level] += np.sum(np.maximum(contribution, 0), dtype=dtype)
            absolute[level] += np.sum(np.abs(contribution), dtype=dtype)
        types += len(n01)
        if progress_every and (
            n11 % progress_every == 0 or n11 == message_bits
        ):
            print(
                f"progress,n11,{n11},{message_bits},pair_types,{types},elapsed_seconds,{time.perf_counter()-started:.3f}",
                flush=True,
            )
    return {
        "message_bits": message_bits,
        "output_bits": output_bits,
        "right_degree": right_degree,
        "enumerated_pair_types": types,
        "entries": [
            {
                "level": level,
                "scaled_orbit_grouped_sum": float(signed[level]),
                "scaled_orbit_grouped_positive_part": float(positive[level]),
                "scaled_orbit_grouped_absolute_sum": float(absolute[level]),
            }
            for level in levels
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--right-degree", type=int, default=33)
    parser.add_argument("--level", type=int, action="append", required=True)
    parser.add_argument("--progress-every", type=int, default=16)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = evaluate(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        sorted(set(args.level)),
        args.progress_every,
    )
    payload = {
        "schema": "pure-ea-sector-zero-complement-orbit-probe-v1",
        "status": "NONDIRECTED_FLOAT_DIAGNOSTIC",
        "result": result,
        "scope": [
            "The four-term complement grouping is an exact identity for odd sparse-row degree.",
            "The displayed sums use nondirected floating-point arithmetic and are not certificate bounds.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
