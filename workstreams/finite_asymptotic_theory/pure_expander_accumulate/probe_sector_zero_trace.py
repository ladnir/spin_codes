#!/usr/bin/env python3
"""Probe the scaled trace of the sector-zero covariance block.

The normalized Walsh transform preserves trace.  Since the sector-zero block
is positive semidefinite, every scaled primal diagonal is at most this trace.
This program evaluates the trace through the exact radial covariance formula
in extended precision.  It is a diagnostic, not an outward certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "sector_zero_trace_probe.json"


def radial_walk_table(
    message_bits: int, right_degree: int, output_bits: int
) -> np.ndarray:
    """Return per-vector XOR probabilities v_t(m) as long doubles."""
    dtype = np.longdouble
    denominator = dtype(math.comb(message_bits, right_degree))
    shell = np.zeros(message_bits + 1, dtype=dtype)
    shell[0] = 1
    result = np.zeros((output_bits + 1, message_bits + 1), dtype=dtype)
    result[0, 0] = 1
    multiplicities = np.asarray(
        [math.comb(message_bits, weight) for weight in range(message_bits + 1)],
        dtype=dtype,
    )
    for step in range(1, output_bits + 1):
        following = np.zeros_like(shell)
        for weight, mass in enumerate(shell):
            if mass == 0:
                continue
            low = max(0, right_degree - (message_bits - weight))
            high = min(weight, right_degree)
            for overlap in range(low, high + 1):
                following[weight + right_degree - 2 * overlap] += (
                    mass
                    * dtype(math.comb(weight, overlap))
                    * dtype(
                        math.comb(
                            message_bits - weight, right_degree - overlap
                        )
                    )
                    / denominator
                )
        shell = following
        result[step] = shell / multiplicities
    return result


def evaluate(message_bits: int, output_bits: int, right_degree: int) -> dict:
    dtype = np.longdouble
    walk = radial_walk_table(message_bits, right_degree, output_bits)
    message_multiplicities = np.asarray(
        [math.comb(message_bits, weight) for weight in range(message_bits + 1)],
        dtype=dtype,
    )
    diagonal_rows = []
    trace = dtype(0)
    minimum_diagonal = dtype(math.inf)
    for subset_weight in range(output_bits + 1):
        value = dtype(0)
        for intersection in range(subset_weight + 1):
            if subset_weight - intersection > output_bits - subset_weight:
                continue
            only = subset_weight - intersection
            joint = np.sum(
                message_multiplicities
                * walk[only]
                * walk[only]
                * walk[intersection],
                dtype=dtype,
            )
            covariance = joint - walk[subset_weight, 0] ** 2
            ways = dtype(math.comb(subset_weight, intersection)) * dtype(
                math.comb(output_bits - subset_weight, only)
            )
            value += ways * covariance
        scaled = dtype(2) ** message_bits * value
        trace += scaled
        minimum_diagonal = min(minimum_diagonal, scaled)
        diagonal_rows.append(
            {
                "dual_level": subset_weight,
                "scaled_sector_zero_diagonal": str(scaled),
            }
        )
    numerically_valid = bool(minimum_diagonal >= 0 and trace >= 0)
    return {
        "schema": "pure-ea-sector-zero-trace-probe-v1",
        "status": (
            "NONRIGOROUS_EXTENDED_PRECISION_DIAGNOSTIC"
            if numerically_valid
            else "NUMERICALLY_UNSTABLE_DIAGNOSTIC"
        ),
        "parameters": {
            "message_bits": message_bits,
            "output_bits": output_bits,
            "right_degree": right_degree,
        },
        "scaled_sector_zero_trace": str(trace),
        "minimum_scaled_dual_diagonal": str(minimum_diagonal),
        "all_computed_diagonals_nonnegative": numerically_valid,
        "trace_below_512_diagnostic": bool(
            numerically_valid and trace < 512
        ),
        "diagonals": diagonal_rows,
        "scope": [
            "The radial covariance and sector-zero trace formulas are exact.",
            "The numerical evaluation uses nondirected long-double arithmetic.",
            "A negative diagonal marks catastrophic cancellation and invalidates the reported trace comparison.",
            "A rigorous trace certificate requires outward interval arithmetic.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--right-degree", type=int, default=33)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(
        args.message_bits, args.output_bits, args.right_degree
    )
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
