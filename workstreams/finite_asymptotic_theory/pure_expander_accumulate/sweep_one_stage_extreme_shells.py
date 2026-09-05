#!/usr/bin/env python3
"""Diagnose extreme-shell failure for one independent-right-regular EA stage.

For each requested odd right degree, this script evaluates the binary64 first
moment of nonzero messages whose accumulated output has weight at most the
low cutoff or at least the high cutoff.  Markov's inequality turns an exact
evaluation of this moment into a spectrum-event failure bound.  The present
binary64 calculation is a parameter diagnostic, not an outward certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "one_stage_extreme_shell_sweep_B512.json"


def krawtchouk(length: int, degree: int, weight: int) -> int:
    return sum(
        (-1) ** j
        * math.comb(weight, j)
        * math.comb(length - weight, degree - j)
        for j in range(max(0, degree - (length - weight)), min(degree, weight) + 1)
    )


def activation_probabilities(message_bits: int, right_degree: int) -> np.ndarray:
    denominator = math.comb(message_bits, right_degree)
    return np.array(
        [
            0.5
            * (
                1.0
                - krawtchouk(message_bits, right_degree, weight) / denominator
            )
            for weight in range(1, message_bits + 1)
        ],
        dtype=np.float64,
    )


def extreme_tail_probabilities(
    activations: np.ndarray, output_bits: int, tail_width: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return probabilities of <=tail_width ones and <=tail_width zeros."""
    rows = activations.shape[0]
    shape = (rows, tail_width + 1)
    low_zero = np.zeros(shape, dtype=np.float64)
    low_one = np.zeros(shape, dtype=np.float64)
    high_zero = np.zeros(shape, dtype=np.float64)
    high_one = np.zeros(shape, dtype=np.float64)
    low_zero[:, 0] = 1.0
    high_zero[:, 0] = 1.0
    q = activations[:, None]
    one_minus_q = 1.0 - q

    for _ in range(output_bits):
        next_low_zero = one_minus_q * low_zero + q * low_one
        next_low_one = np.zeros(shape, dtype=np.float64)
        next_low_one[:, 1:] = (
            q * low_zero[:, :-1] + one_minus_q * low_one[:, :-1]
        )

        next_high_zero = np.zeros(shape, dtype=np.float64)
        next_high_zero[:, 1:] = (
            one_minus_q * high_zero[:, :-1] + q * high_one[:, :-1]
        )
        next_high_one = q * high_zero + one_minus_q * high_one

        low_zero, low_one = next_low_zero, next_low_one
        high_zero, high_one = next_high_zero, next_high_one

    low = np.sum(low_zero + low_one, axis=1)
    high = np.sum(high_zero + high_one, axis=1)
    return low, high


def safe_bits(value: float) -> float | None:
    if not math.isfinite(value):
        return None
    if value == 0.0:
        return None
    return -math.log2(value)


def finite_or_none(value: float) -> float | None:
    return value if math.isfinite(value) else None


def evaluate_degree(
    message_bits: int,
    output_bits: int,
    low_cutoff: int,
    high_cutoff: int,
    right_degree: int,
    rank_attempts: int,
) -> dict[str, object]:
    activations = activation_probabilities(message_bits, right_degree)
    low, high = extreme_tail_probabilities(
        activations, output_bits, max(low_cutoff, output_bits - high_cutoff)
    )
    message_counts = np.array(
        [float(math.comb(message_bits, weight)) for weight in range(1, message_bits + 1)],
        dtype=np.float64,
    )
    zero = (1.0 - activations) ** output_bits
    nonzero_low = np.maximum(0.0, low - zero)
    kernel_terms = message_counts * zero
    low_terms = message_counts * nonzero_low
    high_terms = message_counts * high
    kernel_total = float(np.sum(kernel_terms))
    low_total = float(np.sum(low_terms))
    high_total = float(np.sum(high_terms))
    total = low_total + high_total
    rank_failure_upper = min(1.0, kernel_total)
    rank_success_lower = max(0.0, 1.0 - kernel_total)
    accepted_bad_upper = (
        math.inf if rank_success_lower == 0.0 else total / rank_success_lower
    )
    abort_upper = rank_failure_upper**rank_attempts
    tested_setup_upper = accepted_bad_upper + abort_upper
    largest_index = int(np.argmax(low_terms + high_terms))
    return {
        "right_degree": right_degree,
        "expected_left_degree": output_bits * right_degree / message_bits,
        "sparse_map_xors": output_bits * max(0, right_degree - 1),
        "accumulator_xors": output_bits - 1,
        "total_xors_per_constituent": output_bits * max(0, right_degree - 1)
        + output_bits
        - 1,
        "low_expected_count": low_total,
        "high_expected_count": high_total,
        "nonzero_extreme_expected_count": total,
        "unconditional_kernel_expected_count": kernel_total,
        "unconditional_rank_failure_upper": rank_failure_upper,
        "rank_success_lower": rank_success_lower,
        "rank_attempts": rank_attempts,
        "rank_test_abort_upper": abort_upper,
        "accepted_nonzero_extreme_upper": finite_or_none(accepted_bad_upper),
        "rank_tested_setup_failure_upper": finite_or_none(tested_setup_upper),
        "rank_tested_setup_failure_bits": safe_bits(tested_setup_upper),
        "largest_message_weight": largest_index + 1,
        "largest_weight_contribution": float(
            low_terms[largest_index] + high_terms[largest_index]
        ),
        "largest_weight_contribution_bits": safe_bits(
            float(low_terms[largest_index] + high_terms[largest_index])
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--low-cutoff", type=int, default=41)
    parser.add_argument("--high-cutoff", type=int, default=471)
    parser.add_argument("--degree-min", type=int, default=1)
    parser.add_argument("--degree-max", type=int, default=63)
    parser.add_argument("--rank-attempts", type=int, default=16)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not (0 <= args.low_cutoff < args.high_cutoff <= args.output_bits):
        parser.error("cutoffs must satisfy 0 <= low < high <= output bits")
    if not (1 <= args.degree_min <= args.degree_max <= args.message_bits):
        parser.error("degree range must lie in [1,message bits]")
    if args.rank_attempts < 1:
        parser.error("rank attempts must be positive")

    degrees = [
        degree
        for degree in range(args.degree_min, args.degree_max + 1)
        if degree % 2 == 1
    ]
    rows = [
        evaluate_degree(
            args.message_bits,
            args.output_bits,
            args.low_cutoff,
            args.high_cutoff,
            degree,
            args.rank_attempts,
        )
        for degree in degrees
    ]
    passing = [
        row
        for row in rows
        if row["rank_tested_setup_failure_bits"] is not None
        and row["rank_tested_setup_failure_bits"] > 42.2615
    ]
    payload = {
        "schema": "pure-ea-one-stage-extreme-shell-sweep-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "low_cutoff": args.low_cutoff,
            "high_cutoff": args.high_cutoff,
            "degree_min": args.degree_min,
            "degree_max": args.degree_max,
            "rank_attempts": args.rank_attempts,
            "target_outer_bits": 42.2615,
        },
        "first_passing_odd_degree": None if not passing else passing[0]["right_degree"],
        "rows": rows,
        "scope": [
            "The kernel first moment upper-bounds rejection of one rank-tested draw.",
            "The accepted bad-spectrum bound divides the nonzero-tail first moment by the rank-success lower bound.",
            "Binary64 values select parameters but do not certify inequalities.",
            "This receipt controls only the two extreme weight ranges.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "first_passing_odd_degree": payload["first_passing_odd_degree"],
                "rows": [
                    {
                        "right_degree": row["right_degree"],
                        "total_xors_per_constituent": row[
                            "total_xors_per_constituent"
                        ],
                        "kernel_expected_count": row[
                            "unconditional_kernel_expected_count"
                        ],
                        "rank_tested_setup_failure_bits": row[
                            "rank_tested_setup_failure_bits"
                        ],
                        "largest_message_weight": row["largest_message_weight"],
                    }
                    for row in rows
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
