#!/usr/bin/env python3
"""Diagnose extreme-shell failure for two sparse-mixer--accumulator stages.

The first map has shape ``message_bits`` by ``output_bits``.  The second is
square on ``output_bits`` coordinates.  Each output coordinate independently
samples a fixed-size input neighborhood.  Setup samples both maps once.

The calculation propagates the exact one-word weight law in binary64.  It
then applies a first-moment rank and extreme-shell bound.  The receipt is a
parameter diagnostic, not an outward certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from sweep_one_stage_extreme_shells import (
    activation_probabilities,
    extreme_tail_probabilities,
    finite_or_none,
    safe_bits,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "two_stage_extreme_shell_sweep_B512.json"


def full_accumulator_shells(
    activations: np.ndarray, output_bits: int
) -> np.ndarray:
    rows = activations.shape[0]
    state_zero = np.zeros((rows, output_bits + 1), dtype=np.float64)
    state_one = np.zeros((rows, output_bits + 1), dtype=np.float64)
    state_zero[:, 0] = 1.0
    q = activations[:, None]
    one_minus_q = 1.0 - q
    for _ in range(output_bits):
        next_zero = one_minus_q * state_zero + q * state_one
        next_one = np.zeros_like(state_one)
        next_one[:, 1:] = (
            q * state_zero[:, :-1] + one_minus_q * state_one[:, :-1]
        )
        state_zero, state_one = next_zero, next_one
    result = state_zero + state_one
    residual = np.max(np.abs(np.sum(result, axis=1) - 1.0))
    if residual > 2e-13:
        raise AssertionError(f"weight-law mass residual {residual} is too large")
    return result


def first_stage_expected_spectrum(
    message_bits: int, output_bits: int, right_degree: int
) -> np.ndarray:
    activations = activation_probabilities(message_bits, right_degree)
    shells = full_accumulator_shells(activations, output_bits)
    counts = np.array(
        [float(math.comb(message_bits, weight)) for weight in range(1, message_bits + 1)],
        dtype=np.float64,
    )
    return counts @ shells


def second_stage_tail_laws(
    output_bits: int,
    low_cutoff: int,
    high_cutoff: int,
    right_degree: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    positive_activations = activation_probabilities(output_bits, right_degree)
    activations = np.concatenate((np.array([0.0]), positive_activations))
    low, high = extreme_tail_probabilities(
        activations, output_bits, max(low_cutoff, output_bits - high_cutoff)
    )
    zero = (1.0 - activations) ** output_bits
    nonzero_low = np.maximum(0.0, low - zero)
    return zero, nonzero_low, high


def evaluate_pair(
    first_spectrum: np.ndarray,
    output_bits: int,
    first_degree: int,
    second_degree: int,
    second_laws: tuple[np.ndarray, np.ndarray, np.ndarray],
    rank_attempts: int,
) -> dict[str, object]:
    zero, nonzero_low, high = second_laws
    kernel_expected = float(first_spectrum @ zero)
    low_expected = float(first_spectrum @ nonzero_low)
    high_expected = float(first_spectrum @ high)
    nonzero_extreme = low_expected + high_expected
    rank_failure_upper = min(1.0, kernel_expected)
    rank_success_lower = max(0.0, 1.0 - kernel_expected)
    accepted_bad_upper = (
        math.inf
        if rank_success_lower == 0.0
        else nonzero_extreme / rank_success_lower
    )
    abort_upper = rank_failure_upper**rank_attempts
    setup_upper = accepted_bad_upper + abort_upper
    return {
        "first_right_degree": first_degree,
        "second_right_degree": second_degree,
        "degree_sum": first_degree + second_degree,
        "first_expected_left_degree": output_bits * first_degree / 256,
        "second_expected_left_degree": second_degree,
        "total_xors_per_constituent": output_bits * (first_degree + second_degree) - 2,
        "kernel_expected_count": kernel_expected,
        "nonzero_low_expected_count": low_expected,
        "high_expected_count": high_expected,
        "nonzero_extreme_expected_count": nonzero_extreme,
        "rank_test_abort_upper": abort_upper,
        "accepted_nonzero_extreme_upper": finite_or_none(accepted_bad_upper),
        "rank_tested_setup_failure_upper": finite_or_none(setup_upper),
        "rank_tested_setup_failure_bits": safe_bits(setup_upper),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--low-cutoff", type=int, default=41)
    parser.add_argument("--high-cutoff", type=int, default=471)
    parser.add_argument("--first-degree-min", type=int, default=3)
    parser.add_argument("--first-degree-max", type=int, default=31)
    parser.add_argument("--second-degree-min", type=int, default=1)
    parser.add_argument("--second-degree-max", type=int, default=31)
    parser.add_argument("--rank-attempts", type=int, default=16)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.message_bits != 256:
        parser.error("the current XOR report assumes message_bits=256")
    if not (0 <= args.low_cutoff < args.high_cutoff <= args.output_bits):
        parser.error("invalid extreme-shell cutoffs")
    if args.rank_attempts < 1:
        parser.error("rank attempts must be positive")

    first_degrees = [
        degree
        for degree in range(args.first_degree_min, args.first_degree_max + 1)
        if degree % 2 == 1
    ]
    # The rectangular first map needs odd right degree because an even degree
    # kills the all-one message deterministically.  A later square map may use
    # even degree: the preceding code need not contain the all-one word, and
    # the complete composite rank test catches an actual kernel intersection.
    second_degrees = list(
        range(args.second_degree_min, args.second_degree_max + 1)
    )
    second_laws = {
        degree: second_stage_tail_laws(
            args.output_bits,
            args.low_cutoff,
            args.high_cutoff,
            degree,
        )
        for degree in second_degrees
    }

    rows = []
    for first_degree in first_degrees:
        first_spectrum = first_stage_expected_spectrum(
            args.message_bits, args.output_bits, first_degree
        )
        for second_degree in second_degrees:
            rows.append(
                evaluate_pair(
                    first_spectrum,
                    args.output_bits,
                    first_degree,
                    second_degree,
                    second_laws[second_degree],
                    args.rank_attempts,
                )
            )

    rows.sort(key=lambda row: (row["degree_sum"], row["first_right_degree"]))
    passing = [
        row
        for row in rows
        if row["rank_tested_setup_failure_bits"] is not None
        and row["rank_tested_setup_failure_bits"] > 42.2615
    ]
    best = None if not passing else passing[0]
    payload = {
        "schema": "pure-ea-two-stage-extreme-shell-sweep-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": vars(args) | {"output": str(args.output)},
        "first_passing_by_degree_sum": best,
        "rows": rows,
        "scope": [
            "Both sparse maps are sampled once and reused in every outer row.",
            "The first moment includes rank loss of the complete two-stage map.",
            "Rank-tested resampling accepts only a full-rank composite map.",
            "Binary64 values select parameters but do not certify inequalities.",
            "This receipt controls only the two extreme weight ranges.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"first_passing_by_degree_sum": best}, indent=2))


if __name__ == "__main__":
    main()
