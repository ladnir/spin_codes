#!/usr/bin/env python3
"""Diagnose extreme-shell failure for three sparse-mixer--accumulator stages.

This extends the one- and two-stage first-moment diagnostics.  The first map
is rectangular from 256 to 512 coordinates.  The next two maps are square.
All maps use independent fixed-size right neighborhoods and are sampled once.
The resulting composite is rank-tested before reuse in every outer row.

The binary64 receipt selects candidate parameters.  It is not an outward
certificate and does not address central-shell concentration.
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
from sweep_two_stage_extreme_shells import (
    first_stage_expected_spectrum,
    full_accumulator_shells,
    second_stage_tail_laws,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "three_stage_extreme_shell_sweep_B512.json"


def square_weight_kernel(output_bits: int, right_degree: int) -> np.ndarray:
    positive = activation_probabilities(output_bits, right_degree)
    activations = np.concatenate((np.array([0.0]), positive))
    return full_accumulator_shells(activations, output_bits)


def evaluate_candidate(
    spectrum_after_two: np.ndarray,
    output_bits: int,
    degrees: tuple[int, int, int],
    final_laws: tuple[np.ndarray, np.ndarray, np.ndarray],
    rank_attempts: int,
) -> dict[str, object]:
    zero, nonzero_low, high = final_laws
    kernel_expected = float(spectrum_after_two @ zero)
    low_expected = float(spectrum_after_two @ nonzero_low)
    high_expected = float(spectrum_after_two @ high)
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
    first_degree, second_degree, third_degree = degrees
    return {
        "first_right_degree": first_degree,
        "second_right_degree": second_degree,
        "third_right_degree": third_degree,
        "degree_sum": sum(degrees),
        "total_xors_per_constituent": output_bits * sum(degrees) - 3,
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
    parser.add_argument("--first-degree-max", type=int, default=19)
    parser.add_argument("--second-degree-min", type=int, default=1)
    parser.add_argument("--second-degree-max", type=int, default=9)
    parser.add_argument("--third-degree-min", type=int, default=1)
    parser.add_argument("--third-degree-max", type=int, default=11)
    parser.add_argument("--rank-attempts", type=int, default=16)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.message_bits != 256 or args.output_bits != 512:
        parser.error("the current diagnostic is fixed to the [512,256] target")
    if not (0 <= args.low_cutoff < args.high_cutoff <= args.output_bits):
        parser.error("invalid extreme-shell cutoffs")
    if args.rank_attempts < 1:
        parser.error("rank attempts must be positive")

    first_degrees = [
        degree
        for degree in range(args.first_degree_min, args.first_degree_max + 1)
        if degree % 2 == 1
    ]
    second_degrees = list(range(args.second_degree_min, args.second_degree_max + 1))
    third_degrees = list(range(args.third_degree_min, args.third_degree_max + 1))

    middle_kernels = {
        degree: square_weight_kernel(args.output_bits, degree)
        for degree in second_degrees
    }
    final_laws = {
        degree: second_stage_tail_laws(
            args.output_bits,
            args.low_cutoff,
            args.high_cutoff,
            degree,
        )
        for degree in third_degrees
    }

    rows = []
    for first_degree in first_degrees:
        first_spectrum = first_stage_expected_spectrum(
            args.message_bits, args.output_bits, first_degree
        )
        for second_degree in second_degrees:
            second_spectrum = first_spectrum @ middle_kernels[second_degree]
            for third_degree in third_degrees:
                rows.append(
                    evaluate_candidate(
                        second_spectrum,
                        args.output_bits,
                        (first_degree, second_degree, third_degree),
                        final_laws[third_degree],
                        args.rank_attempts,
                    )
                )

    rows.sort(
        key=lambda row: (
            row["degree_sum"],
            row["first_right_degree"],
            row["second_right_degree"],
        )
    )
    passing = [
        row
        for row in rows
        if row["rank_tested_setup_failure_bits"] is not None
        and row["rank_tested_setup_failure_bits"] > 42.2615
    ]
    best = None if not passing else passing[0]
    payload = {
        "schema": "pure-ea-three-stage-extreme-shell-sweep-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": vars(args) | {"output": str(args.output)},
        "first_passing_by_degree_sum": best,
        "rows": rows,
        "scope": [
            "All sparse maps are sampled once and reused in every outer row.",
            "The first moment includes rank loss of the complete composite map.",
            "Rank-tested resampling accepts only a full-rank composite map.",
            "Binary64 values select parameters but do not certify inequalities.",
            "This receipt controls only the two extreme weight ranges.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"first_passing_by_degree_sum": best}, indent=2))


if __name__ == "__main__":
    main()
