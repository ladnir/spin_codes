#!/usr/bin/env python3
"""Certify when the diagonal Cauchy--Schwarz relaxation cannot reach its target.

This script does not lower-bound the true spectrum variance.  It lower-bounds
the right-hand side produced by the existing diagonal Cauchy--Schwarz route.
If the lower endpoint exceeds the target, that proof route cannot certify the
target without retaining additional signed cross-level information.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_dual_walk_and_accumulator_energy import accumulator_joint_counts
from certify_primal_schur_diagonal_outward import (
    add_lower_scalar,
    coefficient_lower,
    down,
    mean_shell_upper,
    positive_multiply_lower,
)


HERE = Path(__file__).resolve().parent
DEFAULT_RECEIPT = (
    HERE / "primal_schur_sector2_K256_B512_r33_p2_159_mass_rust20_outward.json"
)
DEFAULT_OUTPUT = (
    HERE / "diagonal_route_obstruction_K256_B512_r33_w42_79_outward.json"
)


def positive_sqrt_lower(value: np.float64) -> np.float64:
    if value < 0:
        raise ArithmeticError("cannot take the square root of a negative endpoint")
    if value == 0:
        return np.float64(0)
    return np.float64(max(0.0, down(np.float64(math.sqrt(float(value))))))


def positive_divide_lower(
    numerator: np.float64, denominator: np.float64
) -> np.float64:
    if denominator <= 0:
        raise ArithmeticError("denominator upper endpoint must be positive")
    return np.float64(max(0.0, down(np.float64(numerator / denominator))))


def obstruction_lower(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
    sector_two_lower: dict[int, np.float64],
) -> tuple[np.float64, np.float64, int]:
    counts = accumulator_joint_counts(output_bits)[shell_weight]
    total = np.float64(0)
    used = 0
    for level, count in sorted(counts.items()):
        if level < 2:
            continue
        if level not in sector_two_lower:
            raise ValueError(f"missing sector-two lower endpoint at level {level}")
        product = positive_multiply_lower(
            coefficient_lower(count), sector_two_lower[level]
        )
        total = add_lower_scalar(total, positive_sqrt_lower(product))
        used += 1
    quadratic = positive_multiply_lower(total, total)
    numerator = np.float64(math.ldexp(float(quadratic), -message_bits))
    numerator = np.float64(max(0.0, down(numerator)))
    mean_upper = mean_shell_upper(
        message_bits, output_bits, right_degree, shell_weight
    )
    return positive_divide_lower(numerator, mean_upper), mean_upper, used


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--shell-min", type=int, default=42)
    parser.add_argument("--shell-max", type=int, default=79)
    parser.add_argument("--target", type=float, default=512.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    payload = json.loads(args.receipt.read_text(encoding="utf-8"))
    if payload.get("status") != "OUTWARD_BINARY64_UPPER_BOUND":
        raise ValueError("input is not an outward orbit receipt")
    parameters = payload["parameters"]
    message_bits = int(parameters["message_bits"])
    output_bits = int(parameters["output_bits"])
    right_degree = int(parameters["right_degree"])
    lowers = {
        int(row["level"]): np.float64(row["positive_orbit_lower"])
        for row in payload["entries"]
        if int(row["sector"]) == 2
    }

    rows = []
    for shell_weight in range(args.shell_min, args.shell_max + 1):
        lower, mean_upper, used = obstruction_lower(
            message_bits,
            output_bits,
            right_degree,
            shell_weight,
            lowers,
        )
        rows.append(
            {
                "shell_weight": shell_weight,
                "sector_two_levels_used": used,
                "mean_upper": float(mean_upper),
                "mean_log2_upper": math.log2(float(mean_upper)),
                "diagonal_relaxation_ratio_lower": float(lower),
                "diagonal_relaxation_ratio_log2_lower": math.log2(float(lower)),
                "proved_above_target": bool(lower > args.target),
            }
        )

    first_failure = next(
        (row["shell_weight"] for row in rows if row["proved_above_target"]), None
    )
    result = {
        "schema": "pure-ea-diagonal-route-obstruction-outward-v1",
        "status": "OUTWARD_BINARY64_LOWER_BOUND",
        "parameters": {
            "message_bits": message_bits,
            "output_bits": output_bits,
            "right_degree": right_degree,
            "target": args.target,
        },
        "source_receipt": str(args.receipt.resolve()),
        "shells": rows,
        "claim": {
            "first_shell_proved_above_target": first_failure,
            "some_requested_shell_proved_above_target": first_failure is not None,
        },
        "scope": [
            "The lower bound applies to the diagonal Cauchy--Schwarz proof expression, not to the true variance.",
            "Only certified sector-two diagonal lower endpoints contribute; omitted sectors and levels contribute zero.",
            "A value above target proves that this relaxation cannot establish the requested factor without a sharper signed or correlated contraction.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
