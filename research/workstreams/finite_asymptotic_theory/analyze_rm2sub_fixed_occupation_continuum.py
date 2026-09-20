#!/usr/bin/env python3
"""Evaluate the exact common-tilt continuum for fixed outer occupation.

The matrix formula is proved in RM2SUB_FIXED_OCCUPATION_CONTINUUM.md.  The
optimizer and special-function evaluations use binary64 arithmetic and are
diagnostic rather than outward-rounded certificates.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import hyp1f1, logsumexp


def log_beta_laplace(*, intervals: int, gamma: float) -> np.ndarray:
    """Return log E exp(-gamma X) for X~Beta(k,n-k), 0<=k<=n."""
    indices = np.arange(intervals + 1, dtype=np.float64)
    result = np.empty(intervals + 1, dtype=np.float64)
    result[0] = 0.0
    result[intervals] = -gamma
    if intervals > 1:
        values = hyp1f1(indices[1:intervals], intervals, -gamma)
        if np.any(values <= 0.0) or not np.all(np.isfinite(values)):
            raise FloatingPointError(
                "hyp1f1 left the positive finite range; lower the occupation cap"
            )
        result[1:intervals] = np.log(values)
    return result


def log_transfer_matrix(
    *, occupation: int, theta: float, state_bits: int
) -> np.ndarray:
    """Return entrywise natural logs of T_Q(theta)."""
    if occupation < 1:
        raise ValueError("occupation must be positive")
    state_count = (1 << state_bits) - 1
    reset = 1.0 / state_count
    live_one_probability = (1 << (state_bits - 1)) / state_count
    gamma = live_one_probability * theta
    intervals = occupation + 1
    log_integrals = log_beta_laplace(intervals=intervals, gamma=gamma)

    # At a potential row position, the zero bit applies I and the one bit
    # applies the reset transition P.  Summing both choices gives H=I+P.
    log_h = np.log(
        np.asarray(
            [[1.0, 1.0], [reset, 2.0 - reset]], dtype=np.float64
        )
    )
    output = np.empty((2, 2), dtype=np.float64)
    for initial_state in range(2):
        coefficients = np.full((2, intervals + 1), -math.inf)
        coefficients[initial_state, int(initial_state == 1)] = 0.0
        maximum_live_intervals = 1
        for _ in range(occupation):
            updated = np.full_like(coefficients, -math.inf)
            source = slice(0, maximum_live_intervals + 1)
            updated[0, source] = np.logaddexp(
                coefficients[0, source] + log_h[0, 0],
                coefficients[1, source] + log_h[1, 0],
            )
            destination = slice(1, maximum_live_intervals + 2)
            updated[1, destination] = np.logaddexp(
                coefficients[0, source] + log_h[0, 1],
                coefficients[1, source] + log_h[1, 1],
            )
            coefficients = updated
            maximum_live_intervals += 1
        for final_state in range(2):
            output[initial_state, final_state] = float(
                logsumexp(coefficients[final_state] + log_integrals)
            )
    return output


def log_perron(log_matrix: np.ndarray) -> float:
    scale = float(np.max(log_matrix))
    matrix = np.exp(log_matrix - scale)
    root = (
        matrix[0, 0]
        + matrix[1, 1]
        + math.sqrt(
            (matrix[0, 0] - matrix[1, 1]) ** 2
            + 4.0 * matrix[0, 1] * matrix[1, 0]
        )
    ) / 2.0
    return scale + math.log(root)


def optimize_occupation(
    *, occupation: int, delta: float, state_bits: int, block_constant: float
) -> dict[str, object]:
    state_count = (1 << state_bits) - 1
    live_one_probability = (1 << (state_bits - 1)) / state_count
    natural_tilt = occupation * math.log(2.0) / live_one_probability

    def objective(theta: float) -> float:
        log_radius = log_perron(
            log_transfer_matrix(
                occupation=occupation,
                theta=theta,
                state_bits=state_bits,
            )
        )
        return log_radius / math.log(2.0) + delta * theta / math.log(2.0)

    result = minimize_scalar(
        objective,
        bounds=(0.5 * natural_tilt, 3.0 * natural_tilt),
        method="bounded",
        options={"xatol": 1e-8, "maxiter": 100},
    )
    theta = float(result.x)
    objective_value = float(result.fun)
    gap = occupation / 2.0 - objective_value
    return {
        "occupation": occupation,
        "theta": theta,
        "theta_per_active_row": theta / occupation,
        "objective_bits_per_outer_coordinate": objective_value,
        "objective_per_active_row": objective_value / occupation,
        "gap_bits_per_outer_coordinate": gap,
        "gap_per_active_row": gap / occupation,
        "strict_log2_block_threshold": occupation / gap,
        "selected_block_constant": block_constant,
        "selected_schedule_decay_exponent": block_constant * gap - occupation,
        "selected_schedule_decay_exponent_per_active_row": (
            block_constant * gap / occupation - 1.0
        ),
        "optimizer_success": bool(result.success),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-bits", type=int, default=19)
    parser.add_argument("--delta", type=float, default=0.11)
    parser.add_argument("--block-constant", type=float, default=3.6)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=[1, 2, 3, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512],
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = [
        optimize_occupation(
            occupation=occupation,
            delta=args.delta,
            state_bits=args.state_bits,
            block_constant=args.block_constant,
        )
        for occupation in args.occupations
    ]
    first_failure = next(
        (
            row["occupation"]
            for row in rows
            if row["selected_schedule_decay_exponent"] <= 0.0
        ),
        None,
    )
    payload = {
        "schema": "rm2sub-fixed-occupation-continuum-v1",
        "status": "EXACT_LIMIT_FORMULA_BINARY64_COMMON_TILT_DIAGNOSTIC",
        "state_bits": args.state_bits,
        "delta": args.delta,
        "selected_block_constant": args.block_constant,
        "occupations": rows,
        "first_sampled_failure_of_selected_constant": first_failure,
        "maximum_sampled_threshold": max(
            float(row["strict_log2_block_threshold"]) for row in rows
        ),
        "limitations": [
            "The fixed-occupation continuum theorem is exact in form.",
            "The optimizer and hyp1f1 evaluations use binary64 arithmetic.",
            "Failure of this common-tilt sufficient condition is not a distance counterexample.",
            "The fixed-Q theorem is not uniform when Q grows with L.",
        ],
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
