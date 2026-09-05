#!/usr/bin/env python3
"""All-message finite diagnostic for a random outer and RandomStepConv g=1."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from analyze_riffle_striped_random_outer import log_matrix_power_moment  # noqa: E402


DEFAULT_OUTPUT = WORKSTREAM / "random_outer_randomstepconv_g1_k20_d11.json"
LOG2 = math.log(2.0)


def evaluate_memory(
    *, output_bits: int, message_bits: int, distance: int, memory_bits: int
) -> dict[str, float | int | bool]:
    state_zero = math.ldexp(1.0, -memory_bits)

    def objective(log_surprisal: float) -> float:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        live_moment = (1.0 + z) / 2.0
        terminate = state_zero * live_moment
        survive = (1.0 - state_zero) * live_moment

        # The random outer makes the current input bit independent and fair.
        # From state zero, input zero leaves the process at zero. Input one
        # invokes a random linear step on a nonzero vector.
        transition = (
            0.5 + 0.5 * terminate,
            0.5 * survive,
            terminate,
            survive,
        )
        inner_log = log_matrix_power_moment(transition, output_bits)
        outer_log = message_bits * LOG2 + math.log1p(-math.ldexp(1.0, -message_bits))
        return outer_log + inner_log + distance * surprisal

    grid = np.linspace(-4.0, 2.0, 97)
    values = np.asarray([objective(float(point)) for point in grid])
    index = int(np.argmin(values))
    lower = float(grid[max(0, index - 1)])
    upper = float(grid[min(len(grid) - 1, index + 1)])
    result = minimize_scalar(
        objective,
        bounds=(lower, upper),
        method="bounded",
        options={"xatol": 1.0e-13, "maxiter": 500},
    )
    log2_upper = float(result.fun) / LOG2
    return {
        "memory_bits": memory_bits,
        "log_surprisal": float(result.x),
        "surprisal": math.exp(float(result.x)),
        "z": math.exp(-math.exp(float(result.x))),
        "log2_expected_bad_upper": log2_upper,
        "margin_bits": -log2_upper,
        "closes_one_bit": log2_upper < -1.0,
        "closes_40_bits": log2_upper < -40.0,
        "optimizer_success": bool(result.success),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--output-bits", type=int, default=240 * 8832)
    parser.add_argument("--relative-distance", type=float, default=0.11)
    parser.add_argument(
        "--memory-bits",
        type=int,
        nargs="+",
        default=(8, 9, 10, 12, 14, 16, 18, 19, 20, 21, 22, 23),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    distance = math.floor(args.relative_distance * args.output_bits)
    rows = [
        evaluate_memory(
            output_bits=args.output_bits,
            message_bits=args.message_bits,
            distance=distance,
            memory_bits=memory_bits,
        )
        for memory_bits in args.memory_bits
    ]
    payload = {
        "schema": "random-outer-randomstepconv-g1-finite-v1",
        "status": "BINARY64_DIAGNOSTIC_EXACT_TWO_STATE_TRANSFER",
        "construction": {
            "outer": (
                "one uniform random binary linear map from k input bits to N "
                "outer bits"
            ),
            "inner": (
                "independent uniform unrestricted binary linear (M+1)-by-(M+1) "
                "map at every bit position"
            ),
            "recurrence": "(y_t,s_{t+1}) = M_t (x_t,s_t)",
            "initial_state": 0,
            "terminal_state": "discarded",
            "shared_setup": "all sampled maps are fixed and shared by every message",
        },
        "probability_space": (
            "the random outer generator and all independent RandomStepConv step maps"
        ),
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "relative_distance": args.relative_distance,
            "distance_cutoff": distance,
        },
        "rows": rows,
        "limitations": [
            "The random outer is a full random linear map, not the repeated BA outer.",
            "Arithmetic and tilt optimization use nearest binary64.",
            "The calculation proves no implementation-cost bound.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(rows, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
