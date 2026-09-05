#!/usr/bin/env python3
"""Find the BA-2 variance factor allowed by the Q=1 finite transfer.

The calculation reuses one RandomStepConv Q=1 kernel for many hypothetical
shell-variance factors and forced-zero cutoffs. It is a nearest-binary64
diagnostic. It identifies necessary, not sufficient, BA-2 proof targets.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients
from evaluate_ebch128x4_ba_variance_caps import caps_from_variance, load_means


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch128x4_ba2_variance_threshold_q1_m64.json"
B = 512
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000
LOG2 = math.log(2.0)


def parse_float_list(value: str) -> list[float]:
    return [float(item) for item in value.split(",")]


def parse_int_list(value: str) -> list[int]:
    return [int(item) for item in value.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=64)
    parser.add_argument("--grid-min", type=float, default=-11.0)
    parser.add_argument("--grid-max", type=float, default=-7.0)
    parser.add_argument("--grid-step", type=float, default=0.125)
    parser.add_argument("--cutoffs", default="20,21,22,23,24,25,26,27,28")
    parser.add_argument(
        "--variance-bits",
        default="0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4",
    )
    parser.add_argument("--per-shell-failure-bits", type=float, default=51.0)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    q1_best = np.full(B + 1, math.inf)
    q1_witness = np.full(B + 1, math.nan)
    u_values = np.arange(
        args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step
    )
    for index, u_raw in enumerate(u_values):
        u = float(u_raw)
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        regions = uniform_coefficients(
            transfer.log_entries(zero), transfer.log_entries(active), L, 2
        )
        coordinates = uniform_coefficients(regions[0], regions[1], B, B)
        moments = np.logaddexp(coordinates[:, 0, 0], coordinates[:, 0, 1])
        values = np.minimum(0.0, moments + D * surprisal)
        improved = values < q1_best
        q1_best[improved] = values[improved]
        q1_witness[improved] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    mp.mp.dps = 200
    means = load_means(2)
    delta = mp.power(2, -args.per_shell_failure_bits)
    rows = []
    for cutoff in parse_int_list(args.cutoffs):
        for inflation_bits in parse_float_list(args.variance_bits):
            inflation = mp.power(2, inflation_bits)
            caps, failure, methods = caps_from_variance(
                means, delta, inflation, cutoff
            )
            log_caps = np.full(B + 1, -math.inf)
            for weight in range(1, B + 1):
                if int(caps[weight]):
                    log_caps[weight] = math.log(int(caps[weight]))
            terms = math.log(L) + log_caps + q1_best
            finite = np.isfinite(terms)
            q1_log = float(logsumexp(terms[finite]) / LOG2)
            q1_weight = int(np.argmax(terms))
            event_log = float(mp.log(failure, 2))
            combined_log = float(
                np.logaddexp(event_log * LOG2, q1_log * LOG2) / LOG2
            )
            rows.append(
                {
                    "forced_zero_through": cutoff,
                    "variance_inflation_bits": inflation_bits,
                    "variance_inflation": float(inflation),
                    "event_log2_upper": event_log,
                    "event_margin_bits": -event_log,
                    "q1_log2_upper": q1_log,
                    "q1_margin_bits": -q1_log,
                    "q1_dominant_weight": q1_weight,
                    "q1_dominant_log_surprisal": float(q1_witness[q1_weight]),
                    "event_plus_q1_log2_upper": combined_log,
                    "event_plus_q1_margin_bits": -combined_log,
                    "passes_40_bit_necessary_gate": combined_log < -40.0,
                    "cap_method_shell_counts": methods,
                }
            )

    passing = [row for row in rows if row["passes_40_bit_necessary_gate"]]
    best = max(rows, key=lambda row: row["event_plus_q1_margin_bits"])
    largest_factor = (
        max(passing, key=lambda row: row["variance_inflation_bits"])
        if passing
        else None
    )
    result = {
        "schema": "ebch128x4-ba2-variance-threshold-q1-v1",
        "status": "CONDITIONAL_BINARY64_NECESSARY_GATE_DIAGNOSTIC",
        "parameters": {
            "accumulator_stages": 2,
            "memory_bits": args.memory_bits,
            "distance_cutoff": D,
            "per_shell_failure_bits": args.per_shell_failure_bits,
            "grid": [args.grid_min, args.grid_max, args.grid_step],
        },
        "best": best,
        "largest_tested_variance_factor_passing": largest_factor,
        "rows": rows,
        "limitations": [
            "The shell-variance inequality is a hypothesis, not a proved BA-2 fact.",
            "Passing this event-plus-Q1 gate is necessary but not sufficient for the all-Q certificate.",
            "All transfer arithmetic is nearest binary64.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"best": best, "largest_passing": largest_factor}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
