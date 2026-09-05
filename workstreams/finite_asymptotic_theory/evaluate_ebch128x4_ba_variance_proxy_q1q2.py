#!/usr/bin/env python3
"""Exact Q=1 and Q=2 transfer under hypothetical BA shell caps."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients
from evaluate_single_random_constituent_highprob_renyi import LOG2
from evaluate_single_random_constituent_q2 import pair_support_coefficients


WORKSTREAM = Path(__file__).resolve().parent
CAPS = WORKSTREAM / "ebch128x4_ba_variance_cap_requirements.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba_variance_proxy_q1q2.json"
B = 512
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accumulator-stages", type=int, default=3)
    parser.add_argument("--variance-inflation-bits", type=float, default=4.0)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-7.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--q1-only", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--caps", type=Path, default=CAPS)
    args = parser.parse_args()

    receipt = json.loads(args.caps.read_text(encoding="utf-8"))
    source = next(
        row
        for row in receipt["rows"]
        if int(row["accumulator_stages"]) == args.accumulator_stages
        and float(row["variance_inflation_bits"]) == args.variance_inflation_bits
    )
    caps = np.asarray([int(value) for value in source["caps"]], dtype=object)
    log_caps = np.full(B + 1, -math.inf)
    for weight in range(1, B + 1):
        if int(caps[weight]):
            log_caps[weight] = math.log(int(caps[weight]))

    q1_best = np.full(B + 1, math.inf)
    q2_best = None if args.q1_only else np.full((B + 1, B + 1), math.inf)
    q1_witness = np.full(B + 1, math.nan)
    q2_witness = None if args.q1_only else np.full_like(q2_best, math.nan)
    u_values = np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step)
    for index, u_raw in enumerate(u_values):
        u = float(u_raw)
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        regions = uniform_coefficients(
            transfer.log_entries(zero), transfer.log_entries(active), L, 2
        )

        q1_coordinates = uniform_coefficients(regions[0], regions[1], B, B)
        q1_moments = np.logaddexp(q1_coordinates[:, 0, 0], q1_coordinates[:, 0, 1])
        q1_values = np.minimum(0.0, q1_moments + D * surprisal)
        q1_improved = q1_values < q1_best
        q1_best[q1_improved] = q1_values[q1_improved]
        q1_witness[q1_improved] = u

        if not args.q1_only:
            assert q2_best is not None and q2_witness is not None
            q2_coordinates = pair_support_coefficients(regions, B)
            q2_moments = np.logaddexp(q2_coordinates[..., 0, 0], q2_coordinates[..., 0, 1])
            q2_values = np.minimum(0.0, q2_moments + D * surprisal)
            q2_improved = q2_values < q2_best
            q2_best[q2_improved] = q2_values[q2_improved]
            q2_witness[q2_improved] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    q1_terms = math.log(L) + log_caps + q1_best
    q1_finite = np.isfinite(q1_terms)
    q1_log = float(logsumexp(q1_terms[q1_finite]))
    q1_weight = int(np.argmax(q1_terms))

    q2 = None
    if not args.q1_only:
        assert q2_best is not None and q2_witness is not None
        row_pairs = math.log(math.comb(L, 2))
        q2_terms = row_pairs + log_caps[:, None] + log_caps[None, :] + q2_best
        q2_finite = np.isfinite(q2_terms)
        q2_log = float(logsumexp(q2_terms[q2_finite]))
        q2_flat = int(np.argmax(q2_terms))
        q2_first, q2_second = divmod(q2_flat, B + 1)
        q2 = {
            "log2_upper": q2_log / LOG2,
            "margin_bits": -q2_log / LOG2,
            "dominant_weights": [q2_first, q2_second],
            "dominant_log_surprisal": float(q2_witness[q2_first, q2_second]),
        }

    result = {
        "schema": "ebch128x4-ba-variance-proxy-q1q2-v1",
        "status": "CONDITIONAL_BINARY64_EXACT_SHELL_DIAGNOSTIC",
        "parameters": {
            "accumulator_stages": args.accumulator_stages,
            "variance_inflation_bits": args.variance_inflation_bits,
            "variance_hypothesis": source["variance_hypothesis"],
            "memory_bits": args.memory_bits,
            "distance_cutoff": D,
        },
        "cap_event_failure_log2_upper": source["event_failure_log2_upper"],
        "q1": {
            "log2_upper": q1_log / LOG2,
            "margin_bits": -q1_log / LOG2,
            "dominant_weight": q1_weight,
            "dominant_log_surprisal": float(q1_witness[q1_weight]),
        },
        "q2": q2,
        "limitations": [
            "The shell-variance hypothesis is not proved.",
            "Nearest binary64 arithmetic is not an outward certificate.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"q1": result["q1"], "q2": result["q2"]}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
