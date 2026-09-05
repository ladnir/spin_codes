#!/usr/bin/env python3
"""Sparse transfer with defect, shoulder, and central-body categories."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128x4_ba_variance_proxy_split_sparse import (
    compositions,
    exact_shared_regions,
    merged_band,
)
from evaluate_single_random_constituent_highprob_bands import band_majorant, log_multinomial
from evaluate_single_random_constituent_highprob_renyi import LOG2


WORKSTREAM = Path(__file__).resolve().parent
CAPS = WORKSTREAM / "ebch128x4_ba_variance_cap_requirements.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba_variance_proxy_shoulder_sparse.json"
B = 512
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accumulator-stages", type=int, default=3)
    parser.add_argument("--variance-inflation-bits", type=float, default=12.0)
    parser.add_argument("--shoulder-upper", type=int, default=111)
    parser.add_argument("--minimum-occupation", type=int, default=3)
    parser.add_argument("--maximum-occupation", type=int, default=8)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-3.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if not 80 <= args.shoulder_upper < B // 2:
        parser.error("shoulder upper weight must lie from 80 through 255")

    receipt = json.loads(CAPS.read_text(encoding="utf-8"))
    source = next(
        row
        for row in receipt["rows"]
        if int(row["accumulator_stages"]) == args.accumulator_stages
        and float(row["variance_inflation_bits"]) == args.variance_inflation_bits
    )
    caps = np.asarray([int(value) for value in source["caps"]], dtype=object)
    defect = merged_band(caps, (42, 79), (433, 470))
    shoulder = merged_band(
        caps,
        (80, args.shoulder_upper),
        (B - args.shoulder_upper, 432),
    )
    body = band_majorant(
        caps,
        args.shoulder_upper + 1,
        B - args.shoulder_upper - 1,
        block_bits=B,
    )
    probabilities = [
        float(defect["value_probability"]),
        float(shoulder["value_probability"]),
        float(body["value_probability"]),
    ]
    majorants = [
        float(defect["log_majorant"]),
        float(shoulder["log_majorant"]),
        float(body["log_majorant"]),
    ]
    targets = [
        composition
        for occupation in range(args.minimum_occupation, args.maximum_occupation + 1)
        for composition in compositions(3, occupation)
    ]
    best = {target: math.inf for target in targets}
    witnesses = {target: math.nan for target in targets}
    u_values = np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step)
    for tilt_index, u_raw in enumerate(u_values):
        u = float(u_raw)
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        matrices = [(1.0 - probability) * zero + probability * active for probability in probabilities]
        regions = exact_shared_regions(zero, matrices, args.maximum_occupation)
        for target in targets:
            powered = transfer.log_power(regions[target], B)
            moment = float(np.logaddexp(powered[0, 0], powered[0, 1]))
            outer = log_multinomial((L - sum(target),) + target) + sum(
                count * majorant for count, majorant in zip(target, majorants)
            )
            value = outer + min(0.0, moment + D * surprisal)
            if value < best[target]:
                best[target] = value
                witnesses[target] = u
        print(f"tilt,{tilt_index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    rows = [
        {
            "occupation": sum(target),
            "defect_rows": target[0],
            "shoulder_rows": target[1],
            "body_rows": target[2],
            "pointwise_log2_upper": best[target] / LOG2,
            "margin_bits": -best[target] / LOG2,
            "log_surprisal": witnesses[target],
        }
        for target in targets
    ]
    aggregate = float(logsumexp(list(best.values())))
    result = {
        "schema": "ebch128x4-ba-variance-proxy-shoulder-sparse-v1",
        "status": "CONDITIONAL_BINARY64_EXACT_SHARED_CATEGORY_DIAGNOSTIC",
        "parameters": {
            "accumulator_stages": args.accumulator_stages,
            "variance_inflation_bits": args.variance_inflation_bits,
            "variance_hypothesis": source["variance_hypothesis"],
            "memory_bits": args.memory_bits,
            "occupations": [args.minimum_occupation, args.maximum_occupation],
            "distance_cutoff": D,
            "shoulder_upper": args.shoulder_upper,
        },
        "cap_event_failure_log2_upper": source["event_failure_log2_upper"],
        "bands": {"defect": defect, "shoulder": shoulder, "body": body},
        "claim": {
            "aggregate_log2_upper": aggregate / LOG2,
            "aggregate_margin_bits": -aggregate / LOG2,
            "worst": max(rows, key=lambda row: float(row["pointwise_log2_upper"])),
        },
        "rows": rows,
        "limitations": [
            "The shell-variance hypothesis is not proved.",
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only the stated sparse occupation interval is covered.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
