#!/usr/bin/env python3
"""Sparse SPIN transfer under a hypothetical BA shell-variance lemma.

The input cap receipt assumes Var(A_w)<=F E[A_w].  This program treats those
caps as a fixed-code spectrum envelope and runs the exact shared two-category
region recurrence used by the closed random-constituent proof.

The result is conditional on the stated variance lemma and uses nearest
binary64 arithmetic.  It quantifies the lemma strength needed for transfer.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_single_random_constituent_highprob_bands import band_majorant, log_multinomial
from evaluate_single_random_constituent_highprob_renyi import LOG2
from evaluate_single_random_constituent_shared_two_band import (
    exact_two_category_regions,
)


WORKSTREAM = Path(__file__).resolve().parent
CAPS = WORKSTREAM / "ebch128x4_ba_variance_cap_requirements.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba_variance_proxy_sparse.json"
B = 512
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accumulator-stages", type=int, default=3)
    parser.add_argument("--variance-inflation-bits", type=float, default=12.0)
    parser.add_argument("--minimum-occupation", type=int, default=3)
    parser.add_argument("--maximum-occupation", type=int, default=16)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--low-upper", type=int)
    parser.add_argument("--grid-min", type=float, default=-8.0)
    parser.add_argument("--grid-max", type=float, default=-3.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
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
    if args.low_upper is None:
        low = source["low"]
        central = source["central"]
        high = source["high"]
    else:
        caps = np.asarray([int(value) for value in source["caps"]], dtype=object)
        support = [weight for weight in range(1, B + 1) if int(caps[weight])]
        if not min(support) <= args.low_upper < B // 2:
            parser.error("low upper boundary does not split the spectrum support")
        central_upper = B - args.low_upper - 1
        low = band_majorant(caps, min(support), args.low_upper, block_bits=B)
        central = band_majorant(caps, args.low_upper + 1, central_upper, block_bits=B)
        high = band_majorant(caps, central_upper + 1, max(support), block_bits=B)
    defect_probability = min(
        float(low["value_probability"]),
        1.0 - float(high["value_probability"]),
    )
    defect_log_majorant = math.log(2.0) + max(
        float(low["log_majorant"]), float(high["log_majorant"])
    )
    central_probability = float(central["value_probability"])
    central_log_majorant = float(central["log_majorant"])

    compositions = [
        (defects, occupation - defects)
        for occupation in range(args.minimum_occupation, args.maximum_occupation + 1)
        for defects in range(occupation + 1)
    ]
    best = {composition: math.inf for composition in compositions}
    witnesses = {composition: math.nan for composition in compositions}
    u_values = np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step)
    for index, u_raw in enumerate(u_values):
        u = float(u_raw)
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        defect = (1.0 - defect_probability) * zero + defect_probability * active
        central_matrix = (1.0 - central_probability) * zero + central_probability * active
        regions = exact_two_category_regions(
            zero, defect, central_matrix, args.maximum_occupation, L
        )
        for composition in compositions:
            powered = transfer.log_power(regions[composition], B)
            moment = float(np.logaddexp(powered[0, 0], powered[0, 1]))
            defects, centrals = composition
            outer = (
                log_multinomial((L - defects - centrals, defects, centrals))
                + defects * defect_log_majorant
                + centrals * central_log_majorant
            )
            value = outer + min(0.0, moment + D * surprisal)
            if value < best[composition]:
                best[composition] = value
                witnesses[composition] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    rows = [
        {
            "occupation": sum(composition),
            "defect_rows": composition[0],
            "central_rows": composition[1],
            "pointwise_log2_upper": best[composition] / LOG2,
            "margin_bits": -best[composition] / LOG2,
            "log_surprisal": witnesses[composition],
        }
        for composition in compositions
    ]
    aggregate = float(logsumexp(list(best.values())))
    result = {
        "schema": "ebch128x4-ba-variance-proxy-sparse-v1",
        "status": "CONDITIONAL_BINARY64_EXACT_SHARED_CATEGORY_DIAGNOSTIC",
        "parameters": {
            "accumulator_stages": args.accumulator_stages,
            "variance_inflation_bits": args.variance_inflation_bits,
            "variance_hypothesis": source["variance_hypothesis"],
            "memory_bits": args.memory_bits,
            "occupations": [args.minimum_occupation, args.maximum_occupation],
            "distance_cutoff": D,
            "low_upper": int(low["upper_weight"]),
        },
        "cap_event_failure_log2_upper": source["event_failure_log2_upper"],
        "bands": {"low": low, "central": central, "high": high},
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
