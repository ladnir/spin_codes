#!/usr/bin/env python3
"""Sample dense SPIN compositions under hypothetical BA spectrum caps."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from evaluate_single_random_constituent_highprob_bands import band_majorant
from evaluate_single_random_constituent_shared_two_band_dense import (
    LOG2,
    optimize_counts,
)


WORKSTREAM = Path(__file__).resolve().parent
CAPS = WORKSTREAM / "ebch128x4_ba_variance_cap_requirements.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba_variance_proxy_dense.json"
B = 512
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000


def defect_samples(occupation: int, grid: int) -> list[int]:
    return sorted(
        {
            min(occupation, max(0, int(round(occupation * index / grid))))
            for index in range(grid + 1)
        }
        | {0, min(1, occupation), min(2, occupation), max(0, occupation - 1), occupation}
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accumulator-stages", type=int, default=3)
    parser.add_argument("--variance-inflation-bits", type=float, default=12.0)
    parser.add_argument(
        "--occupations", type=int, nargs="+", default=[17, 32, 64, 128, 160, 256, 512, 1024, 2048, 4096]
    )
    parser.add_argument("--composition-grid", type=int, default=8)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--low-upper", type=int)
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
    defect_majorant = math.log(2.0) + max(
        float(low["log_majorant"]), float(high["log_majorant"])
    )
    central_probability = float(central["value_probability"])
    central_majorant = float(central["log_majorant"])
    values = np.asarray([0.0, defect_probability, central_probability])
    majorants = np.asarray([0.0, defect_majorant, central_majorant])

    rows = []
    for occupation in args.occupations:
        for defects in defect_samples(occupation, args.composition_grid):
            counts = (L - occupation, defects, occupation - defects)
            row = optimize_counts(
                counts,
                values=values,
                log_majorants=majorants,
                block_bits=B,
                distance_cutoff=D,
                memory_bits=args.memory_bits,
                output_bits=N,
                temperatures=(0.6, 1.0),
                scales=(1.4, 3.0),
                maximum_iterations=700,
            )
            rows.append(row)
            print(
                f"q,{occupation},defects,{defects},margin,{row['margin_bits']:.6f}",
                flush=True,
            )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    result = {
        "schema": "ebch128x4-ba-variance-proxy-dense-v1",
        "status": "CONDITIONAL_BINARY64_SAMPLED_COMPOSITION_DIAGNOSTIC",
        "parameters": {
            "accumulator_stages": args.accumulator_stages,
            "variance_inflation_bits": args.variance_inflation_bits,
            "variance_hypothesis": source["variance_hypothesis"],
            "memory_bits": args.memory_bits,
            "occupations": args.occupations,
            "composition_grid": args.composition_grid,
            "distance_cutoff": D,
            "low_upper": int(low["upper_weight"]),
        },
        "claim": {"worst_sampled": worst},
        "rows": rows,
        "limitations": [
            "The shell-variance hypothesis is not proved.",
            "Occupations and compositions are sampled, not a complete cover.",
            "Nearest binary64 arithmetic is not an outward certificate.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
