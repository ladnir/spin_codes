#!/usr/bin/env python3
"""Merge disjoint Chernoff-tilt sweeps by occupation."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


LOG2 = math.log(2.0)


def logsumexp(values: np.ndarray) -> float:
    maximum = float(np.max(values))
    if maximum == -math.inf:
        return -math.inf
    return maximum + math.log(float(np.sum(np.exp(values - maximum))))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in args.inputs]
    best: dict[int, dict[str, object]] = {}
    for source, payload in zip(args.inputs, payloads):
        for raw in payload["occupation_rows"]:
            row = dict(raw)
            occupation = int(row["active_regular_outer_blocks"])
            row["source_receipt"] = str(source)
            current = best.get(occupation)
            if current is None or float(row["pointwise_log2_upper"]) < float(
                current["pointwise_log2_upper"]
            ):
                best[occupation] = row

    rows = [best[key] for key in sorted(best)]
    pointwise = np.asarray(
        [float(row["pointwise_log2_upper"]) * LOG2 for row in rows]
    )
    aggregate = float(logsumexp(pointwise)) / LOG2
    dominant = sorted(
        rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]
    result = {
        "schema": "riffle-splitstate-merged-tilt-receipts-v1",
        "candidate": payloads[0]["candidate"],
        "input_receipts": [str(path) for path in args.inputs],
        "parameters": payloads[0]["parameters"],
        "covered_minimum_occupation": int(rows[0]["active_regular_outer_blocks"]),
        "covered_maximum_occupation": int(rows[-1]["active_regular_outer_blocks"]),
        "partial_log2_upper": aggregate,
        "partial_margin_bits": -aggregate,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "scope": (
            "Pointwise minimum across the supplied floating-point Chernoff "
            "grids. This does not change the proof scope of the input receipts."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"partial_margin,{result['partial_margin_bits']:.6f}")
    for row in dominant[:5]:
        print(
            f"dominant,{row['active_regular_outer_blocks']},"
            f"tilt,{row['best_log_surprisal']:.6f},"
            f"pointwise_log2,{row['pointwise_log2_upper']:.6f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
