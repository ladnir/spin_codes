#!/usr/bin/env python3
"""Evaluate the proved RM(5,11) shells below weight 128 at Q=1."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


COUNTS = HERE / "rm511_low_weight_exact.json"
OUTPUT = HERE / "rm511_q1_low_randomstepconv_M22_d109_diagnostic.json"
LOG2 = math.log(2.0)


def main() -> int:
    receipt = json.loads(COUNTS.read_text(encoding="utf-8"))
    if receipt.get("status") != "EXACT_INTEGER_FORMULA_CHECK":
        raise AssertionError("low-weight receipt is not accepting")
    counts = {int(weight): int(count) for weight, count in receipt["coefficients"].items()}

    block_bits = 2048
    outer_rows = 1024
    output_bits = 1 << 21
    distance_cutoff = (109 * output_bits + 999) // 1000
    memory_bits = 22
    maximum_weight = 128
    weights = sorted([*counts, maximum_weight])
    best = {weight: (math.inf, math.nan) for weight in weights}

    for index, log_surprisal in enumerate(np.arange(-10.0, -6.999, 0.05)):
        surprisal = math.exp(float(log_surprisal))
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, memory_bits)
        regions = uniform_coefficients(
            transfer.log_entries(zero), transfer.log_entries(active), outer_rows, 1
        )
        coordinates = uniform_coefficients(
            regions[0], regions[1], block_bits, maximum_weight
        )
        moments = np.logaddexp(coordinates[:, 0, 0], coordinates[:, 0, 1])
        for weight in weights:
            multiplicity = counts.get(weight, 1 << 1024)
            value = (
                math.log(outer_rows)
                + math.log(multiplicity)
                + min(0.0, float(moments[weight]) + distance_cutoff * surprisal)
            ) / LOG2
            if value < best[weight][0]:
                best[weight] = (value, float(log_surprisal))
        if (index + 1) % 10 == 0:
            print(f"tilt,{index + 1},61", flush=True)

    rows = [
        {
            "weight": weight,
            "multiplicity": str(counts[weight]),
            "pointwise_log2_upper": best[weight][0],
            "margin_bits": -best[weight][0],
            "log_surprisal": best[weight][1],
        }
        for weight in sorted(counts)
    ]
    aggregate = float(logsumexp([row["pointwise_log2_upper"] * LOG2 for row in rows]) / LOG2)
    payload = {
        "schema": "rm511-exact-low-q1-diagnostic-v1",
        "status": "BINARY64_WITNESS_DIAGNOSTIC",
        "parameters": {
            "outer": "one fixed RM(5,11) [2048,1024,64] constituent repeated 1024 times",
            "output_bits": output_bits,
            "distance_cutoff": distance_cutoff,
            "relative_distance_lower": distance_cutoff / output_bits,
            "inner_memory_bits": memory_bits,
        },
        "claim": {
            "exact_shell_range": "nonzero weights below 128",
            "aggregate_log2_upper_diagnostic": aggregate,
            "aggregate_margin_bits_diagnostic": -aggregate,
            "dominant": max(rows, key=lambda row: row["pointwise_log2_upper"]),
            "weight_128_all_codewords_log2_diagnostic": best[128][0],
        },
        "weight_rows": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only the five exact shells below weight 128 are aggregated.",
            "The weight-128 row deliberately assigns all 2^1024 codewords to one shell and is only a rejected coarse bound.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
