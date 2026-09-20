#!/usr/bin/env python3
"""Measure whether distance-only shell packing can cover the RM(5,11) gap."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


OUTPUT = HERE / "rm511_q1_packing_gap_diagnostic.json"
BLOCK_BITS = 2048
OUTER_DIMENSION = 1024
OUTER_ROWS = 1024
MINIMUM_DISTANCE = 64
OUTPUT_BITS = 1 << 21
DISTANCE_CUTOFF = (109 * OUTPUT_BITS + 999) // 1000
MEMORY_BITS = 22
WEIGHTS = (128, 160, 192, 224, 256, 320, 384)


def packing_upper(weight: int) -> int:
    """Johnson-space radius-15 sphere-packing upper bound."""

    johnson_distance = (MINIMUM_DISTANCE + 1) // 2
    radius = (johnson_distance - 1) // 2
    ball = sum(
        math.comb(weight, index) * math.comb(BLOCK_BITS - weight, index)
        for index in range(radius + 1)
    )
    return min((1 << OUTER_DIMENSION) - 1, math.comb(BLOCK_BITS, weight) // ball)


def main() -> int:
    best_inner = {weight: (math.inf, math.nan) for weight in WEIGHTS}
    tilts = np.arange(-10.0, -3.999, 0.05)
    for index, log_surprisal in enumerate(tilts):
        surprisal = math.exp(float(log_surprisal))
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, MEMORY_BITS)
        regions = uniform_coefficients(
            transfer.log_entries(zero), transfer.log_entries(active), OUTER_ROWS, 1
        )
        coordinates = uniform_coefficients(
            regions[0], regions[1], BLOCK_BITS, max(WEIGHTS)
        )
        moments = np.logaddexp(coordinates[:, 0, 0], coordinates[:, 0, 1])
        for weight in WEIGHTS:
            inner = min(0.0, float(moments[weight]) + DISTANCE_CUTOFF * surprisal) / math.log(2.0)
            if inner < best_inner[weight][0]:
                best_inner[weight] = (inner, float(log_surprisal))
        if (index + 1) % 10 == 0:
            print(f"tilt,{index + 1},{len(tilts)}", flush=True)

    rows = []
    for weight in WEIGHTS:
        inner, tilt = best_inner[weight]
        cap = packing_upper(weight)
        cap_log2 = math.log2(cap)
        packed = math.log2(OUTER_ROWS) + cap_log2 + inner
        all_mass = math.log2(OUTER_ROWS) + OUTER_DIMENSION + inner
        required_cap_log2_for_40 = -40.0 - math.log2(OUTER_ROWS) - inner
        rows.append(
            {
                "weight": weight,
                "inner_log2_upper_diagnostic": inner,
                "log_surprisal": tilt,
                "johnson_packing_cap_log2": cap_log2,
                "packing_pointwise_log2_upper_diagnostic": packed,
                "all_mass_pointwise_log2_upper_diagnostic": all_mass,
                "required_shell_log2_cap_for_40_bits": required_cap_log2_for_40,
                "packing_cap_excess_bits": cap_log2 - required_cap_log2_for_40,
            }
        )

    payload = {
        "schema": "rm511-q1-packing-gap-diagnostic-v1",
        "status": "BINARY64_WITNESS_DIAGNOSTIC",
        "parameters": {
            "outer": "RM(5,11) [2048,1024,64]",
            "outer_rows": OUTER_ROWS,
            "output_bits": OUTPUT_BITS,
            "distance_cutoff": DISTANCE_CUTOFF,
            "inner_memory_bits": MEMORY_BITS,
            "johnson_packing_radius": 15,
        },
        "rows": rows,
        "conclusion": (
            "Distance-only constant-weight packing is insufficient wherever "
            "packing_cap_excess_bits is positive. The required cap column is "
            "the per-shell threshold for a 40-bit pointwise contribution."
        ),
        "limitations": [
            "All logarithmic transfer values use nearest binary64 arithmetic.",
            "This diagnostic tests selected weights; it is not a spectrum envelope or certificate.",
            "The per-shell 40-bit thresholds do not allocate a union budget across shells.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(rows, indent=2))
    print(f"output={OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
