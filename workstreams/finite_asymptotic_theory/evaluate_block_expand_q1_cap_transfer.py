#!/usr/bin/env python3
"""Test Block Expand-3 recentered caps in the exact Q=1 SPIN transfer.

The shell means and resulting caps are diagnostic.  Conditional on those
caps, the occupation-one RandomStepConv-M22 transfer is the same exact-shell
calculation used by the one-shot random-constituent proof.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from certify_single_random_constituent_dense_outward import exact_caps
from evaluate_block_expand_cap_budget import caps_for_variance_factor
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


WORKSTREAM = Path(__file__).resolve().parent
INPUT = WORKSTREAM / "block_expand_accumulate_512_256_d14_t0_5_spectrum.json"
OUTPUT = WORKSTREAM / "block_expand3_q1_cap_transfer_diagnostic.json"
B = 512
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000
MEMORY = 22
LOG2 = math.log(2.0)


def main() -> None:
    mp.mp.dps = 100
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    stage = next(row for row in payload["stages"] if row["accumulator_stages"] == 3)
    log_means = stage["expected_spectrum_log2_by_weight"]
    frozen_caps, _ = exact_caps()
    factors = (1, 2, 4, 16, 512)
    caps_by_factor = {
        factor: caps_for_variance_factor(log_means, frozen_caps, factor)
        for factor in factors
    }
    weights = [weight for weight in range(1, B + 1) if frozen_caps[weight]]
    best = {
        factor: {weight: math.inf for weight in weights} for factor in factors
    }
    witnesses = {
        factor: {weight: math.nan for weight in weights} for factor in factors
    }
    coarse = np.arange(-12.0, 1.0 + 0.25, 0.5)
    fine = np.arange(-10.0, -7.0 + 0.025, 0.05)
    u_values = sorted(set(map(float, coarse)) | set(map(float, fine)))

    for index, u in enumerate(u_values):
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, MEMORY)
        region_coefficients = uniform_coefficients(
            transfer.log_entries(zero),
            transfer.log_entries(active),
            L,
            1,
        )
        coordinate_coefficients = uniform_coefficients(
            region_coefficients[0],
            region_coefficients[1],
            B,
            B,
        )
        moments = np.logaddexp(
            coordinate_coefficients[1:, 0, 0],
            coordinate_coefficients[1:, 0, 1],
        )
        for factor in factors:
            caps = caps_by_factor[factor]
            for weight in weights:
                value = (
                    math.log(L)
                    + math.log(caps[weight])
                    + min(0.0, float(moments[weight - 1]) + D * surprisal)
                )
                if value < best[factor][weight]:
                    best[factor][weight] = value
                    witnesses[factor][weight] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    rows = []
    for factor in factors:
        aggregate = float(logsumexp(list(best[factor].values()))) / LOG2
        dominant_weight = max(weights, key=lambda weight: best[factor][weight])
        rows.append(
            {
                "assumed_variance_factor": factor,
                "q1_log2_upper": aggregate,
                "q1_margin_bits": -aggregate,
                "dominant_weight": dominant_weight,
                "dominant_pointwise_log2_upper": best[factor][dominant_weight] / LOG2,
                "dominant_log_surprisal": witnesses[factor][dominant_weight],
            }
        )

    result = {
        "schema": "block-expand3-q1-cap-transfer-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "construction": "degree-14 regional expander followed by three independently interleaved accumulators",
        "outer_event": "zero shells outside 42..470 and per-shell Cantelli caps inside that interval under the stated variance factor",
        "inner": "RandomStepConv-M22 with the certified Structured SPIN routing",
        "rows": rows,
        "limitations": [
            "The Block Expand means are binary64 and the variance bounds are hypotheses.",
            "Only active-row occupation Q=1 is evaluated.",
            "The calculation is a parameter diagnostic, not outward rounded.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
