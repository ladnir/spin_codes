#!/usr/bin/env python3
"""Occupation-one RandomStepConv diagnostic for EBCH128x4 + BA-t.

The spectrum input is the ensemble expectation over the BA interleavers.
Consequently, this calculation is an unconditional occupation-one first
moment over the BA sampler.  It is not a fixed-code spectrum certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


WORKSTREAM = Path(__file__).resolve().parent
SPECTRA = WORKSTREAM / "ebch128x4_ba0_6_B512_expected_spectra.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba0_6_B512_randomstepconv_m22_q1.json"
B = 512
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000
LOG2 = math.log(2.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--minimum-accumulators", type=int, default=0)
    parser.add_argument("--maximum-accumulators", type=int, default=6)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=-3.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    source = json.loads(SPECTRA.read_text(encoding="utf-8"))
    stages = {
        int(stage["accumulators"]): stage
        for stage in source["stages"]
        if args.minimum_accumulators
        <= int(stage["accumulators"])
        <= args.maximum_accumulators
    }
    if set(stages) != set(range(args.minimum_accumulators, args.maximum_accumulators + 1)):
        raise AssertionError("expected-spectrum receipt lacks a requested BA stage")

    spectra = {}
    best = {}
    witnesses = {}
    for stage, payload in stages.items():
        spectrum = np.full(B + 1, -math.inf)
        for row in payload["spectrum"]:
            value = row["log2_expected_multiplicity"]
            if value is not None:
                spectrum[int(row["weight"])] = float(value) * LOG2
        spectra[stage] = spectrum
        best[stage] = np.full(B + 1, math.inf)
        witnesses[stage] = np.full(B + 1, math.nan)

    u_values = np.arange(
        args.grid_min,
        args.grid_max + args.grid_step / 2,
        args.grid_step,
    )
    for index, u_raw in enumerate(u_values):
        u = float(u_raw)
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        regions = uniform_coefficients(
            transfer.log_entries(zero), transfer.log_entries(active), L, 1
        )
        coordinates = uniform_coefficients(regions[0], regions[1], B, B)
        moments = np.logaddexp(coordinates[:, 0, 0], coordinates[:, 0, 1])
        inner = np.minimum(0.0, moments + D * surprisal)
        for stage, spectrum in spectra.items():
            values = math.log(L) + spectrum + inner
            values[0] = math.inf
            improved = values < best[stage]
            best[stage][improved] = values[improved]
            witnesses[stage][improved] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    result_rows = []
    for stage in sorted(stages):
        finite = np.isfinite(spectra[stage][1:])
        weights = np.arange(1, B + 1)[finite]
        values = best[stage][weights]
        aggregate = float(logsumexp(values))
        dominant_index = int(np.argmax(values))
        dominant_weight = int(weights[dominant_index])
        result_rows.append(
            {
                "accumulators": stage,
                "q1_log2_expected_bad_upper": aggregate / LOG2,
                "q1_margin_bits": -aggregate / LOG2,
                "dominant_weight": dominant_weight,
                "dominant_pointwise_log2_upper": float(values[dominant_index] / LOG2),
                "dominant_log_surprisal": float(witnesses[stage][dominant_weight]),
            }
        )

    payload = {
        "schema": "ebch128x4-ba-randomstepconv-q1-v1",
        "status": "BINARY64_ENSEMBLE_EXPECTATION_DIAGNOSTIC",
        "parameters": {
            "base": "direct sum of four extended BCH [128,64,22] codes",
            "outer_block_bits": B,
            "outer_dimension": 256,
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "accumulator_range": [args.minimum_accumulators, args.maximum_accumulators],
        },
        "rows": result_rows,
        "limitations": [
            "Nearest binary64 arithmetic makes the result diagnostic.",
            "Only messages with one active outer row are covered.",
            "The spectrum is averaged over the BA sampler; reuse across two or more active rows requires joint moments or a high-probability fixed-code event.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result_rows, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
