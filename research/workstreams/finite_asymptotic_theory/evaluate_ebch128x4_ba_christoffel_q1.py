#!/usr/bin/env python3
"""Evaluate the Q=1 gate for primal/dual-distance Christoffel caps."""

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
CAPS = WORKSTREAM / "ebch128x4_ba12_christoffel_caps_probe.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba12_christoffel_q1_m64_probe.json"
LENGTH = 512
OUTER_ROWS = 4096
OUTPUT_BITS = 1 << 21
DISTANCE = (109 * OUTPUT_BITS + 999) // 1000
LOG2 = math.log(2.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caps", type=Path, default=CAPS)
    parser.add_argument("--memory-bits", type=int, default=64)
    parser.add_argument("--grid-min", type=float, default=-11.0)
    parser.add_argument("--grid-max", type=float, default=-7.0)
    parser.add_argument("--grid-step", type=float, default=0.125)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    source = json.loads(args.caps.read_text(encoding="utf-8"))
    caps = np.zeros(LENGTH + 1, dtype=object)
    for row in source["caps"]:
        caps[int(row["weight"])] = int(row["cap"])
    log_caps = np.full(LENGTH + 1, -math.inf)
    for weight in range(1, LENGTH + 1):
        if int(caps[weight]):
            log_caps[weight] = math.log(int(caps[weight]))

    best = np.full(LENGTH + 1, math.inf)
    witness = np.full(LENGTH + 1, math.nan)
    tilts = np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step)
    for index, raw_tilt in enumerate(tilts):
        log_surprisal = float(raw_tilt)
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        regions = uniform_coefficients(
            transfer.log_entries(zero), transfer.log_entries(active), OUTER_ROWS, 2
        )
        coordinates = uniform_coefficients(regions[0], regions[1], LENGTH, LENGTH)
        moments = np.logaddexp(coordinates[:, 0, 0], coordinates[:, 0, 1])
        values = np.minimum(0.0, moments + DISTANCE * surprisal)
        improved = values < best
        best[improved] = values[improved]
        witness[improved] = log_surprisal
        print(f"tilt,{index + 1},{len(tilts)},u,{log_surprisal:.6f}", flush=True)

    terms = math.log(OUTER_ROWS) + log_caps + best
    finite = np.isfinite(terms)
    aggregate = float(logsumexp(terms[finite]) / LOG2)
    dominant = int(np.argmax(terms))
    event_log2 = float(source["event"]["union_failure_log2_upper"])
    combined = float(np.logaddexp(event_log2 * LOG2, aggregate * LOG2) / LOG2)
    result = {
        "schema": "ebch128x4-ba-christoffel-q1-probe-v1",
        "status": "CONDITIONAL_BINARY64_EXACT_SHELL_DIAGNOSTIC",
        "parameters": {
            "accumulator_stages": source["construction"]["accumulator_stages"],
            "memory_bits": args.memory_bits,
            "distance_cutoff": DISTANCE,
            "caps": str(args.caps),
        },
        "event_failure_log2_upper": event_log2,
        "q1": {
            "log2_upper": aggregate,
            "margin_bits": -aggregate,
            "dominant_weight": dominant,
            "dominant_log_surprisal": float(witness[dominant]),
        },
        "event_plus_q1": {
            "log2_upper": combined,
            "margin_bits": -combined,
        },
        "limitations": [
            "The Christoffel shell caps are exact conditional on the primal and dual distance event.",
            "The BA expected tail calculations and this transfer use nearest binary64 arithmetic.",
            "Only occupation Q=1 is evaluated.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"q1": result["q1"], "event_plus_q1": result["event_plus_q1"]}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
