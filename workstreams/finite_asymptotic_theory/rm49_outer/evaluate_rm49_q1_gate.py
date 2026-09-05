#!/usr/bin/env python3
"""Binary64 witness search for the exact RM(4,9) occupation-one gate.

This script uses the same uniform-routing and RandomStepConv transfer as the
closed sparse-EA theorem, but replaces spectrum caps by the exact RM(4,9)
weight multiplicities.  The output selects tilts; it is not an outward
certificate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
sys.path.insert(0, str(WORKSTREAM))

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


LOG2 = math.log(2.0)
DEFAULT_SPECTRUM = REPOSITORY / "scripts" / "rm512_256_spectrum.csv"
DEFAULT_OUTPUT = HERE / "rm49_q1_randomstepconv_M22_d109_diagnostic.json"


def read_spectrum(path: Path) -> dict[int, int]:
    with path.open(newline="", encoding="utf-8") as source:
        spectrum = {
            int(row["weight"]): int(row["count"])
            for row in csv.DictReader(source)
        }
    if sum(spectrum.values()) != 1 << 256:
        raise AssertionError("RM spectrum mass is not 2^256")
    if spectrum.get(0) != 1 or spectrum.get(512) != 1:
        raise AssertionError("RM spectrum endpoints changed")
    if min(weight for weight, count in spectrum.items() if weight and count) != 32:
        raise AssertionError("RM minimum distance changed")
    return spectrum


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--fine-min", type=float, default=-10.0)
    parser.add_argument("--fine-max", type=float, default=-5.0)
    parser.add_argument("--fine-step", type=float, default=0.05)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    block_bits = 512
    outer_rows = 4096
    output_bits = 1 << 21
    distance_cutoff = (109 * output_bits + 999) // 1000
    spectrum = read_spectrum(args.spectrum)
    weights = sorted(weight for weight, count in spectrum.items() if weight and count)
    best = {weight: math.inf for weight in weights}
    witnesses = {weight: math.nan for weight in weights}
    coarse = np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step)
    fine = np.arange(args.fine_min, args.fine_max + args.fine_step / 2, args.fine_step)
    tilts = sorted(set(map(float, coarse)) | set(map(float, fine)))

    for index, log_surprisal in enumerate(tilts):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        regions = uniform_coefficients(
            transfer.log_entries(zero), transfer.log_entries(active), outer_rows, 1
        )
        coordinates = uniform_coefficients(regions[0], regions[1], block_bits, block_bits)
        moments = np.logaddexp(coordinates[1:, 0, 0], coordinates[1:, 0, 1])
        for weight in weights:
            log_bound = (
                math.log(outer_rows)
                + math.log(spectrum[weight])
                + min(0.0, float(moments[weight - 1]) + distance_cutoff * surprisal)
            )
            if log_bound < best[weight]:
                best[weight] = log_bound
                witnesses[weight] = log_surprisal
        print(f"tilt,{index + 1},{len(tilts)},u,{log_surprisal:.6f}", flush=True)

    rows = [
        {
            "weight": weight,
            "multiplicity": str(spectrum[weight]),
            "pointwise_log2_upper": best[weight] / LOG2,
            "log_surprisal": witnesses[weight],
        }
        for weight in weights
    ]
    aggregate = float(logsumexp([best[weight] for weight in weights]) / LOG2)
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    payload = {
        "schema": "rm49-exact-q1-randomstepconv-diagnostic-v1",
        "status": "BINARY64_WITNESS_DIAGNOSTIC",
        "parameters": {
            "outer": "one fixed RM(4,9) [512,256,32] constituent repeated 4096 times",
            "outer_rows": outer_rows,
            "output_bits": output_bits,
            "distance_cutoff": distance_cutoff,
            "relative_distance_lower": distance_cutoff / output_bits,
            "memory_bits": args.memory_bits,
            "routing": "independent uniform row-coordinate and transposed-region permutations",
        },
        "spectrum": {
            "file": str(args.spectrum.relative_to(REPOSITORY)).replace("\\", "/"),
            "sha256": hashlib.sha256(args.spectrum.read_bytes()).hexdigest(),
            "total_mass": str(1 << 256),
            "nonzero_supported_weights": len(weights),
        },
        "claim": {
            "q1_log2_upper_diagnostic": aggregate,
            "q1_margin_bits_diagnostic": -aggregate,
            "dominant": dominant,
        },
        "weight_rows": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "The tilt grid supplies witnesses and does not assert optimality.",
            "Only occupation one is covered.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
