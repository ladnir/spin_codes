#!/usr/bin/env python3
"""Column-Hölder diagnostic for any repeated full-support outer code.

Fix a binary [B,K] linear code whose every coordinate is nonzero.  For a
uniform nonzero codeword, each coordinate equals one with the same probability
2^(K-1)/(2^K-1).  Independent row-coordinate and region permutations give an
exact one-region law.  Hölder across the B dependent regions then removes all
remaining dependence on the code's weight spectrum.

This nearest-binary64 calculation tests whether the spectrum-free argument can
cover dense occupations of the EBCH128x4--BA construction.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "full_support_B512_randomstepconv_m22_column_holder.json"
B = 512
K = 256
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000
LOG2 = math.log(2.0)
NONZERO = (1 << K) - 1
ONE_PROBABILITY = (1 << (K - 1)) / NONZERO


def log_choose(total: int, selected: np.ndarray | int) -> np.ndarray | float:
    return gammaln(total + 1.0) - gammaln(np.asarray(selected) + 1.0) - gammaln(total - np.asarray(selected) + 1.0)


def log_weighted_norms(regions: np.ndarray, log_ratio: float) -> np.ndarray:
    row_zero = np.logaddexp(regions[:, 0, 0], regions[:, 0, 1] + log_ratio)
    row_live = np.logaddexp(regions[:, 1, 0] - log_ratio, regions[:, 1, 1])
    return np.maximum(row_zero, row_live)


def evaluate_occupation(
    regions: np.ndarray,
    occupation: int,
    surprisal: float,
) -> dict[str, object]:
    weights = np.arange(occupation + 1, dtype=np.float64)
    log_binomial = (
        log_choose(occupation, weights)
        + weights * math.log(ONE_PROBABILITY)
        + (occupation - weights) * math.log1p(-ONE_PROBABILITY)
    )

    def objective(log_ratio: float) -> float:
        log_norms = log_weighted_norms(regions[: occupation + 1], log_ratio)
        region_holder = float(logsumexp(log_binomial + B * log_norms))
        terminal = max(0.0, -log_ratio)
        outer = float(log_choose(L, occupation)) + occupation * math.log(NONZERO)
        return outer + terminal + region_holder + D * surprisal

    result = minimize_scalar(
        objective,
        bounds=(-20.0, 20.0),
        method="bounded",
        options={"xatol": 1e-10},
    )
    raw_value = float(result.fun)
    value = min(0.0, raw_value)
    return {
        "occupation": occupation,
        "pointwise_log2_upper": value / LOG2,
        "margin_bits": -value / LOG2,
        "raw_log2_upper": raw_value / LOG2,
        "log_norm_ratio": float(result.x),
        "optimizer_success": bool(result.success),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--occupations", type=int, nargs="+", default=[128, 160, 256, 512, 1024, 2048, 3072, 4096]
    )
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument(
        "--log-surprisals", type=float, nargs="+", default=[0.4, 0.55, 0.7, 0.85, 1.0]
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    maximum_q = max(args.occupations)
    best = {occupation: None for occupation in args.occupations}
    for index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        regions = uniform_coefficients(
            transfer.log_entries(zero),
            transfer.log_entries(active),
            L,
            maximum_q,
        )
        for occupation in args.occupations:
            row = evaluate_occupation(regions, occupation, surprisal)
            row["log_surprisal"] = log_surprisal
            current = best[occupation]
            if current is None or float(row["pointwise_log2_upper"]) < float(current["pointwise_log2_upper"]):
                best[occupation] = row
        print(f"tilt,{index + 1},{len(args.log_surprisals)},u,{log_surprisal:.6f}", flush=True)

    rows = [best[occupation] for occupation in args.occupations]
    result = {
        "schema": "full-support-b512-randomstepconv-column-holder-v1",
        "status": "BINARY64_SPECTRUM_FREE_DIAGNOSTIC",
        "parameters": {
            "outer_class": "any full-support binary [512,256] linear code repeated in every row",
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "one_coordinate_probability": ONE_PROBABILITY,
            "occupations": args.occupations,
        },
        "rows": rows,
        "limitations": [
            "Only the listed occupations and tilts are evaluated.",
            "Nearest binary64 arithmetic is not an outward certificate.",
            "The matrix-norm and column-Hölder reduction requires a written proof audit.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(rows, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
