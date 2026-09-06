#!/usr/bin/env python3
"""Exact-shell Q=1 diagnostic for repeated EBCH128 + RandomStepConv."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as base


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch128_randomstepconv_g1_s30_q1_exact_d11_diagnostic.json"


def uniform_coefficients(
    zero: np.ndarray,
    candidate: np.ndarray,
    positions: int,
    maximum_degree: int,
) -> np.ndarray:
    current = np.full((maximum_degree + 1, 2, 2), -math.inf)
    current[0] = base.log_identity()
    for completed in range(positions):
        maximum = min(completed + 1, maximum_degree)
        old_maximum = min(completed, maximum_degree)
        updated = np.full_like(current, -math.inf)
        zero_products = base.log_matmul(current[: old_maximum + 1], zero)
        degrees = np.arange(old_maximum + 1, dtype=np.float64)
        updated[: old_maximum + 1] = zero_products + np.log(
            ((completed + 1.0) - degrees) / (completed + 1.0)
        )[:, None, None]
        candidate_products = base.log_matmul(current[:maximum], candidate)
        selected = np.arange(1, maximum + 1, dtype=np.float64)
        candidate_terms = candidate_products + np.log(
            selected / (completed + 1.0)
        )[:, None, None]
        updated[1 : maximum + 1] = np.logaddexp(
            updated[1 : maximum + 1], candidate_terms
        )
        current = updated
    return current


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    spectrum = base.load_spectrum(args.spectrum)
    weights = sorted(weight for weight in spectrum if weight)
    best = {weight: math.inf for weight in weights}
    witnesses = {weight: math.nan for weight in weights}
    u_values = sorted(
        set(
            float(value)
            for value in np.arange(-12.0, 1.0 + 0.25, 0.5)
        )
        | set(
            float(value)
            for value in np.arange(0.65, 0.80 + 0.0125, 0.025)
        )
        | set(
            float(value)
            for value in np.arange(-9.50, -8.00 + 0.025, 0.05)
        )
    )
    for index, u in enumerate(u_values):
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = base.step_matrices(z, args.memory_bits)
        bit_zero = base.log_entries(zero)
        bit_one = base.log_entries(active)
        region_coefficients = uniform_coefficients(
            bit_zero, bit_one, base.L, 1
        )
        zero_region = region_coefficients[0]
        one_region = region_coefficients[1]
        coordinate_coefficients = uniform_coefficients(
            zero_region, one_region, base.B, base.B
        )
        moments = np.logaddexp(
            coordinate_coefficients[:, 0, 0],
            coordinate_coefficients[:, 0, 1],
        )
        correction = base.D * surprisal
        for weight in weights:
            value = (
                math.log(base.ACTIVE_ROWS)
                + math.log(spectrum[weight])
                + min(0.0, float(moments[weight]) + correction)
            )
            if value < best[weight]:
                best[weight] = value
                witnesses[weight] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    aggregate = float(logsumexp(list(best.values())))
    rows = [
        {
            "outer_weight": weight,
            "multiplicity": spectrum[weight],
            "pointwise_log2_upper": best[weight] / base.LOG2,
            "log_surprisal": witnesses[weight],
        }
        for weight in weights
    ]
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    return {
        "schema": "ebch128-randomstepconv-g1-q1-exact-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "claim": {
            "log2_expected_bad_upper": aggregate / base.LOG2,
            "margin_bits": -aggregate / base.LOG2,
            "closes_40_bits": aggregate < -40.0 * base.LOG2,
            "dominant": dominant,
        },
        "parameters": {
            "outer_code": "extended BCH [128,64,22]",
            "active_outer_rows": base.ACTIVE_ROWS,
            "outer_rows": base.L,
            "output_bits": base.N,
            "distance_cutoff": base.D,
            "memory_bits": args.memory_bits,
            "occupation": 1,
        },
        "probability_space": (
            "the EBCH constituent is fixed and repeated; one active row is "
            "chosen; its coordinate permutation and every region permutation "
            "are uniform; RandomStepConv maps are independent and shared"
        ),
        "shell_rows": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "This receipt covers only occupation one.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=base.SPECTRUM)
    parser.add_argument("--outer-rows", type=int, default=base.L)
    parser.add_argument("--memory-bits", type=int, default=base.MEMORY)
    parser.add_argument("--distance-numerator", type=int, default=11)
    parser.add_argument("--distance-denominator", type=int, default=100)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    base.configure_outer_rows(
        args.outer_rows,
        args.distance_numerator,
        args.distance_denominator,
    )
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
