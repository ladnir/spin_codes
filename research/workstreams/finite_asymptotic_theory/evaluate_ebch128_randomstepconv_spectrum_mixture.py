#!/usr/bin/env python3
"""Diagnostic positive-defect mixture for full EBCH128 occupation."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as base
from evaluate_ebch128_randomstepconv_dense_complement import (
    log_matrix_sum,
    upper_pivot_matrix,
)
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch128_randomstepconv_s30_pow2_spectrum_mixture.json"


def positive_defect_mass() -> Fraction:
    spectrum = base.load_spectrum(base.SPECTRUM)
    baseline = 1 << base.K
    denominator = 1 << (base.B - 1)
    defect = Fraction(0)
    for weight in range(0, base.B + 1, 2):
        code_mass = Fraction(spectrum.get(weight, 0)) if weight else Fraction(0)
        reference_mass = Fraction(baseline * math.comb(base.B, weight), denominator)
        if code_mass > reference_mass:
            defect += code_mass - reference_mass
    return defect


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    base.configure_outer_rows(base.ACTIVE_ROWS)
    baseline_mass = 1 << base.K
    defect_mass = positive_defect_mass()
    surprisal = math.exp(args.log_surprisal)
    z = math.exp(-surprisal)
    zero, active = base.step_matrices(z, args.memory_bits)
    candidate = 0.5 * (zero + active)
    region_by_zeros = uniform_coefficients(
        base.log_entries(candidate),
        base.log_entries(zero),
        base.L,
        args.maximum_region_zeros,
    )
    candidate_row_norm = float(np.max(np.sum(candidate, axis=1)))
    log_baseline = math.log(baseline_mass)
    log_defect = math.log(defect_mass.numerator) - math.log(defect_mass.denominator)

    rows: list[dict[str, object]] = []
    terms: list[float] = []
    for original_inactive in range(args.maximum_original_inactive_rows + 1):
        original_active = base.L - original_inactive
        maximum_defects = min(args.maximum_defect_rows, original_active)
        for defect_rows in range(maximum_defects + 1):
            baseline_rows = original_active - defect_rows
            fixed_inactive = original_inactive + defect_rows
            best = math.inf
            best_log_r = math.nan
            for log_r in np.arange(
                args.log_r_min,
                args.log_r_max + args.log_r_step / 2.0,
                args.log_r_step,
            ):
                pivot_matrix = upper_pivot_matrix(
                    region_by_zeros,
                    inactive_rows=fixed_inactive,
                    active_rows=baseline_rows,
                    log_r=float(log_r),
                    candidate_row_norm=candidate_row_norm,
                )
                all_regions = base.log_power(pivot_matrix, base.B)
                conditional = (
                    math.lgamma(baseline_rows + 1)
                    - baseline_rows * math.log(base.B)
                    - baseline_rows * float(log_r)
                    + log_matrix_sum(all_regions)
                )
                inner = min(0.0, conditional + base.D * surprisal)
                outer = (
                    math.lgamma(base.L + 1)
                    - math.lgamma(original_active + 1)
                    - math.lgamma(original_inactive + 1)
                    + math.lgamma(original_active + 1)
                    - math.lgamma(defect_rows + 1)
                    - math.lgamma(baseline_rows + 1)
                    + baseline_rows * log_baseline
                    + defect_rows * log_defect
                )
                pointwise = (outer + inner) / base.LOG2
                if pointwise < best:
                    best = pointwise
                    best_log_r = float(log_r)
            terms.append(best * base.LOG2)
            rows.append(
                {
                    "original_active_rows": original_active,
                    "original_inactive_rows": original_inactive,
                    "defect_rows": defect_rows,
                    "baseline_rows": baseline_rows,
                    "pointwise_log2_upper": best,
                    "log_pivot_fugacity": best_log_r,
                }
            )

    aggregate = float(logsumexp(terms)) / base.LOG2
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    return {
        "schema": "ebch128-randomstepconv-spectrum-mixture-diagnostic-v1",
        "status": "PARTIAL_BINARY64_DIAGNOSTIC",
        "claim": {
            "covered_defect_rows": [0, args.maximum_defect_rows],
            "covered_original_inactive_rows": [
                0,
                args.maximum_original_inactive_rows,
            ],
            "partial_sum_log2_upper": aggregate,
            "dominant_covered_term": dominant,
        },
        "parameters": {
            "outer_rows": base.L,
            "output_bits": base.N,
            "distance_cutoff": base.D,
            "memory_bits": args.memory_bits,
            "log_surprisal": args.log_surprisal,
            "maximum_exact_region_zeros": args.maximum_region_zeros,
        },
        "measure_decomposition": {
            "uniform_even_baseline_mass": baseline_mass,
            "positive_defect_mass_numerator": defect_mass.numerator,
            "positive_defect_mass_denominator": defect_mass.denominator,
            "defect_to_baseline_ratio": float(defect_mass / baseline_mass),
        },
        "occupation_rows": rows,
        "limitations": [
            "The sum covers only the displayed defect-row range.",
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only the displayed original-occupation and defect-row ranges are included.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=30)
    parser.add_argument("--maximum-original-inactive-rows", type=int, default=0)
    parser.add_argument("--maximum-defect-rows", type=int, default=192)
    parser.add_argument("--maximum-region-zeros", type=int, default=4096)
    parser.add_argument("--log-surprisal", type=float, default=0.74)
    parser.add_argument("--log-r-min", type=float, default=4.25)
    parser.add_argument("--log-r-max", type=float, default=5.25)
    parser.add_argument("--log-r-step", type=float, default=0.03125)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(json.dumps(payload["measure_decomposition"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
