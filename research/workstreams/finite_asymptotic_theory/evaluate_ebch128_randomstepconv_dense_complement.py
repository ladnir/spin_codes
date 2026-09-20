#!/usr/bin/env python3
"""Near-full-occupation diagnostic using exact complement coefficients."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp

import evaluate_ebch128_randomstepconv_g1 as base
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch128_randomstepconv_s30_pow2_dense_complement.json"


def log_matrix_sum(log_matrix: np.ndarray) -> float:
    return float(np.logaddexp(log_matrix[0, 0], log_matrix[0, 1]))


def upper_pivot_matrix(
    region_by_zeros: np.ndarray,
    *,
    inactive_rows: int,
    active_rows: int,
    log_r: float,
    candidate_row_norm: float,
) -> np.ndarray:
    maximum_zeros = len(region_by_zeros) - 1
    maximum_pivot = min(active_rows, maximum_zeros - inactive_rows)
    pivots = np.arange(maximum_pivot + 1, dtype=np.float64)
    weights = pivots * log_r - gammaln(pivots + 1.0)
    exact_terms = (
        region_by_zeros[inactive_rows : inactive_rows + maximum_pivot + 1]
        + weights[:, None, None]
    )
    result = logsumexp(exact_terms, axis=0)

    if maximum_pivot < active_rows:
        tail_pivots = np.arange(
            maximum_pivot + 1, active_rows + 1, dtype=np.float64
        )
        tail_terms = (
            active_rows * math.log(candidate_row_norm)
            + tail_pivots * (log_r - math.log(candidate_row_norm))
            - gammaln(tail_pivots + 1.0)
        )
        tail = float(logsumexp(tail_terms))
        result = np.logaddexp(result, tail)
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    base.configure_outer_rows(
        base.ACTIVE_ROWS,
        args.distance_numerator,
        args.distance_denominator,
    )
    spectrum = base.load_spectrum(base.SPECTRUM)
    envelope, _ = base.body_envelope(spectrum)
    effective_mass_log = math.log(envelope + 1)
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

    rows: list[dict[str, object]] = []
    for inactive in range(args.maximum_inactive_rows + 1):
        q = base.L - inactive
        best = math.inf
        best_log_r = math.nan
        for log_r in np.arange(
            args.log_r_min,
            args.log_r_max + args.log_r_step / 2.0,
            args.log_r_step,
        ):
            pivot_matrix = upper_pivot_matrix(
                region_by_zeros,
                inactive_rows=inactive,
                active_rows=q,
                log_r=float(log_r),
                candidate_row_norm=candidate_row_norm,
            )
            all_regions = base.log_power(pivot_matrix, base.B)
            conditional = (
                math.lgamma(q + 1)
                - q * math.log(base.B)
                - q * float(log_r)
                + log_matrix_sum(all_regions)
            )
            inner = min(0.0, conditional + base.D * surprisal)
            outer = (
                math.lgamma(base.L + 1)
                - math.lgamma(q + 1)
                - math.lgamma(base.L - q + 1)
                + q * effective_mass_log
            )
            pointwise = (outer + inner) / base.LOG2
            if pointwise < best:
                best = pointwise
                best_log_r = float(log_r)
        rows.append(
            {
                "active_outer_rows": q,
                "inactive_outer_rows": inactive,
                "pointwise_log2_upper": best,
                "log_pivot_fugacity": best_log_r,
                "surprisal_hex": float(surprisal).hex(),
                "pivot_fugacity_hex": float(math.exp(best_log_r)).hex(),
            }
        )

    return {
        "schema": "ebch128-randomstepconv-dense-complement-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC_WITH_ANALYTIC_TAIL_CAP",
        "parameters": {
            "outer_rows": base.L,
            "active_message_rows": base.ACTIVE_ROWS,
            "output_bits": base.N,
            "distance_cutoff": base.D,
            "memory_bits": args.memory_bits,
            "maximum_inactive_rows": args.maximum_inactive_rows,
            "maximum_exact_region_zeros": args.maximum_region_zeros,
            "log_surprisal": args.log_surprisal,
        },
        "method": {
            "region_transfer": "exact normalized coefficient indexed by the number of noncandidate positions",
            "pivot_average": "one positive coefficient fugacity after the exact region transfer",
            "tail": "entrywise row-sum cap for pivot loads beyond the exact complement table",
            "outer_spectrum": "existing worst-shell uniform-even pointwise majorant",
        },
        "occupation_rows": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "The Chernoff witness is fixed across the reported occupations.",
            "The outer spectrum still uses the pointwise body majorant.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=30)
    parser.add_argument("--distance-numerator", type=int, default=11)
    parser.add_argument("--distance-denominator", type=int, default=100)
    parser.add_argument("--maximum-inactive-rows", type=int, default=97)
    parser.add_argument("--maximum-region-zeros", type=int, default=512)
    parser.add_argument("--log-surprisal", type=float, default=0.75)
    parser.add_argument("--log-r-min", type=float, default=6.0)
    parser.add_argument("--log-r-max", type=float, default=11.0)
    parser.add_argument("--log-r-step", type=float, default=0.125)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    rows = payload["occupation_rows"]
    print(json.dumps({"first": rows[0], "last": rows[-1]}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
