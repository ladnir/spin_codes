#!/usr/bin/env python3
"""Evaluate occupation one for BA-3 followed by random Toeplitz convolution.

The inner emits zero before the first nonzero routed input bit. Its activation
output is one. Later output bits for a fixed nonzero difference are
independent fair bits. This is the invertible no-reset, packet-width-one
random-convolution model.

The calculation uses the exact expected BA-3 spectrum and the exact
first-activation distribution. Binomial tails are evaluated by SciPy in
nearest binary64. The result is a diagnostic, not an outward-rounded
certificate, and it covers only occupation one.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import bdtr, gammaln, logsumexp

from analyze_golay_ba_rm2sub_joint import expected_ba_log_spectrum


WORKSTREAM = Path(__file__).resolve().parent


def log_binomial(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def region_average_tables(
    outer_bits: int,
    region_length: int,
    distance: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return exact-tail and forced-late averages for each first region.

    Entry j-1 averages over the uniform position of the first nonzero bit in
    region j. The suffix includes the activation position.
    """
    output_bits = outer_bits * region_length
    suffix_bad = np.ones(output_bits, dtype=np.float64)
    if distance < output_bits:
        suffix_lengths = np.arange(distance + 1, output_bits + 1, dtype=np.int64)
        suffix_bad[distance:] = bdtr(distance - 1, suffix_lengths - 1, 0.5)

    forced_late = np.zeros(output_bits, dtype=np.float64)
    forced_late[: min(distance, output_bits)] = 1.0

    # Rows are indexed by increasing suffix length. First region j uses row
    # outer_bits-j, so reverse the row means.
    exact_ascending = suffix_bad.reshape(outer_bits, region_length).mean(axis=1)
    late_ascending = forced_late.reshape(outer_bits, region_length).mean(axis=1)
    return exact_ascending[::-1].copy(), late_ascending[::-1].copy()


def log_average_for_weight(
    weight: int,
    outer_bits: int,
    region_averages: np.ndarray,
) -> float:
    """Average a first-region quantity over a uniform weight-w support."""
    denominator = log_binomial(outer_bits, weight)
    terms = []
    for first_region in range(1, outer_bits - weight + 2):
        average = float(region_averages[first_region - 1])
        if average == 0.0:
            continue
        log_first_region = (
            log_binomial(outer_bits - first_region, weight - 1) - denominator
        )
        terms.append(log_first_region + math.log(average))
    return float(logsumexp(terms)) if terms else -math.inf


def log_capped_row_union(log_single_row: float, outer_rows: int) -> float:
    """Bound any row placement for one reused word by min(1,L*p)."""
    return min(0.0, math.log(outer_rows) + log_single_row)


def margin_bits(log_probability: float) -> float:
    return -log_probability / math.log(2.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-bits", type=int, default=720)
    parser.add_argument("--outer-rows", type=int, default=2944)
    parser.add_argument("--relative-distance", type=float, default=0.11)
    parser.add_argument(
        "--output",
        type=Path,
        default=WORKSTREAM / "ba3_B720_ideal_causal_g1_q1_d11.json",
    )
    args = parser.parse_args()

    if args.outer_bits % 24:
        raise ValueError("outer bits must be divisible by 24")
    output_bits = args.outer_bits * args.outer_rows
    parent_message_bits = (args.outer_bits // 2) * args.outer_rows
    target_message_bits = 2**20
    if parent_message_bits < target_message_bits:
        raise ValueError("parent dimension is smaller than 2^20")
    distance = math.floor(args.relative_distance * output_bits)
    log_spectrum = expected_ba_log_spectrum(args.outer_bits)
    exact_region, late_region = region_average_tables(
        args.outer_bits,
        args.outer_rows,
        distance,
    )

    rows = []
    naive_logs = []
    reuse_aware_logs = []
    for weight in range(1, args.outer_bits + 1):
        log_multiplicity = float(log_spectrum[weight])
        if not math.isfinite(log_multiplicity):
            continue

        log_exact = log_average_for_weight(weight, args.outer_bits, exact_region)
        log_late = log_average_for_weight(weight, args.outer_bits, late_region)
        log_naive = math.log(args.outer_rows) + log_multiplicity + log_exact
        log_reuse_aware = (
            log_multiplicity + log_capped_row_union(log_exact, args.outer_rows)
        )
        naive_logs.append(log_naive)
        reuse_aware_logs.append(log_reuse_aware)
        rows.append(
            {
                "outer_weight": weight,
                "single_row_bad_log2": log_exact / math.log(2.0),
                "single_row_forced_late_log2": log_late / math.log(2.0),
                "expected_outer_multiplicity_log2": (
                    log_multiplicity / math.log(2.0)
                ),
                "naive_row_union_contribution_log2": (
                    log_naive / math.log(2.0)
                ),
                "reuse_aware_contribution_log2": (
                    log_reuse_aware / math.log(2.0)
                ),
            }
        )

    naive_aggregate = float(logsumexp(naive_logs))
    reuse_aware_aggregate = float(logsumexp(reuse_aware_logs))
    rows.sort(key=lambda row: row["reuse_aware_contribution_log2"], reverse=True)

    # This table is not a conditional-ensemble claim. It records only how the
    # same first moment changes when low outer-weight terms are deleted.
    log_by_weight = {
        int(row["outer_weight"]): (
            float(row["reuse_aware_contribution_log2"]) * math.log(2.0)
        )
        for row in rows
    }
    deleted_low_shells = []
    for minimum_weight in range(1, args.outer_bits + 2):
        tail = [
            value for weight, value in log_by_weight.items()
            if weight >= minimum_weight
        ]
        log_tail = float(logsumexp(tail)) if tail else -math.inf
        deleted_low_shells.append(
            {
                "retain_outer_weights_at_least": minimum_weight,
                "reuse_aware_q1_margin_bits": margin_bits(log_tail),
            }
        )
        if margin_bits(log_tail) >= 64.0:
            break

    result = {
        "schema": "ba3-b720-ideal-causal-g1-q1-d11-v2",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "outer": "one sampled Golay--BA-3 code reused in every row",
            "outer_bits": args.outer_bits,
            "outer_dimension": args.outer_bits // 2,
            "outer_rows": args.outer_rows,
            "parent_message_bits": parent_message_bits,
            "target_message_bits": target_message_bits,
            "zero_shortened_input_bits": parent_message_bits - target_message_bits,
            "output_bits": output_bits,
            "bad_output_weight_at_most": distance,
            "relative_distance_target": args.relative_distance,
            "inner": "invertible random lower-triangular Toeplitz convolution",
            "packet_bits": 1,
        },
        "probability_space": [
            "the two BA-3 interleavers sampled once",
            "one independent local coordinate permutation per outer row",
            "one independent permutation within each transposed region",
            "the random Toeplitz inner; for each fixed nonzero difference, its activation output is one and its later outputs are independent fair bits",
        ],
        "q1_bounds": {
            "naive_row_union_log2": naive_aggregate / math.log(2.0),
            "naive_row_union_margin_bits": margin_bits(naive_aggregate),
            "reuse_aware_log2": reuse_aware_aggregate / math.log(2.0),
            "reuse_aware_margin_bits": margin_bits(reuse_aware_aggregate),
        },
        "dominant_rows_reuse_aware": rows[:24],
        "deleted_low_shell_diagnostic": deleted_low_shells,
        "interpretation": [
            "The reuse-aware bound groups all outer-row placements of the same outer codeword before averaging over the one sampled BA code.",
            "The forced-late column is a lower bound on the single-row bad probability for every causal inner because a suffix no longer than the distance threshold cannot have larger weight.",
            "Deleting low shells is an algebraic diagnostic, not conditioning on a verified minimum-distance event.",
        ],
        "limitations": [
            "SciPy binomial tails and all arithmetic use nearest binary64.",
            "The result covers occupation one only.",
            "The no-reset ideal inner is a proof-of-concept benchmark, not the current RM2Sub-S19 inner.",
            "A full repeated-code theorem still requires every occupation and cannot multiply expected BA spectra across rows.",
            "An ensemble-average upper bound above the target does not prove that a typical sampled code fails.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "naive_q1_margin_bits": result["q1_bounds"]["naive_row_union_margin_bits"],
        "reuse_aware_q1_margin_bits": result["q1_bounds"]["reuse_aware_margin_bits"],
        "dominant_outer_weight": rows[0]["outer_weight"],
        "first_deleted_shell_threshold_reaching_40_bits": next(
            (
                row["retain_outer_weights_at_least"]
                for row in deleted_low_shells
                if row["reuse_aware_q1_margin_bits"] >= 40.0
            ),
            None,
        ),
    }, indent=2))


if __name__ == "__main__":
    main()
