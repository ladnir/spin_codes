#!/usr/bin/env python3
"""Tail/body Riffle diagnostic under a shift-robust body hypothesis.

Endpoint-tail BA words are counted exactly. Central-body words use their full
spectrum through an optimized Holder moment. The calculation asks what follows
if the uniform-reference inner bound for b body words remains valid after
adding an arbitrary sum of tail words. This is a diagnostic for the missing
lemma, not a proof that the lemma holds.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    load_outer_spectrum_logs,
    log_choose,
    logsumexp,
    log_two_power_minus_one,
)


def log_expm1(value: float) -> float:
    if value > 50.0:
        return value + math.log1p(-math.exp(-value))
    return math.log(math.expm1(value))


def log_subtract(left: float, right: float) -> float:
    if right == -math.inf:
        return left
    if right >= left:
        if right - left < 1e-10:
            return -math.inf
        raise ArithmeticError("log subtraction would be negative")
    return left + math.log1p(-math.exp(right - left))


def evaluate(
    bulk_path: Path,
    spectrum_path: Path,
    one_active_path: Path,
    tail_endpoint_width: int,
    holder_grid_points: int,
) -> dict[str, object]:
    bulk = json.loads(bulk_path.read_text(encoding="utf-8"))
    one_active = json.loads(one_active_path.read_text(encoding="utf-8"))
    parameters = bulk["parameters"]
    outer_bits = int(parameters["outer_bits"])
    dimension = outer_bits // 2
    outer_blocks = int(parameters["outer_blocks"])
    if not 1 <= tail_endpoint_width < outer_bits // 2:
        raise ValueError("tail endpoint width must be below half the block")

    spectrum = load_outer_spectrum_logs(spectrum_path, outer_bits)
    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - outer_bits * LOG2
    tail_logs = []
    body_log_reference = []
    body_log_likelihood = []
    body_eta = -math.inf
    body_eta_weight = -1
    for weight in range(1, outer_bits + 1):
        multiplicity = float(spectrum[weight])
        if not math.isfinite(multiplicity):
            continue
        endpoint_distance = min(weight, outer_bits - weight)
        if endpoint_distance <= tail_endpoint_width:
            tail_logs.append(multiplicity)
        else:
            ratio = (
                multiplicity
                - log_choose(outer_bits, weight)
                - reference_log_density
            )
            body_log_reference.append(
                log_choose(outer_bits, weight) - outer_bits * LOG2
            )
            body_log_likelihood.append(ratio)
            if ratio > body_eta:
                body_eta = ratio
                body_eta_weight = weight
    if not tail_logs or body_eta_weight < 0:
        raise ValueError("split must leave nonempty tail and body")
    tail_log_count = float(logsumexp(np.asarray(tail_logs)))
    tail_count = math.exp(tail_log_count)
    log_one_plus_tail = math.log1p(tail_count)
    if holder_grid_points < 3:
        raise ValueError("Holder grid needs at least three points")
    p_values = 1.0 + np.exp2(
        np.linspace(-16.0, 16.0, holder_grid_points)
    )
    inverse_p = 1.0 / p_values
    event_exponents = 1.0 - inverse_p
    body_reference = np.asarray(body_log_reference, dtype=np.float64)
    body_likelihood = np.asarray(body_log_likelihood, dtype=np.float64)
    body_log_moments = np.asarray(
        [
            logsumexp(body_reference + float(p) * body_likelihood)
            for p in p_values
        ],
        dtype=np.float64,
    )

    inner_by_body = {
        int(row["active_regular_outer_blocks"]):
        float(row["inner_log2_upper"]) * LOG2
        for row in bulk["occupation_rows"]
    }
    mixed_rows = []
    mixed_logs = []
    body_only_logs = []
    for body_blocks in range(1, outer_blocks + 1):
        inner_log = inner_by_body[body_blocks]
        corrections = (
            body_blocks * body_log_moments * inverse_p
            + inner_log * event_exponents
        )
        best_index = int(np.argmin(corrections))
        best_correction = float(corrections[best_index])
        body_only_contribution = (
            log_choose(outer_blocks, body_blocks)
            + body_blocks * log_mass
            + best_correction
        )
        contribution = (
            body_only_contribution
            + (outer_blocks - body_blocks) * log_one_plus_tail
        )
        body_only_logs.append(body_only_contribution)
        mixed_logs.append(contribution)
        mixed_rows.append(
            {
                "body_blocks": body_blocks,
                "reference_inner_log2_upper": inner_log / LOG2,
                "best_holder_p": float(p_values[best_index]),
                "body_holder_correction_log2": best_correction / LOG2,
                "body_only_pointwise_log2_upper": (
                    body_only_contribution / LOG2
                ),
                "body_only_pointwise_margin_bits": (
                    -body_only_contribution / LOG2
                ),
                "pointwise_log2_upper": contribution / LOG2,
                "pointwise_margin_bits": -contribution / LOG2,
            }
        )
    mixed_log = float(logsumexp(np.asarray(mixed_logs)))
    body_only_log = float(logsumexp(np.asarray(body_only_logs)))

    # Exact supplied one-active contributions restricted to the endpoint tail.
    one_tail_logs = []
    for row in one_active["weight_rows"]:
        weight = int(row["outer_weight"])
        if min(weight, outer_bits - weight) <= tail_endpoint_width:
            one_tail_logs.append(float(row["pointwise_log2_upper"]) * LOG2)
    one_tail_log = (
        float(logsumexp(np.asarray(one_tail_logs)))
        if one_tail_logs
        else -math.inf
    )

    # Count-only bound for two or more tail blocks. It is deliberately crude
    # but isolates the remaining small-tail obligation.
    all_nonempty_tail_log = log_expm1(outer_blocks * log_one_plus_tail)
    one_tail_count_log = math.log(outer_blocks) + tail_log_count
    two_tail_count_log = log_choose(outer_blocks, 2) + 2.0 * tail_log_count
    two_plus_tail_count_log = log_subtract(
        all_nonempty_tail_log, one_tail_count_log
    )
    one_or_two_count_log = float(
        logsumexp(np.asarray([one_tail_count_log, two_tail_count_log]))
    )
    three_plus_tail_count_log = log_subtract(
        all_nonempty_tail_log, one_or_two_count_log
    )
    tail_only_log = float(
        logsumexp(np.asarray([one_tail_log, two_plus_tail_count_log]))
    )
    combined_diagnostic_log = float(
        logsumexp(np.asarray([mixed_log, tail_only_log]))
    )
    dominant_mixed = sorted(
        mixed_rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]
    return {
        "schema": "riffle-bulk-tail-shift-diagnostic-v1",
        "inputs": {
            "bulk_uniform_reference": str(bulk_path),
            "outer_spectrum": str(spectrum_path),
            "one_active": str(one_active_path),
        },
        "parameters": {
            "outer_bits": outer_bits,
            "outer_dimension": dimension,
            "outer_blocks": outer_blocks,
            "tail_endpoint_width": tail_endpoint_width,
            "holder_grid_size": holder_grid_points,
        },
        "split": {
            "tail_log2_expected_count_per_outer_block": tail_log_count / LOG2,
            "tail_log2_probability_given_nonzero_message": (
                tail_log_count - log_mass
            ) / LOG2,
            "body_log2_eta": body_eta / LOG2,
            "body_eta_weight": body_eta_weight,
        },
        "mixed_under_shift_hypothesis_log2_upper": mixed_log / LOG2,
        "mixed_under_shift_hypothesis_margin_bits": -mixed_log / LOG2,
        "body_only_log2_upper": body_only_log / LOG2,
        "body_only_margin_bits": -body_only_log / LOG2,
        "tail_only": {
            "exact_one_active_log2_upper": one_tail_log / LOG2,
            "exact_one_active_margin_bits": -one_tail_log / LOG2,
            "two_plus_count_only_log2_upper": two_plus_tail_count_log / LOG2,
            "two_plus_count_only_margin_bits": -two_plus_tail_count_log / LOG2,
            "exact_two_count_only_log2_upper": two_tail_count_log / LOG2,
            "exact_two_count_only_margin_bits": -two_tail_count_log / LOG2,
            "three_plus_count_only_log2_upper": three_plus_tail_count_log / LOG2,
            "three_plus_count_only_margin_bits": -three_plus_tail_count_log / LOG2,
            "combined_log2_upper": tail_only_log / LOG2,
            "combined_margin_bits": -tail_only_log / LOG2,
        },
        "combined_diagnostic_log2_upper": combined_diagnostic_log / LOG2,
        "combined_diagnostic_margin_bits": -combined_diagnostic_log / LOG2,
        "dominant_mixed_rows": dominant_mixed,
        "scope": (
            "The body-only bound uses complete central-body spectrum moments "
            "and the recorded uniform-reference inner bound directly. The "
            "mixed bound uses the same calculation and "
            "assumes its uniform-reference inner bound is unchanged by an "
            "arbitrary additive endpoint-tail shift. "
            "That shift-robustness statement is not established. Tail-only "
            "one-active terms use the supplied exact spectrum calculation; "
            "two-or-more tail blocks use count only. Nearest-binary64 diagnostic."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, required=True)
    parser.add_argument("--outer-spectrum", type=Path, required=True)
    parser.add_argument("--one-active", type=Path, required=True)
    parser.add_argument("--tail-endpoint-width", type=int, required=True)
    parser.add_argument("--holder-grid-points", type=int, default=257)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(
        args.bulk,
        args.outer_spectrum,
        args.one_active,
        args.tail_endpoint_width,
        args.holder_grid_points,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    split = result["split"]
    tail = result["tail_only"]
    print(
        "split,"
        f"tail_count_log2,{split['tail_log2_expected_count_per_outer_block']:.9f},"
        f"body_eta_log2,{split['body_log2_eta']:.9f},"
        f"body_eta_weight,{split['body_eta_weight']}"
    )
    print(
        "mixed_shift_hypothesis_margin,"
        f"{result['mixed_under_shift_hypothesis_margin_bits']:.9f}"
    )
    print(
        "tail_one_margin,"
        f"{tail['exact_one_active_margin_bits']:.9f},"
        "tail_two_plus_count_only_margin,"
        f"{tail['two_plus_count_only_margin_bits']:.9f}"
    )
    print(
        "combined_diagnostic_margin,"
        f"{result['combined_diagnostic_margin_bits']:.9f}"
    )
    for row in result["dominant_mixed_rows"][:3]:
        print(
            "dominant_mixed,"
            f"body_blocks,{row['body_blocks']},"
            f"p,{row['best_holder_p']:.9f},"
            f"margin,{row['pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
