#!/usr/bin/env python3
"""Bound one endpoint-tail outer block mixed with central-body blocks.

The tail block is uniform on its weight-h shell after the per-block output
permutation.  Compare that shell with one ambient-uniform reference block by
conditioning on its weight.  The body blocks are then changed from ambient
uniform words to the supplied central-body spectrum with Holder.  This avoids
an additive-shift lemma and does not compose a forced tail through entrywise
inner transfer envelopes.
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


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    bulk = json.loads(args.bulk.read_text(encoding="utf-8"))
    parameters = bulk["parameters"]
    outer_bits = int(parameters["outer_bits"])
    outer_blocks = int(parameters["outer_blocks"])
    dimension = outer_bits // 2
    if not 1 <= args.tail_endpoint_width < outer_bits // 2:
        raise ValueError("tail endpoint width must be below half the block")

    inner_by_occupation = {
        int(row["active_regular_outer_blocks"]):
        float(row["inner_log2_upper"]) * LOG2
        for row in bulk["occupation_rows"]
    }
    if any(
        occupation not in inner_by_occupation
        for occupation in range(2, outer_blocks + 1)
    ):
        raise ValueError("bulk receipt must cover occupations 2..outer_blocks")

    spectrum = load_outer_spectrum_logs(args.outer_spectrum, outer_bits)
    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - outer_bits * LOG2
    tail_weights: list[int] = []
    tail_log_counts: list[float] = []
    tail_shell_penalties: list[float] = []
    body_log_reference: list[float] = []
    body_log_likelihood: list[float] = []
    for weight in range(1, outer_bits + 1):
        multiplicity = float(spectrum[weight])
        if not math.isfinite(multiplicity):
            continue
        low_tail = weight <= args.tail_endpoint_width
        high_tail = weight >= outer_bits - args.tail_endpoint_width
        selected_tail = (
            (args.tail_side == "both" and (low_tail or high_tail))
            or (args.tail_side == "low" and low_tail)
            or (args.tail_side == "high" and high_tail)
        )
        shell = log_choose(outer_bits, weight)
        if selected_tail:
            tail_weights.append(weight)
            tail_log_counts.append(multiplicity)
            # -log Pr[Wt(U)=h] for U uniform in {0,1}^outer_bits.
            tail_shell_penalties.append(outer_bits * LOG2 - shell)
        else:
            body_log_reference.append(shell - outer_bits * LOG2)
            body_log_likelihood.append(
                multiplicity - shell - reference_log_density
            )
    if not tail_weights or not body_log_reference:
        raise ValueError("split must leave nonempty tail and body")

    weights = np.asarray(tail_weights, dtype=np.int64)
    log_counts = np.asarray(tail_log_counts, dtype=np.float64)
    shell_penalties = np.asarray(tail_shell_penalties, dtype=np.float64)
    body_reference = np.asarray(body_log_reference, dtype=np.float64)
    body_likelihood = np.asarray(body_log_likelihood, dtype=np.float64)

    p_values = 1.0 + np.exp2(
        np.linspace(-16.0, 16.0, args.holder_grid_points)
    )
    inverse_p = 1.0 / p_values
    event_exponents = 1.0 - inverse_p
    body_log_moments = np.asarray(
        [
            logsumexp(body_reference + float(p) * body_likelihood)
            for p in p_values
        ],
        dtype=np.float64,
    )

    contribution_logs: list[float] = []
    occupation_rows: list[dict[str, object]] = []
    for body_count in range(1, outer_blocks):
        reference_inner = inner_by_occupation[body_count + 1]
        # Conditioning an ambient-uniform tail word on Wt=h costs the exact
        # inverse shell probability.  Cap the resulting event bound at one.
        shell_event = np.minimum(0.0, reference_inner + shell_penalties)
        corrections = (
            body_count * body_log_moments[None, :] * inverse_p[None, :]
            + shell_event[:, None] * event_exponents[None, :]
        )
        holder_indices = np.argmin(corrections, axis=1)
        best_corrections = corrections[
            np.arange(len(weights)), holder_indices
        ]
        location_log = (
            log_choose(outer_blocks, body_count)
            + math.log(outer_blocks - body_count)
        )
        values = (
            location_log
            + body_count * log_mass
            + log_counts
            + best_corrections
        )
        contribution_logs.extend(values.tolist())
        dominant_index = int(np.argmax(values))
        occupation_rows.append(
            {
                "body_blocks": body_count,
                "dominant_tail_weight": int(weights[dominant_index]),
                "tail_spectrum_log2": float(
                    log_counts[dominant_index] / LOG2
                ),
                "tail_shell_conditioning_bits": float(
                    shell_penalties[dominant_index] / LOG2
                ),
                "reference_inner_log2_upper": float(reference_inner / LOG2),
                "shell_conditioned_event_log2_upper": float(
                    shell_event[dominant_index] / LOG2
                ),
                "best_holder_p": float(
                    p_values[holder_indices[dominant_index]]
                ),
                "holder_correction_log2": float(
                    best_corrections[dominant_index] / LOG2
                ),
                "pointwise_log2_upper": float(values[dominant_index] / LOG2),
                "pointwise_margin_bits": float(-values[dominant_index] / LOG2),
                "tail_weight_sum_log2_upper": float(
                    logsumexp(values) / LOG2
                ),
            }
        )

    aggregate = float(logsumexp(np.asarray(contribution_logs)))
    dominant_rows = sorted(
        occupation_rows,
        key=lambda row: float(row["tail_weight_sum_log2_upper"]),
        reverse=True,
    )[:20]
    return {
        "schema": "riffle-one-tail-body-shell-conditioned-v1",
        "inputs": {
            "bulk_uniform_reference": str(args.bulk),
            "outer_spectrum": str(args.outer_spectrum),
        },
        "parameters": {
            "outer_bits": outer_bits,
            "outer_dimension": dimension,
            "outer_blocks": outer_blocks,
            "tail_endpoint_width": args.tail_endpoint_width,
            "tail_side": args.tail_side,
            "tail_weight_count": len(weights),
            "holder_grid_points": args.holder_grid_points,
        },
        "aggregate_log2_upper": aggregate / LOG2,
        "aggregate_margin_bits": -aggregate / LOG2,
        "dominant_rows": dominant_rows,
        "body_occupation_rows": occupation_rows,
        "scope": (
            "Exactly one selected endpoint-tail block and at least one "
            "central-body block. The per-block output permutation makes a "
            "weight-h tail word uniform on its shell. Its inner bad-event "
            "probability is compared to the recorded ambient-uniform "
            "occupation-(b+1) bound by exact shell conditioning; central-body "
            "blocks use their complete spectrum through finite-p Holder. "
            "Inherits the supplied inner-model assumptions and nearest-"
            "binary64 arithmetic; not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, required=True)
    parser.add_argument("--outer-spectrum", type=Path, required=True)
    parser.add_argument("--tail-endpoint-width", type=int, required=True)
    parser.add_argument(
        "--tail-side", choices=("low", "high", "both"), default="both"
    )
    parser.add_argument("--holder-grid-points", type=int, default=257)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"margin_bits,{result['aggregate_margin_bits']:.9f}")
    for row in result["dominant_rows"][:5]:
        print(
            "dominant,"
            f"body_blocks,{row['body_blocks']},"
            f"tail_weight,{row['dominant_tail_weight']},"
            f"p,{row['best_holder_p']:.9f},"
            f"margin,{row['pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
