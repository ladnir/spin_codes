#!/usr/bin/env python3
"""Bound every mixed central-body/endpoint-tail Riffle occupation.

For b central-body blocks and k endpoint-tail blocks, compare each tail shell
with an ambient-uniform reference block by exact shell conditioning.  At a
fixed Holder exponent the sum over all k tail weights factorizes into the k-th
power of one scalar tail moment.  This covers all b>=1, k>=1 without tuple
enumeration.  The shell-conditioned event bound is deliberately left uncapped;
that is weaker but preserves the factorization and remains a valid upper bound.
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

    spectrum = load_outer_spectrum_logs(args.outer_spectrum, outer_bits)
    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - outer_bits * LOG2
    tail_log_counts: list[float] = []
    tail_shell_penalties: list[float] = []
    body_log_reference: list[float] = []
    body_log_likelihood: list[float] = []
    for weight in range(1, outer_bits + 1):
        multiplicity = float(spectrum[weight])
        if not math.isfinite(multiplicity):
            continue
        endpoint = min(weight, outer_bits - weight)
        shell = log_choose(outer_bits, weight)
        if endpoint <= args.tail_endpoint_width:
            tail_log_counts.append(multiplicity)
            tail_shell_penalties.append(outer_bits * LOG2 - shell)
        else:
            body_log_reference.append(shell - outer_bits * LOG2)
            body_log_likelihood.append(
                multiplicity - shell - reference_log_density
            )
    if not tail_log_counts or not body_log_reference:
        raise ValueError("split must leave nonempty tail and body")

    tail_counts = np.asarray(tail_log_counts, dtype=np.float64)
    shell_penalties = np.asarray(tail_shell_penalties, dtype=np.float64)
    body_reference = np.asarray(body_log_reference, dtype=np.float64)
    body_likelihood = np.asarray(body_log_likelihood, dtype=np.float64)
    p_values = 1.0 + np.exp2(
        np.linspace(-16.0, 16.0, args.holder_grid_points)
    )
    inverse_p = 1.0 / p_values
    q_values = 1.0 - inverse_p
    body_log_moments = np.asarray(
        [
            logsumexp(body_reference + float(p) * body_likelihood)
            for p in p_values
        ],
        dtype=np.float64,
    )
    tail_log_moments = np.asarray(
        [
            logsumexp(tail_counts + float(q) * shell_penalties)
            for q in q_values
        ],
        dtype=np.float64,
    )

    all_contributions: list[float] = []
    rows: list[dict[str, object]] = []
    for active_count in range(2, outer_blocks + 1):
        if active_count not in inner_by_occupation:
            raise ValueError(f"bulk receipt lacks occupation {active_count}")
        inner_log = inner_by_occupation[active_count]
        tail_counts_k = np.arange(1, active_count, dtype=np.int64)
        body_counts_b = active_count - tail_counts_k
        corrections = (
            body_counts_b[:, None]
            * body_log_moments[None, :]
            * inverse_p[None, :]
            + q_values[None, :] * inner_log
            + tail_counts_k[:, None] * tail_log_moments[None, :]
        )
        best_indices = np.argmin(corrections, axis=1)
        best_corrections = corrections[
            np.arange(len(tail_counts_k)), best_indices
        ]
        location_logs = np.asarray(
            [
                log_choose(outer_blocks, int(active_count))
                + log_choose(active_count, int(k))
                for k in tail_counts_k
            ],
            dtype=np.float64,
        )
        contributions = (
            location_logs
            + body_counts_b * log_mass
            + best_corrections
        )
        all_contributions.extend(contributions.tolist())
        dominant_index = int(np.argmax(contributions))
        rows.append(
            {
                "active_blocks": active_count,
                "dominant_body_blocks": int(body_counts_b[dominant_index]),
                "dominant_tail_blocks": int(tail_counts_k[dominant_index]),
                "best_holder_p": float(p_values[best_indices[dominant_index]]),
                "reference_inner_log2_upper": float(inner_log / LOG2),
                "holder_and_tail_correction_log2": float(
                    best_corrections[dominant_index] / LOG2
                ),
                "pointwise_log2_upper": float(
                    contributions[dominant_index] / LOG2
                ),
                "pointwise_margin_bits": float(
                    -contributions[dominant_index] / LOG2
                ),
                "tail_count_sum_log2_upper": float(
                    logsumexp(contributions) / LOG2
                ),
            }
        )

    aggregate = float(logsumexp(np.asarray(all_contributions)))
    dominant_rows = sorted(
        rows,
        key=lambda row: float(row["tail_count_sum_log2_upper"]),
        reverse=True,
    )[:20]
    return {
        "schema": "riffle-multi-tail-body-shell-holder-v1",
        "inputs": {
            "bulk_uniform_reference": str(args.bulk),
            "outer_spectrum": str(args.outer_spectrum),
        },
        "parameters": {
            "outer_bits": outer_bits,
            "outer_dimension": dimension,
            "outer_blocks": outer_blocks,
            "tail_endpoint_width": args.tail_endpoint_width,
            "holder_grid_points": args.holder_grid_points,
        },
        "aggregate_log2_upper": aggregate / LOG2,
        "aggregate_margin_bits": -aggregate / LOG2,
        "dominant_rows": dominant_rows,
        "occupation_rows": rows,
        "scope": (
            "All mixed occupations with at least one central-body block and "
            "at least one endpoint-tail block. Each tail word is uniform on "
            "its shell after the per-block permutation and is compared with "
            "an ambient-uniform reference block by exact shell conditioning. "
            "For each finite Holder exponent, the complete tail-weight sum "
            "factorizes and the central-body blocks use their complete "
            "spectrum moment. The conditional event bound is not capped at "
            "one, which is conservative. Inherits the supplied inner-model "
            "assumptions and nearest-binary64 arithmetic; not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, required=True)
    parser.add_argument("--outer-spectrum", type=Path, required=True)
    parser.add_argument("--tail-endpoint-width", type=int, required=True)
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
            f"active,{row['active_blocks']},"
            f"body,{row['dominant_body_blocks']},"
            f"tail,{row['dominant_tail_blocks']},"
            f"p,{row['best_holder_p']:.9f},"
            f"margin,{row['pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
