#!/usr/bin/env python3
"""One endpoint-tail BA block plus central-body blocks in Riffle.

The body blocks are uniform-reference words corrected by the complete central
body spectrum through Holder moments. The one tail block is represented by iid
region marks and converted to its exact fixed-weight support distribution by
the binomial conditioning penalty. This removes an arbitrary-shift assumption
for the exactly-one-tail mixed class.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    binomial_transform,
    candidate_epoch_logs,
    load_nonzero_spectrum,
    load_outer_spectrum_logs,
    load_uniform_nonactivation,
    log_choose,
    logsumexp,
    log_two_power_minus_one,
    regular_region_log_matrices,
    splitstate_impulse_matrices,
)
from analyze_riffle_two_tail_iid_regions import log_matrix_power_pairwise


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    region_bits = output_bits // args.outer_bits
    epochs = region_bits // args.step_bits
    if region_bits % args.step_bits:
        raise ValueError("step size must divide one region")

    spectrum = load_outer_spectrum_logs(args.outer_spectrum, args.outer_bits)
    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - args.outer_bits * LOG2
    tail_weights = []
    tail_log_counts = []
    body_log_reference = []
    body_log_likelihood = []
    for weight in range(1, args.outer_bits + 1):
        multiplicity = float(spectrum[weight])
        if not math.isfinite(multiplicity):
            continue
        shell = log_choose(args.outer_bits, weight)
        is_low_tail = weight <= args.tail_endpoint_width
        is_high_tail = weight >= args.outer_bits - args.tail_endpoint_width
        selected_tail = (
            (args.tail_side == "both" and (is_low_tail or is_high_tail))
            or (args.tail_side == "low" and is_low_tail)
            or (args.tail_side == "high" and is_high_tail)
        )
        if selected_tail:
            tail_weights.append(weight)
            tail_log_counts.append(multiplicity)
        else:
            body_log_reference.append(shell - args.outer_bits * LOG2)
            body_log_likelihood.append(
                multiplicity - shell - reference_log_density
            )
    weights = np.asarray(tail_weights, dtype=np.int64)
    log_counts = np.asarray(tail_log_counts, dtype=np.float64)
    if not len(weights) or not body_log_reference:
        raise ValueError("split must leave tail and body rows")

    mark_probabilities = weights.astype(np.float64) / args.outer_bits
    with np.errstate(divide="ignore"):
        log_mark = np.log(mark_probabilities)
        log_no_mark = np.log1p(-mark_probabilities)
    fixed_weight_log_probability = np.empty(len(weights), dtype=np.float64)
    for index, weight in enumerate(weights):
        probability = float(weight) / args.outer_bits
        if weight == args.outer_bits:
            fixed_weight_log_probability[index] = 0.0
        else:
            fixed_weight_log_probability[index] = (
                log_choose(args.outer_bits, int(weight))
                + weight * math.log(probability)
                + (args.outer_bits - weight) * math.log1p(-probability)
            )
    # In addition to conditioning the iid region marks on total support h,
    # fairize the h nonzero tail values and condition all of them back to one.
    # The latter event has probability 2^-h and avoids subtracting entrywise
    # transfer envelopes to manufacture a forced candidate.
    conditioning_penalty = -fixed_weight_log_probability + weights * LOG2

    p_values = 1.0 + np.exp2(
        np.linspace(-16.0, 16.0, args.holder_grid_points)
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

    nonactivation = load_uniform_nonactivation(args.activation, args.step_bits)
    live_spectrum = load_nonzero_spectrum(
        args.live_spectrum, args.step_bits, args.state_bits
    )
    transform = binomial_transform(args.step_bits)
    body_counts = np.arange(1, outer_blocks, dtype=np.int64)
    best_inner = np.full(
        (len(body_counts), len(weights)), math.inf, dtype=np.float64
    )
    best_tilt = np.full_like(best_inner, math.nan)
    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        impulses = splitstate_impulse_matrices(
            z=z,
            step_bits=args.step_bits,
            state_bits=args.state_bits,
            constituent_distance=args.constituent_distance,
            live_moment_order=args.live_moment_order,
            live_model="support-averaged-preaddmul",
            nonactivation=nonactivation,
            live_spectrum=live_spectrum,
        )
        fair_candidate_logs = candidate_epoch_logs(impulses, transform)
        regions = regular_region_log_matrices(
            fair_candidate_logs,
            args.step_bits,
            epochs,
            outer_blocks,
        )
        zero_rows = regions[1:outer_blocks]
        # A marked region contains one additional fair tail candidate. Actual
        # tail ones are recovered by the 2^h conditioning penalty above.
        marked_rows = regions[2 : outer_blocks + 1]
        mixed = np.logaddexp(
            zero_rows[:, None, :, :] + log_no_mark[None, :, None, None],
            marked_rows[:, None, :, :] + log_mark[None, :, None, None],
        ).reshape(-1, 2, 2)
        powered = log_matrix_power_pairwise(mixed, args.outer_bits)
        moments = np.logaddexp(
            powered[:, 0, 0], powered[:, 0, 1]
        ).reshape(len(body_counts), len(weights))
        bounds = (
            moments
            + conditioning_penalty[None, :]
            + distance * surprisal
        )
        improved = bounds < best_inner
        best_inner[improved] = bounds[improved]
        best_tilt[improved] = log_surprisal
        print(f"tilt,{log_surprisal:.6f}", flush=True)

    body_holder = np.empty_like(best_inner)
    best_holder_p = np.empty_like(best_inner)
    for body_index, body_count in enumerate(body_counts):
        corrections = (
            body_count * body_log_moments[None, :] * inverse_p[None, :]
            + np.minimum(0.0, best_inner[body_index, :])[:, None]
            * event_exponents[None, :]
        )
        indices = np.argmin(corrections, axis=1)
        body_holder[body_index, :] = corrections[
            np.arange(len(weights)), indices
        ]
        best_holder_p[body_index, :] = p_values[indices]

    contribution_rows = []
    contribution_logs = []
    for body_index, body_count in enumerate(body_counts):
        location_log = (
            log_choose(outer_blocks, int(body_count))
            + math.log(outer_blocks - int(body_count))
        )
        values = (
            location_log
            + body_count * log_mass
            + log_counts
            + body_holder[body_index, :]
        )
        contribution_logs.extend(values.tolist())
        dominant_index = int(np.argmax(values))
        contribution_rows.append(
            {
                "body_blocks": int(body_count),
                "dominant_tail_weight": int(weights[dominant_index]),
                "best_log_surprisal": float(
                    best_tilt[body_index, dominant_index]
                ),
                "best_holder_p": float(
                    best_holder_p[body_index, dominant_index]
                ),
                "reference_inner_log2_upper": float(
                    min(0.0, best_inner[body_index, dominant_index]) / LOG2
                ),
                "conditioning_penalty_bits": float(
                    conditioning_penalty[dominant_index] / LOG2
                ),
                "body_holder_correction_log2": float(
                    body_holder[body_index, dominant_index] / LOG2
                ),
                "tail_spectrum_log2": float(
                    log_counts[dominant_index] / LOG2
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
        contribution_rows,
        key=lambda row: float(row["tail_weight_sum_log2_upper"]),
        reverse=True,
    )[:20]
    return {
        "schema": "riffle-one-tail-body-conditioned-v1",
        "inputs": {
            "outer_spectrum": str(args.outer_spectrum),
            "activation": str(args.activation),
            "live_spectrum": str(args.live_spectrum),
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "tail_endpoint_width": args.tail_endpoint_width,
            "tail_side": args.tail_side,
            "tail_weight_count": int(len(weights)),
            "holder_grid_points": args.holder_grid_points,
            "log_surprisals": list(args.log_surprisals),
        },
        "aggregate_log2_upper": aggregate / LOG2,
        "aggregate_margin_bits": -aggregate / LOG2,
        "dominant_rows": dominant_rows,
        "body_occupation_rows": contribution_rows,
        "scope": (
            "Exactly one endpoint-tail block and at least one central-body "
            "block. Tail support is bounded by exact binomial conditioning of "
            "iid region marks; tail values are fairized and conditioned back "
            "to all ones at an exact h-bit cost. Body blocks use the complete "
            "body spectrum through Holder moments. Inherits the recorded "
            "RM2Sub transfer assumptions; nearest binary64, not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=1024)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--tail-endpoint-width", type=int, required=True)
    parser.add_argument(
        "--tail-side", choices=("low", "high", "both"), default="both"
    )
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--state-bits", type=int, default=16)
    parser.add_argument("--constituent-distance", type=int, default=48)
    parser.add_argument("--live-moment-order", type=int, default=3)
    parser.add_argument("--holder-grid-points", type=int, default=129)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0),
    )
    parser.add_argument("--outer-spectrum", type=Path, required=True)
    parser.add_argument("--activation", type=Path, required=True)
    parser.add_argument("--live-spectrum", type=Path, required=True)
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
            f"margin,{row['pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
