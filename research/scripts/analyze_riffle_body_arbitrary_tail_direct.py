#!/usr/bin/env python3
"""Bound central-body Riffle blocks under arbitrary endpoint-tail inputs.

For one transposed region, let F[w] be the inner transfer bound conditioned on
a uniformly located input support of weight w.  With b ambient-uniform body
bits and j forced tail ones in distinct positions, the total support has size
j+Bin(b,1/2) and is uniform conditional on that size.  Hence

    G[b,j] = E[F[j+Bin(b,1/2)]].

The recurrence G[b+1,j]=(G[b,j]+G[b,j+1])/2 computes every row.  Taking the
entrywise maximum over j gives a valid transfer envelope for an arbitrary tail
pattern in each region.  The complete central-body spectrum is then inserted
through Holder, while endpoint-tail codewords are counted exactly in total.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    load_nonzero_spectrum,
    load_outer_spectrum_logs,
    load_uniform_nonactivation,
    log_choose,
    log_matrix_power_moments_batch,
    logsumexp,
    log_two_power_minus_one,
    regular_region_log_matrices,
    splitstate_impulse_matrices,
)


def log_expm1(value: float) -> float:
    if value > 50.0:
        return value + math.log1p(-math.exp(-value))
    return math.log(math.expm1(value))


def arbitrary_tail_region_envelopes(
    fixed_weight_logs: np.ndarray,
    maximum_body_blocks: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return entrywise max_j E[F[j+Bin(b,1/2)]] for b=1..max."""
    current = fixed_weight_logs.copy()
    envelopes = np.empty((maximum_body_blocks, 2, 2), dtype=np.float64)
    maximizers = np.empty((maximum_body_blocks, 4), dtype=np.int64)
    for body_count in range(1, maximum_body_blocks + 1):
        current = np.logaddexp(current[:-1], current[1:]) - LOG2
        envelopes[body_count - 1] = np.max(current, axis=0)
        maximizers[body_count - 1] = np.argmax(
            current.reshape(len(current), 4), axis=0
        )
    return envelopes, maximizers


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    region_bits = output_bits // args.outer_bits
    if region_bits != outer_blocks:
        raise ValueError("this analyzer expects one region position per block")
    if region_bits % args.step_bits:
        raise ValueError("step size must divide one region")
    epochs = region_bits // args.step_bits

    spectrum = load_outer_spectrum_logs(args.outer_spectrum, args.outer_bits)
    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - args.outer_bits * LOG2
    tail_logs: list[float] = []
    body_log_reference: list[float] = []
    body_log_likelihood: list[float] = []
    for weight in range(1, args.outer_bits + 1):
        multiplicity = float(spectrum[weight])
        if not math.isfinite(multiplicity):
            continue
        shell = log_choose(args.outer_bits, weight)
        if min(weight, args.outer_bits - weight) <= args.tail_endpoint_width:
            tail_logs.append(multiplicity)
        else:
            body_log_reference.append(shell - args.outer_bits * LOG2)
            body_log_likelihood.append(
                multiplicity - shell - reference_log_density
            )
    if not tail_logs or not body_log_reference:
        raise ValueError("split must leave nonempty tail and body")
    tail_log_count = float(logsumexp(np.asarray(tail_logs)))
    tail_count = math.exp(tail_log_count)
    log_one_plus_tail = math.log1p(tail_count)

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
    body_counts = np.arange(1, outer_blocks, dtype=np.int64)
    best_inner = np.full(len(body_counts), math.inf, dtype=np.float64)
    best_tilt = np.full(len(body_counts), math.nan, dtype=np.float64)
    best_worst_j = np.full((len(body_counts), 4), -1, dtype=np.int64)
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
        with np.errstate(divide="ignore"):
            impulse_logs = np.log(impulses)
        fixed_weight_logs = regular_region_log_matrices(
            impulse_logs,
            args.step_bits,
            epochs,
            outer_blocks,
        )
        envelopes, maximizers = arbitrary_tail_region_envelopes(
            fixed_weight_logs, outer_blocks - 1
        )
        moments = log_matrix_power_moments_batch(envelopes, args.outer_bits)
        bounds = moments + distance * surprisal
        improved = bounds < best_inner
        best_inner[improved] = bounds[improved]
        best_tilt[improved] = log_surprisal
        best_worst_j[improved] = maximizers[improved]
        print(f"tilt,{log_surprisal:.6f}", flush=True)

    rows: list[dict[str, object]] = []
    contribution_logs: list[float] = []
    for body_index, body_count in enumerate(body_counts):
        inner_log = min(0.0, float(best_inner[body_index]))
        corrections = (
            body_count * body_log_moments * inverse_p
            + inner_log * event_exponents
        )
        holder_index = int(np.argmin(corrections))
        best_correction = float(corrections[holder_index])
        available_tail_blocks = outer_blocks - int(body_count)
        nonempty_tail_choices = log_expm1(
            available_tail_blocks * log_one_plus_tail
        )
        contribution = (
            log_choose(outer_blocks, int(body_count))
            + body_count * log_mass
            + nonempty_tail_choices
            + best_correction
        )
        contribution_logs.append(contribution)
        rows.append(
            {
                "body_blocks": int(body_count),
                "available_tail_blocks": available_tail_blocks,
                "best_log_surprisal": float(best_tilt[body_index]),
                "worst_tail_counts_by_transition": {
                    "zero_to_zero": int(best_worst_j[body_index, 0]),
                    "zero_to_live": int(best_worst_j[body_index, 1]),
                    "live_to_zero": int(best_worst_j[body_index, 2]),
                    "live_to_live": int(best_worst_j[body_index, 3]),
                },
                "arbitrary_tail_inner_log2_upper": inner_log / LOG2,
                "best_holder_p": float(p_values[holder_index]),
                "holder_correction_log2": best_correction / LOG2,
                "tail_choice_log2": nonempty_tail_choices / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
                "pointwise_margin_bits": -contribution / LOG2,
            }
        )

    aggregate = float(logsumexp(np.asarray(contribution_logs)))
    dominant_rows = sorted(
        rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]
    return {
        "schema": "riffle-body-arbitrary-tail-direct-v1",
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
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "holder_grid_points": args.holder_grid_points,
            "log_surprisals": list(args.log_surprisals),
        },
        "tail_log2_expected_count_per_outer_block": tail_log_count / LOG2,
        "aggregate_log2_upper": aggregate / LOG2,
        "aggregate_margin_bits": -aggregate / LOG2,
        "dominant_rows": dominant_rows,
        "body_occupation_rows": rows,
        "scope": (
            "All messages with at least one central-body outer block and at "
            "least one endpoint-tail outer block. The direct region envelope "
            "maximizes over the number of forced tail ones separately in each "
            "region, so it permits arbitrary tail supports. It retains the "
            "independent random permutation within each transposed region and "
            "the ambient-uniform body values. The complete central-body "
            "spectrum is inserted through finite-p Holder, and all nonempty "
            "tail choices are counted. Inherits the recorded RM2Sub transfer "
            "assumptions; nearest binary64, not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=1024)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--tail-endpoint-width", type=int, required=True)
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--state-bits", type=int, default=16)
    parser.add_argument("--constituent-distance", type=int, default=48)
    parser.add_argument("--live-moment-order", type=int, default=3)
    parser.add_argument("--holder-grid-points", type=int, default=257)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(
            -10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
            -5.0, -4.0, -3.0, -2.5, -2.25, -2.0, -1.75, -1.5,
            -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
            0.75, 1.0, 1.25, 1.5,
        ),
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
            f"body,{row['body_blocks']},"
            f"tilt,{row['best_log_surprisal']:.6f},"
            f"p,{row['best_holder_p']:.9f},"
            f"margin,{row['pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
