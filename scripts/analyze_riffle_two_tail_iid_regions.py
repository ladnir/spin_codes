#!/usr/bin/env python3
"""Two-tail Riffle bound via iid regions plus exact conditioning penalties.

For fixed outer weights h1,h2, replace each independently permuted fixed-size
region support by Bernoulli marks with probabilities h1/B and h2/B. Mix the
exact zero/one/two-forced-candidate region transfer matrices, raise the result
across B regions, and sum the complete endpoint-tail BA spectrum. A uniform
fixed-weight support is the iid support conditioned on its total weight. For a
nonnegative Chernoff moment, division by the two binomial point probabilities
therefore gives a valid fixed-support upper bound without requiring the region
transfer matrices to commute.
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
    logsumexp,
    regular_region_log_matrices,
    splitstate_impulse_matrices,
)


def log_matmul_pairwise(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty_like(left)
    result[:, 0, 0] = np.logaddexp(
        left[:, 0, 0] + right[:, 0, 0],
        left[:, 0, 1] + right[:, 1, 0],
    )
    result[:, 0, 1] = np.logaddexp(
        left[:, 0, 0] + right[:, 0, 1],
        left[:, 0, 1] + right[:, 1, 1],
    )
    result[:, 1, 0] = np.logaddexp(
        left[:, 1, 0] + right[:, 0, 0],
        left[:, 1, 1] + right[:, 1, 0],
    )
    result[:, 1, 1] = np.logaddexp(
        left[:, 1, 0] + right[:, 0, 1],
        left[:, 1, 1] + right[:, 1, 1],
    )
    return result


def log_matrix_power_pairwise(matrices: np.ndarray, power: int) -> np.ndarray:
    if power <= 0 or power & (power - 1):
        raise ValueError("this diagnostic expects a positive power of two")
    result = matrices.copy()
    remaining = power.bit_length() - 1
    for _ in range(remaining):
        result = log_matmul_pairwise(result, result)
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    region_bits = output_bits // args.outer_bits
    if region_bits % args.step_bits:
        raise ValueError("step size must divide one transposed region")
    epochs = region_bits // args.step_bits
    if args.outer_bits & (args.outer_bits - 1):
        raise ValueError("fast matrix powering expects a power-of-two outer size")

    spectrum = load_outer_spectrum_logs(args.outer_spectrum, args.outer_bits)
    weights = np.asarray(
        [
            weight
            for weight in range(1, args.outer_bits + 1)
            if min(weight, args.outer_bits - weight) <= args.tail_endpoint_width
            and math.isfinite(float(spectrum[weight]))
        ],
        dtype=np.int64,
    )
    if not len(weights):
        raise ValueError("tail contains no nonzero spectrum rows")
    log_counts = spectrum[weights]
    first = np.repeat(weights, len(weights))
    second = np.tile(weights, len(weights))
    pair_log_counts = np.repeat(log_counts, len(weights)) + np.tile(
        log_counts, len(weights)
    )
    p1 = first.astype(np.float64) / args.outer_bits
    p2 = second.astype(np.float64) / args.outer_bits
    probabilities = np.stack(
        (
            (1.0 - p1) * (1.0 - p2),
            p1 * (1.0 - p2) + (1.0 - p1) * p2,
            p1 * p2,
        ),
        axis=1,
    )
    with np.errstate(divide="ignore"):
        log_probabilities = np.log(probabilities)
    fixed_weight_log_probabilities = np.empty(len(weights), dtype=np.float64)
    for index, weight in enumerate(weights):
        probability = float(weight) / args.outer_bits
        if weight == args.outer_bits:
            fixed_weight_log_probabilities[index] = 0.0
        else:
            fixed_weight_log_probabilities[index] = (
                log_choose(args.outer_bits, int(weight))
                + weight * math.log(probability)
                + (args.outer_bits - weight) * math.log1p(-probability)
            )
    pair_conditioning_penalty = -(
        np.repeat(fixed_weight_log_probabilities, len(weights))
        + np.tile(fixed_weight_log_probabilities, len(weights))
    )

    nonactivation = load_uniform_nonactivation(args.activation, args.step_bits)
    live_spectrum = load_nonzero_spectrum(
        args.live_spectrum, args.step_bits, args.state_bits
    )
    best_inner = np.full(len(first), math.inf, dtype=np.float64)
    best_tilt = np.full(len(first), math.nan, dtype=np.float64)
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
            candidate_logs = np.log(impulses)
        region_logs = regular_region_log_matrices(
            candidate_logs,
            args.step_bits,
            epochs,
            2,
        )
        mixed = np.full((len(first), 2, 2), -math.inf, dtype=np.float64)
        for candidates in range(3):
            mixed = np.logaddexp(
                mixed,
                region_logs[candidates][None, :, :]
                + log_probabilities[:, candidates, None, None],
            )
        powered = log_matrix_power_pairwise(mixed, args.outer_bits)
        moments = np.logaddexp(powered[:, 0, 0], powered[:, 0, 1])
        bounds = moments + distance * surprisal
        improved = bounds < best_inner
        best_inner[improved] = bounds[improved]
        best_tilt[improved] = log_surprisal
        print(f"tilt,{log_surprisal:.6f}", flush=True)

    location_log = log_choose(outer_blocks, 2)
    contributions = (
        location_log
        + pair_log_counts
        + np.minimum(0.0, best_inner)
        + pair_conditioning_penalty
    )
    aggregate = float(logsumexp(contributions))
    dominant_index = int(np.argmax(contributions))
    return {
        "schema": "riffle-two-tail-conditioned-iid-regions-v2",
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
            "tail_weight_count": int(len(weights)),
            "ordered_weight_pairs": int(len(first)),
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "constituent_distance": args.constituent_distance,
            "epochs_per_region": epochs,
            "log_surprisals": list(args.log_surprisals),
        },
        "aggregate_log2_upper": aggregate / LOG2,
        "aggregate_margin_bits": -aggregate / LOG2,
        "dominant_pair": {
            "first_weight": int(first[dominant_index]),
            "second_weight": int(second[dominant_index]),
            "pair_spectrum_log2": float(pair_log_counts[dominant_index] / LOG2),
            "inner_log2_upper": float(
                min(0.0, best_inner[dominant_index]) / LOG2
            ),
            "conditioning_penalty_bits": float(
                pair_conditioning_penalty[dominant_index] / LOG2
            ),
            "best_log_surprisal": float(best_tilt[dominant_index]),
            "pointwise_log2_upper": float(contributions[dominant_index] / LOG2),
            "pointwise_margin_bits": float(-contributions[dominant_index] / LOG2),
        },
        "scope": (
            "Complete ordered two-tail outer-weight sum. The iid-region moment "
            "is converted to a fixed-support upper bound by dividing by the "
            "exact binomial probability of each conditioned outer weight. The "
            "region transfer is the recorded RM2Sub pre-add-multiply model. "
            "Nearest-binary64 and not outward rounded."
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
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0),
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
    dominant = result["dominant_pair"]
    print(
        "dominant_pair,"
        f"{dominant['first_weight']},{dominant['second_weight']},"
        f"inner_log2,{dominant['inner_log2_upper']:.9f},"
        f"conditioning_bits,{dominant['conditioning_penalty_bits']:.9f},"
        f"margin,{dominant['pointwise_margin_bits']:.9f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
