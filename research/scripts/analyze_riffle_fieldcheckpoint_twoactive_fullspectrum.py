#!/usr/bin/env python3
"""Complete two-active modeled-spectrum diagnostic for FieldCheckpointAccumulate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_fieldcheckpoint_twoactive_w38 import (
    category_log_moments,
    region_matrices_up_to_two,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2, log_choose


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal05_two_active_fullspectrum.json"
)


def intersection_terms(
    outer_bits: int,
    first_weight: int,
    second_weight: int,
    category_logs: np.ndarray,
) -> tuple[float, int, float]:
    lower = max(0, first_weight + second_weight - outer_bits)
    upper = min(first_weight, second_weight)
    denominator = log_choose(outer_bits, second_weight)
    terms = []
    dominant_intersection = lower
    dominant_term = -math.inf
    for intersection in range(lower, upper + 1):
        singletons = first_weight + second_weight - 2 * intersection
        log_probability = (
            log_choose(first_weight, intersection)
            + log_choose(
                outer_bits - first_weight,
                second_weight - intersection,
            )
            - denominator
        )
        log_moment = float(category_logs[singletons, intersection])
        term = log_probability + log_moment
        terms.append(term)
        if term > dominant_term:
            dominant_term = term
            dominant_intersection = intersection
    return float(logsumexp(np.asarray(terms))), dominant_intersection, dominant_term


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    weights = [
        weight
        for weight in range(1, args.outer_bits + 1)
        if math.isfinite(float(spectrum[weight]))
    ]
    pairs = [
        (first_weight, second_weight)
        for first_index, first_weight in enumerate(weights)
        for second_weight in weights[first_index:]
    ]
    best = {pair: math.inf for pair in pairs}
    best_tilt = {pair: math.nan for pair in pairs}
    best_intersection = {pair: -1 for pair in pairs}
    best_intersection_term = {pair: -math.inf for pair in pairs}

    grid_count = int(
        math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
    ) + 1
    for grid_index in range(grid_count):
        log_surprisal = args.grid_min + grid_index * args.grid_step
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        regions = region_matrices_up_to_two(
            args.state_bits,
            args.epoch_bits,
            args.epochs_per_region,
            z,
        )
        categories = category_log_moments(
            regions,
            args.outer_bits,
            args.outer_bits,
            args.outer_bits,
        )
        for pair in pairs:
            log_moment, intersection, intersection_term = intersection_terms(
                args.outer_bits, pair[0], pair[1], categories
            )
            candidate = log_moment + distance * surprisal
            if candidate < best[pair]:
                best[pair] = candidate
                best_tilt[pair] = log_surprisal
                best_intersection[pair] = intersection
                best_intersection_term[pair] = intersection_term
        print(
            f"grid,{grid_index + 1},{grid_count},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    block_pair_log = math.log(math.comb(outer_blocks, 2))
    rows = []
    contributions = []
    for pair in pairs:
        first_weight, second_weight = pair
        assignment_log = 0.0 if first_weight == second_weight else math.log(2.0)
        outer_log = (
            block_pair_log
            + assignment_log
            + float(spectrum[first_weight])
            + float(spectrum[second_weight])
        )
        inner_log = min(0.0, best[pair])
        contribution = outer_log + inner_log
        contributions.append(contribution)
        rows.append(
            {
                "first_outer_weight": first_weight,
                "second_outer_weight": second_weight,
                "best_log_surprisal": best_tilt[pair],
                "inner_log2_upper": inner_log / LOG2,
                "outer_log2_multiplicity": outer_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
                "dominant_intersection": best_intersection[pair],
                "dominant_intersection_summand_log2": (
                    best_intersection_term[pair] / LOG2
                ),
            }
        )
    total_log = float(logsumexp(np.asarray(contributions)))
    dominant = sorted(
        rows, key=lambda row: float(row["pointwise_log2_upper"]), reverse=True
    )[:20]
    return {
        "schema": "riffle-fieldcheckpoint-two-active-fullspectrum-v1",
        "candidate": "Riffle FieldCheckpointAccumulate t=32 s=64 K=32",
        "probability_space": {
            "outer_spectrum": (
                "modeled real-valued complement-symmetric even "
                f"[{args.outer_bits},{dimension},"
                f"{args.modeled_minimum_distance}]-shaped spectrum"
            ),
            "active_outer_blocks": "exactly two distinct outer-block positions",
            "coordinate_permutations": "independent by outer block",
            "region_permutations": "independent by region",
            "checkpoint_multipliers": "independent uniform nonzero GF(2^64) elements",
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "epochs_per_region": args.epochs_per_region,
        },
        "weight_pair_count": len(rows),
        "two_active_log2_upper": total_log / LOG2,
        "two_active_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_weight_pairs": dominant,
        "weight_pair_rows": rows,
        "scope": (
            "Complete floating-point Chernoff diagnostic for exactly two active "
            "outer blocks under the modeled spectrum. Direct binary64 "
            "coefficient scaling can underflow negligible categories; the final "
            "certificate must replace it with outward-rounded scaled arithmetic."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=1024)
    parser.add_argument("--epochs-per-region", type=int, default=8)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--grid-min", type=float, default=-8.5)
    parser.add_argument("--grid-max", type=float, default=-6.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "two_active_fullspectrum_lambda_bits,"
        f"{payload['two_active_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
