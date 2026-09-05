#!/usr/bin/env python3
"""Exact floating-point two-active weight-38 diagnostic.

The checker builds R_0(z), R_1(z), and R_2(z) from the exact epoch occupancy
transfers.  It then averages the ordered 256-region product over the exact
intersection law of two independently permuted weight-38 outer words.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_fieldcheckpoint_accumulate_oneblock import (
    occupancy_epoch_matrices,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2, log_choose


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal04_two_active_w38.json"
)


def region_matrices_up_to_two(
    state_bits: int,
    epoch_bits: int,
    epochs_per_region: int,
    z: float,
) -> list[np.ndarray]:
    """Average the ordered epoch product for region weights zero through two."""
    epoch = occupancy_epoch_matrices(state_bits, epoch_bits, z)
    weighted_epoch = [
        math.comb(epoch_bits, weight) * epoch[weight]
        for weight in range(3)
    ]
    coefficients = [np.eye(2), np.zeros((2, 2)), np.zeros((2, 2))]
    for _ in range(epochs_per_region):
        updated = [np.zeros((2, 2)) for _ in range(3)]
        for total_weight in range(3):
            for epoch_weight in range(total_weight + 1):
                updated[total_weight] += (
                    coefficients[total_weight - epoch_weight]
                    @ weighted_epoch[epoch_weight]
                )
        coefficients = updated
    region_bits = epoch_bits * epochs_per_region
    return [
        coefficients[weight] / math.comb(region_bits, weight)
        for weight in range(3)
    ]


def category_log_moments(
    region_matrices: list[np.ndarray],
    outer_bits: int,
    maximum_singletons: int,
    maximum_overlaps: int,
) -> np.ndarray:
    """Average ordered products for every singleton/overlap region count."""
    coefficients = np.zeros(
        (maximum_singletons + 1, maximum_overlaps + 1, 2, 2),
        dtype=np.float64,
    )
    coefficients[0, 0] = np.eye(2)
    inverse_choices = 1.0 / 3.0
    for _ in range(outer_bits):
        updated = np.einsum(
            "...ij,jk->...ik", coefficients, region_matrices[0]
        ) * inverse_choices
        updated[1:] += np.einsum(
            "...ij,jk->...ik", coefficients[:-1], region_matrices[1]
        ) * inverse_choices
        updated[:, 1:] += np.einsum(
            "...ij,jk->...ik", coefficients[:, :-1], region_matrices[2]
        ) * inverse_choices
        coefficients = updated

    log_moments = np.full(
        (maximum_singletons + 1, maximum_overlaps + 1),
        -math.inf,
    )
    for singletons in range(maximum_singletons + 1):
        for overlaps in range(maximum_overlaps + 1):
            empty = outer_bits - singletons - overlaps
            if empty < 0:
                continue
            row_sum = coefficients[singletons, overlaps, 0].sum()
            if row_sum <= 0.0:
                continue
            log_multinomial = (
                math.lgamma(outer_bits + 1)
                - math.lgamma(empty + 1)
                - math.lgamma(singletons + 1)
                - math.lgamma(overlaps + 1)
            )
            log_moments[singletons, overlaps] = (
                math.log(row_sum)
                + outer_bits * math.log(3.0)
                - log_multinomial
            )
    return log_moments


def intersection_log_probabilities(outer_bits: int, weight: int) -> np.ndarray:
    values = np.full(weight + 1, -math.inf)
    denominator = log_choose(outer_bits, weight)
    for intersection in range(weight + 1):
        if weight - intersection > outer_bits - weight:
            continue
        values[intersection] = (
            log_choose(weight, intersection)
            + log_choose(outer_bits - weight, weight - intersection)
            - denominator
        )
    return values


def self_test(args: argparse.Namespace) -> dict[str, float]:
    regions = region_matrices_up_to_two(3, 6, 2, 1.0)
    region_error = max(
        float(np.max(np.abs(matrix.sum(axis=1) - 1.0)))
        for matrix in regions
    )
    categories = category_log_moments(regions, 6, 4, 2)
    category_error = 0.0
    for singletons in range(5):
        for overlaps in range(3):
            if singletons + overlaps <= 6:
                category_error = max(
                    category_error,
                    abs(math.exp(categories[singletons, overlaps]) - 1.0),
                )
    intersection_error = abs(
        float(np.exp(intersection_log_probabilities(8, 3)).sum()) - 1.0
    )
    if max(region_error, category_error, intersection_error) > 3e-13:
        raise AssertionError("two-active transfer normalization failed")
    return {
        "maximum_region_stochastic_error": region_error,
        "maximum_category_stochastic_error": category_error,
        "intersection_mass_error": intersection_error,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    spectrum_log = float(spectrum[args.outer_weight])
    intersection_logs = intersection_log_probabilities(
        args.outer_bits, args.outer_weight
    )

    best_log = math.inf
    best_log_surprisal = math.nan
    best_intersection_rows: list[dict[str, float | int]] = []
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
        category_logs = category_log_moments(
            regions,
            args.outer_bits,
            2 * args.outer_weight,
            args.outer_weight,
        )
        terms = []
        intersection_rows = []
        for intersection in range(args.outer_weight + 1):
            if not math.isfinite(float(intersection_logs[intersection])):
                continue
            singletons = 2 * args.outer_weight - 2 * intersection
            log_moment = float(category_logs[singletons, intersection])
            term = float(intersection_logs[intersection]) + log_moment
            terms.append(term)
            intersection_rows.append(
                {
                    "intersection": intersection,
                    "singleton_regions": singletons,
                    "intersection_log2_probability": (
                        float(intersection_logs[intersection]) / LOG2
                    ),
                    "conditional_log2_moment": log_moment / LOG2,
                    "summand_log2": term / LOG2,
                }
            )
        log_inner_moment = float(logsumexp(np.asarray(terms)))
        candidate = log_inner_moment + distance * surprisal
        if candidate < best_log:
            best_log = candidate
            best_log_surprisal = log_surprisal
            best_intersection_rows = intersection_rows
        print(
            f"grid,{grid_index + 1},{grid_count},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    best_inner = min(0.0, best_log)
    outer_log = (
        math.log(math.comb(outer_blocks, 2)) + 2.0 * spectrum_log
    )
    pointwise = outer_log + best_inner
    dominant_intersection = max(
        best_intersection_rows, key=lambda row: float(row["summand_log2"])
    )
    return {
        "schema": "riffle-fieldcheckpoint-two-active-w38-v1",
        "candidate": "Riffle FieldCheckpointAccumulate t=32 s=64 K=32",
        "probability_space": {
            "outer_words": (
                f"two modeled-spectrum words, each conditioned on weight "
                f"{args.outer_weight}"
            ),
            "coordinate_permutations": "independent for the two outer blocks",
            "region_permutations": "independent for all 256 regions",
            "checkpoint_multipliers": "independent uniform nonzero GF(2^64) elements",
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "outer_weight": args.outer_weight,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "epochs_per_region": args.epochs_per_region,
        },
        "self_test": self_test(args),
        "best_log_surprisal": best_log_surprisal,
        "inner_log2_upper": best_inner / LOG2,
        "modeled_pair_multiplicity_log2": outer_log / LOG2,
        "pointwise_log2_upper": pointwise / LOG2,
        "lambda_bits_lower_for_this_shell": -pointwise / LOG2,
        "dominant_intersection": dominant_intersection,
        "intersection_rows_at_best_tilt": best_intersection_rows,
        "scope": (
            "Complete floating-point Chernoff diagnostic for the shell with "
            "exactly two active outer blocks and both outer weights equal to 38. "
            "The outer spectrum is modeled and the arithmetic is not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--outer-weight", type=int, default=38)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=1024)
    parser.add_argument("--epochs-per-region", type=int, default=8)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-5.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "two_active_w38_lambda_bits,"
        f"{payload['lambda_bits_lower_for_this_shell']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
