#!/usr/bin/env python3
"""Low-occupation regular/all-one envelope for FieldCheckpointAccumulate."""

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
from analyze_riffle_fieldcheckpoint_regular_envelope import (
    log_row_power_moment,
    spectrum_density_envelope_log,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2, log_choose, log_two_power_minus_one


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal07_mixed_envelope_total64.json"
)


def active_region_matrices(
    state_bits: int,
    epoch_bits: int,
    epochs_per_region: int,
    maximum_weight: int,
    z: float,
) -> list[np.ndarray]:
    """Return R_h(z) for uniform active supports through maximum_weight."""
    epoch = occupancy_epoch_matrices(state_bits, epoch_bits, z)
    weighted_epoch = [
        math.comb(epoch_bits, weight) * epoch[weight]
        for weight in range(maximum_weight + 1)
    ]
    coefficients = [np.zeros((2, 2)) for _ in range(maximum_weight + 1)]
    coefficients[0] = np.eye(2)
    for _ in range(epochs_per_region):
        updated = [np.zeros((2, 2)) for _ in range(maximum_weight + 1)]
        for total in range(maximum_weight + 1):
            for current in range(total + 1):
                updated[total] += coefficients[total - current] @ weighted_epoch[current]
        coefficients = updated
    region_bits = epoch_bits * epochs_per_region
    return [
        coefficients[weight] / math.comb(region_bits, weight)
        for weight in range(maximum_weight + 1)
    ]


def mixed_region_matrix(
    active_regions: list[np.ndarray], regular_blocks: int, all_one_blocks: int
) -> np.ndarray:
    """Average b forced impulses and a independently fair regular bits."""
    result = np.zeros((2, 2), dtype=np.float64)
    scale = math.ldexp(1.0, -regular_blocks)
    for regular_ones in range(regular_blocks + 1):
        result += (
            math.comb(regular_blocks, regular_ones)
            * scale
            * active_regions[all_one_blocks + regular_ones]
        )
    return result


def self_test() -> dict[str, float]:
    active = active_region_matrices(3, 6, 2, 6, 1.0)
    maximum_error = 0.0
    for regular in range(4):
        for all_one in range(4 - regular):
            matrix = mixed_region_matrix(active, regular, all_one)
            maximum_error = max(
                maximum_error,
                float(np.max(np.abs(matrix.sum(axis=1) - 1.0))),
            )
    if maximum_error > 4e-13:
        raise AssertionError("mixed region transfer is not stochastic")
    return {"maximum_small_stochastic_error": maximum_error}


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta

    pairs = [
        (regular, all_one)
        for regular in range(args.maximum_total_blocks + 1)
        for all_one in range(args.maximum_total_blocks - regular + 1)
        if regular + all_one > 0
    ]
    best = {pair: math.inf for pair in pairs}
    best_tilt = {pair: math.nan for pair in pairs}
    grid_count = int(
        math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
    ) + 1
    for grid_index in range(grid_count):
        log_surprisal = args.grid_min + grid_index * args.grid_step
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        active = active_region_matrices(
            args.state_bits,
            args.epoch_bits,
            args.epochs_per_region,
            args.maximum_total_blocks,
            z,
        )
        for pair in pairs:
            matrix = mixed_region_matrix(active, pair[0], pair[1])
            log_moment = log_row_power_moment(matrix, args.outer_bits)
            candidate = log_moment + distance * surprisal
            if candidate < best[pair]:
                best[pair] = candidate
                best_tilt[pair] = log_surprisal
        print(
            f"grid,{grid_index + 1},{grid_count},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    contributions = []
    for regular, all_one in pairs:
        outer_log = (
            log_choose(outer_blocks, regular)
            + log_choose(outer_blocks - regular, all_one)
            + regular * regular_log_mass
        )
        inner_log = min(0.0, best[(regular, all_one)])
        contribution = outer_log + inner_log
        contributions.append(contribution)
        rows.append(
            {
                "regular_active_blocks": regular,
                "all_one_active_blocks": all_one,
                "best_log_surprisal": best_tilt[(regular, all_one)],
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
            }
        )
    total_log = float(logsumexp(np.asarray(contributions)))
    dominant = sorted(
        rows, key=lambda row: float(row["pointwise_log2_upper"]), reverse=True
    )[:30]
    all_one_dominant = sorted(
        (row for row in rows if int(row["all_one_active_blocks"]) > 0),
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]
    return {
        "schema": "riffle-fieldcheckpoint-mixed-envelope-v1",
        "candidate": "Riffle FieldCheckpointAccumulate t=32 s=64 K=32",
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "regular_supports": "pointwise random-support density envelope",
            "all_one_support": "retained exactly with multiplicity one per block",
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
            "maximum_total_blocks": args.maximum_total_blocks,
        },
        "self_test": self_test(),
        "mixed_partial_log2_upper": total_log / LOG2,
        "mixed_partial_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_rows": dominant,
        "dominant_rows_with_all_one": all_one_dominant,
        "occupation_rows": rows,
        "scope": (
            "Floating-point envelope diagnostic for every regular/all-one "
            f"occupation with total active blocks at most {args.maximum_total_blocks}. "
            "Larger occupations and outward rounding remain open."
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
    parser.add_argument("--maximum-total-blocks", type=int, default=64)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=0.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "mixed_partial_lambda_bits,"
        f"{payload['mixed_partial_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
