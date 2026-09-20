#!/usr/bin/env python3
"""Exact-row diagnostic for TransposeBitShuffle-RandomStepConv.

Condition on a active outer blocks. In each transposed row, their candidate
bits occupy a uniform a-subset of the L bit positions. A packet containing r
candidate bits is zero with probability 2^-r. The row transition is therefore
the degree-a coefficient of a degree-g matrix polynomial raised through L/g
ordered packets, divided by C(L,a).
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_striped_random_outer import (
    LOG2,
    log_choose,
    log_matmul_batch,
    log_matrix_entries,
    log_matrix_power_moment,
    log_two_power_minus_one,
    tilted_transition_profile,
)
from analyze_riffle_transpose_packetshuffle import transition_by_active_rank


DEFAULT_OUTPUT = Path(
    "constructions/riffle_transpose_bitshuffle_randomstepconv/"
    "receipts/initial.json"
)


def bitshuffle_row_coefficient_logs(
    *, packet_bits: int, sigma: int, blocks: int, z: float
) -> np.ndarray:
    if blocks % packet_bits:
        raise ValueError("packet width must divide the outer-block count")
    positions = blocks // packet_bits
    log_matrices = [
        log_matrix_entries(
            transition_by_active_rank(packet_bits, sigma, z, rank)
        )
        for rank in range(packet_bits + 1)
    ]
    log_multiplicities = [log_choose(packet_bits, rank) for rank in range(packet_bits + 1)]

    coefficients = np.full((blocks + 1, 2, 2), -math.inf)
    coefficients[0, 0, 0] = 0.0
    coefficients[0, 1, 1] = 0.0
    for position in range(positions):
        old_count = position * packet_bits + 1
        active = coefficients[:old_count]
        updated = np.full_like(coefficients, -math.inf)
        for rank in range(packet_bits + 1):
            product = log_matmul_batch(active, log_matrices[rank])
            product += log_multiplicities[rank]
            updated[rank : rank + old_count] = np.logaddexp(
                updated[rank : rank + old_count], product
            )
        coefficients = updated
    return coefficients


def numeric_subset_average(
    *, packet_bits: int, sigma: int, blocks: int, active_blocks: int, z: float
) -> np.ndarray:
    matrices = [
        np.asarray(
            transition_by_active_rank(packet_bits, sigma, z, rank)
        ).reshape(2, 2)
        for rank in range(packet_bits + 1)
    ]
    total = np.zeros((2, 2))
    count = 0
    for subset in itertools.combinations(range(blocks), active_blocks):
        selected = set(subset)
        product = np.eye(2)
        for start in range(0, blocks, packet_bits):
            rank = sum(index in selected for index in range(start, start + packet_bits))
            product = product @ matrices[rank]
        total += product
        count += 1
    return total / count


def self_tests() -> list[dict[str, float | int]]:
    rows = []
    for packet_bits, blocks, sigma, z in ((2, 6, 3, 0.37), (3, 6, 4, 0.61)):
        coefficients = bitshuffle_row_coefficient_logs(
            packet_bits=packet_bits, sigma=sigma, blocks=blocks, z=z
        )
        maximum_error = 0.0
        for active_blocks in range(blocks + 1):
            recovered = np.exp(
                coefficients[active_blocks] - log_choose(blocks, active_blocks)
            )
            brute = numeric_subset_average(
                packet_bits=packet_bits,
                sigma=sigma,
                blocks=blocks,
                active_blocks=active_blocks,
                z=z,
            )
            maximum_error = max(
                maximum_error, float(np.max(np.abs(recovered - brute)))
            )
        if maximum_error > 3e-12:
            raise AssertionError("bit-shuffle coefficient DP changed the row average")
        rows.append(
            {
                "packet_bits": packet_bits,
                "blocks": blocks,
                "maximum_error": maximum_error,
            }
        )
    return rows


def evaluate_case(
    *,
    message_bits: int,
    outer_bits: int,
    packet_bits: int,
    sigma: int,
    relative_distance: float,
    grid_min: float,
    grid_max: float,
    grid_step: float,
) -> dict[str, object]:
    outer_dimension = outer_bits // 2
    if message_bits % outer_dimension:
        raise ValueError("outer dimension must divide the message length")
    blocks = message_bits // outer_dimension
    if blocks % packet_bits:
        raise ValueError("packet width must divide the outer-block count")
    packet_positions = outer_bits * blocks // packet_bits
    output_bits = packet_bits * packet_positions
    distance = math.floor(relative_distance * output_bits)
    log_combinations = np.asarray(
        [log_choose(blocks, active) for active in range(blocks + 1)]
    )
    best_inner = np.full(blocks + 1, math.inf)
    best_surprisal = np.full(blocks + 1, math.nan)

    grid_count = int(math.floor((grid_max - grid_min) / grid_step + 0.5)) + 1
    for grid_index in range(grid_count):
        log_surprisal = grid_min + grid_index * grid_step
        z = math.exp(-math.exp(log_surprisal))
        coefficients = bitshuffle_row_coefficient_logs(
            packet_bits=packet_bits,
            sigma=sigma,
            blocks=blocks,
            z=z,
        )
        correction = distance * math.exp(log_surprisal)
        for active_blocks in range(1, blocks + 1):
            entry_logs = (
                coefficients[active_blocks].reshape(4)
                - log_combinations[active_blocks]
            )
            scale = float(np.max(entry_logs))
            row_matrix = tuple(
                math.exp(float(entry - scale)) for entry in entry_logs
            )
            moment = (
                log_matrix_power_moment(row_matrix, outer_bits)
                + outer_bits * scale
                + correction
            )
            if moment < best_inner[active_blocks]:
                best_inner[active_blocks] = moment
                best_surprisal[active_blocks] = log_surprisal
        if grid_index % max(1, grid_count // 10) == 0:
            print(
                f"grid,{grid_index + 1},{grid_count},"
                f"log_surprisal,{log_surprisal:.6f}",
                flush=True,
            )

    messages_per_block = log_two_power_minus_one(outer_dimension)
    conditioning_penalty = -math.log1p(-math.ldexp(1.0, -outer_bits))
    rows: list[dict[str, object]] = []
    for active_blocks in range(1, blocks + 1):
        outer_log = log_combinations[active_blocks] + active_blocks * (
            messages_per_block + conditioning_penalty
        )
        inner_log2 = min(0.0, float(best_inner[active_blocks]) / LOG2)
        rows.append(
            {
                "active_outer_blocks": active_blocks,
                "outer_message_log2": outer_log / LOG2,
                "inner_log2_upper": inner_log2,
                "pointwise_log2_upper": outer_log / LOG2 + inner_log2,
                "log_output_surprisal": float(best_surprisal[active_blocks]),
            }
        )

    pointwise = np.asarray(
        [float(row["pointwise_log2_upper"]) for row in rows]
    )
    total_log2 = float(logsumexp(pointwise * LOG2) / LOG2)
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    dominant_active = int(dominant["active_outer_blocks"])
    dominant_u = float(dominant["log_output_surprisal"])
    dominant_z = math.exp(-math.exp(dominant_u))
    dominant_coefficients = bitshuffle_row_coefficient_logs(
        packet_bits=packet_bits,
        sigma=sigma,
        blocks=blocks,
        z=dominant_z,
    )
    dominant_logs = (
        dominant_coefficients[dominant_active].reshape(4)
        - log_combinations[dominant_active]
    )
    dominant_scale = float(np.max(dominant_logs))
    dominant_matrix = tuple(
        math.exp(float(entry - dominant_scale)) for entry in dominant_logs
    )
    dominant["tilted_row_transition_profile"] = tilted_transition_profile(
        dominant_matrix, outer_bits
    )

    return {
        "outer_bits": outer_bits,
        "outer_dimension_bits": outer_dimension,
        "outer_blocks": blocks,
        "packet_bits": packet_bits,
        "packet_positions_per_row": blocks // packet_bits,
        "transposed_rows": outer_bits,
        "packet_positions": packet_positions,
        "sigma": sigma,
        "distance": distance,
        "coefficient_method": "exact_bitshuffle_row_log_polynomial_dp",
        "output_grid": {
            "minimum": grid_min,
            "maximum": grid_max,
            "step": grid_step,
            "count": grid_count,
        },
        "log2_expected_bad_upper": total_log2,
        "lambda_bits_lower_float": -total_log2,
        "dominant_profile": dominant,
        "occupation_rows": rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=1024)
    parser.add_argument("--g", type=int, default=4)
    parser.add_argument("--sigma", type=int, default=13)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    checks = self_tests()
    if args.self_test_only:
        print(json.dumps({"bitshuffle_row_checks": checks}, indent=2))
        return
    case = evaluate_case(
        message_bits=args.message_bits,
        outer_bits=args.outer_bits,
        packet_bits=args.g,
        sigma=args.sigma,
        relative_distance=args.relative_distance,
        grid_min=args.grid_min,
        grid_max=args.grid_max,
        grid_step=args.grid_step,
    )
    payload = {
        "schema": "riffle-transpose-bitshuffle-diagnostic-v1",
        "candidate": "Riffle TransposeBitShuffle-RandomStepConv g",
        "message_bits": args.message_bits,
        "relative_distance": args.relative_distance,
        "target_lambda_bits": 40,
        "small_exact_checks": checks,
        "case": case,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "outer_bits": case["outer_bits"],
                "packet_bits": case["packet_bits"],
                "sigma": case["sigma"],
                "lambda_bits_lower_float": case["lambda_bits_lower_float"],
                "dominant_profile": case["dominant_profile"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()

