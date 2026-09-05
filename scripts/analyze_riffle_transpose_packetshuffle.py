#!/usr/bin/env python3
"""Packed-group diagnostic for TransposePacketShuffle-RandomStepConv.

For each active-outer-block count, the evaluator packs the active blocks into
the fewest fixed groups of size g.  It computes the exact packet-order average
for one transposed row by a log-domain matrix-polynomial dynamic program, then
composes B independent rows. The candidate proof folder proves packing
domination. The resulting number remains provisional until numerical
rounding is certified.
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


DEFAULT_MESSAGE_BITS = 1 << 20
DEFAULT_OUTPUT = Path(
    "constructions/riffle_transpose_packetshuffle_randomstepconv/"
    "receipts/initial_g4_b1024_sigma15.json"
)


def transition_by_active_rank(
    packet_bits: int, sigma: int, z: float, active_rank: int
) -> tuple[float, float, float, float]:
    """Return the averaged transition moment for one packet group."""
    if not 0 <= active_rank <= packet_bits:
        raise ValueError("active rank must lie in [0,g]")
    state_zero = math.ldexp(1.0, -sigma)
    output_moment = math.ldexp((1.0 + z) ** packet_bits, -packet_bits)
    off = state_zero * output_moment
    live = (1.0 - state_zero) * output_moment
    zero_input = (1.0, 0.0, off, live)
    if active_rank == 0:
        return zero_input
    nonzero_input = (off, live, off, live)
    packet_zero = math.ldexp(1.0, -active_rank)
    return tuple(
        packet_zero * zero_value + (1.0 - packet_zero) * nonzero_value
        for zero_value, nonzero_value in zip(zero_input, nonzero_input)
    )


def packed_coefficient_logs(
    *, packet_bits: int, sigma: int, positions: int, z: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return coefficients for full groups and for one partial group.

    The first result has shape (P+1,2,2) and stores the coefficient with f
    full groups.  The second has shape (g-1,P+1,2,2) and stores the
    coefficient with f full groups and exactly one group of rank r+1.
    """
    log_zero = log_matrix_entries(
        transition_by_active_rank(packet_bits, sigma, z, 0)
    )
    log_full = log_matrix_entries(
        transition_by_active_rank(packet_bits, sigma, z, packet_bits)
    )
    log_partial = np.stack(
        [
            log_matrix_entries(
                transition_by_active_rank(packet_bits, sigma, z, rank)
            )
            for rank in range(1, packet_bits)
        ]
    )

    no_partial = np.full((positions + 1, 2, 2), -math.inf)
    no_partial[0, 0, 0] = 0.0
    no_partial[0, 1, 1] = 0.0
    one_partial = np.full(
        (packet_bits - 1, positions + 1, 2, 2), -math.inf
    )

    for index in range(positions):
        degree_count = index + 1
        active_zero = no_partial[:degree_count]
        next_zero = log_matmul_batch(active_zero, log_zero)
        next_full = log_matmul_batch(active_zero, log_full)

        updated_zero = np.full_like(no_partial, -math.inf)
        updated_zero[:degree_count] = next_zero
        updated_zero[1 : degree_count + 1] = np.logaddexp(
            updated_zero[1 : degree_count + 1], next_full
        )

        updated_partial = np.full_like(one_partial, -math.inf)
        if packet_bits > 1:
            old = one_partial[:, :degree_count].reshape(-1, 2, 2)
            old_zero = log_matmul_batch(old, log_zero).reshape(
                packet_bits - 1, degree_count, 2, 2
            )
            old_full = log_matmul_batch(old, log_full).reshape(
                packet_bits - 1, degree_count, 2, 2
            )
            updated_partial[:, :degree_count] = old_zero
            updated_partial[:, 1 : degree_count + 1] = np.logaddexp(
                updated_partial[:, 1 : degree_count + 1], old_full
            )
            for rank_index in range(packet_bits - 1):
                inserted = log_matmul_batch(
                    active_zero, log_partial[rank_index]
                )
                updated_partial[rank_index, :degree_count] = np.logaddexp(
                    updated_partial[rank_index, :degree_count], inserted
                )

        no_partial = updated_zero
        one_partial = updated_partial

    return no_partial, one_partial


def normalized_packed_row_logs(
    no_partial: np.ndarray,
    one_partial: np.ndarray,
    *,
    active_blocks: int,
    packet_bits: int,
    positions: int,
) -> np.ndarray:
    """Return the log entries of the averaged packed-row transition."""
    full_groups, remainder = divmod(active_blocks, packet_bits)
    if remainder == 0:
        coefficient = no_partial[full_groups]
        normalizer = log_choose(positions, full_groups)
    else:
        coefficient = one_partial[remainder - 1, full_groups]
        normalizer = log_choose(positions, full_groups) + math.log(
            positions - full_groups
        )
    return coefficient - normalizer


def numeric_permutation_average(
    matrices: list[np.ndarray], labels: tuple[int, ...]
) -> np.ndarray:
    """Brute-force a distinct-type packet-order average for a self-test."""
    orders = sorted(set(itertools.permutations(labels)))
    total = np.zeros((2, 2))
    for order in orders:
        product = np.eye(2)
        for label in order:
            product = product @ matrices[label]
        total += product
    return total / len(orders)


def self_tests() -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for packet_bits, sigma, positions, z in (
        (2, 3, 4, 0.37),
        (4, 5, 5, 0.61),
    ):
        no_partial, one_partial = packed_coefficient_logs(
            packet_bits=packet_bits,
            sigma=sigma,
            positions=positions,
            z=z,
        )
        matrices = [
            np.asarray(
                transition_by_active_rank(packet_bits, sigma, z, rank)
            ).reshape(2, 2)
            for rank in range(packet_bits + 1)
        ]
        maximum_error = 0.0
        for active_blocks in range(1, positions * packet_bits + 1):
            full_groups, remainder = divmod(active_blocks, packet_bits)
            if remainder:
                labels = (
                    (packet_bits,) * full_groups
                    + (remainder,)
                    + (0,) * (positions - full_groups - 1)
                )
            else:
                labels = (
                    (packet_bits,) * full_groups
                    + (0,) * (positions - full_groups)
                )
            brute = numeric_permutation_average(matrices, labels)
            log_average = normalized_packed_row_logs(
                no_partial,
                one_partial,
                active_blocks=active_blocks,
                packet_bits=packet_bits,
                positions=positions,
            )
            recovered = np.exp(log_average)
            maximum_error = max(
                maximum_error, float(np.max(np.abs(brute - recovered)))
            )
        if maximum_error > 2e-12:
            raise AssertionError("packed coefficient DP changed the row average")
        rows.append(
            {
                "packet_bits": packet_bits,
                "positions": positions,
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
    if outer_bits % 2:
        raise ValueError("outer block length must be even")
    outer_dimension = outer_bits // 2
    if message_bits % outer_dimension:
        raise ValueError("outer dimension must divide the message length")
    blocks = message_bits // outer_dimension
    if blocks % packet_bits:
        raise ValueError("packet width must divide the outer-block count")
    packet_groups = blocks // packet_bits
    packet_positions = outer_bits * packet_groups
    distance = math.floor(relative_distance * packet_bits * packet_positions)

    best_inner = np.full(blocks + 1, math.inf)
    best_surprisal = np.full(blocks + 1, math.nan)
    grid_count = int(math.floor((grid_max - grid_min) / grid_step + 0.5)) + 1
    grid = [grid_min + index * grid_step for index in range(grid_count)]

    for grid_index, log_surprisal in enumerate(grid):
        z = math.exp(-math.exp(log_surprisal))
        no_partial, one_partial = packed_coefficient_logs(
            packet_bits=packet_bits,
            sigma=sigma,
            positions=packet_groups,
            z=z,
        )
        correction = distance * math.exp(log_surprisal)
        for active_blocks in range(1, blocks + 1):
            entry_logs = normalized_packed_row_logs(
                no_partial,
                one_partial,
                active_blocks=active_blocks,
                packet_bits=packet_bits,
                positions=packet_groups,
            ).reshape(4)
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

    log_messages = log_two_power_minus_one(outer_dimension)
    conditioning_penalty = -math.log1p(-math.ldexp(1.0, -outer_bits))
    rows: list[dict[str, object]] = []
    for active_blocks in range(1, blocks + 1):
        outer_log = log_choose(blocks, active_blocks) + active_blocks * (
            log_messages + conditioning_penalty
        )
        inner_log2 = min(0.0, float(best_inner[active_blocks]) / LOG2)
        full_groups, remainder = divmod(active_blocks, packet_bits)
        rows.append(
            {
                "active_outer_blocks": active_blocks,
                "packed_full_groups": full_groups,
                "packed_remainder_rank": remainder,
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
    dominant_log_surprisal = float(dominant["log_output_surprisal"])
    dominant_z = math.exp(-math.exp(dominant_log_surprisal))
    no_partial, one_partial = packed_coefficient_logs(
        packet_bits=packet_bits,
        sigma=sigma,
        positions=packet_groups,
        z=dominant_z,
    )
    dominant_logs = normalized_packed_row_logs(
        no_partial,
        one_partial,
        active_blocks=dominant_active,
        packet_bits=packet_bits,
        positions=packet_groups,
    ).reshape(4)
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
        "fixed_block_groups": packet_groups,
        "transposed_rows": outer_bits,
        "packet_positions": packet_positions,
        "sigma": sigma,
        "distance": distance,
        "coefficient_method": "exact_packed_group_row_log_dp_output_grid",
        "packing_domination_status": "proved_in_proof/PACKING_DOMINATION.md",
        "output_grid": {
            "minimum": grid_min,
            "maximum": grid_max,
            "step": grid_step,
            "count": grid_count,
        },
        "log2_expected_bad_upper_if_packing_lemma": total_log2,
        "lambda_bits_if_packing_lemma": -total_log2,
        "lambda_bits_lower_float": -total_log2,
        "dominant_profile": dominant,
        "occupation_rows": rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--message-bits", type=int, default=DEFAULT_MESSAGE_BITS)
    parser.add_argument("--outer-bits", type=int, default=1024)
    parser.add_argument("--g", type=int, default=4)
    parser.add_argument("--sigma", type=int, default=15)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=2.0)
    parser.add_argument("--grid-step", type=float, default=0.1)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    checks = self_tests()
    if args.self_test_only:
        print(json.dumps({"packed_row_checks": checks}, indent=2))
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
        "schema": "riffle-transpose-packetshuffle-diagnostic-v1",
        "candidate": "Riffle TransposePacketShuffle-RandomStepConv g",
        "message_bits": args.message_bits,
        "relative_distance": args.relative_distance,
        "target_lambda_bits": 40,
        "small_exact_checks": checks,
        "scope": (
            "Exact packed-row packet-order coefficients and complete active-"
            "block sum. Packing domination is proved in the candidate proof "
            "folder. The complete calculation is not outward-rounded."
        ),
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
                "lambda_bits_if_packing_lemma": case[
                    "lambda_bits_if_packing_lemma"
                ],
                "dominant_profile": case["dominant_profile"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(f"wrote,{args.output}", flush=True)


if __name__ == "__main__":
    main()
