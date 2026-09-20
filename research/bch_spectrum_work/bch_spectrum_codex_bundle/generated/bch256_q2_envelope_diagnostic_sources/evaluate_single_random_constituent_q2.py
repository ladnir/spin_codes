#!/usr/bin/env python3
"""Exact-shell Q=2 diagnostic for one repeated random constituent.

Sample one uniform linear injection G:F_2^K -> F_2^B and reuse G in all
outer rows.  This program separates two active local messages by their rank.
Equal nonzero messages have rank one.  Distinct nonzero messages have rank
two over F_2.  Independent row-coordinate permutations make the routed
supports independent conditional on their two output weights.

The program computes every conditional support-pair moment.  Nearest
binary64 log arithmetic makes the result diagnostic, not outward certified.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = WORKSTREAM / "single_random_constituent_B256_q2_s22_d109.json"
LOG2 = math.log(2.0)


def log_two_power_minus_one(exponent: int) -> float:
    return exponent * LOG2 + math.log1p(-math.ldexp(1.0, -exponent))


def log_two_power_minus_two(exponent: int) -> float:
    return exponent * LOG2 + math.log1p(-math.ldexp(2.0, -exponent))


def log_choose_rows(length: int) -> np.ndarray:
    return np.asarray(
        [
            math.lgamma(length + 1)
            - math.lgamma(weight + 1)
            - math.lgamma(length - weight + 1)
            for weight in range(length + 1)
        ]
    )


def pair_support_coefficients(
    region_matrices: np.ndarray,
    block_bits: int,
) -> np.ndarray:
    """Average ordered products for two independent fixed-weight supports.

    Entry ``result[a,b]`` averages over an independent uniform a-subset and
    b-subset of the B outer coordinates.  ``region_matrices[j]`` is the
    transfer matrix for a region containing j active input bits.
    """

    current = np.full((block_bits + 1, block_bits + 1, 2, 2), -math.inf)
    current[0, 0] = transfer.log_identity()
    for completed in range(block_bits):
        size = completed + 1
        old = current[:size, :size]
        updated = np.full((size + 1, size + 1, 2, 2), -math.inf)

        unselected = np.log(
            (size - np.arange(size, dtype=np.float64)) / float(size)
        )
        selected = np.log(
            np.arange(1, size + 1, dtype=np.float64) / float(size)
        )

        term = transfer.log_matmul(old, region_matrices[0])
        term += unselected[:, None, None, None]
        term += unselected[None, :, None, None]
        updated[:size, :size] = term

        term = transfer.log_matmul(old, region_matrices[1])
        term += selected[:, None, None, None]
        term += unselected[None, :, None, None]
        updated[1:, :size] = np.logaddexp(updated[1:, :size], term)

        term = transfer.log_matmul(old, region_matrices[1])
        term += unselected[:, None, None, None]
        term += selected[None, :, None, None]
        updated[:size, 1:] = np.logaddexp(updated[:size, 1:], term)

        term = transfer.log_matmul(old, region_matrices[2])
        term += selected[:, None, None, None]
        term += selected[None, :, None, None]
        updated[1:, 1:] = np.logaddexp(updated[1:, 1:], term)
        current[: size + 1, : size + 1] = updated
    return current


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    block_bits = args.block_bits
    dimension = args.dimension
    output_bits = args.output_bits
    if not 0 < dimension < block_bits:
        raise ValueError("dimension must lie strictly between zero and block length")
    if output_bits % block_bits:
        raise ValueError("block length must divide the output length")
    outer_rows = output_bits // block_bits
    distance_cutoff = (
        args.distance_numerator * output_bits
        + args.distance_denominator
        - 1
    ) // args.distance_denominator

    coarse = np.arange(
        args.grid_min,
        args.grid_max + args.grid_step / 2.0,
        args.grid_step,
    )
    fine = np.arange(
        args.fine_min,
        args.fine_max + args.fine_step / 2.0,
        args.fine_step,
    )
    u_values = sorted(
        set(float(value) for value in coarse)
        | set(float(value) for value in fine)
    )

    best = np.full((block_bits + 1, block_bits + 1), math.inf)
    witnesses = np.full((block_bits + 1, block_bits + 1), math.nan)
    for index, u in enumerate(u_values):
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        region_matrices = uniform_coefficients(
            transfer.log_entries(zero),
            transfer.log_entries(active),
            outer_rows,
            2,
        )
        coefficients = pair_support_coefficients(region_matrices, block_bits)
        moments = np.logaddexp(coefficients[..., 0, 0], coefficients[..., 0, 1])
        values = np.minimum(0.0, moments + distance_cutoff * surprisal)
        improved = values < best
        best[improved] = values[improved]
        witnesses[improved] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    log_binomial = log_choose_rows(block_bits)
    nonzero = slice(1, block_bits + 1)

    message_nonzero = log_two_power_minus_one(dimension)
    ambient_nonzero = log_two_power_minus_one(block_bits)
    row_pair_count = (
        math.lgamma(outer_rows + 1)
        - math.lgamma(3)
        - math.lgamma(outer_rows - 1)
    )

    diagonal = np.arange(1, block_bits + 1)
    rank_one_shells = (
        row_pair_count
        + message_nonzero
        - ambient_nonzero
        + log_binomial[nonzero]
        + best[diagonal, diagonal]
    )
    rank_one_aggregate = float(logsumexp(rank_one_shells))

    rank_two_image_logs = (
        log_binomial[:, None] + log_binomial[None, :]
    )
    same_weight = diagonal
    diagonal_correction = np.full(block_bits, -math.inf)
    nontrivial_shells = log_binomial[same_weight] > 0.0
    diagonal_correction[nontrivial_shells] = np.log1p(
        -np.exp(-log_binomial[same_weight][nontrivial_shells])
    )
    rank_two_image_logs[same_weight, same_weight] += diagonal_correction

    image_pair_mass = float(logsumexp(rank_two_image_logs[nonzero, nonzero]))
    expected_image_pair_mass = ambient_nonzero + log_two_power_minus_two(
        block_bits
    )
    if abs(image_pair_mass - expected_image_pair_mass) > 1e-10:
        raise ArithmeticError("rank-two image-pair mass check failed")
    message_pair_mass = float(
        np.logaddexp(
            message_nonzero,
            message_nonzero + log_two_power_minus_two(dimension),
        )
    )
    expected_message_pair_mass = 2.0 * message_nonzero
    if abs(message_pair_mass - expected_message_pair_mass) > 1e-10:
        raise ArithmeticError("rank split does not cover every message pair")
    rank_two_shells = (
        row_pair_count
        + message_nonzero
        + log_two_power_minus_two(dimension)
        - ambient_nonzero
        - log_two_power_minus_two(block_bits)
        + rank_two_image_logs[nonzero, nonzero]
        + best[nonzero, nonzero]
    )
    rank_two_aggregate = float(logsumexp(rank_two_shells))
    total = float(np.logaddexp(rank_one_aggregate, rank_two_aggregate))

    rank_one_top_indices = np.argsort(rank_one_shells)[-10:][::-1]
    flat_rank_two = rank_two_shells.ravel()
    rank_two_top_indices = np.argsort(flat_rank_two)[-20:][::-1]
    rank_two_width = block_bits

    return {
        "schema": "single-random-constituent-randomstepconv-q2-v1",
        "status": "BINARY64_DIAGNOSTIC_EXACT_Q2_SHELL_LAW",
        "claim": {
            "q2_log2_expected_bad_upper": total / LOG2,
            "q2_margin_bits": -total / LOG2,
            "closes_40_bits": total < -40.0 * LOG2,
            "rank_one_log2_expected_bad_upper": rank_one_aggregate / LOG2,
            "rank_one_margin_bits": -rank_one_aggregate / LOG2,
            "rank_two_log2_expected_bad_upper": rank_two_aggregate / LOG2,
            "rank_two_margin_bits": -rank_two_aggregate / LOG2,
        },
        "parameters": {
            "outer_code": f"one random [{block_bits},{dimension}] linear constituent",
            "outer_rows": outer_rows,
            "message_bits": dimension * outer_rows,
            "output_bits": output_bits,
            "distance_cutoff": distance_cutoff,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "memory_bits": args.memory_bits,
            "occupation": 2,
        },
        "probability_space": (
            "sample one uniform linear injection and reuse it in every outer row; "
            "independently sample each row-coordinate permutation, each region "
            "permutation, and every RandomStepConv map"
        ),
        "rank_split": {
            "rank_one": (
                "the two nonzero local messages are equal; their common image is "
                "uniform over the nonzero ambient words"
            ),
            "rank_two": (
                "the two nonzero local messages are distinct; their images are a "
                "uniform ordered pair of distinct nonzero ambient words"
            ),
            "conditional_routing": (
                "conditional on the two image weights, the independent row-coordinate "
                "permutations produce independent uniform supports of those weights"
            ),
        },
        "mass_checks": {
            "rank_one_message_count_log2": message_nonzero / LOG2,
            "rank_two_message_count_log2": (
                message_nonzero + log_two_power_minus_two(dimension)
            ) / LOG2,
            "all_message_pair_count_log2": message_pair_mass / LOG2,
            "all_message_pair_count_expected_log2": (
                expected_message_pair_mass / LOG2
            ),
            "rank_two_image_pair_count_log2": image_pair_mass / LOG2,
            "rank_two_image_pair_count_expected_log2": (
                expected_image_pair_mass / LOG2
            ),
        },
        "top_rank_one_shells": [
            {
                "weight": int(index + 1),
                "pointwise_log2_upper": float(rank_one_shells[index] / LOG2),
                "log_surprisal": float(witnesses[index + 1, index + 1]),
            }
            for index in rank_one_top_indices
        ],
        "top_rank_two_shells": [
            {
                "first_weight": int(index // rank_two_width + 1),
                "second_weight": int(index % rank_two_width + 1),
                "pointwise_log2_upper": float(flat_rank_two[index] / LOG2),
                "log_surprisal": float(
                    witnesses[index // rank_two_width + 1, index % rank_two_width + 1]
                ),
            }
            for index in rank_two_top_indices
        ],
        "grid": {
            "coarse_minimum": args.grid_min,
            "coarse_maximum": args.grid_max,
            "coarse_step": args.grid_step,
            "fine_minimum": args.fine_min,
            "fine_maximum": args.fine_max,
            "fine_step": args.fine_step,
            "count": len(u_values),
        },
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only occupation two is covered.",
            "The calculation averages over the one sampled constituent; it does not identify a fixed constituent.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, default=256)
    parser.add_argument("--dimension", type=int, default=128)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--fine-min", type=float, default=-10.0)
    parser.add_argument("--fine-max", type=float, default=-7.0)
    parser.add_argument("--fine-step", type=float, default=0.05)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
