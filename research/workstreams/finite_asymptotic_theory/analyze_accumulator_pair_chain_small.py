#!/usr/bin/env python3
"""Analyze exact pair-type Markov chains at small block lengths.

The chain is induced by one common uniform coordinate permutation followed by
prefix accumulation of both words. This program materializes the rank-two
pair-type chain for one or two RM(1,3) blocks. It verifies stationarity and
propagates the exact pair distribution of the direct-sum code. The result is
a diagnostic for a length-512 contraction proof, not an extrapolation.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from verify_accumulator_pair_type_kernel import kernel_counts, multinomial, pair_type


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "accumulator_pair_chain_B8_exact.json"
LOCAL_BITS = 8


Type4 = tuple[int, int, int, int]


def compositions4(total: int) -> list[Type4]:
    return [
        (a, b, c, total - a - b - c)
        for a in range(total + 1)
        for b in range(total - a + 1)
        for c in range(total - a - b + 1)
    ]


def rank_two(kind: Type4) -> bool:
    # Symbol order is 00, 01, 10, 11.
    first_nonzero = kind[2] + kind[3] > 0
    second_nonzero = kind[1] + kind[3] > 0
    distinct = kind[1] + kind[2] > 0
    return first_nonzero and second_nonzero and distinct


def rm13_words() -> list[int]:
    points = [
        (value >> 2 & 1, value >> 1 & 1, value & 1)
        for value in range(LOCAL_BITS)
    ]
    basis = []
    for coordinate in range(4):
        word = 0
        for index, point in enumerate(points):
            bit = 1 if coordinate == 0 else point[coordinate - 1]
            word |= bit << index
        basis.append(word)
    words = []
    for message in range(1 << 4):
        word = 0
        for index, row in enumerate(basis):
            if message >> index & 1:
                word ^= row
        words.append(word)
    return words


def direct_sum_words(blocks: int) -> list[int]:
    local = rm13_words()
    return [
        sum(word << (LOCAL_BITS * block) for block, word in enumerate(parts))
        for parts in itertools.product(local, repeat=blocks)
    ]


def code_pair_distribution(
    types: list[Type4], words: list[int], block_bits: int
) -> np.ndarray:
    index = {kind: position for position, kind in enumerate(types)}
    nonzero_words = words[1:]
    counts = np.zeros(len(types), dtype=np.float64)
    for first in nonzero_words:
        for second in nonzero_words:
            if first != second:
                counts[index[pair_type(first, second, block_bits)]] += 1.0
    expected = len(nonzero_words) * (len(nonzero_words) - 1)
    if counts.sum() != expected:
        raise AssertionError("RM(1,3) pair distribution has the wrong mass")
    return counts / expected


def shell_statistics(
    distribution: np.ndarray,
    types: list[Type4],
    nonzero_words: int,
    block_bits: int,
) -> list[dict[str, float | int]]:
    ordered_pairs = nonzero_words * (nonzero_words - 1)
    rows = []
    for weight in range(1, block_bits + 1):
        marginal = sum(
            distribution[index]
            for index, kind in enumerate(types)
            if kind[2] + kind[3] == weight
        )
        joint = sum(
            distribution[index]
            for index, kind in enumerate(types)
            if kind[2] + kind[3] == weight and kind[1] + kind[3] == weight
        )
        mean = nonzero_words * marginal
        factorial_second = ordered_pairs * joint
        variance = max(0.0, factorial_second + mean - mean * mean)
        rows.append(
            {
                "weight": weight,
                "mean": mean,
                "variance": variance,
                "variance_over_mean": variance / mean if mean else 0.0,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-stages", type=int, default=6)
    parser.add_argument("--blocks", type=int, choices=(1, 2), default=1)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    block_bits = LOCAL_BITS * args.blocks
    words = direct_sum_words(args.blocks)
    nonzero_words = len(words) - 1
    types = [kind for kind in compositions4(block_bits) if rank_two(kind)]
    index = {kind: position for position, kind in enumerate(types)}
    transition = np.zeros((len(types), len(types)), dtype=np.float64)
    for row, input_type in enumerate(types):
        counts = kernel_counts(input_type)
        denominator = multinomial(input_type)
        for output_type, count in counts.items():
            if output_type not in index:
                raise AssertionError("invertible accumulator changed pair rank")
            transition[row, index[output_type]] = count / denominator
    if not np.allclose(transition.sum(axis=1), 1.0, atol=2e-15):
        raise AssertionError("transition rows are not stochastic")

    stationary_counts = np.asarray([multinomial(kind) for kind in types], dtype=np.float64)
    stationary = stationary_counts / stationary_counts.sum()
    expected_rank_two_pairs = (1 << block_bits) - 1
    expected_rank_two_pairs *= (1 << block_bits) - 2
    if int(stationary_counts.sum()) != expected_rank_two_pairs:
        raise AssertionError("rank-two orbit masses have the wrong total")
    stationary_error = float(np.max(np.abs(stationary @ transition - stationary)))

    weighted = (
        np.sqrt(stationary)[:, None]
        * transition
        / np.sqrt(stationary)[None, :]
    )
    singular_values = np.linalg.svd(weighted, compute_uv=False)

    distribution = code_pair_distribution(types, words, block_bits)
    stages = []
    transition_power = np.eye(len(types))
    for stage in range(args.maximum_stages + 1):
        shell_rows = shell_statistics(
            distribution, types, nonzero_words, block_bits
        )
        active = [row for row in shell_rows if float(row["mean"]) > 1e-14]
        worst = max(active, key=lambda row: float(row["variance_over_mean"]))
        chi_square = float(np.sum((distribution - stationary) ** 2 / stationary))
        maximum_ratio = float(np.max(distribution / stationary))
        worst_point_chi_square = 0.0
        worst_point_ratio = 0.0
        for row in transition_power:
            worst_point_chi_square = max(
                worst_point_chi_square,
                float(np.sum((row - stationary) ** 2 / stationary)),
            )
            worst_point_ratio = max(worst_point_ratio, float(np.max(row / stationary)))
        stages.append(
            {
                "accumulator_stages": stage,
                "rm13_pair_chi_square_from_stationary": chi_square,
                "rm13_pair_maximum_density_ratio": maximum_ratio,
                "worst_start_chi_square_from_stationary": worst_point_chi_square,
                "worst_start_maximum_density_ratio": worst_point_ratio,
                "worst_rm13_shell_variance_ratio": worst,
                "rm13_shells": shell_rows,
            }
        )
        print(
            f"stage,{stage},rm_chi2,{chi_square:.9g},"
            f"rm_density,{maximum_ratio:.9g},shell_var_ratio,{worst['variance_over_mean']:.9g}",
            flush=True,
        )
        distribution = distribution @ transition
        transition_power = transition_power @ transition

    result = {
        "schema": "accumulator-pair-chain-small-exact-v2",
        "status": "EXACT_KERNEL_WITH_BINARY64_LINEAR_ALGEBRA_DIAGNOSTIC",
        "parameters": {
            "length": block_bits,
            "rank_two_pair_types": len(types),
            "base_code": (
                f"direct sum of {args.blocks} RM(1,3) [8,4,4] blocks"
            ),
            "maximum_accumulator_stages": args.maximum_stages,
        },
        "checks": {
            "maximum_stationarity_error": stationary_error,
            "largest_weighted_singular_values": singular_values[:10].tolist(),
        },
        "stages": stages,
        "limitations": [
            "Kernel entries are exact integer ratios but spectral diagnostics use binary64.",
            "Small-block contraction factors do not prove a length-512 bound.",
            "Worst-start mixing is stronger than the shell variance needed by the theorem target.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
