#!/usr/bin/env python3
"""Evaluate the radial covariance-block relaxation at moderate dimensions.

This diagnostic combines the radial rank-(k+1) moment factorization, the
exact row-permutation block formula, and exact accumulator level energies.
It uses nondirected binary64 arithmetic for probabilities and block sums.
"""

from __future__ import annotations

import argparse
import functools
import json
import math
from pathlib import Path

import numpy as np

from analyze_dual_walk_and_accumulator_energy import (
    accumulator_joint_counts,
    krawtchouk_row,
)
from probe_one_stage_variance_scaling import accumulator_shell_law
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "radial_block_bound_K32_B64_r7_w40_probe.json"


def radial_walk_table(
    message_bits: int, right_degree: int, maximum_steps: int
) -> np.ndarray:
    """Return per-vector probabilities indexed by step and output weight."""
    denominator = math.comb(message_bits, right_degree)
    transition = np.zeros((message_bits + 1, message_bits + 1))
    for weight in range(message_bits + 1):
        low = max(0, right_degree - (message_bits - weight))
        high = min(weight, right_degree)
        for overlap in range(low, high + 1):
            following = weight + right_degree - 2 * overlap
            transition[weight, following] += (
                math.comb(weight, overlap)
                * math.comb(message_bits - weight, right_degree - overlap)
                / denominator
            )
    if np.max(np.abs(transition.sum(axis=1) - 1)) > 3e-15:
        raise AssertionError("radial transition is not stochastic")

    shell = np.zeros(message_bits + 1)
    shell[0] = 1
    table = np.zeros((maximum_steps + 1, message_bits + 1))
    multiplicities = np.asarray(
        [math.comb(message_bits, weight) for weight in range(message_bits + 1)],
        dtype=float,
    )
    table[0] = shell / multiplicities
    maximum_mass_residual = 0.0
    for step in range(1, maximum_steps + 1):
        shell = shell @ transition
        maximum_mass_residual = max(maximum_mass_residual, abs(shell.sum() - 1))
        table[step] = shell / multiplicities
    if maximum_mass_residual > 2e-13:
        raise AssertionError("radial walk mass residual is too large")
    return table


def exact_level_energies(output_bits: int, shell_weight: int) -> list[int]:
    joint = accumulator_joint_counts(output_bits)
    kraw = krawtchouk_row(output_bits, shell_weight)
    return [
        sum(
            count * kraw[derivative_weight] ** 2
            for derivative_weight, count in row.items()
        )
        for row in joint
    ]


def marginal_mean(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
) -> float:
    denominator = math.comb(message_bits, right_degree)
    mean = 0.0
    for weight in range(1, message_bits + 1):
        bias = krawtchouk(message_bits, right_degree, weight) / denominator
        law = accumulator_shell_law(np.longdouble((1 - bias) / 2), output_bits)
        mean += math.comb(message_bits, weight) * float(law[shell_weight])
    return mean


def evaluate(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
) -> dict[str, object]:
    walk = radial_walk_table(message_bits, right_degree, output_bits)
    multiplicities = np.asarray(
        [math.comb(message_bits, weight) for weight in range(message_bits + 1)],
        dtype=float,
    )

    @functools.lru_cache(maxsize=None)
    def covariance(first_weight: int, second_weight: int, intersection: int) -> float:
        first_only = first_weight - intersection
        second_only = second_weight - intersection
        joint = float(
            np.dot(
                multiplicities
                * walk[first_only]
                * walk[second_only],
                walk[intersection],
            )
        )
        return joint - walk[first_weight, 0] * walk[second_weight, 0]

    entry_norm = np.zeros((output_bits + 1, output_bits + 1))
    block_entries = 0
    maximum_symmetry_residual = 0.0
    for first_weight in range(output_bits + 1):
        for second_weight in range(first_weight, output_bits + 1):
            maximum_sector = min(
                first_weight,
                second_weight,
                output_bits - first_weight,
                output_bits - second_weight,
            )
            maximum_entry = 0.0
            for sector in range(maximum_sector + 1):
                outside = output_bits - 2 * sector
                first_outside = first_weight - sector
                second_outside = second_weight - sector
                first_count = math.comb(outside, first_outside)
                second_count = math.comb(outside, second_outside)
                terms = []
                for intersection in range(
                    max(0, first_weight + second_weight - output_bits),
                    min(first_weight, second_weight) + 1,
                ):
                    coefficient = 0
                    for same_pairs in range(sector + 1):
                        outside_intersection = intersection - same_pairs
                        if not 0 <= outside_intersection <= first_outside:
                            continue
                        remainder = second_outside - outside_intersection
                        if not 0 <= remainder <= outside - first_outside:
                            continue
                        coefficient += (
                            (-1) ** (sector - same_pairs)
                            * math.comb(sector, same_pairs)
                            * math.comb(first_outside, outside_intersection)
                            * math.comb(outside - first_outside, remainder)
                        )
                    if coefficient:
                        terms.append(
                            coefficient
                            * covariance(
                                first_weight, second_weight, intersection
                            )
                        )
                entry = math.sqrt(first_count / second_count) * math.fsum(terms)
                maximum_entry = max(maximum_entry, abs(entry))
                block_entries += 1
            entry_norm[first_weight, second_weight] = maximum_entry
            entry_norm[second_weight, first_weight] = maximum_entry
            maximum_symmetry_residual = max(
                maximum_symmetry_residual,
                abs(
                    entry_norm[first_weight, second_weight]
                    - entry_norm[second_weight, first_weight]
                ),
            )

    energies = exact_level_energies(output_bits, shell_weight)
    norms = np.sqrt(np.asarray(energies, dtype=float))
    quadratic_bound = float(norms @ entry_norm @ norms)
    mean = marginal_mean(
        message_bits, output_bits, right_degree, shell_weight
    )
    ratio_bound = (
        2.0 ** (2 * (message_bits - output_bits))
        * quadratic_bound
        / mean
    )
    return {
        "message_bits": message_bits,
        "output_bits": output_bits,
        "right_degree": right_degree,
        "shell_weight": shell_weight,
        "shell_mean": mean,
        "weight_level_block_variance_to_mean_upper": ratio_bound,
        "cached_covariance_entries": covariance.cache_info().currsize,
        "evaluated_block_entries": block_entries,
        "maximum_symmetry_residual": maximum_symmetry_residual,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=32)
    parser.add_argument("--output-bits", type=int, default=64)
    parser.add_argument("--right-degree", type=int, default=7)
    parser.add_argument("--shell-weight", type=int, default=40)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits > 128:
        parser.error("this binary64 probe is intentionally limited to output_bits <= 128")
    if args.output_bits != 2 * args.message_bits:
        parser.error("the present probe expects output_bits = 2 * message_bits")

    result = evaluate(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        args.shell_weight,
    )
    payload = {
        "schema": "pure-ea-radial-weight-level-block-bound-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "result": result,
        "scope": [
            "The radial and block formulas are exact identities.",
            "All probabilities and block sums in this receipt use nondirected binary64 arithmetic.",
            "The displayed upper bound is not rigorous without outward rounding.",
            "The receipt does not extrapolate to length 512.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
