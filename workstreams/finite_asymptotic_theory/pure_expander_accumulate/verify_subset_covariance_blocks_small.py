#!/usr/bin/env python3
"""Verify the dual covariance formula and its symmetric-group blocks.

The default instance is small enough to enumerate every sparse map.  Kernel
entries and the shell-variance identity use exact rational arithmetic.  The
comparison of the full and block eigenvalue multisets uses binary64.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

from verify_pair_kernel_small import accumulator, krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "subset_covariance_blocks_K3_B4_r1_exact.json"


def region_counts(first: int, second: int, length: int) -> tuple[int, int, int]:
    mask = (1 << length) - 1
    first_only = (first & (mask ^ second)).bit_count()
    second_only = (second & (mask ^ first)).bit_count()
    intersection = (first & second).bit_count()
    return first_only, second_only, intersection


def row_biases(message_bits: int, right_degree: int) -> list[Fraction]:
    denominator = math.comb(message_bits, right_degree)
    return [
        Fraction(krawtchouk(message_bits, right_degree, weight), denominator)
        for weight in range(message_bits + 1)
    ]


def covariance_kernel(
    message_bits: int, output_bits: int, right_degree: int
) -> tuple[list[Fraction], dict[tuple[int, int, int], Fraction]]:
    biases = row_biases(message_bits, right_degree)
    messages = range(1 << message_bits)
    returns = []
    for steps in range(output_bits + 1):
        returns.append(
            sum(biases[message.bit_count()] ** steps for message in messages)
            / (1 << message_bits)
        )

    kernel = {}
    for first_only in range(output_bits + 1):
        for second_only in range(output_bits - first_only + 1):
            for intersection in range(
                output_bits - first_only - second_only + 1
            ):
                joint = sum(
                    biases[first.bit_count()] ** first_only
                    * biases[second.bit_count()] ** second_only
                    * biases[(first ^ second).bit_count()] ** intersection
                    for first in messages
                    for second in messages
                ) / (1 << (2 * message_bits))
                kernel[(first_only, second_only, intersection)] = joint - (
                    returns[first_only + intersection]
                    * returns[second_only + intersection]
                )
    return returns, kernel


def full_covariance(output_bits: int, kernel: dict[tuple[int, int, int], Fraction]):
    subsets = range(1 << output_bits)
    return [
        [kernel[region_counts(first, second, output_bits)] for second in subsets]
        for first in subsets
    ]


def covariance_blocks(
    output_bits: int, kernel: dict[tuple[int, int, int], Fraction]
) -> list[dict[str, object]]:
    blocks = []
    for sector in range(output_bits // 2 + 1):
        outside = output_bits - 2 * sector
        levels = list(range(sector, output_bits - sector + 1))
        matrix = np.zeros((len(levels), len(levels)), dtype=float)
        for row, first_weight in enumerate(levels):
            first_outside = first_weight - sector
            first_count = math.comb(outside, first_outside)
            for column, second_weight in enumerate(levels):
                second_outside = second_weight - sector
                second_count = math.comb(outside, second_outside)
                value = Fraction(0)
                for same_pairs in range(sector + 1):
                    pair_factor = (
                        (-1) ** (sector - same_pairs)
                        * math.comb(sector, same_pairs)
                    )
                    low = max(
                        0, second_outside - (outside - first_outside)
                    )
                    high = min(first_outside, second_outside)
                    for outside_intersection in range(low, high + 1):
                        ways = math.comb(
                            first_outside, outside_intersection
                        ) * math.comb(
                            outside - first_outside,
                            second_outside - outside_intersection,
                        )
                        total_intersection = same_pairs + outside_intersection
                        value += pair_factor * ways * kernel[
                            (
                                first_weight - total_intersection,
                                second_weight - total_intersection,
                                total_intersection,
                            )
                        ]
                matrix[row, column] = float(value) * math.sqrt(
                    first_count / second_count
                )
        symmetry_residual = float(np.max(np.abs(matrix - matrix.T)))
        if symmetry_residual > 2e-12:
            raise AssertionError("covariance block is not symmetric")
        dimension = math.comb(output_bits, sector) - (
            math.comb(output_bits, sector - 1) if sector else 0
        )
        blocks.append(
            {
                "sector": sector,
                "levels": levels,
                "irreducible_dimension": dimension,
                "matrix": matrix,
                "symmetry_residual": symmetry_residual,
            }
        )
    return blocks


def derivative_weight(subset: int) -> int:
    return (subset ^ (subset >> 1)).bit_count()


def shell_coefficient(output_bits: int, shell_weight: int, subset: int) -> int:
    return krawtchouk(
        output_bits, shell_weight, derivative_weight(subset)
    )


def apply_sparse_map(message: int, rows: tuple[int, ...]) -> int:
    output = 0
    for index, row in enumerate(rows):
        output |= ((message & row).bit_count() & 1) << index
    return output


def exhaustive_shell_variance(
    message_bits: int, output_bits: int, right_degree: int, shell_weight: int
) -> tuple[Fraction, Fraction, int]:
    rows = tuple(
        sum(1 << index for index in support)
        for support in itertools.combinations(range(message_bits), right_degree)
    )
    counts = []
    for matrix_rows in itertools.product(rows, repeat=output_bits):
        shell_count = 0
        for message in range(1 << message_bits):
            encoded = accumulator(
                apply_sparse_map(message, matrix_rows), output_bits
            )
            shell_count += encoded.bit_count() == shell_weight
        counts.append(shell_count)
    mean = Fraction(sum(counts), len(counts))
    second = Fraction(sum(count * count for count in counts), len(counts))
    return mean, second - mean * mean, len(counts)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=3)
    parser.add_argument("--output-bits", type=int, default=4)
    parser.add_argument("--right-degree", type=int, default=1)
    parser.add_argument("--shell-weight", type=int, default=2)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 1 <= args.right_degree <= args.message_bits:
        parser.error("invalid right degree")
    if not 0 <= args.shell_weight <= args.output_bits:
        parser.error("invalid shell weight")
    if args.output_bits > 8:
        parser.error("this exhaustive verifier is limited to output_bits <= 8")

    returns, kernel = covariance_kernel(
        args.message_bits, args.output_bits, args.right_degree
    )
    full_exact = full_covariance(args.output_bits, kernel)
    full_float = np.asarray(
        [[float(value) for value in row] for row in full_exact], dtype=float
    )
    blocks = covariance_blocks(args.output_bits, kernel)

    full_eigenvalues = sorted(np.linalg.eigvalsh(full_float).tolist())
    block_eigenvalues = []
    for block in blocks:
        values = np.linalg.eigvalsh(block["matrix"]).tolist()
        block_eigenvalues.extend(values * int(block["irreducible_dimension"]))
    block_eigenvalues.sort()
    eigenvalue_residual = max(
        abs(first - second)
        for first, second in zip(full_eigenvalues, block_eigenvalues, strict=True)
    )
    if eigenvalue_residual > 3e-12:
        raise AssertionError("block and full covariance spectra differ")

    coefficients = [
        shell_coefficient(args.output_bits, args.shell_weight, subset)
        for subset in range(1 << args.output_bits)
    ]
    quadratic = sum(
        coefficients[first]
        * full_exact[first][second]
        * coefficients[second]
        for first in range(1 << args.output_bits)
        for second in range(1 << args.output_bits)
    )
    dual_variance = Fraction(2) ** (
        2 * (args.message_bits - args.output_bits)
    ) * quadratic
    exhaustive_mean, exhaustive_variance, map_count = exhaustive_shell_variance(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        args.shell_weight,
    )
    if dual_variance != exhaustive_variance:
        raise AssertionError("dual and exhaustive shell variances differ")

    block_rows = []
    for block in blocks:
        eigenvalues = np.linalg.eigvalsh(block["matrix"])
        block_rows.append(
            {
                "sector": block["sector"],
                "levels": block["levels"],
                "irreducible_dimension": block["irreducible_dimension"],
                "maximum_eigenvalue": float(eigenvalues[-1]),
                "minimum_eigenvalue": float(eigenvalues[0]),
                "symmetry_residual": block["symmetry_residual"],
            }
        )

    payload = {
        "schema": "pure-ea-subset-covariance-blocks-small-v1",
        "status": "EXACT_RATIONAL_IDENTITIES_WITH_BINARY64_EIGENVALUE_CHECK",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "shell_weight": args.shell_weight,
            "enumerated_maps": map_count,
        },
        "return_probabilities": [fraction_text(value) for value in returns],
        "shell_mean": fraction_text(exhaustive_mean),
        "exhaustive_shell_variance": fraction_text(exhaustive_variance),
        "dual_quadratic_shell_variance": fraction_text(dual_variance),
        "full_vs_block_eigenvalue_maximum_residual": eigenvalue_residual,
        "blocks": block_rows,
        "verified_identities": [
            "dual covariance quadratic form equals exhaustive ensemble shell variance",
            "symmetric-group block eigenvalues reproduce the full covariance spectrum",
        ],
        "scope": [
            "Rational probability and variance comparisons are exact.",
            "The eigenvalue comparison uses nondirected binary64 arithmetic.",
            "The result validates only the displayed small instance.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
