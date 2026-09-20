#!/usr/bin/env python3
"""Verify the dual covariance formula for a two-stage sparse-EA chain.

The default model exhaustively enumerates all degree-one maps with shape
2-to-3 followed by all degree-one square maps on three coordinates. It
checks direct shell variances, the law of total variance, and the composed
dual return-kernel formula in exact rational arithmetic.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "two_stage_dual_covariance_K2_B3_r1_1_exact.json"


def parity(value: int) -> int:
    return value.bit_count() & 1


def accumulator(value: int, length: int) -> int:
    state = 0
    result = 0
    for index in range(length):
        state ^= (value >> index) & 1
        result |= state << index
    return result


def accumulator_transpose(value: int, length: int) -> int:
    state = 0
    result = 0
    for index in range(length - 1, -1, -1):
        state ^= (value >> index) & 1
        result |= state << index
    return result


def apply_rows(rows: tuple[int, ...], value: int) -> int:
    result = 0
    for index, row in enumerate(rows):
        result |= parity(row & value) << index
    return result


def apply_rows_transpose(rows: tuple[int, ...], value: int) -> int:
    result = 0
    for index, row in enumerate(rows):
        if (value >> index) & 1:
            result ^= row
    return result


def fixed_weight_vectors(length: int, weight: int) -> tuple[int, ...]:
    return tuple(
        sum(1 << index for index in support)
        for support in itertools.combinations(range(length), weight)
    )


def maps(input_length: int, output_length: int, degree: int):
    rows = fixed_weight_vectors(input_length, degree)
    yield from itertools.product(rows, repeat=output_length)


def krawtchouk(length: int, shell: int, level: int) -> int:
    return sum(
        (-1) ** index
        * math.comb(level, index)
        * math.comb(length - level, shell - index)
        for index in range(
            max(0, shell - (length - level)), min(shell, level) + 1
        )
    )


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def mean_and_variance(values: list[int | Fraction]) -> tuple[Fraction, Fraction]:
    count = len(values)
    mean = sum(values, Fraction(0)) / count
    second = sum((Fraction(value) ** 2 for value in values), Fraction(0)) / count
    return mean, second - mean * mean


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=2)
    parser.add_argument("--output-bits", type=int, default=3)
    parser.add_argument("--first-degree", type=int, default=1)
    parser.add_argument("--second-degree", type=int, default=1)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    k = args.message_bits
    n = args.output_bits
    if not (1 <= args.first_degree <= k):
        parser.error("invalid first degree")
    if not (1 <= args.second_degree <= n):
        parser.error("invalid second degree")

    first_maps = tuple(maps(k, n, args.first_degree))
    second_maps = tuple(maps(n, n, args.second_degree))

    rho_single = []
    rho_pair = []
    for first_index in range(1 << n):
        single_count = sum(
            apply_rows_transpose(first, first_index) == 0 for first in first_maps
        )
        rho_single.append(Fraction(single_count, len(first_maps)))
        pair_row = []
        for second_index in range(1 << n):
            pair_count = sum(
                apply_rows_transpose(first, first_index) == 0
                and apply_rows_transpose(first, second_index) == 0
                for first in first_maps
            )
            pair_row.append(Fraction(pair_count, len(first_maps)))
        rho_pair.append(pair_row)

    dual_images = []
    for second in second_maps:
        row = []
        for final_character in range(1 << n):
            before_final_accumulator = accumulator_transpose(final_character, n)
            before_second_map = apply_rows_transpose(
                second, before_final_accumulator
            )
            before_first_accumulator = accumulator_transpose(before_second_map, n)
            row.append(before_first_accumulator)
        dual_images.append(row)

    dual_mean = []
    dual_pair = []
    for first_character in range(1 << n):
        dual_mean.append(
            sum(
                (rho_single[row[first_character]] for row in dual_images),
                Fraction(0),
            )
            / len(second_maps)
        )
        pair_row = []
        for second_character in range(1 << n):
            pair_row.append(
                sum(
                    (
                        rho_pair[row[first_character]][row[second_character]]
                        for row in dual_images
                    ),
                    Fraction(0),
                )
                / len(second_maps)
            )
        dual_pair.append(pair_row)

    direct_by_shell = [[] for _ in range(n + 1)]
    conditional_means = [[Fraction(0) for _ in first_maps] for _ in range(n + 1)]
    conditional_variances = [
        [Fraction(0) for _ in first_maps] for _ in range(n + 1)
    ]
    for first_position, first in enumerate(first_maps):
        conditional_counts = [[] for _ in range(n + 1)]
        for second in second_maps:
            shell_counts = [0 for _ in range(n + 1)]
            for message in range(1 << k):
                word = apply_rows(first, message)
                word = accumulator(word, n)
                word = apply_rows(second, word)
                word = accumulator(word, n)
                shell_counts[word.bit_count()] += 1
            for shell, value in enumerate(shell_counts):
                direct_by_shell[shell].append(value)
                conditional_counts[shell].append(value)
        for shell in range(n + 1):
            mean, variance = mean_and_variance(conditional_counts[shell])
            conditional_means[shell][first_position] = mean
            conditional_variances[shell][first_position] = variance

    rows = []
    for shell in range(n + 1):
        direct_mean, direct_variance = mean_and_variance(direct_by_shell[shell])
        mean_conditional_variance = sum(
            conditional_variances[shell], Fraction(0)
        ) / len(first_maps)
        _conditional_mean, variance_conditional_mean = mean_and_variance(
            conditional_means[shell]
        )
        total_variance = mean_conditional_variance + variance_conditional_mean

        coefficients = [
            krawtchouk(n, shell, character.bit_count())
            for character in range(1 << n)
        ]
        dual_variance = Fraction(0)
        for first_character in range(1 << n):
            for second_character in range(1 << n):
                covariance = (
                    dual_pair[first_character][second_character]
                    - dual_mean[first_character] * dual_mean[second_character]
                )
                dual_variance += (
                    coefficients[first_character]
                    * coefficients[second_character]
                    * covariance
                )
        dual_variance *= Fraction(1 << (2 * k), 1 << (2 * n))

        if direct_variance != total_variance:
            raise AssertionError(f"law of total variance failed at shell {shell}")
        if direct_variance != dual_variance:
            raise AssertionError(f"dual covariance failed at shell {shell}")
        rows.append(
            {
                "shell_weight": shell,
                "direct_mean": fraction_text(direct_mean),
                "direct_variance": fraction_text(direct_variance),
                "mean_conditional_variance": fraction_text(
                    mean_conditional_variance
                ),
                "variance_conditional_mean": fraction_text(
                    variance_conditional_mean
                ),
                "dual_variance": fraction_text(dual_variance),
            }
        )

    payload = {
        "schema": "two-stage-sparse-ea-dual-covariance-small-v1",
        "status": "EXACT_RATIONAL_SMALL_MODEL",
        "parameters": {
            "message_bits": k,
            "output_bits": n,
            "first_degree": args.first_degree,
            "second_degree": args.second_degree,
        },
        "ensemble_sizes": {
            "first_maps": len(first_maps),
            "second_maps": len(second_maps),
            "map_pairs": len(first_maps) * len(second_maps),
        },
        "verified": {
            "all_shell_laws_of_total_variance": True,
            "all_shell_dual_covariance_formulas": True,
        },
        "shells": rows,
        "scope": [
            "Every map and every message in the stated small ensemble is enumerated.",
            "All comparisons use exact rational arithmetic.",
            "The result validates the identity but does not extrapolate its numerical bounds to length 512.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), **payload["verified"]}, indent=2))


if __name__ == "__main__":
    main()
