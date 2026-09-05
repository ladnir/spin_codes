#!/usr/bin/env python3
"""Verify the primal Schur-power factorization for a small sparse map.

For a fixed ordered message pair, one sparse-map row has a 2-by-2 output
probability matrix P.  The n-row pair kernel is P tensor n.  Its two-row
Specht block in sector j is det(P)^j times the normalized symmetric power
of degree n-2j.  This program compares the aggregate of those blocks with
the Walsh conjugate of the exact dual covariance blocks.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from probe_walsh_conjugated_level_bound import invariant_blocks, walsh_blocks
from verify_pair_kernel_small import krawtchouk
from verify_subset_covariance_blocks_small import covariance_kernel


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "primal_schur_power_K4_B8_r3_probe.json"


def symmetric_power_block(matrix: np.ndarray, degree: int) -> np.ndarray:
    """Return the normalized weight-basis action of matrix tensor degree."""
    a, b = matrix[0]
    c, d = matrix[1]
    result = np.zeros((degree + 1, degree + 1), dtype=float)
    for output_weight in range(degree + 1):
        output_count = math.comb(degree, output_weight)
        for input_weight in range(degree + 1):
            input_count = math.comb(degree, input_weight)
            terms = []
            low = max(0, output_weight + input_weight - degree)
            high = min(output_weight, input_weight)
            for intersection in range(low, high + 1):
                terms.append(
                    math.comb(output_weight, intersection)
                    * math.comb(degree - output_weight, input_weight - intersection)
                    * a ** (degree - output_weight - input_weight + intersection)
                    * b ** (input_weight - intersection)
                    * c ** (output_weight - intersection)
                    * d ** intersection
                )
            result[output_weight, input_weight] = math.sqrt(
                output_count / input_count
            ) * math.fsum(terms)
    return result


def pair_probability(
    first: int, second: int, biases: list[float]
) -> np.ndarray:
    first_bias = biases[first.bit_count()]
    second_bias = biases[second.bit_count()]
    difference_bias = biases[(first ^ second).bit_count()]
    return np.asarray(
        [
            [
                (1 + first_bias + second_bias + difference_bias) / 4,
                (1 + first_bias - second_bias - difference_bias) / 4,
            ],
            [
                (1 - first_bias + second_bias - difference_bias) / 4,
                (1 - first_bias - second_bias + difference_bias) / 4,
            ],
        ]
    )


def evaluate(message_bits: int, output_bits: int, right_degree: int) -> dict[str, object]:
    if output_bits != 2 * message_bits:
        raise ValueError("the verifier expects rate one half")
    denominator = math.comb(message_bits, right_degree)
    biases = [
        krawtchouk(message_bits, right_degree, weight) / denominator
        for weight in range(message_bits + 1)
    ]

    primal_blocks = [
        np.zeros((output_bits - 2 * sector + 1,) * 2, dtype=float)
        for sector in range(output_bits // 2 + 1)
    ]
    mean_by_word_weight = np.zeros(output_bits + 1, dtype=float)
    for message in range(1 << message_bits):
        q = (1 - biases[message.bit_count()]) / 2
        for weight in range(output_bits + 1):
            mean_by_word_weight[weight] += (
                q**weight * (1 - q) ** (output_bits - weight)
            )

    for first in range(1 << message_bits):
        for second in range(1 << message_bits):
            probability = pair_probability(first, second, biases)
            determinant = float(np.linalg.det(probability))
            for sector, block in enumerate(primal_blocks):
                degree = output_bits - 2 * sector
                block += determinant**sector * symmetric_power_block(
                    probability, degree
                )

    radial_mean = np.asarray(
        [
            math.sqrt(math.comb(output_bits, weight)) * mean_by_word_weight[weight]
            for weight in range(output_bits + 1)
        ]
    )
    primal_blocks[0] -= np.outer(radial_mean, radial_mean)

    _returns, exact_kernel = covariance_kernel(
        message_bits, output_bits, right_degree
    )
    dual_blocks = invariant_blocks(
        output_bits,
        lambda first, second, intersection: float(
            exact_kernel[
                (first - intersection, second - intersection, intersection)
            ]
        ),
    )
    transforms, involution_residual = walsh_blocks(output_bits)
    normalization = math.ldexp(1.0, -output_bits)
    conjugated = [
        normalization * transform @ block @ transform.T
        for transform, block in zip(transforms, dual_blocks, strict=True)
    ]
    scale = math.ldexp(1.0, output_bits - 2 * message_bits)
    residuals = [
        float(np.max(np.abs(left - scale * right)))
        for left, right in zip(conjugated, primal_blocks, strict=True)
    ]
    relative = [
        residual / max(1.0, float(np.max(np.abs(left))))
        for residual, left in zip(residuals, conjugated, strict=True)
    ]
    maximum = max(relative)
    if maximum > 3e-12:
        raise AssertionError("primal Schur-power blocks do not match conjugated blocks")
    return {
        "message_bits": message_bits,
        "output_bits": output_bits,
        "right_degree": right_degree,
        "ordered_message_pairs": 1 << (2 * message_bits),
        "maximum_relative_block_residual": maximum,
        "sector_relative_block_residuals": relative,
        "normalized_walsh_involution_residual": involution_residual,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--output-bits", type=int, default=8)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.message_bits > 5 or args.output_bits > 10:
        parser.error("this verifier is intentionally limited to k <= 5 and n <= 10")
    result = evaluate(args.message_bits, args.output_bits, args.right_degree)
    payload = {
        "schema": "pure-ea-primal-schur-power-factorization-small-v1",
        "status": "BINARY64_IDENTITY_CHECK",
        "result": result,
        "scope": [
            "The Schur-power formula is algebraic; this comparison evaluates it in binary64.",
            "The dual covariance kernel used in the comparison is exact rational arithmetic.",
            "The receipt validates only the displayed small instance.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
