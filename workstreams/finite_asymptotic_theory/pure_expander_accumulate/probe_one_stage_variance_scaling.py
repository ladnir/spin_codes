#!/usr/bin/env python3
"""Probe one-stage pure-EA shell variance at moderate dimensions.

This program evaluates the exact finite probability formulas in binary
floating-point arithmetic.  It averages over the exact multinomial law of
ordered rank-two message pairs.  It does not enumerate sparse maps.

The output is diagnostic: cancellation in the variance is performed with
``numpy.longdouble``, but no outward rounding is provided.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from verify_pair_kernel_small import compositions, krawtchouk, multinomial


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "one_stage_variance_scaling_probe.json"
TYPE4 = tuple[int, int, int, int]


def accumulator_shell_law(q: np.longdouble, length: int) -> np.ndarray:
    state = np.zeros((2, length + 1), dtype=np.longdouble)
    state[0, 0] = 1
    for step in range(length):
        following = np.zeros_like(state)
        following[0, : step + 1] += (1 - q) * state[0, : step + 1]
        following[1, 1 : step + 2] += q * state[0, : step + 1]
        following[1, 1 : step + 2] += (1 - q) * state[1, : step + 1]
        following[0, : step + 1] += q * state[1, : step + 1]
        state = following
    return state.sum(axis=0)


def pair_activation(
    message_bits: int,
    right_degree: int,
    input_type: TYPE4,
    biases: np.ndarray,
) -> np.ndarray:
    _n00, n01, n10, n11 = input_type
    h1 = n10 + n11
    h2 = n01 + n11
    h3 = n01 + n10
    b1, b2, b3 = biases[h1], biases[h2], biases[h3]
    return np.asarray(
        [
            (1 + b1 + b2 + b3) / 4,
            (1 + b1 - b2 - b3) / 4,
            (1 - b1 + b2 - b3) / 4,
            (1 - b1 - b2 + b3) / 4,
        ],
        dtype=np.longdouble,
    )


def pair_accumulator_diagonal_law(p: np.ndarray, length: int) -> np.ndarray:
    state = np.zeros((4, length + 1, length + 1), dtype=np.longdouble)
    state[0, 0, 0] = 1
    for step in range(length):
        following = np.zeros_like(state)
        size = step + 1
        for old_symbol in range(4):
            source = state[old_symbol, :size, :size]
            for input_symbol in range(4):
                next_symbol = old_symbol ^ input_symbol
                first = (next_symbol >> 1) & 1
                second = next_symbol & 1
                following[
                    next_symbol,
                    first : first + size,
                    second : second + size,
                ] += p[input_symbol] * source
        state = following
    diagonal = np.zeros(length + 1, dtype=np.longdouble)
    indices = np.arange(length + 1)
    for symbol in range(4):
        diagonal += state[symbol, indices, indices]
    return diagonal


def evaluate(message_bits: int, output_bits: int, right_degree: int) -> dict[str, object]:
    denominator = math.comb(message_bits, right_degree)
    biases = np.asarray(
        [
            np.longdouble(krawtchouk(message_bits, right_degree, weight))
            / np.longdouble(denominator)
            for weight in range(message_bits + 1)
        ],
        dtype=np.longdouble,
    )

    nonzero_messages = (1 << message_bits) - 1
    rank_two_pairs = nonzero_messages * (nonzero_messages - 1)
    marginal = np.zeros(output_bits + 1, dtype=np.longdouble)
    single_laws = []
    single_laws.append(accumulator_shell_law(np.longdouble(0), output_bits))
    for weight in range(1, message_bits + 1):
        q = (1 - biases[weight]) / 2
        law = accumulator_shell_law(q, output_bits)
        single_laws.append(law)
        marginal += (
            np.longdouble(math.comb(message_bits, weight))
            / np.longdouble(nonzero_messages)
        ) * law

    joint = np.zeros(output_bits + 1, dtype=np.longdouble)
    covariance_by_difference = np.zeros(
        (message_bits + 1, output_bits + 1), dtype=np.longdouble
    )
    retained_types = 0
    for input_type in compositions(message_bits, 4):
        _n00, n01, n10, n11 = input_type
        if min(n10 + n11, n01 + n11, n01 + n10) == 0:
            continue
        retained_types += 1
        type_probability = np.longdouble(multinomial(input_type)) / np.longdouble(
            rank_two_pairs
        )
        p = pair_activation(message_bits, right_degree, input_type, biases)
        if np.min(p) < -np.longdouble("1e-13"):
            raise AssertionError("negative pair activation probability")
        p = np.maximum(p, 0)
        p /= p.sum()
        pair_law = pair_accumulator_diagonal_law(p, output_bits)
        joint += type_probability * pair_law
        h1 = n10 + n11
        h2 = n01 + n11
        h3 = n01 + n10
        covariance_by_difference[h3] += type_probability * (
            pair_law - single_laws[h1] * single_laws[h2]
        )

    joint_equal_weight_mass = joint.sum()
    rows = []
    maximum_ratio = np.longdouble(0)
    maximum_weight = 0
    for weight in range(output_bits + 1):
        pi = marginal[weight]
        if pi == 0:
            ratio = None
        else:
            ratio = (
                1
                + np.longdouble(nonzero_messages - 1) * joint[weight] / pi
                - np.longdouble(nonzero_messages) * pi
            )
            if weight > 0 and ratio > maximum_ratio:
                maximum_ratio = ratio
                maximum_weight = weight
        rows.append(
            {
                "weight": weight,
                "mean_log2": None
                if pi == 0
                else float(np.log2(pi) + np.longdouble(message_bits)),
                "variance_to_mean": None if ratio is None else float(ratio),
            }
        )

    worst_pi = marginal[maximum_weight]
    difference_rows = []
    off_diagonal_ratio_sum = np.longdouble(0)
    for difference_weight in range(1, message_bits + 1):
        contribution = (
            np.longdouble(nonzero_messages - 1)
            * covariance_by_difference[difference_weight, maximum_weight]
            / worst_pi
        )
        off_diagonal_ratio_sum += contribution
        difference_rows.append(
            {
                "difference_weight": difference_weight,
                "variance_to_mean_contribution": float(contribution),
            }
        )
    complement_pairs = []
    for difference_weight in range(1, message_bits // 2 + 1):
        complement = message_bits - difference_weight
        first = np.longdouble(
            difference_rows[difference_weight - 1]["variance_to_mean_contribution"]
        )
        if complement == difference_weight:
            paired = first
        else:
            paired = first + np.longdouble(
                difference_rows[complement - 1]["variance_to_mean_contribution"]
            )
        complement_pairs.append(
            {
                "difference_weights": (
                    [difference_weight]
                    if complement == difference_weight
                    else [difference_weight, complement]
                ),
                "variance_to_mean_contribution": float(paired),
            }
        )
    diagonal_numerator = np.longdouble(0)
    for weight in range(1, message_bits + 1):
        probability = np.longdouble(math.comb(message_bits, weight)) / np.longdouble(
            nonzero_messages
        )
        pi = single_laws[weight][maximum_weight]
        diagonal_numerator += probability * pi * (1 - pi)
    diagonal_ratio = diagonal_numerator / worst_pi

    return {
        "message_bits": message_bits,
        "output_bits": output_bits,
        "right_degree": right_degree,
        "right_degree_fraction": right_degree / message_bits,
        "rank_two_pair_types": retained_types,
        "joint_equal_weight_probability": float(joint_equal_weight_mass),
        "maximum_positive_shell_variance_to_mean": float(maximum_ratio),
        "maximum_positive_shell_weight": maximum_weight,
        "worst_shell_decomposition": {
            "weight": maximum_weight,
            "diagonal_variance_to_mean": float(diagonal_ratio),
            "off_diagonal_variance_to_mean": float(off_diagonal_ratio_sum),
            "reconstructed_variance_to_mean": float(
                diagonal_ratio + off_diagonal_ratio_sum
            ),
            "by_message_difference_weight": difference_rows,
            "by_complement_paired_difference_weights": complement_pairs,
        },
        "shells": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=12)
    parser.add_argument("--output-bits", type=int)
    parser.add_argument("--right-degree", type=int, default=2)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output_bits = args.output_bits or 2 * args.message_bits
    if not 1 <= args.right_degree <= args.message_bits:
        parser.error("right degree must lie in [1,message bits]")
    if output_bits < 1:
        parser.error("output bits must be positive")

    result = evaluate(args.message_bits, output_bits, args.right_degree)
    payload = {
        "schema": "pure-ea-one-stage-variance-scaling-probe-v1",
        "status": "BINARY_FLOAT_DIAGNOSTIC",
        "arithmetic": {
            "dtype": "numpy.longdouble",
            "mantissa_bits": int(np.finfo(np.longdouble).nmant),
            "directed_rounding": False,
        },
        "result": result,
        "scope": [
            "The probability formulas are exact identities for the ensemble.",
            "The numerical evaluation is not an interval certificate.",
            "The result does not extrapolate to other dimensions.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "message_bits",
                    "output_bits",
                    "right_degree",
                    "rank_two_pair_types",
                    "joint_equal_weight_probability",
                    "maximum_positive_shell_variance_to_mean",
                    "maximum_positive_shell_weight",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
