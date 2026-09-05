#!/usr/bin/env python3
"""Probe a positive pair bound that preserves the dominant row character.

For each ordered nonzero distinct message pair, the probe keeps the largest
of the three row-character biases exactly.  It bounds the other two biases
by a pointwise likelihood ratio.  Negative pairwise excess is discarded.
The resulting variance bound is valid algebraically, but this implementation
uses nondirected binary64 aggregation and is therefore only a diagnostic.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

from analyze_dual_walk_and_accumulator_energy import (
    accumulator_joint_counts,
    krawtchouk_row,
)
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "dominant_character_likelihood_bound_probe.json"


@lru_cache(maxsize=None)
def cached_joint_counts(output_bits: int) -> list[dict[int, int]]:
    return accumulator_joint_counts(output_bits)


def polynomial_ratios(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
    precision: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return one-word and unbiased-correlated pair ratios to random."""
    mp.mp.dps = precision
    joint = cached_joint_counts(output_bits)
    krawtchouk_values = krawtchouk_row(output_bits, shell_weight)
    shell_size = math.comb(output_bits, shell_weight)
    one_coefficients = []
    pair_coefficients = []
    for level in range(output_bits + 1):
        one = 0
        pair = 0
        for derivative_weight, count in joint[level].items():
            coefficient = krawtchouk_values[derivative_weight]
            one += count * coefficient
            pair += count * coefficient * coefficient
        one_coefficients.append(mp.mpf(one) / shell_size)
        pair_coefficients.append(mp.mpf(pair) / (shell_size * shell_size))

    denominator = math.comb(message_bits, right_degree)
    one_ratios = np.empty(message_bits + 1, dtype=np.float64)
    pair_ratios = np.empty(message_bits + 1, dtype=np.float64)
    for weight in range(message_bits + 1):
        numerator = krawtchouk(message_bits, right_degree, weight)
        if abs(numerator) == denominator:
            bias = 1 if numerator > 0 else -1
            if shell_weight == 0 and bias == 1:
                one_ratios[weight] = math.ldexp(1.0, output_bits) / shell_size
            elif bias == -1 and shell_weight == (output_bits + 1) // 2:
                one_ratios[weight] = math.ldexp(1.0, output_bits) / shell_size
            else:
                one_ratios[weight] = 0.0
            # Under the correlated reference, beta=-1 fixes the XOR of the
            # two accumulated words to the alternating word of weight n/2.
            if bias == -1 and 2 * shell_weight < output_bits // 2:
                pair_ratios[weight] = 0.0
            else:
                value = mp.polyval(
                    list(reversed(pair_coefficients)), mp.mpf(bias)
                )
                pair_ratios[weight] = float(value)
            continue
        bias = mp.mpf(numerator) / denominator
        one_value = mp.polyval(list(reversed(one_coefficients)), bias)
        pair_value = mp.polyval(list(reversed(pair_coefficients)), bias)
        if one_value < 0 or pair_value < 0:
            raise ArithmeticError(
                "a probability ratio became negative: "
                f"shell={shell_weight}, weight={weight}, "
                f"one={one_value}, pair={pair_value}"
            )
        one_ratios[weight] = float(one_value)
        pair_ratios[weight] = float(pair_value)
    return one_ratios, pair_ratios


def pair_bound(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weights: list[int],
    precision: int,
) -> list[dict[str, object]]:
    denominator = math.comb(message_bits, right_degree)
    biases = np.asarray(
        [
            krawtchouk(message_bits, right_degree, weight) / denominator
            for weight in range(message_bits + 1)
        ],
        dtype=np.float64,
    )
    absolute_biases = np.abs(biases)
    ratio_tables = {
        shell: polynomial_ratios(
            message_bits,
            output_bits,
            right_degree,
            shell,
            precision,
        )
        for shell in shell_weights
    }

    positive_excess = {shell: np.longdouble(0) for shell in shell_weights}
    pair_mass = np.longdouble(0)
    pair_types = 0
    dominant_counts = [0, 0, 0]
    for n11 in range(message_bits + 1):
        first_factor = math.comb(message_bits, n11)
        remaining = message_bits - n11
        for n10 in range(remaining + 1):
            second_factor = math.comb(remaining, n10)
            last = remaining - n10
            n01 = np.arange(last + 1, dtype=np.int64)
            first_weight = n10 + n11
            second_weight = n01 + n11
            difference_weight = n10 + n01
            valid = (
                (first_weight != 0)
                & (second_weight != 0)
                & (difference_weight != 0)
            )
            if not np.any(valid):
                continue
            second_weight = second_weight[valid]
            difference_weight = difference_weight[valid]
            n01_valid = n01[valid]
            multiplicities = np.asarray(
                [
                    float(
                        first_factor
                        * second_factor
                        * math.comb(last, int(value))
                    )
                    for value in n01_valid
                ],
                dtype=np.float64,
            )
            weights = np.ldexp(multiplicities, -2 * message_bits)
            first_weights = np.full_like(second_weight, first_weight)
            triples = np.stack(
                (
                    absolute_biases[first_weights],
                    absolute_biases[second_weight],
                    absolute_biases[difference_weight],
                )
            )
            dominant = np.argmax(triples, axis=0)
            dominant_bias = np.take_along_axis(
                triples, dominant[np.newaxis, :], axis=0
            )[0]
            residual = np.sum(triples, axis=0) - dominant_bias
            with np.errstate(divide="ignore", invalid="ignore"):
                delta = residual / (1.0 - dominant_bias)
            log_likelihood = output_bits * np.log1p(delta)
            likelihood_excess = np.expm1(
                np.minimum(log_likelihood, math.log(np.finfo(float).max))
            )

            pair_mass += np.sum(weights, dtype=np.longdouble)
            pair_types += len(weights)
            for index in range(3):
                dominant_counts[index] += int(np.count_nonzero(dominant == index))

            for shell in shell_weights:
                one, correlated = ratio_tables[shell]
                first_ratio = one[first_weights]
                second_ratio = one[second_weight]
                independent_ratio = (
                    1.0
                    + (first_ratio - 1.0)
                    + (second_ratio - 1.0)
                    + (first_ratio - 1.0) * (second_ratio - 1.0)
                )
                base = np.where(
                    dominant == 0,
                    first_ratio,
                    np.where(
                        dominant == 1,
                        second_ratio,
                        correlated[difference_weight],
                    ),
                )
                with np.errstate(over="ignore", invalid="ignore"):
                    finite_candidate = base + base * likelihood_excess
                candidate = np.where(
                    np.isfinite(likelihood_excess), finite_candidate, np.inf
                )
                random_shell_probability = math.comb(
                    output_bits, shell
                ) / math.ldexp(1.0, output_bits)
                marginal_first = first_ratio / random_shell_probability
                marginal_second = second_ratio / random_shell_probability
                upper_ratio = np.minimum(
                    candidate, np.minimum(marginal_first, marginal_second)
                )
                impossible = (first_ratio == 0.0) | (second_ratio == 0.0)
                upper_ratio[impossible] = 0.0
                excess = np.maximum(upper_ratio - independent_ratio, 0.0)
                positive_excess[shell] += np.sum(
                    weights * excess, dtype=np.longdouble
                )

    expected_pair_mass = (
        (1 - math.ldexp(1.0, -message_bits))
        * (1 - math.ldexp(2.0, -message_bits))
    )
    rows = []
    for shell in shell_weights:
        one, _correlated = ratio_tables[shell]
        message_weights = np.asarray(
            [
                math.ldexp(float(math.comb(message_bits, weight)), -message_bits)
                for weight in range(message_bits + 1)
            ],
            dtype=np.float64,
        )
        # The zero message is excluded.  For positive shells its ratio is zero.
        mean_ratio = np.sum(message_weights[1:] * one[1:], dtype=np.longdouble)
        random_mean = math.ldexp(
            float(math.comb(output_bits, shell)), message_bits - output_bits
        )
        variance_ratio = 1.0 + (
            random_mean * float(positive_excess[shell]) / float(mean_ratio)
        )
        rows.append(
            {
                "shell_weight": shell,
                "mean_ratio_to_random": float(mean_ratio),
                "positive_pair_excess_normalized": float(positive_excess[shell]),
                "variance_to_mean_upper_diagnostic": variance_ratio,
                "passes_factor_512_diagnostic": variance_ratio <= 512,
            }
        )
    if abs(float(pair_mass) - expected_pair_mass) > 2e-13:
        raise AssertionError("ordered-pair mass check failed")
    return rows, pair_types, dominant_counts, float(pair_mass)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--right-degree", type=int, default=33)
    parser.add_argument("--shell-weight", type=int, action="append", default=[])
    parser.add_argument("--precision", type=int, default=180)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    shells = sorted(set(args.shell_weight or [42, 54, 65, 70, 79]))
    rows, pair_types, dominant_counts, pair_mass = pair_bound(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        shells,
        args.precision,
    )
    payload = {
        "schema": "pure-ea-dominant-character-likelihood-probe-v1",
        "status": "NONDIRECTED_BINARY64_DIAGNOSTIC",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "polynomial_precision_decimal_digits": args.precision,
        },
        "ordered_pair_types": pair_types,
        "normalized_ordered_pair_mass": pair_mass,
        "dominant_character_type_counts": {
            "first_message": dominant_counts[0],
            "second_message": dominant_counts[1],
            "difference_message": dominant_counts[2],
        },
        "shells": rows,
        "proved_reduction": [
            "The reference row law preserves the largest of beta(x), beta(y), and beta(x+y) in absolute value.",
            "The pointwise n-row likelihood ratio is at most (1+delta)^n for the displayed delta.",
            "Replacing each pairwise excess over the product of marginals by its positive part gives a valid variance upper bound.",
        ],
        "scope": [
            "Polynomial ratios use high-precision arithmetic before conversion to binary64.",
            "Pair aggregation, likelihood ratios, and the final bound are nondirected binary64 diagnostics.",
            "A passing row is evidence for an outward implementation, not a certificate.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
