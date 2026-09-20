#!/usr/bin/env python3
"""Verify the dominant-character deviation bound in exact small models.

The verifier checks every ordered pair of distinct nonzero messages and every
output shell.  It checks both the likelihood domination and the final
cancellation-free covariance-excess inequality in exact rational arithmetic.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from verify_pair_kernel_small import (
    pair_accumulator_diagonal_shells,
    pair_activation,
    pair_type,
    single_accumulator_shells,
    single_activation,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "dominant_character_deviation_small_exact.json"


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def reference_pair_shells(
    bias: Fraction, output_bits: int
) -> tuple[Fraction, ...]:
    probabilities = tuple(
        Fraction(1, 4)
        * (1 + (-1) ** (((symbol >> 1) & 1) ^ (symbol & 1)) * bias)
        for symbol in range(4)
    )
    return pair_accumulator_diagonal_shells(probabilities, output_bits)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--output-bits", type=int, default=8)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not (1 <= args.right_degree <= args.message_bits):
        parser.error("right degree must lie in [1,message_bits]")

    denominator = math.comb(args.message_bits, args.right_degree)
    biases = []
    one_shells = []
    for message in range(1 << args.message_bits):
        weight = message.bit_count()
        activation = single_activation(
            args.message_bits, args.right_degree, weight
        )
        biases.append(1 - 2 * activation)
        one_shells.append(single_accumulator_shells(activation, args.output_bits))

    reference_by_bias = {
        bias: reference_pair_shells(bias, args.output_bits)
        for bias in set(biases)
    }
    random_probability = tuple(
        Fraction(math.comb(args.output_bits, shell), 1 << args.output_bits)
        for shell in range(args.output_bits + 1)
    )
    likelihood_checks = 0
    covariance_checks = 0
    exact_positive_sum = [Fraction(0) for _ in random_probability]
    bound_positive_sum = [Fraction(0) for _ in random_probability]

    for first in range(1, 1 << args.message_bits):
        for second in range(1, 1 << args.message_bits):
            if first == second:
                continue
            difference = first ^ second
            triple = (biases[first], biases[second], biases[difference])
            dominant = max(range(3), key=lambda index: abs(triple[index]))
            joint_shells = pair_accumulator_diagonal_shells(
                pair_activation(
                    args.message_bits,
                    args.right_degree,
                    pair_type(first, second, args.message_bits),
                ),
                args.output_bits,
            )
            deterministic_marginal = (
                first.bit_count() == args.message_bits
                or second.bit_count() == args.message_bits
            )
            if deterministic_marginal:
                factor = None
            elif dominant == 2 and difference.bit_count() == args.message_bits:
                factor = (1 + abs(triple[0])) ** args.output_bits
            else:
                retained = abs(triple[dominant])
                residual = sum(
                    abs(triple[index]) for index in range(3) if index != dominant
                )
                if retained == 1:
                    raise AssertionError("unhandled unit retained bias")
                factor = (1 + residual / (1 - retained)) ** args.output_bits

            for shell, random in enumerate(random_probability):
                first_probability = one_shells[first][shell]
                second_probability = one_shells[second][shell]
                joint_probability = joint_shells[shell]
                covariance = (
                    joint_probability - first_probability * second_probability
                )
                exact_positive_sum[shell] += max(Fraction(0), covariance)
                if deterministic_marginal:
                    if covariance != 0:
                        raise AssertionError("deterministic marginal has covariance")
                    candidate = Fraction(0)
                else:
                    if dominant == 0:
                        reference = first_probability * random
                    elif dominant == 1:
                        reference = second_probability * random
                    else:
                        reference = reference_by_bias[triple[2]][shell]
                    if joint_probability > reference * factor:
                        raise AssertionError("likelihood domination failed")
                    likelihood_checks += 1

                    if random == 0:
                        if joint_probability != 0:
                            raise AssertionError("zero random shell has positive mass")
                        candidate = Fraction(0)
                    else:
                        ux = first_probability / random
                        uy = second_probability / random
                        dx = abs(ux - 1)
                        dy = abs(uy - 1)
                        ell = factor - 1
                        if dominant == 0:
                            candidate = dy + dx * dy + ell + dx * ell
                        elif dominant == 1:
                            candidate = dx + dy * dx + ell + dy * ell
                        else:
                            bz = reference_by_bias[triple[2]][shell] / (random * random)
                            dz = abs(bz - 1)
                            candidate = (
                                dz + dx + dy + dx * dy + ell + dz * ell
                            )
                        exact_excess = joint_probability / (random * random) - ux * uy
                        if exact_excess > candidate:
                            raise AssertionError("covariance-excess bound failed")
                        covariance_checks += 1
                if random:
                    bound_positive_sum[shell] += candidate * random * random

    shell_rows = []
    for shell, (exact, bound) in enumerate(
        zip(exact_positive_sum, bound_positive_sum, strict=True)
    ):
        shell_rows.append(
            {
                "shell_weight": shell,
                "exact_positive_covariance_sum": fraction_text(exact),
                "deviation_bound_sum": fraction_text(bound),
                "bound_dominates": bound >= exact,
            }
        )
    if not all(row["bound_dominates"] for row in shell_rows):
        raise AssertionError("aggregate deviation bound failed")

    payload = {
        "schema": "pure-ea-dominant-character-deviation-small-exact-v1",
        "status": "EXACT_RATIONAL_SMALL_MODEL",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "row_neighborhood_count": denominator,
        },
        "verified": {
            "likelihood_pair_shell_checks": likelihood_checks,
            "covariance_excess_pair_shell_checks": covariance_checks,
            "all_shell_aggregates_dominated": True,
        },
        "shells": shell_rows,
        "scope": [
            "The receipt proves the displayed small finite model exactly.",
            "The general inequality uses the same algebra but still requires a separate outward target computation.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["verified"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
