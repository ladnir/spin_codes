#!/usr/bin/env python3
"""Verify the sector-zero Hoeffding expansion in exact rational arithmetic."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from verify_pair_kernel_small import krawtchouk
from verify_radial_triple_factorization_small import radial_walk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "sector_zero_hoeffding_K4_B8_r3_exact.json"


def coefficient(n: int, j: int, p: int, q: Fraction) -> Fraction:
    """Return [z^p](1-z)^j(1-q+qz)^(n-j)."""
    result = Fraction(0)
    low = max(0, p - (n - j))
    high = min(j, p)
    for selected in range(low, high + 1):
        remainder = p - selected
        result += (
            (-1) ** selected
            * math.comb(j, selected)
            * math.comb(n - j, remainder)
            * q**remainder
            * (1 - q) ** (n - j - remainder)
        )
    return result


def joint_shell_probability(
    n: int, p: int, probability: tuple[Fraction, Fraction, Fraction, Fraction]
) -> Fraction:
    a, b, c, d = probability
    result = Fraction(0)
    low = max(0, 2 * p - n)
    for intersection in range(low, p + 1):
        cross = p - intersection
        result += (
            math.comb(p, intersection)
            * math.comb(n - p, cross)
            * a ** (n - 2 * p + intersection)
            * (b * c) ** cross
            * d**intersection
        )
    return math.comb(n, p) * result


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--output-bits", type=int, default=8)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.message_bits > 8 or args.output_bits > 16:
        parser.error("the exact verifier is limited to small parameters")

    denominator = math.comb(args.message_bits, args.right_degree)
    biases = [
        Fraction(
            krawtchouk(args.message_bits, args.right_degree, weight),
            denominator,
        )
        for weight in range(args.message_bits + 1)
    ]
    checks = 0
    aggregate_orders = {
        (p, j): Fraction(0)
        for p in range(args.output_bits + 1)
        for j in range(1, args.output_bits + 1)
    }
    samples = []
    for first in range(1 << args.message_bits):
        for second in range(1 << args.message_bits):
            b1 = biases[first.bit_count()]
            b2 = biases[second.bit_count()]
            b3 = biases[(first ^ second).bit_count()]
            probability = (
                (1 + b1 + b2 + b3) / 4,
                (1 + b1 - b2 - b3) / 4,
                (1 - b1 + b2 - b3) / 4,
                (1 - b1 - b2 + b3) / 4,
            )
            q1 = (1 - b1) / 2
            q2 = (1 - b2) / 2
            determinant = (b3 - b1 * b2) / 4
            for p in range(args.output_bits + 1):
                marginal = (
                    math.comb(args.output_bits, p) ** 2
                    * q1**p
                    * (1 - q1) ** (args.output_bits - p)
                    * q2**p
                    * (1 - q2) ** (args.output_bits - p)
                )
                direct = joint_shell_probability(
                    args.output_bits, p, probability
                ) - marginal
                terms = []
                for j in range(1, args.output_bits + 1):
                    term = (
                        math.comb(args.output_bits, j)
                        * determinant**j
                        * coefficient(args.output_bits, j, p, q1)
                        * coefficient(args.output_bits, j, p, q2)
                    )
                    terms.append(term)
                    aggregate_orders[(p, j)] += term
                expanded = sum(terms, Fraction(0))
                if expanded != direct:
                    raise AssertionError("Hoeffding expansion mismatch")
                checks += 1
                if len(samples) < 10:
                    samples.append(
                        {
                            "messages": [first, second],
                            "shell": p,
                            "covariance": fraction_text(direct),
                        }
                    )

    negative = [
        (p, j, value)
        for (p, j), value in aggregate_orders.items()
        if value < 0
    ]
    if negative:
        raise AssertionError(f"negative fixed-order aggregate: {negative[0]}")
    walk = radial_walk(
        args.message_bits, args.right_degree, args.output_bits
    )
    radial_checks = 0
    polynomial_transform_checks = 0
    for p in range(args.output_bits + 1):
        for j in range(1, args.output_bits + 1):
            terms = []
            for convolution_power in range(j + 1):
                radial_function = [
                    coefficient(
                        args.output_bits,
                        j,
                        p,
                        (1 - biases[weight]) / 2,
                    )
                    * biases[weight] ** (j - convolution_power)
                    for weight in range(args.message_bits + 1)
                ]
                transform = [
                    sum(
                        krawtchouk(args.message_bits, weight, frequency_weight)
                        * radial_function[weight]
                        for weight in range(args.message_bits + 1)
                    )
                    for frequency_weight in range(args.message_bits + 1)
                ]
                polynomial_transform = [
                    Fraction(2**args.message_bits, 2 ** (args.output_bits - j))
                    * sum(
                        math.comb(args.output_bits - j, power)
                        * krawtchouk(
                            args.output_bits,
                            p,
                            j + power,
                        )
                        * walk[power + j - convolution_power][
                            frequency_weight
                        ]
                        for power in range(args.output_bits - j + 1)
                    )
                    for frequency_weight in range(args.message_bits + 1)
                ]
                if polynomial_transform != transform:
                    raise AssertionError(
                        "polynomial radial transform mismatch "
                        f"p={p} j={j} ell={convolution_power}"
                    )
                polynomial_transform_checks += len(transform)
                convolution_quadratic = sum(
                    math.comb(args.message_bits, frequency_weight)
                    * walk[convolution_power][frequency_weight]
                    * transform[frequency_weight] ** 2
                    for frequency_weight in range(args.message_bits + 1)
                )
                terms.append(
                    (-1) ** (j - convolution_power)
                    * math.comb(j, convolution_power)
                    * convolution_quadratic
                )
            radial = (
                math.comb(args.output_bits, j)
                * sum(terms, Fraction(0))
                / 4**j
            )
            if radial != aggregate_orders[(p, j)]:
                raise AssertionError(
                    "radial fixed-order expansion mismatch "
                    f"p={p} j={j} radial={radial} "
                    f"direct={aggregate_orders[(p, j)]}"
                )
            radial_checks += 1
    payload = {
        "schema": "pure-ea-sector-zero-hoeffding-small-v1",
        "status": "EXACT_RATIONAL_SMALL_MODEL",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
        },
        "checked_pair_shell_identities": checks,
        "checked_nonnegative_order_aggregates": len(aggregate_orders),
        "checked_radial_order_aggregates": radial_checks,
        "checked_polynomial_radial_transforms": polynomial_transform_checks,
        "samples": samples,
        "verified": [
            "The exact pair covariance equals the fixed-order Hoeffding expansion.",
            "Every fixed-order aggregate over all ordered message pairs is nonnegative in the checked model.",
            "Every fixed-order aggregate equals the rank-(k+1) radial transform formula.",
            "Every radial transform equals the Krawtchouk-polynomial expansion through the row-XOR walk.",
        ],
        "scope": [
            "All arithmetic and comparisons in this receipt are exact.",
            "General fixed-order nonnegativity follows separately from the Schur product theorem.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
