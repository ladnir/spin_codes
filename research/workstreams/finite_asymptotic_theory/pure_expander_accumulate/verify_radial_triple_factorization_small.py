#!/usr/bin/env python3
"""Verify the radial rank-(k+1) factorization of triple bias moments."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "radial_triple_factorization_K4_r3_m8_exact.json"


def radial_walk(
    message_bits: int, right_degree: int, maximum_steps: int
) -> list[list[Fraction]]:
    """Return per-vector probabilities v_h(m), indexed by m then h."""
    denominator = math.comb(message_bits, right_degree)
    shell = [Fraction(0) for _ in range(message_bits + 1)]
    shell[0] = 1
    result = [[Fraction(int(weight == 0)) for weight in range(message_bits + 1)]]
    for _step in range(maximum_steps):
        following = [Fraction(0) for _ in range(message_bits + 1)]
        for weight, mass in enumerate(shell):
            if mass == 0:
                continue
            low = max(0, right_degree - (message_bits - weight))
            high = min(weight, right_degree)
            for overlap in range(low, high + 1):
                next_weight = weight + right_degree - 2 * overlap
                probability = Fraction(
                    math.comb(weight, overlap)
                    * math.comb(message_bits - weight, right_degree - overlap),
                    denominator,
                )
                following[next_weight] += mass * probability
        if sum(following) != 1:
            raise AssertionError("radial shell walk lost mass")
        shell = following
        result.append(
            [
                shell[weight] / math.comb(message_bits, weight)
                for weight in range(message_bits + 1)
            ]
        )
    return result


def biases(message_bits: int, right_degree: int) -> list[Fraction]:
    denominator = math.comb(message_bits, right_degree)
    return [
        Fraction(krawtchouk(message_bits, right_degree, weight), denominator)
        for weight in range(message_bits + 1)
    ]


def direct_triple(
    message_bits: int,
    row_bias: list[Fraction],
    first_power: int,
    second_power: int,
    difference_power: int,
) -> Fraction:
    total = Fraction(0)
    for first in range(1 << message_bits):
        for second in range(1 << message_bits):
            total += (
                row_bias[first.bit_count()] ** first_power
                * row_bias[second.bit_count()] ** second_power
                * row_bias[(first ^ second).bit_count()] ** difference_power
            )
    return total / (1 << (2 * message_bits))


def factored_triple(
    message_bits: int,
    walk: list[list[Fraction]],
    first_power: int,
    second_power: int,
    difference_power: int,
) -> Fraction:
    return sum(
        math.comb(message_bits, weight)
        * walk[first_power][weight]
        * walk[second_power][weight]
        * walk[difference_power][weight]
        for weight in range(message_bits + 1)
    )


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--maximum-power", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 1 <= args.right_degree <= args.message_bits:
        parser.error("invalid right degree")
    if args.message_bits > 8:
        parser.error("the direct verifier is limited to message_bits <= 8")

    walk = radial_walk(
        args.message_bits, args.right_degree, args.maximum_power
    )
    row_bias = biases(args.message_bits, args.right_degree)
    checks = 0
    samples = []
    for first_power in range(args.maximum_power + 1):
        for second_power in range(args.maximum_power - first_power + 1):
            for difference_power in range(
                args.maximum_power - first_power - second_power + 1
            ):
                direct = direct_triple(
                    args.message_bits,
                    row_bias,
                    first_power,
                    second_power,
                    difference_power,
                )
                factored = factored_triple(
                    args.message_bits,
                    walk,
                    first_power,
                    second_power,
                    difference_power,
                )
                if direct != factored:
                    raise AssertionError("radial triple factorization failed")
                checks += 1
                if len(samples) < 12:
                    samples.append(
                        {
                            "powers": [
                                first_power,
                                second_power,
                                difference_power,
                            ],
                            "value": fraction_text(direct),
                        }
                    )

    payload = {
        "schema": "pure-ea-radial-triple-factorization-small-v1",
        "status": "EXACT_RATIONAL_SMALL_MODEL",
        "parameters": {
            "message_bits": args.message_bits,
            "right_degree": args.right_degree,
            "maximum_total_power": args.maximum_power,
        },
        "checked_triples": checks,
        "sample_values": samples,
        "verified_identity": (
            "E[beta(X)^a beta(Y)^b beta(X+Y)^d] equals "
            "sum_h binom(k,h) v_h(a) v_h(b) v_h(d)"
        ),
        "scope": [
            "Every comparison in this receipt is exact rational arithmetic.",
            "The receipt validates only the displayed finite parameter range.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
