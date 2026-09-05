#!/usr/bin/env python3
"""Verify exact first and second moments of a weighted sparse-EA spectrum.

This is an exhaustive small-model check of the transfer-weighted object

    F_f(C) = sum_{x != 0} f(wt(Cx)).

The analytic calculation uses the exact one-word and ordered-pair kernels.
The exhaustive calculation enumerates every right-regular sparse map.  All
arithmetic is rational.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIR_DIR = HERE.parent / "pure_expander_accumulate"
sys.path.insert(0, str(PAIR_DIR))

from verify_pair_kernel_small import (  # noqa: E402
    accumulator,
    apply_sparse_map,
    compositions,
    multinomial,
    pair_activation,
    single_accumulator_shells,
    single_activation,
)


DEFAULT_OUTPUT = HERE / "combined_functional_K4_B5_r3_exact.json"


def pair_functional_expectation(
    activation: tuple[Fraction, Fraction, Fraction, Fraction],
    output_bits: int,
    weights: tuple[Fraction, ...],
) -> Fraction:
    state: dict[tuple[int, int, int], Fraction] = {(0, 0, 0): Fraction(1)}
    for _ in range(output_bits):
        following: defaultdict[tuple[int, int, int], Fraction] = defaultdict(Fraction)
        for (acc_symbol, first_weight, second_weight), mass in state.items():
            for input_symbol, probability in enumerate(activation):
                next_symbol = acc_symbol ^ input_symbol
                first = (next_symbol >> 1) & 1
                second = next_symbol & 1
                following[
                    (next_symbol, first_weight + first, second_weight + second)
                ] += mass * probability
        state = dict(following)
    return sum(
        mass * weights[first_weight] * weights[second_weight]
        for (_symbol, first_weight, second_weight), mass in state.items()
    )


def analytic_moments(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    weights: tuple[Fraction, ...],
) -> tuple[Fraction, Fraction]:
    first = Fraction(0)
    diagonal_second = Fraction(0)
    for input_weight in range(1, message_bits + 1):
        activation = single_activation(message_bits, right_degree, input_weight)
        shells = single_accumulator_shells(activation, output_bits)
        count = math.comb(message_bits, input_weight)
        first += count * sum(p * f for p, f in zip(shells, weights, strict=True))
        diagonal_second += count * sum(
            p * f * f for p, f in zip(shells, weights, strict=True)
        )

    off_diagonal_second = Fraction(0)
    for input_type in compositions(message_bits, 4):
        _n00, n01, n10, n11 = input_type
        h1 = n10 + n11
        h2 = n01 + n11
        h3 = n01 + n10
        if min(h1, h2, h3) == 0:
            continue
        activation = pair_activation(message_bits, right_degree, input_type)
        off_diagonal_second += multinomial(input_type) * pair_functional_expectation(
            activation, output_bits, weights
        )
    return first, diagonal_second + off_diagonal_second


def exhaustive_moments(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    weights: tuple[Fraction, ...],
) -> tuple[Fraction, Fraction, int]:
    masks = tuple(
        sum(1 << index for index in support)
        for support in itertools.combinations(range(message_bits), right_degree)
    )
    map_count = len(masks) ** output_bits
    first = Fraction(0)
    second = Fraction(0)
    for neighborhoods in itertools.product(masks, repeat=output_bits):
        functional = Fraction(0)
        for message in range(1, 1 << message_bits):
            intermediate = apply_sparse_map(message, neighborhoods)
            encoded = accumulator(intermediate, output_bits)
            functional += weights[encoded.bit_count()]
        first += functional
        second += functional * functional
    return first / map_count, second / map_count, map_count


def text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--output-bits", type=int, default=5)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 1 <= args.right_degree <= args.message_bits:
        parser.error("right degree must lie in [1,message_bits]")
    if args.output_bits < 1:
        parser.error("output bits must be positive")

    # A non-exponential positive vector ensures that the verifier checks the
    # full bivariate weight law rather than one accidental matrix identity.
    weights = tuple(
        Fraction(weight + 1, 2 ** (weight + 1))
        for weight in range(args.output_bits + 1)
    )
    analytic_first, analytic_second = analytic_moments(
        args.message_bits, args.output_bits, args.right_degree, weights
    )
    exhaustive_first, exhaustive_second, map_count = exhaustive_moments(
        args.message_bits, args.output_bits, args.right_degree, weights
    )
    if analytic_first != exhaustive_first:
        raise AssertionError("analytic and exhaustive first moments differ")
    if analytic_second != exhaustive_second:
        raise AssertionError("analytic and exhaustive second moments differ")
    variance = analytic_second - analytic_first * analytic_first
    payload = {
        "schema": "ea-transfer-weighted-functional-small-v1",
        "status": "EXACT_RATIONAL_SMALL_MODEL",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "enumerated_maps": map_count,
            "functional_weights": [text(value) for value in weights],
        },
        "claim": {
            "analytic_first_moment": text(analytic_first),
            "exhaustive_first_moment": text(exhaustive_first),
            "analytic_second_moment": text(analytic_second),
            "exhaustive_second_moment": text(exhaustive_second),
            "variance": text(variance),
            "identities_match": True,
        },
        "scope": [
            "This proves the displayed small finite instance exactly.",
            "The same pair-kernel formula applies to any fixed nonnegative weight vector f.",
            "It does not extrapolate a numerical moment bound to K=256,B=512.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
