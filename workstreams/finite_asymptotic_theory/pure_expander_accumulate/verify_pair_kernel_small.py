#!/usr/bin/env python3
"""Verify the exact one-stage sparse-mixer--accumulator pair kernel.

The sparse map has ``message_bits`` left coordinates and ``output_bits``
right coordinates.  Each right coordinate independently samples a uniform
``right_degree``-subset of the left coordinates and outputs its parity.

The verifier compares the analytic one-word and ordered-pair shell moments
with exhaustive enumeration of every sparse map.  Defaults are deliberately
small.  This script proves only the finite identities at the requested small
parameters; it does not extrapolate them to the length-512 target.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "pair_kernel_small_exact.json"
TYPE4 = tuple[int, int, int, int]


def accumulator(word: int, length: int) -> int:
    state = 0
    output = 0
    for index in range(length):
        state ^= (word >> index) & 1
        output |= state << index
    return output


def krawtchouk(length: int, degree: int, weight: int) -> int:
    return sum(
        (-1) ** j
        * math.comb(weight, j)
        * math.comb(length - weight, degree - j)
        for j in range(max(0, degree - (length - weight)), min(degree, weight) + 1)
    )


def pair_type(first: int, second: int, length: int) -> TYPE4:
    counts = [0, 0, 0, 0]
    for index in range(length):
        a = (first >> index) & 1
        b = (second >> index) & 1
        counts[(a << 1) | b] += 1
    return tuple(counts)  # type: ignore[return-value]


def multinomial(parts: Iterable[int]) -> int:
    remaining = sum(parts)
    result = 1
    for part in parts:
        result *= math.comb(remaining, part)
        remaining -= part
    return result


def compositions(total: int, parts: int) -> Iterable[tuple[int, ...]]:
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first, *rest)


def single_activation(message_bits: int, right_degree: int, weight: int) -> Fraction:
    denominator = math.comb(message_bits, right_degree)
    bias = Fraction(krawtchouk(message_bits, right_degree, weight), denominator)
    return (1 - bias) / 2


def pair_activation(
    message_bits: int, right_degree: int, input_type: TYPE4
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    _n00, n01, n10, n11 = input_type
    first_weight = n10 + n11
    second_weight = n01 + n11
    difference_weight = n01 + n10
    denominator = math.comb(message_bits, right_degree)
    first_bias = Fraction(
        krawtchouk(message_bits, right_degree, first_weight), denominator
    )
    second_bias = Fraction(
        krawtchouk(message_bits, right_degree, second_weight), denominator
    )
    difference_bias = Fraction(
        krawtchouk(message_bits, right_degree, difference_weight), denominator
    )
    probabilities = []
    for symbol in range(4):
        first = (symbol >> 1) & 1
        second = symbol & 1
        probability = (
            1
            + (-1) ** first * first_bias
            + (-1) ** second * second_bias
            + (-1) ** (first ^ second) * difference_bias
        ) / 4
        if probability < 0:
            raise AssertionError("pair activation probability is negative")
        probabilities.append(probability)
    if sum(probabilities) != 1:
        raise AssertionError("pair activation probabilities do not sum to one")
    return tuple(probabilities)  # type: ignore[return-value]


def single_accumulator_shells(
    activation: Fraction, output_bits: int
) -> tuple[Fraction, ...]:
    state: dict[tuple[int, int], Fraction] = {(0, 0): Fraction(1)}
    for _ in range(output_bits):
        following: defaultdict[tuple[int, int], Fraction] = defaultdict(Fraction)
        for (acc_state, weight), mass in state.items():
            for input_bit, probability in ((0, 1 - activation), (1, activation)):
                next_state = acc_state ^ input_bit
                following[(next_state, weight + next_state)] += mass * probability
        state = dict(following)
    shells = [Fraction(0) for _ in range(output_bits + 1)]
    for (_acc_state, weight), mass in state.items():
        shells[weight] += mass
    if sum(shells) != 1:
        raise AssertionError("single accumulator law does not preserve mass")
    return tuple(shells)


def pair_accumulator_diagonal_shells(
    activation: tuple[Fraction, Fraction, Fraction, Fraction], output_bits: int
) -> tuple[Fraction, ...]:
    state: dict[tuple[int, int, int], Fraction] = {(0, 0, 0): Fraction(1)}
    for _ in range(output_bits):
        following: defaultdict[tuple[int, int, int], Fraction] = defaultdict(Fraction)
        for (acc_symbol, first_weight, second_weight), mass in state.items():
            for input_symbol, probability in enumerate(activation):
                next_symbol = acc_symbol ^ input_symbol
                next_first = (next_symbol >> 1) & 1
                next_second = next_symbol & 1
                following[
                    (
                        next_symbol,
                        first_weight + next_first,
                        second_weight + next_second,
                    )
                ] += mass * probability
        state = dict(following)
    diagonal = [Fraction(0) for _ in range(output_bits + 1)]
    for (_acc_symbol, first_weight, second_weight), mass in state.items():
        if first_weight == second_weight:
            diagonal[first_weight] += mass
    return tuple(diagonal)


def analytic_moments(
    message_bits: int, output_bits: int, right_degree: int
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    mean = [Fraction(0) for _ in range(output_bits + 1)]
    for weight in range(1, message_bits + 1):
        activation = single_activation(message_bits, right_degree, weight)
        shells = single_accumulator_shells(activation, output_bits)
        count = math.comb(message_bits, weight)
        for output_weight, probability in enumerate(shells):
            mean[output_weight] += count * probability

    factorial = [Fraction(0) for _ in range(output_bits + 1)]
    for input_type in compositions(message_bits, 4):
        n00, n01, n10, n11 = input_type
        first_weight = n10 + n11
        second_weight = n01 + n11
        difference_weight = n01 + n10
        if min(first_weight, second_weight, difference_weight) == 0:
            continue
        activation = pair_activation(message_bits, right_degree, input_type)
        diagonal = pair_accumulator_diagonal_shells(activation, output_bits)
        count = multinomial(input_type)
        for output_weight, probability in enumerate(diagonal):
            factorial[output_weight] += count * probability
    return tuple(mean), tuple(factorial)


def apply_sparse_map(message: int, neighborhoods: tuple[int, ...]) -> int:
    output = 0
    for index, mask in enumerate(neighborhoods):
        output |= ((message & mask).bit_count() & 1) << index
    return output


def exhaustive_moments(
    message_bits: int, output_bits: int, right_degree: int
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...], int]:
    masks = tuple(
        sum(1 << index for index in support)
        for support in itertools.combinations(range(message_bits), right_degree)
    )
    map_count = len(masks) ** output_bits
    mean_total = [0 for _ in range(output_bits + 1)]
    factorial_total = [0 for _ in range(output_bits + 1)]
    for neighborhoods in itertools.product(masks, repeat=output_bits):
        shell_counts: Counter[int] = Counter()
        for message in range(1, 1 << message_bits):
            intermediate = apply_sparse_map(message, neighborhoods)
            encoded = accumulator(intermediate, output_bits)
            shell_counts[encoded.bit_count()] += 1
        for weight in range(output_bits + 1):
            count = shell_counts[weight]
            mean_total[weight] += count
            factorial_total[weight] += count * (count - 1)
    return (
        tuple(Fraction(value, map_count) for value in mean_total),
        tuple(Fraction(value, map_count) for value in factorial_total),
        map_count,
    )


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--output-bits", type=int, default=5)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not (1 <= args.right_degree <= args.message_bits):
        parser.error("right degree must lie in [1,message_bits]")
    if args.output_bits < 1:
        parser.error("output bits must be positive")

    analytic_mean, analytic_factorial = analytic_moments(
        args.message_bits, args.output_bits, args.right_degree
    )
    exhaustive_mean, exhaustive_factorial, map_count = exhaustive_moments(
        args.message_bits, args.output_bits, args.right_degree
    )
    if analytic_mean != exhaustive_mean:
        raise AssertionError("analytic and exhaustive shell means differ")
    if analytic_factorial != exhaustive_factorial:
        raise AssertionError("analytic and exhaustive factorial moments differ")

    shells = []
    maximum_ratio = Fraction(0)
    maximum_weight = 0
    for weight, (mean, factorial) in enumerate(
        zip(analytic_mean, analytic_factorial, strict=True)
    ):
        variance = mean + factorial - mean * mean
        ratio = None if mean == 0 else variance / mean
        if weight > 0 and ratio is not None and ratio > maximum_ratio:
            maximum_ratio = ratio
            maximum_weight = weight
        shells.append(
            {
                "weight": weight,
                "mean": fraction_text(mean),
                "second_factorial_moment": fraction_text(factorial),
                "variance": fraction_text(variance),
                "variance_to_mean": None if ratio is None else fraction_text(ratio),
            }
        )

    payload = {
        "schema": "pure-ea-independent-right-regular-pair-kernel-small-v1",
        "status": "EXACT_RATIONAL_SMALL_MODEL",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "enumerated_maps": map_count,
        },
        "verified_identities": [
            "analytic shell means equal exhaustive ensemble means",
            "analytic shell second factorial moments equal exhaustive ensemble moments",
        ],
        "maximum_positive_shell_variance_to_mean": fraction_text(maximum_ratio),
        "maximum_positive_shell_weight": maximum_weight,
        "shells": shells,
        "scope": [
            "The receipt proves only the displayed small finite instance.",
            "It does not prove rank, spectrum caps, or distance at length 512.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
