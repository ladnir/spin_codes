#!/usr/bin/env python3
"""Compute exact shell variances for small pure EA chains.

The first sparse map is rectangular.  Later maps are square.  Every output
coordinate independently samples a fixed-size input neighborhood.  This
script propagates the complete ordered-pair type measure with rational
arithmetic and reports the exact shell variance after every stage.

The default instance is intentionally small.  It validates multistage moment
composition but does not extrapolate to the length-512 target.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from verify_pair_kernel_small import compositions, multinomial, pair_activation


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "multistage_pair_moments_K4_B8_r3_3_3_exact.json"
TYPE4 = tuple[int, int, int, int]


def increment_type(pair_type: TYPE4, symbol: int) -> TYPE4:
    values = list(pair_type)
    values[symbol] += 1
    return tuple(values)  # type: ignore[return-value]


@lru_cache(maxsize=None)
def stage_row(
    input_length: int,
    output_length: int,
    right_degree: int,
    input_type: TYPE4,
) -> tuple[tuple[TYPE4, Fraction], ...]:
    activation = pair_activation(input_length, right_degree, input_type)
    state: dict[tuple[int, TYPE4], Fraction] = {
        (0, (0, 0, 0, 0)): Fraction(1)
    }
    for _ in range(output_length):
        following: defaultdict[tuple[int, TYPE4], Fraction] = defaultdict(Fraction)
        for (acc_symbol, output_type), mass in state.items():
            for input_symbol, probability in enumerate(activation):
                next_symbol = acc_symbol ^ input_symbol
                following[(next_symbol, increment_type(output_type, next_symbol))] += (
                    mass * probability
                )
        state = dict(following)
    result: defaultdict[TYPE4, Fraction] = defaultdict(Fraction)
    for (_acc_symbol, output_type), mass in state.items():
        result[output_type] += mass
    if sum(result.values()) != 1:
        raise AssertionError("stage row does not preserve probability mass")
    return tuple(sorted(result.items()))


def initial_pair_measure(message_bits: int) -> dict[TYPE4, Fraction]:
    result = {}
    for input_type in compositions(message_bits, 4):
        n00, n01, n10, n11 = input_type
        first_weight = n10 + n11
        second_weight = n01 + n11
        difference_weight = n01 + n10
        if min(first_weight, second_weight, difference_weight) == 0:
            continue
        result[input_type] = Fraction(multinomial(input_type))
    expected = ((1 << message_bits) - 1) * ((1 << message_bits) - 2)
    if sum(result.values()) != expected:
        raise AssertionError("initial ordered-pair measure has the wrong mass")
    return result


def apply_stage(
    measure: dict[TYPE4, Fraction],
    input_length: int,
    output_length: int,
    right_degree: int,
) -> dict[TYPE4, Fraction]:
    result: defaultdict[TYPE4, Fraction] = defaultdict(Fraction)
    for input_type, count in measure.items():
        for output_type, probability in stage_row(
            input_length, output_length, right_degree, input_type
        ):
            result[output_type] += count * probability
    if sum(result.values()) != sum(measure.values()):
        raise AssertionError("pair measure changed mass")
    return dict(result)


def shell_statistics(
    measure: dict[TYPE4, Fraction], message_bits: int, output_bits: int
) -> tuple[list[dict[str, object]], Fraction, int]:
    partners = (1 << message_bits) - 2
    first_marginal: defaultdict[int, Fraction] = defaultdict(Fraction)
    factorial_diagonal: defaultdict[int, Fraction] = defaultdict(Fraction)
    for pair_type, mass in measure.items():
        _n00, n01, n10, n11 = pair_type
        first_weight = n10 + n11
        second_weight = n01 + n11
        first_marginal[first_weight] += mass
        if first_weight == second_weight:
            factorial_diagonal[first_weight] += mass

    rows = []
    maximum_ratio = Fraction(0)
    maximum_weight = 0
    for weight in range(output_bits + 1):
        mean = first_marginal[weight] / partners
        factorial = factorial_diagonal[weight]
        variance = mean + factorial - mean * mean
        ratio = None if mean == 0 else variance / mean
        if weight > 0 and ratio is not None and ratio > maximum_ratio:
            maximum_ratio = ratio
            maximum_weight = weight
        rows.append(
            {
                "weight": weight,
                "mean": f"{mean.numerator}/{mean.denominator}",
                "second_factorial_moment": (
                    f"{factorial.numerator}/{factorial.denominator}"
                ),
                "variance": f"{variance.numerator}/{variance.denominator}",
                "variance_to_mean": (
                    None
                    if ratio is None
                    else f"{ratio.numerator}/{ratio.denominator}"
                ),
                "variance_to_mean_float": None if ratio is None else float(ratio),
            }
        )
    return rows, maximum_ratio, maximum_weight


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--output-bits", type=int, default=8)
    parser.add_argument("--degrees", type=int, nargs="+", default=[3, 3, 3])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.message_bits < 2 or args.output_bits < 1:
        parser.error("message and output lengths must be positive")
    if not args.degrees:
        parser.error("at least one stage is required")
    if not (1 <= args.degrees[0] <= args.message_bits):
        parser.error("first degree must lie in [1,message bits]")
    if any(not 1 <= degree <= args.output_bits for degree in args.degrees[1:]):
        parser.error("later degrees must lie in [1,output bits]")

    measure = initial_pair_measure(args.message_bits)
    stages = []
    input_length = args.message_bits
    for stage_index, degree in enumerate(args.degrees, start=1):
        measure = apply_stage(measure, input_length, args.output_bits, degree)
        input_length = args.output_bits
        shells, maximum_ratio, maximum_weight = shell_statistics(
            measure, args.message_bits, args.output_bits
        )
        stages.append(
            {
                "stage": stage_index,
                "right_degree": degree,
                "pair_type_support": len(measure),
                "maximum_positive_shell_variance_to_mean": (
                    f"{maximum_ratio.numerator}/{maximum_ratio.denominator}"
                ),
                "maximum_positive_shell_variance_to_mean_float": float(
                    maximum_ratio
                ),
                "maximum_positive_shell_weight": maximum_weight,
                "shells": shells,
            }
        )

    payload = {
        "schema": "pure-ea-multistage-pair-moments-small-v1",
        "status": "EXACT_RATIONAL_SMALL_MODEL",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degrees": args.degrees,
        },
        "ordered_pair_mass": ((1 << args.message_bits) - 1)
        * ((1 << args.message_bits) - 2),
        "stages": stages,
        "scope": [
            "Every reported moment is exact over the stated finite ensemble.",
            "The result does not extrapolate to the length-512 target.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            [
                {
                    "stage": stage["stage"],
                    "right_degree": stage["right_degree"],
                    "maximum_positive_shell_variance_to_mean_float": stage[
                        "maximum_positive_shell_variance_to_mean_float"
                    ],
                    "maximum_positive_shell_weight": stage[
                        "maximum_positive_shell_weight"
                    ],
                }
                for stage in stages
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
