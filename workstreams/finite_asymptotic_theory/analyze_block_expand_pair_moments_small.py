#!/usr/bin/env python3
"""Compute exact shell variances for a small Block Expand-t ensemble.

The checked instance has K=4, B=8, and four independent regions of length
two.  It uses the exact regional pair EGF and the exact common-permutation
accumulator pair kernel.  The calculation is a small-model diagnostic, not a
length-512 theorem.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from verify_accumulator_pair_type_kernel import kernel_counts, multinomial
from verify_block_expand_pair_region_small import compositions, egf_count


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "block_expand_pair_moments_K4_B8_exact.json"
K = 4
B = 8
REGIONS = (2, 2, 2, 2)
TYPE4 = tuple[int, int, int, int]


def add_types(first: TYPE4, second: TYPE4) -> TYPE4:
    return tuple(first[index] + second[index] for index in range(4))  # type: ignore[return-value]


@lru_cache(maxsize=None)
def regional_distribution(length: int, pair_type: TYPE4) -> dict[TYPE4, Fraction]:
    inputs = (pair_type[2], pair_type[1], pair_type[3])
    denominator = length ** sum(inputs)
    return {
        output: Fraction(count, denominator)
        for output in compositions(length, 4)
        if (count := egf_count(length, inputs, output))
    }


def expander_pair_measure() -> dict[TYPE4, Fraction]:
    """Return expected ordered distinct-message counts by output pair type."""
    result: defaultdict[TYPE4, Fraction] = defaultdict(Fraction)
    for input_type in compositions(K, 4):
        n00, n01, n10, n11 = input_type
        if not (n10 + n11 and n01 + n11 and n10 + n01):
            continue
        distribution: dict[TYPE4, Fraction] = {(0, 0, 0, 0): Fraction(1)}
        for length in REGIONS:
            following: defaultdict[TYPE4, Fraction] = defaultdict(Fraction)
            for prefix_type, prefix_mass in distribution.items():
                for region_type, probability in regional_distribution(
                    length, input_type
                ).items():
                    following[add_types(prefix_type, region_type)] += (
                        prefix_mass * probability
                    )
            distribution = dict(following)
        input_pairs = multinomial(input_type)
        for output_type, probability in distribution.items():
            result[output_type] += input_pairs * probability
    expected_pairs = ((1 << K) - 1) * ((1 << K) - 2)
    if sum(result.values()) != expected_pairs:
        raise AssertionError("expander pair measure has the wrong mass")
    return dict(result)


@lru_cache(maxsize=None)
def accumulator_row(input_type: TYPE4) -> dict[TYPE4, Fraction]:
    denominator = multinomial(input_type)
    return {
        output_type: Fraction(count, denominator)
        for output_type, count in kernel_counts(input_type).items()
    }


def apply_accumulator(measure: dict[TYPE4, Fraction]) -> dict[TYPE4, Fraction]:
    result: defaultdict[TYPE4, Fraction] = defaultdict(Fraction)
    for input_type, mass in measure.items():
        for output_type, probability in accumulator_row(input_type).items():
            result[output_type] += mass * probability
    if sum(result.values()) != sum(measure.values()):
        raise AssertionError("accumulator pair kernel did not preserve mass")
    return dict(result)


def shell_statistics(measure: dict[TYPE4, Fraction]) -> list[dict[str, object]]:
    messages_other_than_pair_member = (1 << K) - 2
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
    for weight in range(B + 1):
        mean = first_marginal[weight] / messages_other_than_pair_member
        factorial = factorial_diagonal[weight]
        variance = mean + factorial - mean * mean
        ratio = None if mean == 0 else variance / mean
        rows.append(
            {
                "weight": weight,
                "mean": f"{mean.numerator}/{mean.denominator}",
                "second_factorial_moment": f"{factorial.numerator}/{factorial.denominator}",
                "variance": f"{variance.numerator}/{variance.denominator}",
                "variance_to_mean": None if ratio is None else float(ratio),
            }
        )
    return rows


def main() -> None:
    global K, B, REGIONS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=K)
    parser.add_argument("--region-lengths", type=int, nargs="+", default=REGIONS)
    parser.add_argument("--stages", type=int, default=5)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    K = args.message_bits
    REGIONS = tuple(args.region_lengths)
    B = sum(REGIONS)
    if K <= 0 or B <= 0 or args.stages < 0:
        parser.error("message bits, output bits, and stage count must be valid")
    regional_distribution.cache_clear()
    accumulator_row.cache_clear()

    measure = expander_pair_measure()
    stages = []
    expected_mass = ((1 << K) - 1) * ((1 << K) - 2)
    for stage in range(args.stages + 1):
        rows = shell_statistics(measure)
        positive_ratios = [
            (row["variance_to_mean"], row["weight"])
            for row in rows
            if row["variance_to_mean"] is not None and row["weight"] != 0
        ]
        maximum_ratio, maximum_weight = max(positive_ratios)
        stages.append(
            {
                "accumulator_stages": stage,
                "pair_type_support": len(measure),
                "ordered_pair_mass": expected_mass,
                "maximum_positive_shell_variance_to_mean": maximum_ratio,
                "maximum_positive_shell_weight": maximum_weight,
                "shells": rows,
            }
        )
        if stage < args.stages:
            measure = apply_accumulator(measure)

    payload = {
        "schema": "block-expand-pair-moments-small-exact-v1",
        "status": "EXACT_RATIONAL_SMALL_MODEL",
        "parameters": {
            "message_bits": K,
            "output_bits": B,
            "region_lengths": list(REGIONS),
            "left_degree": len(REGIONS),
        },
        "stages": stages,
        "notes": [
            "Every moment is exact over rational numbers.",
            "Ordered input pairs exclude x=0, y=0, and x=y.",
            "The small model does not prove a variance factor at length 512.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            [
                {
                    "accumulator_stages": row["accumulator_stages"],
                    "maximum_positive_shell_variance_to_mean": row[
                        "maximum_positive_shell_variance_to_mean"
                    ],
                    "maximum_positive_shell_weight": row[
                        "maximum_positive_shell_weight"
                    ],
                }
                for row in stages
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
