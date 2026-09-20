#!/usr/bin/env python3
"""Exact expected spectrum of the random-layer ParityMix8 outer ensemble.

A length-B, dimension-B/2 codeword starts as (message, 0).  Each layer
uniformly permutes the B coordinates, divides them into consecutive groups of
eight, and applies y_i = x_i + sum_j x_j within every group.  For even group
size this local map is linear, symmetric, and involutory.

For a group containing r ones, the output weight is r when r is even and
8-r when r is odd.  Consequently the exact joint input/output count for one
random layer is the coefficient table of

    (sum_{r=0}^8 binom(8,r) x^r y^f(r))^(B/8).

The script constructs that polynomial over the integers with python-flint and
composes its normalized transition table using exact rational arithmetic.
Only log2 values are serialized; exact total mass is checked before output.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from time import perf_counter

from flint import fmpq, fmpz_mpoly_ctx, fmpz_poly


GROUP_SIZE = 8


def output_group_weight(weight: int) -> int:
    return weight if weight % 2 == 0 else GROUP_SIZE - weight


def log2_positive_integer(value: int) -> float:
    if value <= 0:
        raise ValueError("log2 input must be positive")
    high_bit = value.bit_length() - 1
    shift = max(0, high_bit - 52)
    leading = value >> shift
    return math.log2(leading) + shift


def log2_rational(value: fmpq) -> float | None:
    if not value:
        return None
    numerator = int(value.numerator)
    denominator = int(value.denominator)
    return log2_positive_integer(numerator) - log2_positive_integer(denominator)


def random_linear_expected_spectrum(block_size: int, dimension: int) -> list[float | None]:
    # Expected homogeneous spectrum of a uniformly random dimension-k subspace.
    nonzero_probability_log2 = (
        log2_positive_integer((1 << dimension) - 1)
        - log2_positive_integer((1 << block_size) - 1)
    )
    result: list[float | None] = [None] * (block_size + 1)
    result[0] = 0.0
    for weight in range(1, block_size + 1):
        result[weight] = (
            log2_positive_integer(math.comb(block_size, weight))
            + nonzero_probability_log2
        )
    return result


def first_at_least(log_spectrum: list[float | None], threshold: float) -> int | None:
    for weight in range(1, len(log_spectrum)):
        value = log_spectrum[weight]
        if value is not None and value >= threshold:
            return weight
    return None


def summarize(log_spectrum: list[float | None]) -> dict[str, int | None]:
    return {
        "first_weight_expected_at_least_2^-40": first_at_least(log_spectrum, -40.0),
        "first_weight_expected_at_least_2^-20": first_at_least(log_spectrum, -20.0),
        "first_weight_expected_at_least_1": first_at_least(log_spectrum, 0.0),
    }


def summarize_exact(spectrum: list[fmpq]) -> dict[str, object]:
    cumulative = fmpq(0)
    cumulative_first: dict[str, int | None] = {
        "2^-40": None,
        "2^-20": None,
        "1": None,
    }
    thresholds = {
        "2^-40": fmpq(1, 1 << 40),
        "2^-20": fmpq(1, 1 << 20),
        "1": fmpq(1),
    }
    selected_cutoffs: dict[str, float | None] = {}
    requested_cutoffs = {31, 37, 63, 95, 127}
    for weight in range(1, len(spectrum)):
        cumulative += spectrum[weight]
        for label, threshold in thresholds.items():
            if cumulative_first[label] is None and cumulative >= threshold:
                cumulative_first[label] = weight
        if weight in requested_cutoffs:
            selected_cutoffs[str(weight)] = log2_rational(cumulative)
    return {
        "first_weight_cumulative_expected_at_least": cumulative_first,
        "cumulative_log2_expected_nonzero_codewords_through_weight": selected_cutoffs,
    }


def safe_comb(top: int, bottom: int) -> int:
    if top < 0 or bottom < 0 or bottom > top:
        return 0
    return math.comb(top, bottom)


def apply_accumulator(
    spectrum: list[fmpq], denominators: list[fmpq]
) -> list[fmpq]:
    """Apply a uniform interleaver followed by the binary accumulator."""
    block_size = len(spectrum) - 1
    result = [fmpq(0) for _ in range(block_size + 1)]
    result[0] = spectrum[0]
    for input_weight in range(1, block_size + 1):
        multiplicity = spectrum[input_weight]
        if not multiplicity:
            continue
        normalized = multiplicity / denominators[input_weight]
        low = (input_weight + 1) // 2
        for output_weight in range(low, block_size + 1):
            count = safe_comb(
                block_size - output_weight, input_weight // 2
            ) * safe_comb(
                output_weight - 1, (input_weight + 1) // 2 - 1
            )
            if count:
                result[output_weight] += normalized * count
    return result


def analyze(
    block_size: int,
    layers: int,
    accumulators: int,
    initial_transform: str,
    constituent_size: int,
    expander_degree: int,
) -> dict[str, object]:
    if block_size <= 0 or block_size % GROUP_SIZE:
        raise ValueError("block size must be a positive multiple of eight")
    if block_size % 2:
        raise ValueError("block size must be even")
    if layers < 0:
        raise ValueError("layers must be nonnegative")
    if accumulators < 0:
        raise ValueError("accumulators must be nonnegative")

    dimension = block_size // 2
    context = fmpz_mpoly_ctx.get(("x", "y"), "lex")
    local = context.from_dict(
        {
            (weight, output_group_weight(weight)): math.comb(GROUP_SIZE, weight)
            for weight in range(GROUP_SIZE + 1)
        }
    )

    begin = perf_counter()
    joint_polynomial = local ** (block_size // GROUP_SIZE)
    terms = list(joint_polynomial.terms())
    polynomial_seconds = perf_counter() - begin

    rows: list[list[tuple[int, object]]] = [[] for _ in range(block_size + 1)]
    for (input_weight, output_weight), coefficient in terms:
        rows[input_weight].append((output_weight, coefficient))

    denominators = [fmpq(math.comb(block_size, weight)) for weight in range(block_size + 1)]
    spectrum = [fmpq(0) for _ in range(block_size + 1)]
    if initial_transform == "zeropad":
        for weight in range(dimension + 1):
            spectrum[weight] = fmpq(math.comb(dimension, weight))
    elif initial_transform == "repeat2":
        for input_weight in range(dimension + 1):
            spectrum[2 * input_weight] = fmpq(math.comb(dimension, input_weight))
    elif initial_transform == "systematic-block-random":
        if constituent_size <= 0 or dimension % constituent_size:
            raise ValueError("constituent size must divide the message dimension")
        local_coefficients = [
            (
                (1 << constituent_size)
                if weight == 0
                else math.comb(2 * constituent_size, weight)
                - (
                    math.comb(constituent_size, weight)
                    if weight <= constituent_size
                    else 0
                )
            )
            for weight in range(2 * constituent_size + 1)
        ]
        numerator = fmpz_poly(local_coefficients) ** (
            dimension // constituent_size
        )
        scale = fmpq(1, 1 << dimension)
        for weight in range(block_size + 1):
            spectrum[weight] = scale * numerator[weight]
    elif initial_transform == "systematic-expander":
        if expander_degree <= 0 or expander_degree > dimension:
            raise ValueError("expander degree must be in [1, message dimension]")
        denominator = math.comb(dimension, expander_degree)
        common_denominator = denominator**dimension
        numerators = [0 for _ in range(block_size + 1)]
        for input_weight in range(dimension + 1):
            odd = sum(
                safe_comb(input_weight, intersection)
                * safe_comb(
                    dimension - input_weight,
                    expander_degree - intersection,
                )
                for intersection in range(1, expander_degree + 1, 2)
            )
            even = denominator - odd
            input_count = math.comb(dimension, input_weight)
            odd_powers = [1]
            even_powers = [1]
            for exponent in range(1, dimension + 1):
                odd_powers.append(odd_powers[-1] * odd)
                even_powers.append(even_powers[-1] * even)
            for parity_weight in range(dimension + 1):
                coefficient = (
                    math.comb(dimension, parity_weight)
                    * odd_powers[parity_weight]
                    * even_powers[dimension - parity_weight]
                )
                if coefficient:
                    numerators[input_weight + parity_weight] += (
                        input_count * coefficient
                    )
        spectrum = [
            fmpq(numerator, common_denominator) for numerator in numerators
        ]
    else:
        raise ValueError(f"unsupported initial transform: {initial_transform}")

    layer_seconds: list[float] = []
    layer_summaries: list[dict[str, object]] = []
    for layer in range(1, layers + 1):
        layer_begin = perf_counter()
        next_spectrum = [fmpq(0) for _ in range(block_size + 1)]
        for input_weight, multiplicity in enumerate(spectrum):
            if not multiplicity:
                continue
            normalized = multiplicity / denominators[input_weight]
            for output_weight, count in rows[input_weight]:
                next_spectrum[output_weight] += normalized * count
        spectrum = next_spectrum
        elapsed = perf_counter() - layer_begin
        layer_seconds.append(elapsed)
        layer_log = [log2_rational(value) for value in spectrum]
        layer_summaries.append(
            {
                "stage": f"paritymix_layer_{layer}",
                **summarize(layer_log),
                **summarize_exact(spectrum),
            }
        )

    accumulator_seconds: list[float] = []
    for accumulator in range(1, accumulators + 1):
        stage_begin = perf_counter()
        spectrum = apply_accumulator(spectrum, denominators)
        elapsed = perf_counter() - stage_begin
        accumulator_seconds.append(elapsed)
        stage_log = [log2_rational(value) for value in spectrum]
        layer_summaries.append(
            {
                "stage": f"accumulator_{accumulator}",
                **summarize(stage_log),
                **summarize_exact(spectrum),
            }
        )

    expected_mass = fmpq(1 << dimension)
    if sum(spectrum) != expected_mass:
        raise AssertionError("expected spectrum mass changed")
    if spectrum[0] != 1:
        raise AssertionError("zero codeword multiplicity changed")

    log_spectrum = [log2_rational(value) for value in spectrum]
    random_log_spectrum = random_linear_expected_spectrum(block_size, dimension)
    serialized_spectrum = [
        {
            "weight": weight,
            "log2_expected_multiplicity": log_spectrum[weight],
            "random_linear_log2_expected_multiplicity": random_log_spectrum[weight],
            "excess_over_random_bits": (
                None
                if log_spectrum[weight] is None or random_log_spectrum[weight] is None
                else log_spectrum[weight] - random_log_spectrum[weight]
            ),
        }
        for weight in range(block_size + 1)
    ]

    return {
        "schema": "riffle-paritymix8-outer-expected-spectrum-v1",
        "ensemble": {
            "block_size": block_size,
            "dimension": dimension,
            "rate": 0.5,
            "initial_transform": initial_transform,
            "constituent_size": (
                constituent_size
                if initial_transform == "systematic-block-random"
                else None
            ),
            "expander_degree": (
                expander_degree
                if initial_transform == "systematic-expander"
                else None
            ),
            "group_size": GROUP_SIZE,
            "layers": layers,
            "accumulators": accumulators,
            "local_weight_map": {
                str(weight): output_group_weight(weight)
                for weight in range(GROUP_SIZE + 1)
            },
            "layer_randomness": "independent uniform coordinate permutation before each fixed grouping",
        },
        "method": {
            "local_polynomial": "sum_r binom(8,r) x^r y^(r if r even else 8-r)",
            "joint_transition_terms": len(terms),
            "arithmetic": "exact integer polynomial and exact rational Markov composition; binary64 only for serialized logarithms",
            "mass_check": "PASS",
        },
        "timing_seconds": {
            "joint_polynomial": polynomial_seconds,
            "layers": layer_seconds,
            "accumulators": accumulator_seconds,
        },
        "summary": {**summarize(log_spectrum), **summarize_exact(spectrum)},
        "random_linear_summary": summarize(random_log_spectrum),
        "layer_summaries": layer_summaries,
        "spectrum": serialized_spectrum,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--block-size", type=int, default=1024)
    parser.add_argument("--layers", type=int, default=4)
    parser.add_argument("--accumulators", type=int, default=0)
    parser.add_argument(
        "--initial-transform",
        choices=(
            "zeropad",
            "repeat2",
            "systematic-block-random",
            "systematic-expander",
        ),
        default="zeropad",
    )
    parser.add_argument("--constituent-size", type=int, default=16)
    parser.add_argument("--expander-degree", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = analyze(
        args.block_size,
        args.layers,
        args.accumulators,
        args.initial_transform,
        args.constituent_size,
        args.expander_degree,
    )
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
