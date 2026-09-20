#!/usr/bin/env python3
"""Stream selected primal Schur-block diagonals over message-pair types.

This diagnostic evaluates selected entries of

    2^k D^(j) = 2^k sum_(x,y) det(P_xy)^j Sym^(n-2j)(P_xy)

at rate one half.  Sector zero also subtracts the radial mean outer product.
The program streams multinomial message-pair types and never materializes a
dense block or the 2^(2k) ordered message pairs.

The formulas are exact identities.  Binary64 vector operations and a
long-double reduction make the output diagnostic rather than certified.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from probe_walsh_conjugated_level_bound import (
    covariance_blocks,
    walsh_blocks,
)
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "primal_schur_diagonal_target_probe.json"


def parse_entry(text: str) -> tuple[int, int]:
    try:
        sector_text, level_text = text.split(":", maxsplit=1)
        sector, level = int(sector_text), int(level_text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("entry must have form SECTOR:LEVEL") from error
    if sector < 0 or level < sector:
        raise argparse.ArgumentTypeError("entry must satisfy 0 <= sector <= level")
    return sector, level


def symmetric_diagonal(
    probability: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    degree: int,
    weight: int,
) -> np.ndarray:
    a, b, c, d = probability
    bc = b * c
    result = np.zeros_like(a)
    low = max(0, 2 * weight - degree)
    for intersection in range(low, weight + 1):
        cross = weight - intersection
        coefficient = math.comb(weight, intersection) * math.comb(
            degree - weight, cross
        )
        result += coefficient * (
            np.power(a, degree - 2 * weight + intersection)
            * np.power(bc, cross)
            * np.power(d, intersection)
        )
    return result


def centered_symmetric_diagonal(
    probability: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    determinant: np.ndarray,
    degree: int,
    weight: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate s(P)-s(P_x tensor P_y) with local cancellation."""
    a, b, c, d = probability
    first_one = c + d
    second_one = b + d
    a0 = (1 - first_one) * (1 - second_one)
    b0 = (1 - first_one) * second_one
    c0 = first_one * (1 - second_one)
    d0 = first_one * second_one
    bc = b * c
    bc0 = b0 * c0
    result = np.zeros_like(a)
    absolute_result = np.zeros_like(a)
    telescoping_result = np.zeros_like(a)
    rational_log_result = np.zeros_like(a)
    low = max(0, 2 * weight - degree)
    for intersection in range(low, weight + 1):
        cross = weight - intersection
        zero_zero = degree - 2 * weight + intersection
        coefficient = math.comb(weight, intersection) * math.comb(
            degree - weight, cross
        )
        baseline = coefficient * (
            np.power(a0, zero_zero)
            * np.power(bc0, cross)
            * np.power(d0, intersection)
        )
        regular = np.ones_like(a, dtype=bool)
        if zero_zero:
            regular &= a0 > 0
        if cross:
            regular &= bc0 > 0
        if intersection:
            regular &= d0 > 0
        current = coefficient * (
            np.power(a, zero_zero)
            * np.power(bc, cross)
            * np.power(d, intersection)
        )
        difference = current - baseline
        stable_indices = np.flatnonzero(regular & (baseline > 0))
        if len(stable_indices):
            exponent = np.zeros(len(stable_indices), dtype=a.dtype)
            with np.errstate(divide="ignore", invalid="raise"):
                if zero_zero:
                    relative_a = np.maximum(
                        determinant[stable_indices] / a0[stable_indices], -1
                    )
                    exponent += zero_zero * np.log1p(relative_a)
                if cross:
                    relative_bc = (
                        -determinant[stable_indices]
                        * (b0[stable_indices] + c0[stable_indices])
                        + determinant[stable_indices] * determinant[stable_indices]
                    ) / bc0[stable_indices]
                    exponent += cross * np.log1p(np.maximum(relative_bc, -1))
                if intersection:
                    relative_d = np.maximum(
                        determinant[stable_indices] / d0[stable_indices], -1
                    )
                    exponent += intersection * np.log1p(relative_d)
            nearby = np.isfinite(exponent) & (np.abs(exponent) <= 0.5)
            if np.any(nearby):
                indices = stable_indices[nearby]
                difference[indices] = baseline[indices] * np.expm1(exponent[nearby])
        result += difference
        absolute_result += np.abs(difference)
        maximum_a = np.maximum(a, a0)
        maximum_bc = np.maximum(bc, bc0)
        maximum_d = np.maximum(d, d0)
        absolute_determinant = np.abs(determinant)
        absolute_bc_difference = absolute_determinant * np.abs(
            determinant - (b0 + c0)
        )
        if zero_zero:
            telescoping_result += coefficient * zero_zero * absolute_determinant * (
                np.power(maximum_a, zero_zero - 1)
                * np.power(maximum_bc, cross)
                * np.power(maximum_d, intersection)
            )
        if cross:
            telescoping_result += coefficient * cross * absolute_bc_difference * (
                np.power(maximum_a, zero_zero)
                * np.power(maximum_bc, cross - 1)
                * np.power(maximum_d, intersection)
            )
        if intersection:
            telescoping_result += coefficient * intersection * absolute_determinant * (
                np.power(maximum_a, zero_zero)
                * np.power(maximum_bc, cross)
                * np.power(maximum_d, intersection - 1)
            )
        rational_log = np.zeros_like(a)
        rational_regular = baseline > 0
        if zero_zero:
            denominator_a = np.minimum(a, a0)
            rational_regular &= denominator_a > 0
        if cross:
            denominator_bc = np.minimum(bc, bc0)
            rational_regular &= denominator_bc > 0
        if intersection:
            denominator_d = np.minimum(d, d0)
            rational_regular &= denominator_d > 0
        if np.any(rational_regular):
            if zero_zero:
                rational_log[rational_regular] += (
                    zero_zero
                    * absolute_determinant[rational_regular]
                    / denominator_a[rational_regular]
                )
            if cross:
                rational_log[rational_regular] += (
                    cross
                    * absolute_bc_difference[rational_regular]
                    / denominator_bc[rational_regular]
                )
            if intersection:
                rational_log[rational_regular] += (
                    intersection
                    * absolute_determinant[rational_regular]
                    / denominator_d[rational_regular]
                )
        small = rational_regular & (rational_log < 1)
        rational_term = current + baseline
        rational_term[small] = (
            baseline[small]
            * rational_log[small]
            / (1 - rational_log[small])
        )
        rational_log_result += rational_term
    return result, absolute_result, telescoping_result, rational_log_result


def mean_word_probability(
    message_bits: int,
    output_bits: int,
    biases: np.ndarray,
    level: int,
) -> np.longdouble:
    terms = []
    for weight in range(message_bits + 1):
        q = np.longdouble((1 - biases[weight]) / 2)
        terms.append(
            np.longdouble(math.comb(message_bits, weight))
            * q**level
            * (1 - q) ** (output_bits - level)
        )
    return sum(terms, np.longdouble(0))


def streamed_entries(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    entries: list[tuple[int, int]],
    extended: bool,
) -> tuple[
    dict[tuple[int, int], np.longdouble],
    dict[tuple[int, int], np.longdouble],
    dict[tuple[int, int], np.longdouble],
    dict[tuple[int, int], np.longdouble],
    dict[tuple[int, int], np.longdouble],
    dict[tuple[int, int], np.longdouble],
    dict[tuple[int, int], np.longdouble],
    dict[tuple[int, int], np.longdouble],
    int,
]:
    dtype = np.longdouble if extended else np.float64
    denominator = math.comb(message_bits, right_degree)
    biases = np.asarray(
        [
            dtype(krawtchouk(message_bits, right_degree, weight))
            / dtype(denominator)
            for weight in range(message_bits + 1)
        ],
        dtype=dtype,
    )
    log_factorial = np.asarray(
        [math.lgamma(value + 1) for value in range(message_bits + 1)]
    )
    scaled = {entry: np.longdouble(0) for entry in entries}
    positive_parts = {entry: np.longdouble(0) for entry in entries}
    centered_positive_parts = {entry: np.longdouble(0) for entry in entries}
    centered_sums = {entry: np.longdouble(0) for entry in entries}
    sector_zero_tv_bounds = {entry: np.longdouble(0) for entry in entries}
    centered_termwise_absolute_parts = {
        entry: np.longdouble(0) for entry in entries
    }
    centered_telescoping_bounds = {
        entry: np.longdouble(0) for entry in entries
    }
    centered_rational_log_bounds = {
        entry: np.longdouble(0) for entry in entries
    }
    type_count = 0
    scale_log = message_bits * math.log(2)
    scale_exact = dtype(2) ** message_bits

    # Fix n11.  The remaining three counts form one vectorized chunk.
    for n11 in range(message_bits + 1):
        n00_rows: list[int] = []
        n01_rows: list[int] = []
        n10_rows: list[int] = []
        for n01 in range(message_bits - n11 + 1):
            for n10 in range(message_bits - n11 - n01 + 1):
                n00_rows.append(message_bits - n11 - n01 - n10)
                n01_rows.append(n01)
                n10_rows.append(n10)
        n00 = np.asarray(n00_rows, dtype=np.int16)
        n01 = np.asarray(n01_rows, dtype=np.int16)
        n10 = np.asarray(n10_rows, dtype=np.int16)
        first_weight = n10 + n11
        second_weight = n01 + n11
        difference_weight = n01 + n10
        first_bias = biases[first_weight]
        second_bias = biases[second_weight]
        difference_bias = biases[difference_weight]
        probability = (
            (1 + first_bias + second_bias + difference_bias) / 4,
            (1 + first_bias - second_bias - difference_bias) / 4,
            (1 - first_bias + second_bias - difference_bias) / 4,
            (1 - first_bias - second_bias + difference_bias) / 4,
        )
        determinant = (
            difference_bias - first_bias * second_bias
        ) / 4
        if extended:
            exact_multiplicities = [
                math.comb(message_bits, int(d11))
                * math.comb(message_bits - int(d11), int(d01))
                * math.comb(
                    message_bits - int(d11) - int(d01), int(d10)
                )
                for d01, d10, d11 in zip(
                    n01, n10, np.full_like(n01, n11), strict=True
                )
            ]
            type_scale = np.asarray(exact_multiplicities, dtype=dtype) * scale_exact
        else:
            log_multiplicity = (
                log_factorial[message_bits]
                - log_factorial[n00]
                - log_factorial[n01]
                - log_factorial[n10]
                - log_factorial[n11]
            )
            type_scale = np.exp(log_multiplicity + scale_log)
        for sector, level in entries:
            degree = output_bits - 2 * sector
            weight = level - sector
            diagonal = symmetric_diagonal(probability, degree, weight)
            contribution = type_scale * np.power(determinant, sector) * diagonal
            scaled[(sector, level)] += np.sum(
                contribution, dtype=np.longdouble
            )
            if sector % 2:
                positive_parts[(sector, level)] += np.sum(
                    np.maximum(contribution, 0), dtype=np.longdouble
                )
            elif sector == 0:
                (
                    centered,
                    centered_absolute,
                    centered_telescoping,
                    centered_rational_log,
                ) = centered_symmetric_diagonal(
                    probability, determinant, degree, weight
                )
                centered *= type_scale
                centered_absolute *= type_scale
                centered_telescoping *= type_scale
                centered_rational_log *= type_scale
                centered_sums[(sector, level)] += np.sum(
                    centered, dtype=np.longdouble
                )
                centered_positive_parts[(sector, level)] += np.sum(
                    np.maximum(centered, 0), dtype=np.longdouble
                )
                centered_termwise_absolute_parts[(sector, level)] += np.sum(
                    centered_absolute, dtype=np.longdouble
                )
                centered_telescoping_bounds[(sector, level)] += np.sum(
                    centered_telescoping, dtype=np.longdouble
                )
                centered_rational_log_bounds[(sector, level)] += np.sum(
                    centered_rational_log, dtype=np.longdouble
                )
                first_one = probability[2] + probability[3]
                first_pattern = (
                    np.power(first_one, weight)
                    * np.power(1 - first_one, degree - weight)
                )
                interior = (first_one > 0) & (first_one < 1)
                sensitivity = np.zeros_like(first_one)
                sensitivity[interior] = (
                    weight / first_one[interior]
                    + (degree - weight) / (1 - first_one[interior])
                )
                conditional_tv = np.minimum(
                    1, np.abs(determinant) * sensitivity
                )
                tv_contribution = type_scale * first_pattern * conditional_tv
                sector_zero_tv_bounds[(sector, level)] += np.sum(
                    tv_contribution, dtype=np.longdouble
                )
        type_count += len(n00)

    uncentered = scaled.copy()
    for sector, level in entries:
        if sector != 0:
            continue
        mean = mean_word_probability(
            message_bits, output_bits, biases, level
        )
        subtraction = (
            np.longdouble(2) ** message_bits
            * np.longdouble(math.comb(output_bits, level))
            * mean
            * mean
        )
        scaled[(sector, level)] -= subtraction
    return (
        scaled,
        uncentered,
        positive_parts,
        centered_sums,
        centered_positive_parts,
        centered_termwise_absolute_parts,
        centered_telescoping_bounds,
        centered_rational_log_bounds,
        sector_zero_tv_bounds,
        type_count,
    )


def full_block_comparison(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    entries: list[tuple[int, int]],
) -> dict[tuple[int, int], float]:
    covariance, _cached = covariance_blocks(
        message_bits, output_bits, right_degree
    )
    walsh, _residual = walsh_blocks(output_bits)
    normalization = math.ldexp(1.0, -output_bits)
    result = {}
    for sector, level in entries:
        block = normalization * walsh[sector] @ covariance[sector] @ walsh[sector].T
        result[(sector, level)] = math.ldexp(
            float(block[level - sector, level - sector]), message_bits
        )
    return result


def evaluate(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    entries: list[tuple[int, int]],
    validate_full: bool,
    extended: bool,
) -> dict[str, object]:
    (
        scaled,
        uncentered,
        positive_parts,
        centered_sums,
        centered_positive_parts,
        centered_termwise_absolute_parts,
        centered_telescoping_bounds,
        centered_rational_log_bounds,
        sector_zero_tv_bounds,
        type_count,
    ) = streamed_entries(
        message_bits, output_bits, right_degree, entries, extended
    )
    reference = (
        full_block_comparison(
            message_bits, output_bits, right_degree, entries
        )
        if validate_full
        else {}
    )
    rows = []
    maximum_relative_residual = 0.0
    for entry in entries:
        value = float(scaled[entry])
        expected = reference.get(entry)
        relative = None
        if expected is not None:
            relative = abs(value - expected) / max(1.0, abs(expected))
            maximum_relative_residual = max(maximum_relative_residual, relative)
        rows.append(
            {
                "sector": entry[0],
                "level": entry[1],
                "scaled_diagonal": value,
                "scaled_uncentered_pair_sum": (
                    float(uncentered[entry]) if entry[0] == 0 else None
                ),
                "scaled_positive_part": (
                    float(positive_parts[entry]) if entry[0] % 2 else None
                ),
                "scaled_centered_positive_part": (
                    float(centered_positive_parts[entry])
                    if entry[0] == 0
                    else None
                ),
                "scaled_locally_centered_sum": (
                    float(centered_sums[entry]) if entry[0] == 0 else None
                ),
                "scaled_centered_termwise_absolute_upper": (
                    float(centered_termwise_absolute_parts[entry])
                    if entry[0] == 0
                    else None
                ),
                "scaled_centered_telescoping_upper": (
                    float(centered_telescoping_bounds[entry])
                    if entry[0] == 0
                    else None
                ),
                "scaled_centered_rational_log_upper": (
                    float(centered_rational_log_bounds[entry])
                    if entry[0] == 0
                    else None
                ),
                "scaled_sector_zero_conditional_tv_upper": (
                    float(sector_zero_tv_bounds[entry])
                    if entry[0] == 0
                    else None
                ),
                "positive_part_over_value": (
                    float(positive_parts[entry] / scaled[entry])
                    if entry[0] % 2 and scaled[entry] > 0
                    else None
                ),
                "unscaled_log2": (
                    None if value <= 0 else math.log2(value) - message_bits
                ),
                "full_block_scaled_diagonal": expected,
                "relative_residual": relative,
            }
        )
    return {
        "message_bits": message_bits,
        "output_bits": output_bits,
        "right_degree": right_degree,
        "enumerated_pair_types": type_count,
        "arithmetic": (
            "numpy.longdouble with exact integer multinomial conversion"
            if extended
            else "binary64 with long-double reduction"
        ),
        "entries": rows,
        "maximum_relative_validation_residual": (
            maximum_relative_residual if validate_full else None
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=32)
    parser.add_argument("--output-bits", type=int, default=64)
    parser.add_argument("--right-degree", type=int, default=7)
    parser.add_argument(
        "--entry", type=parse_entry, action="append", default=[]
    )
    parser.add_argument("--all-through", type=int)
    parser.add_argument("--maximum-sector", type=int, default=2)
    parser.add_argument("--validate-full", action="store_true")
    parser.add_argument("--extended", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits != 2 * args.message_bits:
        parser.error("the present probe expects output_bits = 2 * message_bits")
    if args.all_through is not None:
        if args.entry:
            parser.error("use either --entry or --all-through")
        entries = [
            (sector, level)
            for sector in range(args.maximum_sector + 1)
            for level in range(max(1, sector), args.all_through + 1)
        ]
    else:
        entries = args.entry or [(1, 1), (1, 4), (2, 8)]
    if any(level > args.output_bits - sector for sector, level in entries):
        parser.error("entry is outside its sector block")
    if args.validate_full and args.output_bits > 128:
        parser.error("full validation is limited to output_bits <= 128")
    result = evaluate(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        entries,
        args.validate_full,
        args.extended,
    )
    payload = {
        "schema": "pure-ea-primal-schur-diagonal-target-probe-v1",
        "status": (
            "EXTENDED_FLOAT_WITH_EXACT_MULTINOMIAL_WEIGHTS_DIAGNOSTIC"
            if args.extended
            else "BINARY64_WITH_LONGDOUBLE_REDUCTION_DIAGNOSTIC"
        ),
        "result": result,
        "scope": [
            "The streamed pair-type and Schur-power formulas are exact identities.",
            (
                "Vector operations and reductions use nondirected numpy.longdouble. Multinomial weights are converted from exact integers."
                if args.extended
                else "Vector operations use nondirected binary64; reductions use numpy.longdouble."
            ),
            "The output is not a rigorous interval certificate.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
