#!/usr/bin/env python3
"""Outward upper bounds for the three controlling primal Schur diagonals.

The calculation uses only positive one-sided sums.  Sector zero uses the
rational-log majorant for each locally centered monomial.  Sector one drops
negative determinant contributions.  Sector two is nonnegative.  Binary64
interval endpoints are expanded with ``nextafter`` after every basic
operation; integer powers use repeated multiplication.

The program is intentionally a selected-entry verifier.  Use ``--entry`` to
certify isolated gates, or ``--all-through`` for a complete prefix of levels.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import time

import numpy as np

from analyze_dual_walk_and_accumulator_energy import accumulator_joint_counts
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "primal_schur_diagonal_outward.json"
POSITIVE_INFINITY = np.float64(math.inf)
NEGATIVE_INFINITY = np.float64(-math.inf)
UNIT_ROUNDOFF = 2.0**-53
MIN_SUBNORMAL = np.nextafter(np.float64(0), POSITIVE_INFINITY)


def parse_entry(text: str) -> tuple[int, int]:
    try:
        sector_text, level_text = text.split(":", maxsplit=1)
        sector, level = int(sector_text), int(level_text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("entry must have form SECTOR:LEVEL") from error
    if sector not in (0, 1, 2) or level < sector:
        raise argparse.ArgumentTypeError("entry must satisfy sector in {0,1,2} and level >= sector")
    return sector, level


def down(value: np.ndarray | np.float64) -> np.ndarray | np.float64:
    return np.nextafter(value, NEGATIVE_INFINITY)


def up(value: np.ndarray | np.float64) -> np.ndarray | np.float64:
    return np.nextafter(value, POSITIVE_INFINITY)


def add_interval(
    left: tuple[np.ndarray, np.ndarray], right: tuple[np.ndarray, np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    return down(left[0] + right[0]), up(left[1] + right[1])


def subtract_interval(
    left: tuple[np.ndarray, np.ndarray], right: tuple[np.ndarray, np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    return down(left[0] - right[1]), up(left[1] - right[0])


def multiply_interval(
    left: tuple[np.ndarray, np.ndarray], right: tuple[np.ndarray, np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return down(np.minimum.reduce(products)), up(np.maximum.reduce(products))


def positive_multiply_upper(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return up(left * right)


def positive_divide_upper(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return up(numerator / denominator)


def positive_power_upper(base: np.ndarray, exponent: int) -> np.ndarray:
    result = np.ones_like(base)
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = positive_multiply_upper(result, factor)
        power >>= 1
        if power:
            factor = positive_multiply_upper(factor, factor)
    return result


def power_prefix_upper(base: np.ndarray, maximum: int) -> list[np.ndarray]:
    """Return outward powers 0 through ``maximum`` by sequential products."""
    result = [np.ones_like(base)]
    for _exponent in range(maximum):
        result.append(positive_multiply_upper(result[-1], base))
    return result


def positive_sum_upper(values: np.ndarray) -> np.float64:
    """Upper-bound a nonnegative reduction without assuming its summation tree."""
    if np.any(values < 0) or np.any(~np.isfinite(values)):
        raise ArithmeticError("positive reduction received an invalid term")
    count = int(values.size)
    if count == 0:
        return np.float64(0)
    rounded = np.float64(np.sum(values, dtype=np.float64))
    operations = max(0, count - 1)
    relative_allowance = 2 * operations * UNIT_ROUNDOFF
    # Any binary summation tree uses at most count-1 additions.  The additive
    # term covers gradual-underflow errors in those additions.
    corrected = (
        rounded / (1 - relative_allowance) + operations * MIN_SUBNORMAL
    )
    return np.float64(up(np.float64(corrected)))


def add_upper_scalar(left: np.float64, right: np.float64) -> np.float64:
    return np.float64(up(np.float64(left + right)))


def rational_interval(numerator: int, denominator: int) -> tuple[float, float]:
    exact = Fraction(numerator, denominator)
    nearest = float(exact)
    represented = Fraction.from_float(nearest)
    low = nearest if represented <= exact else math.nextafter(nearest, -math.inf)
    high = nearest if represented >= exact else math.nextafter(nearest, math.inf)
    return low, high


def coefficient_upper(value: int) -> np.float64:
    return np.float64(math.nextafter(float(value), math.inf))


def nonnegative_interval(
    interval: tuple[np.ndarray, np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    low, high = interval
    if np.any(high < 0):
        raise ArithmeticError("probability interval lies below zero")
    return np.maximum(low, 0), np.maximum(high, 0)


def symmetric_diagonal_upper(
    probability: tuple[tuple[np.ndarray, np.ndarray], ...],
    degree: int,
    weight: int,
) -> np.ndarray:
    a, b, c, d = probability
    bc = multiply_interval(b, c)
    result = np.zeros_like(a[1])
    low = max(0, 2 * weight - degree)
    cross_powers = power_prefix_upper(bc[1], weight - low)
    a_power = positive_power_upper(a[1], degree - 2 * weight + low)
    d_power = positive_power_upper(d[1], low)
    for intersection in range(low, weight + 1):
        cross = weight - intersection
        zero_zero = degree - 2 * weight + intersection
        term = np.full_like(result, coefficient_upper(
            math.comb(weight, intersection) * math.comb(degree - weight, cross)
        ))
        term = positive_multiply_upper(term, a_power)
        term = positive_multiply_upper(term, cross_powers[cross])
        term = positive_multiply_upper(term, d_power)
        result = up(result + term)
        if intersection < weight:
            a_power = positive_multiply_upper(a_power, a[1])
            d_power = positive_multiply_upper(d_power, d[1])
    return result


def centered_rational_log_upper(
    probability: tuple[tuple[np.ndarray, np.ndarray], ...],
    determinant: tuple[np.ndarray, np.ndarray],
    degree: int,
    weight: int,
) -> np.ndarray:
    """Bound the absolute locally centered diagonal by positive operations."""
    a, b, c, d = probability
    one = (np.ones_like(a[0]), np.ones_like(a[1]))
    first_one = add_interval(c, d)
    second_one = add_interval(b, d)
    first_zero = subtract_interval(one, first_one)
    second_zero = subtract_interval(one, second_one)
    a0 = nonnegative_interval(multiply_interval(first_zero, second_zero))
    b0 = nonnegative_interval(multiply_interval(first_zero, second_one))
    c0 = nonnegative_interval(multiply_interval(first_one, second_zero))
    d0 = nonnegative_interval(multiply_interval(first_one, second_one))
    bc = nonnegative_interval(multiply_interval(b, c))
    bc0 = nonnegative_interval(multiply_interval(b0, c0))
    absolute_determinant = np.maximum(np.abs(determinant[0]), np.abs(determinant[1]))
    b0_plus_c0 = add_interval(b0, c0)
    det_minus_sum = subtract_interval(determinant, b0_plus_c0)
    absolute_det_minus_sum = np.maximum(
        np.abs(det_minus_sum[0]), np.abs(det_minus_sum[1])
    )
    absolute_bc_difference = positive_multiply_upper(
        absolute_determinant, absolute_det_minus_sum
    )

    result = np.zeros_like(a[0])
    low = max(0, 2 * weight - degree)
    current_cross_powers = power_prefix_upper(bc[1], weight - low)
    baseline_cross_powers = power_prefix_upper(bc0[1], weight - low)
    current_a_power = positive_power_upper(a[1], degree - 2 * weight + low)
    baseline_a_power = positive_power_upper(a0[1], degree - 2 * weight + low)
    current_d_power = positive_power_upper(d[1], low)
    baseline_d_power = positive_power_upper(d0[1], low)
    for intersection in range(low, weight + 1):
        cross = weight - intersection
        zero_zero = degree - 2 * weight + intersection
        coefficient = coefficient_upper(
            math.comb(weight, intersection) * math.comb(degree - weight, cross)
        )
        baseline = np.full_like(result, coefficient)
        baseline = positive_multiply_upper(baseline, baseline_a_power)
        baseline = positive_multiply_upper(baseline, baseline_cross_powers[cross])
        baseline = positive_multiply_upper(baseline, baseline_d_power)
        current = np.full_like(result, coefficient)
        current = positive_multiply_upper(current, current_a_power)
        current = positive_multiply_upper(current, current_cross_powers[cross])
        current = positive_multiply_upper(current, current_d_power)

        regular = np.ones_like(result, dtype=bool)
        logarithm = np.zeros_like(result)
        if zero_zero:
            denominator = np.minimum(a[0], a0[0])
            regular &= denominator > 0
            indices = denominator > 0
            quotient = positive_divide_upper(
                absolute_determinant[indices], denominator[indices]
            )
            logarithm[indices] = positive_multiply_upper(
                np.full_like(quotient, coefficient_upper(zero_zero)), quotient
            )
        if cross:
            denominator = np.minimum(bc[0], bc0[0])
            regular &= denominator > 0
            indices = denominator > 0
            quotient = positive_divide_upper(
                absolute_bc_difference[indices], denominator[indices]
            )
            addition = positive_multiply_upper(
                np.full_like(quotient, coefficient_upper(cross)), quotient
            )
            logarithm[indices] = up(logarithm[indices] + addition)
        if intersection:
            denominator = np.minimum(d[0], d0[0])
            regular &= denominator > 0
            indices = denominator > 0
            quotient = positive_divide_upper(
                absolute_determinant[indices], denominator[indices]
            )
            addition = positive_multiply_upper(
                np.full_like(quotient, coefficient_upper(intersection)), quotient
            )
            logarithm[indices] = up(logarithm[indices] + addition)

        small = regular & (logarithm < 1)
        term = up(current + baseline)
        if np.any(small):
            denominator = down(np.float64(1) - logarithm[small])
            controlled = positive_divide_upper(logarithm[small], denominator)
            term[small] = positive_multiply_upper(baseline[small], controlled)

        result = up(result + term)
        if intersection < weight:
            current_a_power = positive_multiply_upper(current_a_power, a[1])
            baseline_a_power = positive_multiply_upper(baseline_a_power, a0[1])
            current_d_power = positive_multiply_upper(current_d_power, d[1])
            baseline_d_power = positive_multiply_upper(baseline_d_power, d0[1])
    return result


def bias_intervals(message_bits: int, right_degree: int) -> tuple[np.ndarray, np.ndarray]:
    denominator = math.comb(message_bits, right_degree)
    rows = [
        rational_interval(
            krawtchouk(message_bits, right_degree, weight), denominator
        )
        for weight in range(message_bits + 1)
    ]
    return (
        np.asarray([row[0] for row in rows], dtype=np.float64),
        np.asarray([row[1] for row in rows], dtype=np.float64),
    )


def streamed_upper_entries(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    entries: list[tuple[int, int]],
    progress_every: int,
) -> tuple[dict[tuple[int, int], np.float64], int, float]:
    bias_low, bias_high = bias_intervals(message_bits, right_degree)
    totals = {entry: np.float64(0) for entry in entries}
    type_count = 0
    started = time.perf_counter()
    one = None
    for n11 in range(message_bits + 1):
        n00_rows: list[int] = []
        n01_rows: list[int] = []
        n10_rows: list[int] = []
        for n01 in range(message_bits - n11 + 1):
            for n10 in range(message_bits - n11 - n01 + 1):
                n00_rows.append(message_bits - n11 - n01 - n10)
                n01_rows.append(n01)
                n10_rows.append(n10)
        n01_array = np.asarray(n01_rows, dtype=np.int16)
        n10_array = np.asarray(n10_rows, dtype=np.int16)
        first_weight = n10_array + n11
        second_weight = n01_array + n11
        difference_weight = n01_array + n10_array
        first_bias = (bias_low[first_weight], bias_high[first_weight])
        second_bias = (bias_low[second_weight], bias_high[second_weight])
        difference_bias = (bias_low[difference_weight], bias_high[difference_weight])
        if one is None or len(one[0]) != len(n01_array):
            one = (np.ones(len(n01_array)), np.ones(len(n01_array)))
        p00 = add_interval(add_interval(add_interval(one, first_bias), second_bias), difference_bias)
        p01 = subtract_interval(add_interval(one, first_bias), add_interval(second_bias, difference_bias))
        p10 = subtract_interval(add_interval(one, second_bias), add_interval(first_bias, difference_bias))
        p11 = add_interval(subtract_interval(one, add_interval(first_bias, second_bias)), difference_bias)
        probability = tuple(
            nonnegative_interval((np.ldexp(row[0], -2), np.ldexp(row[1], -2)))
            for row in (p00, p01, p10, p11)
        )
        determinant = (
            np.ldexp(
                subtract_interval(
                    difference_bias, multiply_interval(first_bias, second_bias)
                )[0],
                -2,
            ),
            np.ldexp(
                subtract_interval(
                    difference_bias, multiply_interval(first_bias, second_bias)
                )[1],
                -2,
            ),
        )
        exact_multiplicities = [
            math.comb(message_bits, n11)
            * math.comb(message_bits - n11, int(n01))
            * math.comb(message_bits - n11 - int(n01), int(n10))
            for n01, n10 in zip(n01_array, n10_array, strict=True)
        ]
        multiplicity_upper = up(np.asarray(exact_multiplicities, dtype=np.float64))
        type_scale_upper = np.ldexp(multiplicity_upper, message_bits)
        absolute_determinant = np.maximum(
            np.abs(determinant[0]), np.abs(determinant[1])
        )
        for sector, level in entries:
            degree = output_bits - 2 * sector
            weight = level - sector
            if sector == 0:
                diagonal = centered_rational_log_upper(
                    probability, determinant, degree, weight
                )
                contribution = positive_multiply_upper(type_scale_upper, diagonal)
            else:
                diagonal = symmetric_diagonal_upper(probability, degree, weight)
                if sector == 1:
                    determinant_factor = np.maximum(determinant[1], 0)
                else:
                    determinant_factor = positive_multiply_upper(
                        absolute_determinant, absolute_determinant
                    )
                contribution = positive_multiply_upper(
                    type_scale_upper,
                    positive_multiply_upper(determinant_factor, diagonal),
                )
            totals[(sector, level)] = add_upper_scalar(
                totals[(sector, level)], positive_sum_upper(contribution)
            )
        type_count += len(n01_array)
        if progress_every and (
            n11 % progress_every == 0 or n11 == message_bits
        ):
            elapsed = time.perf_counter() - started
            print(
                f"progress,n11,{n11},{message_bits},pair_types,{type_count},elapsed_seconds,{elapsed:.3f}",
                flush=True,
            )
    return totals, type_count, time.perf_counter() - started


def mean_shell_lower(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
) -> np.float64:
    bias_low, bias_high = bias_intervals(message_bits, right_degree)
    run_counts = accumulator_joint_counts(output_bits)[shell_weight]
    total = np.float64(0)
    for message_weight in range(1, message_bits + 1):
        q_low = max(0.0, math.nextafter((1 - bias_high[message_weight]) / 2, -math.inf))
        q_high = min(1.0, math.nextafter((1 - bias_low[message_weight]) / 2, math.inf))
        one_minus_q_low = max(0.0, math.nextafter(1 - q_high, -math.inf))
        probability = np.float64(0)
        for transition_weight, count in run_counts.items():
            first = positive_power_lower(np.float64(q_low), transition_weight)
            second = positive_power_lower(
                np.float64(one_minus_q_low), output_bits - transition_weight
            )
            coefficient = coefficient_lower(count)
            term = positive_multiply_lower(
                coefficient, positive_multiply_lower(first, second)
            )
            probability = add_lower_scalar(probability, term)
        multiplicity = coefficient_lower(math.comb(message_bits, message_weight))
        total = add_lower_scalar(
            total, positive_multiply_lower(multiplicity, probability)
        )
    return total


def mean_shell_upper(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
) -> np.float64:
    """Return an outward binary64 upper bound on the nonzero-message mean."""
    bias_low, bias_high = bias_intervals(message_bits, right_degree)
    run_counts = accumulator_joint_counts(output_bits)[shell_weight]
    total = np.float64(0)
    for message_weight in range(1, message_bits + 1):
        q_high = min(
            1.0,
            math.nextafter((1 - bias_low[message_weight]) / 2, math.inf),
        )
        one_minus_q_high = min(
            1.0,
            math.nextafter((1 + bias_high[message_weight]) / 2, math.inf),
        )
        probability = np.float64(0)
        for transition_weight, count in run_counts.items():
            first = positive_power_upper(
                np.asarray([q_high], dtype=np.float64), transition_weight
            )[0]
            second = positive_power_upper(
                np.asarray([one_minus_q_high], dtype=np.float64),
                output_bits - transition_weight,
            )[0]
            term = positive_multiply_upper(
                np.asarray([coefficient_upper(count)]),
                positive_multiply_upper(
                    np.asarray([first]), np.asarray([second])
                ),
            )[0]
            probability = add_upper_scalar(probability, term)
        term = positive_multiply_upper(
            np.asarray([coefficient_upper(math.comb(message_bits, message_weight))]),
            np.asarray([probability]),
        )[0]
        total = add_upper_scalar(total, term)
    return total


def positive_multiply_lower(left: np.float64, right: np.float64) -> np.float64:
    return np.float64(max(0.0, down(np.float64(left * right))))


def positive_power_lower(base: np.float64, exponent: int) -> np.float64:
    result = np.float64(1)
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = positive_multiply_lower(result, factor)
        power >>= 1
        if power:
            factor = positive_multiply_lower(factor, factor)
    return result


def coefficient_lower(value: int) -> np.float64:
    nearest = float(value)
    represented = Fraction.from_float(nearest)
    return np.float64(
        nearest if represented <= value else math.nextafter(nearest, -math.inf)
    )


def add_lower_scalar(left: np.float64, right: np.float64) -> np.float64:
    return np.float64(max(0.0, down(np.float64(left + right))))


def shell_ratio_upper(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
    diagonal_upper: dict[tuple[int, int], np.float64],
) -> tuple[np.float64, np.float64]:
    counts = accumulator_joint_counts(output_bits)[shell_weight]
    total = np.float64(0)
    for level, count in sorted(counts.items()):
        candidates = [
            diagonal_upper[(sector, level)]
            for sector in range(min(2, level) + 1)
        ]
        diagonal = max(candidates)
        product = positive_multiply_upper(
            np.asarray([coefficient_upper(count)]),
            np.asarray([diagonal]),
        )[0]
        root = np.float64(up(np.float64(math.sqrt(product))))
        total = add_upper_scalar(total, root)
    quadratic = positive_multiply_upper(
        np.asarray([total]), np.asarray([total])
    )[0]
    variance_upper = np.ldexp(quadratic, -message_bits)
    mean_lower = mean_shell_lower(
        message_bits, output_bits, right_degree, shell_weight
    )
    if mean_lower <= 0:
        raise ArithmeticError("mean lower endpoint is not positive")
    ratio = positive_divide_upper(
        np.asarray([variance_upper]), np.asarray([mean_lower])
    )[0]
    return np.float64(ratio), mean_lower


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--right-degree", type=int, default=33)
    parser.add_argument("--entry", type=parse_entry, action="append", default=[])
    parser.add_argument("--all-through", type=int)
    parser.add_argument("--level-min", type=int, default=0)
    parser.add_argument("--shell-weight", type=int, action="append", default=[])
    parser.add_argument("--progress-every", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    entries = list(args.entry)
    if args.all_through is not None:
        if not 0 <= args.level_min <= args.all_through:
            parser.error("level-min must lie between zero and all-through")
        entries.extend(
            (sector, level)
            for level in range(args.level_min, args.all_through + 1)
            for sector in range(min(2, level) + 1)
        )
    entries = sorted(set(entries))
    if not entries:
        parser.error("supply --entry or --all-through")
    if any(level > args.output_bits // 2 for _sector, level in entries):
        parser.error("the three-sector reduction is used only through n/2")

    values, type_count, elapsed = streamed_upper_entries(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        entries,
        args.progress_every,
    )
    shell_rows = []
    for shell_weight in args.shell_weight:
        ratio, mean_lower = shell_ratio_upper(
            args.message_bits,
            args.output_bits,
            args.right_degree,
            shell_weight,
            values,
        )
        shell_rows.append(
            {
                "shell_weight": shell_weight,
                "mean_lower": float(mean_lower),
                "mean_log2_lower": math.log2(float(mean_lower)),
                "variance_to_mean_upper": float(ratio),
                "variance_to_mean_log2_upper": math.log2(float(ratio)),
                "passes_factor_512": bool(ratio <= 512),
            }
        )
    payload = {
        "schema": "pure-ea-primal-schur-diagonal-outward-v1",
        "status": "OUTWARD_BINARY64_UPPER_BOUND",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
        },
        "enumerated_pair_types": type_count,
        "elapsed_seconds": elapsed,
        "entries": [
            {
                "sector": sector,
                "level": level,
                "scaled_diagonal_upper": float(values[(sector, level)]),
            }
            for sector, level in entries
        ],
        "shells": shell_rows,
        "proof_scope": [
            "Every reported diagonal is an upper bound under IEEE-754 binary64 round-to-nearest arithmetic with gradual underflow.",
            "Sector zero uses a termwise rational-log bound on the absolute locally centered contribution.",
            "Sector one discards negative determinant contributions; sector two is nonnegative.",
            "The exact diagonal insertion lemma reduces all higher sectors to sector two for levels at most output_bits/2.",
            "The shell contraction uses exact accumulator run counts and an outward lower bound on the shell mean.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
