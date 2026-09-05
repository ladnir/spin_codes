#!/usr/bin/env python3
"""Outward upper bounds from four-message complement orbits.

For odd row degree, complementing either input message complements that
message's sparse-map output.  Averaging a centered pair contribution over
the resulting four-message orbit leaves the total sector-zero diagonal
unchanged.  This verifier encloses every orbit with directed binary64
intervals and sums the positive upper endpoints.  It does not assume that
the orbit contributions are nonnegative.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from certify_primal_schur_diagonal_outward import (
    add_interval,
    add_upper_scalar,
    bias_intervals,
    down,
    multiply_interval,
    nonnegative_interval,
    positive_multiply_upper,
    positive_sum_upper,
    rational_interval,
    subtract_interval,
    up,
)
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "sector_zero_complement_orbits_outward.json"
Interval = tuple[np.ndarray, np.ndarray]


def pair_type_orbit(row: tuple[int, int, int, int]) -> set[tuple[int, int, int, int]]:
    """Return the distinct complement-and-swap images of (00,01,10,11)."""
    a, b, c, d = row
    complements = {
        (a, b, c, d),
        (c, d, a, b),
        (b, a, d, c),
        (d, c, b, a),
    }
    return complements | {(x[0], x[2], x[1], x[3]) for x in complements}


def positive_power_interval(base: Interval, exponent: int) -> Interval:
    one = np.ones_like(base[0])
    result = (one, one)
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = multiply_interval(result, factor)
        power >>= 1
        if power:
            factor = multiply_interval(factor, factor)
    return result


def divide_by_positive(numerator: Interval, denominator: Interval) -> Interval:
    """Divide a signed interval by an interval with positive lower bound."""
    valid = denominator[0] > 0
    safe_low = np.where(valid, denominator[0], 1)
    safe_high = np.where(valid, denominator[1], 1)
    low_denominator = np.where(numerator[0] >= 0, safe_high, safe_low)
    high_denominator = np.where(numerator[1] >= 0, safe_low, safe_high)
    low = down(numerator[0] / low_denominator)
    high = up(numerator[1] / high_denominator)
    return np.where(valid, low, 0), np.where(valid, high, 0)


def delta_product(left: Interval, right: Interval) -> Interval:
    """Given intervals for x and y, enclose (1+x)(1+y)-1."""
    return add_interval(add_interval(left, right), multiply_interval(left, right))


def delta_power(delta: Interval, exponent: int) -> Interval:
    """Enclose (1+delta)^exponent-1 without subtracting nearby numbers."""
    zero = np.zeros_like(delta[0])
    result = (zero, zero)
    factor = delta
    power = exponent
    while power:
        if power & 1:
            result = delta_product(result, factor)
        power >>= 1
        if power:
            factor = delta_product(factor, factor)
    return result


def centered_diagonal_interval(
    first_one: Interval,
    second_one: Interval,
    determinant: Interval,
    degree: int,
    weight: int,
) -> Interval:
    """Enclose s(P)-s(P_x tensor P_y) with correlation-scaled error."""
    one = (np.ones_like(first_one[0]), np.ones_like(first_one[1]))
    first_zero = subtract_interval(one, first_one)
    second_zero = subtract_interval(one, second_one)
    a0 = nonnegative_interval(multiply_interval(first_zero, second_zero))
    b0 = nonnegative_interval(multiply_interval(first_zero, second_one))
    c0 = nonnegative_interval(multiply_interval(first_one, second_zero))
    d0 = nonnegative_interval(multiply_interval(first_one, second_one))
    a = nonnegative_interval(add_interval(a0, determinant))
    b = nonnegative_interval(subtract_interval(b0, determinant))
    c = nonnegative_interval(subtract_interval(c0, determinant))
    d = nonnegative_interval(add_interval(d0, determinant))
    bc = nonnegative_interval(multiply_interval(b, c))
    bc0 = nonnegative_interval(multiply_interval(b0, c0))
    bc_difference = multiply_interval(
        determinant,
        subtract_interval(determinant, add_interval(b0, c0)),
    )

    zero = np.zeros_like(a[0])
    result = (zero, zero)
    low_intersection = max(0, 2 * weight - degree)
    for intersection in range(low_intersection, weight + 1):
        cross = weight - intersection
        zero_zero = degree - 2 * weight + intersection
        coefficient = math.comb(weight, intersection) * math.comb(
            degree - weight, cross
        )
        coefficient_low, coefficient_high = rational_interval(coefficient, 1)
        baseline = (
            np.full_like(a[0], coefficient_low),
            np.full_like(a[1], coefficient_high),
        )
        baseline = multiply_interval(
            baseline, positive_power_interval(a0, zero_zero)
        )
        baseline = multiply_interval(
            baseline, positive_power_interval(bc0, cross)
        )
        baseline = multiply_interval(
            baseline, positive_power_interval(d0, intersection)
        )

        regular = np.ones_like(a[0], dtype=bool)
        if zero_zero:
            regular &= a0[0] > 0
        if cross:
            regular &= bc0[0] > 0
        if intersection:
            regular &= d0[0] > 0

        term = subtract_interval(
            multiply_interval(
                multiply_interval(
                    (
                        np.full_like(a[0], coefficient_low),
                        np.full_like(a[1], coefficient_high),
                    ),
                    positive_power_interval(a, zero_zero),
                ),
                positive_power_interval(bc, cross),
            ),
            baseline,
        )
        term = multiply_interval(
            term,
            positive_power_interval(d, intersection),
        ) if intersection else term
        # The direct path above has already subtracted a baseline without the
        # d0 factor when intersection is nonzero.  Recompute those uncommon
        # fallback entries below; regular entries use the stable ratio path.
        if intersection:
            current = (
                np.full_like(a[0], coefficient_low),
                np.full_like(a[1], coefficient_high),
            )
            current = multiply_interval(current, positive_power_interval(a, zero_zero))
            current = multiply_interval(current, positive_power_interval(bc, cross))
            current = multiply_interval(current, positive_power_interval(d, intersection))
            term = subtract_interval(current, baseline)

        if np.any(regular):
            local_zero = np.zeros_like(a[0])
            ratio_delta = (local_zero, local_zero)
            ratio_rows: list[tuple[Interval, int]] = []
            if zero_zero:
                ratio_rows.append((divide_by_positive(determinant, a0), zero_zero))
            if cross:
                ratio_rows.append((divide_by_positive(bc_difference, bc0), cross))
            if intersection:
                ratio_rows.append((divide_by_positive(determinant, d0), intersection))
            stable_regular = regular.copy()
            for ratio, _exponent in ratio_rows:
                stable_regular &= np.isfinite(ratio[0]) & np.isfinite(ratio[1])
                stable_regular &= np.maximum(np.abs(ratio[0]), np.abs(ratio[1])) <= 0.5
            for ratio, exponent in ratio_rows:
                masked_ratio = (
                    np.where(stable_regular, ratio[0], 0),
                    np.where(stable_regular, ratio[1], 0),
                )
                ratio_delta = delta_product(
                    ratio_delta, delta_power(masked_ratio, exponent)
                )
            stable = multiply_interval(baseline, ratio_delta)
            term = (
                np.where(stable_regular, stable[0], term[0]),
                np.where(stable_regular, stable[1], term[1]),
            )
        result = add_interval(result, term)
    return result


def evaluate(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    levels: list[int],
    progress_every: int,
    canonical_orbits: bool,
) -> dict:
    if right_degree % 2 != 1:
        raise ValueError("complement grouping requires odd row degree")
    denominator = math.comb(message_bits, right_degree)
    bias_numerators = [
        krawtchouk(message_bits, right_degree, weight)
        for weight in range(message_bits + 1)
    ]
    bias_low, bias_high = bias_intervals(message_bits, right_degree)
    totals = {level: np.float64(0) for level in levels}
    type_count = 0
    covered_type_count = 0
    started = time.perf_counter()
    one = None
    for n11 in range(message_bits + 1):
        n01_rows: list[int] = []
        n10_rows: list[int] = []
        orbit_sizes: list[int] = []
        for n01 in range(message_bits - n11 + 1):
            for n10 in range(message_bits - n11 - n01 + 1):
                row = (message_bits - n11 - n01 - n10, n01, n10, n11)
                orbit = pair_type_orbit(row)
                if canonical_orbits and row != min(orbit):
                    continue
                n01_rows.append(n01)
                n10_rows.append(n10)
                orbit_sizes.append(len(orbit) if canonical_orbits else 1)
        if not n01_rows:
            continue
        n01 = np.asarray(n01_rows, dtype=np.int16)
        n10 = np.asarray(n10_rows, dtype=np.int16)
        first_weight = n10 + n11
        second_weight = n01 + n11
        difference_weight = n01 + n10
        first_bias = (bias_low[first_weight], bias_high[first_weight])
        second_bias = (bias_low[second_weight], bias_high[second_weight])
        if one is None or len(one[0]) != len(n01):
            one = (np.ones(len(n01)), np.ones(len(n01)))
        first_one = (
            np.ldexp(subtract_interval(one, first_bias)[0], -1),
            np.ldexp(subtract_interval(one, first_bias)[1], -1),
        )
        second_one = (
            np.ldexp(subtract_interval(one, second_bias)[0], -1),
            np.ldexp(subtract_interval(one, second_bias)[1], -1),
        )
        determinant_rows = [
            rational_interval(
                bias_numerators[int(difference)] * denominator
                - bias_numerators[int(first)] * bias_numerators[int(second)],
                4 * denominator * denominator,
            )
            for first, second, difference in zip(
                first_weight, second_weight, difference_weight, strict=True
            )
        ]
        determinant = (
            np.asarray([row[0] for row in determinant_rows]),
            np.asarray([row[1] for row in determinant_rows]),
        )
        negative_determinant = (-determinant[1], -determinant[0])
        first_zero = subtract_interval(one, first_one)
        second_zero = subtract_interval(one, second_one)
        transforms = (
            (first_one, second_one, determinant),
            (first_zero, second_one, negative_determinant),
            (first_one, second_zero, negative_determinant),
            (first_zero, second_zero, determinant),
        )
        multiplicities = [
            orbit_size
            * math.comb(message_bits, n11)
            * math.comb(message_bits - n11, int(d01))
            * math.comb(message_bits - n11 - int(d01), int(d10))
            for d01, d10, orbit_size in zip(
                n01, n10, orbit_sizes, strict=True
            )
        ]
        multiplicity_upper = up(np.asarray(multiplicities, dtype=np.float64))
        type_scale_upper = np.ldexp(multiplicity_upper, message_bits - 2)
        for level in levels:
            orbit = (
                np.zeros_like(first_one[0]),
                np.zeros_like(first_one[1]),
            )
            for transformed_first, transformed_second, transformed_determinant in transforms:
                orbit = add_interval(
                    orbit,
                    centered_diagonal_interval(
                        transformed_first,
                        transformed_second,
                        transformed_determinant,
                        output_bits,
                        level,
                    ),
                )
            positive_upper = np.maximum(orbit[1], 0)
            contribution_upper = positive_multiply_upper(
                type_scale_upper, positive_upper
            )
            totals[level] = add_upper_scalar(
                totals[level], positive_sum_upper(contribution_upper)
            )
        type_count += len(n01)
        covered_type_count += sum(orbit_sizes)
        if progress_every and (
            n11 % progress_every == 0 or n11 == message_bits
        ):
            print(
                f"progress,n11,{n11},{message_bits},pair_types,{type_count},elapsed_seconds,{time.perf_counter()-started:.3f}",
                flush=True,
            )
    return {
        "message_bits": message_bits,
        "output_bits": output_bits,
        "right_degree": right_degree,
        "enumerated_pair_types": type_count,
        "covered_pair_types": covered_type_count,
        "canonical_complement_swap_orbits": canonical_orbits,
        "entries": [
            {
                "sector": 0,
                "level": level,
                "scaled_diagonal_upper": float(totals[level]),
            }
            for level in levels
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--right-degree", type=int, default=33)
    parser.add_argument("--level", type=int, action="append", required=True)
    parser.add_argument("--progress-every", type=int, default=16)
    parser.add_argument("--no-canonical-orbits", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    levels = sorted(set(args.level))
    if any(level < 0 or level > args.output_bits for level in levels):
        parser.error("level is outside the output range")
    result = evaluate(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        levels,
        args.progress_every,
        not args.no_canonical_orbits,
    )
    payload = {
        "schema": "pure-ea-sector-zero-complement-orbits-outward-v1",
        "status": "OUTWARD_BINARY64_UPPER_BOUND",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
        },
        "entries": result["entries"],
        "enumerated_pair_types": result["enumerated_pair_types"],
        "covered_pair_types": result["covered_pair_types"],
        "canonical_complement_swap_orbits": result[
            "canonical_complement_swap_orbits"
        ],
        "elapsed_seconds": result["elapsed_seconds"],
        "scope": [
            "Complement-orbit averaging is exact for odd row degree.",
            "Every centered orbit is enclosed with directed binary64 intervals.",
            "The bound sums positive orbit upper endpoints and does not assume orbit positivity.",
            "The receipt bounds only the listed sector-zero diagonals.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
