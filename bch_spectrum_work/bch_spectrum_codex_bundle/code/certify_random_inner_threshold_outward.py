#!/usr/bin/env python3
"""Certify the RandomStepConv application threshold with outward rounding.

The transfer recurrences use positive arithmetic.  Every binary64 operation is
rounded upward by one ulp, while transcendental operations are enclosed by
Arb.  The final aggregation is exact rational arithmetic.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import arb, ctx


ROOT = Path(__file__).resolve().parents[1]
DIAGNOSTIC_PATH = ROOT / "generated" / "random_inner_oa15_majorant_diagnostic.json"
OUTPUT_PATH = ROOT / "generated" / "random_inner_threshold_outward.json"
JOHNSON_PATHS = {
    52: ROOT / "generated" / "johnson_n256_w52_d38.json",
    54: ROOT / "generated" / "johnson_n256_w54_d38.json",
}

MESSAGE_BITS = 1 << 20
OUTER_ROWS = 8192
OUTPUT_BITS = 1 << 21
DISTANCE = (OUTPUT_BITS + 9) // 10
MEMORY_BITS = 22
INNER_POSITIONS = 256
LOW_WEIGHTS = tuple(range(38, 52, 2))
HIGH_WEIGHTS = tuple(range(52, 130, 2))
INTERIOR_WEIGHTS = tuple(range(38, 220, 2))
TARGET = Fraction(1, 1 << 40)
LOCATOR_SCALE = Fraction(3968, 19)
RANDOM_RANK_LOCATOR_REFERENCE = Fraction(
    125188836519193760215181565,
    41658296553177088,
)
DEGENERATE_LOCATOR_CAP = 87550900

ctx.prec = 192


def up(value: float) -> float:
    if value == 0.0 or not math.isfinite(value):
        return value
    return math.nextafter(float(value), math.inf)


def up_array(values: np.ndarray) -> np.ndarray:
    return np.where(values == 0.0, 0.0, np.nextafter(values, np.inf))


def arb_of_float(value: float) -> arb:
    numerator, denominator = float(value).as_integer_ratio()
    return arb(numerator) / denominator


def upper_float(value: arb) -> float:
    _, upper = value.lower(), value.upper()
    result = float(upper)
    while arb_of_float(result) < upper:
        result = math.nextafter(result, math.inf)
    return result


def mul_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    product = left * right
    rounded = up_array(product)
    positive_underflow = (product == 0.0) & (left > 0.0) & (right > 0.0)
    return np.where(positive_underflow, np.nextafter(0.0, np.inf), rounded)


def add_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return up_array(left + right)


def ldexp_up(values: np.ndarray, exponents: np.ndarray) -> np.ndarray:
    scaled = np.ldexp(values, exponents)
    rounded = up_array(scaled)
    positive_underflow = (scaled == 0.0) & (values > 0.0)
    return np.where(positive_underflow, np.nextafter(0.0, np.inf), rounded)


def matmul_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Upper-bound a nonnegative 2-by-2 product, including batches."""
    result = np.empty(np.broadcast_shapes(left.shape, right.shape), dtype=np.float64)
    result[..., 0, 0] = add_up(
        mul_up(left[..., 0, 0], right[..., 0, 0]),
        mul_up(left[..., 0, 1], right[..., 1, 0]),
    )
    result[..., 0, 1] = add_up(
        mul_up(left[..., 0, 0], right[..., 0, 1]),
        mul_up(left[..., 0, 1], right[..., 1, 1]),
    )
    result[..., 1, 0] = add_up(
        mul_up(left[..., 1, 0], right[..., 0, 0]),
        mul_up(left[..., 1, 1], right[..., 1, 0]),
    )
    result[..., 1, 1] = add_up(
        mul_up(left[..., 1, 0], right[..., 0, 1]),
        mul_up(left[..., 1, 1], right[..., 1, 1]),
    )
    return result


def transition_upper(surprisal: float, memory_bits: int) -> tuple[np.ndarray, np.ndarray]:
    s = arb_of_float(surprisal)
    z = (-s).exp()
    collision = arb(1) / (1 << memory_bits)
    bit_moment = (1 + z) / 2
    terminate = upper_float(collision * bit_moment)
    survive = upper_float((1 - collision) * bit_moment)
    zero = np.asarray(((1.0, 0.0), (terminate, survive)), dtype=np.float64)
    active = np.asarray(
        ((terminate, survive), (terminate, survive)), dtype=np.float64
    )
    return zero, active


def uniform_coefficients_upper(
    zero: np.ndarray,
    selected: np.ndarray,
    positions: int,
    maximum_degree: int,
) -> tuple[np.ndarray, np.ndarray]:
    state_count = zero.shape[0]
    mantissas = np.zeros((maximum_degree + 1, state_count, state_count), dtype=np.float64)
    exponents = np.full(maximum_degree + 1, np.iinfo(np.int64).min, dtype=np.int64)
    mantissas[0] = np.eye(state_count, dtype=np.float64)
    exponents[0] = 0

    for completed in range(positions):
        maximum_next_degree = min(completed + 1, maximum_degree)
        maximum_current_degree = min(completed, maximum_degree)
        count = maximum_next_degree + 1
        zero_terms = np.zeros((count, state_count, state_count), dtype=np.float64)
        selected_terms = np.zeros_like(zero_terms)
        zero_exponents = np.full(count, np.iinfo(np.int64).min, dtype=np.int64)
        selected_exponents = np.full_like(zero_exponents, np.iinfo(np.int64).min)

        zero_count = maximum_current_degree + 1
        zero_products = matmul_up(mantissas[:zero_count], zero)
        zero_weights = up_array(
            (completed + 1 - np.arange(zero_count, dtype=np.float64)) / (completed + 1)
        )
        zero_terms[:zero_count] = mul_up(zero_products, zero_weights[:, None, None])
        zero_exponents[:zero_count] = exponents[:zero_count]

        selected_count = min(maximum_current_degree + 1, maximum_next_degree)
        if selected_count:
            selected_products = matmul_up(mantissas[:selected_count], selected)
            selected_weights = up_array(
                np.arange(1, selected_count + 1, dtype=np.float64) / (completed + 1)
            )
            selected_terms[1:selected_count + 1] = mul_up(
                selected_products, selected_weights[:, None, None]
            )
            selected_exponents[1:selected_count + 1] = exponents[:selected_count]

        common_exponents = np.maximum(zero_exponents, selected_exponents)
        zero_shifts = np.zeros(count, dtype=np.int64)
        selected_shifts = np.zeros(count, dtype=np.int64)
        zero_valid = zero_exponents != np.iinfo(np.int64).min
        selected_valid = selected_exponents != np.iinfo(np.int64).min
        zero_shifts[zero_valid] = (
            zero_exponents[zero_valid] - common_exponents[zero_valid]
        )
        selected_shifts[selected_valid] = (
            selected_exponents[selected_valid] - common_exponents[selected_valid]
        )
        total = add_up(
            ldexp_up(zero_terms, zero_shifts[:, None, None]),
            ldexp_up(selected_terms, selected_shifts[:, None, None]),
        )
        maxima = np.max(total, axis=(1, 2))
        _, shifts = np.frexp(maxima)
        mantissas[:count] = ldexp_up(total, -shifts[:, None, None])
        exponents[:count] = common_exponents + shifts

    return mantissas, exponents


def row_sum_log2_upper(matrix: np.ndarray, exponent: int) -> float:
    total = up(float(matrix[0, 0]) + float(matrix[0, 1]))
    logarithm = upper_float(arb_of_float(total).log() / arb(2).log())
    return up(logarithm + exponent)


def coefficient_logs_for_tilt(log_tilt: float) -> list[float]:
    surprisal = math.exp(log_tilt)
    zero, active = transition_upper(surprisal, MEMORY_BITS)
    bit_mantissas, bit_exponents = uniform_coefficients_upper(
        zero, active, OUTER_ROWS, 1
    )
    assert all(-900 < int(exponent) < 900 for exponent in bit_exponents)
    bit_zero = np.ldexp(bit_mantissas[0], int(bit_exponents[0]))
    bit_one = np.ldexp(bit_mantissas[1], int(bit_exponents[1]))
    inner_mantissas, inner_exponents = uniform_coefficients_upper(
        bit_zero, bit_one, INNER_POSITIONS, INNER_POSITIONS
    )
    correction = upper_float(arb(DISTANCE) * arb_of_float(surprisal) / arb(2).log())
    logs: list[float] = []
    for weight in range(INNER_POSITIONS + 1):
        moment = row_sum_log2_upper(inner_mantissas[weight], int(inner_exponents[weight]))
        tail = up(moment + correction)
        logs.append(up(13.0 + min(0.0, tail)))
    return logs


def pow2_upper(logarithm: float) -> float:
    value = (arb(2).log() * arb_of_float(logarithm)).exp()
    return upper_float(value)


def fraction_log2(value: Fraction) -> float:
    return math.log2(value.numerator) - math.log2(value.denominator)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    diagnostic = json.loads(DIAGNOSTIC_PATH.read_text())
    diagnostic_tilts = diagnostic["randomstepconv"]["best_tilt_by_weight"]
    tilt_by_weight = {
        weight: float(diagnostic_tilts[weight])
        for weight in INTERIOR_WEIGHTS
    }
    unique_tilts = sorted(set(tilt_by_weight.values()))
    tables = {str(tilt): coefficient_logs_for_tilt(tilt) for tilt in unique_tilts}

    coefficient_log2_upper: dict[int, float] = {}
    coefficient_upper: dict[int, Fraction] = {}
    for weight, tilt in tilt_by_weight.items():
        logarithm = tables[str(tilt)][weight]
        coefficient_log2_upper[weight] = logarithm
        coefficient_upper[weight] = Fraction.from_float(pow2_upper(logarithm))

    johnson_caps = {
        weight: int(json.loads(path.read_text())["integer_shell_upper"])
        for weight, path in JOHNSON_PATHS.items()
    }

    high_functional_upper = Fraction(0)
    high_rows = []
    for weight in HIGH_WEIGHTS:
        subset_size = weight - 18
        packing_upper = math.comb(256, subset_size) // math.comb(weight, subset_size)
        multiplicity_upper = min(1 << 128, packing_upper)
        if weight in johnson_caps:
            multiplicity_upper = min(multiplicity_upper, johnson_caps[weight])
        pair_coefficient = coefficient_upper[weight]
        if weight != 128:
            pair_coefficient += coefficient_upper[256 - weight]
        contribution = multiplicity_upper * pair_coefficient
        high_functional_upper += contribution
        high_rows.append(
            {
                "weight": weight,
                "multiplicity_upper": str(multiplicity_upper),
                "pair_coefficient_log2_upper": fraction_log2(pair_coefficient),
                "contribution_log2_upper": fraction_log2(contribution),
            }
        )

    residual_lower = TARGET - high_functional_upper
    if residual_lower <= 0:
        raise RuntimeError(
            "outward high-weight bound exceeds target: "
            f"log2(high)={fraction_log2(high_functional_upper):.12f}, "
            f"target=-40"
        )
    remaining_codewords = (1 << 128) - 2
    normalization = sum(math.comb(256, weight) for weight in INTERIOR_WEIGHTS)
    low_null_functional_upper = Fraction(0)
    low_rows = []
    null_count_38 = None
    for weight in LOW_WEIGHTS:
        null_count = Fraction(remaining_codewords * math.comb(256, weight), normalization)
        if weight == 38:
            null_count_38 = null_count
        pair_coefficient = coefficient_upper[weight] + coefficient_upper[256 - weight]
        contribution = null_count * pair_coefficient
        low_null_functional_upper += contribution
        low_rows.append(
            {
                "weight": weight,
                "null_count_log2": fraction_log2(null_count),
                "pair_coefficient_log2_upper": fraction_log2(pair_coefficient),
                "contribution_log2_upper": fraction_log2(contribution),
            }
        )
    assert null_count_38 is not None

    joint_multiplier_lower = residual_lower / low_null_functional_upper
    a38_cap_lower = null_count_38 * joint_multiplier_lower
    integer_a38_cap = a38_cap_lower.numerator // a38_cap_lower.denominator
    locator_cap_lower = Fraction(integer_a38_cap, 1) / LOCATOR_SCALE
    r6_threshold_lower = locator_cap_lower / RANDOM_RANK_LOCATOR_REFERENCE
    smooth_r6_threshold_lower = (
        locator_cap_lower - DEGENERATE_LOCATOR_CAP
    ) / RANDOM_RANK_LOCATOR_REFERENCE

    constants = {}
    for constant in (5, 6):
        candidate_locator = constant * RANDOM_RANK_LOCATOR_REFERENCE
        candidate_a38 = LOCATOR_SCALE * candidate_locator
        constants[str(constant)] = {
            "candidate_a38_log2": fraction_log2(candidate_a38),
            "certified_below_integer_cap": candidate_a38 <= integer_a38_cap,
            "margin_bits": fraction_log2(Fraction(integer_a38_cap, 1) / candidate_a38),
        }

    receipt = {
        "classification": (
            "directed outward certificate: Arb transcendental enclosures, "
            "one-ULP positive binary64 recurrences, and exact rational aggregation"
        ),
        "parameters": {
            "message_bits": MESSAGE_BITS,
            "outer_rows": OUTER_ROWS,
            "output_bits": OUTPUT_BITS,
            "distance": DISTANCE,
            "memory_bits": MEMORY_BITS,
            "target": "1/1099511627776",
            "unique_diagnostic_tilts": unique_tilts,
            "tilt_interpretation": (
                "Each exp(best_log_tilt) binary64 value is treated as an exact "
                "dyadic Chernoff witness; optimality of the diagnostic grid is "
                "not assumed."
            ),
        },
        "source_sha256": {
            "diagnostic": sha256(DIAGNOSTIC_PATH),
            "johnson_w52": sha256(JOHNSON_PATHS[52]),
            "johnson_w54": sha256(JOHNSON_PATHS[54]),
        },
        "high_functional_upper": {
            "numerator": str(high_functional_upper.numerator),
            "denominator": str(high_functional_upper.denominator),
            "log2": fraction_log2(high_functional_upper),
            "rows": high_rows,
        },
        "residual_lower": {
            "numerator": str(residual_lower.numerator),
            "denominator": str(residual_lower.denominator),
            "log2": fraction_log2(residual_lower),
        },
        "low_null_functional_upper": {
            "numerator": str(low_null_functional_upper.numerator),
            "denominator": str(low_null_functional_upper.denominator),
            "log2": fraction_log2(low_null_functional_upper),
            "rows": low_rows,
        },
        "joint_multiplier_lower": {
            "numerator": str(joint_multiplier_lower.numerator),
            "denominator": str(joint_multiplier_lower.denominator),
            "log2": fraction_log2(joint_multiplier_lower),
        },
        "a38_sufficient_cap": {
            "strict_real_lower_log2": fraction_log2(a38_cap_lower),
            "certified_integer_cap": str(integer_a38_cap),
            "integer_cap_log2": math.log2(integer_a38_cap),
        },
        "locator_cap_lower": {
            "numerator": str(locator_cap_lower.numerator),
            "denominator": str(locator_cap_lower.denominator),
            "log2": fraction_log2(locator_cap_lower),
        },
        "r6_factor_threshold_lower": {
            "numerator": str(r6_threshold_lower.numerator),
            "denominator": str(r6_threshold_lower.denominator),
            "decimal": float(r6_threshold_lower),
            "log2": fraction_log2(r6_threshold_lower),
        },
        "smooth_r6_factor_threshold_after_degenerate_cap": {
            "degenerate_locator_cap": DEGENERATE_LOCATOR_CAP,
            "numerator": str(smooth_r6_threshold_lower.numerator),
            "denominator": str(smooth_r6_threshold_lower.denominator),
            "decimal": float(smooth_r6_threshold_lower),
        },
        "constant_checks": constants,
        "coefficient_log2_upper": {
            str(weight): coefficient_log2_upper[weight]
            for weight in sorted(coefficient_log2_upper)
        },
    }
    OUTPUT_PATH.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({
        "output": str(OUTPUT_PATH),
        "unique_tilts": len(unique_tilts),
        "integer_a38_cap": integer_a38_cap,
        "r6_factor_threshold_lower": float(r6_threshold_lower),
        "constant_6_certified": constants["6"]["certified_below_integer_cap"],
        "constant_6_margin_bits": constants["6"]["margin_bits"],
    }, indent=2))


if __name__ == "__main__":
    main()
