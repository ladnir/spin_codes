#!/usr/bin/env python3
"""Outward-rounded certificate for the explicit-XBCH Riffle code.

Arb encloses transcendental values and exact combinatorial ratios.  Batched
matrix recurrences use binary64 upper endpoints.  Every nonnegative binary64
addition, multiplication, and division is followed by nextafter(+infinity).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from flint import arb, ctx, fmpz_poly


OUTER_BITS = 1024
OUTER_DIMENSION = 512
OUTER_BLOCKS = 2048
STEP_BITS = 128
STATE_BITS = 16
EPOCHS = 16
OUTPUT_BITS = 1 << 21
DISTANCE = math.floor(0.09 * OUTPUT_BITS)
TAIL_ENDPOINT = 97
LOG_SURPRISALS = (-10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0)
BODY_LOG_SURPRISALS = (
    -7.0,
    -6.5,
    -6.0,
    -5.0,
    -4.0,
    -3.0,
    -2.5,
    -2.25,
    -2.0,
    -1.75,
    -1.5,
    -1.25,
    -1.0,
    -0.75,
    -0.5,
    -0.25,
    0.0,
    0.25,
    0.5,
    0.75,
)


def up(value: np.ndarray | float) -> np.ndarray | float:
    return np.nextafter(value, math.inf)


def add_up(left: np.ndarray | float, right: np.ndarray | float):
    raw = left + right
    rounded = up(raw)
    return np.where((np.asarray(left) == 0.0) & (np.asarray(right) == 0.0), 0.0, rounded)


def mul_up(left: np.ndarray | float, right: np.ndarray | float):
    raw = left * right
    rounded = up(raw)
    return np.where((np.asarray(left) == 0.0) | (np.asarray(right) == 0.0), 0.0, rounded)


def div_up(left: np.ndarray | float, right: np.ndarray | float):
    raw = left / right
    rounded = up(raw)
    return np.where(np.asarray(left) == 0.0, 0.0, rounded)


def arb_upper_float(value: arb) -> float:
    return float(np.nextafter(float(value.upper()), math.inf))


def arb_lower_float(value: arb) -> float:
    return float(np.nextafter(float(value.lower()), -math.inf))


def arb_from_float(value: float) -> arb:
    # Python-flint converts the binary64 value exactly at the selected
    # precision because every binary64 significand has at most 53 bits.
    return arb(value)


def log2_upper(value: arb) -> float:
    return arb_upper_float(value.log() / arb(2).log())


def read_local_spectrum(path: Path) -> list[int]:
    values = [0] * 65
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            values[int(row["weight"])] = int(row["count"])
    if sum(values) != 1 << 32 or values[0] != 1:
        raise ValueError("constituent spectrum has the wrong mass")
    return values


def transition_upper() -> np.ndarray:
    """Return an entrywise upper accumulator transition matrix."""
    transition = np.zeros((OUTER_BITS + 1, OUTER_BITS + 1))
    transition[0, 0] = 1.0
    for input_weight in range(1, OUTER_BITS + 1):
        half_down = input_weight // 2
        half_up = (input_weight + 1) // 2
        denominator = math.comb(OUTER_BITS, input_weight)
        low = half_up
        high = OUTER_BITS - half_down
        for output_weight in range(low, high + 1):
            numerator = (
                math.comb(OUTER_BITS - output_weight, half_down)
                * math.comb(output_weight - 1, half_up - 1)
            )
            transition[input_weight, output_weight] = arb_upper_float(
                arb(numerator) / denominator
            )
    return transition


def apply_transition_upper(
    spectrum: np.ndarray, transition: np.ndarray
) -> np.ndarray:
    result = np.zeros_like(spectrum)
    for input_weight in range(OUTER_BITS + 1):
        if spectrum[input_weight] == 0.0:
            continue
        contribution = mul_up(
            spectrum[input_weight], transition[input_weight]
        )
        result = add_up(result, contribution)
    return result


def outer_spectrum_upper(local_path: Path) -> tuple[np.ndarray, dict[str, object]]:
    local = read_local_spectrum(local_path)
    direct = fmpz_poly(local) ** 16
    spectrum = np.zeros(OUTER_BITS + 1)
    for weight in range(OUTER_BITS + 1):
        coefficient = int(direct[weight])
        if coefficient:
            spectrum[weight] = arb_upper_float(arb(coefficient))
    transition = transition_upper()
    first = apply_transition_upper(spectrum, transition)
    second = apply_transition_upper(first, transition)
    mass_upper = arb(0)
    for value in second:
        mass_upper += arb_from_float(float(value))
    return second, {
        "direct_sum": "exact fmpz polynomial power",
        "accumulator_transition": (
            "exact integer ratios enclosed by 192-bit Arb; upper binary64 "
            "endpoints"
        ),
        "accumulator_composition": (
            "nonnegative binary64 products and sums, each rounded upward"
        ),
        "log2_mass_upper": log2_upper(mass_upper),
    }


def load_activation(path: Path) -> list[arb]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    values: list[arb | None] = [None] * (STEP_BITS + 1)
    values[0] = arb(1)
    for row in payload["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight > STEP_BITS:
            continue
        key = "uniform_support_average_distinct_upper_bound"
        if key not in row:
            key = "uniform_256_support_average_distinct_upper_bound"
        # The serialized decimal is itself the proved upper bound supplied by
        # the activation certificate.
        values[weight] = arb(str(row[key]))
    if any(value is None for value in values):
        raise ValueError("activation certificate is incomplete")
    return [value for value in values if value is not None]


def load_live_counts(path: Path) -> list[int]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    counts = [0] * (STEP_BITS + 1)
    for row in payload["spectrum"]:
        counts[int(row["weight"])] = int(row["count"])
    if sum(counts) != 1 << STATE_BITS or counts[0] != 1:
        raise ValueError("live-state spectrum has the wrong mass")
    counts[0] = 0
    return counts


def support_moment(
    z: arb, input_weight: int, codeword_weight: int
) -> arb:
    denominator = math.comb(STEP_BITS, input_weight)
    result = arb(0)
    low = max(0, input_weight - (STEP_BITS - codeword_weight))
    high = min(input_weight, codeword_weight)
    for intersection in range(low, high + 1):
        numerator = (
            math.comb(codeword_weight, intersection)
            * math.comb(
                STEP_BITS - codeword_weight,
                input_weight - intersection,
            )
        )
        exponent = input_weight + codeword_weight - 2 * intersection
        result += (arb(numerator) / denominator) * z**exponent
    return result


def impulse_upper(
    z: arb,
    activation: list[arb],
    live_counts: list[int],
    maximum_input_weight: int,
) -> np.ndarray:
    denominator = (1 << STATE_BITS) - 1
    punctured = arb(denominator) / (denominator - 1)
    moments: list[arb] = []
    for input_weight in range(maximum_input_weight + 1):
        moment = arb(0)
        for codeword_weight, count in enumerate(live_counts):
            if count:
                moment += (
                    arb(count)
                    / denominator
                    * support_moment(z, input_weight, codeword_weight)
                )
        moments.append(moment)

    matrices = np.zeros((maximum_input_weight + 1, 2, 2))
    matrices[0, 0, 0] = 1.0
    matrices[0, 1, 1] = arb_upper_float(punctured * moments[0])
    for weight in range(1, maximum_input_weight + 1):
        output = z**weight
        averaged = punctured * moments[weight]
        matrices[weight, 0, 0] = arb_upper_float(
            activation[weight] * output
        )
        matrices[weight, 0, 1] = arb_upper_float(output)
        matrices[weight, 1, 0] = arb_upper_float(
            averaged / denominator
        )
        matrices[weight, 1, 1] = arb_upper_float(averaged)
    return matrices


def matmul_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty_like(left)
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


def marked_region_upper(
    zero_epoch: np.ndarray, active_epoch: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    powers = [np.eye(2)]
    for _ in range(EPOCHS):
        powers.append(matmul_up(powers[-1], zero_epoch))
    active = np.zeros((2, 2))
    for epoch in range(EPOCHS):
        term = matmul_up(
            matmul_up(powers[epoch], active_epoch),
            powers[EPOCHS - 1 - epoch],
        )
        active = add_up(active, term)
    active = div_up(active, float(EPOCHS))
    return powers[EPOCHS], active


def subset_moments_upper(
    zero_row: np.ndarray, active_row: np.ndarray
) -> np.ndarray:
    """Average products over uniform subsets, with upward rounding.

    The normalized recurrence avoids binomial coefficients and underflow from
    dividing a large coefficient at the end.
    """
    current = np.zeros((1, 2, 2))
    current[0] = np.eye(2)
    for old_rows in range(OUTER_BITS):
        next_rows = old_rows + 1
        updated = np.zeros((next_rows + 1, 2, 2))
        zero_products = matmul_up(current, zero_row)
        active_products = matmul_up(current, active_row)
        weights = np.arange(old_rows + 1, dtype=np.float64)
        zero_prob = up((next_rows - weights) / next_rows)
        active_prob = up((weights + 1.0) / next_rows)
        updated[: old_rows + 1] = mul_up(
            zero_products, zero_prob[:, None, None]
        )
        active_terms = mul_up(
            active_products, active_prob[:, None, None]
        )
        updated[1:] = add_up(updated[1:], active_terms)
        current = updated
    return add_up(current[:, 0, 0], current[:, 0, 1])


def chernoff_log2_upper(moment: float, surprisal: arb) -> float:
    if moment <= 0.0:
        raise ArithmeticError("outward recurrence produced a zero moment")
    value = arb_from_float(moment).log() + DISTANCE * surprisal
    upper = arb_upper_float(value / arb(2).log())
    return min(0.0, upper)


def arb_logsum2_upper(log2_values: list[float]) -> float:
    pivot = max(log2_values)
    total = arb(0)
    two = arb(2)
    for value in log2_values:
        total += two ** (arb(value) - arb(pivot))
    return arb_upper_float(arb(pivot) + total.log() / two.log())


def one_tail_upper(
    outer: np.ndarray,
    activation: list[arb],
    live_counts: list[int],
) -> tuple[float, list[dict[str, float]]]:
    best = np.full(OUTER_BITS + 1, math.inf)
    best_tilt = np.full(OUTER_BITS + 1, math.nan)
    for log_surprisal in LOG_SURPRISALS:
        surprisal = arb(log_surprisal).exp()
        z = (-surprisal).exp()
        impulse = impulse_upper(z, activation, live_counts, 1)
        zero_row, active_row = marked_region_upper(impulse[0], impulse[1])
        moments = subset_moments_upper(zero_row, active_row)
        for weight, moment in enumerate(moments):
            candidate = chernoff_log2_upper(float(moment), surprisal)
            if candidate < best[weight]:
                best[weight] = candidate
                best_tilt[weight] = log_surprisal

    rows: list[dict[str, float]] = []
    contribution_logs: list[float] = []
    for weight in range(1, OUTER_BITS + 1):
        if min(weight, OUTER_BITS - weight) > TAIL_ENDPOINT:
            continue
        if outer[weight] == 0.0:
            continue
        value = (
            arb(OUTER_BLOCKS)
            * arb_from_float(float(outer[weight]))
            * arb(2) ** arb(float(best[weight]))
        )
        log_value = log2_upper(value)
        contribution_logs.append(log_value)
        rows.append(
            {
                "outer_weight": weight,
                "best_log_surprisal": float(best_tilt[weight]),
                "inner_log2_upper": float(best[weight]),
                "pointwise_log2_upper": log_value,
                "pointwise_margin_bits": -log_value,
            }
        )
    aggregate = arb_logsum2_upper(contribution_logs)
    rows.sort(key=lambda row: row["pointwise_log2_upper"], reverse=True)
    return aggregate, rows


def region_fixed_count_upper(impulse: np.ndarray) -> np.ndarray:
    current = impulse[:3].copy()
    for completed in range(1, EPOCHS):
        updated = np.zeros((3, 2, 2))
        total_bits = (completed + 1) * STEP_BITS
        previous_bits = completed * STEP_BITS
        for total_count in range(3):
            denominator = math.comb(total_bits, total_count)
            for next_count in range(total_count + 1):
                previous_count = total_count - next_count
                probability = arb_upper_float(
                    arb(
                        math.comb(STEP_BITS, next_count)
                        * math.comb(previous_bits, previous_count)
                    )
                    / denominator
                )
                term = mul_up(
                    matmul_up(
                        current[previous_count], impulse[next_count]
                    ),
                    probability,
                )
                updated[total_count] = add_up(
                    updated[total_count], term
                )
        current = updated
    return current


def matrix_power_1024_upper(matrices: np.ndarray) -> np.ndarray:
    result = matrices
    for _ in range(10):
        result = matmul_up(result, result)
    return result


def fixed_iid_probability_lower(weight: int) -> float:
    if weight == OUTER_BITS:
        return 1.0
    probability = arb(weight) / OUTER_BITS
    value = (
        arb(math.comb(OUTER_BITS, weight))
        * probability**weight
        * (1 - probability) ** (OUTER_BITS - weight)
    )
    lower = arb_lower_float(value)
    if lower <= 0.0:
        raise ArithmeticError("conditioning probability lost its lower bound")
    return lower


def two_tail_upper(
    outer: np.ndarray,
    activation: list[arb],
    live_counts: list[int],
) -> tuple[float, dict[str, float]]:
    weights = np.asarray(
        [
            weight
            for weight in range(1, OUTER_BITS + 1)
            if min(weight, OUTER_BITS - weight) <= TAIL_ENDPOINT
            and outer[weight] > 0.0
        ],
        dtype=np.int64,
    )
    first = np.repeat(weights, len(weights))
    second = np.tile(weights, len(weights))
    p1 = first.astype(np.float64) / OUTER_BITS
    p2 = second.astype(np.float64) / OUTER_BITS
    q0 = (1.0 - p1) * (1.0 - p2)
    q1 = p1 * (1.0 - p2) + (1.0 - p1) * p2
    q2 = p1 * p2

    surprisal = arb(-10).exp()
    z = (-surprisal).exp()
    impulse = impulse_upper(z, activation, live_counts, 2)
    regions = region_fixed_count_upper(impulse)
    mixed = np.zeros((len(first), 2, 2))
    for probability, region in zip((q0, q1, q2), regions):
        mixed = add_up(
            mixed,
            mul_up(probability[:, None, None], region),
        )
    powered = matrix_power_1024_upper(mixed)
    moments = add_up(powered[:, 0, 0], powered[:, 0, 1])
    conditioning = {
        int(weight): fixed_iid_probability_lower(int(weight))
        for weight in weights
    }
    location = math.comb(OUTER_BLOCKS, 2)
    contribution_logs: list[float] = []
    dominant: dict[str, float] | None = None
    for index, (weight1, weight2) in enumerate(zip(first, second)):
        inner = chernoff_log2_upper(float(moments[index]), surprisal)
        value = (
            arb(location)
            * arb_from_float(float(outer[int(weight1)]))
            * arb_from_float(float(outer[int(weight2)]))
            * arb(2) ** arb(inner)
            / arb_from_float(conditioning[int(weight1)])
            / arb_from_float(conditioning[int(weight2)])
        )
        log_value = log2_upper(value)
        contribution_logs.append(log_value)
        if dominant is None or log_value > dominant["pointwise_log2_upper"]:
            dominant = {
                "first_weight": int(weight1),
                "second_weight": int(weight2),
                "inner_log2_upper": inner,
                "pointwise_log2_upper": log_value,
                "pointwise_margin_bits": -log_value,
            }
    if dominant is None:
        raise AssertionError("two-tail sum is empty")
    return arb_logsum2_upper(contribution_logs), dominant


def three_plus_tail_upper(outer: np.ndarray) -> float:
    tail_count = arb(0)
    for weight in range(1, OUTER_BITS + 1):
        if min(weight, OUTER_BITS - weight) <= TAIL_ENDPOINT:
            tail_count += arb_from_float(float(outer[weight]))
    # For k >= 3,
    # C(N,k) <= C(N,3) C(N-3,k-3).  This avoids cancellation in the exact
    # binomial tail and gives a close nonnegative upper bound.
    value = (
        arb(math.comb(OUTER_BLOCKS, 3))
        * tail_count**3
        * (1 + tail_count) ** (OUTER_BLOCKS - 3)
    )
    return log2_upper(value)


def region_combination_weights_upper() -> list[list[np.ndarray]]:
    """Precompute hypergeometric epoch-combination probabilities."""
    stages: list[list[np.ndarray]] = []
    for completed in range(1, EPOCHS):
        previous_bits = completed * STEP_BITS
        total_bits = (completed + 1) * STEP_BITS
        by_next_count: list[np.ndarray] = []
        for next_count in range(STEP_BITS + 1):
            previous_counts = np.arange(previous_bits + 1, dtype=np.int64)
            values = np.zeros(previous_bits + 1)
            epoch_choices = math.comb(STEP_BITS, next_count)
            for index, previous_count in enumerate(previous_counts):
                total_count = int(previous_count) + next_count
                values[index] = arb_upper_float(
                    arb(
                        epoch_choices
                        * math.comb(previous_bits, int(previous_count))
                    )
                    / math.comb(total_bits, total_count)
                )
            by_next_count.append(values)
        stages.append(by_next_count)
    return stages


def fixed_region_matrices_upper(
    impulse: np.ndarray,
    combination_weights: list[list[np.ndarray]],
) -> np.ndarray:
    """Condition on every possible number of active bits in one region."""
    current = impulse.copy()
    for completed, by_next_count in enumerate(combination_weights, start=1):
        previous_bits = completed * STEP_BITS
        total_bits = (completed + 1) * STEP_BITS
        updated = np.zeros((total_bits + 1, 2, 2))
        for next_count, probabilities in enumerate(by_next_count):
            products = matmul_up(current, impulse[next_count])
            terms = mul_up(products, probabilities[:, None, None])
            destination = slice(next_count, next_count + previous_bits + 1)
            updated[destination] = add_up(updated[destination], terms)
        current = updated
    return current


def normalize_scaled_upper(
    matrices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    maximum = np.max(matrices, axis=(-2, -1))
    nonzero = maximum > 0.0
    exponents = np.zeros(maximum.shape, dtype=np.int64)
    if np.any(nonzero):
        _, selected = np.frexp(maximum[nonzero])
        exponents[nonzero] = selected.astype(np.int64)
    scaled = np.zeros_like(matrices)
    if matrices.ndim == 2:
        if bool(nonzero):
            raw = np.ldexp(matrices, -int(exponents))
            scaled = np.where(matrices == 0.0, 0.0, up(raw))
    else:
        raw = np.ldexp(matrices, -exponents[..., None, None])
        scaled = np.where(matrices == 0.0, 0.0, up(raw))
    return scaled, exponents


def ldexp_up(values: np.ndarray, shifts: np.ndarray) -> np.ndarray:
    raw = np.ldexp(values, shifts[..., None, None])
    return np.where(values == 0.0, 0.0, up(raw))


def scaled_matmul_upper(
    left: np.ndarray,
    left_exp: np.ndarray,
    right: np.ndarray,
    right_exp: np.ndarray | np.int64,
) -> tuple[np.ndarray, np.ndarray]:
    product = matmul_up(left, right)
    scaled, extra = normalize_scaled_upper(product)
    return scaled, left_exp + right_exp + extra


def scaled_add_upper(
    left: np.ndarray,
    left_exp: np.ndarray,
    right: np.ndarray,
    right_exp: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    left_zero = np.max(left, axis=(-2, -1)) == 0.0
    right_zero = np.max(right, axis=(-2, -1)) == 0.0
    target = np.maximum(left_exp, right_exp)
    left_aligned = ldexp_up(left, left_exp - target)
    right_aligned = ldexp_up(right, right_exp - target)
    combined = add_up(left_aligned, right_aligned)
    scaled, extra = normalize_scaled_upper(combined)
    exponent = target + extra
    only_left = ~left_zero & right_zero
    only_right = left_zero & ~right_zero
    both_zero = left_zero & right_zero
    if np.any(only_left):
        scaled[only_left] = left[only_left]
        exponent[only_left] = left_exp[only_left]
    if np.any(only_right):
        scaled[only_right] = right[only_right]
        exponent[only_right] = right_exp[only_right]
    if np.any(both_zero):
        scaled[both_zero] = 0.0
        exponent[both_zero] = 0
    return scaled, exponent


def multiply_scaled_probability_upper(
    matrices: np.ndarray,
    exponents: np.ndarray,
    probabilities: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mantissas, probability_exponents = np.frexp(probabilities)
    products = mul_up(matrices, mantissas[:, None, None])
    scaled, extra = normalize_scaled_upper(products)
    return (
        scaled,
        exponents + probability_exponents.astype(np.int64) + extra,
    )


def fixed_region_matrices_scaled_upper(
    impulse: np.ndarray,
    combination_weights: list[list[np.ndarray]],
) -> tuple[np.ndarray, np.ndarray]:
    current, current_exp = normalize_scaled_upper(impulse)
    impulse_scaled = current.copy()
    impulse_exp = current_exp.copy()
    for completed, by_next_count in enumerate(combination_weights, start=1):
        previous_bits = completed * STEP_BITS
        total_bits = (completed + 1) * STEP_BITS
        updated = np.zeros((total_bits + 1, 2, 2))
        updated_exp = np.zeros(total_bits + 1, dtype=np.int64)
        for next_count, probabilities in enumerate(by_next_count):
            products, product_exp = scaled_matmul_upper(
                current,
                current_exp,
                impulse_scaled[next_count],
                impulse_exp[next_count],
            )
            terms, term_exp = multiply_scaled_probability_upper(
                products, product_exp, probabilities
            )
            destination = slice(next_count, next_count + previous_bits + 1)
            summed, summed_exp = scaled_add_upper(
                updated[destination],
                updated_exp[destination],
                terms,
                term_exp,
            )
            updated[destination] = summed
            updated_exp[destination] = summed_exp
        current = updated
        current_exp = updated_exp
    return current, current_exp


def scaled_entrywise_max_upper(
    matrices: np.ndarray, exponents: np.ndarray
) -> tuple[np.ndarray, np.int64]:
    target = int(np.max(exponents))
    aligned = ldexp_up(matrices, exponents - target)
    maximum = np.max(aligned, axis=0)
    scaled, extra = normalize_scaled_upper(maximum)
    return scaled, np.int64(target + int(extra))


def fair_body_and_envelope_matrices_scaled_upper(
    fixed_regions: np.ndarray, fixed_exp: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    body = np.zeros((OUTER_BLOCKS, 2, 2))
    body_exp = np.zeros(OUTER_BLOCKS, dtype=np.int64)
    mixed = np.zeros((OUTER_BLOCKS - 1, 2, 2))
    mixed_exp = np.zeros(OUTER_BLOCKS - 1, dtype=np.int64)
    current = fixed_regions
    current_exp = fixed_exp
    for body_count in range(1, OUTER_BLOCKS + 1):
        current, current_exp = scaled_add_upper(
            current[:-1],
            current_exp[:-1],
            current[1:],
            current_exp[1:],
        )
        current_exp -= 1  # exact division by two
        body[body_count - 1] = current[0]
        body_exp[body_count - 1] = current_exp[0]
        if body_count < OUTER_BLOCKS:
            envelope, envelope_exp = scaled_entrywise_max_upper(
                current, current_exp
            )
            mixed[body_count - 1] = envelope
            mixed_exp[body_count - 1] = envelope_exp
    return body, body_exp, mixed, mixed_exp


def fair_body_and_envelope_matrices_upper(
    fixed_regions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return G[b,0] and max_j G[b,j] for all feasible b."""
    body = np.zeros((OUTER_BLOCKS, 2, 2))
    mixed = np.zeros((OUTER_BLOCKS - 1, 2, 2))
    current = fixed_regions
    for body_count in range(1, OUTER_BLOCKS + 1):
        current = mul_up(add_up(current[:-1], current[1:]), 0.5)
        body[body_count - 1] = current[0]
        if body_count < OUTER_BLOCKS:
            mixed[body_count - 1] = np.max(current, axis=0)
    return body, mixed


def normalize_power_of_two_upper(
    matrices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    maximum = np.max(matrices, axis=(1, 2))
    if np.any(maximum <= 0.0):
        raise ArithmeticError("matrix normalization encountered zero")
    _, exponents = np.frexp(maximum)
    scaled = np.ldexp(matrices, -exponents[:, None, None])
    # Scaling by a power of two is exact unless a subnormal is encountered.
    # The outward step also covers that case.
    rounded = up(scaled)
    rounded = np.where(matrices == 0.0, 0.0, rounded)
    return rounded, exponents.astype(np.int64)


def matrix_power_1024_log2_upper(
    matrices: np.ndarray, input_exponents: np.ndarray | None = None
) -> np.ndarray:
    scaled, exponents = normalize_power_of_two_upper(matrices)
    if input_exponents is not None:
        exponents += input_exponents
    for _ in range(10):
        scaled = matmul_up(scaled, scaled)
        exponents *= 2
        scaled, extra = normalize_power_of_two_upper(scaled)
        exponents += extra
    row_sum = add_up(scaled[:, 0, 0], scaled[:, 0, 1])
    result = np.empty(len(row_sum))
    for index, value in enumerate(row_sum):
        result[index] = float(
            up(
                log2_upper(arb_from_float(float(value)))
                + int(exponents[index])
            )
        )
    return result


def all_body_inner_upper(
    activation: list[arb], live_counts: list[int]
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    combination_weights = region_combination_weights_upper()
    best_body = np.full(OUTER_BLOCKS, math.inf)
    best_mixed = np.full(OUTER_BLOCKS - 1, math.inf)
    best_body_tilt = np.full(OUTER_BLOCKS, math.nan)
    best_mixed_tilt = np.full(OUTER_BLOCKS - 1, math.nan)
    for log_surprisal in BODY_LOG_SURPRISALS:
        surprisal = arb(log_surprisal).exp()
        z = (-surprisal).exp()
        impulse = impulse_upper(
            z, activation, live_counts, STEP_BITS
        )
        fixed_regions, fixed_exp = fixed_region_matrices_scaled_upper(
            impulse, combination_weights
        )
        body_matrices, body_exp, mixed_matrices, mixed_exp = (
            fair_body_and_envelope_matrices_scaled_upper(
                fixed_regions, fixed_exp
            )
        )
        body_moments = matrix_power_1024_log2_upper(
            body_matrices, body_exp
        )
        mixed_moments = matrix_power_1024_log2_upper(
            mixed_matrices, mixed_exp
        )
        chernoff = arb_upper_float(
            arb(DISTANCE) * surprisal / arb(2).log()
        )
        body_candidate = np.minimum(0.0, up(body_moments + chernoff))
        mixed_candidate = np.minimum(0.0, up(mixed_moments + chernoff))
        improved = body_candidate < best_body
        best_body[improved] = body_candidate[improved]
        best_body_tilt[improved] = log_surprisal
        improved = mixed_candidate < best_mixed
        best_mixed[improved] = mixed_candidate[improved]
        best_mixed_tilt[improved] = log_surprisal
        print(
            f"body_tilt,{log_surprisal:.2f},"
            f"body1,{best_body[0]:.9f},mixed1,{best_mixed[0]:.9f}",
            flush=True,
        )
    return best_body, best_mixed, {
        "log_surprisals": list(BODY_LOG_SURPRISALS),
        "body_best_tilts": best_body_tilt.tolist(),
        "mixed_best_tilts": best_mixed_tilt.tolist(),
        "region_recurrence": (
            "exact hypergeometric coefficients enclosed by Arb; all matrix "
            "operations rounded upward"
        ),
        "fair_body_recurrence": "G[b+1,j]=(G[b,j]+G[b,j+1])/2",
        "matrix_power": (
            "ten outward-rounded squarings with exact power-of-two scaling"
        ),
    }


def holder_p_values() -> list[float]:
    return list(
        1.0
        + np.exp2(np.linspace(-16.0, 16.0, 257, dtype=np.float64))
    )


def holder_body_moments_log2_upper(
    outer: np.ndarray,
) -> tuple[list[float], list[float]]:
    mass = (1 << OUTER_DIMENSION) - 1
    body_weights = [
        weight
        for weight in range(1, OUTER_BITS + 1)
        if min(weight, OUTER_BITS - weight) > TAIL_ENDPOINT
        and outer[weight] > 0.0
    ]
    references: list[arb] = []
    likelihoods: list[arb] = []
    ambient = arb(2) ** OUTER_BITS
    density = arb(mass) / ambient
    for weight in body_weights:
        shell = math.comb(OUTER_BITS, weight)
        references.append(arb(shell) / ambient)
        likelihoods.append(
            arb_from_float(float(outer[weight])) / (arb(shell) * density)
        )
    p_values = holder_p_values()
    moments: list[float] = []
    for p in p_values:
        exponent = arb_from_float(p)
        moment = arb(0)
        for reference, likelihood in zip(references, likelihoods):
            moment += reference * likelihood**exponent
        moments.append(log2_upper(moment))
    return p_values, moments


def tail_count_upper_arb(outer: np.ndarray) -> arb:
    result = arb(0)
    for weight in range(1, OUTER_BITS + 1):
        if min(weight, OUTER_BITS - weight) <= TAIL_ENDPOINT:
            result += arb_from_float(float(outer[weight]))
    return result


def body_components_upper(
    outer: np.ndarray,
    body_inner: np.ndarray,
    mixed_inner: np.ndarray,
) -> tuple[float, float, dict[str, object]]:
    p_values, moment_logs = holder_body_moments_log2_upper(outer)
    mass_log = (arb((1 << OUTER_DIMENSION) - 1).log() / arb(2).log())
    tail_count = tail_count_upper_arb(outer)
    body_logs: list[float] = []
    mixed_logs: list[float] = []
    body_rows: list[dict[str, float]] = []
    mixed_rows: list[dict[str, float]] = []
    for body_count in range(1, OUTER_BLOCKS + 1):
        best_body = math.inf
        best_mixed = math.inf
        best_body_p = math.nan
        best_mixed_p = math.nan
        for p, moment_log in zip(p_values, moment_logs):
            p_arb = arb_from_float(p)
            inverse = 1 / p_arb
            event = 1 - inverse
            body_correction = (
                arb(body_count) * arb(moment_log) * inverse
                + arb(float(body_inner[body_count - 1])) * event
            )
            body_upper = arb_upper_float(body_correction)
            if body_upper < best_body:
                best_body = body_upper
                best_body_p = p
            if body_count < OUTER_BLOCKS:
                mixed_correction = (
                    arb(body_count) * arb(moment_log) * inverse
                    + arb(float(mixed_inner[body_count - 1])) * event
                )
                mixed_upper = arb_upper_float(mixed_correction)
                if mixed_upper < best_mixed:
                    best_mixed = mixed_upper
                    best_mixed_p = p
        location_log = (
            arb(math.comb(OUTER_BLOCKS, body_count)).log()
            / arb(2).log()
        )
        body_value = (
            location_log + body_count * mass_log + arb(best_body)
        )
        body_log = arb_upper_float(body_value)
        body_logs.append(body_log)
        body_rows.append(
            {
                "body_blocks": body_count,
                "best_holder_p": best_body_p,
                "inner_log2_upper": float(body_inner[body_count - 1]),
                "pointwise_log2_upper": body_log,
                "pointwise_margin_bits": -body_log,
            }
        )
        if body_count < OUTER_BLOCKS:
            available = OUTER_BLOCKS - body_count
            # (1+T)^m-1 <= m T (1+T)^(m-1), avoiding subtraction.
            tail_choices = (
                arb(available)
                * tail_count
                * (1 + tail_count) ** (available - 1)
            )
            mixed_value = (
                location_log
                + body_count * mass_log
                + tail_choices.log() / arb(2).log()
                + arb(best_mixed)
            )
            mixed_log = arb_upper_float(mixed_value)
            mixed_logs.append(mixed_log)
            mixed_rows.append(
                {
                    "body_blocks": body_count,
                    "best_holder_p": best_mixed_p,
                    "inner_log2_upper": float(
                        mixed_inner[body_count - 1]
                    ),
                    "pointwise_log2_upper": mixed_log,
                    "pointwise_margin_bits": -mixed_log,
                }
            )
    body_rows.sort(key=lambda row: row["pointwise_log2_upper"], reverse=True)
    mixed_rows.sort(key=lambda row: row["pointwise_log2_upper"], reverse=True)
    return (
        arb_logsum2_upper(body_logs),
        arb_logsum2_upper(mixed_logs),
        {
            "body_dominant_rows": body_rows[:20],
            "mixed_dominant_rows": mixed_rows[:20],
            "holder_grid_points": len(p_values),
            "tail_choice_bound": "m T (1+T)^(m-1)",
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constituent-spectrum", type=Path, required=True)
    parser.add_argument("--activation", type=Path, required=True)
    parser.add_argument("--live-spectrum", type=Path, required=True)
    parser.add_argument("--nearest-one-tail", type=Path, required=True)
    parser.add_argument("--nearest-two-tail", type=Path, required=True)
    parser.add_argument("--nearest-tail-body", type=Path, required=True)
    parser.add_argument("--nearest-mixed", type=Path, required=True)
    parser.add_argument("--arb-precision", type=int, default=192)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.arb_precision < 128:
        raise ValueError("Arb precision must be at least 128 bits")
    ctx.prec = args.arb_precision
    outer, outer_method = outer_spectrum_upper(args.constituent_spectrum)
    activation = load_activation(args.activation)
    live_counts = load_live_counts(args.live_spectrum)
    one_log, one_rows = one_tail_upper(outer, activation, live_counts)
    print(f"one_tail_margin_bits,{-one_log:.12f}", flush=True)
    two_log, two_dominant = two_tail_upper(outer, activation, live_counts)
    print(f"two_tail_margin_bits,{-two_log:.12f}", flush=True)
    three_log = three_plus_tail_upper(outer)
    print(f"three_plus_tail_margin_bits,{-three_log:.12f}", flush=True)
    tail_log = arb_logsum2_upper([one_log, two_log, three_log])
    body_inner, mixed_inner, inner_method = all_body_inner_upper(
        activation, live_counts
    )
    body_log, mixed_log, body_method = body_components_upper(
        outer, body_inner, mixed_inner
    )
    print(f"body_only_margin_bits,{-body_log:.12f}", flush=True)
    print(f"body_and_tail_margin_bits,{-mixed_log:.12f}", flush=True)
    aggregate_log = arb_logsum2_upper(
        [one_log, two_log, three_log, body_log, mixed_log]
    )

    nearest_one = json.loads(args.nearest_one_tail.read_text(encoding="utf-8"))
    nearest_two = json.loads(args.nearest_two_tail.read_text(encoding="utf-8"))
    nearest_split = json.loads(args.nearest_tail_body.read_text(encoding="utf-8"))
    nearest_mixed = json.loads(args.nearest_mixed.read_text(encoding="utf-8"))
    payload = {
        "schema": "riffle-xbch-outward-certificate-v1",
        "parameters": {
            "outer_bits": OUTER_BITS,
            "outer_dimension": OUTER_DIMENSION,
            "outer_blocks": OUTER_BLOCKS,
            "step_bits": STEP_BITS,
            "state_bits": STATE_BITS,
            "epochs": EPOCHS,
            "output_bits": OUTPUT_BITS,
            "distance": DISTANCE,
            "tail_endpoint_width": TAIL_ENDPOINT,
            "arb_precision_bits": ctx.prec,
        },
        "outer_spectrum_method": outer_method,
        "components": {
            "exactly_one_tail_no_body": {
                "log2_upper": one_log,
                "margin_bits_lower": -one_log,
                "nearest_margin_bits": nearest_one["aggregate_margin_bits"],
                "dominant_rows": one_rows[:20],
            },
            "exactly_two_tail_no_body": {
                "log2_upper": two_log,
                "margin_bits_lower": -two_log,
                "nearest_margin_bits": nearest_two["aggregate_margin_bits"],
                "dominant_pair": two_dominant,
                "tilt": -10.0,
            },
            "three_or_more_tail_no_body": {
                "log2_upper": three_log,
                "margin_bits_lower": -three_log,
                "nearest_margin_bits": nearest_split["tail_only"][
                    "three_plus_count_only_margin_bits"
                ],
                "bound": "C(N,3) T^3 (1+T)^(N-3)",
            },
            "body_no_tail": {
                "log2_upper": body_log,
                "margin_bits_lower": -body_log,
                "nearest_margin_bits": nearest_split[
                    "body_only_margin_bits"
                ],
            },
            "body_and_tail": {
                "log2_upper": mixed_log,
                "margin_bits_lower": -mixed_log,
                "nearest_margin_bits": nearest_mixed[
                    "aggregate_margin_bits"
                ],
            },
        },
        "tail_only_log2_upper": tail_log,
        "tail_only_margin_bits_lower": -tail_log,
        "aggregate_log2_upper": aggregate_log,
        "aggregate_margin_bits_lower": -aggregate_log,
        "body_inner_method": inner_method,
        "body_holder_method": body_method,
        "scope": (
            "Outward-rounded certificate for the five disjoint and exhaustive "
            "tail/body message classes under the recorded SplitState transfer "
            "model."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"tail_only_margin_bits,{-tail_log:.12f}")
    print(f"aggregate_margin_bits,{-aggregate_log:.12f}")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
