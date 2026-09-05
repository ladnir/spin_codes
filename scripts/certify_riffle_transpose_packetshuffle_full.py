#!/usr/bin/env python3
"""Outward-rounded occupation certificate for the packet transpose variant.

The checker reads fixed output tilts from an exploratory receipt. It does not
optimize them. Every nonnegative floating-point addition and multiplication in
the row dynamic program is rounded upward with ``nextafter``. Matrix powers
use exact binary rescaling, and scalar logarithms and the final occupation sum
use interval arithmetic.

The arithmetic argument assumes IEEE-754 binary64 operations with NumPy's
separate multiply and add ufuncs. No fused multiply-add operation is used.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np


DEFAULT_SOURCE = Path(
    "constructions/riffle_transpose_packetshuffle_randomstepconv/"
    "receipts/g4_b1024_sigma15_initial.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_transpose_packetshuffle_randomstepconv/"
    "receipts/g4_b1024_sigma15_full_interval.json"
)
POSITIVE_INFINITY = np.float64(math.inf)


def upward(value):
    value = np.asarray(value, dtype=np.float64)
    return np.nextafter(value, POSITIVE_INFINITY)


def multiply_up(left, right):
    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    product = np.multiply(left_array, right_array)
    rounded = upward(product)
    return np.where((left_array == 0.0) | (right_array == 0.0), 0.0, rounded)


def add_up(left, right):
    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    result = np.add(left_array, right_array)
    rounded = upward(result)
    return np.where((left_array == 0.0) & (right_array == 0.0), 0.0, rounded)


def ratio_up(numerator: int, denominator: int) -> float:
    if numerator == 0:
        return 0.0
    return math.nextafter(numerator / denominator, math.inf)


def matrix_multiply_batch_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty_like(left)
    result[..., 0, 0] = add_up(
        multiply_up(left[..., 0, 0], right[0, 0]),
        multiply_up(left[..., 0, 1], right[1, 0]),
    )
    result[..., 0, 1] = add_up(
        multiply_up(left[..., 0, 0], right[0, 1]),
        multiply_up(left[..., 0, 1], right[1, 1]),
    )
    result[..., 1, 0] = add_up(
        multiply_up(left[..., 1, 0], right[0, 0]),
        multiply_up(left[..., 1, 1], right[1, 0]),
    )
    result[..., 1, 1] = add_up(
        multiply_up(left[..., 1, 0], right[0, 1]),
        multiply_up(left[..., 1, 1], right[1, 1]),
    )
    return result


def transition_upper(packet_bits: int, sigma: int, z_text: str, rank: int):
    z_nearest = float(z_text)
    z_upper = math.nextafter(z_nearest, math.inf)
    half_sum = math.ldexp(float(add_up(1.0, z_upper)), -1)
    output_moment = 1.0
    for _ in range(packet_bits):
        output_moment = float(multiply_up(output_moment, half_sum))

    state_zero = math.ldexp(1.0, -sigma)
    state_live = 1.0 - state_zero
    off = float(multiply_up(state_zero, output_moment))
    live = float(multiply_up(state_live, output_moment))
    zero = np.asarray(((1.0, 0.0), (off, live)), dtype=np.float64)
    if rank == 0:
        return zero

    nonzero = np.asarray(((off, live), (off, live)), dtype=np.float64)
    packet_zero = math.ldexp(1.0, -rank)
    packet_live = 1.0 - packet_zero
    return add_up(
        multiply_up(packet_zero, zero), multiply_up(packet_live, nonzero)
    )


def scaled_normalize(mantissa: np.ndarray, exponent: np.ndarray):
    fractions, shifts = np.frexp(mantissa)
    return fractions, exponent + shifts.astype(np.int64)


def scaled_from_float(values: np.ndarray):
    values = np.asarray(values, dtype=np.float64)
    return scaled_normalize(values, np.zeros(values.shape, dtype=np.int64))


def scaled_multiply(left, right):
    left_mantissa, left_exponent = left
    right_mantissa, right_exponent = right
    mantissa = multiply_up(left_mantissa, right_mantissa)
    exponent = left_exponent + right_exponent
    return scaled_normalize(mantissa, exponent)


def scaled_add(left, right):
    left_mantissa, left_exponent = left
    right_mantissa, right_exponent = right
    left_mantissa, right_mantissa = np.broadcast_arrays(
        left_mantissa, right_mantissa
    )
    left_exponent, right_exponent = np.broadcast_arrays(
        left_exponent, right_exponent
    )
    left_zero = left_mantissa == 0.0
    right_zero = right_mantissa == 0.0
    choose_left = left_exponent >= right_exponent
    large_mantissa = np.where(choose_left, left_mantissa, right_mantissa)
    large_exponent = np.where(choose_left, left_exponent, right_exponent)
    small_mantissa = np.where(choose_left, right_mantissa, left_mantissa)
    small_exponent = np.where(choose_left, right_exponent, left_exponent)
    with np.errstate(under="ignore"):
        aligned_small = np.ldexp(small_mantissa, small_exponent - large_exponent)
    mantissa = add_up(large_mantissa, aligned_small)
    exponent = large_exponent.copy()
    mantissa = np.where(left_zero, right_mantissa, mantissa)
    exponent = np.where(left_zero, right_exponent, exponent)
    mantissa = np.where(right_zero, left_mantissa, mantissa)
    exponent = np.where(right_zero, left_exponent, exponent)
    return scaled_normalize(mantissa, exponent)


def scaled_matrix_multiply_batch(left, right):
    left_mantissa, left_exponent = left
    right_mantissa, right_exponent = right
    shape = np.broadcast_shapes(left_mantissa.shape[:-2], right_mantissa.shape[:-2])
    result_mantissa = np.empty(shape + (2, 2), dtype=np.float64)
    result_exponent = np.empty(shape + (2, 2), dtype=np.int64)
    for row in range(2):
        for column in range(2):
            first = scaled_multiply(
                (left_mantissa[..., row, 0], left_exponent[..., row, 0]),
                (right_mantissa[..., 0, column], right_exponent[..., 0, column]),
            )
            second = scaled_multiply(
                (left_mantissa[..., row, 1], left_exponent[..., row, 1]),
                (right_mantissa[..., 1, column], right_exponent[..., 1, column]),
            )
            value = scaled_add(first, second)
            result_mantissa[..., row, column] = value[0]
            result_exponent[..., row, column] = value[1]
    return result_mantissa, result_exponent


def scaled_weight(values, weight: float):
    weight_scaled = scaled_from_float(np.asarray(weight))
    return scaled_multiply(values, weight_scaled)


def normalized_packed_rows_scaled_upper(
    *, packet_bits: int, sigma: int, positions: int, z_text: str
):
    zero_matrix = scaled_from_float(
        transition_upper(packet_bits, sigma, z_text, 0)
    )
    full_matrix = scaled_from_float(
        transition_upper(packet_bits, sigma, z_text, packet_bits)
    )
    partial_matrices = [
        scaled_from_float(transition_upper(packet_bits, sigma, z_text, rank))
        for rank in range(1, packet_bits)
    ]

    no_mantissa = np.zeros((positions + 1, 2, 2), dtype=np.float64)
    no_exponent = np.zeros((positions + 1, 2, 2), dtype=np.int64)
    no_mantissa[0] = np.eye(2)
    no_partial = scaled_normalize(no_mantissa, no_exponent)
    one_partial = (
        np.zeros((packet_bits - 1, positions + 1, 2, 2), dtype=np.float64),
        np.zeros((packet_bits - 1, positions + 1, 2, 2), dtype=np.int64),
    )

    for length in range(positions):
        denominator = length + 1
        old_zero = (
            no_partial[0][: length + 1],
            no_partial[1][: length + 1],
        )
        zero_products = scaled_matrix_multiply_batch(old_zero, zero_matrix)
        full_products = scaled_matrix_multiply_batch(old_zero, full_matrix)

        updated_zero = (
            np.zeros_like(no_partial[0]), np.zeros_like(no_partial[1])
        )
        zero_weights = scaled_from_float(
            np.asarray(
                [
                    ratio_up(denominator - degree, denominator)
                    for degree in range(length + 1)
                ]
            )[:, None, None]
        )
        full_weights = scaled_from_float(
            np.asarray(
                [
                    ratio_up(degree + 1, denominator)
                    for degree in range(length + 1)
                ]
            )[:, None, None]
        )
        zero_terms = scaled_multiply(zero_products, zero_weights)
        full_terms = scaled_multiply(full_products, full_weights)
        updated_zero[0][: length + 1] = zero_terms[0]
        updated_zero[1][: length + 1] = zero_terms[1]
        combined = scaled_add(
            (
                updated_zero[0][1 : length + 2],
                updated_zero[1][1 : length + 2],
            ),
            full_terms,
        )
        updated_zero[0][1 : length + 2] = combined[0]
        updated_zero[1][1 : length + 2] = combined[1]

        updated_partial = (
            np.zeros_like(one_partial[0]), np.zeros_like(one_partial[1])
        )
        if length:
            old_partial = (
                one_partial[0][:, :length].reshape(-1, 2, 2),
                one_partial[1][:, :length].reshape(-1, 2, 2),
            )
            partial_zero_products = scaled_matrix_multiply_batch(
                old_partial, zero_matrix
            )
            partial_full_products = scaled_matrix_multiply_batch(
                old_partial, full_matrix
            )
            partial_zero_products = (
                partial_zero_products[0].reshape(packet_bits - 1, length, 2, 2),
                partial_zero_products[1].reshape(packet_bits - 1, length, 2, 2),
            )
            partial_full_products = (
                partial_full_products[0].reshape(packet_bits - 1, length, 2, 2),
                partial_full_products[1].reshape(packet_bits - 1, length, 2, 2),
            )
            partial_zero_weights = scaled_from_float(
                np.asarray(
                    [
                        ratio_up(length - degree, denominator)
                        for degree in range(length)
                    ]
                )[None, :, None, None]
            )
            partial_full_weights = scaled_from_float(
                np.asarray(
                    [
                        ratio_up(degree + 1, denominator)
                        for degree in range(length)
                    ]
                )[None, :, None, None]
            )
            zero_terms = scaled_multiply(
                partial_zero_products, partial_zero_weights
            )
            full_terms = scaled_multiply(
                partial_full_products, partial_full_weights
            )
            updated_partial[0][:, :length] = zero_terms[0]
            updated_partial[1][:, :length] = zero_terms[1]
            combined = scaled_add(
                (
                    updated_partial[0][:, 1 : length + 1],
                    updated_partial[1][:, 1 : length + 1],
                ),
                full_terms,
            )
            updated_partial[0][:, 1 : length + 1] = combined[0]
            updated_partial[1][:, 1 : length + 1] = combined[1]

        insertion_weight = ratio_up(1, denominator)
        for rank_index in range(packet_bits - 1):
            inserted = scaled_matrix_multiply_batch(
                old_zero, partial_matrices[rank_index]
            )
            inserted = scaled_weight(inserted, insertion_weight)
            destination = (
                updated_partial[0][rank_index, : length + 1],
                updated_partial[1][rank_index, : length + 1],
            )
            combined = scaled_add(destination, inserted)
            updated_partial[0][rank_index, : length + 1] = combined[0]
            updated_partial[1][rank_index, : length + 1] = combined[1]

        no_partial = updated_zero
        one_partial = updated_partial

    return no_partial, one_partial


def scaled_matrix_power_moment_log2_upper(matrix, exponent: int):
    result = scaled_from_float(np.eye(2, dtype=np.float64))
    power = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = scaled_matrix_multiply_batch(result, power)
        remaining >>= 1
        if remaining:
            power = scaled_matrix_multiply_batch(power, power)
    moment = scaled_add(
        (result[0][0, 0], result[1][0, 0]),
        (result[0][0, 1], result[1][0, 1]),
    )
    mantissa = float(moment[0])
    binary_exponent = int(moment[1])
    numerator, denominator = mantissa.as_integer_ratio()
    mantissa_interval = mp.iv.mpf(numerator) / denominator
    return mp.iv.mpf(binary_exponent) + mp.iv.log(mantissa_interval) / mp.iv.log(2)


def normalized_packed_rows_upper(
    *, packet_bits: int, sigma: int, positions: int, z_text: str
) -> tuple[np.ndarray, np.ndarray]:
    """Return upward bounds on every normalized packed-row matrix."""
    zero_matrix = transition_upper(packet_bits, sigma, z_text, 0)
    full_matrix = transition_upper(packet_bits, sigma, z_text, packet_bits)
    partial_matrices = np.stack(
        [
            transition_upper(packet_bits, sigma, z_text, rank)
            for rank in range(1, packet_bits)
        ]
    )

    no_partial = np.zeros((positions + 1, 2, 2), dtype=np.float64)
    no_partial[0] = np.eye(2)
    one_partial = np.zeros(
        (packet_bits - 1, positions + 1, 2, 2), dtype=np.float64
    )

    for length in range(positions):
        denominator = length + 1
        old_zero = no_partial[: length + 1]
        zero_products = matrix_multiply_batch_up(old_zero, zero_matrix)
        full_products = matrix_multiply_batch_up(old_zero, full_matrix)

        updated_zero = np.zeros_like(no_partial)
        zero_weights = np.asarray(
            [ratio_up(denominator - degree, denominator) for degree in range(length + 1)]
        )[:, None, None]
        full_weights = np.asarray(
            [ratio_up(degree + 1, denominator) for degree in range(length + 1)]
        )[:, None, None]
        updated_zero[: length + 1] = multiply_up(zero_weights, zero_products)
        updated_zero[1 : length + 2] = add_up(
            updated_zero[1 : length + 2],
            multiply_up(full_weights, full_products),
        )

        updated_partial = np.zeros_like(one_partial)
        if length:
            old_partial = one_partial[:, :length].reshape(-1, 2, 2)
            partial_zero_products = matrix_multiply_batch_up(
                old_partial, zero_matrix
            ).reshape(packet_bits - 1, length, 2, 2)
            partial_full_products = matrix_multiply_batch_up(
                old_partial, full_matrix
            ).reshape(packet_bits - 1, length, 2, 2)
            partial_zero_weights = np.asarray(
                [ratio_up(length - degree, denominator) for degree in range(length)]
            )[None, :, None, None]
            partial_full_weights = np.asarray(
                [ratio_up(degree + 1, denominator) for degree in range(length)]
            )[None, :, None, None]
            updated_partial[:, :length] = multiply_up(
                partial_zero_weights, partial_zero_products
            )
            updated_partial[:, 1 : length + 1] = add_up(
                updated_partial[:, 1 : length + 1],
                multiply_up(partial_full_weights, partial_full_products),
            )

        insertion_weight = ratio_up(1, denominator)
        for rank_index in range(packet_bits - 1):
            inserted = matrix_multiply_batch_up(
                old_zero, partial_matrices[rank_index]
            )
            updated_partial[rank_index, : length + 1] = add_up(
                updated_partial[rank_index, : length + 1],
                multiply_up(insertion_weight, inserted),
            )

        no_partial = updated_zero
        one_partial = updated_partial

    return no_partial, one_partial


def matrix_multiply_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return matrix_multiply_batch_up(left[np.newaxis, ...], right)[0]


def normalize_binary(matrix: np.ndarray, exponent: int):
    maximum = float(np.max(matrix))
    if maximum == 0.0:
        return matrix, exponent
    scale_exponent = math.frexp(maximum)[1]
    return np.ldexp(matrix, -scale_exponent), exponent + scale_exponent


def matrix_power_moment_log2_upper(matrix: np.ndarray, exponent: int):
    result = np.eye(2, dtype=np.float64)
    result_exponent = 0
    power, power_exponent = normalize_binary(matrix, 0)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = matrix_multiply_up(result, power)
            result_exponent += power_exponent
            result, result_exponent = normalize_binary(
                result, result_exponent
            )
        remaining >>= 1
        if remaining:
            power = matrix_multiply_up(power, power)
            power_exponent *= 2
            power, power_exponent = normalize_binary(power, power_exponent)
    mantissa = float(add_up(result[0, 0], result[0, 1]))
    numerator, denominator = mantissa.as_integer_ratio()
    mantissa_interval = mp.iv.mpf(numerator) / denominator
    return mp.iv.mpf(result_exponent) + mp.iv.log(mantissa_interval) / mp.iv.log(2)


def interval_upper(value):
    return mp.mpf(value._mpi_[1])


def interval_lower(value):
    return mp.mpf(value._mpi_[0])


def decimal(value, digits: int = 30) -> str:
    return mp.nstr(value, digits)


def fixed_z_text(grid_tenth: int) -> str:
    log_surprisal = mp.mpf(grid_tenth) / 10
    return mp.nstr(mp.exp(-mp.exp(log_surprisal)), 50)


def outer_log2_interval(
    *, blocks: int, outer_dimension: int, outer_bits: int, active_blocks: int
):
    exact_count = (
        math.comb(blocks, active_blocks)
        * ((1 << outer_dimension) - 1) ** active_blocks
    )
    conditioning = (
        mp.iv.mpf(1) - mp.iv.mpf(1) / (1 << outer_bits)
    ) ** (-active_blocks)
    return (
        mp.iv.log(mp.iv.mpf(exact_count)) + mp.iv.log(conditioning)
    ) / mp.iv.log(2)


def self_test() -> dict[str, float]:
    no_partial, one_partial = normalized_packed_rows_upper(
        packet_bits=4,
        sigma=5,
        positions=5,
        z_text="0.61",
    )
    from analyze_riffle_transpose_packetshuffle import (
        normalized_packed_row_logs,
        packed_coefficient_logs,
    )

    exact_no, exact_one = packed_coefficient_logs(
        packet_bits=4, sigma=5, positions=5, z=0.61
    )
    maximum_gap = 0.0
    for active_blocks in range(1, 21):
        exact = np.exp(
            normalized_packed_row_logs(
                exact_no,
                exact_one,
                active_blocks=active_blocks,
                packet_bits=4,
                positions=5,
            )
        )
        full, remainder = divmod(active_blocks, 4)
        upper = no_partial[full] if remainder == 0 else one_partial[remainder - 1, full]
        if np.any(upper < exact):
            raise AssertionError("upward row bound fell below the reference")
        maximum_gap = max(maximum_gap, float(np.max(upper - exact)))
    scaled_no, scaled_one = normalized_packed_rows_scaled_upper(
        packet_bits=4,
        sigma=5,
        positions=5,
        z_text="0.61",
    )
    for active_blocks in range(1, 21):
        full, remainder = divmod(active_blocks, 4)
        scaled = (
            (scaled_no[0][full], scaled_no[1][full])
            if remainder == 0
            else (
                scaled_one[0][remainder - 1, full],
                scaled_one[1][remainder - 1, full],
            )
        )
        recovered = np.ldexp(scaled[0], scaled[1])
        exact = np.exp(
            normalized_packed_row_logs(
                exact_no,
                exact_one,
                active_blocks=active_blocks,
                packet_bits=4,
                positions=5,
            )
        )
        if np.any(recovered < exact):
            raise AssertionError("scaled upward row bound fell below reference")
        maximum_gap = max(maximum_gap, float(np.max(recovered - exact)))
    return {"maximum_small_row_upward_gap": maximum_gap}


def evaluate(source: Path, dps: int, tilt_quantum_tenths: int) -> dict[str, object]:
    mp.mp.dps = dps
    mp.iv.dps = dps
    source_payload = json.loads(source.read_text(encoding="utf-8"))
    case = source_payload["case"]
    message_bits = int(source_payload["message_bits"])
    outer_bits = int(case["outer_bits"])
    outer_dimension = int(case["outer_dimension_bits"])
    blocks = int(case["outer_blocks"])
    packet_bits = int(case["packet_bits"])
    groups = int(case["fixed_block_groups"])
    sigma = int(case["sigma"])
    distance = int(case["distance"])

    tilt_by_occupation: dict[int, int] = {}
    for row in case["occupation_rows"]:
        active_blocks = int(row["active_outer_blocks"])
        source_tenth = int(round(float(row["log_output_surprisal"]) * 10))
        grid_tenth = int(
            round(source_tenth / tilt_quantum_tenths) * tilt_quantum_tenths
        )
        tilt_by_occupation[active_blocks] = grid_tenth

    occupations_by_tilt: dict[int, list[int]] = {}
    for active_blocks, grid_tenth in tilt_by_occupation.items():
        occupations_by_tilt.setdefault(grid_tenth, []).append(active_blocks)

    rows: list[dict[str, object]] = []
    total = mp.iv.mpf(0)
    completed = 0
    for grid_tenth in sorted(occupations_by_tilt):
        z_text = fixed_z_text(grid_tenth)
        no_partial, one_partial = normalized_packed_rows_scaled_upper(
            packet_bits=packet_bits,
            sigma=sigma,
            positions=groups,
            z_text=z_text,
        )
        z_interval = mp.iv.mpf(z_text)
        correction = -distance * mp.iv.log(z_interval) / mp.iv.log(2)
        for active_blocks in occupations_by_tilt[grid_tenth]:
            full, remainder = divmod(active_blocks, packet_bits)
            row_matrix = (
                (no_partial[0][full], no_partial[1][full])
                if remainder == 0
                else (
                    one_partial[0][remainder - 1, full],
                    one_partial[1][remainder - 1, full],
                )
            )
            log_moment = scaled_matrix_power_moment_log2_upper(
                row_matrix, outer_bits
            )
            inner_upper = min(
                mp.mpf(0), interval_upper(log_moment + correction)
            )
            outer_upper = interval_upper(
                outer_log2_interval(
                    blocks=blocks,
                    outer_dimension=outer_dimension,
                    outer_bits=outer_bits,
                    active_blocks=active_blocks,
                )
            )
            point_interval = mp.iv.mpf(outer_upper) + mp.iv.mpf(inner_upper)
            point_upper = interval_upper(point_interval)
            total += mp.iv.power(2, mp.iv.mpf(point_upper))
            rows.append(
                {
                    "active_outer_blocks": active_blocks,
                    "fixed_log_surprisal": decimal(mp.mpf(grid_tenth) / 10, 8),
                    "fixed_z": z_text,
                    "outer_log2_upper": decimal(outer_upper),
                    "inner_log2_upper": decimal(inner_upper),
                    "pointwise_log2_upper": decimal(point_upper),
                }
            )
        completed += len(occupations_by_tilt[grid_tenth])
        print(
            f"tilt,{grid_tenth / 10:.1f},occupations,"
            f"{completed},{blocks}",
            flush=True,
        )

    total_upper = interval_upper(total)
    lambda_interval = -mp.iv.log(mp.iv.mpf(total_upper)) / mp.iv.log(2)
    lambda_lower = interval_lower(lambda_interval)
    rows.sort(key=lambda row: int(row["active_outer_blocks"]))
    dominant = sorted(
        rows, key=lambda row: mp.mpf(str(row["pointwise_log2_upper"])), reverse=True
    )[:10]
    return {
        "schema": "riffle-transpose-packetshuffle-full-interval-v1",
        "candidate": "Riffle TransposePacketShuffle-RandomStepConv g=4",
        "source_tilt_receipt": str(source),
        "arithmetic": {
            "binary_row_dp": "IEEE-754 binary64 mantissa plus integer binary exponent; nextafter after every nonnegative mantissa add and multiply",
            "matrix_power_scaling": "per-entry exact powers of two",
            "scalar_and_sum": f"mpmath interval arithmetic at {dps} decimal digits",
            "optimizer_in_checker": False,
            "tilt_quantum_tenths": tilt_quantum_tenths,
        },
        "parameters": {
            "message_bits": message_bits,
            "output_bits": outer_bits * groups * packet_bits,
            "outer_bits": outer_bits,
            "outer_dimension_bits": outer_dimension,
            "outer_blocks": blocks,
            "packet_bits": packet_bits,
            "fixed_block_groups": groups,
            "sigma": sigma,
            "distance": distance,
            "relative_distance_floor": distance / (outer_bits * groups * packet_bits),
        },
        "packing_domination": "proved in proof/PACKING_DOMINATION.md",
        "self_test": self_test(),
        "occupation_count": len(rows),
        "log2_expected_bad_upper": decimal(-lambda_lower),
        "lambda_bits_lower": decimal(lambda_lower),
        "target_lambda_bits": 40,
        "passes_target": bool(lambda_lower > 40),
        "dominant_certified_rows": dominant,
        "occupation_rows": rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument(
        "--tilt-quantum-tenths",
        type=int,
        default=1,
        help="round source log-surprisal tilts to this many tenths",
    )
    parser.add_argument("--self-test-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mp.mp.dps = args.dps
    mp.iv.dps = args.dps
    if args.self_test_only:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    if args.tilt_quantum_tenths < 1:
        raise ValueError("tilt quantum must be positive")
    payload = evaluate(args.source, args.dps, args.tilt_quantum_tenths)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = {
        "lambda_bits_lower": payload["lambda_bits_lower"],
        "passes_target": payload["passes_target"],
        "dominant_certified_rows": payload["dominant_certified_rows"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"wrote,{args.output}", flush=True)


if __name__ == "__main__":
    main()
