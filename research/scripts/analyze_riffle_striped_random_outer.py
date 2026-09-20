#!/usr/bin/env python3
"""First-moment diagnostic for Riffle S-Stripe RandomStepConv.

The message is divided among independent random rate-half outer injections.
Packets are arranged by outer block and packet coordinate.  The packet
coordinates are divided equally among S consecutive inner regions, and each
region receives an independent uniform permutation of its assigned packets.

The evaluator conditions on the number of active outer blocks.  Its general
path uses entrywise Chernoff bounds for regional candidate counts.  When every
packet coordinate has its own stripe, an optional log-domain dynamic program
computes all regional coefficients exactly.  Both paths use a Chernoff bound
for final output weight.  The resulting values are diagnostics in a
proof-safe direction.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


LOG2 = math.log(2.0)
DEFAULT_MESSAGE_BITS = 1 << 20
DEFAULT_DISTANCE = 0.09
DEFAULT_OUTPUT = Path(
    "constructions/riffle_striped_randomstepconv_g8/receipts/initial.json"
)


def log_two_power_minus_one(exponent: int) -> float:
    return exponent * LOG2 + math.log1p(-math.ldexp(1.0, -exponent))


def log_choose(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def matmul2(
    left: tuple[float, float, float, float],
    right: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    a, b, c, d = left
    e, f, g, h = right
    return (
        a * e + b * g,
        a * f + b * h,
        c * e + d * g,
        c * f + d * h,
    )


def log_matrix_power_moment(
    matrix: tuple[float, float, float, float], exponent: int
) -> float:
    """Return log([1,0] matrix^exponent [1,1]^T) with scaling.

    The explicit two-by-two kernel is an optimizer hot path.  It avoids NumPy
    allocation and dispatch on every objective evaluation.
    """
    power = matrix
    power_scale = 0.0
    vector0 = 1.0
    vector1 = 0.0
    vector_scale = 0.0
    remaining = exponent
    while remaining:
        if remaining & 1:
            a, b, c, d = power
            next0 = vector0 * a + vector1 * c
            next1 = vector0 * b + vector1 * d
            norm = max(next0, next1)
            if not math.isfinite(norm) or norm <= 0.0:
                raise ArithmeticError("matrix moment vector lost positive scale")
            vector0 = next0 / norm
            vector1 = next1 / norm
            vector_scale += power_scale + math.log(norm)
        remaining >>= 1
        if not remaining:
            break
        power = matmul2(power, power)
        norm = max(power)
        if not math.isfinite(norm) or norm <= 0.0:
            raise ArithmeticError("matrix moment power lost positive scale")
        power = tuple(value / norm for value in power)
        power_scale = 2.0 * power_scale + math.log(norm)
    return vector_scale + math.log(vector0 + vector1)


def scaled_matrix_power(
    matrix: tuple[float, float, float, float], exponent: int
) -> tuple[tuple[float, float, float, float], float]:
    """Return a scaled two-by-two matrix power."""
    result = (1.0, 0.0, 0.0, 1.0)
    result_scale = 0.0
    power = matrix
    power_scale = 0.0
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = matmul2(result, power)
            norm = max(result)
            if not math.isfinite(norm) or norm <= 0.0:
                raise ArithmeticError("matrix product lost positive scale")
            result = tuple(value / norm for value in result)
            result_scale += power_scale + math.log(norm)
        remaining >>= 1
        if not remaining:
            break
        power = matmul2(power, power)
        norm = max(power)
        if not math.isfinite(norm) or norm <= 0.0:
            raise ArithmeticError("matrix square lost positive scale")
        power = tuple(value / norm for value in power)
        power_scale = 2.0 * power_scale + math.log(norm)
    return result, result_scale


def log_banded_matrix_moment(
    bands: list[tuple[tuple[float, float, float, float], int]]
) -> float:
    """Return the moment of an ordered product of powered matrices."""
    vector0 = 1.0
    vector1 = 0.0
    scale = 0.0
    for matrix, exponent in bands:
        powered, power_scale = scaled_matrix_power(matrix, exponent)
        a, b, c, d = powered
        next0 = vector0 * a + vector1 * c
        next1 = vector0 * b + vector1 * d
        norm = max(next0, next1)
        if not math.isfinite(norm) or norm <= 0.0:
            raise ArithmeticError("banded moment vector lost positive scale")
        vector0 = next0 / norm
        vector1 = next1 / norm
        scale += power_scale + math.log(norm)
    return scale + math.log(vector0 + vector1)


def contiguous_band_counts(stripes: int, tilt_bands: int) -> list[int]:
    if not 1 <= tilt_bands <= stripes:
        raise ValueError("tilt-band count must lie between one and S")
    quotient, remainder = divmod(stripes, tilt_bands)
    return [quotient + (index < remainder) for index in range(tilt_bands)]


def transition_matrices(
    packet_bits: int, sigma: int, z: float
) -> tuple[
    tuple[float, float, float, float],
    tuple[float, float, float, float],
]:
    """Return the zero-input and uniform-candidate transition moments."""
    state_zero = math.ldexp(1.0, -sigma)
    output_moment = math.ldexp((1.0 + z) ** packet_bits, -packet_bits)
    off = state_zero * output_moment
    live = (1.0 - state_zero) * output_moment

    zero_input = (1.0, 0.0, off, live)
    nonzero_input = (off, live, off, live)
    packet_zero = math.ldexp(1.0, -packet_bits)
    packet_nonzero = 1.0 - packet_zero
    uniform_candidate = tuple(
        packet_zero * zero_value + packet_nonzero * nonzero_value
        for zero_value, nonzero_value in zip(zero_input, nonzero_input)
    )
    return zero_input, uniform_candidate


@dataclass(frozen=True)
class InnerPoint:
    log2_upper: float
    log_candidate_tilt: float | list[float] | None
    log_output_surprisal: float


def optimize_inner_bound(
    *,
    packet_bits: int,
    sigma: int,
    packet_positions: int,
    blocks: int,
    packets_per_block: int,
    stripes: int,
    occupation: int,
    distance: int,
) -> InnerPoint:
    """Bound low output weight at one fixed outer-block occupation."""
    columns_per_stripe = packets_per_block // stripes
    region_positions = blocks * columns_per_stripe
    region_candidates = occupation * columns_per_stripe

    if occupation == blocks:
        def objective_full(parameters: np.ndarray) -> float:
            log_surprisal = float(parameters[0])
            z = math.exp(-math.exp(log_surprisal))
            _zero, candidate = transition_matrices(packet_bits, sigma, z)
            return (
                log_matrix_power_moment(candidate, packet_positions)
                + distance * math.exp(log_surprisal)
            )

        results = [
            minimize(
                objective_full,
                np.asarray([start]),
                method="Nelder-Mead",
                bounds=((-18.0, 3.0),),
                options={"xatol": 2e-8, "fatol": 2e-8, "maxiter": 300},
            )
            for start in (-4.0, -1.0, 1.0)
        ]
        best = min(results, key=lambda result: float(result.fun))
        return InnerPoint(
            log2_upper=min(0.0, float(best.fun) / LOG2),
            log_candidate_tilt=None,
            log_output_surprisal=float(best.x[0]),
        )

    coefficient_normalizer = stripes * log_choose(
        region_positions, region_candidates
    )
    fraction = occupation / blocks
    base_log_x = math.log(fraction / (1.0 - fraction))

    def common_objective(parameters: np.ndarray) -> float:
        log_x = float(parameters[0])
        log_surprisal = float(parameters[1])
        x = math.exp(log_x)
        z = math.exp(-math.exp(log_surprisal))
        zero_input, candidate = transition_matrices(packet_bits, sigma, z)
        marked = tuple(
            zero_value + x * candidate_value
            for zero_value, candidate_value in zip(zero_input, candidate)
        )
        return (
            log_matrix_power_moment(marked, packet_positions)
            - stripes * region_candidates * log_x
            - coefficient_normalizer
            + distance * math.exp(log_surprisal)
        )

    common_results = [
        minimize(
            common_objective,
            np.asarray([base_log_x, output_start]),
            method="L-BFGS-B",
            bounds=((-40.0, 40.0), (-18.0, 3.0)),
            options={"ftol": 1e-11, "gtol": 1e-8, "maxiter": 400},
        )
        for output_start in (-4.0, -1.0, 1.0)
    ]
    common_best = min(common_results, key=lambda result: float(result.fun))

    def objective(parameters: np.ndarray) -> float:
        log_x_values = parameters[:4]
        log_surprisal = float(parameters[-1])
        z = math.exp(-math.exp(log_surprisal))
        zero_input, candidate = transition_matrices(packet_bits, sigma, z)
        entry_logs: list[float] = []
        region_log_choose = coefficient_normalizer / stripes
        for entry, log_x in enumerate(log_x_values):
            x = math.exp(float(log_x))
            marked = tuple(
                zero_value + x * candidate_value
                for zero_value, candidate_value in zip(zero_input, candidate)
            )
            powered, power_scale = scaled_matrix_power(marked, region_positions)
            entry_value = powered[entry]
            entry_logs.append(
                power_scale
                + math.log(max(entry_value, np.finfo(float).tiny))
                - region_candidates * float(log_x)
                - region_log_choose
            )
        if not all(math.isfinite(value) for value in entry_logs):
            return 1e300
        entry_scale = max(entry_logs)
        region_bound = tuple(math.exp(value - entry_scale) for value in entry_logs)
        try:
            region_moment = log_matrix_power_moment(region_bound, stripes)
        except ArithmeticError:
            return 1e300
        return region_moment + stripes * entry_scale + distance * math.exp(
            log_surprisal
        )

    starts: list[np.ndarray] = [
        np.asarray([base_log_x] * 4 + [output_start])
        for output_start in (-4.0, -1.0, 1.0)
    ]
    # The off-to-off coefficient is dominated by candidate packets that
    # realize as zero.  Its natural marker saddle is shifted by g*ln(2).
    zero_shift = packet_bits * LOG2
    starts.append(
        np.asarray(
            [float(common_best.x[0])] * 4 + [float(common_best.x[1])]
        )
    )
    starts.append(
        np.asarray(
            [base_log_x + zero_shift, base_log_x, base_log_x, base_log_x, -1.0]
        )
    )
    bounds = [(-40.0, 40.0)] * 4 + [(-18.0, 3.0)]
    results = [
        minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"ftol": 1e-11, "gtol": 1e-8, "maxiter": 400},
        )
        for start in starts
    ]
    best = min(results, key=lambda result: float(result.fun))
    if float(common_best.fun) < float(best.fun):
        best_value = float(common_best.fun)
        best_x = np.asarray(
            [float(common_best.x[0])] * 4 + [float(common_best.x[1])]
        )
    else:
        best_value = float(best.fun)
        best_x = best.x
    return InnerPoint(
        log2_upper=min(0.0, best_value / LOG2),
        log_candidate_tilt=[float(value) for value in best_x[:-1]],
        log_output_surprisal=float(best_x[-1]),
    )


def exact_region_matrix(
    zero_input: np.ndarray,
    candidate: np.ndarray,
    positions: int,
    candidates: int,
) -> np.ndarray:
    """Exact candidate-position coefficient for a small self-test."""
    coefficients = np.zeros((candidates + 1, 2, 2), dtype=np.float64)
    coefficients[0] = np.eye(2)
    for index in range(positions):
        stop = min(candidates, index + 1)
        updated = np.zeros_like(coefficients)
        for degree in range(stop + 1):
            if degree <= index:
                updated[degree] += coefficients[degree] @ zero_input
            if degree and degree - 1 <= index:
                updated[degree] += coefficients[degree - 1] @ candidate
        coefficients = updated
    return coefficients[candidates] / math.comb(positions, candidates)


def log_matrix_entries(
    matrix: tuple[float, float, float, float]
) -> np.ndarray:
    return np.asarray(
        [math.log(value) if value > 0.0 else -math.inf for value in matrix],
        dtype=np.float64,
    ).reshape(2, 2)


def log_matmul_batch(coefficients: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Multiply a batch of positive two-by-two matrices in the log domain."""
    result = np.empty_like(coefficients)
    result[:, 0, 0] = np.logaddexp(
        coefficients[:, 0, 0] + right[0, 0],
        coefficients[:, 0, 1] + right[1, 0],
    )
    result[:, 0, 1] = np.logaddexp(
        coefficients[:, 0, 0] + right[0, 1],
        coefficients[:, 0, 1] + right[1, 1],
    )
    result[:, 1, 0] = np.logaddexp(
        coefficients[:, 1, 0] + right[0, 0],
        coefficients[:, 1, 1] + right[1, 0],
    )
    result[:, 1, 1] = np.logaddexp(
        coefficients[:, 1, 0] + right[0, 1],
        coefficients[:, 1, 1] + right[1, 1],
    )
    return result


def exact_region_coefficient_logs(
    zero_input: tuple[float, float, float, float],
    candidate: tuple[float, float, float, float],
    positions: int,
) -> np.ndarray:
    """Return log coefficient matrices for every candidate count."""
    log_zero = log_matrix_entries(zero_input)
    log_candidate = log_matrix_entries(candidate)
    coefficients = np.full((positions + 1, 2, 2), -math.inf, dtype=np.float64)
    coefficients[0, 0, 0] = 0.0
    coefficients[0, 1, 1] = 0.0
    for index in range(positions):
        active = coefficients[: index + 1]
        zero_product = log_matmul_batch(active, log_zero)
        candidate_product = log_matmul_batch(active, log_candidate)
        updated = np.full_like(coefficients, -math.inf)
        updated[: index + 1] = zero_product
        updated[1 : index + 2] = np.logaddexp(
            updated[1 : index + 2], candidate_product
        )
        coefficients = updated
    return coefficients


def tilted_transition_profile(
    region: tuple[float, float, float, float], stripes: int
) -> dict[str, float]:
    """Return expected stripe-level state transitions under a tilted moment."""
    matrix = np.asarray(region, dtype=np.float64).reshape(2, 2)
    backward = np.empty((stripes + 1, 2), dtype=np.float64)
    backward[stripes] = 1.0
    backward[stripes] /= np.max(backward[stripes])
    for index in range(stripes - 1, -1, -1):
        backward[index] = matrix @ backward[index + 1]
        backward[index] /= np.max(backward[index])

    distribution = np.asarray([1.0, 0.0], dtype=np.float64)
    counts = np.zeros((2, 2), dtype=np.float64)
    live_probability_sum = 0.0
    for index in range(stripes):
        live_probability_sum += distribution[1]
        beta = backward[index + 1]
        row_weights = matrix * beta[np.newaxis, :]
        row_sums = row_weights.sum(axis=1)
        transitions = row_weights / row_sums[:, np.newaxis]
        edge_mass = distribution[:, np.newaxis] * transitions
        counts += edge_mass
        distribution = edge_mass.sum(axis=0)
    return {
        "off_to_off": float(counts[0, 0]),
        "off_to_live": float(counts[0, 1]),
        "live_to_off": float(counts[1, 0]),
        "live_to_live": float(counts[1, 1]),
        "expected_live_stripe_starts": live_probability_sum,
        "terminal_live_probability": float(distribution[1]),
    }


def evaluate_fullstripe_exact_grid(
    *,
    message_bits: int,
    outer_bits: int,
    packet_bits: int,
    sigma: int,
    relative_distance: float,
    grid_min: float,
    grid_max: float,
    grid_step: float,
) -> dict[str, object]:
    """Evaluate full striping with exact regional candidate coefficients."""
    outer_dimension = outer_bits // 2
    blocks = message_bits // outer_dimension
    packets_per_block = outer_bits // packet_bits
    stripes = packets_per_block
    packet_positions = blocks * packets_per_block
    distance = math.floor(relative_distance * packet_bits * packet_positions)
    log_combinations = np.asarray(
        [log_choose(blocks, occupation) for occupation in range(blocks + 1)]
    )
    best_inner = np.full(blocks + 1, math.inf, dtype=np.float64)
    best_surprisal = np.full(blocks + 1, math.nan, dtype=np.float64)

    grid_count = int(math.floor((grid_max - grid_min) / grid_step + 0.5)) + 1
    grid = [grid_min + index * grid_step for index in range(grid_count)]
    for grid_index, log_surprisal in enumerate(grid):
        z = math.exp(-math.exp(log_surprisal))
        zero_input, candidate = transition_matrices(packet_bits, sigma, z)
        coefficient_logs = exact_region_coefficient_logs(
            zero_input, candidate, blocks
        )
        output_correction = distance * math.exp(log_surprisal)
        for occupation in range(1, blocks + 1):
            entry_logs = coefficient_logs[occupation].reshape(4) - log_combinations[
                occupation
            ]
            entry_scale = float(np.max(entry_logs))
            region = tuple(math.exp(float(value - entry_scale)) for value in entry_logs)
            moment = (
                log_matrix_power_moment(region, stripes)
                + stripes * entry_scale
                + output_correction
            )
            if moment < best_inner[occupation]:
                best_inner[occupation] = moment
                best_surprisal[occupation] = log_surprisal
        if grid_index % max(1, grid_count // 10) == 0:
            print(
                f"grid,{grid_index + 1},{grid_count},"
                f"log_surprisal,{log_surprisal:.6f}",
                flush=True,
            )

    log_messages_per_active_block = log_two_power_minus_one(outer_dimension)
    conditioning_penalty = -math.log1p(-math.ldexp(1.0, -outer_bits))
    rows: list[dict[str, object]] = []
    for occupation in range(1, blocks + 1):
        outer_log = log_combinations[occupation] + occupation * (
            log_messages_per_active_block + conditioning_penalty
        )
        inner_log2 = min(0.0, float(best_inner[occupation]) / LOG2)
        rows.append(
            {
                "active_outer_blocks": occupation,
                "candidate_packet_fraction": occupation / blocks,
                "outer_message_log2": outer_log / LOG2,
                "inner_log2_upper": inner_log2,
                "pointwise_log2_upper": outer_log / LOG2 + inner_log2,
                "log_output_surprisal": float(best_surprisal[occupation]),
            }
        )
    total_log2 = float(
        logsumexp(np.asarray([row["pointwise_log2_upper"] for row in rows]) * LOG2)
        / LOG2
    )
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    dominant_occupation = int(dominant["active_outer_blocks"])
    dominant_surprisal = float(dominant["log_output_surprisal"])
    dominant_z = math.exp(-math.exp(dominant_surprisal))
    dominant_zero, dominant_candidate = transition_matrices(
        packet_bits, sigma, dominant_z
    )
    dominant_coefficients = exact_region_coefficient_logs(
        dominant_zero, dominant_candidate, blocks
    )
    dominant_entry_logs = (
        dominant_coefficients[dominant_occupation].reshape(4)
        - log_combinations[dominant_occupation]
    )
    dominant_scale = float(np.max(dominant_entry_logs))
    dominant_region = tuple(
        math.exp(float(value - dominant_scale)) for value in dominant_entry_logs
    )
    dominant["tilted_stripe_transition_profile"] = tilted_transition_profile(
        dominant_region, stripes
    )
    return {
        "outer_bits": outer_bits,
        "outer_dimension_bits": outer_dimension,
        "outer_blocks": blocks,
        "packets_per_outer_block": packets_per_block,
        "packet_positions": packet_positions,
        "stripe_count": stripes,
        "columns_per_stripe": 1,
        "positions_per_stripe": blocks,
        "sigma": sigma,
        "distance": distance,
        "evaluated_occupation_range": [1, blocks],
        "omitted_occupation_range": None,
        "coefficient_method": "exact_fullstripe_log_dp_output_grid",
        "output_grid": {
            "minimum": grid_min,
            "maximum": grid_max,
            "step": grid_step,
            "count": grid_count,
        },
        "log2_expected_bad_upper": total_log2,
        "lambda_bits_lower": -total_log2,
        "dominant_profile": dominant,
        "occupation_rows": rows,
    }


def fixed_z_coefficient_bound(
    *,
    packet_bits: int,
    sigma: int,
    blocks: int,
    packets_per_block: int,
    stripes: int,
    occupation: int,
    z: float,
) -> float:
    columns_per_stripe = packets_per_block // stripes
    region_positions = blocks * columns_per_stripe
    region_candidates = occupation * columns_per_stripe
    packet_positions = blocks * packets_per_block
    normalizer = stripes * log_choose(region_positions, region_candidates)
    base = math.log(occupation / (blocks - occupation))
    zero_input, candidate = transition_matrices(packet_bits, sigma, z)

    def objective(parameters: np.ndarray) -> float:
        entry_logs: list[float] = []
        for entry, log_x in enumerate(parameters):
            x = math.exp(float(log_x))
            marked = tuple(
                zero_value + x * candidate_value
                for zero_value, candidate_value in zip(zero_input, candidate)
            )
            powered, power_scale = scaled_matrix_power(marked, region_positions)
            entry_logs.append(
                power_scale
                + math.log(max(powered[entry], np.finfo(float).tiny))
                - region_candidates * float(log_x)
                - normalizer / stripes
            )
        entry_scale = max(entry_logs)
        region_bound = tuple(math.exp(value - entry_scale) for value in entry_logs)
        return log_matrix_power_moment(region_bound, stripes) + stripes * entry_scale

    starts = (
        np.asarray([base] * 4),
        np.asarray([base + packet_bits * LOG2, base, base, base]),
    )
    results = [
        minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=[(-40.0, 40.0)] * 4,
            options={"ftol": 1e-13, "gtol": 1e-10, "maxiter": 500},
        )
        for start in starts
    ]
    return math.exp(min(float(result.fun) for result in results))


def self_tests() -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for packet_bits, sigma, blocks, packets_per_block, stripes, occupation, z in (
        (2, 2, 3, 4, 2, 1, 0.70),
        (1, 3, 4, 4, 4, 2, 0.55),
    ):
        columns_per_stripe = packets_per_block // stripes
        region_positions = blocks * columns_per_stripe
        region_candidates = occupation * columns_per_stripe
        zero_tuple, candidate_tuple = transition_matrices(packet_bits, sigma, z)
        zero_input = np.asarray(zero_tuple).reshape(2, 2)
        candidate = np.asarray(candidate_tuple).reshape(2, 2)
        region = exact_region_matrix(
            zero_input, candidate, region_positions, region_candidates
        )
        coefficient_logs = exact_region_coefficient_logs(
            zero_tuple, candidate_tuple, region_positions
        )
        region_from_logs = np.exp(
            coefficient_logs[region_candidates] - log_choose(
                region_positions, region_candidates
            )
        )
        coefficient_error = float(np.max(np.abs(region - region_from_logs)))
        if coefficient_error > 1e-11:
            raise AssertionError("log-domain regional coefficient changed value")
        total = np.linalg.matrix_power(region, stripes)
        exact = float(total[0].sum())
        bound = fixed_z_coefficient_bound(
            packet_bits=packet_bits,
            sigma=sigma,
            blocks=blocks,
            packets_per_block=packets_per_block,
            stripes=stripes,
            occupation=occupation,
            z=z,
        )
        if bound + 1e-11 < exact:
            raise AssertionError("striped coefficient bound fell below exact moment")
        rows.append(
            {
                "packet_bits": packet_bits,
                "sigma": sigma,
                "blocks": blocks,
                "packets_per_block": packets_per_block,
                "stripes": stripes,
                "occupation": occupation,
                "z": z,
                "exact_moment": exact,
                "chernoff_upper": bound,
                "log_coefficient_max_error": coefficient_error,
            }
        )
    return rows


def occupation_schedule(blocks: int, maximum: int | None) -> list[int]:
    stop = blocks if maximum is None else min(blocks, maximum)
    return list(range(1, stop + 1))


def evaluate_case(
    *,
    message_bits: int,
    outer_bits: int,
    packet_bits: int,
    sigma: int,
    stripes: int,
    relative_distance: float,
    max_occupation: int | None,
) -> dict[str, object]:
    if outer_bits <= 0 or outer_bits % 2:
        raise ValueError("outer length must be positive and even")
    if outer_bits % packet_bits:
        raise ValueError("packet size must divide outer length")
    outer_dimension = outer_bits // 2
    if message_bits % outer_dimension:
        raise ValueError("outer dimension must divide message length")
    packets_per_block = outer_bits // packet_bits
    if stripes <= 0 or packets_per_block % stripes:
        raise ValueError("stripe count must divide packets per outer block")

    blocks = message_bits // outer_dimension
    packet_positions = blocks * packets_per_block
    distance = math.floor(relative_distance * packet_bits * packet_positions)
    log_messages_per_active_block = log_two_power_minus_one(outer_dimension)
    # A fixed nonzero input to a uniform random injection has a uniform
    # nonzero output.  Dropping the nonzero condition and dividing by its
    # probability gives an upper bound with independent uniform output bits.
    nonzero_condition_penalty = -math.log1p(-math.ldexp(1.0, -outer_bits))

    rows: list[dict[str, object]] = []
    for occupation in occupation_schedule(blocks, max_occupation):
        inner = optimize_inner_bound(
            packet_bits=packet_bits,
            sigma=sigma,
            packet_positions=packet_positions,
            blocks=blocks,
            packets_per_block=packets_per_block,
            stripes=stripes,
            occupation=occupation,
            distance=distance,
        )
        outer_log = (
            log_choose(blocks, occupation)
            + occupation
            * (log_messages_per_active_block + nonzero_condition_penalty)
        )
        contribution = outer_log / LOG2 + inner.log2_upper
        rows.append(
            {
                "active_outer_blocks": occupation,
                "candidate_packet_fraction": occupation / blocks,
                "outer_message_log2": outer_log / LOG2,
                "inner_log2_upper": inner.log2_upper,
                "pointwise_log2_upper": contribution,
                "log_candidate_tilt": inner.log_candidate_tilt,
                "log_output_surprisal": inner.log_output_surprisal,
            }
        )

    total_log2 = float(
        logsumexp(np.asarray([row["pointwise_log2_upper"] for row in rows]) * LOG2)
        / LOG2
    )
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    omitted = None
    if max_occupation is not None and max_occupation < blocks:
        omitted = [max_occupation + 1, blocks]
    return {
        "outer_bits": outer_bits,
        "outer_dimension_bits": outer_dimension,
        "outer_blocks": blocks,
        "packets_per_outer_block": packets_per_block,
        "packet_positions": packet_positions,
        "stripe_count": stripes,
        "coefficient_method": "entrywise_region_chernoff",
        "columns_per_stripe": packets_per_block // stripes,
        "positions_per_stripe": packet_positions // stripes,
        "sigma": sigma,
        "distance": distance,
        "evaluated_occupation_range": [1, rows[-1]["active_outer_blocks"]],
        "omitted_occupation_range": omitted,
        "log2_expected_bad_upper": total_log2,
        "lambda_bits_lower": -total_log2,
        "dominant_profile": dominant,
        "occupation_rows": rows,
    }


def parse_stripes(values: list[str], packets_per_block: int) -> list[int]:
    result: list[int] = []
    for value in values:
        stripes = packets_per_block if value.lower() == "full" else int(value)
        if stripes not in result:
            result.append(stripes)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--g", type=int, choices=(1, 2, 4, 8), default=8)
    parser.add_argument("--outer-bits", type=int, default=4096)
    parser.add_argument("--sigmas", nargs="+", type=int, default=[40])
    parser.add_argument("--stripes", nargs="+", default=["1", "2", "4", "full"])
    parser.add_argument("--message-bits", type=int, default=DEFAULT_MESSAGE_BITS)
    parser.add_argument("--relative-distance", type=float, default=DEFAULT_DISTANCE)
    parser.add_argument("--max-occupation", type=int)
    parser.add_argument(
        "--fullstripe-exact-grid",
        action="store_true",
        help="use exact regional coefficients when S=B/g",
    )
    parser.add_argument("--grid-min", type=float, default=-4.0)
    parser.add_argument("--grid-max", type=float, default=2.0)
    parser.add_argument("--grid-step", type=float, default=0.05)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    checks = self_tests()
    if args.self_test_only:
        print(json.dumps({"striped_moment_checks": checks}, indent=2, sort_keys=True))
        return

    packets_per_block = args.outer_bits // args.g
    stripe_counts = parse_stripes(args.stripes, packets_per_block)
    cases: list[dict[str, object]] = []
    for stripes in stripe_counts:
        for sigma in args.sigmas:
            print(
                f"case,g,{args.g},B,{args.outer_bits},S,{stripes},sigma,{sigma}",
                flush=True,
            )
            if args.fullstripe_exact_grid and stripes == packets_per_block:
                if args.max_occupation is not None:
                    raise ValueError(
                        "exact full-stripe grid always covers every occupation"
                    )
                case = evaluate_fullstripe_exact_grid(
                    message_bits=args.message_bits,
                    outer_bits=args.outer_bits,
                    packet_bits=args.g,
                    sigma=sigma,
                    relative_distance=args.relative_distance,
                    grid_min=args.grid_min,
                    grid_max=args.grid_max,
                    grid_step=args.grid_step,
                )
            else:
                case = evaluate_case(
                    message_bits=args.message_bits,
                    outer_bits=args.outer_bits,
                    packet_bits=args.g,
                    sigma=sigma,
                    stripes=stripes,
                    relative_distance=args.relative_distance,
                    max_occupation=args.max_occupation,
                )
            cases.append(case)
            dominant = case["dominant_profile"]
            print(
                f"result,S,{stripes},sigma,{sigma},"
                f"lambda_lower,{case['lambda_bits_lower']:.6f},"
                f"occupation,{dominant['active_outer_blocks']},"
                f"pointwise,{dominant['pointwise_log2_upper']:.6f}",
                flush=True,
            )

    payload = {
        "schema": "riffle-striped-random-outer-diagnostic-v1",
        "model": (
            "independent random rate-half outer injections; packet columns "
            "split equally among S consecutive regions; independent uniform "
            "packet permutation inside each region; one-lap RandomStepConv"
        ),
        "message_bits": args.message_bits,
        "packet_bits": args.g,
        "relative_distance": args.relative_distance,
        "target_lambda_bits": 40,
        "small_exact_checks": checks,
        "cases": cases,
        "scope": (
            "The occupation sum is exact over its displayed range. Each case "
            "states whether candidate-position coefficients are exact or use "
            "entrywise Chernoff bounds. "
            "Uniform-nonzero outer outputs are relaxed to independent uniform "
            "bits with the exact conditioning penalty. A complete occupation "
            "range therefore gives a first-moment upper bound, subject to "
            "floating-point validation and outward rounding."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"receipt,{args.output}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
