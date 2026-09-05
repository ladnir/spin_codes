#!/usr/bin/env python3
"""Hybrid exact/saddle diagnostic for the binary biregular EC bridge.

The uniform-slice computation can underflow in a few very small regional
matrix entries at moderate support.  This experiment retains every nonzero
floating-point exact entry and fills only zero entries with an entrywise
positive-coefficient bound.  It is a parameter-selection diagnostic, not an
interval certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from binary_biregular_activation_diagnostic import (
    activation_group_matrix,
    log_matrix_power_entries,
)
from binary_biregular_diagnostic import (
    biregular_parity_shell_counts,
    degree_three_parity_shell_counts,
)
from expander_bounds import log_binom, log_weight_mgf
from regular_ec_diagnostic import wrapping_input_matrices


ScaledMatrix = tuple[np.ndarray, float]


def normalize_scaled(matrix: np.ndarray, log_scale: float = 0.0) -> ScaledMatrix:
    """Normalize a nonnegative matrix and carry its logarithmic scale."""
    maximum = float(np.max(matrix))
    if maximum == 0.0:
        return matrix, -math.inf
    return matrix / maximum, log_scale + math.log(maximum)


def combine_scaled_uniform_slices(
    left: tuple[int, list[ScaledMatrix]],
    right: tuple[int, list[ScaledMatrix]],
    max_weight: int,
) -> tuple[int, list[ScaledMatrix]]:
    """Concatenate normalized conditional slice matrices."""
    left_length, left_slices = left
    right_length, right_slices = right
    total_length = left_length + right_length
    size = left_slices[0][0].shape[0]
    result: list[ScaledMatrix] = []
    for weight in range(min(max_weight, total_length) + 1):
        accumulator = np.zeros((size, size), dtype=np.float64)
        accumulator_scale = -math.inf
        lo = max(0, weight - right_length)
        hi = min(weight, left_length, len(left_slices) - 1)
        for left_weight in range(lo, hi + 1):
            right_weight = weight - left_weight
            if right_weight >= len(right_slices):
                continue
            left_matrix, left_scale = left_slices[left_weight]
            right_matrix, right_scale = right_slices[right_weight]
            if not math.isfinite(left_scale) or not math.isfinite(right_scale):
                continue
            product, product_scale = normalize_scaled(left_matrix @ right_matrix)
            term_scale = (
                left_scale
                + right_scale
                + product_scale
                + log_binom(left_length, left_weight)
                + log_binom(right_length, right_weight)
                - log_binom(total_length, weight)
            )
            if not math.isfinite(accumulator_scale):
                accumulator = product
                accumulator_scale = term_scale
            elif term_scale > accumulator_scale:
                accumulator *= math.exp(accumulator_scale - term_scale)
                accumulator += product
                accumulator_scale = term_scale
            else:
                accumulator += math.exp(term_scale - accumulator_scale) * product
        result.append(normalize_scaled(accumulator, accumulator_scale))
    return total_length, result


def scaled_uniform_slice_transfers(
    *, length: int, max_weight: int, output_marker: float, memory: int,
) -> list[ScaledMatrix]:
    """Compute conditional slice transfers with one scale per slice."""
    zero, one = wrapping_input_matrices(output_marker, memory)
    identity = np.eye(memory + 1, dtype=np.float64)
    result: tuple[int, list[ScaledMatrix]] = (0, [(identity, 0.0)])
    power: tuple[int, list[ScaledMatrix]] = (
        1,
        [normalize_scaled(zero), normalize_scaled(one)],
    )
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_scaled_uniform_slices(result, power, max_weight)
        remaining >>= 1
        if remaining:
            power = combine_scaled_uniform_slices(power, power, max_weight)
    return result[1]


def scaled_exact_region_logs(
    *, k: int, right_degree: int, memory: int, message_weight: int,
    output_marker: float,
) -> np.ndarray:
    """Return entrywise logs of the exact regional uniform-slice matrix."""
    slices = scaled_uniform_slice_transfers(
        length=k // right_degree,
        max_weight=message_weight,
        output_marker=output_marker,
        memory=memory,
    )
    return scaled_exact_region_logs_from_slices(
        k=k,
        right_degree=right_degree,
        message_weight=message_weight,
        slices=slices,
        memory=memory,
    )


def scaled_exact_region_logs_from_slices(
    *, k: int, right_degree: int, message_weight: int,
    slices: list[ScaledMatrix], memory: int,
) -> np.ndarray:
    """Mix one parity shell into an existing scaled slice family."""
    if right_degree == 3:
        denominator, shell_counts = degree_three_parity_shell_counts(
            left_vertices=k,
            support_size=message_weight,
        )
    else:
        denominator, shell_counts = biregular_parity_shell_counts(
            left_vertices=k,
            right_degree=right_degree,
            support_size=message_weight,
        )
    accumulator = np.zeros((memory + 1, memory + 1), dtype=np.float64)
    accumulator_scale = -math.inf
    log_denominator = math.log(denominator)
    for weight, count in shell_counts.items():
        matrix, matrix_scale = slices[weight]
        term_scale = matrix_scale + math.log(count) - log_denominator
        if not math.isfinite(accumulator_scale):
            accumulator = matrix.copy()
            accumulator_scale = term_scale
        elif term_scale > accumulator_scale:
            accumulator *= math.exp(accumulator_scale - term_scale)
            accumulator += matrix
            accumulator_scale = term_scale
        else:
            accumulator += math.exp(term_scale - accumulator_scale) * matrix
    with np.errstate(divide="ignore"):
        return np.log(accumulator) + accumulator_scale


def scaled_exact_fixed_marker_block(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, support_start: int, support_limit: int,
    output_marker: float,
) -> tuple[float, int, float, tuple[float, ...]]:
    """Evaluate a support block using one scaled slice family and marker."""
    slices = scaled_uniform_slice_transfers(
        length=k // right_degree,
        max_weight=support_limit,
        output_marker=output_marker,
        memory=memory,
    )
    terms: list[float] = []
    for support in range(support_start, support_limit + 1):
        regional_logs = scaled_exact_region_logs_from_slices(
            k=k,
            right_degree=right_degree,
            message_weight=support,
            slices=slices,
            memory=memory,
        )
        scale = float(np.max(regional_logs))
        regional = np.exp(regional_logs - scale)
        log_term = (
            log_binom(k, support)
            + left_degree * scale
            + log_weight_mgf(regional, left_degree, memory)
            - cutoff * math.log(output_marker)
        ) / math.log(2.0)
        terms.append(log_term)
    worst_index = int(np.argmax(terms))
    worst = terms[worst_index]
    summed = worst + math.log2(sum(2.0 ** (term - worst) for term in terms))
    return summed, support_start + worst_index, worst, tuple(terms)


def entrywise_saddle_logs(
    *, k: int, right_degree: int, memory: int, message_weight: int,
    output_marker: float, log_input_min: float, log_input_max: float,
    input_steps: int,
) -> np.ndarray:
    """Return gridded entrywise coefficient bounds for one region."""
    region_length = k // right_degree
    best = np.full((memory + 1, memory + 1), math.inf, dtype=np.float64)
    for log_x in np.linspace(log_input_min, log_input_max, input_steps):
        matrix = activation_group_matrix(
            right_degree=right_degree,
            memory=memory,
            input_marker=math.exp(float(log_x)),
            output_marker=output_marker,
            activation_marker=1.0,
        )
        powered = log_matrix_power_entries(matrix, region_length)
        best = np.minimum(best, powered - message_weight * log_x)
    if np.any(np.isposinf(best)):
        raise ArithmeticError("input-marker grid did not bound every entry")
    return best - log_binom(k, message_weight)


def hybrid_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
    log_input_min: float, log_input_max: float, input_steps: int,
) -> tuple[float, int, int, float]:
    """Return the hybrid first-moment term and regional diagnostics."""
    exact_logs = scaled_exact_region_logs(
        k=k,
        right_degree=right_degree,
        message_weight=message_weight,
        output_marker=output_marker,
        memory=memory,
    )
    exact_nonzero = np.isfinite(exact_logs)
    saddle_logs = entrywise_saddle_logs(
        k=k,
        right_degree=right_degree,
        memory=memory,
        message_weight=message_weight,
        output_marker=output_marker,
        log_input_min=log_input_min,
        log_input_max=log_input_max,
        input_steps=input_steps,
    )
    hybrid_logs = saddle_logs.copy()
    hybrid_logs[exact_nonzero] = exact_logs[exact_nonzero]
    scale = float(np.max(hybrid_logs))
    hybrid = np.exp(hybrid_logs - scale)
    log_term = (
        log_binom(k, message_weight)
        + left_degree * scale
        + log_weight_mgf(hybrid, left_degree, memory)
        - cutoff * math.log(output_marker)
    ) / math.log(2.0)
    smallest_log2 = (
        float(np.min(exact_logs[exact_nonzero])) / math.log(2.0)
        if np.any(exact_nonzero)
        else -math.inf
    )
    return (
        log_term,
        int(np.count_nonzero(exact_nonzero)),
        int(np.count_nonzero(~exact_nonzero)),
        smallest_log2,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=1_048_575)
    parser.add_argument("--left-degree", type=int, default=6)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--cutoff", type=int, default=230_729)
    parser.add_argument("--memory", type=int, default=80)
    parser.add_argument("--message-weight", type=int, default=256)
    parser.add_argument("--support-start", type=int)
    parser.add_argument("--support-limit", type=int)
    parser.add_argument("--output-marker", type=float, default=0.99339)
    parser.add_argument("--log-input-min", type=float, default=-12.0)
    parser.add_argument("--log-input-max", type=float, default=-5.0)
    parser.add_argument("--input-steps", type=int, default=113)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if (args.support_start is None) != (args.support_limit is None):
        raise ValueError("support-start and support-limit must be used together")
    if args.support_start is not None:
        summed, worst_support, worst, _ = scaled_exact_fixed_marker_block(
            k=args.k,
            left_degree=args.left_degree,
            right_degree=args.right_degree,
            cutoff=args.cutoff,
            memory=args.memory,
            support_start=args.support_start,
            support_limit=args.support_limit,
            output_marker=args.output_marker,
        )
        print(f"support_block={args.support_start}..{args.support_limit}")
        print(f"worst_support={worst_support}")
        print(f"worst_log2_term={worst:.12f}")
        print(f"summed_log2_bound={summed:.12f}")
        return
    result = hybrid_logterm(
        k=args.k,
        left_degree=args.left_degree,
        right_degree=args.right_degree,
        cutoff=args.cutoff,
        memory=args.memory,
        message_weight=args.message_weight,
        output_marker=args.output_marker,
        log_input_min=args.log_input_min,
        log_input_max=args.log_input_max,
        input_steps=args.input_steps,
    )
    print(f"hybrid_log2_term={result[0]:.12f}")
    print(f"nonzero_exact_entries={result[1]}")
    print(f"zero_exact_entries={result[2]}")
    print(f"smallest_nonzero_exact_log2={result[3]:.12f}")


if __name__ == "__main__":
    main()
