#!/usr/bin/env python3
"""Activation-refined coefficient diagnostic for binary biregular EC.

This experiment adds a marker to transitions that leave the all-zero
convolution state.  Separating a few activation counts can substantially
tighten the positive-coefficient saddle for sparse supports.  The calculation
uses a finite marker grid and floating-point matrix powers, so it is a
diagnostic rather than an interval certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.special import logsumexp

from binary_biregular_diagnostic import (
    _exact_region_matrix,
    biregular_parity_shell_distribution,
)
from expander_bounds import log_binom, log_weight_mgf, normalized


def activation_group_matrix(
    *, right_degree: int, memory: int, input_marker: float,
    output_marker: float, activation_marker: float,
) -> np.ndarray:
    """Return one right-group transfer with input and activation markers."""
    even = sum(
        math.comb(right_degree, j) * input_marker**j
        for j in range(0, right_degree + 1, 2)
    )
    odd = sum(
        math.comb(right_degree, j) * input_marker**j
        for j in range(1, right_degree + 1, 2)
    )
    matrix = np.zeros((memory + 1, memory + 1), dtype=np.float64)
    mixed = 0.5 * (even + odd)
    for state in range(memory - 1):
        matrix[state, 0] = mixed * output_marker
        matrix[state, state + 1] = mixed
    matrix[memory - 1, 0] = even * output_marker
    matrix[memory - 1, memory] = odd
    matrix[memory, 0] = odd * output_marker * activation_marker
    matrix[memory, memory] = even
    return matrix


def log_matrix_power_entries(matrix: np.ndarray, exponent: int) -> np.ndarray:
    """Return entrywise logs of a nonnegative matrix power."""
    result = np.eye(matrix.shape[0], dtype=np.float64)
    result_log_scale = 0.0
    power, power_log_scale = normalized(matrix)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result, local_scale = normalized(result @ power)
            result_log_scale += power_log_scale + local_scale
        remaining >>= 1
        if remaining:
            power, local_scale = normalized(power @ power)
            power_log_scale = 2.0 * power_log_scale + local_scale
    with np.errstate(divide="ignore"):
        return np.log(result) + result_log_scale


def activation_split_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
    activation_cutoff: int, log_input_min: float, log_input_max: float,
    input_steps: int, log_activation_limit: float, activation_steps: int,
) -> tuple[float, int, np.ndarray]:
    """Return a gridded activation-split first-moment bound and unknown count."""
    if k % right_degree:
        raise ValueError("right_degree must divide k")
    region_length = k // right_degree
    size = memory + 1
    exact_counts = np.full(
        (activation_cutoff, size, size), math.inf, dtype=np.float64
    )
    tail = np.full((size, size), math.inf, dtype=np.float64)

    input_grid = np.linspace(log_input_min, log_input_max, input_steps)
    activation_grid = np.linspace(
        -log_activation_limit, log_activation_limit, activation_steps
    )
    for log_x in input_grid:
        x = math.exp(float(log_x))
        for log_v in activation_grid:
            v = math.exp(float(log_v))
            matrix = activation_group_matrix(
                right_degree=right_degree,
                memory=memory,
                input_marker=x,
                output_marker=output_marker,
                activation_marker=v,
            )
            powered = log_matrix_power_entries(matrix, region_length)
            base = powered - message_weight * log_x
            for activations in range(activation_cutoff):
                exact_counts[activations] = np.minimum(
                    exact_counts[activations], base - activations * log_v
                )
            if log_v >= 0.0:
                tail = np.minimum(
                    tail, base - activation_cutoff * log_v
                )

    # Starting in the inactive state, zero activations can only end inactive.
    # These are mathematical zeros, not floating-point failures.
    exact_counts[0, memory, :memory] = -math.inf

    unknown = np.any(np.isposinf(exact_counts), axis=0) | np.isposinf(tail)
    regional_logs = logsumexp(
        np.concatenate((exact_counts, tail[np.newaxis, :, :]), axis=0),
        axis=0,
    ) - log_binom(k, message_weight)
    unknown_count = int(np.count_nonzero(unknown))
    if unknown_count:
        raise ArithmeticError(
            f"marker grid did not bound {unknown_count} regional entries"
        )

    scale = float(np.max(regional_logs))
    regional = np.exp(regional_logs - scale)
    log_term = (
        log_binom(k, message_weight)
        + left_degree * scale
        + log_weight_mgf(regional, left_degree, memory)
        - cutoff * math.log(output_marker)
    )
    return log_term / math.log(2.0), unknown_count, regional_logs


def activation_bands(message_weight: int, width: int) -> list[tuple[int, int]]:
    """Partition feasible activation counts, keeping zero separate."""
    if message_weight < 0 or width < 1:
        raise ValueError("invalid activation-band parameters")
    if message_weight == 0:
        return [(0, 0)]
    result = [(0, 0)]
    lo = 1
    while lo <= message_weight:
        hi = min(message_weight, lo + width - 1)
        result.append((lo, hi))
        lo = hi + 1
    return result


def activation_band_split_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
    band_width: int, log_input_min: float, log_input_max: float,
    input_steps: int, log_activation_limit: float, activation_steps: int,
    coupled_input: bool = False,
) -> tuple[float, list[tuple[int, int]], np.ndarray]:
    """Return a bivariate coefficient bound split into activation bands."""
    if k % right_degree:
        raise ValueError("right_degree must divide k")
    if input_steps < 1 or activation_steps < 1:
        raise ValueError("marker grids must be nonempty")
    region_length = k // right_degree
    size = memory + 1
    bands = activation_bands(message_weight, band_width)
    band_logs = np.full(
        (len(bands), size, size), math.inf, dtype=np.float64
    )
    input_grid = np.linspace(log_input_min, log_input_max, input_steps)
    activation_grid = np.linspace(
        -log_activation_limit, log_activation_limit, activation_steps
    )
    for input_coordinate in input_grid:
        for log_v in activation_grid:
            log_x = (
                float(input_coordinate) - float(log_v)
                if coupled_input
                else float(input_coordinate)
            )
            x = math.exp(log_x)
            v = math.exp(float(log_v))
            matrix = activation_group_matrix(
                right_degree=right_degree,
                memory=memory,
                input_marker=x,
                output_marker=output_marker,
                activation_marker=v,
            )
            base = log_matrix_power_entries(matrix, region_length)
            base -= message_weight * log_x
            for index, (lo, hi) in enumerate(bands):
                endpoint = max(-lo * log_v, -hi * log_v)
                candidate = base + math.log(hi - lo + 1) + endpoint
                band_logs[index] = np.minimum(band_logs[index], candidate)

    # From the inactive state, zero activations can only end inactive.
    band_logs[0, memory, :memory] = -math.inf
    if np.any(np.isposinf(band_logs)):
        raise ArithmeticError("marker grid left an activation band unbounded")
    regional_logs = logsumexp(band_logs, axis=0) - log_binom(
        k, message_weight
    )
    scale = float(np.max(regional_logs))
    regional = np.exp(regional_logs - scale)
    log_term = (
        log_binom(k, message_weight)
        + left_degree * scale
        + log_weight_mgf(regional, left_degree, memory)
        - cutoff * math.log(output_marker)
    )
    return log_term / math.log(2.0), bands, regional_logs


def exact_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
) -> float:
    """Return the exact floating-point first-moment term at one marker."""
    region_length = k // right_degree
    shell = biregular_parity_shell_distribution(
        left_vertices=k,
        right_degree=right_degree,
        support_size=message_weight,
    )
    region = _exact_region_matrix(
        code="ec",
        region_length=region_length,
        support_size=message_weight,
        shell=shell,
        output_marker=output_marker,
        memory=memory,
    )
    return (
        log_binom(k, message_weight)
        + log_weight_mgf(region, left_degree, memory)
        - cutoff * math.log(output_marker)
    ) / math.log(2.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=1_048_575)
    parser.add_argument("--left-degree", type=int, default=6)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--cutoff", type=int, default=230_729)
    parser.add_argument("--memory", type=int, default=60)
    parser.add_argument("--message-weight", type=int, default=64)
    parser.add_argument("--output-marker", type=float, default=0.9983461751813475)
    parser.add_argument("--activation-cutoff", type=int, default=6)
    parser.add_argument("--log-input-min", type=float, default=-13.0)
    parser.add_argument("--log-input-max", type=float, default=-7.0)
    parser.add_argument("--input-steps", type=int, default=49)
    parser.add_argument("--log-activation-limit", type=float, default=35.0)
    parser.add_argument("--activation-steps", type=int, default=71)
    parser.add_argument(
        "--band-width",
        type=int,
        default=0,
        help="if positive, replace the low-count split by full activation bands",
    )
    parser.add_argument(
        "--coupled-input",
        action="store_true",
        help="interpret the input grid as log(x*v), following the activation ridge",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    common = dict(
        k=args.k,
        left_degree=args.left_degree,
        right_degree=args.right_degree,
        cutoff=args.cutoff,
        memory=args.memory,
        message_weight=args.message_weight,
        output_marker=args.output_marker,
    )
    exact = exact_logterm(**common)
    if args.band_width:
        split, bands, regional_logs = activation_band_split_logterm(
            **common,
            band_width=args.band_width,
            log_input_min=args.log_input_min,
            log_input_max=args.log_input_max,
            input_steps=args.input_steps,
            log_activation_limit=args.log_activation_limit,
            activation_steps=args.activation_steps,
            coupled_input=args.coupled_input,
        )
        unknown = 0
        print(f"activation_bands={len(bands)}")
    else:
        split, unknown, regional_logs = activation_split_logterm(
            **common,
            activation_cutoff=args.activation_cutoff,
            log_input_min=args.log_input_min,
            log_input_max=args.log_input_max,
            input_steps=args.input_steps,
            log_activation_limit=args.log_activation_limit,
            activation_steps=args.activation_steps,
        )
    print(f"exact_log2_term={exact:.12f}")
    print(f"activation_split_log2_bound={split:.12f}")
    print(f"gap_bits={split - exact:.12f}")
    print(f"unknown_entries={unknown}")

    shell = biregular_parity_shell_distribution(
        left_vertices=args.k,
        right_degree=args.right_degree,
        support_size=args.message_weight,
    )
    exact_region = _exact_region_matrix(
        code="ec",
        region_length=args.k // args.right_degree,
        support_size=args.message_weight,
        shell=shell,
        output_marker=args.output_marker,
        memory=args.memory,
    )
    with np.errstate(divide="ignore"):
        exact_logs = np.log(exact_region)
    finite = np.isfinite(exact_logs)
    gap = (regional_logs - exact_logs) / math.log(2.0)
    print(f"regional_gap_min_bits={float(np.min(gap[finite])):.12f}")
    print(f"regional_gap_median_bits={float(np.median(gap[finite])):.12f}")
    print(f"regional_gap_max_bits={float(np.max(gap[finite])):.12f}")
    inactive_gaps = gap[args.memory, :]
    for endpoint, value in enumerate(inactive_gaps):
        if math.isfinite(float(value)):
            print(f"inactive_row_gap[{endpoint}]={float(value):.12f}")


if __name__ == "__main__":
    main()
