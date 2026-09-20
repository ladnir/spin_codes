#!/usr/bin/env python3
"""Floating-point diagnostics for region-stratified left-regular EA."""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import minimize_scalar

from expander_bounds import ea_matrix, log_binom, log_weight_mgf
from regular_ec_diagnostic import (
    RegularECTerm,
    combine_uniform_slice_transfers,
    dense_regular_logterm,
    diagnostic_weights,
    log_poisson_conditioning_penalty,
    regular_activation_probability,
    regional_shell_distribution,
)


def accumulator_input_matrices(output_marker: float) -> tuple[np.ndarray, np.ndarray]:
    """Return the input-zero and input-one accumulator transfer matrices."""
    return (
        np.array([[1.0, 0.0], [0.0, output_marker]], dtype=np.float64),
        np.array([[0.0, output_marker], [1.0, 0.0]], dtype=np.float64),
    )


def accumulator_uniform_slice_transfers(
    *, length: int, max_weight: int, output_marker: float
) -> list[np.ndarray]:
    zero, one = accumulator_input_matrices(output_marker)
    result: tuple[int, list[np.ndarray]] = (0, [np.eye(2, dtype=np.float64)])
    power: tuple[int, list[np.ndarray]] = (1, [zero, one])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers(result, power, max_weight)
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers(power, power, max_weight)
    return result[1]


def exact_regular_ea_logterm(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    message_weight: int,
) -> RegularECTerm:
    shell = regional_shell_distribution(region_length, message_weight)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        slices = accumulator_uniform_slice_transfers(
            length=region_length,
            max_weight=message_weight,
            output_marker=z,
        )
        region = np.zeros((2, 2), dtype=np.float64)
        for weight, probability in shell.items():
            region += probability * slices[weight]
        return log_weight_mgf(region, region_count, 0) - cutoff * log_z

    optimum = minimize_scalar(
        objective,
        bounds=(-8.0, 0.0),
        method="bounded",
        options={"xatol": 2e-9, "maxiter": 100},
    )
    log_bound = log_binom(k, message_weight) + float(optimum.fun)
    return RegularECTerm(
        message_weight=message_weight,
        activation_probability=regular_activation_probability(
            message_weight, region_length
        ),
        output_marker=math.exp(float(optimum.x)),
        log2_bound=log_bound / math.log(2.0),
    )


def poisson_regular_ea_logterm(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    message_weight: int,
) -> RegularECTerm:
    n = region_count * region_length
    q = regular_activation_probability(message_weight, region_length)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        return log_weight_mgf(ea_matrix(q, z), n, 0) - cutoff * log_z

    optimum = minimize_scalar(
        objective,
        bounds=(-8.0, 0.0),
        method="bounded",
        options={"xatol": 2e-10, "maxiter": 160},
    )
    log_bound = (
        log_binom(k, message_weight)
        + region_count * log_poisson_conditioning_penalty(message_weight)
        + float(optimum.fun)
    )
    return RegularECTerm(
        message_weight=message_weight,
        activation_probability=q,
        output_marker=math.exp(float(optimum.x)),
        log2_bound=log_bound / math.log(2.0),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=1_048_572)
    parser.add_argument("--regions", type=int, default=28)
    parser.add_argument("--region-length", type=int, default=74_898)
    parser.add_argument("--cutoff", type=int, default=230_730)
    parser.add_argument("--exact-limit", type=int, default=10)
    args = parser.parse_args()
    n = args.regions * args.region_length

    print(
        f"k={args.k} n={n} rate={args.k/n:.9f} regions={args.regions} "
        f"region_length={args.region_length} cutoff={args.cutoff} "
        f"relative_cutoff={args.cutoff/n:.12f}"
    )
    print("exact regional terms")
    for r in range(1, args.exact_limit + 1):
        term = exact_regular_ea_logterm(
            k=args.k,
            region_length=args.region_length,
            region_count=args.regions,
            cutoff=args.cutoff,
            message_weight=r,
        )
        print(f"{r}\t{term.output_marker:.9g}\t{term.log2_bound:.6f}")

    terms = [
        poisson_regular_ea_logterm(
            k=args.k,
            region_length=args.region_length,
            region_count=args.regions,
            cutoff=args.cutoff,
            message_weight=r,
        )
        for r in diagnostic_weights(args.k)
        if r > args.exact_limit
    ]
    worst = max(terms, key=lambda term: term.log2_bound)
    print(
        f"sampled Poisson worst: r={worst.message_weight}, "
        f"log2_bound={worst.log2_bound:.6f}"
    )
    dense = [
        (
            r,
            dense_regular_logterm(
                k=args.k,
                region_length=args.region_length,
                region_count=args.regions,
                cutoff=args.cutoff,
                message_weight=r,
            ),
        )
        for r in diagnostic_weights(args.k)
        if r >= args.k // 4
    ]
    dense_worst = max(dense, key=lambda item: item[1])
    print(
        f"sampled dense Fourier worst: r={dense_worst[0]}, "
        f"log2_bound={dense_worst[1]:.6f}"
    )


if __name__ == "__main__":
    main()
