#!/usr/bin/env python3
"""Floating-point diagnostic for labeled prime-field regular EA.

This evaluates representative exact-shell, Poissonized, and dense Fourier
terms.  It is not an interval certificate and does not cover every message
support unless the caller performs a complete block scan.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import minimize_scalar

from expander_bounds import log_binom, log_weight_mgf
from prime_field_ea_diagnostic import gv_distance, regular_shell_distribution
from regular_ec_diagnostic import (
    RegularECTerm,
    combine_uniform_slice_transfers,
    diagnostic_weights,
    log_poisson_conditioning_penalty,
)


def log_projective_support_count(prime: int, k: int, weight: int) -> float:
    return log_binom(k, weight) + (weight - 1) * math.log(prime - 1)


def log_support_grouped_bound(
    prime: int, k: int, weight: int, log_line_probability_bound: float
) -> float:
    """Union over supports after capping each projective family at one."""
    log_projective_family = (weight - 1) * math.log(prime - 1)
    return log_binom(k, weight) + min(
        0.0, log_projective_family + log_line_probability_bound
    )


def combine_line_bound(
    *, prime: int, k: int, weight: int, log_line_probability_bound: float,
    group_projective_family: bool,
) -> float:
    if group_projective_family:
        return log_support_grouped_bound(
            prime, k, weight, log_line_probability_bound
        )
    return (
        log_projective_support_count(prime, k, weight)
        + log_line_probability_bound
    )


def activation_probability(prime: int, message_weight: int, region_length: int) -> float:
    s = prime - 1.0
    return s / prime * (
        -math.expm1(-prime * message_weight / (s * region_length))
    )


def accumulator_matrix(prime: int, q: float, z: float) -> np.ndarray:
    s = prime - 1.0
    return np.array(
        ((1.0 - q, q * z), (q / s, (1.0 - q / s) * z)),
        dtype=np.float64,
    )


def accumulator_input_matrices(prime: int, z: float) -> tuple[np.ndarray, np.ndarray]:
    """Input-zero and uniform-nonzero transition matrices."""
    s = prime - 1.0
    return (
        np.array([[1.0, 0.0], [0.0, z]], dtype=np.float64),
        np.array([[0.0, z], [1.0 / s, (s - 1.0) / s * z]], dtype=np.float64),
    )


def uniform_shell_transfers(
    *, prime: int, length: int, max_weight: int, output_marker: float
) -> list[np.ndarray]:
    zero, nonzero = accumulator_input_matrices(prime, output_marker)
    result: tuple[int, list[np.ndarray]] = (0, [np.eye(2, dtype=np.float64)])
    power: tuple[int, list[np.ndarray]] = (1, [zero, nonzero])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers(result, power, max_weight)
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers(power, power, max_weight)
    return result[1]


def exact_logterm(
    *, prime: int, k: int, region_length: int, region_count: int,
    cutoff: int, message_weight: int, group_projective_family: bool = True,
) -> RegularECTerm:
    shell = regular_shell_distribution(prime, region_length, message_weight)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        slices = uniform_shell_transfers(
            prime=prime,
            length=region_length,
            max_weight=message_weight,
            output_marker=z,
        )
        region = np.zeros((2, 2), dtype=np.float64)
        for weight, probability in enumerate(shell):
            if probability:
                region += probability * slices[weight]
        return log_weight_mgf(region, region_count, 0) - cutoff * log_z

    lower_log_z = -max(12.0, math.log(prime) + 8.0)
    optimum = minimize_scalar(
        objective,
        bounds=(lower_log_z, 0.0),
        method="bounded",
        options={"xatol": 2e-9, "maxiter": 140},
    )
    value = combine_line_bound(
        prime=prime,
        k=k,
        weight=message_weight,
        log_line_probability_bound=float(optimum.fun),
        group_projective_family=group_projective_family,
    )
    return RegularECTerm(
        message_weight=message_weight,
        activation_probability=activation_probability(
            prime, message_weight, region_length
        ),
        output_marker=math.exp(float(optimum.x)),
        log2_bound=value / math.log(2.0),
    )


def poisson_logterm(
    *, prime: int, k: int, region_length: int, region_count: int,
    cutoff: int, message_weight: int, group_projective_family: bool = True,
) -> RegularECTerm:
    n = region_length * region_count
    q = activation_probability(prime, message_weight, region_length)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        return log_weight_mgf(accumulator_matrix(prime, q, z), n, 0) - cutoff * log_z

    lower_log_z = -max(12.0, math.log(prime) + 8.0)
    optimum = minimize_scalar(
        objective,
        bounds=(lower_log_z, 0.0),
        method="bounded",
        options={"xatol": 2e-10, "maxiter": 180},
    )
    value = combine_line_bound(
        prime=prime,
        k=k,
        weight=message_weight,
        log_line_probability_bound=(
            region_count * log_poisson_conditioning_penalty(message_weight)
            + float(optimum.fun)
        ),
        group_projective_family=group_projective_family,
    )
    return RegularECTerm(
        message_weight=message_weight,
        activation_probability=q,
        output_marker=math.exp(float(optimum.x)),
        log2_bound=value / math.log(2.0),
    )


def log_qary_hamming_ball_geometric(prime: int, length: int, cutoff: int) -> float:
    s = prime - 1.0
    if not 0 <= cutoff < s * length / prime:
        raise ValueError("cutoff must lie below the modal shell")
    ratio = cutoff / (s * (length - cutoff + 1.0))
    return log_binom(length, cutoff) + cutoff * math.log(s) - math.log1p(-ratio)


def logaddexp(left: float, right: float) -> float:
    high = max(left, right)
    return high + math.log(math.exp(left - high) + math.exp(right - high))


def dense_logterm(
    *, prime: int, k: int, region_length: int, region_count: int,
    cutoff: int, message_weight: int, group_projective_family: bool = True,
) -> float:
    s = prime - 1.0
    n = region_length * region_count
    log_a = region_length * math.log1p(
        s * math.exp(-prime * message_weight / (s * region_length))
    )
    log_b = (
        (region_length - message_weight) * math.log(s)
        + region_length * math.log1p(
            math.exp(-prime * message_weight / region_length) / s
        )
    )
    log_line_probability_bound = (
        log_qary_hamming_ball_geometric(prime, n, cutoff)
        - n * math.log(prime)
        + region_count * logaddexp(log_a, log_b)
    )
    return combine_line_bound(
        prime=prime,
        k=k,
        weight=message_weight,
        log_line_probability_bound=log_line_probability_bound,
        group_projective_family=group_projective_family,
    ) / math.log(2.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, default=3)
    parser.add_argument("--k", type=int)
    parser.add_argument("--regions", type=int, default=48)
    parser.add_argument("--region-length", type=int, default=43_690)
    parser.add_argument("--cutoff", type=int)
    parser.add_argument("--gv-fraction", type=float, default=0.999)
    parser.add_argument("--exact-limit", type=int, default=12)
    parser.add_argument(
        "--ungrouped-projective",
        action="store_false",
        dest="group_projective_family",
        help="report the older word-by-word projective first moment",
    )
    args = parser.parse_args()
    n = args.regions * args.region_length
    k = args.k if args.k is not None else n // 2
    gv = gv_distance(args.prime, k / n)
    cutoff = args.cutoff if args.cutoff is not None else math.floor(args.gv_fraction * gv * n)
    print(
        f"p={args.prime} k={k} n={n} rate={k/n:.9f} regions={args.regions} "
        f"region_length={args.region_length} cutoff={cutoff} "
        f"relative_cutoff={cutoff/n:.12f} gv={gv:.12f} "
        f"support_grouped={args.group_projective_family}"
    )
    print("exact regional terms")
    for r in range(1, args.exact_limit + 1):
        term = exact_logterm(
            prime=args.prime, k=k, region_length=args.region_length,
            region_count=args.regions, cutoff=cutoff, message_weight=r,
            group_projective_family=args.group_projective_family,
        )
        print(f"{r}\t{term.output_marker:.9g}\t{term.log2_bound:.6f}")

    sampled = [
        poisson_logterm(
            prime=args.prime, k=k, region_length=args.region_length,
            region_count=args.regions, cutoff=cutoff, message_weight=r,
            group_projective_family=args.group_projective_family,
        )
        for r in diagnostic_weights(k)
        if r > args.exact_limit
    ]
    worst = max(sampled, key=lambda term: term.log2_bound)
    print(
        f"sampled Poisson worst: r={worst.message_weight}, "
        f"log2_bound={worst.log2_bound:.6f}"
    )
    dense = [
        (r, dense_logterm(
            prime=args.prime, k=k, region_length=args.region_length,
            region_count=args.regions, cutoff=cutoff, message_weight=r,
            group_projective_family=args.group_projective_family,
        ))
        for r in diagnostic_weights(k) if r >= k // 4
    ]
    dense_worst = max(dense, key=lambda item: item[1])
    print(
        f"sampled dense Fourier worst: r={dense_worst[0]}, "
        f"log2_bound={dense_worst[1]:.6f}"
    )


if __name__ == "__main__":
    main()
