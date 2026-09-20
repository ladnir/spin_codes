#!/usr/bin/env python3
"""Floating-point diagnostics for region-stratified left-regular wrapped EC.

The calculation implements the Poissonized first-moment bound from
Lemma ``wrapping-regular-poisson-bound``.  It is an exploration tool, not a
certificate: a final table must recompute selected blocks with interval
arithmetic.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from expander_bounds import ec_matrix, log_binom, log_weight_mgf


@dataclass(frozen=True)
class RegularECTerm:
    message_weight: int
    activation_probability: float
    output_marker: float
    log2_bound: float


def wrapping_input_matrices(output_marker: float, memory: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the input-zero and input-one wrapping transfer matrices."""
    size = memory + 1
    zero = np.zeros((size, size), dtype=np.float64)
    one = np.zeros((size, size), dtype=np.float64)
    for state in range(memory - 1):
        zero[state, 0] = one[state, 0] = output_marker / 2.0
        zero[state, state + 1] = one[state, state + 1] = 0.5
    zero[memory - 1, 0] = output_marker
    one[memory - 1, memory] = 1.0
    zero[memory, memory] = 1.0
    one[memory, 0] = output_marker
    return zero, one


def combine_uniform_slice_transfers(
    left: tuple[int, list[np.ndarray]],
    right: tuple[int, list[np.ndarray]],
    max_weight: int,
) -> tuple[int, list[np.ndarray]]:
    """Concatenate two normalized uniform-slice transfer families."""
    left_length, left_slices = left
    right_length, right_slices = right
    total_length = left_length + right_length
    size = left_slices[0].shape[0]
    result: list[np.ndarray] = []
    for weight in range(min(max_weight, total_length) + 1):
        matrix = np.zeros((size, size), dtype=np.float64)
        denominator = math.comb(total_length, weight)
        lo = max(0, weight - right_length)
        hi = min(weight, left_length, len(left_slices) - 1)
        for left_weight in range(lo, hi + 1):
            right_weight = weight - left_weight
            if right_weight >= len(right_slices):
                continue
            probability = (
                math.comb(left_length, left_weight)
                * math.comb(right_length, right_weight)
                / denominator
            )
            matrix += probability * (
                left_slices[left_weight] @ right_slices[right_weight]
            )
        result.append(matrix)
    return total_length, result


def uniform_slice_transfer_matrices(
    *, length: int, max_weight: int, output_marker: float, memory: int
) -> list[np.ndarray]:
    """Compute [X^u] W(X,z)^length / binom(length,u) for u<=max_weight."""
    if length < 1 or not 0 <= max_weight <= length:
        raise ValueError("invalid slice-transfer dimensions")
    zero, one = wrapping_input_matrices(output_marker, memory)
    identity = np.eye(memory + 1, dtype=np.float64)
    result: tuple[int, list[np.ndarray]] = (0, [identity])
    power: tuple[int, list[np.ndarray]] = (1, [zero, one])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers(result, power, max_weight)
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers(power, power, max_weight)
    return result[1]


def regional_shell_distribution(region_length: int, rows: int) -> dict[int, float]:
    """Weight law of the parity of ``rows`` independent uniform unit vectors."""
    distribution = {0: 1.0}
    for _ in range(rows):
        next_distribution: dict[int, float] = {}
        for weight, probability in distribution.items():
            if weight:
                next_distribution[weight - 1] = next_distribution.get(weight - 1, 0.0) + (
                    probability * weight / region_length
                )
            if weight < region_length:
                next_distribution[weight + 1] = next_distribution.get(weight + 1, 0.0) + (
                    probability * (region_length - weight) / region_length
                )
        distribution = next_distribution
    return distribution


def exact_regular_ec_logterm(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    memory: int,
    message_weight: int,
) -> RegularECTerm:
    """Optimize the exact region-shell transfer for one small message weight."""
    n = region_count * region_length
    shell = regional_shell_distribution(region_length, message_weight)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        slices = uniform_slice_transfer_matrices(
            length=region_length,
            max_weight=message_weight,
            output_marker=z,
            memory=memory,
        )
        region = np.zeros_like(slices[0])
        for weight, probability in shell.items():
            region += probability * slices[weight]
        return log_weight_mgf(region, region_count, memory) - cutoff * log_z

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


def regular_activation_probability(message_weight: int, region_length: int) -> float:
    """Parity probability after throwing ``message_weight`` balls in a region."""
    if message_weight < 0 or region_length < 1:
        raise ValueError("invalid message weight or region length")
    return -0.5 * math.expm1(-2.0 * message_weight / region_length)


def log_poisson_conditioning_penalty(message_weight: int) -> float:
    """Return log(exp(r) r! / r^r), evaluated stably."""
    if message_weight < 1:
        raise ValueError("message weight must be positive")
    r = message_weight
    return math.lgamma(r + 1) + r - r * math.log(r)


def log_hamming_ball_geometric(length: int, cutoff: int) -> float:
    """Geometric-series upper bound on log(sum_{j<=cutoff} binom(length,j))."""
    if not 0 <= cutoff < length / 2:
        raise ValueError("cutoff must lie below half the length")
    return (
        log_binom(length, cutoff)
        + math.log((length - cutoff + 1) / (length - 2 * cutoff + 1))
    )


def dense_regular_logterm(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    message_weight: int,
) -> float:
    """Fourier point-mass upper bound for one dense message-weight shell."""
    n = region_count * region_length
    mixing = math.log1p(math.exp(-2.0 * message_weight / region_length))
    return (
        log_binom(k, message_weight)
        + log_hamming_ball_geometric(n, cutoff)
        + (region_count - n) * math.log(2.0)
        + n * mixing
    ) / math.log(2.0)


def regular_ec_logterm(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    memory: int,
    message_weight: int,
) -> RegularECTerm:
    """Optimize the Poissonized first-moment term for one message weight."""
    if not (1 <= message_weight <= k):
        raise ValueError("message weight must lie in [1,k]")
    if region_count < 1 or region_length < 1 or memory < 1:
        raise ValueError("invalid regular-expander parameters")
    n = region_count * region_length
    if not 0 < cutoff < n // 2:
        raise ValueError("invalid cutoff")

    q = regular_activation_probability(message_weight, region_length)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        return (
            log_weight_mgf(ec_matrix(q, z, memory, True), n, memory)
            - cutoff * log_z
        )

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


def diagnostic_weights(k: int) -> list[int]:
    """Return a compact grid that resolves sparse and dense message weights."""
    weights = set(range(1, min(k, 64) + 1))
    value = 64
    while value < k:
        value = min(k, max(value + 1, int(math.ceil(value * 1.35))))
        weights.add(value)
    for numerator in range(1, 21):
        weights.add(max(1, min(k, numerator * k // 20)))
    return sorted(weights)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=1_048_575)
    parser.add_argument("--regions", type=int, default=5)
    parser.add_argument("--region-length", type=int, default=419_430)
    parser.add_argument("--relative-cutoff", type=float, default=0.11)
    parser.add_argument("--cutoff", type=int)
    parser.add_argument("--memory", type=int, default=9)
    args = parser.parse_args()

    n = args.regions * args.region_length
    cutoff = (
        args.cutoff
        if args.cutoff is not None
        else int(math.floor(args.relative_cutoff * n))
    )
    if args.k * 2 != n:
        print(f"warning: rate is {args.k / n:.12g}, not exactly 1/2")
    print(
        f"k={args.k} n={n} rate={args.k/n:.9f} regions={args.regions} "
        f"region_length={args.region_length} memory={args.memory} "
        f"cutoff={cutoff} relative_cutoff={cutoff/n:.12f}"
    )

    terms = [
        regular_ec_logterm(
            k=args.k,
            region_length=args.region_length,
            region_count=args.regions,
            cutoff=cutoff,
            memory=args.memory,
            message_weight=r,
        )
        for r in diagnostic_weights(args.k)
    ]
    worst = max(terms, key=lambda term: term.log2_bound)
    print("r\tq_r\tz\tlog2(first-moment term)")
    for term in terms:
        if term.message_weight <= 64 or term.log2_bound > -100:
            print(
                f"{term.message_weight}\t{term.activation_probability:.9g}\t"
                f"{term.output_marker:.9g}\t{term.log2_bound:.6f}"
            )
    print(
        f"sampled worst: r={worst.message_weight}, "
        f"log2_bound={worst.log2_bound:.6f}"
    )

    dense_terms = [
        (r, dense_regular_logterm(
            k=args.k,
            region_length=args.region_length,
            region_count=args.regions,
            cutoff=cutoff,
            message_weight=r,
        ))
        for r in diagnostic_weights(args.k)
        if r >= args.k // 4
    ]
    dense_worst = max(dense_terms, key=lambda item: item[1])
    print(
        f"sampled dense Fourier worst: r={dense_worst[0]}, "
        f"log2_bound={dense_worst[1]:.6f}"
    )


if __name__ == "__main__":
    main()
