#!/usr/bin/env python3
"""Arb certificate for region-stratified left-regular Expand--Accumulate."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from flint import arb, arb_mat, ctx
from scipy.optimize import minimize_scalar

from ea_certificate import decimal_marker, hamming_ball_bound_arb, spectral_tail_bound_arb
from expander_bounds import log_binom
from regular_ea_diagnostic import exact_regular_ea_logterm
from regular_ec_certificate import (
    binomial_shell_bound_arb,
    combine_uniform_slice_transfers_arb,
    dense_fourier_bound_arb,
    identity_matrix,
    poisson_penalty_arb,
    regional_shell_distribution_arb,
    regular_q_arb,
    zero_matrix,
)
from regular_ec_diagnostic import log_poisson_conditioning_penalty, regular_activation_probability


@dataclass(frozen=True)
class RegularEAVerificationResult:
    success: bool
    total_bound: arb
    security_bits: arb
    exact_bound: arb
    intermediate_bound: arb
    dense_bound: arb
    largest_exact_r: int
    largest_exact_term: arb
    intermediate_blocks: int
    dense_blocks: int


def accumulator_input_matrices_arb(output_marker: arb) -> tuple[arb_mat, arb_mat]:
    return (
        arb_mat([[arb(1), arb(0)], [arb(0), output_marker]]),
        arb_mat([[arb(0), output_marker], [arb(1), arb(0)]]),
    )


def accumulator_uniform_slice_transfers_arb(
    *, length: int, max_weight: int, output_marker: arb
) -> list[arb_mat]:
    zero, one = accumulator_input_matrices_arb(output_marker)
    result: tuple[int, list[arb_mat]] = (0, [identity_matrix(2)])
    power: tuple[int, list[arb_mat]] = (1, [zero, one])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers_arb(result, power, max_weight)
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers_arb(power, power, max_weight)
    return result[1]


def exact_regular_ea_tail_arb(
    *,
    region_length: int,
    region_count: int,
    cutoff: int,
    message_weight: int,
    output_marker: str,
) -> arb:
    z = arb(output_marker)
    if not (z > 0 and z <= 1):
        raise ValueError("invalid exact output marker")
    slices = accumulator_uniform_slice_transfers_arb(
        length=region_length,
        max_weight=message_weight,
        output_marker=z,
    )
    region = zero_matrix(2)
    for weight, probability in regional_shell_distribution_arb(
        region_length, message_weight
    ).items():
        region += slices[weight] * probability
    power = region**region_count
    transform = power[0, 0] + power[0, 1]
    return transform * z ** (-cutoff)


def exact_small_terms_arb(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    limit: int,
) -> list[arb]:
    terms: list[arb] = []
    for r in range(1, limit + 1):
        selected = exact_regular_ea_logterm(
            k=k,
            region_length=region_length,
            region_count=region_count,
            cutoff=cutoff,
            message_weight=r,
        )
        tail = exact_regular_ea_tail_arb(
            region_length=region_length,
            region_count=region_count,
            cutoff=cutoff,
            message_weight=r,
            output_marker=decimal_marker(selected.output_marker),
        )
        terms.append(arb(math.comb(k, r)) * tail)
    return terms


def spectral_radius_float(q: float, z: float) -> float:
    discriminant = (1.0 - q) ** 2 * (1.0 - z) ** 2 + 4.0 * q * q * z
    return 0.5 * ((1.0 - q) * (1.0 + z) + math.sqrt(discriminant))


def regular_ea_block_marker_and_logterm(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    lo: int,
    hi: int,
) -> tuple[str, float]:
    n = region_count * region_length
    q_lo = regular_activation_probability(lo, region_length)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        return (
            0.5 * math.log1p(z)
            + n * math.log(spectral_radius_float(q_lo, z))
            - cutoff * log_z
        )

    optimum = minimize_scalar(
        objective,
        bounds=(-8.0, 0.0),
        method="bounded",
        options={"xatol": 2e-10, "maxiter": 160},
    )
    log_count = math.log(hi - lo + 1) + log_binom(k, hi)
    log_term = (
        log_count
        + region_count * log_poisson_conditioning_penalty(hi)
        + float(optimum.fun)
    )
    return decimal_marker(math.exp(float(optimum.x))), log_term


def generate_intermediate_blocks(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    start: int,
    stop: int,
    target_log2: float,
) -> list[tuple[int, int, str]]:
    blocks: list[tuple[int, int, str]] = []
    target = target_log2 * math.log(2.0)
    lo = start
    while lo <= stop:
        hi = lo
        marker, value = regular_ea_block_marker_and_logterm(
            k=k,
            region_length=region_length,
            region_count=region_count,
            cutoff=cutoff,
            lo=lo,
            hi=hi,
        )
        if value > target:
            raise ValueError(f"singleton intermediate block r={lo} exceeds target")
        step = 1
        while hi < stop:
            candidate = min(stop, hi + step)
            candidate_marker, candidate_value = regular_ea_block_marker_and_logterm(
                k=k,
                region_length=region_length,
                region_count=region_count,
                cutoff=cutoff,
                lo=lo,
                hi=candidate,
            )
            if candidate_value > target:
                break
            hi, marker, value = candidate, candidate_marker, candidate_value
            step *= 2
        low = hi + 1
        high = min(stop, hi + step - 1)
        while low <= high:
            middle = (low + high) // 2
            candidate_marker, candidate_value = regular_ea_block_marker_and_logterm(
                k=k,
                region_length=region_length,
                region_count=region_count,
                cutoff=cutoff,
                lo=lo,
                hi=middle,
            )
            if candidate_value <= target:
                hi, marker, value = middle, candidate_marker, candidate_value
                low = middle + 1
            else:
                high = middle - 1
        del value
        blocks.append((lo, hi, marker))
        lo = hi + 1
    return blocks


def uniform_poisson_ea_block_tail_arb(
    *,
    n: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    lo: int,
    hi: int,
    output_marker: str,
) -> arb:
    q_lo = regular_q_arb(lo, region_length)
    z = arb(output_marker)
    if not (z > 0 and z <= 1):
        raise ValueError("invalid intermediate output marker")
    return (
        poisson_penalty_arb(hi) ** region_count
        * spectral_tail_bound_arb(q_lo, z, n, cutoff)
    )


def verify_regular_ea_parameters(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    exact_limit: int,
    dense_start: int | None = None,
    target_bits: int = 20,
    target_block_log2: float = -40.0,
    dense_block_width: int = 256,
    precision_bits: int = 192,
) -> RegularEAVerificationResult:
    if dense_start is None:
        dense_start = k // 4
    if not (1 <= exact_limit < dense_start <= k // 2):
        raise ValueError("invalid proof regimes")
    ctx.prec = precision_bits
    n = region_count * region_length

    exact_terms = exact_small_terms_arb(
        k=k,
        region_length=region_length,
        region_count=region_count,
        cutoff=cutoff,
        limit=exact_limit,
    )
    exact = sum(exact_terms, arb(0))
    largest_index = max(range(len(exact_terms)), key=lambda index: exact_terms[index])

    blocks = generate_intermediate_blocks(
        k=k,
        region_length=region_length,
        region_count=region_count,
        cutoff=cutoff,
        start=exact_limit + 1,
        stop=dense_start - 1,
        target_log2=target_block_log2,
    )
    intermediate = arb(0)
    expected = exact_limit + 1
    for lo, hi, marker in blocks:
        if lo != expected:
            raise ValueError(f"intermediate coverage breaks at r={expected}")
        expected = hi + 1
        count = arb(hi - lo + 1) * binomial_shell_bound_arb(k, hi)
        intermediate += count * uniform_poisson_ea_block_tail_arb(
            n=n,
            region_length=region_length,
            region_count=region_count,
            cutoff=cutoff,
            lo=lo,
            hi=hi,
            output_marker=marker,
        )
    if expected != dense_start:
        raise ValueError("intermediate blocks do not reach dense start")

    dense, dense_blocks = dense_fourier_bound_arb(
        k=k,
        region_length=region_length,
        region_count=region_count,
        cutoff=cutoff,
        start=dense_start,
        block_width=dense_block_width,
    )
    total = exact + intermediate + dense
    return RegularEAVerificationResult(
        success=bool(total < arb(2) ** (-target_bits)),
        total_bound=total,
        security_bits=-total.log() / arb(2).log(),
        exact_bound=exact,
        intermediate_bound=intermediate,
        dense_bound=dense,
        largest_exact_r=largest_index + 1,
        largest_exact_term=exact_terms[largest_index],
        intermediate_blocks=len(blocks),
        dense_blocks=dense_blocks,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=1_048_576)
    parser.add_argument("--regions", type=int, default=32)
    parser.add_argument("--region-length", type=int, default=163_840)
    parser.add_argument("--cutoff", type=int, default=262_144)
    parser.add_argument("--exact-limit", type=int, default=16)
    parser.add_argument("--dense-start", type=int)
    parser.add_argument("--target-bits", type=int, default=20)
    parser.add_argument("--target-block-log2", type=float, default=-40.0)
    parser.add_argument("--dense-block-width", type=int, default=256)
    parser.add_argument("--precision-bits", type=int, default=192)
    args = parser.parse_args()
    result = verify_regular_ea_parameters(
        k=args.k,
        region_length=args.region_length,
        region_count=args.regions,
        cutoff=args.cutoff,
        exact_limit=args.exact_limit,
        dense_start=args.dense_start,
        target_bits=args.target_bits,
        target_block_log2=args.target_block_log2,
        dense_block_width=args.dense_block_width,
        precision_bits=args.precision_bits,
    )
    print(f"success: {result.success}")
    print(f"total bound: {result.total_bound}")
    print(f"security bits: {result.security_bits}")
    print(f"exact contribution: {result.exact_bound}")
    print(f"intermediate contribution: {result.intermediate_bound}")
    print(f"dense contribution: {result.dense_bound}")
    print(f"largest exact term: r={result.largest_exact_r}, {result.largest_exact_term}")
    print(f"intermediate blocks: {result.intermediate_blocks}")
    print(f"dense blocks: {result.dense_blocks}")


if __name__ == "__main__":
    main()
