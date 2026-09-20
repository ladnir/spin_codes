#!/usr/bin/env python3
"""Arb certificate for region-stratified left-regular wrapped EC.

Small message weights use exact normalized regional slice matrices.
Intermediate weights use the Poissonized wrapping-transfer bound.  Dense
weights use the Fourier point-mass bound with a positive block partition.
Floating point selects markers and block endpoints only; Arb recomputes every
probability bound with outward rounding.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from flint import arb, arb_mat, ctx
from scipy.optimize import minimize_scalar

from ea_certificate import decimal_marker, hamming_ball_bound_arb
from exact_row_ec_certificate import wrapping_enumerator_matrix_arb
from expander_bounds import ec_enumerator_matrix, log_binom, log_weight_mgf
from regular_ec_diagnostic import (
    exact_regular_ec_logterm,
    log_poisson_conditioning_penalty,
    regular_activation_probability,
)


@dataclass(frozen=True)
class RegularECVerificationResult:
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


def wrapping_input_matrices_arb(
    output_marker: arb, memory: int
) -> tuple[arb_mat, arb_mat]:
    size = memory + 1
    zero = [[arb(0) for _ in range(size)] for _ in range(size)]
    one = [[arb(0) for _ in range(size)] for _ in range(size)]
    for state in range(memory - 1):
        zero[state][0] = one[state][0] = output_marker / 2
        zero[state][state + 1] = one[state][state + 1] = arb(1) / 2
    zero[memory - 1][0] = output_marker
    one[memory - 1][memory] = arb(1)
    zero[memory][memory] = arb(1)
    one[memory][0] = output_marker
    return arb_mat(zero), arb_mat(one)


def zero_matrix(size: int) -> arb_mat:
    return arb_mat([[arb(0) for _ in range(size)] for _ in range(size)])


def identity_matrix(size: int) -> arb_mat:
    return arb_mat(
        [[arb(1 if row == column else 0) for column in range(size)] for row in range(size)]
    )


def combine_uniform_slice_transfers_arb(
    left: tuple[int, list[arb_mat]],
    right: tuple[int, list[arb_mat]],
    max_weight: int,
) -> tuple[int, list[arb_mat]]:
    left_length, left_slices = left
    right_length, right_slices = right
    total_length = left_length + right_length
    size = left_slices[0].nrows()
    result: list[arb_mat] = []
    for weight in range(min(max_weight, total_length) + 1):
        matrix = zero_matrix(size)
        denominator = math.comb(total_length, weight)
        lo = max(0, weight - right_length)
        hi = min(weight, left_length, len(left_slices) - 1)
        for left_weight in range(lo, hi + 1):
            right_weight = weight - left_weight
            if right_weight >= len(right_slices):
                continue
            probability = (
                arb(math.comb(left_length, left_weight))
                * math.comb(right_length, right_weight)
                / denominator
            )
            matrix += (left_slices[left_weight] * right_slices[right_weight]) * probability
        result.append(matrix)
    return total_length, result


def uniform_slice_transfer_matrices_arb(
    *, length: int, max_weight: int, output_marker: arb, memory: int
) -> list[arb_mat]:
    zero, one = wrapping_input_matrices_arb(output_marker, memory)
    result: tuple[int, list[arb_mat]] = (0, [identity_matrix(memory + 1)])
    power: tuple[int, list[arb_mat]] = (1, [zero, one])
    remaining = length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers_arb(result, power, max_weight)
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers_arb(power, power, max_weight)
    return result[1]


def regional_shell_distribution_arb(region_length: int, rows: int) -> dict[int, arb]:
    distribution = {0: arb(1)}
    for _ in range(rows):
        next_distribution: dict[int, arb] = {}
        for weight, probability in distribution.items():
            if weight:
                next_distribution[weight - 1] = next_distribution.get(weight - 1, arb(0)) + (
                    probability * weight / region_length
                )
            if weight < region_length:
                next_distribution[weight + 1] = next_distribution.get(weight + 1, arb(0)) + (
                    probability * (region_length - weight) / region_length
                )
        distribution = next_distribution
    return distribution


def exact_regular_tail_arb(
    *,
    region_length: int,
    region_count: int,
    cutoff: int,
    memory: int,
    message_weight: int,
    output_marker: str,
) -> arb:
    z = arb(output_marker)
    if not (z > 0 and z <= 1):
        raise ValueError("invalid exact output marker")
    slices = uniform_slice_transfer_matrices_arb(
        length=region_length,
        max_weight=message_weight,
        output_marker=z,
        memory=memory,
    )
    region = zero_matrix(memory + 1)
    for weight, probability in regional_shell_distribution_arb(
        region_length, message_weight
    ).items():
        region += slices[weight] * probability
    power = region**region_count
    transform = sum((power[memory, column] for column in range(memory + 1)), arb(0))
    return transform * z ** (-cutoff)


def exact_small_terms_arb(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    memory: int,
    limit: int,
) -> list[arb]:
    terms: list[arb] = []
    for r in range(1, limit + 1):
        selected = exact_regular_ec_logterm(
            k=k,
            region_length=region_length,
            region_count=region_count,
            cutoff=cutoff,
            memory=memory,
            message_weight=r,
        )
        tail = exact_regular_tail_arb(
            region_length=region_length,
            region_count=region_count,
            cutoff=cutoff,
            memory=memory,
            message_weight=r,
            output_marker=decimal_marker(selected.output_marker),
        )
        terms.append(arb(math.comb(k, r)) * tail)
    return terms


def regular_block_marker_and_logterm(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    memory: int,
    lo: int,
    hi: int,
) -> tuple[str, float]:
    n = region_count * region_length
    q_lo = regular_activation_probability(lo, region_length)
    q_hi = regular_activation_probability(hi, region_length)
    input_marker = q_hi / (1.0 - q_hi)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        matrix = (1.0 - q_lo) * ec_enumerator_matrix(
            input_marker, z, memory, True
        )
        return log_weight_mgf(matrix, n, memory) - cutoff * log_z

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
    memory: int,
    start: int,
    stop: int,
    target_log2: float,
) -> list[tuple[int, int, str]]:
    blocks: list[tuple[int, int, str]] = []
    target = target_log2 * math.log(2.0)
    lo = start
    while lo <= stop:
        hi = lo
        marker, value = regular_block_marker_and_logterm(
            k=k,
            region_length=region_length,
            region_count=region_count,
            cutoff=cutoff,
            memory=memory,
            lo=lo,
            hi=hi,
        )
        if value > target:
            raise ValueError(f"singleton intermediate block r={lo} exceeds target")
        step = 1
        while hi < stop:
            candidate = min(stop, hi + step)
            candidate_marker, candidate_value = regular_block_marker_and_logterm(
                k=k,
                region_length=region_length,
                region_count=region_count,
                cutoff=cutoff,
                memory=memory,
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
            candidate_marker, candidate_value = regular_block_marker_and_logterm(
                k=k,
                region_length=region_length,
                region_count=region_count,
                cutoff=cutoff,
                memory=memory,
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


def regular_q_arb(message_weight: int, region_length: int) -> arb:
    return (1 - (-arb(2 * message_weight) / region_length).exp()) / 2


def poisson_penalty_arb(message_weight: int) -> arb:
    r = message_weight
    return arb(r).exp() * arb(math.factorial(r)) / arb(r) ** r


def uniform_poisson_block_tail_arb(
    *,
    n: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    memory: int,
    lo: int,
    hi: int,
    output_marker: str,
) -> arb:
    q_lo = regular_q_arb(lo, region_length)
    q_hi = regular_q_arb(hi, region_length)
    input_marker = q_hi / (1 - q_hi)
    z = arb(output_marker)
    matrix = wrapping_enumerator_matrix_arb(input_marker, z, memory) * (1 - q_lo)
    power = matrix**n
    transform = sum((power[memory, column] for column in range(memory + 1)), arb(0))
    return (
        poisson_penalty_arb(hi) ** region_count
        * transform
        * z ** (-cutoff)
    )


def binomial_shell_bound_arb(length: int, weight: int) -> arb:
    return (
        arb(length + 1).lgamma()
        - arb(weight + 1).lgamma()
        - arb(length - weight + 1).lgamma()
    ).exp()


def dense_fourier_bound_arb(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    start: int,
    block_width: int,
) -> tuple[arb, int]:
    n = region_count * region_length
    common = hamming_ball_bound_arb(n, cutoff) * arb(2) ** (region_count - n)
    total = arb(0)
    blocks = 0
    lo = start
    while lo <= k:
        hi = min(k, lo + block_width - 1)
        if hi < k // 2:
            peak = hi
        elif lo > k // 2:
            peak = lo
        else:
            peak = k // 2
        count = arb(hi - lo + 1) * binomial_shell_bound_arb(k, peak)
        mixing = (1 + (-arb(2 * lo) / region_length).exp()) ** n
        total += common * count * mixing
        blocks += 1
        lo = hi + 1
    return total, blocks


def verify_regular_ec_parameters(
    *,
    k: int,
    region_length: int,
    region_count: int,
    cutoff: int,
    memory: int,
    exact_limit: int = 10,
    dense_start: int | None = None,
    target_bits: int = 20,
    target_block_log2: float = -40.0,
    dense_block_width: int = 256,
    precision_bits: int = 192,
) -> RegularECVerificationResult:
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
        memory=memory,
        limit=exact_limit,
    )
    exact = sum(exact_terms, arb(0))
    largest_index = max(range(len(exact_terms)), key=lambda index: exact_terms[index])

    blocks = generate_intermediate_blocks(
        k=k,
        region_length=region_length,
        region_count=region_count,
        cutoff=cutoff,
        memory=memory,
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
        intermediate += count * uniform_poisson_block_tail_arb(
            n=n,
            region_length=region_length,
            region_count=region_count,
            cutoff=cutoff,
            memory=memory,
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
    return RegularECVerificationResult(
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
    parser.add_argument("--k", type=int, default=1_048_572)
    parser.add_argument("--regions", type=int, default=28)
    parser.add_argument("--region-length", type=int, default=74_898)
    parser.add_argument("--cutoff", type=int, default=230_730)
    parser.add_argument("--memory", type=int, default=9)
    parser.add_argument("--exact-limit", type=int, default=10)
    parser.add_argument("--dense-start", type=int)
    parser.add_argument("--target-bits", type=int, default=20)
    parser.add_argument("--target-block-log2", type=float, default=-40.0)
    parser.add_argument("--dense-block-width", type=int, default=256)
    parser.add_argument("--precision-bits", type=int, default=192)
    args = parser.parse_args()
    result = verify_regular_ec_parameters(
        k=args.k,
        region_length=args.region_length,
        region_count=args.regions,
        cutoff=args.cutoff,
        memory=args.memory,
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
