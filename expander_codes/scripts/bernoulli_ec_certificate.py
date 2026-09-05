#!/usr/bin/env python3
"""Arb verifier for Bernoulli-expander wrapped Expand--Convolute."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from flint import arb, arb_mat, ctx
from scipy.optimize import minimize_scalar

from ea_certificate import (
    decimal_marker,
    entropy_count_bound_arb,
    hamming_ball_bound_arb,
)
from exact_row_ec_certificate import wrapping_enumerator_matrix_arb
from expander_bounds import (
    activation_probability,
    ec_enumerator_matrix,
    log_binom,
    log_tail_bound,
    log_weight_mgf,
)


@dataclass(frozen=True)
class BernoulliECVerificationResult:
    total_bound: arb
    security_bits: arb
    exact_bound: arb
    intermediate_bound: arb
    dense_bound: arb
    largest_exact_r: int
    largest_exact_term: arb
    intermediate_blocks: int
    dense_blocks: int


def density_arb(expected_row_coefficient: str, n: int) -> arb:
    """Return p=C ln(n)/n with outward rounding."""
    return arb(expected_row_coefficient) * arb(n).log() / n


def q_from_density_arb(density: arb, message_weight: int) -> arb:
    """Return the parity activation probability for one output coordinate."""
    return (arb(1) - (arb(1) - 2 * density) ** message_weight) / 2


def bernoulli_wrapping_matrix_arb(q: arb, z: arb, memory: int) -> arb_mat:
    """Return the exact Bernoulli-input wrapping matrix N_{memory,q}(z)."""
    size = memory + 1
    entries = [[arb(0) for _ in range(size)] for _ in range(size)]
    for state in range(memory - 1):
        entries[state][0] = z / 2
        entries[state][state + 1] = arb(1) / 2
    entries[memory - 1][0] = (1 - q) * z
    entries[memory - 1][memory] = q
    entries[memory][0] = q * z
    entries[memory][memory] = 1 - q
    return arb_mat(entries)


def exact_tail_bound_arb(*, q: arb, z: arb, n: int, cutoff: int, memory: int) -> arb:
    """Verify one fixed-message Chernoff bound."""
    power = bernoulli_wrapping_matrix_arb(q, z, memory) ** n
    transform = sum((power[memory, column] for column in range(memory + 1)), arb(0))
    return transform * z ** (-cutoff)


def exact_terms_arb(
    *,
    k: int,
    n: int,
    cutoff: int,
    expected_row_coefficient: str,
    memory: int,
    limit: int,
) -> list[arb]:
    """Verify message weights one through limit individually."""
    density_float = float(expected_row_coefficient) * math.log(n) / n
    density = density_arb(expected_row_coefficient, n)
    terms: list[arb] = []
    for message_weight in range(1, limit + 1):
        q_float = activation_probability(density_float, message_weight)
        _, marker = log_tail_bound(
            code="ec",
            q=q_float,
            length=n,
            cutoff=cutoff,
            memory=memory,
            wrapping=True,
        )
        z = arb(decimal_marker(marker))
        q = q_from_density_arb(density, message_weight)
        tail = exact_tail_bound_arb(q=q, z=z, n=n, cutoff=cutoff, memory=memory)
        terms.append(arb(math.comb(k, message_weight)) * tail)
    return terms


def block_marker_and_logterm(
    *,
    k: int,
    n: int,
    cutoff: int,
    expected_row_coefficient: float,
    memory: int,
    lo: int,
    hi: int,
) -> tuple[str, float]:
    """Select a decimal marker and estimate one uniform transfer block."""
    density = expected_row_coefficient * math.log(n) / n
    q_lo = activation_probability(density, lo)
    q_hi = activation_probability(density, hi)
    input_marker = q_hi / (1.0 - q_hi)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        matrix = (1.0 - q_lo) * ec_enumerator_matrix(
            input_marker, z, memory, True
        )
        return log_weight_mgf(matrix, n, memory) - cutoff * log_z

    optimum = minimize_scalar(
        objective,
        bounds=(-1.0, 0.0),
        method="bounded",
        options={"xatol": 1e-10, "maxiter": 160},
    )
    log_count = math.log(hi - lo + 1) + log_binom(k, hi)
    return decimal_marker(math.exp(float(optimum.x))), log_count + float(optimum.fun)


def generate_blocks(
    *,
    k: int,
    n: int,
    cutoff: int,
    expected_row_coefficient: float,
    memory: int,
    start: int,
    stop: int,
    target_log2: float,
) -> list[tuple[int, int, str]]:
    """Greedily partition the intermediate message weights."""
    blocks: list[tuple[int, int, str]] = []
    lo = start
    target = target_log2 * math.log(2.0)
    while lo <= stop:
        hi = lo
        marker, value = block_marker_and_logterm(
            k=k,
            n=n,
            cutoff=cutoff,
            expected_row_coefficient=expected_row_coefficient,
            memory=memory,
            lo=lo,
            hi=hi,
        )
        step = 1
        while hi < stop:
            candidate = min(stop, hi + step)
            candidate_marker, candidate_value = block_marker_and_logterm(
                k=k,
                n=n,
                cutoff=cutoff,
                expected_row_coefficient=expected_row_coefficient,
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
            candidate_marker, candidate_value = block_marker_and_logterm(
                k=k,
                n=n,
                cutoff=cutoff,
                expected_row_coefficient=expected_row_coefficient,
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


def uniform_block_tail_arb(
    *,
    n: int,
    cutoff: int,
    density: arb,
    memory: int,
    lo: int,
    hi: int,
    output_marker: str,
) -> arb:
    """Verify the transfer tail uniformly for message weights lo through hi."""
    q_lo = q_from_density_arb(density, lo)
    q_hi = q_from_density_arb(density, hi)
    input_marker = q_hi / (1 - q_hi)
    z = arb(output_marker)
    matrix = wrapping_enumerator_matrix_arb(input_marker, z, memory) * (1 - q_lo)
    power = matrix**n
    transform = sum((power[memory, column] for column in range(memory + 1)), arb(0))
    return transform * z ** (-cutoff)


def dense_bound_arb(
    *, k: int, n: int, cutoff: int, density: arb, start: int
) -> tuple[arb, int]:
    """Bound dense messages with a short activation-probability partition."""
    if not (1 <= start <= k // 2):
        raise ValueError("dense start must lie in [1,k/2]")

    output_ball = hamming_ball_bound_arb(n, cutoff)
    two_to_k = arb(2) ** k
    # The points become finer near k/2, where the message count grows quickly.
    # Every endpoint is an integer function of k, so certificate coverage is
    # reproducible and independent of floating-point rounding.
    candidate_stops = (
        3 * k // 10,
        7 * k // 20,
        2 * k // 5,
        9 * k // 20,
        47 * k // 100,
        12 * k // 25,
        49 * k // 100,
        99 * k // 200,
        k // 2,
    )
    first_fixed_stop = 3 * k // 10
    stops: list[int] = []

    # When the dense regime starts below 0.3k, select a few additional stops.
    # The floating-point calculation only chooses endpoints. Arb verifies the
    # contribution of every resulting block below.
    if start < first_fixed_stop and k >= 1024:
        density_float = float(density)
        output_log_ball = (
            math.lgamma(n + 1)
            - math.lgamma(cutoff + 1)
            - math.lgamma(n - cutoff + 1)
            + math.log((n - cutoff + 1) / (n - 2 * cutoff + 1))
        )
        target = -40.0 * math.log(2.0)

        def estimated_log_block(lo: int, stop: int) -> float:
            hi = stop - 1
            log_count = (
                math.lgamma(k + 1)
                - math.lgamma(hi + 1)
                - math.lgamma(k - hi + 1)
                + math.log((k - hi + 1) / (k - 2 * hi + 1))
            )
            q_lo = activation_probability(density_float, lo)
            return log_count + output_log_ball + n * math.log1p(-q_lo)

        lo = start
        while lo < first_fixed_stop:
            if estimated_log_block(lo, lo + 1) > target:
                raise ValueError(
                    "dense start is too early for the point-mass partition"
                )
            stop = lo + 1
            step = 1
            while stop < first_fixed_stop:
                candidate = min(first_fixed_stop, stop + step)
                if estimated_log_block(lo, candidate) > target:
                    break
                stop = candidate
                step *= 2
            low = stop + 1
            high = min(first_fixed_stop, stop + step - 1)
            while low <= high:
                middle = (low + high) // 2
                if estimated_log_block(lo, middle) <= target:
                    stop = middle
                    low = middle + 1
                else:
                    high = middle - 1
            stops.append(stop)
            lo = stop

    stops.extend(
        stop for stop in candidate_stops if start < stop <= k // 2
    )
    stops = sorted(set(stops))

    total = arb(0)
    lo = start
    block_count = 0
    for stop in stops:
        hi = stop - 1
        count = hamming_ball_bound_arb(k, hi)
        if count > two_to_k:
            count = two_to_k
        q_lo = q_from_density_arb(density, lo)
        total += count * output_ball * (1 - q_lo) ** n
        lo = stop
        block_count += 1

    # The final block includes the upper half of the message weights. Its
    # cardinality is at most 2^k, and q_r is nondecreasing in r.
    q_lo = q_from_density_arb(density, lo)
    total += two_to_k * output_ball * (1 - q_lo) ** n
    return total, block_count + 1


def verify_parameters(
    *,
    k: int,
    n: int,
    cutoff: int,
    expected_row_coefficient: str,
    memory: int,
    exact_limit: int = 32,
    dense_start: int | None = None,
    target_block_log2: float = -40.0,
    precision_bits: int = 192,
) -> BernoulliECVerificationResult:
    """Verify all nonzero message weights for one Bernoulli wrapped-EC set."""
    if dense_start is None:
        dense_start = k // 4
    if not (1 <= exact_limit < dense_start <= k // 2):
        raise ValueError("require 1 <= exact_limit < dense_start <= k/2")
    if not (0 < cutoff < n // 2 and memory >= 1):
        raise ValueError("invalid cutoff or memory")
    ctx.prec = precision_bits

    exact_terms = exact_terms_arb(
        k=k,
        n=n,
        cutoff=cutoff,
        expected_row_coefficient=expected_row_coefficient,
        memory=memory,
        limit=exact_limit,
    )
    exact = sum(exact_terms, arb(0))
    largest_index = max(range(len(exact_terms)), key=lambda index: exact_terms[index])

    coefficient_float = float(expected_row_coefficient)
    blocks = generate_blocks(
        k=k,
        n=n,
        cutoff=cutoff,
        expected_row_coefficient=coefficient_float,
        memory=memory,
        start=exact_limit + 1,
        stop=dense_start - 1,
        target_log2=target_block_log2,
    )
    density = density_arb(expected_row_coefficient, n)
    intermediate = arb(0)
    expected = exact_limit + 1
    for lo, hi, marker in blocks:
        if lo != expected or hi < lo:
            raise ValueError(f"intermediate coverage breaks at r={expected}")
        expected = hi + 1
        count = entropy_count_bound_arb(k, lo, hi)
        intermediate += count * uniform_block_tail_arb(
            n=n,
            cutoff=cutoff,
            density=density,
            memory=memory,
            lo=lo,
            hi=hi,
            output_marker=marker,
        )
    if expected != dense_start:
        raise ValueError("intermediate blocks do not reach dense_start")

    dense, dense_blocks = dense_bound_arb(
        k=k,
        n=n,
        cutoff=cutoff,
        density=density,
        start=dense_start,
    )
    total = exact + intermediate + dense
    return BernoulliECVerificationResult(
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
    parser.add_argument("--k-log2", type=int, default=20)
    parser.add_argument("--rate-denominator", type=int, default=5)
    parser.add_argument("--relative-cutoff", type=float, default=0.05)
    parser.add_argument("--cutoff", type=int)
    parser.add_argument("--coefficient", required=True)
    parser.add_argument("--memory", type=int, default=5)
    parser.add_argument("--exact-limit", type=int, default=32)
    parser.add_argument("--dense-fraction-denominator", type=int, default=4)
    parser.add_argument("--dense-start", type=int)
    parser.add_argument("--target-block-log2", type=float, default=-40.0)
    parser.add_argument("--precision-bits", type=int, default=192)
    args = parser.parse_args()

    k = 1 << args.k_log2
    n = args.rate_denominator * k
    cutoff = (
        args.cutoff
        if args.cutoff is not None
        else math.floor(args.relative_cutoff * n)
    )
    result = verify_parameters(
        k=k,
        n=n,
        cutoff=cutoff,
        expected_row_coefficient=args.coefficient,
        memory=args.memory,
        exact_limit=args.exact_limit,
        dense_start=(
            args.dense_start
            if args.dense_start is not None
            else k // args.dense_fraction_denominator
        ),
        target_block_log2=args.target_block_log2,
        precision_bits=args.precision_bits,
    )
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
