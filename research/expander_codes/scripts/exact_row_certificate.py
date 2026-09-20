#!/usr/bin/env python3
"""Rigorous finite certificates for exact-row Expand--Accumulate codes.

The verifier uses three complementary bounds. Small message weights use the
positive shell recurrence exactly. Intermediate weights condition a Bernoulli
expander on every selected row having the prescribed weight. Dense weights use
an L2 Fourier bound for the exact-row walk. Floating point selects Chernoff
markers only; all inequalities are checked with outward-rounded Arb balls.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from flint import arb, ctx

from ea_certificate import (
    decimal_marker,
    entropy_count_bound_arb,
    geometric_blocks,
    optimize_spectral_marker,
    q_arb,
    spectral_tail_bound_arb,
)


@dataclass(frozen=True)
class ExactRowVerificationResult:
    success: bool
    total_bound: arb
    security_bits: arb
    exact_bound: arb
    conditioning_bound: arb
    dense_bound: arb
    largest_exact_r: int
    largest_exact_term: arb
    conditioning_blocks: int


def accumulator_slice_tail_bound_arb(n: int, cutoff: int, shell_weight: int) -> arb:
    """Upper-bound the positive hypergeometric slice tail with Arb.

    The summand ratio decreases above the lower endpoint. Bounding the finite
    tail by the resulting infinite geometric series avoids forming hundreds
    of enormous binomial coefficients for every reachable shell.
    """
    u = shell_weight
    if u == 0:
        return arb(1)
    threshold = (u + 1) // 2
    upper = min(u, cutoff)
    if threshold > upper:
        return arb(0)
    first = (
        arb(math.comb(cutoff, threshold))
        * math.comb(n - cutoff, u - threshold)
        / math.comb(n, u)
    )
    if threshold == upper:
        return first
    ratio = (
        arb(cutoff - threshold)
        * (u - threshold)
        / ((threshold + 1) * (n - cutoff - u + threshold + 1))
    )
    if not (ratio >= 0 and ratio < 1):
        raise ValueError("geometric slice-tail ratio is not in [0,1)")
    return first / (1 - ratio)


def exact_small_weight_terms_arb(
    *, k: int, n: int, cutoff: int, row_weight: int, limit: int
) -> list[arb]:
    """Evaluate first-moment terms 1..limit by the positive shell chain."""
    denominator = math.comb(n, row_weight)
    distribution: dict[int, arb] = {0: arb(1)}
    slice_cache: dict[int, arb] = {0: arb(1)}
    transition_cache: dict[int, tuple[tuple[int, arb], ...]] = {}
    terms: list[arb] = []

    for r in range(1, limit + 1):
        next_distribution: dict[int, arb] = {}
        for u, probability in distribution.items():
            if u not in transition_cache:
                j_lo = max(0, row_weight - (n - u))
                j_hi = min(row_weight, u)
                transitions: list[tuple[int, arb]] = []
                for overlap in range(j_lo, j_hi + 1):
                    v = u + row_weight - 2 * overlap
                    numerator = math.comb(u, overlap) * math.comb(
                        n - u, row_weight - overlap
                    )
                    transitions.append((v, arb(numerator) / denominator))
                transition_cache[u] = tuple(transitions)
            for v, transition in transition_cache[u]:
                next_distribution[v] = next_distribution.get(v, arb(0)) + (
                    probability * transition
                )
        distribution = next_distribution

        fixed_tail = arb(0)
        for u, probability in distribution.items():
            if u not in slice_cache:
                slice_cache[u] = accumulator_slice_tail_bound_arb(n, cutoff, u)
            fixed_tail += probability * slice_cache[u]
        terms.append(arb(math.comb(k, r)) * fixed_tail)
    return terms


def row_condition_probability_arb(n: int, row_weight: int) -> arb:
    """Return Pr[Bin(n,row_weight/n)=row_weight]."""
    d = row_weight
    p = arb(d) / n
    return arb(math.comb(n, d)) * p**d * (1 - p) ** (n - d)


def dense_l2_bound_arb(
    *, k: int, n: int, cutoff: int, row_weight: int, start: int
) -> arb:
    """Bound all message weights at least start by exact-row L2 mixing."""
    x = arb(cutoff) / n
    entropy = -x * x.log() - (1 - x) * (1 - x).log()
    relative_ball_bound = (arb(n) * entropy).exp() / arb(2) ** n
    exponent = (
        arb(2)
        * start
        * row_weight
        * (n - row_weight)
        / (arb(n) * (n - 1))
    )
    spectral_sum_bound = arb(2) * (1 + (-exponent).exp()) ** n
    return arb(2) ** k * (relative_ball_bound * spectral_sum_bound).sqrt()


def verify_exact_row_parameters(
    *,
    k: int,
    n: int,
    cutoff: int,
    row_weight: int,
    target_bits: int = 20,
    exact_limit: int = 7,
    dense_start: int | None = None,
    block_growth: float = 1.04,
    precision_bits: int = 192,
) -> ExactRowVerificationResult:
    if dense_start is None:
        dense_start = k // 8
    if not (1 <= exact_limit < dense_start <= k):
        raise ValueError("require 1 <= exact_limit < dense_start <= k")
    if not (0 < row_weight < n // 2 and 0 < cutoff < n // 2):
        raise ValueError("invalid row weight or cutoff")

    ctx.prec = precision_bits
    exact_terms = exact_small_weight_terms_arb(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=row_weight,
        limit=exact_limit,
    )
    exact_total = sum(exact_terms, arb(0))
    largest_exact_r = max(range(1, exact_limit + 1), key=lambda r: exact_terms[r - 1])
    largest_exact_term = exact_terms[largest_exact_r - 1]

    condition_probability = row_condition_probability_arb(n, row_weight)
    conditioning_total = arb(0)
    blocks = geometric_blocks(exact_limit + 1, dense_start - 1, block_growth)
    for lo, hi in blocks:
        q_float = (1.0 - (1.0 - 2.0 * row_weight / n) ** lo) / 2.0
        marker = arb(decimal_marker(optimize_spectral_marker(q_float, n, cutoff)))
        if not (marker > 0 and marker <= 1):
            raise ValueError(f"invalid marker for conditioning block [{lo},{hi}]")
        count_bound = entropy_count_bound_arb(k, lo, hi)
        bernoulli_tail = spectral_tail_bound_arb(
            q_arb(row_weight, n, lo), marker, n, cutoff
        )
        conditioning_total += (
            count_bound * bernoulli_tail * condition_probability ** (-hi)
        )

    dense_bound = dense_l2_bound_arb(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=row_weight,
        start=dense_start,
    )
    total = exact_total + conditioning_total + dense_bound
    target = arb(2) ** (-target_bits)
    return ExactRowVerificationResult(
        success=bool(total < target),
        total_bound=total,
        security_bits=-total.log() / arb(2).log(),
        exact_bound=exact_total,
        conditioning_bound=conditioning_total,
        dense_bound=dense_bound,
        largest_exact_r=largest_exact_r,
        largest_exact_term=largest_exact_term,
        conditioning_blocks=len(blocks),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-log2", type=int, required=True)
    parser.add_argument("--rate-denominator", type=int, default=5)
    parser.add_argument("--relative-cutoff", type=float, default=0.05)
    parser.add_argument("--row-weight", type=int, required=True)
    parser.add_argument("--target-bits", type=int, default=20)
    parser.add_argument("--exact-limit", type=int, default=7)
    parser.add_argument("--dense-fraction-denominator", type=int, default=8)
    parser.add_argument("--block-growth", type=float, default=1.04)
    parser.add_argument("--precision-bits", type=int, default=192)
    args = parser.parse_args()

    k = 1 << args.k_log2
    n = args.rate_denominator * k
    cutoff = math.floor(args.relative_cutoff * n)
    result = verify_exact_row_parameters(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=args.row_weight,
        target_bits=args.target_bits,
        exact_limit=args.exact_limit,
        dense_start=k // args.dense_fraction_denominator,
        block_growth=args.block_growth,
        precision_bits=args.precision_bits,
    )
    print(f"success: {result.success}")
    print(f"total bound: {result.total_bound}")
    print(f"security bits: {result.security_bits}")
    print(f"exact small-weight bound: {result.exact_bound}")
    print(f"conditioning bound: {result.conditioning_bound}")
    print(f"dense L2 bound: {result.dense_bound}")
    print(
        f"largest exact term: r={result.largest_exact_r}, "
        f"value={result.largest_exact_term}"
    )
    print(f"conditioning blocks: {result.conditioning_blocks}")
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
