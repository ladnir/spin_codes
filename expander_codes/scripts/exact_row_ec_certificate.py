#!/usr/bin/env python3
"""All-weight Arb verifier for wrapped exact-row Expand--Convolute.

The verifier composes the positive exact-row shell recurrence with the exact
two-marker wrapping transfer matrix. A floating-point heuristic selects marker
values, but the trusted calculation parses them as fixed decimals and
recomputes every bound with outward-rounded Arb balls.

Small weights use the positive shell recurrence. Intermediate weights use a
uniform conditioned-Bernoulli transfer bound. Dense weights use the L2 bound
that applies to every invertible inner map.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from flint import arb, arb_mat, ctx
from scipy.optimize import minimize_scalar

from ea_certificate import decimal_marker, entropy_count_bound_arb, q_arb
from exact_row_certificate import dense_l2_bound_arb, row_condition_probability_arb
from expander_bounds import (
    activation_probability,
    ec_enumerator_matrix,
    log_binom,
    log_weight_mgf,
)


@dataclass(frozen=True)
class SmallWeightECResult:
    small_weight_bound: arb
    security_bits: arb
    largest_r: int
    largest_term: arb
    terms: tuple[arb, ...]
    slice_count: int


@dataclass(frozen=True)
class ExactRowECVerificationResult:
    success: bool
    total_bound: arb
    security_bits: arb
    small_weight_bound: arb
    intermediate_bound: arb
    dense_bound: arb
    largest_small_r: int
    largest_small_term: arb
    slice_count: int
    intermediate_blocks: int


def wrapping_enumerator_matrix_arb(
    input_marker: arb, output_marker: arb, memory: int
) -> arb_mat:
    """Return the exact bivariate wrapping matrix over Arb balls."""
    size = memory + 1
    entries = [[arb(0) for _ in range(size)] for _ in range(size)]
    for ell in range(memory - 1):
        entries[ell][0] = (1 + input_marker) * output_marker / 2
        entries[ell][ell + 1] = (1 + input_marker) / 2
    entries[memory - 1][0] = output_marker
    entries[memory - 1][memory] = input_marker
    entries[memory][0] = input_marker * output_marker
    entries[memory][memory] = arb(1)
    return arb_mat(entries)


def wrapping_slice_tail_bound_arb(
    *,
    n: int,
    cutoff: int,
    shell_weight: int,
    memory: int,
    input_marker: str,
    output_marker: str,
) -> arb:
    """Check one fixed-shell wrapping tail bound with Arb."""
    u = shell_weight
    if u == 0:
        return arb(1)
    x = arb(input_marker)
    z = arb(output_marker)
    if not (x > 0 and z > 0 and z <= 1):
        raise ValueError("invalid two-marker certificate")
    power = wrapping_enumerator_matrix_arb(x, z, memory) ** n
    transform = sum((power[memory, j] for j in range(memory + 1)), arb(0))
    marker_bound = transform * x ** (-u) * z ** (-cutoff) / math.comb(n, u)
    # Probability one is always a valid slice bound. This case matters when a
    # rare shell is too small for the two-marker optimization to improve it.
    if marker_bound > 1:
        return arb(1)
    return marker_bound


def selected_markers(
    *, n: int, cutoff: int, shell_weight: int, memory: int
) -> tuple[str, str]:
    """Select decimal markers outside the trusted Arb calculation."""
    del memory
    # These formulas are close to the numerical saddle at relative cutoff
    # 0.05. Optimality is irrelevant: any positive x and z<=1 give a valid
    # bound when the Arb verifier substitutes their fixed decimal values.
    input_marker = math.exp(math.log(shell_weight / n) + 2.238)
    output_marker = math.exp(-shell_weight / cutoff)
    return decimal_marker(input_marker), decimal_marker(output_marker)


def exact_small_weight_terms_arb(
    *,
    k: int,
    n: int,
    cutoff: int,
    row_weight: int,
    memory: int,
    limit: int,
) -> tuple[list[arb], int]:
    """Verify first-moment terms for message weights one through ``limit``."""
    denominator = math.comb(n, row_weight)
    distribution: dict[int, arb] = {0: arb(1)}
    transition_cache: dict[int, tuple[tuple[int, arb], ...]] = {}
    slice_cache: dict[int, arb] = {0: arb(1)}
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
                input_marker, output_marker = selected_markers(
                    n=n,
                    cutoff=cutoff,
                    shell_weight=u,
                    memory=memory,
                )
                slice_cache[u] = wrapping_slice_tail_bound_arb(
                    n=n,
                    cutoff=cutoff,
                    shell_weight=u,
                    memory=memory,
                    input_marker=input_marker,
                    output_marker=output_marker,
                )
            fixed_tail += probability * slice_cache[u]
        terms.append(arb(math.comb(k, r)) * fixed_tail)
    return terms, len(slice_cache)


def verify_small_weight_parameters(
    *,
    k: int,
    n: int,
    cutoff: int,
    row_weight: int,
    memory: int,
    limit: int,
    precision_bits: int = 192,
) -> SmallWeightECResult:
    """Return the certified contribution of message weights ``1..limit``."""
    if not (1 <= limit <= k):
        raise ValueError("limit must lie in [1,k]")
    if not (0 < row_weight < n // 2 and 0 < cutoff < n // 2 and memory >= 1):
        raise ValueError("invalid EC parameters")
    ctx.prec = precision_bits
    terms, slice_count = exact_small_weight_terms_arb(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=row_weight,
        memory=memory,
        limit=limit,
    )
    total = sum(terms, arb(0))
    largest_r = max(range(1, limit + 1), key=lambda r: terms[r - 1])
    return SmallWeightECResult(
        small_weight_bound=total,
        security_bits=-total.log() / arb(2).log(),
        largest_r=largest_r,
        largest_term=terms[largest_r - 1],
        terms=tuple(terms),
        slice_count=slice_count,
    )


def uniform_block_marker_and_logterm(
    *,
    k: int,
    n: int,
    cutoff: int,
    row_weight: int,
    memory: int,
    lo: int,
    hi: int,
) -> tuple[str, float]:
    """Select a marker and estimate one conditioned-Bernoulli block."""
    density = row_weight / n
    q_lo = activation_probability(density, lo)
    q_hi = activation_probability(density, hi)
    input_marker = q_hi / (1.0 - q_hi)

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        # Form the scaled matrix before exponentiation. Evaluating
        # n*log(1-q_lo)+log(F) separately loses many bits at large n.
        matrix = (1.0 - q_lo) * ec_enumerator_matrix(
            input_marker, z, memory, True
        )
        transform = log_weight_mgf(
            matrix, n, memory
        )
        return transform - cutoff * log_z

    optimum = minimize_scalar(
        objective,
        bounds=(-1.0, 0.0),
        method="bounded",
        options={"xatol": 1e-10, "maxiter": 160},
    )
    log_count = math.log(hi - lo + 1) + log_binom(k, hi)
    p = row_weight / n
    log_condition = (
        log_binom(n, row_weight)
        + row_weight * math.log(p)
        + (n - row_weight) * math.log1p(-p)
    )
    log_term = log_count - hi * log_condition + float(optimum.fun)
    return decimal_marker(math.exp(float(optimum.x))), log_term


def generate_intermediate_blocks(
    *,
    k: int,
    n: int,
    cutoff: int,
    row_weight: int,
    memory: int,
    start: int,
    stop: int,
    target_log2: float = -34.0,
) -> list[tuple[int, int, str]]:
    """Greedily select a contiguous partition for the middle range."""
    blocks: list[tuple[int, int, str]] = []
    lo = start
    target = target_log2 * math.log(2.0)
    while lo <= stop:
        hi = lo
        marker, value = uniform_block_marker_and_logterm(
            k=k,
            n=n,
            cutoff=cutoff,
            row_weight=row_weight,
            memory=memory,
            lo=lo,
            hi=hi,
        )
        step = 1
        while hi < stop:
            candidate = min(stop, hi + step)
            candidate_marker, candidate_value = uniform_block_marker_and_logterm(
                k=k,
                n=n,
                cutoff=cutoff,
                row_weight=row_weight,
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
            candidate_marker, candidate_value = uniform_block_marker_and_logterm(
                k=k,
                n=n,
                cutoff=cutoff,
                row_weight=row_weight,
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


def uniform_bernoulli_block_tail_arb(
    *,
    n: int,
    cutoff: int,
    row_weight: int,
    memory: int,
    lo: int,
    hi: int,
    output_marker: str,
) -> arb:
    """Verify one tail bound uniformly for message weights ``lo..hi``."""
    q_lo = q_arb(row_weight, n, lo)
    q_hi = q_arb(row_weight, n, hi)
    input_marker = q_hi / (1 - q_hi)
    z = arb(output_marker)
    if not (z > 0 and z <= 1):
        raise ValueError("invalid output marker")
    matrix = wrapping_enumerator_matrix_arb(input_marker, z, memory) * (1 - q_lo)
    power = matrix ** n
    transform = sum((power[memory, j] for j in range(memory + 1)), arb(0))
    return transform * z ** (-cutoff)


def verify_exact_row_ec_parameters(
    *,
    k: int,
    n: int,
    cutoff: int,
    row_weight: int,
    memory: int,
    target_bits: int = 20,
    exact_limit: int = 59,
    dense_start: int | None = None,
    precision_bits: int = 192,
) -> ExactRowECVerificationResult:
    """Verify all message weights for one exact-row wrapped-EC parameter set."""
    if dense_start is None:
        dense_start = k // 4
    if not (1 <= exact_limit < dense_start <= k // 2):
        raise ValueError("require 1 <= exact_limit < dense_start <= k/2")
    ctx.prec = precision_bits

    small = verify_small_weight_parameters(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=row_weight,
        memory=memory,
        limit=exact_limit,
        precision_bits=precision_bits,
    )
    blocks = generate_intermediate_blocks(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=row_weight,
        memory=memory,
        start=exact_limit + 1,
        stop=dense_start - 1,
    )
    expected = exact_limit + 1
    condition_probability = row_condition_probability_arb(n, row_weight)
    intermediate = arb(0)
    for lo, hi, marker in blocks:
        if lo != expected or hi < lo:
            raise ValueError(f"intermediate coverage breaks at r={expected}")
        expected = hi + 1
        count = entropy_count_bound_arb(k, lo, hi)
        tail = uniform_bernoulli_block_tail_arb(
            n=n,
            cutoff=cutoff,
            row_weight=row_weight,
            memory=memory,
            lo=lo,
            hi=hi,
            output_marker=marker,
        )
        intermediate += count * tail * condition_probability ** (-hi)
    if expected != dense_start:
        raise ValueError("intermediate blocks do not reach dense_start")

    dense = dense_l2_bound_arb(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=row_weight,
        start=dense_start,
    )
    total = small.small_weight_bound + intermediate + dense
    return ExactRowECVerificationResult(
        success=bool(total < arb(2) ** (-target_bits)),
        total_bound=total,
        security_bits=-total.log() / arb(2).log(),
        small_weight_bound=small.small_weight_bound,
        intermediate_bound=intermediate,
        dense_bound=dense,
        largest_small_r=small.largest_r,
        largest_small_term=small.largest_term,
        slice_count=small.slice_count,
        intermediate_blocks=len(blocks),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-log2", type=int, required=True)
    parser.add_argument("--rate-denominator", type=int, default=5)
    parser.add_argument("--relative-cutoff", type=float, default=0.05)
    parser.add_argument("--row-weight", type=int, required=True)
    parser.add_argument("--memory", type=int, default=21)
    parser.add_argument("--target-bits", type=int, default=20)
    parser.add_argument("--exact-limit", type=int, default=59)
    parser.add_argument("--dense-fraction-denominator", type=int, default=4)
    parser.add_argument("--precision-bits", type=int, default=192)
    args = parser.parse_args()

    k = 1 << args.k_log2
    n = args.rate_denominator * k
    cutoff = math.floor(args.relative_cutoff * n)
    result = verify_exact_row_ec_parameters(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=args.row_weight,
        memory=args.memory,
        target_bits=args.target_bits,
        exact_limit=args.exact_limit,
        dense_start=k // args.dense_fraction_denominator,
        precision_bits=args.precision_bits,
    )
    print(f"success: {result.success}")
    print(f"total bound: {result.total_bound}")
    print(f"security bits: {result.security_bits}")
    print(f"small-weight bound: {result.small_weight_bound}")
    print(f"intermediate bound: {result.intermediate_bound}")
    print(f"dense bound: {result.dense_bound}")
    print(f"verified convolution shells: {result.slice_count}")
    print(
        f"largest small-weight term: r={result.largest_small_r}, "
        f"value={result.largest_small_term}"
    )
    print(f"intermediate blocks: {result.intermediate_blocks}")
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
