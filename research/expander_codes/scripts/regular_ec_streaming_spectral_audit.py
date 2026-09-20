#!/usr/bin/env python3
"""Audit the exact libOTe signed streaming topology at cyclic frequencies.

The script replays the public-seed topology sampler used by
RegularEcStreamingFieldCode.  It checks all complex roots of X^q-1
numerically.  It also checks the order-1, order-79, and order-1021 components
exactly over the Goldilocks field.  The primitive order-80659 component is
reported only by the numerical audit because its six Goldilocks factors have
degree 13260.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass

import numpy as np
from flint import fmpz_mod_poly_ctx, fq_default_ctx, nmod_poly


MASK64 = (1 << 64) - 1
PERMUTATION_DOMAIN = 0x177DEE2D21EA3BC5
ATTEMPT_STEP = 0xD1B54A32D192ED03
INDEX_STEP = 0x9E3779B97F4A7C15
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
GOLDILOCKS = (1 << 64) - (1 << 32) + 1
AUXILIARY_PRIME = 2_150_530_259
DEPLOYED_REGION_SIZE = 80_659


@dataclass(frozen=True)
class Topology:
    offsets: np.ndarray
    signs: np.ndarray
    attempts: int


def rotl64(value: int, shift: int) -> int:
    return ((value << shift) | (value >> (64 - shift))) & MASK64


def mix64(value: int) -> int:
    value &= MASK64
    value ^= value >> 30
    value = (value * 0xBF58476D1CE4E5B9) & MASK64
    value ^= value >> 27
    value = (value * 0x94D049BB133111EB) & MASK64
    value ^= value >> 31
    return value & MASK64


def has_distinct_differences(offsets: np.ndarray, region_size: int) -> bool:
    left_degree, right_degree = offsets.shape
    for first in range(right_degree):
        for second in range(first + 1, right_degree):
            seen: set[int] = set()
            for region in range(left_degree):
                difference = int(
                    (int(offsets[region, first]) -
                     int(offsets[region, second])) % region_size
                )
                if difference in seen:
                    return False
                seen.add(difference)
    return True


def sample_topology(
    *,
    left_degree: int,
    right_degree: int,
    region_size: int,
    seed_low: int,
    seed_high: int,
) -> Topology:
    collision_budget = (
        right_degree * (right_degree - 1) // 2
        * left_degree * (left_degree - 1) // 2
    )
    require_distinct = region_size >= 2 * collision_budget
    seed_material = seed_low ^ rotl64(seed_high, 23)

    for attempt in range(128):
        offsets = np.zeros((left_degree, right_degree), dtype=np.int64)
        signs = np.ones((left_degree, right_degree), dtype=np.int8)
        for region in range(1, left_degree):
            for slot in range(right_degree):
                index = (region - 1) * right_degree + slot
                sample = mix64(
                    seed_material
                    + PERMUTATION_DOMAIN
                    + attempt * ATTEMPT_STEP
                    + index * INDEX_STEP
                )
                offsets[region, slot] = sample % region_size
                negative = (
                    slot == region if region < right_degree
                    else bool(sample >> 63)
                )
                signs[region, slot] = -1 if negative else 1
        if not require_distinct or has_distinct_differences(
            offsets, region_size
        ):
            return Topology(offsets, signs, attempt + 1)
    raise RuntimeError("libOTe topology sampler exhausted 128 attempts")


def topology_digest(topology: Topology) -> int:
    digest = FNV_OFFSET
    for offset, sign in zip(
        topology.offsets.flat, topology.signs.flat, strict=True
    ):
        encoded = int(offset) | ((int(sign) < 0) << 32)
        digest ^= encoded
        digest = (digest * FNV_PRIME) & MASK64
    return digest


def offset_digest(topology: Topology) -> int:
    digest = FNV_OFFSET
    for offset in topology.offsets.flat:
        digest ^= int(offset)
        digest = (digest * FNV_PRIME) & MASK64
    return digest


def matrix_rank(elements: list[list[object]], zero: object, one: object) -> int:
    work = [row.copy() for row in elements]
    rows = len(work)
    columns = len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, rows) if work[row][column] != zero),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = one / work[rank][column]
        work[rank] = [value * inverse for value in work[rank]]
        for row in range(rows):
            if row == rank or work[row][column] == zero:
                continue
            factor = work[row][column]
            work[row] = [
                lhs - factor * rhs
                for lhs, rhs in zip(work[row], work[rank], strict=True)
            ]
        rank += 1
        if rank == columns:
            break
    return rank


def evaluate_rank(
    topology: Topology, context: object, root: object, order: int
) -> int:
    powers = [context.one()]
    for _ in range(1, order):
        powers.append(powers[-1] * root)
    matrix: list[list[object]] = []
    for region in range(topology.offsets.shape[0]):
        row = []
        for slot in range(topology.offsets.shape[1]):
            value = powers[int(topology.offsets[region, slot]) % order]
            row.append(value if topology.signs[region, slot] > 0 else -value)
        matrix.append(row)
    return matrix_rank(matrix, context.zero(), context.one())


def goldilocks_component_ranks(topology: Topology) -> dict[str, object]:
    ranks: dict[str, object] = {}

    integer_matrix = [
        [int(value) % GOLDILOCKS for value in row]
        for row in topology.signs.tolist()
    ]
    # Use nmod values to reuse the generic exact elimination.
    base_context = fq_default_ctx(GOLDILOCKS, 1, "b", fq_type="FQ_NMOD")
    base_matrix = [
        [base_context(value) for value in row] for row in integer_matrix
    ]
    ranks["order_1"] = matrix_rank(
        base_matrix, base_context.zero(), base_context.one()
    )

    phi79 = nmod_poly([1] * 79, GOLDILOCKS)
    factors79 = phi79.factor()[1]
    large_modulus_polynomials = fmpz_mod_poly_ctx(GOLDILOCKS)
    ranks79 = []
    for factor, exponent in factors79:
        if exponent != 1 or factor.degree() != 39:
            raise RuntimeError("unexpected factorization of Phi_79")
        modulus = large_modulus_polynomials([int(value) for value in factor])
        context = fq_default_ctx(modulus=modulus, var="r79", fq_type="FQ")
        ranks79.append(evaluate_rank(topology, context, context.gen(), 79))
    ranks["order_79"] = ranks79

    phi1021 = large_modulus_polynomials([1] * 1021)
    context1021 = fq_default_ctx(
        modulus=phi1021, var="r1021", fq_type="FQ"
    )
    ranks["order_1021"] = [
        evaluate_rank(topology, context1021, context1021.gen(), 1021)
    ]
    return ranks


def complex_frequency_audit(
    topology: Topology, region_size: int, batch_size: int
) -> dict[str, object]:
    minimum = math.inf
    minimum_frequency = -1
    below = {"1e-8": 0, "1e-6": 0, "1e-4": 0}
    offsets = topology.offsets.astype(np.float64)
    signs = topology.signs.astype(np.float64)
    scale = 2j * np.pi / region_size
    for start in range(0, region_size, batch_size):
        frequencies = np.arange(
            start, min(start + batch_size, region_size), dtype=np.float64
        )
        matrices = signs[None, :, :] * np.exp(
            scale * frequencies[:, None, None] * offsets[None, :, :]
        )
        smallest = np.linalg.svd(matrices, compute_uv=False)[:, -1]
        local = int(np.argmin(smallest))
        if float(smallest[local]) < minimum:
            minimum = float(smallest[local])
            minimum_frequency = start + local
        below["1e-8"] += int(np.count_nonzero(smallest < 1e-8))
        below["1e-6"] += int(np.count_nonzero(smallest < 1e-6))
        below["1e-4"] += int(np.count_nonzero(smallest < 1e-4))
    return {
        "frequencies_checked": region_size,
        "minimum_singular_value": minimum,
        "minimum_frequency": minimum_frequency,
        "counts_below": below,
    }


def primitive_root_of_order(order: int, prime: int) -> int:
    if (prime - 1) % order:
        raise ValueError("the requested order does not divide prime - 1")
    prime_factors: list[int] = []
    remaining = order
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            prime_factors.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        prime_factors.append(remaining)
    exponent = (prime - 1) // order
    for candidate in range(2, 1000):
        root = pow(candidate, exponent, prime)
        if root != 1 and all(
            pow(root, order // factor, prime) != 1
            for factor in prime_factors
        ):
            return root
    raise RuntimeError("could not find a root of the requested order")


def auxiliary_prime_frequency_audit(
    topology: Topology,
    region_size: int,
    batch_size: int,
    prime: int = AUXILIARY_PRIME,
) -> dict[str, object]:
    root = primitive_root_of_order(region_size, prime)
    powers = np.empty(region_size, dtype=np.int64)
    powers[0] = 1
    for exponent in range(1, region_size):
        powers[exponent] = (int(powers[exponent - 1]) * root) % prime

    minimum_rank = topology.offsets.shape[1]
    deficient = 0
    first_deficient = None
    offsets = topology.offsets.astype(np.int64)
    negative = topology.signs < 0
    for start in range(0, region_size, batch_size):
        frequencies = np.arange(
            start, min(start + batch_size, region_size), dtype=np.int64
        )
        indices = (
            frequencies[:, None, None] * offsets[None, :, :]
        ) % region_size
        matrices = powers[indices]
        matrices = np.where(negative[None, :, :], prime - matrices, matrices)
        ranks = np.zeros(len(frequencies), dtype=np.int16)
        for column in range(offsets.shape[1]):
            candidates = matrices[:, column:, column] != 0
            has_pivot = np.any(candidates, axis=1)
            pivot_rows = column + np.argmax(candidates, axis=1)
            batch_indices = np.arange(len(frequencies))
            saved = matrices[batch_indices, column].copy()
            matrices[batch_indices, column] = matrices[
                batch_indices, pivot_rows
            ]
            matrices[batch_indices, pivot_rows] = saved

            pivots = matrices[:, column, column]
            inverses = np.array(
                [
                    pow(int(value), prime - 2, prime) if present else 1
                    for value, present in zip(
                        pivots, has_pivot, strict=True
                    )
                ],
                dtype=np.int64,
            )
            matrices[:, column, column:] = (
                matrices[:, column, column:] * inverses[:, None]
            ) % prime
            for row in range(offsets.shape[0]):
                if row == column:
                    continue
                factors = matrices[:, row, column].copy()
                products = (
                    factors[:, None] * matrices[:, column, column:]
                ) % prime
                matrices[:, row, column:] = (
                    matrices[:, row, column:] - products
                ) % prime
            ranks += has_pivot

        minimum_rank = min(minimum_rank, int(ranks.min()))
        batch_deficient = ranks < offsets.shape[1]
        deficient += int(np.count_nonzero(batch_deficient))
        if first_deficient is None and np.any(batch_deficient):
            first_deficient = start + int(np.flatnonzero(batch_deficient)[0])
    return {
        "prime": prime,
        "root": root,
        "frequencies_checked": region_size,
        "minimum_rank": minimum_rank,
        "deficient_frequencies": deficient,
        "first_deficient_frequency": first_deficient,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-degree", type=int, default=26)
    parser.add_argument("--right-degree", type=int, default=13)
    parser.add_argument("--region-size", type=int, default=80659)
    parser.add_argument("--seed-low", type=int, default=0)
    parser.add_argument("--seed-high", type=int, default=MASK64)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--skip-exact", action="store_true")
    parser.add_argument("--skip-complex", action="store_true")
    parser.add_argument("--skip-auxiliary", action="store_true")
    args = parser.parse_args()

    topology = sample_topology(
        left_degree=args.left_degree,
        right_degree=args.right_degree,
        region_size=args.region_size,
        seed_low=args.seed_low & MASK64,
        seed_high=args.seed_high & MASK64,
    )
    scope = ["exact libOTe topology replay"]
    if not args.skip_exact:
        scope.append(
            "exact Goldilocks audit excluding primitive order-80659 factors"
        )
    if not args.skip_complex:
        scope.append("numerical complex-frequency audit")
    if not args.skip_auxiliary:
        scope.append("exact auxiliary-prime frequency audit")
    result: dict[str, object] = {
        "scope": "; ".join(scope),
        "parameters": {
            "left_degree": args.left_degree,
            "right_degree": args.right_degree,
            "region_size": args.region_size,
            "seed_low": args.seed_low & MASK64,
            "seed_high": args.seed_high & MASK64,
        },
        "sampling_attempts": topology.attempts,
        "distinct_slot_differences": has_distinct_differences(
            topology.offsets, args.region_size
        ),
        "topology_digest": f"0x{topology_digest(topology):016x}",
        "offset_digest": f"0x{offset_digest(topology):016x}",
    }
    if not args.skip_exact:
        if args.region_size != DEPLOYED_REGION_SIZE:
            parser.error(
                "the exact Goldilocks component audit is specialized to "
                f"region size {DEPLOYED_REGION_SIZE}; pass --skip-exact"
            )
        result["exact_goldilocks_ranks"] = goldilocks_component_ranks(topology)
    if not args.skip_complex:
        result["complex_frequency_audit"] = complex_frequency_audit(
            topology, args.region_size, args.batch_size
        )
    if not args.skip_auxiliary:
        result["auxiliary_prime_frequency_audit"] = (
            auxiliary_prime_frequency_audit(
                topology, args.region_size, args.batch_size
            )
        )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
