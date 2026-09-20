#!/usr/bin/env python3
"""Audit compact exact-regular expand--convolute heuristics.

This tool studies the expander before convolution.  The convolution matrix and
the per-coordinate labels are invertible, so neither can repair a rank defect
in the expander.  The tool compares three exact-regular constructions:

* unsigned stripes, which have an explicit constant-slot kernel;
* signed stripes, the cache-friendly odd-characteristic implementation; and
* independent full-index permutations, the characteristic-two fallback.

The sampled matrices are reduced-size ensemble models.  They do not reproduce
libOTe's public-seed generator bit for bit and they do not certify distance.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class AuditRow:
    region_size: int
    message_size: int
    code_size: int
    unsigned_rank: int
    signed_min_rank: int
    signed_rank_failures: int
    permutation_min_rank: int
    permutation_rank_failures: int


def modular_rank(matrix: np.ndarray, prime: int) -> int:
    """Return matrix rank over F_prime."""
    work = np.asarray(matrix, dtype=np.int64).copy() % prime
    rows, columns = work.shape
    rank = 0
    for column in range(columns):
        pivots = np.flatnonzero(work[rank:, column])
        if not len(pivots):
            continue
        pivot = rank + int(pivots[0])
        if pivot != rank:
            work[[rank, pivot]] = work[[pivot, rank]]
        work[rank] = (
            work[rank] * pow(int(work[rank, column]), prime - 2, prime)
        ) % prime
        active = np.flatnonzero(work[rank + 1 :, column]) + rank + 1
        if len(active):
            factors = work[active, column].copy()
            work[active] = (
                work[active] - factors[:, None] * work[rank]
            ) % prime
        rank += 1
        if rank == rows:
            break
    return rank


def striped_expander(
    *,
    left_degree: int,
    right_degree: int,
    region_size: int,
    seed: int,
    signed: bool,
) -> np.ndarray:
    """Return the row generator of a striped exact-regular expander."""
    if not (1 <= right_degree < left_degree and region_size >= 1):
        raise ValueError("invalid regular-expander dimensions")
    message_size = right_degree * region_size
    code_size = left_degree * region_size
    matrix = np.zeros((message_size, code_size), dtype=np.int64)
    generator = np.random.Generator(np.random.PCG64(seed))
    offsets = generator.integers(
        0,
        region_size,
        size=(left_degree - 1, right_degree),
        dtype=np.int64,
    )
    signs = generator.integers(
        0,
        2,
        size=(left_degree - 1, right_degree),
        dtype=np.int8,
    )

    for base in range(region_size):
        for slot in range(right_degree):
            left = base * right_degree + slot
            matrix[left, base] = 1
            for region in range(1, left_degree):
                right = (base + int(offsets[region - 1, slot])) % region_size
                coefficient = 1
                if signed:
                    if region < right_degree:
                        negative = slot == region
                    else:
                        negative = bool(signs[region - 1, slot])
                    coefficient = -1 if negative else 1
                matrix[left, region * region_size + right] = coefficient
    return matrix


def permutation_expander(
    *,
    left_degree: int,
    right_degree: int,
    region_size: int,
    seed: int,
) -> np.ndarray:
    """Return the row generator using one full left permutation per region."""
    if not (1 <= right_degree < left_degree and region_size >= 1):
        raise ValueError("invalid regular-expander dimensions")
    message_size = right_degree * region_size
    code_size = left_degree * region_size
    matrix = np.zeros((message_size, code_size), dtype=np.int64)
    generator = np.random.Generator(np.random.PCG64(seed))
    for left in range(message_size):
        matrix[left, left // right_degree] = 1
    for region in range(1, left_degree):
        permutation = generator.permutation(message_size)
        for left, value in enumerate(permutation):
            right = int(value) // right_degree
            matrix[left, region * region_size + right] = 1
    return matrix


def constant_slot_witness(
    *, right_degree: int, region_size: int, prime: int
) -> np.ndarray:
    """Return a nonzero message killed by every unsigned striped expander."""
    if right_degree < 2:
        raise ValueError("the witness needs at least two slots")
    message = np.zeros(right_degree * region_size, dtype=np.int64)
    for base in range(region_size):
        message[base * right_degree] = 1
        message[base * right_degree + 1] = prime - 1
    return message


def audit(
    *,
    left_degree: int,
    right_degree: int,
    region_sizes: list[int],
    prime: int,
    seeds: int,
) -> list[AuditRow]:
    if prime < 3:
        raise ValueError("signed stripes require an odd prime")
    if seeds < 1:
        raise ValueError("seed count must be positive")

    rows: list[AuditRow] = []
    for region_size in region_sizes:
        unsigned = striped_expander(
            left_degree=left_degree,
            right_degree=right_degree,
            region_size=region_size,
            seed=0,
            signed=False,
        )
        witness = constant_slot_witness(
            right_degree=right_degree,
            region_size=region_size,
            prime=prime,
        )
        if np.any((witness @ unsigned) % prime):
            raise AssertionError("the unsigned constant-slot witness did not vanish")
        unsigned_rank = modular_rank(unsigned, prime)

        signed_ranks: list[int] = []
        permutation_ranks: list[int] = []
        for seed in range(seeds):
            signed_ranks.append(
                modular_rank(
                    striped_expander(
                        left_degree=left_degree,
                        right_degree=right_degree,
                        region_size=region_size,
                        seed=seed,
                        signed=True,
                    ),
                    prime,
                )
            )
            permutation_ranks.append(
                modular_rank(
                    permutation_expander(
                        left_degree=left_degree,
                        right_degree=right_degree,
                        region_size=region_size,
                        seed=seed,
                    ),
                    prime,
                )
            )

        message_size = right_degree * region_size
        rows.append(
            AuditRow(
                region_size=region_size,
                message_size=message_size,
                code_size=left_degree * region_size,
                unsigned_rank=unsigned_rank,
                signed_min_rank=min(signed_ranks),
                signed_rank_failures=sum(
                    rank != message_size for rank in signed_ranks
                ),
                permutation_min_rank=min(permutation_ranks),
                permutation_rank_failures=sum(
                    rank != message_size for rank in permutation_ranks
                ),
            )
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-degree", type=int, default=26)
    parser.add_argument("--right-degree", type=int, default=13)
    parser.add_argument(
        "--region-sizes",
        type=int,
        nargs="+",
        default=[1, 2, 3, 5, 7, 11],
    )
    parser.add_argument("--prime", type=int, default=127)
    parser.add_argument("--seeds", type=int, default=32)
    args = parser.parse_args()
    result = audit(
        left_degree=args.left_degree,
        right_degree=args.right_degree,
        region_sizes=args.region_sizes,
        prime=args.prime,
        seeds=args.seeds,
    )
    print(
        json.dumps(
            {
                "model": "reduced-size ensemble; not libOTe seed replay",
                "left_degree": args.left_degree,
                "right_degree": args.right_degree,
                "prime": args.prime,
                "seeds": args.seeds,
                "rows": [asdict(row) for row in result],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
