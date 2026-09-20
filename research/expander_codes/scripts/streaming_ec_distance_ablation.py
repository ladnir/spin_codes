#!/usr/bin/env python3
"""Reduced-distance ablations for streaming field Expand--Convolute.

The experiments vary topology, edge signs, column labels, and convolution
taps.  They materialize reduced generator matrices over a small prime field
and apply deterministic low-weight searches.  The output is experimental
evidence, not a distance certificate for the full Goldilocks instance.
"""

from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np

from streaming_ec_heuristic_audit import modular_rank


TOPOLOGY_TAG = 0x544F504F4C4F4759
EDGE_LABEL_TAG = 0x454447454C414245
COLUMN_LABEL_TAG = 0x434F4C554D4E4C42
TAP_TAG = 0x434F4E5654415053
EDGE_SIGN_TAG = 0x454447455349474E
MASK64 = (1 << 64) - 1


@dataclass(frozen=True)
class Variant:
    name: str
    topology: str
    labels: str
    taps: str


@dataclass(frozen=True)
class SearchResult:
    variant: str
    rank: int
    row_weight: int
    pair_weight: int
    prefix_weight: int
    suffix_weight: int
    interval_weight: int
    interval_start: int
    short_interval_weight: int | None
    short_interval_start: int | None
    phase_weight: int
    phase_count: int


VARIANTS = (
    Variant("proved_like", "regular", "edge", "independent"),
    Variant("striped_topology", "striped", "edge", "independent"),
    Variant("on_the_fly_heuristic", "striped", "column", "independent"),
    Variant("periodic_tap_control", "striped", "column", "periodic"),
)

COLLISION_FREE_VARIANTS = (
    Variant(
        "collision_free_topology",
        "striped_collision_free",
        "edge",
        "independent",
    ),
    Variant(
        "collision_free_heuristic",
        "striped_collision_free",
        "column",
        "independent",
    ),
    Variant(
        "collision_free_edge_signs",
        "striped_collision_free_edge_signs",
        "column",
        "independent",
    ),
    Variant(
        "collision_free_splitmix_labels",
        "striped_collision_free",
        "column_splitmix",
        "independent",
    ),
    Variant(
        "collision_free_affine_labels",
        "striped_collision_free",
        "column_affine",
        "independent",
    ),
    Variant(
        "collision_free_unit_labels",
        "striped_collision_free",
        "column_unit",
        "independent",
    ),
    Variant(
        "collision_free_sign_labels",
        "striped_collision_free",
        "column_sign",
        "independent",
    ),
    Variant(
        "collision_free_periodic",
        "striped_collision_free",
        "column",
        "periodic",
    ),
)

LABEL_STRESS_VARIANTS = (
    Variant(
        "splitmix_label_control",
        "striped",
        "column_splitmix",
        "independent",
    ),
    Variant(
        "affine_label_control",
        "striped",
        "column_affine",
        "independent",
    ),
    Variant("unit_label_control", "striped", "column_unit", "independent"),
    Variant("sign_label_control", "striped", "column_sign", "independent"),
)

ALL_VARIANTS = VARIANTS + COLLISION_FREE_VARIANTS + LABEL_STRESS_VARIANTS


def rng(seed: int, tag: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, tag])))


def nonzero_samples(
    generator: np.random.Generator, shape: tuple[int, ...], prime: int
) -> np.ndarray:
    return generator.integers(1, prime, size=shape, dtype=np.int64)


def mix64(value: int) -> int:
    value &= MASK64
    value ^= value >> 30
    value = (value * 0xBF58476D1CE4E5B9) & MASK64
    value ^= value >> 27
    value = (value * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def structured_column_labels(
    *, mode: str, size: int, prime: int, seed: int
) -> np.ndarray:
    if mode == "column_unit":
        return np.ones(size, dtype=np.int64)

    material = mix64(seed + COLUMN_LABEL_TAG)
    if mode == "column_sign":
        return np.fromiter(
            (
                prime - 1 if mix64(
                    material + index * 0x9E3779B97F4A7C15
                ) >> 63 else 1
                for index in range(size)
            ),
            dtype=np.int64,
            count=size,
        )
    if mode == "column_affine":
        multiplier = mix64(material + 0x9E3779B97F4A7C15) | 1
        addend = mix64(material + 0xD1B54A32D192ED03)
        words = (
            (multiplier * index + addend) & MASK64
            for index in range(size)
        )
    elif mode == "column_splitmix":
        words = (
            mix64(material + index * 0x9E3779B97F4A7C15)
            for index in range(size)
        )
    else:
        raise ValueError(f"unknown structured label mode: {mode}")

    labels = np.fromiter(
        (word % prime for word in words), dtype=np.int64, count=size
    )
    labels[labels == 0] = 1
    return labels


def regular_neighbors(
    *, left_degree: int, right_degree: int, region_size: int, seed: int
) -> np.ndarray:
    message_size = right_degree * region_size
    neighbors = np.empty((message_size, left_degree), dtype=np.int64)
    neighbors[:, 0] = np.arange(message_size, dtype=np.int64) // right_degree
    generator = rng(seed, TOPOLOGY_TAG)
    for region in range(1, left_degree):
        permutation = generator.permutation(message_size)
        neighbors[:, region] = (
            region * region_size + permutation // right_degree
        )
    return neighbors


def striped_neighbors_and_signs(
    *,
    left_degree: int,
    right_degree: int,
    region_size: int,
    seed: int,
    collision_free: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    message_size = right_degree * region_size
    neighbors = np.empty((message_size, left_degree), dtype=np.int64)
    signs = np.ones((message_size, left_degree), dtype=np.int64)
    generator = rng(seed, TOPOLOGY_TAG)

    def is_prime(value: int) -> bool:
        if value < 2:
            return False
        divisor = 2
        while divisor * divisor <= value:
            if value % divisor == 0:
                return False
            divisor += 1
        return True

    def differences_are_distinct(offsets: np.ndarray) -> bool:
        extended = np.vstack((
            np.zeros((1, right_degree), dtype=np.int64), offsets
        ))
        for first in range(right_degree):
            for second in range(first + 1, right_degree):
                differences = (
                    extended[:, first] - extended[:, second]
                ) % region_size
                if len(np.unique(differences)) != left_degree:
                    return False
        return True

    if (
        collision_free
        and is_prime(region_size)
        and left_degree <= region_size
        and right_degree <= region_size
    ):
        slopes = generator.choice(
            np.arange(1, region_size, dtype=np.int64),
            size=left_degree - 1,
            replace=False,
        )
        slots = generator.choice(
            np.arange(region_size, dtype=np.int64),
            size=right_degree,
            replace=False,
        )
        intercepts = generator.integers(
            0,
            region_size,
            size=(left_degree - 1, 1),
            dtype=np.int64,
        )
        offsets = (slopes[:, None] * slots[None, :] + intercepts) % region_size
        random_signs = generator.integers(
            0,
            2,
            size=(left_degree - 1, right_degree),
            dtype=np.int8,
        )
    else:
        for _ in range(10_000):
            offsets = generator.integers(
                0,
                region_size,
                size=(left_degree - 1, right_degree),
                dtype=np.int64,
            )
            random_signs = generator.integers(
                0,
                2,
                size=(left_degree - 1, right_degree),
                dtype=np.int8,
            )
            if not collision_free or differences_are_distinct(offsets):
                break
        else:
            raise ValueError("could not sample collision-free striped offsets")
    for base in range(region_size):
        for slot in range(right_degree):
            left = base * right_degree + slot
            neighbors[left, 0] = base
            for region in range(1, left_degree):
                right = (base + int(offsets[region - 1, slot])) % region_size
                neighbors[left, region] = region * region_size + right
                if region < right_degree:
                    negative = slot == region
                else:
                    negative = bool(random_signs[region - 1, slot])
                signs[left, region] = -1 if negative else 1
    return neighbors, signs


def expander_generator(
    *,
    variant: Variant,
    left_degree: int,
    right_degree: int,
    region_size: int,
    prime: int,
    seed: int,
) -> np.ndarray:
    message_size = right_degree * region_size
    code_size = left_degree * region_size
    if variant.topology == "regular":
        neighbors = regular_neighbors(
            left_degree=left_degree,
            right_degree=right_degree,
            region_size=region_size,
            seed=seed,
        )
        signs = np.ones_like(neighbors)
    elif variant.topology in (
        "striped",
        "striped_collision_free",
        "striped_collision_free_edge_signs",
    ):
        neighbors, signs = striped_neighbors_and_signs(
            left_degree=left_degree,
            right_degree=right_degree,
            region_size=region_size,
            seed=seed,
            collision_free=variant.topology != "striped",
        )
        if variant.topology == "striped_collision_free_edge_signs":
            edge_sign_bits = rng(seed, EDGE_SIGN_TAG).integers(
                0,
                2,
                size=(message_size, left_degree - 1),
                dtype=np.int8,
            )
            signs[:, 1:] = 1 - 2 * edge_sign_bits
    else:
        raise ValueError(f"unknown topology: {variant.topology}")

    matrix = np.zeros((message_size, code_size), dtype=np.int64)
    if variant.labels == "edge":
        labels = nonzero_samples(
            rng(seed, EDGE_LABEL_TAG), neighbors.shape, prime
        )
    elif variant.labels == "column":
        column_labels = nonzero_samples(
            rng(seed, COLUMN_LABEL_TAG), (code_size,), prime
        )
        labels = column_labels[neighbors]
    elif variant.labels in (
        "column_splitmix",
        "column_affine",
        "column_sign",
        "column_unit",
    ):
        column_labels = structured_column_labels(
            mode=variant.labels,
            size=code_size,
            prime=prime,
            seed=seed,
        )
        labels = column_labels[neighbors]
    else:
        raise ValueError(f"unknown label mode: {variant.labels}")

    row_indices = np.arange(message_size)
    for region in range(left_degree):
        matrix[row_indices, neighbors[:, region]] = (
            signs[:, region] * labels[:, region]
        ) % prime
    return matrix


def convolution_taps(
    *,
    code_size: int,
    memory: int,
    period: int,
    prime: int,
    seed: int,
    mode: str,
) -> np.ndarray:
    generator = rng(seed, TAP_TAG)
    if mode == "independent":
        return generator.integers(
            0, prime, size=(code_size, memory), dtype=np.int64
        )
    if mode == "periodic":
        base = generator.integers(
            0, prime, size=(period, memory), dtype=np.int64
        )
        return base[np.arange(code_size) % period]
    raise ValueError(f"unknown tap mode: {mode}")


def convolve_generator(
    matrix: np.ndarray, taps: np.ndarray, prime: int
) -> np.ndarray:
    result = matrix.copy()
    memory = taps.shape[1]
    for time in range(result.shape[1]):
        for lag in range(1, min(memory, time) + 1):
            result[:, time] = (
                result[:, time] + taps[time, lag - 1] * result[:, time - lag]
            ) % prime
    return result


def null_vector(matrix: np.ndarray, prime: int) -> np.ndarray | None:
    """Return nonzero x with matrix @ x = 0, or None if the kernel is zero."""
    work = np.asarray(matrix, dtype=np.int64).copy() % prime
    row_count, column_count = work.shape
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        candidates = np.flatnonzero(work[pivot_row:, column])
        if not len(candidates):
            continue
        selected = pivot_row + int(candidates[0])
        if selected != pivot_row:
            work[[pivot_row, selected]] = work[[selected, pivot_row]]
        work[pivot_row] = (
            work[pivot_row]
            * pow(int(work[pivot_row, column]), prime - 2, prime)
        ) % prime
        active = np.flatnonzero(work[:, column])
        active = active[active != pivot_row]
        if len(active):
            factors = work[active, column].copy()
            work[active] = (
                work[active] - factors[:, None] * work[pivot_row]
            ) % prime
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break

    if len(pivot_columns) == column_count:
        return None
    pivot_set = set(pivot_columns)
    free = next(
        column for column in range(column_count) if column not in pivot_set
    )
    vector = np.zeros(column_count, dtype=np.int64)
    vector[free] = 1
    for row, pivot in reversed(list(enumerate(pivot_columns))):
        vector[pivot] = -int(work[row, free]) % prime
    if np.any((matrix @ vector) % prime):
        raise AssertionError("nullspace solver returned an invalid vector")
    return vector


def supported_word(
    generator: np.ndarray, support: np.ndarray, prime: int
) -> tuple[int, np.ndarray] | None:
    outside = np.ones(generator.shape[1], dtype=bool)
    outside[support] = False
    message = null_vector(generator[:, outside].T, prime)
    if message is None:
        return None
    word = (message @ generator) % prime
    if np.any(word[outside]):
        raise AssertionError("shortening witness escaped its support")
    return int(np.count_nonzero(word)), message


def projective_pair_weight(generator: np.ndarray, prime: int) -> int:
    inverses = np.zeros(prime, dtype=np.int64)
    for value in range(1, prime):
        inverses[value] = pow(value, prime - 2, prime)
    best = int(np.min(np.count_nonzero(generator, axis=1)))
    for first, second in itertools.combinations(range(generator.shape[0]), 2):
        lhs = generator[first]
        rhs = generator[second]
        lhs_nonzero = lhs != 0
        rhs_nonzero = rhs != 0
        both = lhs_nonzero & rhs_nonzero
        union_weight = int(np.count_nonzero(lhs_nonzero | rhs_nonzero))
        if np.any(both):
            ratios = (-lhs[both] * inverses[rhs[both]]) % prime
            cancellation = int(np.max(np.bincount(ratios, minlength=prime)))
            best = min(best, union_weight - cancellation)
    return best


def interval_starts(code_size: int, region_size: int, period: int) -> list[int]:
    starts = {0, code_size - (code_size // 2 + 1)}
    starts.update(range(0, code_size, region_size))
    # Keep the diagnostic bounded when a deliberately tiny tap period is used.
    tap_step = max(period, max(1, code_size // 32))
    starts.update(range(0, code_size, tap_step))
    return sorted(start for start in starts if 0 <= start < code_size)


def best_interval(
    generator: np.ndarray,
    *,
    length: int,
    starts: Iterable[int],
    prime: int,
) -> tuple[int, int] | None:
    best: tuple[int, int] | None = None
    code_size = generator.shape[1]
    for start in starts:
        if start + length > code_size:
            continue
        support = np.arange(start, start + length)
        witness = supported_word(generator, support, prime)
        if witness is None:
            continue
        candidate = (witness[0], start)
        if best is None or candidate < best:
            best = candidate
    return best


def best_phase_support(
    generator: np.ndarray, *, period: int, prime: int
) -> tuple[int, int]:
    code_size = generator.shape[1]
    lo, hi = 1, period
    best: tuple[int, int] | None = None
    while lo <= hi:
        count = (lo + hi) // 2
        support = np.flatnonzero((np.arange(code_size) % period) < count)
        witness = supported_word(generator, support, prime)
        if witness is None:
            lo = count + 1
        else:
            best = (witness[0], count)
            hi = count - 1
    if best is None:
        raise AssertionError("all tap phases did not support a codeword")
    return best


def run_variant(
    *,
    variant: Variant,
    left_degree: int,
    right_degree: int,
    region_size: int,
    memory: int,
    period: int,
    prime: int,
    seed: int,
) -> SearchResult:
    expander = expander_generator(
        variant=variant,
        left_degree=left_degree,
        right_degree=right_degree,
        region_size=region_size,
        prime=prime,
        seed=seed,
    )
    taps = convolution_taps(
        code_size=expander.shape[1],
        memory=memory,
        period=period,
        prime=prime,
        seed=seed,
        mode=variant.taps,
    )
    generator = convolve_generator(expander, taps, prime)
    rank = modular_rank(generator, prime)
    if rank != generator.shape[0]:
        return SearchResult(
            variant=variant.name,
            rank=rank,
            row_weight=0,
            pair_weight=0,
            prefix_weight=0,
            suffix_weight=0,
            interval_weight=0,
            interval_start=0,
            short_interval_weight=0,
            short_interval_start=0,
            phase_weight=0,
            phase_count=0,
        )

    message_size, code_size = generator.shape
    singleton_length = code_size - message_size + 1
    prefix = supported_word(
        generator, np.arange(singleton_length), prime
    )
    suffix = supported_word(
        generator,
        np.arange(code_size - singleton_length, code_size),
        prime,
    )
    if prefix is None or suffix is None:
        raise AssertionError("Singleton-length interval had no codeword")

    starts = interval_starts(code_size, region_size, period)
    interval = best_interval(
        generator,
        length=singleton_length,
        starts=starts,
        prime=prime,
    )
    if interval is None:
        raise AssertionError("no Singleton-length interval witness")
    short = best_interval(
        generator,
        length=singleton_length - 1,
        starts=starts,
        prime=prime,
    )
    phase_weight, phase_count = best_phase_support(
        generator, period=period, prime=prime
    )

    return SearchResult(
        variant=variant.name,
        rank=rank,
        row_weight=int(np.min(np.count_nonzero(generator, axis=1))),
        pair_weight=projective_pair_weight(generator, prime),
        prefix_weight=prefix[0],
        suffix_weight=suffix[0],
        interval_weight=interval[0],
        interval_start=interval[1],
        short_interval_weight=None if short is None else short[0],
        short_interval_start=None if short is None else short[1],
        phase_weight=phase_weight,
        phase_count=phase_count,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-degree", type=int, default=26)
    parser.add_argument("--right-degree", type=int, default=13)
    parser.add_argument("--region-size", type=int, default=21)
    parser.add_argument("--memory", type=int, default=4)
    parser.add_argument("--period", type=int, default=256)
    parser.add_argument("--prime", type=int, default=127)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=[variant.name for variant in ALL_VARIANTS],
        default=[variant.name for variant in VARIANTS],
    )
    args = parser.parse_args()
    if args.period < 1:
        raise SystemExit("period must be positive")

    selected = [
        variant for variant in ALL_VARIANTS if variant.name in args.variants
    ]
    results = [
        run_variant(
            variant=variant,
            left_degree=args.left_degree,
            right_degree=args.right_degree,
            region_size=args.region_size,
            memory=args.memory,
            period=args.period,
            prime=args.prime,
            seed=args.seed,
        )
        for variant in selected
    ]
    print(
        json.dumps(
            {
                "scope": "reduced experimental search; not a distance certificate",
                "parameters": {
                    "left_degree": args.left_degree,
                    "right_degree": args.right_degree,
                    "region_size": args.region_size,
                    "message_size": args.right_degree * args.region_size,
                    "code_size": args.left_degree * args.region_size,
                    "memory": args.memory,
                    "period": args.period,
                    "prime": args.prime,
                    "seed": args.seed,
                },
                "results": [asdict(result) for result in results],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
