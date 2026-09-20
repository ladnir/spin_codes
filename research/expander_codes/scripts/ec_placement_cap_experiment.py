#!/usr/bin/env python3
"""Exhaustive small-instance audit of placement-level EC union-bound slack.

The certified EC trace caps each output trace separately and then averages over
the random regular-expander placement.  For one fixed placement, however, all
low-weight traces describe alternatives for the same bad-codeword event.  This
script measures the possible gain from first taking their union bound and then
capping that conditional probability at one.

Only tiny instances are practical.  The calculation is a diagnostic, not a
replacement for the large-instance certificate.
"""

from __future__ import annotations

import argparse
import itertools
import math
from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class PlacementCapResult:
    uncapped_mean: float
    capped_mean: float
    improvement_bits: float
    saturated_probability: float
    placement_count: int


def regional_pattern_distribution(
    *, region_length: int, group_size: int, message_weight: int
) -> dict[tuple[int, ...], float]:
    """Distribution of occupied groups from a uniform support-slot subset."""
    slot_count = region_length * group_size
    if not 0 <= message_weight <= slot_count:
        raise ValueError("message_weight exceeds the regional slot count")
    counts: Counter[tuple[int, ...]] = Counter()
    for slots in itertools.combinations(range(slot_count), message_weight):
        occupied = {slot // group_size for slot in slots}
        pattern = tuple(int(index in occupied) for index in range(region_length))
        counts[pattern] += 1
    denominator = math.comb(slot_count, message_weight)
    return {pattern: count / denominator for pattern, count in counts.items()}


def trace_counts_by_equations(
    *, occupancy: tuple[int, ...], memory: int, cutoff: int
) -> list[int]:
    """Count low-weight trace paths by their number of fresh equations."""
    if memory < 1 or cutoff < 0:
        raise ValueError("memory must be positive and cutoff nonnegative")
    length = len(occupancy)
    # (state, output weight, equation count) -> path multiplicity.
    current = {(memory, 0, 0): 1}
    for is_occupied in occupancy:
        following: dict[tuple[int, int, int], int] = {}
        for (state, output_weight, equations), multiplicity in current.items():
            if state < memory:
                if output_weight < cutoff:
                    key = (0, output_weight + 1, equations)
                    following[key] = following.get(key, 0) + multiplicity
                key = (state + 1, output_weight, equations + 1)
                following[key] = following.get(key, 0) + multiplicity
            elif is_occupied:
                if output_weight < cutoff:
                    key = (0, output_weight + 1, equations)
                    following[key] = following.get(key, 0) + multiplicity
                key = (memory, output_weight, equations + 1)
                following[key] = following.get(key, 0) + multiplicity
            else:
                key = (memory, output_weight, equations)
                following[key] = following.get(key, 0) + multiplicity
        current = following
    counts = [0] * (length + 1)
    for (_, _, equations), multiplicity in current.items():
        counts[equations] += multiplicity
    return counts


def projective_trace_union(
    *, equation_counts: list[int], prime: int, message_weight: int
) -> float:
    """Sum the per-trace projective caps for one fixed placement."""
    if prime < 2 or message_weight < 1:
        raise ValueError("invalid projective parameters")
    scale = prime - 1
    return math.fsum(
        multiplicity * min(1.0, scale ** (message_weight - 1 - equations))
        for equations, multiplicity in enumerate(equation_counts)
    )


def exhaustive_placement_cap(
    *, prime: int, region_count: int, region_length: int, group_size: int,
    memory: int, message_weight: int, cutoff: int,
) -> PlacementCapResult:
    """Average before and after the cap conditional on expander placement."""
    regional = regional_pattern_distribution(
        region_length=region_length,
        group_size=group_size,
        message_weight=message_weight,
    )
    items = tuple(regional.items())
    uncapped = 0.0
    capped = 0.0
    saturated_probability = 0.0
    placement_count = 0
    for regions in itertools.product(items, repeat=region_count):
        occupancy = tuple(bit for (pattern, _) in regions for bit in pattern)
        probability = math.prod(probability for _, probability in regions)
        conditional_union = projective_trace_union(
            equation_counts=trace_counts_by_equations(
                occupancy=occupancy,
                memory=memory,
                cutoff=cutoff,
            ),
            prime=prime,
            message_weight=message_weight,
        )
        uncapped += probability * conditional_union
        capped += probability * min(1.0, conditional_union)
        if conditional_union >= 1.0:
            saturated_probability += probability
        placement_count += 1
    improvement = (
        math.log2(uncapped / capped) if capped > 0.0 else math.inf
    )
    return PlacementCapResult(
        uncapped_mean=uncapped,
        capped_mean=capped,
        improvement_bits=improvement,
        saturated_probability=saturated_probability,
        placement_count=placement_count,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, default=127)
    parser.add_argument("--degree", type=int, default=4)
    parser.add_argument("--region-length", type=int, default=4)
    parser.add_argument("--memory", type=int, default=3)
    parser.add_argument("--message-weight", type=int, default=2)
    parser.add_argument("--cutoff", type=int, default=7)
    args = parser.parse_args()
    if args.degree % 2:
        raise ValueError("rate-half experiment requires even degree")
    result = exhaustive_placement_cap(
        prime=args.prime,
        region_count=args.degree,
        region_length=args.region_length,
        group_size=args.degree // 2,
        memory=args.memory,
        message_weight=args.message_weight,
        cutoff=args.cutoff,
    )
    print(result)


if __name__ == "__main__":
    main()
