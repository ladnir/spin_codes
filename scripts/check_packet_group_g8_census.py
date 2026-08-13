#!/usr/bin/env python3
"""Verify the exact support census for the g=8 profile domain."""

from __future__ import annotations

import math

from packet_group_drive_stratified import feasible_profile_count
from packet_group_outer_profile import N, atom_count


GROUP_BITS = 8
CLASSES = GROUP_BITS + 1
MINIMUM_PHYSICAL_WEIGHT = 21
EXPECTED_FEASIBLE_SUPPORTS = 510
EXPECTED_DIMENSION_COUNTS = (8, 36, 84, 126, 126, 84, 36, 9, 1)
EXPECTED_EXCLUDED_BY_SUPPORT_SIZE = (1, 52, 437, 957, 582, 70, 0, 0, 0)
EXPECTED_PROFILE_COUNT = 553169839211945865258921061892182603726


def bounded_partition_count(weights: tuple[int, ...], threshold: int) -> int:
    """Count nonnegative combinations with weighted sum below threshold."""

    if threshold <= 0:
        return 0
    counts = [0] * threshold
    counts[0] = 1
    for weight in weights:
        for total in range(weight, threshold):
            counts[total] += counts[total - weight]
    return sum(counts)


def support_profile_count(mask: int) -> tuple[int, int]:
    support = tuple(index for index in range(CLASSES) if mask & (1 << index))
    size = len(support)
    if support == (0,):
        return 0, 1
    total = math.comb(atom_count(GROUP_BITS) - 1, size - 1)
    if 0 not in support:
        return total, 0
    positive = tuple(index for index in support if index)
    threshold = MINIMUM_PHYSICAL_WEIGHT - sum(positive)
    excluded = bounded_partition_count(positive, threshold)
    return total - excluded, excluded


def main() -> None:
    dimension_counts = [0] * CLASSES
    excluded_by_size = [0] * CLASSES
    support_counts = []
    for mask in range(1, 1 << CLASSES):
        count, excluded = support_profile_count(mask)
        size = mask.bit_count()
        excluded_by_size[size - 1] += excluded
        if count:
            dimension_counts[size - 1] += 1
            support_counts.append((mask, count))

    total = sum(count for _mask, count in support_counts)
    if len(support_counts) != EXPECTED_FEASIBLE_SUPPORTS:
        raise SystemExit("g=8 census: feasible support count changed")
    if tuple(dimension_counts) != EXPECTED_DIMENSION_COUNTS:
        raise SystemExit("g=8 census: dimension distribution changed")
    if tuple(excluded_by_size) != EXPECTED_EXCLUDED_BY_SUPPORT_SIZE:
        raise SystemExit("g=8 census: exclusion distribution changed")
    if total != EXPECTED_PROFILE_COUNT:
        raise SystemExit("g=8 census: support counts do not sum to expected total")
    if total != feasible_profile_count(GROUP_BITS, N, MINIMUM_PHYSICAL_WEIGHT):
        raise SystemExit("g=8 census: independent profile counts disagree")

    print(f"feasible_supports={len(support_counts)}")
    print("dimension_counts=" + ",".join(map(str, dimension_counts)))
    print("excluded_by_support_size=" + ",".join(map(str, excluded_by_size)))
    print(f"feasible_profiles={total}")
    print(f"log2_feasible_profiles={math.log2(total):.14f}")
    print("g=8 exact-support census: PASS")


if __name__ == "__main__":
    main()
