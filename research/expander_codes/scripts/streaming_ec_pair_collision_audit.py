#!/usr/bin/env python3
"""Audit projective pair collisions in the signed striped expander.

Two message rows can share a code coordinate only within the same region.  A
shared column label cancels from their coefficient ratio.  Thus one projective
combination cancels every shared edge with the same sign ratio.  This tool
counts the maximum such multiplicity without materializing the full expander.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter

import numpy as np

from streaming_ec_distance_ablation import TOPOLOGY_TAG, rng


def sample_offsets_and_signs(
    *, left_degree: int, right_degree: int, region_size: int, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    generator = rng(seed, TOPOLOGY_TAG)
    offsets = np.zeros((left_degree, right_degree), dtype=np.int64)
    offsets[1:] = generator.integers(
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
    signs = np.ones((left_degree, right_degree), dtype=np.int8)
    for region in range(1, left_degree):
        for slot in range(right_degree):
            if region < right_degree:
                negative = slot == region
            else:
                negative = bool(random_signs[region - 1, slot])
            signs[region, slot] = -1 if negative else 1
    return offsets, signs


def collision_profile(
    *, left_degree: int, right_degree: int, region_size: int, seed: int
) -> dict[str, int]:
    offsets, signs = sample_offsets_and_signs(
        left_degree=left_degree,
        right_degree=right_degree,
        region_size=region_size,
        seed=seed,
    )
    multiplicities: Counter[int] = Counter()
    maximum = 0
    repeated_keys = 0
    for first in range(right_degree):
        for second in range(first + 1, right_degree):
            keys: Counter[tuple[int, int]] = Counter()
            for region in range(left_degree):
                difference = int(
                    (offsets[region, first] - offsets[region, second])
                    % region_size
                )
                sign_ratio = int(signs[region, first] * signs[region, second])
                keys[(difference, sign_ratio)] += 1
            for count in keys.values():
                multiplicities[count] += 1
                maximum = max(maximum, count)
                repeated_keys += count > 1
    return {
        "maximum_cancellable_shared_edges": maximum,
        "repeated_projective_keys": repeated_keys,
        "keys_of_multiplicity_2": multiplicities[2],
        "keys_of_multiplicity_3": multiplicities[3],
        "keys_of_multiplicity_4_or_more": sum(
            count for multiplicity, count in multiplicities.items()
            if multiplicity >= 4
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-degree", type=int, default=26)
    parser.add_argument("--right-degree", type=int, default=13)
    parser.add_argument("--region-size", type=int, default=80659)
    parser.add_argument("--seeds", type=int, default=1000)
    args = parser.parse_args()
    profiles = [
        collision_profile(
            left_degree=args.left_degree,
            right_degree=args.right_degree,
            region_size=args.region_size,
            seed=seed,
        )
        for seed in range(args.seeds)
    ]
    maxima = Counter(
        profile["maximum_cancellable_shared_edges"] for profile in profiles
    )
    print(
        json.dumps(
            {
                "left_degree": args.left_degree,
                "right_degree": args.right_degree,
                "region_size": args.region_size,
                "seeds": args.seeds,
                "maximum_histogram": {
                    str(key): value for key, value in sorted(maxima.items())
                },
                "seeds_with_repeated_projective_key": sum(
                    profile["repeated_projective_keys"] > 0
                    for profile in profiles
                ),
                "total_multiplicity_2_keys": sum(
                    profile["keys_of_multiplicity_2"] for profile in profiles
                ),
                "total_multiplicity_3_keys": sum(
                    profile["keys_of_multiplicity_3"] for profile in profiles
                ),
                "total_multiplicity_4_or_more_keys": sum(
                    profile["keys_of_multiplicity_4_or_more"]
                    for profile in profiles
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
