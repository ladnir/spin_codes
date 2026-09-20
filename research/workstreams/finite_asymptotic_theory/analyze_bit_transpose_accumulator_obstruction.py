#!/usr/bin/env python3
"""Evaluate the pure bit-transpose accumulator obstruction.

The calculation applies only to independent random half-rate block subspaces
and a deterministic transpose with persistent row adjacency.
"""

from __future__ import annotations

import argparse
import json

import mpmath as mp


def trivial_intersection_probability(half_dimension: int) -> mp.mpf:
    numerator = mp.mpf(1)
    denominator = mp.mpf(1)
    for index in range(1, half_dimension + 1):
        numerator *= 1 - mp.power(2, -index)
    for index in range(half_dimension + 1, 2 * half_dimension + 1):
        denominator *= 1 - mp.power(2, -index)
    return numerator / denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--block-size", type=int, default=32)
    parser.add_argument("--outer-blocks", type=int, default=32768)
    args = parser.parse_args()

    if args.block_size <= 0 or args.block_size % 2:
        raise SystemExit("block size must be a positive even integer")
    if args.outer_blocks <= 0:
        raise SystemExit("outer block count must be positive")

    mp.mp.dps = 80
    half_dimension = args.block_size // 2
    pair_probability = trivial_intersection_probability(half_dimension)
    disjoint_pairs = args.outer_blocks // 2
    log2_no_collision = disjoint_pairs * mp.log(pair_probability, 2)
    total_length = args.block_size * args.outer_blocks

    payload = {
        "schema": "bit-transpose-accumulator-obstruction-v1",
        "scope": "pure deterministic transpose; no region shuffles",
        "parameters": {
            "block_size": args.block_size,
            "half_dimension": half_dimension,
            "outer_blocks": args.outer_blocks,
            "total_length": total_length,
        },
        "exact_consequences": {
            "trivial_intersection_probability_per_disjoint_pair": str(
                pair_probability
            ),
            "log2_upper_bound_on_no_collision_probability": str(
                log2_no_collision
            ),
            "bad_codeword_weight_upper_bound": args.block_size,
            "bad_codeword_relative_weight_upper_bound": (
                args.block_size / total_length
            ),
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
