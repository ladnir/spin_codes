#!/usr/bin/env python3
"""Exact small-instance distance ablation for field Expand--Convolute."""

from __future__ import annotations

import argparse
import json
from collections import Counter

import numpy as np

from streaming_ec_distance_ablation import (
    ALL_VARIANTS,
    VARIANTS,
    Variant,
    convolve_generator,
    convolution_taps,
    expander_generator,
)


def all_messages(message_size: int, prime: int) -> np.ndarray:
    count = prime**message_size
    if count > 10_000_000:
        raise ValueError("exact message space is too large")
    values = np.arange(1, count, dtype=np.int64)
    messages = np.empty((count - 1, message_size), dtype=np.int64)
    for column in range(message_size):
        messages[:, column] = values % prime
        values //= prime
    return messages


def exact_minimum_distance(
    generator: np.ndarray,
    messages: np.ndarray,
    prime: int,
    *,
    chunk_size: int = 1 << 16,
) -> int:
    best = generator.shape[1] + 1
    for start in range(0, len(messages), chunk_size):
        words = (messages[start : start + chunk_size] @ generator) % prime
        best = min(best, int(np.min(np.count_nonzero(words, axis=1))))
    return best


def run_exact_sweep(
    *,
    left_degree: int,
    right_degree: int,
    region_size: int,
    memory: int,
    period: int,
    prime: int,
    seeds: int,
    variants: tuple[Variant, ...] = VARIANTS,
) -> dict[str, object]:
    message_size = right_degree * region_size
    code_size = left_degree * region_size
    messages = all_messages(message_size, prime)
    distances: dict[str, list[int]] = {
        variant.name: [] for variant in variants
    }
    for seed in range(seeds):
        for variant in variants:
            expander = expander_generator(
                variant=variant,
                left_degree=left_degree,
                right_degree=right_degree,
                region_size=region_size,
                prime=prime,
                seed=seed,
            )
            taps = convolution_taps(
                code_size=code_size,
                memory=memory,
                period=period,
                prime=prime,
                seed=seed,
                mode=variant.taps,
            )
            generator = convolve_generator(expander, taps, prime)
            distances[variant.name].append(
                exact_minimum_distance(generator, messages, prime)
            )

    return {
        "scope": "exact distance for the stated small instances only",
        "parameters": {
            "left_degree": left_degree,
            "right_degree": right_degree,
            "region_size": region_size,
            "message_size": message_size,
            "code_size": code_size,
            "memory": memory,
            "period": period,
            "prime": prime,
            "seeds": seeds,
        },
        "variants": {
            name: {
                "minimum": min(values),
                "mean": sum(values) / len(values),
                "histogram": {
                    str(distance): count
                    for distance, count in sorted(Counter(values).items())
                },
            }
            for name, values in distances.items()
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-degree", type=int, default=10)
    parser.add_argument("--right-degree", type=int, default=5)
    parser.add_argument("--region-size", type=int, default=2)
    parser.add_argument("--memory", type=int, default=4)
    parser.add_argument("--period", type=int, default=4)
    parser.add_argument("--prime", type=int, default=3)
    parser.add_argument("--seeds", type=int, default=32)
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=[variant.name for variant in ALL_VARIANTS],
        default=[variant.name for variant in VARIANTS],
    )
    args = parser.parse_args()
    selected = tuple(
        variant for variant in ALL_VARIANTS if variant.name in args.variants
    )
    print(
        json.dumps(
            run_exact_sweep(
                left_degree=args.left_degree,
                right_degree=args.right_degree,
                region_size=args.region_size,
                memory=args.memory,
                period=args.period,
                prime=args.prime,
                seeds=args.seeds,
                variants=selected,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
