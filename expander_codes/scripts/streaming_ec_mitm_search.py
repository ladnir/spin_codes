#!/usr/bin/env python3
"""Meet-in-the-middle low-weight search for reduced field EC codes.

The search splits the message rows into two halves and enumerates every linear
combination in each half.  For each trial, it hashes the two tables on a random
set of code coordinates.  A matching hash forces those coordinates to zero.
The search then measures the complete codeword for every collision.

This is an experimental finder for reduced instances.  It is neither an exact
minimum-distance algorithm nor evidence for the full Goldilocks parameters.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from dataclasses import asdict, dataclass

import numpy as np

from streaming_ec_distance_ablation import (
    ALL_VARIANTS,
    VARIANTS,
    Variant,
    convolve_generator,
    convolution_taps,
    expander_generator,
)
from streaming_ec_exact_ablation import all_messages
from streaming_ec_heuristic_audit import modular_rank


@dataclass(frozen=True)
class MitmResult:
    variant: str
    rank: int
    best_weight: int
    relative_weight: float
    candidate_pairs: int
    nonempty_trials: int
    best_trial: int


def combinations_with_zero(size: int, prime: int) -> np.ndarray:
    nonzero = all_messages(size, prime)
    return np.vstack((np.zeros((1, size), dtype=np.int64), nonzero))


def partial_codewords(
    messages: np.ndarray, generator: np.ndarray, prime: int
) -> np.ndarray:
    values = (messages @ generator) % prime
    if prime <= np.iinfo(np.uint8).max:
        return values.astype(np.uint8)
    if prime <= np.iinfo(np.uint16).max:
        return values.astype(np.uint16)
    return values.astype(np.uint32)


def projection_keys(
    words: np.ndarray,
    coordinates: np.ndarray,
    prime: int,
    *,
    negate: bool,
) -> np.ndarray:
    if prime ** len(coordinates) > np.iinfo(np.uint64).max:
        raise ValueError("projection key does not fit in u64")
    digits = words[:, coordinates].astype(np.uint64)
    if negate:
        digits = (prime - digits) % prime
    powers = np.array(
        [prime**index for index in range(len(coordinates))],
        dtype=np.uint64,
    )
    return np.sum(digits * powers, axis=1, dtype=np.uint64)


def score_cross_product(
    left: np.ndarray,
    right: np.ndarray,
    left_indices: np.ndarray,
    right_indices: np.ndarray,
    prime: int,
    current_best: int,
) -> tuple[int, tuple[int, int] | None, int]:
    best = current_best
    witness: tuple[int, int] | None = None
    tested = 0
    for left_index in left_indices:
        candidates = right[right_indices].astype(np.int64)
        words = (candidates + left[int(left_index)].astype(np.int64)) % prime
        weights = np.count_nonzero(words, axis=1)
        if int(left_index) == 0:
            weights = weights.copy()
            weights[right_indices == 0] = words.shape[1] + 1
        tested += len(right_indices)
        local = int(np.argmin(weights))
        weight = int(weights[local])
        if weight < best:
            best = weight
            witness = (int(left_index), int(right_indices[local]))
    return best, witness, tested


def mitm_search(
    generator: np.ndarray,
    *,
    prime: int,
    zero_coordinates: int,
    trials: int,
    seed: int,
) -> tuple[int, int, int, int]:
    message_size, code_size = generator.shape
    if trials > math.comb(code_size, zero_coordinates):
        raise ValueError("trial count exceeds the number of coordinate subsets")
    left_size = message_size // 2
    right_size = message_size - left_size
    left_messages = combinations_with_zero(left_size, prime)
    right_messages = combinations_with_zero(right_size, prime)
    left_words = partial_codewords(
        left_messages, generator[:left_size], prime
    )
    right_words = partial_codewords(
        right_messages, generator[left_size:], prime
    )

    best = min(
        int(np.min(np.count_nonzero(left_words[1:], axis=1))),
        int(np.min(np.count_nonzero(right_words[1:], axis=1))),
    )
    best_trial = -1
    candidate_pairs = 0
    nonempty_trials = 0
    subset_rng = np.random.Generator(
        np.random.PCG64(np.random.SeedSequence([seed, 0x4D49544D]))
    )
    seen_subsets: set[tuple[int, ...]] = set()
    trial = 0
    while trial < trials:
        coordinates = np.sort(
            subset_rng.choice(code_size, size=zero_coordinates, replace=False)
        )
        signature = tuple(int(value) for value in coordinates)
        if signature in seen_subsets:
            continue
        seen_subsets.add(signature)

        left_keys = projection_keys(
            left_words, coordinates, prime, negate=False
        )
        right_keys = projection_keys(
            right_words, coordinates, prime, negate=True
        )
        left_order = np.argsort(left_keys, kind="stable")
        right_order = np.argsort(right_keys, kind="stable")
        left_sorted = left_keys[left_order]
        right_sorted = right_keys[right_order]

        left_at = 0
        right_at = 0
        trial_candidates = 0
        while left_at < len(left_order) and right_at < len(right_order):
            left_key = left_sorted[left_at]
            right_key = right_sorted[right_at]
            if left_key < right_key:
                left_at = int(np.searchsorted(
                    left_sorted, left_key, side="right"
                ))
                continue
            if right_key < left_key:
                right_at = int(np.searchsorted(
                    right_sorted, right_key, side="right"
                ))
                continue
            left_end = int(np.searchsorted(
                left_sorted, left_key, side="right"
            ))
            right_end = int(np.searchsorted(
                right_sorted, right_key, side="right"
            ))
            best, witness, tested = score_cross_product(
                left_words,
                right_words,
                left_order[left_at:left_end],
                right_order[right_at:right_end],
                prime,
                best,
            )
            trial_candidates += tested
            if witness is not None:
                best_trial = trial
            left_at = left_end
            right_at = right_end

        candidate_pairs += trial_candidates
        nonempty_trials += trial_candidates > 1  # Exclude the zero/zero match.
        trial += 1
    return best, candidate_pairs, nonempty_trials, best_trial


def build_generator(
    *,
    variant: Variant,
    left_degree: int,
    right_degree: int,
    region_size: int,
    memory: int,
    period: int,
    prime: int,
    seed: int,
) -> np.ndarray:
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
    return convolve_generator(expander, taps, prime)


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
    zero_coordinates: int,
    trials: int,
) -> MitmResult:
    generator = build_generator(
        variant=variant,
        left_degree=left_degree,
        right_degree=right_degree,
        region_size=region_size,
        memory=memory,
        period=period,
        prime=prime,
        seed=seed,
    )
    best, candidates, nonempty, best_trial = mitm_search(
        generator,
        prime=prime,
        zero_coordinates=zero_coordinates,
        trials=trials,
        seed=seed,
    )
    return MitmResult(
        variant=variant.name,
        rank=modular_rank(generator, prime),
        best_weight=best,
        relative_weight=best / generator.shape[1],
        candidate_pairs=candidates,
        nonempty_trials=nonempty,
        best_trial=best_trial,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-degree", type=int, default=10)
    parser.add_argument("--right-degree", type=int, default=5)
    parser.add_argument("--region-size", type=int, default=4)
    parser.add_argument("--memory", type=int, default=4)
    parser.add_argument("--period", type=int, default=8)
    parser.add_argument("--prime", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--zero-coordinates", type=int, default=15)
    parser.add_argument("--trials", type=int, default=32)
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=[variant.name for variant in ALL_VARIANTS],
        default=[variant.name for variant in VARIANTS],
    )
    args = parser.parse_args()
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
            zero_coordinates=args.zero_coordinates,
            trials=args.trials,
        )
        for variant in selected
    ]
    print(
        json.dumps(
            {
                "scope": "probabilistic reduced-instance finder",
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
                    "zero_coordinates": args.zero_coordinates,
                    "trials": args.trials,
                },
                "results": [asdict(result) for result in results],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
