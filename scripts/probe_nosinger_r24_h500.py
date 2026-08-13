#!/usr/bin/env python3
"""Conservative floating probe for the no-Singer, graph-r=24 h<=500 ledger.

This is a design diagnostic, not the theorem-facing rational certificate.  It
keeps episode terms through e=32, bounds the remaining binomial tail with a
2^-35 restart cap, and refines every placement bucket in 250-block pieces.
The restart cap combines the exact termination cap 2^-36 with the
geometric-envelope restart
factor ``prefactor/rho < 2``.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from certify_fullsplit_h500_rational import load_outer_coefficients
from punctured_ebch_outer import heterogeneous_coefficients


ROOT = Path(__file__).resolve().parent
LOCAL_SPECTRUM = ROOT / "ebch128_64_spectrum.csv"
N = 2_097_408
DISTANCE = 9 * N // 100
H_MAX = 500
LATE_BLOCKS = 7527
GRAPH_CODIMENSION = 24
EPISODE_MAX = 32
TERMINATION_SHIFT = 35
CHERNOFF_CROSSING = 7526
LOG2_RHO = -0.615274494049
LOG2_PREFACTOR = 0.140751026335
LOG2_POLE = math.log2(2333 / 2373)
OUTER_BLOCKS = 16_386
OUTER_WEIGHT_FLOOR = 44
FIRST_SUBBUCKETS = tuple(
    (LATE_BLOCKS + 250 * index, LATE_BLOCKS + 250 * index + 249)
    for index in range(16)
)
SECOND_SUBBUCKETS = tuple(
    (11_527 + 250 * index, 11_527 + 250 * index + 249)
    for index in range(16)
)
TAIL_SUBBUCKETS = tuple(
    (low, min(low + 249, 32_772)) for low in range(15_527, 32_773, 250)
)


def log2_sum(values: list[float] | np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    if not array.size:
        return -math.inf
    maximum = float(np.max(array))
    if not math.isfinite(maximum):
        return -math.inf
    return maximum + math.log2(float(np.sum(np.exp2(array - maximum))))


def log2_binomial(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return (
        math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    ) / math.log(2)


def survival_log2(live: int) -> float:
    return min(
        0.0,
        -DISTANCE * LOG2_POLE
        + LOG2_PREFACTOR
        + (live - 1) * LOG2_RHO,
    )


def peak_live(
    total: int,
    nonempty: int,
    episodes: int,
    low: int,
    high: int,
    multiplier_log2: float,
) -> int:
    while low < high:
        live = (low + high) // 2
        numerator = (total - live) * (live + 1 - episodes)
        denominator = (total - live + episodes - 1) * (live + 1 - nonempty)
        if math.log2(numerator) - math.log2(denominator) + multiplier_log2 > 0:
            low = live + 1
        else:
            high = live
    return low


def a_term_log2(total: int, nonempty: int, episodes: int, live: int) -> float:
    return log2_binomial(total - live + episodes - 1, episodes - 1) + log2_binomial(
        live - episodes, nonempty - episodes
    )


def suffix_log2(total: int, nonempty: int, episodes: int) -> float:
    if episodes == nonempty + 1:
        return log2_binomial(total, nonempty)
    if not 1 <= episodes <= nonempty:
        return -math.inf

    pieces: list[float] = []
    low_end = min(total, CHERNOFF_CROSSING)
    if nonempty <= low_end:
        peak = peak_live(total, nonempty, episodes, nonempty, low_end, 0.0)
        pieces.append(
            math.log2(low_end - nonempty + 1)
            + a_term_log2(total, nonempty, episodes, peak)
        )

    high_start = max(nonempty, CHERNOFF_CROSSING + 1)
    if high_start <= total:
        peak = peak_live(
            total,
            nonempty,
            episodes,
            high_start,
            total,
            LOG2_RHO,
        )
        pieces.append(
            math.log2(total - high_start + 1)
            + a_term_log2(total, nonempty, episodes, peak)
            + survival_log2(peak)
        )
    return log2_sum(pieces)


def nonempty_log_coefficients(remaining_ones_max: int) -> np.ndarray:
    """Log2 coefficients of ((1+z)^64-1)^x, vectorized over H."""

    block = np.array([log2_binomial(64, weight) for weight in range(65)])
    rows = np.full((remaining_ones_max + 1, remaining_ones_max + 1), -np.inf)
    rows[0, 0] = 0.0
    for nonempty in range(1, remaining_ones_max + 1):
        row = np.full(remaining_ones_max + 1, -np.inf)
        previous = rows[nonempty - 1]
        for weight in range(1, 65):
            row[weight:] = np.logaddexp2(
                row[weight:], previous[:-weight] + block[weight]
            )
        rows[nonempty] = row
    return rows


def inner_log_bounds(total: int, coefficients: np.ndarray) -> np.ndarray:
    remaining_ones_max = coefficients.shape[1] - 1
    suffix = np.full((EPISODE_MAX + 1, remaining_ones_max + 1), -np.inf)
    choose_episodes = np.full_like(suffix, -np.inf)
    for episodes in range(1, EPISODE_MAX + 1):
        for nonempty in range(episodes - 1, remaining_ones_max + 1):
            suffix[episodes, nonempty] = suffix_log2(
                total, nonempty, episodes
            )
            choose_episodes[episodes, nonempty] = log2_binomial(
                nonempty + 1, episodes
            )

    result = np.empty(remaining_ones_max + 1)
    denominators = np.array(
        [
            log2_binomial(64 * total, remaining_ones)
            for remaining_ones in range(remaining_ones_max + 1)
        ]
    )
    for remaining_ones in range(remaining_ones_max + 1):
        terms = [survival_log2(total)]
        denominator = denominators[remaining_ones]
        for episodes in range(1, min(EPISODE_MAX, remaining_ones + 1) + 1):
            indices = np.arange(episodes - 1, remaining_ones + 1)
            terms.append(
                log2_sum(
                    coefficients[indices, remaining_ones]
                    + choose_episodes[episodes, indices]
                    + suffix[episodes, indices]
                )
                - denominator
                - TERMINATION_SHIFT * episodes
            )
        if remaining_ones >= EPISODE_MAX:
            terms.append(
                log2_sum(
                    [
                        log2_binomial(remaining_ones + 1, episodes)
                        - TERMINATION_SHIFT * episodes
                        for episodes in range(
                            EPISODE_MAX + 1, remaining_ones + 2
                        )
                    ]
                )
            )
        result[remaining_ones] = min(0.0, log2_sum(terms))
    return result


def placement_bucket_log2(
    low: int,
    high: int,
    inner: np.ndarray,
    outer: list[int],
    h_min: int,
    h_max: int,
) -> float:
    rows: list[float] = []
    for weight in range(max(OUTER_WEIGHT_FLOOR, h_min), h_max + 1):
        if not outer[weight]:
            continue
        first_weights = [
            log2_binomial(64, first_weight)
            + math.log2(high - low + 1)
            + log2_binomial(64 * (high - 1), weight - first_weight)
            + inner[weight - first_weight]
            for first_weight in range(1, min(64, weight) + 1)
        ]
        rows.append(
            math.log2(outer[weight])
            + log2_sum(first_weights)
            - log2_binomial(N, weight)
        )
    return log2_sum(rows)


def main() -> None:
    global N, DISTANCE, LATE_BLOCKS, CHERNOFF_CROSSING, OUTER_BLOCKS
    global OUTER_WEIGHT_FLOOR, FIRST_SUBBUCKETS, SECOND_SUBBUCKETS, TAIL_SUBBUCKETS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h-min", type=int, default=1)
    parser.add_argument("--h-max", type=int, default=H_MAX)
    parser.add_argument(
        "--power2-punctured",
        action="store_true",
        help="use 256 punctured EBCH blocks and N=2^21",
    )
    args = parser.parse_args()
    if not 1 <= args.h_min <= args.h_max:
        raise SystemExit("no-Singer r24 episode probe: invalid weight interval")

    if args.power2_punctured:
        N = 1 << 21
        DISTANCE = 9 * N // 100
        LATE_BLOCKS = 7526
        CHERNOFF_CROSSING = 7525
        OUTER_BLOCKS = 16_386
        OUTER_WEIGHT_FLOOR = 43
        first_start = LATE_BLOCKS
        second_start = first_start + 4_000
        tail_start = second_start + 4_000
        FIRST_SUBBUCKETS = tuple(
            (first_start + 250 * index, first_start + 250 * index + 249)
            for index in range(16)
        )
        SECOND_SUBBUCKETS = tuple(
            (second_start + 250 * index, second_start + 250 * index + 249)
            for index in range(16)
        )
        TAIL_SUBBUCKETS = tuple(
            (low, min(low + 249, 32_768))
            for low in range(tail_start, 32_769, 250)
        )

    coefficients = nonempty_log_coefficients(args.h_max - 1)
    if args.power2_punctured:
        outer = heterogeneous_coefficients(
            full_blocks=16_130,
            punctured_blocks=256,
            h_max=args.h_max,
            spectrum_path=LOCAL_SPECTRUM,
        )
    else:
        outer = load_outer_coefficients(
            LOCAL_SPECTRUM, blocks=OUTER_BLOCKS, h_max=args.h_max
        )
    for weight in range(OUTER_WEIGHT_FLOOR):
        outer[weight] = 0

    ultra_late = log2_sum(
        [
            math.log2(outer[weight])
            + log2_binomial(64 * LATE_BLOCKS, weight)
            - log2_binomial(N, weight)
            for weight in range(max(OUTER_WEIGHT_FLOOR, args.h_min), args.h_max + 1)
            if outer[weight]
        ]
    )

    first_parts = []
    for low, high in FIRST_SUBBUCKETS:
        first_parts.append(
            placement_bucket_log2(
                low, high, inner_log_bounds(low, coefficients), outer,
                args.h_min, args.h_max
            )
        )
    first_bucket = log2_sum(first_parts)

    second_parts = []
    for low, high in SECOND_SUBBUCKETS:
        second_parts.append(
            placement_bucket_log2(
                low, high, inner_log_bounds(low, coefficients), outer,
                args.h_min, args.h_max
            )
        )
    second_bucket = log2_sum(second_parts)

    tail_parts = []
    for low, high in TAIL_SUBBUCKETS:
        tail_parts.append(
            placement_bucket_log2(
                low, high, inner_log_bounds(low, coefficients), outer,
                args.h_min, args.h_max
            )
        )

    tail_bucket = log2_sum(tail_parts)
    ambient_terms = [ultra_late, first_bucket, second_bucket, tail_bucket]
    ambient = log2_sum(ambient_terms)
    graph = ambient - GRAPH_CODIMENSION
    print(f"no-Singer r=24 h={args.h_min}..{args.h_max} floating ledger probe")
    print(
        f"construction={'power2-punctured' if args.power2_punctured else 'unpunctured'} "
        f"N={N} distance={DISTANCE} late_blocks={LATE_BLOCKS}"
    )
    print(f"ambient_terms_log2={ambient_terms}")
    print(f"ambient_total_log2_upper_approx={ambient:.12f}")
    print(f"graph_total_log2_upper_approx={graph:.12f}")
    print(f"diagnostic_margin_beyond_40_bits={-40.0 - graph:.12f}")
    print("status=DIAGNOSTIC_ONLY_NOT_RATIONAL_CERTIFICATE")


if __name__ == "__main__":
    main()
