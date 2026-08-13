#!/usr/bin/env python3
"""Exhaust the 9! ordered packet-profile simplex chambers.

For a permutation ``pi``, the chamber

    a_pi[0] >= a_pi[1] >= ... >= a_pi[8] >= 0

is a simplex.  Its nine vertices are the uniform profiles on the nested prefix
sets of ``pi``.  A fixed witness covers the entire chamber if it is below the
per-profile target at all those vertices.  When class zero is first, the pure
zero vertex is excluded by total outer weight at least 21; clipping that
simplex adds the eight intersections from the pure-zero vertex to the other
prefix vertices.

Only 511 uniform-subset vertices and 255 clipped vertices exist globally.
This script precomputes the safe-witness bitmask for each one, then checks all
362880 chambers by integer mask intersection.  A complete PASS is exhaustive
and sample-independent, apart from the cached witness and binary64 arithmetic.
"""

from __future__ import annotations

import argparse
import itertools
import math
from collections import Counter
from pathlib import Path

import numpy as np

from probe_packet8_adaptive_simplex_atlas import load_witness_cache
from probe_packet8_profile_simplex_landscape import (
    ANCHORS,
    D,
    K,
    M,
    N,
    PROFILE_COUNT_LOG2,
    FixedWitness,
    load_anchor_file,
    normalization_log2,
)
from probe_packet8_shared_witness_io import load_shared_witness_arrays


def add_full_bijection(witnesses: list[FixedWitness]) -> list[FixedWitness]:
    constant = N * math.log2(11) - (N - D) * math.log2(10) + K
    return witnesses + [FixedWitness("full_bijection", np.zeros(9), constant)]


def uniform_subset(mask: int) -> np.ndarray:
    point = np.zeros(9, dtype=np.float64)
    active = [index for index in range(9) if mask >> index & 1]
    point[active] = M / len(active)
    return point


def clipped_from_zero(point: np.ndarray, minimum_weight: int) -> np.ndarray:
    zero = np.zeros(9, dtype=np.float64)
    zero[0] = M
    weight = float(np.dot(np.arange(9), point))
    fraction = minimum_weight / weight
    return zero + fraction * (point - zero)


def point_values(
    point: np.ndarray, constants: np.ndarray, charges: np.ndarray
) -> np.ndarray:
    return (
        constants
        - point @ charges.T
        - float(normalization_log2(point[None, :])[0])
    )


def safe_mask(values: np.ndarray, target: float) -> int:
    result = 0
    for index in np.flatnonzero(values <= target):
        result |= 1 << int(index)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--minimum-outer-weight", type=int, default=21)
    parser.add_argument("--show-uncovered", type=int, default=20)
    parser.add_argument("--upgraded-cache", type=Path)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    args = parser.parse_args()
    if args.show_uncovered <= 0:
        raise SystemExit("ordered chambers: invalid uncovered-row limit")

    if args.upgraded_cache:
        names, constants, charges = load_shared_witness_arrays(
            args.upgraded_cache, args.extra_shared_report
        )
        full = add_full_bijection([])[0]
        names.append(full.name)
        constants = np.append(constants, full.constant_log2)
        charges = np.vstack((charges, full.linear_charge))
    else:
        if args.extra_shared_report:
            raise SystemExit(
                "ordered chambers: --extra-shared-report requires --upgraded-cache"
            )
        extras = load_anchor_file(args.atlas)
        anchors = ANCHORS + extras
        witnesses = load_witness_cache(
            args.atlas.with_suffix(".witnesses.npz"), anchors
        )
        if witnesses is None:
            raise SystemExit("ordered chambers: matching witness cache missing")
        witnesses = add_full_bijection(witnesses)
        names = [row.name for row in witnesses]
        constants = np.asarray([row.constant_log2 for row in witnesses])
        charges = np.vstack([row.linear_charge for row in witnesses])
    target = -40.0 - PROFILE_COUNT_LOG2
    all_witnesses = (1 << len(constants)) - 1

    uniform_values: dict[int, np.ndarray] = {}
    uniform_masks: dict[int, int] = {}
    clipped_values: dict[int, np.ndarray] = {}
    clipped_masks: dict[int, int] = {}
    for subset in range(1, 1 << 9):
        point = uniform_subset(subset)
        values = point_values(point, constants, charges)
        uniform_values[subset] = values
        uniform_masks[subset] = safe_mask(values, target)
        if subset & 1 and subset != 1:
            clipped = clipped_from_zero(point, args.minimum_outer_weight)
            values = point_values(clipped, constants, charges)
            clipped_values[subset] = values
            clipped_masks[subset] = safe_mask(values, target)

    uncovered_uniform = [
        (
            float(np.min(uniform_values[subset])),
            subset,
            int(np.argmin(uniform_values[subset])),
        )
        for subset in uniform_values
        if not uniform_masks[subset]
    ]
    uncovered_uniform.sort(reverse=True)
    uncovered_clipped = [
        (
            float(np.min(clipped_values[subset])),
            subset,
            int(np.argmin(clipped_values[subset])),
        )
        for subset in clipped_values
        if not clipped_masks[subset]
    ]
    uncovered_clipped.sort(reverse=True)

    total_uniform_edges = 0
    covered_uniform_edges = 0
    uncovered_edge_rows = []
    for subset in range(1, 1 << 9):
        for packet_class in range(9):
            if subset >> packet_class & 1:
                continue
            larger = subset | (1 << packet_class)
            total_uniform_edges += 1
            shared = uniform_masks[subset] & uniform_masks[larger]
            if shared:
                covered_uniform_edges += 1
            elif len(uncovered_edge_rows) < args.show_uncovered:
                stacked = np.vstack(
                    (uniform_values[subset], uniform_values[larger])
                )
                maxima = np.max(stacked, axis=0)
                leader_index = int(np.argmin(maxima))
                uncovered_edge_rows.append(
                    (
                        float(maxima[leader_index]),
                        subset,
                        larger,
                        leader_index,
                    )
                )

    total_comparable_edges = 0
    covered_comparable_edges = 0
    for subset in range(2, 1 << 9):
        complement = ((1 << 9) - 1) ^ subset
        addition = complement
        while addition:
            larger = subset | addition
            total_comparable_edges += 1
            if uniform_masks[subset] & uniform_masks[larger]:
                covered_comparable_edges += 1
            addition = (addition - 1) & complement

    covered = 0
    uncovered = 0
    leaders: Counter[str] = Counter()
    uncovered_rows = []
    for permutation in itertools.permutations(range(9)):
        mask = all_witnesses
        subset = 0
        value_rows = []
        labels = []
        starts_zero = permutation[0] == 0
        for depth, packet_class in enumerate(permutation):
            subset |= 1 << packet_class
            if not (starts_zero and depth == 0):
                mask &= uniform_masks[subset]
                if len(uncovered_rows) < args.show_uncovered:
                    value_rows.append(uniform_values[subset])
                    labels.append(f"uniform_{subset:03x}")
            if starts_zero and depth:
                mask &= clipped_masks[subset]
                if len(uncovered_rows) < args.show_uncovered:
                    value_rows.append(clipped_values[subset])
                    labels.append(f"clipped_{subset:03x}")
            if not mask and len(uncovered_rows) >= args.show_uncovered:
                break
        if mask:
            covered += 1
            leader_index = (mask & -mask).bit_length() - 1
            leaders[names[leader_index]] += 1
            continue

        uncovered += 1
        if len(uncovered_rows) < args.show_uncovered:
            matrix = np.vstack(value_rows)
            maxima = np.max(matrix, axis=0)
            leader_index = int(np.argmin(maxima))
            worst_index = int(np.argmax(matrix[:, leader_index]))
            uncovered_rows.append(
                (
                    float(maxima[leader_index]),
                    permutation,
                    names[leader_index],
                    labels[worst_index],
                )
            )

    print("packet-8 ordered simplex-chamber cover probe")
    print(
        f"atlas={args.atlas} witnesses={len(constants)} target={target:.12f} "
        f"minimum_outer_weight={args.minimum_outer_weight}"
    )
    print(f"uniform_subset_vertices={len(uniform_values)}")
    print(f"clipped_zero_vertices={len(clipped_values)}")
    print(
        f"covered_uniform_vertices={len(uniform_values)-len(uncovered_uniform)} "
        f"uncovered_uniform_vertices={len(uncovered_uniform)}"
    )
    print(
        f"covered_clipped_vertices={len(clipped_values)-len(uncovered_clipped)} "
        f"uncovered_clipped_vertices={len(uncovered_clipped)}"
    )
    print(
        f"covered_uniform_hasse_edges={covered_uniform_edges} "
        f"uncovered_uniform_hasse_edges={total_uniform_edges-covered_uniform_edges} "
        f"total_uniform_hasse_edges={total_uniform_edges}"
    )
    print(
        f"covered_uniform_comparable_edges={covered_comparable_edges} "
        f"uncovered_uniform_comparable_edges="
        f"{total_comparable_edges-covered_comparable_edges} "
        f"total_uniform_comparable_edges={total_comparable_edges}"
    )
    for rank, (score, subset, leader_index) in enumerate(
        uncovered_uniform[: args.show_uncovered], 1
    ):
        print(
            f"uncovered_uniform_rank={rank} subset={subset:03x} "
            f"score={score:.9f} leader={names[leader_index]}"
        )
    for rank, (score, subset, larger, leader_index) in enumerate(
        uncovered_edge_rows, 1
    ):
        print(
            f"uncovered_edge_rank={rank} edge={subset:03x}->{larger:03x} "
            f"score={score:.9f} leader={names[leader_index]}"
        )
    for rank, (score, subset, leader_index) in enumerate(
        uncovered_clipped[: args.show_uncovered], 1
    ):
        print(
            f"uncovered_clipped_rank={rank} subset={subset:03x} "
            f"score={score:.9f} leader={names[leader_index]}"
        )
    print(f"total_chambers={math.factorial(9)}")
    print(f"covered_chambers={covered}")
    print(f"uncovered_chambers={uncovered}")
    for rank, (score, permutation, leader, worst_label) in enumerate(
        uncovered_rows, 1
    ):
        print(
            f"rank={rank} score={score:.9f} leader={leader} "
            f"worst_vertex={worst_label} order="
            + ",".join(map(str, permutation))
        )
    for leader, count in leaders.most_common(20):
        print(f"leader={leader} covered_chambers={count}")
    print(f"complete_ordered_chamber_cover={'PASS' if not uncovered else 'NO'}")
    print("status=DIAGNOSTIC_BINARY64_EXHAUSTIVE_ORDERED_CHAMBERS")


if __name__ == "__main__":
    main()
