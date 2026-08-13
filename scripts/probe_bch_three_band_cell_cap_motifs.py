#!/usr/bin/env python3
"""Solver-free cell-cap bounds for fixed three-band collision motifs.

For each unknown fixed-band split count N[a,b,c], every exact aggregate row
containing that cell is an individual upper bound.  On each exact total-weight
diagonal, greedily filling the largest objective cells up to the smallest of
those caps is therefore a rigorous relaxation of the missing joint split
enumerator.

Collision clusters are separated into one-codeword factors as follows.
Singletons contribute q^a.  More generally the coordinate OR kernel has the
nonnegative rank-one decomposition

    q^OR(x_1,...,x_d) = q + (1-q) product_i 1[x_i=0].

Tensoring over coordinates and averaging fixed-weight slices preserves a
completely positive symmetric tensor.  Generalized Hölder therefore gives

    K_d(a_1,...,a_d) <= product_i K_d(a_i,...,a_i)^(1/d).

The resulting factor for each active block depends only on its own three band
weights.  This script enumerates compatible colored partition patterns and
applies the greedy fixed-band cell-cap bound to each distinct vertex type.
Floating arithmetic makes this a design diagnostic; all combinatorial counts
and relaxations are exact, but the final certificate needs outward rounding.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from probe_bch_three_band_transport_lp import (
    BAND_SIZES,
    CODE_SIZE,
    exact_rows,
    hypergeometric_pair_diagonal,
    remove_trivial_codewords,
    triples,
)
from probe_three_band_spreading import (
    incidence_components,
    pairwise_compatible,
)
from probe_two_band_spreading import set_partitions


BLOCKS = 16384
LANES = 64
GRAPH_CODIMENSION = 24


def cell_caps(
    states: list[tuple[int, int, int]],
    rows: list[tuple[str, int, list[int]]],
) -> tuple[list[int], dict[int, int], dict[int, list[int]]]:
    caps = [CODE_SIZE] * len(states)
    diagonal_mass: dict[int, int] = {}
    by_total: dict[int, list[int]] = defaultdict(list)
    for index, state in enumerate(states):
        by_total[sum(state)].append(index)
    for name, rhs, indices in rows:
        if name.startswith("total_"):
            diagonal_mass[int(name.removeprefix("total_"))] = rhs
        for index in indices:
            caps[index] = min(caps[index], rhs)
    return caps, diagonal_mass, by_total


def factor_table(columns: int, cluster_size: int, pole: float) -> np.ndarray:
    if cluster_size == 1:
        return np.array([pole**weight for weight in range(columns + 1)])
    return np.array(
        [
            cluster_diagonal(columns, weight, cluster_size, pole)
            ** (1.0 / cluster_size)
            for weight in range(columns + 1)
        ]
    )


def cluster_diagonal(
    columns: int, weight: int, cluster_size: int, pole: float
) -> float:
    if cluster_size == 2:
        return hypergeometric_pair_diagonal(columns, weight, pole)
    distribution = np.zeros(columns + 1)
    distribution[0] = 1.0
    denominator = math.comb(columns, weight)
    for _ in range(cluster_size):
        following = np.zeros(columns + 1)
        for union, mass in enumerate(distribution):
            if not mass:
                continue
            for intersection in range(
                max(0, union + weight - columns), min(union, weight) + 1
            ):
                new_union = union + weight - intersection
                following[new_union] += (
                    mass
                    * math.comb(union, intersection)
                    * math.comb(columns - union, weight - intersection)
                    / denominator
                )
        distribution = following
    return sum(mass * pole**union for union, mass in enumerate(distribution))


def greedy_factor_bound(
    states: list[tuple[int, int, int]],
    caps: list[int],
    diagonal_mass: dict[int, int],
    by_total: dict[int, list[int]],
    factors: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> float:
    total = 0.0
    for weight, mass in diagonal_mass.items():
        if not mass:
            continue
        scored = sorted(
            (
                factors[0][states[index][0]]
                * factors[1][states[index][1]]
                * factors[2][states[index][2]],
                caps[index],
            )
            for index in by_total[weight]
            if caps[index]
        )
        remaining = mass
        diagonal_value = 0.0
        for value, cap in reversed(scored):
            take = min(remaining, cap)
            diagonal_value += take * value
            remaining -= take
            if not remaining:
                break
        if remaining:
            raise SystemExit(
                f"cell-cap motif: caps fail to cover weight {weight} by {remaining}"
            )
        total += diagonal_value
    return total


def cluster_sizes(partition: tuple[int, ...]) -> tuple[int, ...]:
    counts = [0] * (max(partition) + 1)
    for label in partition:
        counts[label] += 1
    return tuple(counts[label] for label in partition)


def logaddexp2(left: float, right: float) -> float:
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    maximum = max(left, right)
    return maximum + math.log2(2 ** (left - maximum) + 2 ** (right - maximum))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group-pole", type=float, default=0.181)
    parser.add_argument("--blocks", type=int, default=4)
    parser.add_argument(
        "--max-cluster",
        type=int,
        default=None,
        help="omit patterns containing a larger same-tile cluster",
    )
    parser.add_argument(
        "--min-collision-bands",
        type=int,
        default=0,
        help="omit patterns colliding in fewer bands",
    )
    parser.add_argument(
        "--show-factors",
        action="store_true",
        help="print every cached one-codeword factor",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="also write the compact result to this generated artifact",
    )
    args = parser.parse_args()
    if not 0 < args.group_pole < 1 or not 2 <= args.blocks <= 7:
        raise SystemExit("cell-cap motif: invalid arguments")

    states = triples()
    rows = remove_trivial_codewords(states, exact_rows(states))
    caps, diagonal_mass, by_total = cell_caps(states, rows)
    factor_tables = {
        (band, cluster_size): factor_table(
            BAND_SIZES[band], cluster_size, args.group_pole
        )
        for band in range(3)
        for cluster_size in range(1, args.blocks + 1)
    }
    bound_cache: dict[tuple[int, int, int], float] = {}

    def vertex_bound(vertex_type: tuple[int, int, int]) -> float:
        if vertex_type not in bound_cache:
            factors = tuple(
                factor_tables[band, vertex_type[band]] for band in range(3)
            )
            bound_cache[vertex_type] = greedy_factor_bound(
                states, caps, diagonal_mass, by_total, factors
            )
        return bound_cache[vertex_type]

    partitions = list(set_partitions(args.blocks))
    compatibility_masks = []
    for left in range(len(partitions)):
        mask = 0
        for right in range(len(partitions)):
            if pairwise_compatible(partitions[left], partitions[right]):
                mask |= 1 << right
        compatibility_masks.append(mask)
    sizes = [cluster_sizes(partition) for partition in partitions]
    total_log2 = -math.inf
    connected_log2 = -math.inf
    patterns = 0
    connected_patterns = 0
    worst: list[
        tuple[float, tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]
    ] = []
    for first, first_partition in enumerate(partitions):
        second_mask = compatibility_masks[first]
        while second_mask:
            second_low = second_mask & -second_mask
            second = second_low.bit_length() - 1
            second_mask ^= second_low
            second_partition = partitions[second]
            third_mask = compatibility_masks[first] & compatibility_masks[second]
            while third_mask:
                third_low = third_mask & -third_mask
                third = third_low.bit_length() - 1
                third_mask ^= third_low
                third_partition = partitions[third]
                partition_triple = (
                    first_partition,
                    second_partition,
                    third_partition,
                )
                components = incidence_components(partition_triple)
                vertex_types = tuple(
                    (sizes[first][vertex], sizes[second][vertex], sizes[third][vertex])
                    for vertex in range(args.blocks)
                )
                if args.max_cluster is not None and any(
                    cluster_size > args.max_cluster
                    for vertex_type in vertex_types
                    for cluster_size in vertex_type
                ):
                    continue
                collision_bands = sum(
                    any(cluster_size > 1 for cluster_size in sizes[index])
                    for index in (first, second, third)
                )
                if collision_bands < args.min_collision_bands:
                    continue
                value_log2 = (
                    components * math.log2(BLOCKS)
                    + (args.blocks - components) * math.log2(LANES)
                    - math.lgamma(args.blocks + 1) / math.log(2)
                    - GRAPH_CODIMENSION
                    + sum(math.log2(vertex_bound(kind)) for kind in vertex_types)
                )
                total_log2 = logaddexp2(total_log2, value_log2)
                patterns += 1
                if components == 1:
                    connected_log2 = logaddexp2(connected_log2, value_log2)
                    connected_patterns += 1
                    item = (value_log2, partition_triple, vertex_types)
                    if len(worst) < 10:
                        heapq.heappush(worst, item)
                    elif item[0] > worst[0][0]:
                        heapq.heapreplace(worst, item)

    worst.sort(reverse=True)
    print("fixed three-band cell-cap motif ledger")
    print(
        f"blocks={args.blocks} group_pole={args.group_pole:.12f} "
        f"patterns={patterns} connected_patterns={connected_patterns}"
    )
    print(f"all_patterns_log2={total_log2:.12f}")
    print(f"connected_patterns_log2={connected_log2:.12f}")
    for rank, (value, partition_triple, vertex_types) in enumerate(worst[:10], 1):
        print(
            f"worst_rank={rank} log2={value:.12f} "
            f"partitions={partition_triple} vertex_types={vertex_types}"
        )
    if args.show_factors:
        for kind, bound in sorted(bound_cache.items()):
            print(f"vertex_type={kind} raw_factor_log2={math.log2(bound):.12f}")
    if args.summary_json is not None:
        args.summary_json.write_text(
            json.dumps(
                {
                    "blocks": args.blocks,
                    "group_pole": args.group_pole,
                    "patterns": patterns,
                    "connected_patterns": connected_patterns,
                    "all_patterns_log2": total_log2,
                    "connected_patterns_log2": connected_log2,
                    "worst": [
                        {
                            "log2": value,
                            "partitions": partition_triple,
                            "vertex_types": vertex_types,
                        }
                        for value, partition_triple, vertex_types in worst
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    print("status=RIGOROUS_RELAXATION_FLOATING_ARITHMETIC")


if __name__ == "__main__":
    main()
