#!/usr/bin/env python3
"""Three-band pairwise-linear outer-spreading diagnostics.

The 128 coordinates of every EBCH word are split into bands of sizes
42, 43, and 43.  Each band has 256 tiles of 64 blocks, and the proposed
block-to-tile map has pairwise cell capacity one.  This script checks:

* the exact joint three-band weight law and its density over the product of
  its marginals;
* exact-joint star moments (one colliding band, two dispersed bands); and
* all pairwise-linear partition triples through a requested small block
  count, using a rigorous layout-embedding upper bound M^c 64^(s-c) for an
  incidence hypergraph with c connected components.

Long-double arithmetic is used.  The result is a design diagnostic and must
be replaced by outward-rounded arithmetic in the final certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from probe_two_band_spreading import cluster_moments, set_partitions
from punctured_ebch_outer import full_spectrum


BLOCKS = 16384
TILES = 256
LANES = 64
BAND_SIZES = (42, 43, 43)
NONTRIVIAL_CARDINALITY = (1 << 64) - 2
GRAPH_CODIMENSION = 24


def joint_weight_law() -> np.ndarray:
    spectrum = full_spectrum()
    joint = np.zeros(tuple(size + 1 for size in BAND_SIZES), dtype=np.longdouble)
    for first in range(BAND_SIZES[0] + 1):
        for second in range(BAND_SIZES[1] + 1):
            for third in range(BAND_SIZES[2] + 1):
                total = first + second + third
                if not 0 < total < 128 or not spectrum[total]:
                    continue
                joint[first, second, third] = (
                    np.longdouble(spectrum[total])
                    / NONTRIVIAL_CARDINALITY
                    * np.longdouble(math.comb(BAND_SIZES[0], first))
                    * np.longdouble(math.comb(BAND_SIZES[1], second))
                    * np.longdouble(math.comb(BAND_SIZES[2], third))
                    / np.longdouble(math.comb(128, total))
                )
    if abs(float(np.sum(joint)) - 1.0) > 1e-15:
        raise SystemExit("three-band spreading: joint law does not sum to one")
    return joint


def product_density(
    joint: np.ndarray, marginals: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> tuple[float, tuple[int, int, int]]:
    maximum = 0.0
    argument = (0, 0, 0)
    for index in zip(*np.nonzero(joint)):
        density = float(
            joint[index]
            / (
                marginals[0][index[0]]
                * marginals[1][index[1]]
                * marginals[2][index[2]]
            )
        )
        if density > maximum:
            maximum = density
            argument = index
    return maximum, argument


def union_transitions(columns: int) -> np.ndarray:
    transition = np.zeros((columns + 1, columns + 1, columns + 1))
    for union in range(columns + 1):
        for weight in range(columns + 1):
            denominator = math.comb(columns, weight)
            for intersection in range(
                max(0, union + weight - columns), min(union, weight) + 1
            ):
                transition[union, weight, union + weight - intersection] = (
                    math.comb(union, intersection)
                    * math.comb(columns - union, weight - intersection)
                    / denominator
                )
    return transition


def exact_star_rows(
    joint: np.ndarray, pole: float, maximum_blocks: int, clustered_band: int
) -> list[tuple[int, float]]:
    columns = BAND_SIZES[clustered_band]
    axes = tuple(axis for axis in range(3) if axis != clustered_band)
    outside_weighted = np.zeros(columns + 1)
    for index in zip(*np.nonzero(joint)):
        outside = index[axes[0]] + index[axes[1]]
        outside_weighted[index[clustered_band]] += float(joint[index]) * pole**outside
    transition = union_transitions(columns)
    kernel = np.einsum("uav,a->uv", transition, outside_weighted, optimize=True)
    powers = np.array([pole**weight for weight in range(columns + 1)])
    distribution = np.zeros(columns + 1)
    distribution[0] = 1.0
    scale_log2 = 0.0
    rows = []
    for blocks in range(1, maximum_blocks + 1):
        distribution = distribution @ kernel
        scale = float(np.max(distribution))
        distribution /= scale
        scale_log2 += math.log2(scale)
        moment_log2 = scale_log2 + math.log2(float(distribution @ powers))
        # One common tile in the clustered band.  Pairwise capacity one makes
        # every block use a distinct tile in each of the other two bands.
        assignment_log2 = (1 - blocks) * math.log2(TILES)
        for _ in range(2):
            assignment_log2 += sum(
                math.log2(TILES - offset) for offset in range(blocks)
            ) - blocks * math.log2(TILES)
        outer = (
            blocks
            * (math.log2(BLOCKS) + math.log2(NONTRIVIAL_CARDINALITY))
            - math.lgamma(blocks + 1) / math.log(2)
            - GRAPH_CODIMENSION
            + assignment_log2
            + moment_log2
        )
        rows.append((blocks, outer))
    return rows


def pairwise_compatible(first: tuple[int, ...], second: tuple[int, ...]) -> bool:
    return len(set(zip(first, second))) == len(first)


def incidence_components(partitions: tuple[tuple[int, ...], ...]) -> int:
    blocks = len(partitions[0])
    parent = list(range(blocks))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def merge(left: int, right: int) -> None:
        left = find(left)
        right = find(right)
        if left != right:
            parent[right] = left

    for partition in partitions:
        first: dict[int, int] = {}
        for block, label in enumerate(partition):
            if label in first:
                merge(block, first[label])
            else:
                first[label] = block
    return len({find(block) for block in range(blocks)})


def partition_moment(partition: tuple[int, ...], cluster: np.ndarray) -> float:
    degrees = [0] * (max(partition) + 1)
    for label in partition:
        degrees[label] += 1
    return math.prod(float(cluster[degree]) for degree in degrees)


def enumerate_partition_triples(
    blocks: int, clusters: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> tuple[float, int]:
    partitions = list(set_partitions(blocks))
    count = len(partitions)
    compatible = [
        [pairwise_compatible(partitions[left], partitions[right]) for right in range(count)]
        for left in range(count)
    ]
    weights = [
        [partition_moment(partition, clusters[band]) for partition in partitions]
        for band in range(3)
    ]
    total = 0.0
    patterns = 0
    for first, first_partition in enumerate(partitions):
        for second, second_partition in enumerate(partitions):
            if not compatible[first][second]:
                continue
            for third, third_partition in enumerate(partitions):
                if not compatible[first][third] or not compatible[second][third]:
                    continue
                components = incidence_components(
                    (first_partition, second_partition, third_partition)
                )
                # Embed the first edge of every component in at most M ways.
                # Every later connected edge shares a tile and has at most 64
                # choices.  Relative to M^s this is 256^{-(s-c)}.
                total += (
                    weights[0][first]
                    * weights[1][second]
                    * weights[2][third]
                    * TILES ** (-(blocks - components))
                )
                patterns += 1
    return total, patterns


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group-pole", type=float, default=0.181)
    parser.add_argument("--enumerate-max", type=int, default=6)
    parser.add_argument("--star-max", type=int, default=64)
    args = parser.parse_args()
    if not 0 < args.group_pole < 1 or not 1 <= args.enumerate_max <= 7:
        raise SystemExit("three-band spreading: invalid arguments")

    joint = joint_weight_law()
    marginals = (
        np.sum(joint, axis=(1, 2)),
        np.sum(joint, axis=(0, 2)),
        np.sum(joint, axis=(0, 1)),
    )
    density, density_argument = product_density(joint, marginals)
    clusters = tuple(
        cluster_moments(args.enumerate_max, args.group_pole, marginal)
        for marginal in marginals
    )

    print("three-band pairwise-linear spreading")
    print(f"band_sizes={BAND_SIZES}")
    print(f"group_pole={args.group_pole:.12f}")
    print(
        f"product_density_log2={math.log2(density):.12f} "
        f"at={density_argument}"
    )
    for blocks in range(1, args.enumerate_max + 1):
        moment, patterns = enumerate_partition_triples(blocks, clusters)
        outer = (
            blocks
            * (
                math.log2(BLOCKS)
                + math.log2(NONTRIVIAL_CARDINALITY)
                + math.log2(density)
            )
            - math.lgamma(blocks + 1) / math.log(2)
            - GRAPH_CODIMENSION
            + math.log2(moment)
        )
        print(
            f"active_blocks={blocks:2d} linear_patterns={patterns} "
            f"outer_weighted_log2={outer:.12f}"
        )

    for clustered_band in range(3):
        rows = exact_star_rows(
            joint, args.group_pole, args.star_max, clustered_band
        )
        worst = max(rows, key=lambda row: row[1])
        nontrivial_worst = max(rows[1:], key=lambda row: row[1]) if len(rows) > 1 else worst
        print(
            f"exact_star_clustered_band={clustered_band} "
            f"worst_log2={worst[1]:.12f}@s={worst[0]} "
            f"worst_s_ge_2_log2={nontrivial_worst[1]:.12f}@s={nontrivial_worst[0]}"
        )
    print("status=DIAGNOSTIC_NEEDS_LARGE_COMPONENT_BOUND_AND_OUTWARD_ARITHMETIC")


if __name__ == "__main__":
    main()
