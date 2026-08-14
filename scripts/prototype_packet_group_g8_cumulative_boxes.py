#!/usr/bin/env python3
"""Prototype exact cumulative-mass boxes for one g=8 support.

This script studies geometry and counting only.  It does not evaluate a
witness and does not emit certificate shards.  The output binds the frozen
manifest and records exact cells for a balanced dyadic cumulative grid.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import certify_packet_group_g8_support_shard as verifier


SCHEMA = "permute-conv.packet-group-g8-cumulative-box-prototype.v1"
MANIFEST_SHA256 = "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
MASS = 262144
CUTOFF = 21
CLASSES = 9
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "G8_SUPPORT_MANIFEST.json"

Profile = tuple[Fraction, ...]


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def support_of(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(CLASSES) if mask & (1 << index))


def balanced_bins(total: int, bins: int) -> tuple[tuple[int, int], ...]:
    """Partition the integers 0,...,total into contiguous balanced bins."""

    if bins < 1 or bins > total + 1:
        raise ValueError("cumulative grid has an invalid bin count")
    points = total + 1
    result = []
    for index in range(bins):
        low = index * points // bins
        high = (index + 1) * points // bins - 1
        result.append((low, high))
    if result[0][0] != 0 or result[-1][1] != total:
        raise RuntimeError("balanced cumulative bins lost an endpoint")
    if any(left[1] + 1 != right[0] for left, right in zip(result, result[1:])):
        raise RuntimeError("balanced cumulative bins are not contiguous")
    return tuple(result)


def shifted_cumulatives(profile: Iterable[int], support: tuple[int, ...]) -> tuple[int, ...]:
    values = tuple(int(value) for value in profile)
    running = 0
    result = []
    for index in support[:-1]:
        running += values[index] - 1
        result.append(running)
    return tuple(result)


def bin_owner(cumulatives: tuple[int, ...], intervals: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    result = []
    for value in cumulatives:
        owner = next(
            (index for index, (low, high) in enumerate(intervals) if low <= value <= high),
            None,
        )
        if owner is None:
            raise RuntimeError("cumulative value has no grid owner")
        result.append(owner)
    if any(left > right for left, right in zip(result, result[1:])):
        raise RuntimeError("cumulative bin owner is not nondecreasing")
    return tuple(result)


def runs(indices: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    result = []
    start = 0
    while start < len(indices):
        end = start + 1
        while end < len(indices) and indices[end] == indices[start]:
            end += 1
        result.append((start, end, indices[start]))
        start = end
    return tuple(result)


def unconstrained_cell_count(
    indices: tuple[int, ...], intervals: tuple[tuple[int, int], ...]
) -> int:
    """Count nondecreasing cumulative sequences in one bin tuple."""

    result = 1
    for start, end, bin_index in runs(indices):
        width = intervals[bin_index][1] - intervals[bin_index][0] + 1
        length = end - start
        result *= math.comb(width + length - 1, length)
    return result


def profile_from_cumulatives(
    cumulatives: tuple[int, ...], support: tuple[int, ...], free_mass: int
) -> tuple[int, ...]:
    shifted = []
    previous = 0
    for value in cumulatives:
        shifted.append(value - previous)
        previous = value
    shifted.append(free_mass - previous)
    profile = [0] * CLASSES
    for index, value in zip(support, shifted):
        profile[index] = value + 1
    return tuple(profile)


def explicit_vertices(
    support: tuple[int, ...],
    indices: tuple[int, ...],
    intervals: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, ...], ...]:
    """Return product-of-order-simplex vertices when the weight cut is redundant."""

    choices = []
    for start, end, bin_index in runs(indices):
        low, high = intervals[bin_index]
        length = end - start
        choices.append(
            tuple((low,) * split + (high,) * (length - split) for split in range(length + 1))
        )
    free_mass = MASS - len(support)
    vertices = set()
    for selected in itertools.product(*choices):
        cumulative = tuple(value for block in selected for value in block)
        profile = profile_from_cumulatives(cumulative, support, free_mass)
        if min(profile[index] for index in support) < 1:
            raise RuntimeError("explicit cumulative vertex left exact support")
        vertices.add(profile)
    return tuple(sorted(vertices))


def weight_cut_redundant(support: tuple[int, ...]) -> bool:
    if 0 not in support:
        return True
    return sum(support) >= CUTOFF


def clipped_vertices(
    support: tuple[int, ...],
    indices: tuple[int, ...],
    intervals: tuple[tuple[int, int], ...],
    maximum_systems: int,
) -> tuple[Profile, ...]:
    """Use the verifier's exact chart when the physical-weight cut is active."""

    constraints = list(verifier.root_constraints(support, MASS, CUTOFF))
    dimension = len(support) - 1
    for coordinate in range(dimension):
        low, high = intervals[indices[coordinate]]
        # y_k=sum_{i<=k}(a_i-1).
        left = [Fraction(0)] * dimension
        for index in range(coordinate + 1):
            left[index] = 1
        offset = coordinate + 1
        constraints.append((tuple(left), Fraction(high + offset)))
        constraints.append((tuple(-value for value in left), Fraction(-(low + offset))))
    return verifier.enumerate_vertices(
        support, tuple(constraints), maximum_systems
    )


def excluded_profiles(support: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Enumerate the globally tiny weight-below-21 part of one support."""

    if 0 not in support:
        return ()
    positive = tuple(index for index in support if index)
    profile = [0] * CLASSES
    result = []

    def visit(position: int, weight: int, used: int) -> None:
        if position == len(positive):
            profile[0] = MASS - used
            if profile[0] >= 1:
                result.append(tuple(profile))
            profile[0] = 0
            return
        index = positive[position]
        for count in range(1, (CUTOFF - 1 - weight) // index + 1):
            profile[index] = count
            visit(position + 1, weight + index * count, used + count)
        profile[index] = 0

    visit(0, 0, 0)
    return tuple(result)


def exclusion_counts_by_cell(
    support: tuple[int, ...], intervals: tuple[tuple[int, int], ...]
) -> Counter:
    result: Counter = Counter()
    for profile in excluded_profiles(support):
        owner = bin_owner(shifted_cumulatives(profile, support), intervals)
        result[owner] += 1
    return result


def prefix_split(support: tuple[int, ...], coordinate: int, upper: int) -> dict:
    """Return the canonical split y_coordinate<=upper in profile coordinates."""

    coefficients = [0] * CLASSES
    for index in support[: coordinate + 1]:
        coefficients[index] = 1
    threshold = upper + coordinate + 1
    parsed = verifier.validate_split(
        {"coefficients": coefficients, "threshold": str(threshold)}, support
    )
    if parsed != (tuple(coefficients), threshold):
        raise RuntimeError("cumulative prefix split is not verifier-canonical")
    return {"coefficients": coefficients, "threshold": str(threshold)}


def scaling_row(dimension: int, levels: int) -> dict:
    bins = 1 << levels
    cells = math.comb(bins + dimension - 1, dimension)
    # Sum and maximize product(run_length+1) over all nondecreasing bin tuples.
    total_vertices = 0
    maximum_vertices = 0
    for groups in range(1, min(bins, dimension) + 1):
        for composition in positive_compositions(dimension, groups):
            product = math.prod(length + 1 for length in composition)
            placements = math.comb(bins, groups)
            total_vertices += placements * product
            maximum_vertices = max(maximum_vertices, product)
    if dimension == 0:
        total_vertices = maximum_vertices = 1
    return {
        "levels": levels,
        "bins": bins,
        "nonempty_cells": cells,
        "maximum_vertices_per_cell": maximum_vertices,
        "total_vertex_incidences": total_vertices,
    }


def positive_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(1, total - parts + 2):
        for tail in positive_compositions(total - first, parts - 1):
            yield (first, *tail)


def build(support: tuple[int, ...], levels: int, maximum_systems: int) -> dict:
    if not support or support == (0,):
        raise ValueError("cumulative prototype needs a feasible support")
    dimension = len(support) - 1
    bins = 1 << levels
    free_mass = MASS - len(support)
    intervals = balanced_bins(free_mass, bins)
    exclusions = exclusion_counts_by_cell(support, intervals)
    cells = []
    for indices in itertools.combinations_with_replacement(range(bins), dimension):
        raw_count = unconstrained_cell_count(indices, intervals)
        count = raw_count - exclusions[indices]
        if count < 0:
            raise RuntimeError("physical-weight exclusions exceeded a cell count")
        if weight_cut_redundant(support):
            vertices = tuple(
                tuple(Fraction(value) for value in profile)
                for profile in explicit_vertices(support, indices, intervals)
            )
            geometry = "explicit-product-of-chain-simplices-v1"
        else:
            if dimension > 2:
                raise ValueError(
                    "prototype restricts active physical cuts to dimension at most two"
                )
            vertices = clipped_vertices(support, indices, intervals, maximum_systems)
            geometry = "exact-rational-verifier-intersection-v1"
        if count and not vertices:
            raise RuntimeError("nonempty integer cumulative cell has no real vertices")
        cells.append(
            {
                "bin_indices": list(indices),
                "bin_intervals": [list(intervals[index]) for index in indices],
                "exact_count": str(count),
                "unconstrained_count": str(raw_count),
                "weight_cut_exclusions": str(exclusions[indices]),
                "geometry": geometry,
                "vertices": [
                    [fraction_text(value) for value in profile] for profile in vertices
                ],
            }
        )
    expected, expected_excluded = verifier.exact_support_count(support, MASS, CUTOFF)
    if sum(int(row["exact_count"]) for row in cells) != expected:
        raise RuntimeError("cumulative cells do not sum to the exact root count")
    if sum(int(row["weight_cut_exclusions"]) for row in cells) != expected_excluded:
        raise RuntimeError("cumulative cells do not recover the weight exclusions")
    prefix_splits = [
        prefix_split(support, coordinate, intervals[boundary][1])
        for coordinate in range(dimension)
        for boundary in range(bins - 1)
    ]
    return {
        "schema": SCHEMA,
        "status": "STRUCTURAL_PROTOTYPE_NO_WITNESS_CLAIM",
        "manifest_sha256": MANIFEST_SHA256,
        "support_mask": verifier.mask_text(sum(1 << index for index in support)),
        "active_classes": list(support),
        "dimension": dimension,
        "levels": levels,
        "bins": bins,
        "free_mass": str(free_mass),
        "integer_bin_intervals": [list(interval) for interval in intervals],
        "ownership": (
            "unique balanced bin of each shifted cumulative mass; bin tuple is "
            "nondecreasing"
        ),
        "canonical_prefix_splits": prefix_splits,
        "exact_root_count": str(expected),
        "weight_cut_exclusions": str(expected_excluded),
        "cells": cells,
        "summary": {
            "nonempty_cells": len(cells),
            "maximum_vertices_per_cell": max(len(row["vertices"]) for row in cells),
            "total_vertex_incidences": sum(len(row["vertices"]) for row in cells),
            "scaling_full_support_dimension_8": [
                scaling_row(8, local_levels) for local_levels in range(1, 5)
            ],
        },
        "scope_limit": (
            "geometry and exact counting only; no fixed witness, outward replay, "
            "shard completion, or probability claim"
        ),
    }


def run_self_test() -> None:
    if file_sha256(DEFAULT_MANIFEST) != MANIFEST_SHA256:
        raise SystemExit("cumulative prototype: manifest digest changed")
    line = build((0, 1), 2, 1000)
    triangle = build((0, 1, 2), 2, 10000)
    full = build(tuple(range(9)), 1, 10000)
    if (
        len(line["cells"]) != 4
        or sum(int(row["exact_count"]) for row in line["cells"]) != 262123
        or len(triangle["cells"]) != 10
        or int(triangle["exact_root_count"]) != 34359345063
        or len(full["cells"]) != 9
        or full["summary"]["maximum_vertices_per_cell"] != 25
        or sum(int(row["exact_count"]) for row in full["cells"])
        != math.comb(MASS - 1, 8)
    ):
        raise SystemExit("cumulative prototype self-test failed")
    intervals = balanced_bins(MASS - 2, 4)
    boundary = intervals[1][1] + 1
    owner = bin_owner((boundary,), intervals)
    if owner != (2,):
        raise SystemExit("cumulative boundary ownership self-test failed")
    # Exhaust a toy chain and compare each cell with the run-product formula.
    toy_intervals = balanced_bins(9, 4)
    for dimension in (1, 2, 3):
        observed: Counter = Counter()
        for chain in itertools.combinations_with_replacement(range(10), dimension):
            observed[bin_owner(chain, toy_intervals)] += 1
        expected = {
            indices: unconstrained_cell_count(indices, toy_intervals)
            for indices in itertools.combinations_with_replacement(range(4), dimension)
        }
        if dict(observed) != expected or sum(observed.values()) != math.comb(9 + dimension, dimension):
            raise SystemExit("cumulative toy-chain count self-test failed")
    print("dimension1_clipped_cells=4")
    print("dimension2_clipped_cells=10")
    print("full_support_level1_cells=9")
    print("full_support_level1_max_vertices=25")
    print("exact_count_conservation=PASS")
    print("toy_chain_product_counts_d1_d3=PASS")
    print("canonical_prefix_splits=PASS")
    print("status=PASS_G8_CUMULATIVE_BOX_PROTOTYPE_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--support-mask", type=lambda value: int(value, 0))
    parser.add_argument("--levels", type=int, default=2)
    parser.add_argument("--maximum-vertex-systems", type=int, default=10000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.support_mask is None or args.output is None:
        parser.error("--support-mask and --output are required")
    if file_sha256(args.manifest) != MANIFEST_SHA256:
        parser.error("manifest digest does not match the frozen structured-geometry input")
    if args.levels < 0 or args.levels > 8 or args.maximum_vertex_systems < 1:
        parser.error("invalid cumulative prototype limit")
    support = support_of(args.support_mask)
    report = build(support, args.levels, args.maximum_vertex_systems)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"mask={report['support_mask']} cells={report['summary']['nonempty_cells']} "
        f"max_vertices={report['summary']['maximum_vertices_per_cell']} "
        f"vertex_incidences={report['summary']['total_vertex_incidences']}"
    )
    print(f"output={args.output}")
    print("status=STRUCTURAL_PROTOTYPE_NO_WITNESS_CLAIM")


if __name__ == "__main__":
    main()
