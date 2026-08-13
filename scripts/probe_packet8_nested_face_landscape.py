#!/usr/bin/env python3
"""Audit all nested-prefix triangle faces under the fixed witness atlas.

An ordered packet-profile chamber is an eight-simplex with uniform-prefix
vertices.  Its distinct two-faces correspond to strict chains
``S subset T subset U`` of nonempty packet-weight-class subsets.  This probe
exhausts those chains (excluding the forbidden pure-zero vertex), counts
triangles covered by one convex witness at all three vertices, and evaluates
every triangle barycenter under the complete fixed atlas plus the convex
inverse-orbit branch.

The enumeration is complete and sample-independent, but barycenter coverage
does not certify a whole triangle.  Arithmetic remains binary64.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from probe_packet8_adaptive_simplex_atlas import load_witness_cache
from probe_packet8_profile_edge_interval_cover import orbit_value
from probe_packet8_profile_ordered_chambers import (
    add_full_bijection,
    point_values,
    safe_mask,
    uniform_subset,
)
from probe_packet8_profile_simplex_landscape import (
    ANCHORS,
    PROFILE_COUNT_LOG2,
    integral_profile,
    load_anchor_file,
    normalization_log2,
)


ALL_CLASSES = (1 << 9) - 1


def nested_triples() -> list[tuple[int, int, int]]:
    result = []
    for left in range(2, 1 << 9):
        left_complement = ALL_CLASSES ^ left
        middle_addition = left_complement
        while middle_addition:
            middle = left | middle_addition
            right_complement = ALL_CLASSES ^ middle
            right_addition = right_complement
            while right_addition:
                result.append((left, middle, middle | right_addition))
                right_addition = (right_addition - 1) & right_complement
            middle_addition = (middle_addition - 1) & left_complement
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--show-worst", type=int, default=20)
    args = parser.parse_args()
    if args.batch_size <= 0 or args.show_worst <= 0:
        raise SystemExit("nested face landscape: invalid batch or row limit")

    anchors = ANCHORS + load_anchor_file(args.atlas)
    witnesses = load_witness_cache(
        args.atlas.with_suffix(".witnesses.npz"), anchors
    )
    if witnesses is None:
        raise SystemExit("nested face landscape: matching witness cache missing")
    witnesses = add_full_bijection(witnesses)
    names = [row.name for row in witnesses]
    constants = np.asarray([row.constant_log2 for row in witnesses])
    charges = np.vstack([row.linear_charge for row in witnesses])
    target = -40.0 - PROFILE_COUNT_LOG2

    points = {subset: uniform_subset(subset) for subset in range(2, 1 << 9)}
    values = {
        subset: point_values(point, constants, charges)
        for subset, point in points.items()
    }
    masks = {subset: safe_mask(row, target) for subset, row in values.items()}
    orbit_safe = {
        subset: orbit_value(0.0, point, np.zeros(9), target) <= 0.0
        for subset, point in points.items()
    }

    triples = nested_triples()
    common_witness = 0
    common_orbit = 0
    barycenter_covered = 0
    worst_rows = []
    for start in range(0, len(triples), args.batch_size):
        batch = triples[start : start + args.batch_size]
        profiles = np.vstack(
            [
                (points[left] + points[middle] + points[right]) / 3.0
                for left, middle, right in batch
            ]
        )
        normalizations = normalization_log2(profiles)
        best = np.full(len(batch), np.inf)
        leaders = np.full(len(batch), -1, dtype=np.int64)
        for witness_start in range(0, len(witnesses), 128):
            witness_end = min(witness_start + 128, len(witnesses))
            candidates = (
                constants[None, witness_start:witness_end]
                - profiles @ charges[witness_start:witness_end].T
                - normalizations[:, None]
            )
            local = np.argmin(candidates, axis=1)
            local_values = candidates[np.arange(len(batch)), local]
            improve = local_values < best
            best[improve] = local_values[improve]
            leaders[improve] = witness_start + local[improve]

        orbit = np.asarray(
            [orbit_value(0.0, point, np.zeros(9), 0.0) for point in profiles]
        )
        improve = orbit < best
        best[improve] = orbit[improve]
        leaders[improve] = -2
        barycenter_covered += int(np.sum(best <= target))

        for index, (left, middle, right) in enumerate(batch):
            shared = masks[left] & masks[middle] & masks[right]
            if shared:
                common_witness += 1
            elif orbit_safe[left] and orbit_safe[middle] and orbit_safe[right]:
                common_orbit += 1
            row = (
                float(best[index]),
                left,
                middle,
                right,
                int(leaders[index]),
            )
            if len(worst_rows) < args.show_worst:
                worst_rows.append(row)
                worst_rows.sort(reverse=True)
            elif row[0] > worst_rows[-1][0]:
                worst_rows[-1] = row
                worst_rows.sort(reverse=True)

    total = len(triples)
    print("packet-8 nested-prefix triangle-face landscape")
    print(
        f"atlas={args.atlas} witnesses={len(witnesses)} "
        f"target={target:.12f}"
    )
    print(f"total_nested_triangles={total}")
    print(f"single_fixed_witness_triangles={common_witness}")
    print(f"single_inverse_orbit_triangles={common_orbit}")
    print(
        f"single_branch_triangles={common_witness+common_orbit} "
        f"multi_branch_triangles={total-common_witness-common_orbit}"
    )
    print(f"covered_barycenters={barycenter_covered}")
    print(f"uncovered_barycenters={total-barycenter_covered}")
    for rank, (score, left, middle, right, leader) in enumerate(worst_rows, 1):
        leader_name = names[leader] if leader >= 0 else "inverse_orbit"
        print(
            f"rank={rank} score={score:.9f} "
            f"face={left:03x}<{middle:03x}<{right:03x} "
            f"leader={leader_name} profile="
            + ",".join(
                map(
                    str,
                    integral_profile(
                        (points[left] + points[middle] + points[right]) / 3.0
                    ),
                )
            )
        )
    print("status=DIAGNOSTIC_BINARY64_EXHAUSTIVE_TRIANGLE_BARYCENTERS")


if __name__ == "__main__":
    main()
