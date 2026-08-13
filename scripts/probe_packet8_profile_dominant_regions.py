#!/usr/bin/env python3
"""Test the nine dominant-class packet-profile polytopes.

The simplex region in which class ``j`` is at least every other coordinate is
convex.  Its vertices are the uniform profiles on every subset containing
``j``.  Thus it has only ``2^8`` vertices.  For dominant class zero, the pure
zero vertex is removed by the exact outer minimum-weight constraint; clipping
at total bit weight 21 adds the intersections of its incident rays with that
weight boundary.

For each region this diagnostic asks whether one cached fixed witness is below
the per-profile union target at every vertex.  A PASS is sample-independent
for the whole real polytope, apart from current binary64 arithmetic.
"""

from __future__ import annotations

import argparse
import math
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


def dominant_vertices(dominant: int, minimum_weight: int) -> np.ndarray:
    others = [index for index in range(9) if index != dominant]
    rows = []
    zero = np.zeros(9, dtype=np.float64)
    zero[0] = M
    for mask in range(1 << 8):
        support = [dominant] + [
            index for bit, index in enumerate(others) if mask >> bit & 1
        ]
        point = np.zeros(9, dtype=np.float64)
        point[support] = M / len(support)
        weight = float(np.dot(np.arange(9), point))
        if weight >= minimum_weight - 1e-12:
            rows.append(point)
        elif dominant == 0 and len(support) == 1:
            # The pure-zero endpoint itself is the excluded zero message.  Its
            # clipped replacements are added below along every other vertex.
            continue

    if dominant == 0:
        # The pure-zero vertex is adjacent to every nonzero subset vertex in
        # this star representation.  Intersecting all such rays with the
        # minimum-weight plane safely includes every new clipped vertex (and
        # may include redundant boundary points, which only strengthens a PASS).
        base = list(rows)
        for point in base:
            weight = float(np.dot(np.arange(9), point))
            if weight <= minimum_weight:
                continue
            fraction = minimum_weight / weight
            rows.append(zero + fraction * (point - zero))

    unique = {
        tuple(round(float(value), 12) for value in point): point
        for point in rows
    }
    return np.vstack(list(unique.values()))


def add_full_bijection(witnesses: list[FixedWitness]) -> list[FixedWitness]:
    constant = N * math.log2(11) - (N - D) * math.log2(10) + K
    return witnesses + [FixedWitness("full_bijection", np.zeros(9), constant)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--minimum-outer-weight", type=int, default=21)
    args = parser.parse_args()

    extras = load_anchor_file(args.atlas)
    anchors = ANCHORS + extras
    witnesses = load_witness_cache(
        args.atlas.with_suffix(".witnesses.npz"), anchors
    )
    if witnesses is None:
        raise SystemExit("dominant regions: matching witness cache missing")
    witnesses = add_full_bijection(witnesses)
    constants = np.asarray([row.constant_log2 for row in witnesses])
    charges = np.vstack([row.linear_charge for row in witnesses])
    names = [row.name for row in witnesses]
    target = -40.0 - PROFILE_COUNT_LOG2

    print("packet-8 dominant packet-class region probe")
    print(
        f"atlas={args.atlas} witnesses={len(witnesses)} "
        f"target={target:.12f} minimum_outer_weight={args.minimum_outer_weight}"
    )
    passes = 0
    for dominant in range(9):
        vertices = dominant_vertices(dominant, args.minimum_outer_weight)
        normalizations = normalization_log2(vertices)
        values = (
            constants[None, :]
            - vertices @ charges.T
            - normalizations[:, None]
        )
        maxima = np.max(values, axis=0)
        leader_index = int(np.argmin(maxima))
        worst_index = int(np.argmax(values[:, leader_index]))
        score = float(maxima[leader_index])
        status = "PASS" if score <= target else "NO"
        passes += status == "PASS"
        print(
            f"dominant_class={dominant} vertices={len(vertices)} "
            f"score={score:.12f} leader={names[leader_index]} status={status} "
            f"worst_vertex="
            + ",".join(f"{value:.6f}" for value in vertices[worst_index])
        )
    print(f"covered_dominant_regions={passes}/9")
    print(f"complete_dominant_partition={'PASS' if passes == 9 else 'NO'}")
    print("status=DIAGNOSTIC_BINARY64_DOMINANT_REGION_VERTEX_CHECK")


if __name__ == "__main__":
    main()
