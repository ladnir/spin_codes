#!/usr/bin/env python3
"""Branch-and-bound convex-cell cover for packet-weight profiles.

For every fixed outer/inner witness the logged profile bound is

    C - <a,q> - log2 Q(a).

Since ``-log Q(a)`` is convex on the real simplex, one fixed witness covers a
box/simplex cell whenever it is below the target at every vertex of that
cell.  This diagnostic recursively splits unresolved integer-bound cells and
reports their exact box bounds.  Unlike random landscape sampling, a reported
covered cell contains every real profile in that polytope and hence every
integer profile in it (subject to current binary64 arithmetic).

The root also imposes total binary weight at least 21, the exact minimum
allowed after one possible outer puncture.  Vertices are enumerated both with
the weight inequality slack and on its active boundary.

This is not yet an outward certificate: fixed witness constants, gamma
evaluation, and vertex comparisons are binary64 diagnostics.
"""

from __future__ import annotations

import argparse
import heapq
import math
from dataclasses import dataclass, field
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


MINIMUM_OUTER_WEIGHT = 21


@dataclass(frozen=True)
class Cell:
    lower: tuple[int, ...]
    upper: tuple[int, ...]
    depth: int = 0


@dataclass(order=True)
class QueueItem:
    priority: float
    serial: int
    cell: Cell = field(compare=False)
    score: float = field(compare=False)
    leader: str = field(compare=False)
    vertices: int = field(compare=False)
    worst_vertex: tuple[float, ...] = field(compare=False)


def cell_vertices(cell: Cell, minimum_weight: int) -> np.ndarray:
    lower = np.asarray(cell.lower, dtype=np.float64)
    upper = np.asarray(cell.upper, dtype=np.float64)
    rows: dict[tuple[float, ...], np.ndarray] = {}

    # Vertices at which only box bounds are active: eight coordinates are at
    # bounds and the ninth is determined by sum a_j=M.
    for free in range(9):
        fixed = [index for index in range(9) if index != free]
        for mask in range(1 << 8):
            point = lower.copy()
            for bit, index in enumerate(fixed):
                if mask >> bit & 1:
                    point[index] = upper[index]
            point[free] = M - float(np.sum(point[fixed]))
            if point[free] < lower[free] - 1e-9 or point[free] > upper[free] + 1e-9:
                continue
            if float(np.dot(np.arange(9), point)) < minimum_weight - 1e-9:
                continue
            key = tuple(round(float(value), 12) for value in point)
            rows[key] = point

    # Vertices on the active total-weight boundary.  Seven coordinates are at
    # bounds; the remaining two solve the sum and weight equalities.
    for first in range(9):
        for second in range(first + 1, 9):
            fixed = [
                index for index in range(9) if index not in (first, second)
            ]
            for mask in range(1 << 7):
                point = lower.copy()
                for bit, index in enumerate(fixed):
                    if mask >> bit & 1:
                        point[index] = upper[index]
                remainder = M - float(np.sum(point[fixed]))
                weight_remainder = minimum_weight - float(
                    np.dot(np.asarray(fixed), point[fixed])
                )
                point[first] = (
                    second * remainder - weight_remainder
                ) / (second - first)
                point[second] = remainder - point[first]
                if (
                    point[first] < lower[first] - 1e-9
                    or point[first] > upper[first] + 1e-9
                    or point[second] < lower[second] - 1e-9
                    or point[second] > upper[second] + 1e-9
                ):
                    continue
                key = tuple(round(float(value), 12) for value in point)
                rows[key] = point

    if not rows:
        return np.empty((0, 9), dtype=np.float64)
    return np.vstack(list(rows.values()))


def add_full_bijection(witnesses: list[FixedWitness]) -> list[FixedWitness]:
    constant = N * math.log2(11) - (N - D) * math.log2(10) + K
    return witnesses + [
        FixedWitness("full_bijection", np.zeros(9), constant)
    ]


def score_cell(
    cell: Cell,
    constants: np.ndarray,
    charges: np.ndarray,
    names: list[str],
    minimum_weight: int,
):
    vertices = cell_vertices(cell, minimum_weight)
    if not len(vertices):
        return None
    normalizations = normalization_log2(vertices)
    candidates = (
        constants[None, :]
        - vertices @ charges.T
        - normalizations[:, None]
    )
    witness_maxima = np.max(candidates, axis=0)
    leader_index = int(np.argmin(witness_maxima))
    worst_index = int(np.argmax(candidates[:, leader_index]))
    return (
        float(witness_maxima[leader_index]),
        names[leader_index],
        len(vertices),
        tuple(float(value) for value in vertices[worst_index]),
    )


def choose_split(
    cell: Cell, rare_cutoff: int, worst_vertex: tuple[float, ...]
) -> int | None:
    widths = np.asarray(cell.upper, dtype=np.int64) - np.asarray(
        cell.lower, dtype=np.int64
    )
    if int(np.max(widths)) <= 0:
        return None
    # Isolate exact support faces before bisecting positive count ranges.  A
    # sharp witness with negligible inactive-class fugacities is useful only
    # after the corresponding coordinates have been fixed to zero; midpoint
    # splitting postpones that event for roughly log2(M) levels per class.
    support_candidates = [
        index
        for index, (low, high) in enumerate(zip(cell.lower, cell.upper))
        if low == 0 and high > 0
    ]
    if support_candidates:
        return max(
            support_candidates,
            key=lambda index: (
                abs(worst_vertex[index] - (cell.lower[index] + cell.upper[index]) / 2),
                widths[index],
            ),
        )
    rare_candidates = [
        index
        for index, (low, high) in enumerate(zip(cell.lower, cell.upper))
        if low == 1 and high > rare_cutoff
    ]
    if rare_candidates:
        return max(
            rare_candidates,
            key=lambda index: (worst_vertex[index], widths[index]),
        )
    candidates = np.flatnonzero(widths > 0)
    return int(
        max(
            candidates,
            key=lambda index: (
                worst_vertex[int(index)],
                widths[int(index)],
            ),
        )
    )


def child_cells(cell: Cell, coordinate: int, rare_cutoff: int):
    low = cell.lower[coordinate]
    high = cell.upper[coordinate]
    if low == 0 < high:
        midpoint = 0
    elif low == 1 and high > rare_cutoff:
        midpoint = rare_cutoff
    else:
        midpoint = (low + high) // 2
    left_upper = list(cell.upper)
    left_upper[coordinate] = midpoint
    right_lower = list(cell.lower)
    right_lower[coordinate] = midpoint + 1
    return (
        Cell(cell.lower, tuple(left_upper), cell.depth + 1),
        Cell(tuple(right_lower), cell.upper, cell.depth + 1),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--max-nodes", type=int, default=10000)
    parser.add_argument("--max-depth", type=int, default=96)
    parser.add_argument("--minimum-outer-weight", type=int, default=MINIMUM_OUTER_WEIGHT)
    parser.add_argument("--show-unresolved", type=int, default=20)
    parser.add_argument("--progress-every", type=int, default=1000)
    parser.add_argument(
        "--queue-order", choices=("depth", "score", "hybrid"), default="hybrid"
    )
    parser.add_argument("--rare-cutoff", type=int, default=1024)
    parser.add_argument("--upgraded-cache", type=Path)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    args = parser.parse_args()
    if args.max_nodes <= 0 or args.max_depth <= 0 or args.show_unresolved <= 0:
        raise SystemExit("convex-cell cover: invalid limit")

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
                "convex-cell cover: --extra-shared-report requires --upgraded-cache"
            )
        extras = load_anchor_file(args.atlas)
        all_anchors = ANCHORS + extras
        cache_path = args.atlas.with_suffix(".witnesses.npz")
        witnesses = load_witness_cache(cache_path, all_anchors)
        if witnesses is None:
            raise SystemExit(
                "convex-cell cover: matching witness cache missing; resume the adaptive atlas once"
            )
        witnesses = add_full_bijection(witnesses)
        constants = np.asarray(
            [fixed.constant_log2 for fixed in witnesses], dtype=np.float64
        )
        charges = np.vstack([fixed.linear_charge for fixed in witnesses])
        names = [fixed.name for fixed in witnesses]
    target = -40.0 - PROFILE_COUNT_LOG2

    root = Cell((0,) * 9, (M,) * 9)
    root_row = score_cell(
        root, constants, charges, names, args.minimum_outer_weight
    )
    if root_row is None:
        raise SystemExit("convex-cell cover: root unexpectedly infeasible")

    serial = 0
    def priority(cell: Cell, score: float) -> float:
        if args.queue_order == "depth":
            return float(cell.depth)
        if args.queue_order == "score":
            return -score
        # Complete the nine zero-versus-positive support decisions first;
        # then pursue the numerically worst rare/common cells.
        return (
            float(cell.depth)
            if cell.depth < 9
            else 9.0 - score / 1e9
        )

    queue = [
        QueueItem(
            priority(root, root_row[0]),
            serial,
            root,
            root_row[0],
            root_row[1],
            root_row[2],
            root_row[3],
        )
    ]
    covered = 0
    infeasible = 0
    leaves: list[QueueItem] = []
    processed = 0
    maximum_queue = 1

    print("packet-8 convex profile-cell cover probe")
    print(
        f"atlas={args.atlas} fixed_witnesses={len(constants)} "
        f"target={target:.12f} minimum_outer_weight={args.minimum_outer_weight}"
    )
    print(
        f"root_score={root_row[0]:.12f} root_leader={root_row[1]} "
        f"root_vertices={root_row[2]}"
    )
    while queue and processed < args.max_nodes:
        item = heapq.heappop(queue)
        processed += 1
        if item.score <= target:
            covered += 1
            continue
        coordinate = choose_split(
            item.cell, args.rare_cutoff, item.worst_vertex
        )
        if coordinate is None or item.cell.depth >= args.max_depth:
            leaves.append(item)
            continue
        for child in child_cells(item.cell, coordinate, args.rare_cutoff):
            row = score_cell(
                child, constants, charges, names, args.minimum_outer_weight
            )
            if row is None:
                infeasible += 1
                continue
            serial += 1
            child_item = QueueItem(
                priority(child, row[0]),
                serial,
                child,
                row[0],
                row[1],
                row[2],
                row[3],
            )
            if row[0] <= target:
                covered += 1
            else:
                heapq.heappush(queue, child_item)
        maximum_queue = max(maximum_queue, len(queue))
        if args.progress_every and processed % args.progress_every == 0:
            worst = max((queued.score for queued in queue), default=-math.inf)
            print(
                f"progress_nodes={processed} covered_cells={covered} "
                f"open_cells={len(queue)} unresolved_leaves={len(leaves)} "
                f"worst_open_score={worst:.9f}",
                flush=True,
            )

    leaves.extend(queue)
    leaves.sort(key=lambda item: item.score, reverse=True)
    print(f"processed_nodes={processed}")
    print(f"covered_cells={covered}")
    print(f"infeasible_cells={infeasible}")
    print(f"maximum_queue={maximum_queue}")
    print(f"unresolved_cells={len(leaves)}")
    if leaves:
        print(f"worst_unresolved_score={leaves[0].score:.12f}")
        print(f"minimum_unresolved_depth={min(item.cell.depth for item in leaves)}")
        print(f"maximum_unresolved_depth={max(item.cell.depth for item in leaves)}")
    for rank, item in enumerate(leaves[: args.show_unresolved], 1):
        print(
            f"rank={rank} score={item.score:.9f} leader={item.leader} "
            f"depth={item.cell.depth} vertices={item.vertices} "
            f"worst_vertex={','.join(f'{value:.6f}' for value in item.worst_vertex)} "
            f"lower={','.join(map(str,item.cell.lower))} "
            f"upper={','.join(map(str,item.cell.upper))}"
        )
    print(
        "sample_independent_cell_cover="
        + ("PASS" if not leaves else "NO")
    )
    print("status=DIAGNOSTIC_BINARY64_CONVEX_CELL_BRANCH_AND_BOUND")


if __name__ == "__main__":
    main()
