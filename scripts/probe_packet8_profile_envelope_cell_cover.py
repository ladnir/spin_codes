#!/usr/bin/env python3
"""Branch-and-bound profile cover using secants and the full witness envelope.

The packet-profile normalization is separable and concave.  On an integer box,
the secant of each coordinate term is a rigorous affine lower bound.  Replacing
the normalization by those secants turns the maximum of the pointwise minimum
of every fixed witness into a ten-variable LP.  A cell is covered whenever
that LP upper bound is below the per-profile target.

Children ``[lo,mid]`` and ``[mid+1,hi]`` partition integer profiles exactly;
the LP safely relaxes each child to the corresponding real box.  This script
uses binary64/HiGHS and is a diagnostic precursor to an outward certificate.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.special import gammaln

from probe_packet8_profile_ordered_chambers import add_full_bijection
from probe_packet8_profile_simplex_landscape import M, PROFILE_COUNT_LOG2
from probe_packet8_shared_witness_io import load_shared_witness_arrays
from probe_packet8_weight_profile_scalar import CLASSES


LOG2 = math.log(2.0)
LOG_CLASSES = np.log(np.asarray(CLASSES, dtype=np.float64))
NORMALIZATION_CONSTANT = float(gammaln(M + 1)) / LOG2
WEIGHTS = np.arange(9, dtype=np.float64)


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
    maximizer: tuple[float, ...] = field(compare=False)
    secant_gaps: tuple[float, ...] = field(compare=False)


def coordinate_term(index: int, value: np.ndarray | float):
    return (-gammaln(np.asarray(value) + 1.0) + np.asarray(value) * LOG_CLASSES[index]) / LOG2


def secant(cell: Cell) -> tuple[float, np.ndarray]:
    slopes = np.empty(9)
    intercept = NORMALIZATION_CONSTANT
    for index, (low, high) in enumerate(zip(cell.lower, cell.upper)):
        low_value = float(coordinate_term(index, low))
        if low == high:
            slopes[index] = 0.0
            intercept += low_value
        else:
            high_value = float(coordinate_term(index, high))
            slopes[index] = (high_value - low_value) / (high - low)
            intercept += low_value - slopes[index] * low
    return intercept, slopes


def score_cell(
    cell: Cell,
    constants: np.ndarray,
    charges: np.ndarray,
    minimum_outer_weight: int,
):
    lower = np.asarray(cell.lower, dtype=np.float64)
    upper = np.asarray(cell.upper, dtype=np.float64)
    if np.sum(lower) > M or np.sum(upper) < M:
        return None
    if WEIGHTS @ upper < minimum_outer_weight:
        return None

    intercept, slopes = secant(cell)
    # Variables are the nine profile counts followed by the affine-envelope t.
    # Each fixed witness imposes t <= C_i - <a,q_i>.
    objective = np.append(slopes, -1.0)
    inequalities = np.hstack((charges, np.ones((len(constants), 1))))
    weight_row = np.append(-WEIGHTS, 0.0)[None, :]
    result = linprog(
        objective,
        A_ub=np.vstack((inequalities, weight_row)),
        b_ub=np.append(constants, -minimum_outer_weight),
        A_eq=np.asarray([[1.0] * 9 + [0.0]]),
        b_eq=np.asarray([M]),
        bounds=[(float(low), float(high)) for low, high in zip(lower, upper)]
        + [(None, None)],
        method="highs",
    )
    if result.status == 2:
        return None
    if not result.success:
        raise RuntimeError(result.message)
    point = result.x[:9]
    score = -float(result.fun) - intercept
    exact_terms = np.asarray(
        [coordinate_term(index, point[index]) for index in range(9)]
    )
    secant_terms = np.asarray(
        [
            coordinate_term(index, lower[index])
            + slopes[index] * (point[index] - lower[index])
            for index in range(9)
        ]
    )
    gaps = exact_terms - secant_terms
    return score, tuple(float(value) for value in point), tuple(float(value) for value in gaps)


def choose_split(item: QueueItem) -> int | None:
    cell = item.cell
    widths = np.asarray(cell.upper) - np.asarray(cell.lower)
    candidates = np.flatnonzero(widths > 0)
    if not len(candidates):
        return None
    zero_candidates = [
        int(index)
        for index in candidates
        if cell.lower[int(index)] == 0 < cell.upper[int(index)]
    ]
    if zero_candidates:
        return max(zero_candidates, key=lambda index: widths[index])
    return max(
        (int(index) for index in candidates),
        key=lambda index: (item.secant_gaps[index], widths[index]),
    )


def split_cell(item: QueueItem, coordinate: int) -> tuple[Cell, Cell]:
    cell = item.cell
    low = cell.lower[coordinate]
    high = cell.upper[coordinate]
    if low == 0:
        midpoint = 0
    else:
        candidate = int(math.floor(item.maximizer[coordinate]))
        midpoint = min(max(candidate, low), high - 1)
        # Avoid repeatedly shaving one count from a huge interval when the LP
        # maximizer lies numerically on a boundary.
        if midpoint == low and high - low > 16:
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
    parser.add_argument("--upgraded-cache", type=Path, required=True)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    parser.add_argument("--minimum-outer-weight", type=int, default=21)
    parser.add_argument("--max-nodes", type=int, default=10000)
    parser.add_argument("--max-depth", type=int, default=256)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--show-unresolved", type=int, default=20)
    parser.add_argument("--queue-order", choices=("worst", "depth"), default="worst")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    names, constants, charges = load_shared_witness_arrays(
        args.upgraded_cache, args.extra_shared_report
    )
    full = add_full_bijection([])[0]
    names.append(full.name)
    constants = np.append(constants, full.constant_log2)
    charges = np.vstack((charges, full.linear_charge))
    target = -40.0 - PROFILE_COUNT_LOG2

    root = Cell((0,) * 9, (M,) * 9)
    root_score = score_cell(root, constants, charges, args.minimum_outer_weight)
    if root_score is None:
        raise SystemExit("envelope cell cover: infeasible root")
    serial = 0
    def priority(cell: Cell, score: float) -> float:
        return -float(cell.depth) if args.queue_order == "depth" else -score

    queue = [
        QueueItem(
            priority(root, root_score[0]),
            serial,
            root,
            root_score[0],
            root_score[1],
            root_score[2],
        )
    ]
    covered = 0
    infeasible = 0
    leaves: list[QueueItem] = []
    processed = 0
    print("packet-8 secant/envelope integer-cell cover", flush=True)
    print(
        f"witnesses={len(constants)} target={target:.12f} root_score={root_score[0]:.12f}",
        flush=True,
    )
    while queue and processed < args.max_nodes:
        item = heapq.heappop(queue)
        processed += 1
        if item.score <= target:
            covered += 1
            continue
        coordinate = choose_split(item)
        if coordinate is None or item.cell.depth >= args.max_depth:
            leaves.append(item)
            continue
        for child in split_cell(item, coordinate):
            row = score_cell(child, constants, charges, args.minimum_outer_weight)
            if row is None:
                infeasible += 1
                continue
            serial += 1
            child_item = QueueItem(
                priority(child, row[0]), serial, child, row[0], row[1], row[2]
            )
            if row[0] <= target:
                covered += 1
            else:
                heapq.heappush(queue, child_item)
        if args.progress_every and processed % args.progress_every == 0:
            worst = max((row.score for row in queue), default=-math.inf)
            print(
                f"processed={processed} covered={covered} open={len(queue)} "
                f"leaves={len(leaves)} worst_open={worst:.9f}",
                flush=True,
            )

    leaves.extend(queue)
    leaves.sort(key=lambda row: row.score, reverse=True)
    report = {
        "status": "DIAGNOSTIC_BINARY64_SECANT_ENVELOPE_INTEGER_CELL_COVER",
        "target_log2": target,
        "processed_nodes": processed,
        "covered_cells": covered,
        "infeasible_cells": infeasible,
        "unresolved_cells": len(leaves),
        "complete_cover": not leaves,
        "worst_unresolved": [
            {
                "score_log2": row.score,
                "depth": row.cell.depth,
                "lower": list(row.cell.lower),
                "upper": list(row.cell.upper),
                "maximizer": list(row.maximizer),
                "secant_gaps": list(row.secant_gaps),
            }
            for row in leaves[: args.show_unresolved]
        ],
    }
    print(f"processed_nodes={processed}")
    print(f"covered_cells={covered}")
    print(f"infeasible_cells={infeasible}")
    print(f"unresolved_cells={len(leaves)}")
    if leaves:
        print(f"worst_unresolved_score={leaves[0].score:.12f}")
    for row in leaves[: args.show_unresolved]:
        print(
            f"score={row.score:.9f} depth={row.cell.depth} "
            f"lower={','.join(map(str,row.cell.lower))} "
            f"upper={','.join(map(str,row.cell.upper))}"
        )
    print("complete_envelope_cell_cover=" + ("PASS" if not leaves else "NO"))
    print("status=DIAGNOSTIC_BINARY64_SECANT_ENVELOPE_INTEGER_CELL_COVER")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
