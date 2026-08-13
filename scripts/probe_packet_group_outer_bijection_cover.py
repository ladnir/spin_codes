#!/usr/bin/env python3
"""Integer-cell cover by outer-profile and all-message bijection branches."""

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

from packet_group_drive_stratified import profile_classes, profile_count
from packet_group_outer_profile import D, K, N, atom_count


LOG2 = math.log(2.0)
MINIMUM_WEIGHT = 21


@dataclass(frozen=True)
class Cell:
    lower: tuple[int, ...]
    upper: tuple[int, ...]
    depth: int = 0


@dataclass(order=True)
class Item:
    priority: float
    serial: int
    cell: Cell = field(compare=False)
    score: float = field(compare=False)
    maximizer: tuple[float, ...] = field(compare=False)
    gaps: tuple[float, ...] = field(compare=False)


def term(classes: np.ndarray, index: int, value):
    value = np.asarray(value)
    return (-gammaln(value + 1.0) + value * math.log(classes[index])) / LOG2


def secant(group_bits: int, cell: Cell):
    classes = np.asarray(profile_classes(group_bits), dtype=np.float64)
    intercept = float(gammaln(atom_count(group_bits) + 1)) / LOG2
    slopes = np.zeros(group_bits + 1)
    for index, (low, high) in enumerate(zip(cell.lower, cell.upper)):
        low_value = float(term(classes, index, low))
        if low == high:
            intercept += low_value
        else:
            high_value = float(term(classes, index, high))
            slopes[index] = (high_value - low_value) / (high - low)
            intercept += low_value - slopes[index] * low
    return intercept, slopes


def score_cell(
    group_bits: int,
    cell: Cell,
    outer_constants: np.ndarray,
    outer_charges: np.ndarray,
    combined_rows: list[dict],
):
    atoms = atom_count(group_bits)
    lower = np.asarray(cell.lower, dtype=np.float64)
    upper = np.asarray(cell.upper, dtype=np.float64)
    weights = np.arange(group_bits + 1, dtype=np.float64)
    if np.sum(lower) > atoms or np.sum(upper) < atoms:
        return None
    if weights @ upper < MINIMUM_WEIGHT:
        return None
    intercept, slopes = secant(group_bits, cell)
    bijection_constant = K + N * math.log2(11) - (N - D) * math.log2(10)
    # Variables are the profile followed by the common upper-envelope value t.
    inequality_rows = [
        np.hstack((outer_charges, np.ones((len(outer_constants), 1))))
    ]
    inequality_bounds = [outer_constants]
    eligible_combined = [
        row
        for row in combined_rows
        if all(
            index in row["support"] or cell.upper[index] == 0
            for index in range(group_bits + 1)
        )
    ]
    if eligible_combined:
        combined_charges = np.asarray([row["charge"] for row in eligible_combined])
        inequality_rows.append(
            np.hstack(
                (
                    combined_charges + slopes[None, :],
                    np.ones((len(eligible_combined), 1)),
                )
            )
        )
        inequality_bounds.append(
            np.asarray([row["constant_log2"] for row in eligible_combined])
            - intercept
        )
    inequality_rows.extend(
        (
            np.append(slopes, 1.0)[None, :],
            np.append(-weights, 0.0)[None, :],
        )
    )
    inequality_bounds.extend(
        (np.asarray([bijection_constant - intercept]), np.asarray([-MINIMUM_WEIGHT]))
    )
    rows = np.vstack(inequality_rows)
    bounds = [(float(low), float(high)) for low, high in zip(lower, upper)]
    result = linprog(
        np.append(np.zeros(group_bits + 1), -1.0),
        A_ub=rows,
        b_ub=np.concatenate(inequality_bounds),
        A_eq=np.asarray([[1.0] * (group_bits + 1) + [0.0]]),
        b_eq=np.asarray([atoms]),
        bounds=bounds + [(None, None)],
        method="highs",
    )
    if result.status == 2:
        return None
    if not result.success:
        raise RuntimeError(result.message)
    point = result.x[: group_bits + 1]
    classes = np.asarray(profile_classes(group_bits), dtype=np.float64)
    gaps = []
    for index in range(group_bits + 1):
        low = lower[index]
        secant_value = term(classes, index, low) + slopes[index] * (point[index] - low)
        gaps.append(float(term(classes, index, point[index]) - secant_value))
    return -float(result.fun), tuple(map(float, point)), tuple(gaps)


def choose_split(item: Item) -> int | None:
    widths = np.asarray(item.cell.upper) - np.asarray(item.cell.lower)
    candidates = np.flatnonzero(widths > 0)
    if not len(candidates):
        return None
    zero = [int(index) for index in candidates if item.cell.lower[int(index)] == 0]
    if zero:
        return max(zero, key=lambda index: widths[index])
    return max(
        (int(index) for index in candidates),
        key=lambda index: (item.gaps[index], widths[index]),
    )


def children(item: Item, coordinate: int):
    low = item.cell.lower[coordinate]
    high = item.cell.upper[coordinate]
    if low == 0:
        midpoint = 0
    else:
        candidate = int(math.floor(item.maximizer[coordinate]))
        midpoint = min(max(candidate, low), high - 1)
        if midpoint == low and high - low > 16:
            midpoint = (low + high) // 2
    left_upper = list(item.cell.upper)
    left_upper[coordinate] = midpoint
    right_lower = list(item.cell.lower)
    right_lower[coordinate] = midpoint + 1
    return (
        Cell(item.cell.lower, tuple(left_upper), item.cell.depth + 1),
        Cell(tuple(right_lower), item.cell.upper, item.cell.depth + 1),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--combined-atlas", type=Path, action="append", default=[])
    parser.add_argument("--max-nodes", type=int, default=10000)
    parser.add_argument("--max-depth", type=int, default=256)
    parser.add_argument("--queue-order", choices=("depth", "worst"), default="depth")
    parser.add_argument("--progress-every", type=int, default=250)
    parser.add_argument("--show-unresolved", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    atlas = json.loads(args.atlas.read_text(encoding="utf-8"))
    g = int(atlas["group_bits"])
    constants = np.asarray([row["constant_log2"] for row in atlas["rows"]])
    charges = np.asarray([row["charge"] for row in atlas["rows"]])
    combined_rows = []
    for path in args.combined_atlas:
        combined = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(combined, list):
            combined_group = int(combined[0]["group_bits"]) if combined else g
            rows = combined
        else:
            combined_group = int(combined["group_bits"])
            rows = combined["rows"]
        if combined_group != g:
            raise SystemExit("outer/bijection cover: combined atlas group mismatch")
        combined_rows.extend(rows)
    target = -40.0 - math.log2(profile_count(g, N))
    atoms = atom_count(g)
    root = Cell((0,) * (g + 1), (atoms,) * (g + 1))
    root_row = score_cell(g, root, constants, charges, combined_rows)
    if root_row is None:
        raise SystemExit("outer/bijection cover: infeasible root")

    def priority(cell: Cell, score: float):
        return -cell.depth if args.queue_order == "depth" else -score

    serial = 0
    queue = [Item(priority(root, root_row[0]), serial, root, *root_row)]
    covered = 0
    infeasible = 0
    leaves = []
    processed = 0
    print("generalized packet outer/bijection cell cover", flush=True)
    print(
        f"g={g} outer_witnesses={len(constants)} combined_witnesses={len(combined_rows)} target={target:.12f} "
        f"root_score={root_row[0]:.12f}",
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
        for child in children(item, coordinate):
            row = score_cell(g, child, constants, charges, combined_rows)
            if row is None:
                infeasible += 1
                continue
            serial += 1
            child_item = Item(priority(child, row[0]), serial, child, *row)
            if row[0] <= target:
                covered += 1
            else:
                heapq.heappush(queue, child_item)
        if args.progress_every and processed % args.progress_every == 0:
            print(
                f"processed={processed} covered={covered} open={len(queue)} "
                f"leaves={len(leaves)}",
                flush=True,
            )
    leaves.extend(queue)
    leaves.sort(key=lambda row: row.score, reverse=True)
    report = {
        "status": "DIAGNOSTIC_BINARY64_GENERALIZED_OUTER_BIJECTION_CELL_COVER",
        "group_bits": g,
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
    print("complete_outer_bijection_cover=" + ("PASS" if not leaves else "NO"))
    print("status=DIAGNOSTIC_BINARY64_GENERALIZED_OUTER_BIJECTION_CELL_COVER")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
