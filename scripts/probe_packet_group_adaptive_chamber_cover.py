#!/usr/bin/env python3
"""Adaptive barycentric refinement of generalized ordered profile chambers.

Each ordered chamber is a simplex, except when the all-zero pure vertex is
cut off by physical weight 21; that small truncated chamber is triangulated
once.  A cell passes when one fixed affine witness or fixed convex mixture is
safe at all ``g+1`` vertices.  A failed simplex is split into ``g+1`` children
by its exact barycenter.  Failed terminal cells emit their rounded barycenters
as profiles for the parallel witness tuner.

This is discovery geometry.  The selected final mesh, rational mixtures, and
all witness inequalities require a separate exact/outward verifier.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay
from scipy.special import logsumexp

from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import N, normalization_log2
from probe_packet_group_ordered_chamber_cover import chamber_vertices
from probe_packet_group_support_face_cover import (
    parse_atlas,
    rounded_integer_profile,
    solve_fixed_mixture,
)


Profile = tuple[Fraction, ...]
Cell = tuple[Profile, ...]


def initial_cells(group_bits: int, order: tuple[int, ...]) -> list[Cell]:
    vertices = [profile for _name, profile in chamber_vertices(group_bits, order)]
    dimension = group_bits
    if len(vertices) == dimension + 1:
        return [tuple(vertices)]
    coordinates = np.asarray(
        [[float(value) for value in profile[1:]] for profile in vertices]
    )
    # QJ resolves exact cocircular degeneracies in the one-time diagnostic
    # triangulation.  The final proof artifact must replace this by a rational
    # triangulation and exact volume/boundary checks.
    triangulation = Delaunay(coordinates, qhull_options="Qbb Qc Q12 QJ")
    result = []
    seen = set()
    for indices in triangulation.simplices:
        if any(int(index) >= len(vertices) for index in indices):
            continue
        cell = tuple(vertices[int(index)] for index in indices)
        key = tuple(sorted(cell))
        if key not in seen:
            seen.add(key)
            result.append(cell)
    if not result:
        raise RuntimeError(f"ordered chamber {order} produced no initial simplices")
    return result


def barycenter(cell: Cell) -> Profile:
    count = len(cell)
    return tuple(sum(vertex[index] for vertex in cell) / count for index in range(len(cell[0])))


def subdivide(cell: Cell) -> list[Cell]:
    center = barycenter(cell)
    return [
        tuple(center if index == replaced else vertex for index, vertex in enumerate(cell))
        for replaced in range(len(cell))
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument("--max-depth", type=int, default=0)
    parser.add_argument(
        "--cell-budget-log2",
        type=float,
        default=20.0,
        help="reserve this many union bits for the eventual number of closed cells",
    )
    parser.add_argument("--max-residuals", type=int, default=100000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_depth < 0 or args.max_residuals <= 0:
        raise SystemExit("adaptive chamber cover: invalid refinement limits")

    g = args.group
    witnesses = parse_atlas(args.atlas, g)
    constants = np.asarray([row["constant_log2"] for row in witnesses])
    charges = np.asarray([row["charge"] for row in witnesses])
    total_profiles = profile_count(g, N)
    target = -40.0 - math.log2(total_profiles) - args.cell_budget_log2
    queue = []
    serial = 0
    initial = 0
    for order in itertools.permutations(range(g + 1)):
        for cell in initial_cells(g, order):
            queue.append((order, 0, serial, cell))
            serial += 1
            initial += 1

    covered = []
    residuals = []
    terminal_failures = []
    while queue:
        order, depth, cell_id, cell = queue.pop()
        matrix = np.asarray(
            [[float(value) for value in profile] for profile in cell]
        )
        normalizations = normalization_log2(g, matrix)
        values = constants[:, None] - charges @ matrix.T - normalizations[None, :]
        finite = np.flatnonzero(np.all(np.isfinite(values), axis=1))
        result = solve_fixed_mixture(values[finite])
        if not result.success:
            raise RuntimeError(f"cell {cell_id} LP failed: {result.message}")
        weights = result.x[: len(finite)]
        mixed = weights @ values[finite]
        maximum = float(np.max(mixed))
        if maximum <= target:
            covered.append(
                {
                    "cell_id": cell_id,
                    "order": list(order),
                    "depth": depth,
                    "vertices": [
                        [str(value) for value in profile] for profile in cell
                    ],
                    "maximum_vertex_log2": maximum,
                    "target_log2": target,
                    "margin_bits": target - maximum,
                    "mixture": [
                        {
                            "witness": witnesses[int(finite[index])]["name"],
                            "weight": float(weight),
                        }
                        for index, weight in enumerate(weights)
                        if weight > 1e-10
                    ],
                }
            )
            continue
        if depth < args.max_depth:
            for child in subdivide(cell):
                queue.append((order, depth + 1, serial, child))
                serial += 1
            continue

        center = barycenter(cell)
        profile = rounded_integer_profile(center)
        failure = {
            "cell_id": cell_id,
            "order": list(order),
            "depth": depth,
            "maximum_vertex_log2": maximum,
            "target_log2": target,
            "profile": profile,
        }
        terminal_failures.append(failure)
        if len(residuals) < args.max_residuals:
            residuals.append(
                {
                    "name": f"cell_{cell_id}_d{depth}_centroid",
                    "profile": profile,
                    "value_log2": maximum,
                    "target_log2": target,
                }
            )

    complete = not terminal_failures
    if complete:
        terms = [
            row["maximum_vertex_log2"] + math.log2(total_profiles)
            for row in covered
        ]
        union = float(logsumexp(np.asarray(terms) * math.log(2.0)) / math.log(2.0))
    else:
        union = math.inf
    report = {
        "status": "DIAGNOSTIC_BINARY64_ADAPTIVE_ORDERED_CHAMBER_COVER",
        "group_bits": g,
        "max_depth": args.max_depth,
        "cell_budget_log2": args.cell_budget_log2,
        "target_log2": target,
        "witnesses": len(witnesses),
        "initial_cells": initial,
        "covered_cells": len(covered),
        "terminal_failed_cells": len(terminal_failures),
        "complete_cover": complete,
        "global_union_log2": union,
        "global_security_margin_bits": -union,
        "global_slack_beyond_40_bits": -40.0 - union,
        "passed_global_40_bits": bool(complete and union <= -40.0),
        "covered_cell_ledger": covered,
        "terminal_failures": terminal_failures,
        "uncovered_integer_residuals": residuals,
        "residuals_truncated": len(residuals) < len(terminal_failures),
        "sources": [str(path) for path in args.atlas],
        "assumptions": [
            "binary64 witness evaluation and LP selection are diagnostic",
            "truncated initial chambers use Qhull discovery triangulations",
            "each closed covered cell is charged the full global profile count",
        ],
    }
    rendered = json.dumps(report, indent=2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(f"initial_cells={initial}")
    print(f"covered_cells={len(covered)}")
    print(f"terminal_failed_cells={len(terminal_failures)}")
    print(f"residuals_emitted={len(residuals)}")
    print(f"complete_cover={complete}")
    if complete:
        print(f"global_union_log2={union:.12f}")
        print(f"global_security_margin_bits={-union:.12f}")


if __name__ == "__main__":
    main()
