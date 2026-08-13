#!/usr/bin/env python3
"""Anchor-seeded Delaunay cover diagnostic for generalized packet profiles.

The mesh uses every optimized witness profile as an interior anchor plus the
exact vertices of the simplex clipped by physical weight 21.  This aligns
cells with the regions where witnesses were actually optimized, avoiding the
large diameter mismatch of a fixed ordered-chamber subdivision.

Qhull and binary64 LPs are discovery tools.  A closing mesh must later be
reconstructed with rational vertices and exact orientation/volume checks.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay

from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import N, atom_count, normalization_log2
from probe_packet_group_support_face_cover import (
    MINIMUM_PHYSICAL_WEIGHT,
    parse_atlas,
    rounded_integer_profile,
    solve_fixed_mixture,
)


def clipped_hull_vertices(group_bits: int) -> list[tuple[Fraction, ...]]:
    atoms = atom_count(group_bits)
    result = []
    for weight in range(1, group_bits + 1):
        pure = [Fraction(0)] * (group_bits + 1)
        pure[weight] = atoms
        result.append(tuple(pure))
        boundary = [Fraction(0)] * (group_bits + 1)
        boundary[weight] = Fraction(MINIMUM_PHYSICAL_WEIGHT, weight)
        boundary[0] = atoms - boundary[weight]
        result.append(tuple(boundary))
    return result


def load_anchor_profiles(paths: list[Path], group_bits: int) -> list[tuple[Fraction, ...]]:
    atoms = atom_count(group_bits)
    rows = set(clipped_hull_vertices(group_bits))
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        candidates = payload.get("rows", payload) if isinstance(payload, dict) else payload
        if not isinstance(candidates, list):
            continue
        for row in candidates:
            profile = row.get("profile")
            if not isinstance(profile, list) or len(profile) != group_bits + 1:
                continue
            point = tuple(Fraction(int(value)) for value in profile)
            if sum(point) != atoms:
                continue
            if sum(index * value for index, value in enumerate(point)) < 21:
                continue
            rows.add(point)
    return sorted(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument(
        "--cell-budget-log2",
        type=float,
        default=24.0,
        help="reserve this many union bits for the final number of mesh cells",
    )
    parser.add_argument("--max-residuals", type=int, default=100000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    g = args.group
    witnesses = parse_atlas(args.atlas, g)
    anchors = load_anchor_profiles(args.atlas, g)
    coordinates = np.asarray(
        [[float(value) for value in profile[1:]] for profile in anchors]
    )
    triangulation = Delaunay(coordinates, qhull_options="Qbb Qc Q12 QJ")
    simplices = []
    seen = set()
    for indices in triangulation.simplices:
        if any(int(index) >= len(anchors) for index in indices):
            continue
        key = tuple(sorted(int(index) for index in indices))
        if key not in seen:
            seen.add(key)
            simplices.append(key)

    constants = np.asarray([row["constant_log2"] for row in witnesses])
    charges = np.asarray([row["charge"] for row in witnesses])
    target = (
        -40.0
        - math.log2(profile_count(g, N))
        - args.cell_budget_log2
    )
    covered = []
    failures = []
    residuals = []
    covered_volume = 0.0
    failed_volume = 0.0
    for cell_id, indices in enumerate(simplices):
        cell = tuple(anchors[index] for index in indices)
        cell_coordinates = np.asarray(
            [[float(value) for value in profile[1:]] for profile in cell]
        )
        # The common factor 1/g! cancels in all volume ratios, so retain the
        # absolute determinant directly for a better-conditioned diagnostic.
        volume = abs(
            float(
                np.linalg.det(
                    cell_coordinates[1:] - cell_coordinates[0]
                )
            )
        )
        matrix = np.asarray(
            [[float(value) for value in profile] for profile in cell]
        )
        normalizations = normalization_log2(g, matrix)
        values = constants[:, None] - charges @ matrix.T - normalizations[None, :]
        finite = np.flatnonzero(np.all(np.isfinite(values), axis=1))
        result = solve_fixed_mixture(values[finite])
        if not result.success:
            raise RuntimeError(f"anchor cell {cell_id} LP failed: {result.message}")
        weights = result.x[: len(finite)]
        mixed = weights @ values[finite]
        maximum = float(np.max(mixed))
        record = {
            "cell_id": cell_id,
            "anchor_indices": list(indices),
            "volume_proxy": volume,
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
        if maximum <= target:
            covered.append(record)
            covered_volume += volume
            continue
        center = tuple(sum(profile[index] for profile in cell) / len(cell) for index in range(g + 1))
        profile = rounded_integer_profile(center)
        record["profile"] = profile
        failures.append(record)
        failed_volume += volume
        if len(residuals) < args.max_residuals:
            residuals.append(
                {
                    "name": f"anchor_cell_{cell_id}_centroid",
                    "profile": profile,
                    "value_log2": maximum,
                    "target_log2": target,
                    "volume_proxy": volume,
                }
            )

    report = {
        "status": "DIAGNOSTIC_BINARY64_ANCHOR_DELAUNAY_PROFILE_COVER",
        "group_bits": g,
        "anchors": len(anchors),
        "witnesses": len(witnesses),
        "simplices": len(simplices),
        "cell_budget_log2": args.cell_budget_log2,
        "target_log2": target,
        "covered_cells": len(covered),
        "failed_cells": len(failures),
        "covered_volume_fraction": covered_volume / (covered_volume + failed_volume),
        "failed_volume_fraction": failed_volume / (covered_volume + failed_volume),
        "complete_cover": not failures,
        "covered_cell_ledger": covered,
        "failed_cell_ledger": failures,
        "uncovered_integer_residuals": residuals,
        "residuals_truncated": len(residuals) < len(failures),
        "anchor_profiles": [[str(value) for value in profile] for profile in anchors],
        "sources": [str(path) for path in args.atlas],
        "assumptions": [
            "Qhull Delaunay geometry and binary64 LP selection are diagnostic",
            "the final mesh requires exact rational reconstruction and volume audit",
        ],
    }
    rendered = json.dumps(report, indent=2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(f"anchors={len(anchors)} witnesses={len(witnesses)}")
    print(f"simplices={len(simplices)}")
    print(f"covered_cells={len(covered)}")
    print(f"failed_cells={len(failures)}")
    print(f"covered_volume_fraction={covered_volume / (covered_volume + failed_volume):.12f}")
    print(f"failed_volume_fraction={failed_volume / (covered_volume + failed_volume):.12f}")
    print(f"complete_cover={not failures}")


if __name__ == "__main__":
    main()
