#!/usr/bin/env python3
"""Diagnostic exact-root, recursive-stellar profile cover for ``g=4``.

The producer starts from a deterministic pulling triangulation of the exact
simplex clipped by physical weight 21.  A leaf is retained as soon as one
fixed witness or one fixed convex mixture is safe at all five vertices.  Only
failed leaves are replaced by their five exact barycentric cones.

Geometry and mixture weights are emitted as exact rational strings in the
recursive schema consumed by ``certify_packet_group_g4_anchor_mesh.py``.
Witness evaluation and LP selection remain binary64 discovery arithmetic;
the verifier must outward-harden a complete ledger before theorem use.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.optimize import linprog

import certify_packet_group_g4_anchor_mesh as verifier
from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import D, K, N, atom_count, normalization_log2
from probe_packet_group_support_face_cover import rounded_integer_profile


GROUP_BITS = 4
CLASSES = GROUP_BITS + 1
M = atom_count(GROUP_BITS)
MINIMUM_PHYSICAL_WEIGHT = 21
OWNER_RULE = "minimum_leaf_owner_rank_among_closed_simplices"
ROOT_MESH_SCHEMA = "packet-group-g4-exact-regular-root-mesh-v1"

Profile = tuple[Fraction, ...]


@dataclass(frozen=True)
class Witness:
    reference: str
    constant_log2: float
    charge: tuple[float, ...]
    support: frozenset[int]
    subtract_normalization: bool

    def value(self, profile: Profile, normalization: float) -> float:
        if any(count and index not in self.support for index, count in enumerate(profile)):
            return math.inf
        result = self.constant_log2 - math.fsum(
            float(count) * coefficient
            for count, coefficient in zip(profile, self.charge)
        )
        return result - normalization if self.subtract_normalization else result


@dataclass(frozen=True)
class Cell:
    identifier: str
    anchor_indices: tuple[int, ...]
    depth: int


@dataclass(frozen=True)
class AtlasAnchor:
    profile: Profile
    worst_gap_bits: float
    references: tuple[str, ...]


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact_rows(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload if isinstance(payload, list) else payload.get(
        "rows", payload.get("witnesses", [])
    )
    if isinstance(rows, dict):
        rows = list(rows.values())
    if not isinstance(rows, list):
        raise ValueError(f"g4 stellar cover: atlas {path} has no witness rows")
    return rows


def integer_anchor_profile(row: dict[str, Any]) -> Profile | None:
    raw = row.get("profile")
    if not isinstance(raw, list) or len(raw) != CLASSES:
        return None
    values = []
    for value in raw:
        if isinstance(value, bool) or isinstance(value, float):
            return None
        try:
            parsed = Fraction(str(value))
        except (ValueError, ZeroDivisionError):
            return None
        if parsed.denominator != 1:
            return None
        values.append(parsed)
    profile = tuple(values)
    if (
        any(value <= 0 for value in profile)
        or sum(profile) != M
        or sum(index * value for index, value in enumerate(profile))
        < MINIMUM_PHYSICAL_WEIGHT
    ):
        return None
    return profile


def row_gap_bits(row: dict[str, Any]) -> float:
    if "gap_bits" in row:
        return float(row["gap_bits"])
    if "margin_bits" in row:
        return -float(row["margin_bits"])
    if "combined_log2" in row and "target_log2" in row:
        return float(row["combined_log2"]) - float(row["target_log2"])
    return -math.inf


def load_witnesses(
    paths: list[Path],
) -> tuple[list[Witness], list[dict[str, str]], list[AtlasAnchor], int]:
    witnesses = []
    sources = []
    anchor_rows: dict[Profile, dict[str, Any]] = {}
    duplicate_anchor_rows = 0
    basenames = [path.name for path in paths]
    if len(set(basenames)) != len(basenames):
        raise ValueError("g4 stellar cover: atlas basenames must be unique")
    for path in paths:
        rows = artifact_rows(path)
        for index, row in enumerate(rows):
            if int(row.get("group_bits", GROUP_BITS)) != GROUP_BITS:
                raise ValueError(f"g4 stellar cover: group mismatch in {path}:{index}")
            raw_charge = row.get("charge")
            if not isinstance(raw_charge, list) or len(raw_charge) != CLASSES:
                continue
            charge = tuple(float(value) for value in raw_charge)
            constant = float(row["constant_log2"])
            if not math.isfinite(constant) or not all(map(math.isfinite, charge)):
                continue
            support = frozenset(int(value) for value in row.get("support", range(CLASSES)))
            if not support or not support.issubset(range(CLASSES)):
                raise ValueError(f"g4 stellar cover: invalid support in {path}:{index}")
            witnesses.append(
                Witness(
                    reference=f"{path.name}:{index}",
                    constant_log2=constant,
                    charge=charge,
                    support=support,
                    subtract_normalization=bool(row.get("subtract_normalization", True)),
                )
            )
            anchor = integer_anchor_profile(row)
            if anchor is not None:
                reference = f"{path.name}:{index}"
                existing = anchor_rows.get(anchor)
                if existing is None:
                    anchor_rows[anchor] = {
                        "gap": row_gap_bits(row),
                        "references": [reference],
                    }
                else:
                    duplicate_anchor_rows += 1
                    existing["gap"] = max(float(existing["gap"]), row_gap_bits(row))
                    existing["references"].append(reference)
        sources.append({"path": str(path.resolve()), "sha256": file_sha256(path)})

    witnesses.append(
        Witness(
            reference="full_bijection",
            constant_log2=K + N * math.log2(11.0) - (N - D) * math.log2(10.0),
            charge=(0.0,) * CLASSES,
            support=frozenset(range(CLASSES)),
            subtract_normalization=True,
        )
    )
    if len(witnesses) == 1:
        raise ValueError("g4 stellar cover: supplied atlases contained no usable witnesses")
    atlas_anchors = [
        AtlasAnchor(
            profile=profile,
            worst_gap_bits=float(metadata["gap"]),
            references=tuple(sorted(metadata["references"])),
        )
        for profile, metadata in sorted(anchor_rows.items())
    ]
    return witnesses, sources, atlas_anchors, duplicate_anchor_rows


def clipped_hull_vertices() -> tuple[Profile, ...]:
    vertices = []
    for weight in range(1, CLASSES):
        boundary = [Fraction(0)] * CLASSES
        boundary[0] = M - Fraction(MINIMUM_PHYSICAL_WEIGHT, weight)
        boundary[weight] = Fraction(MINIMUM_PHYSICAL_WEIGHT, weight)
        vertices.append(tuple(boundary))
        pure = [Fraction(0)] * CLASSES
        pure[weight] = M
        vertices.append(tuple(pure))
    return tuple(sorted(vertices))


def load_root_mesh(path: Path) -> tuple[list[Profile], tuple[Cell, ...], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != ROOT_MESH_SCHEMA:
        raise ValueError(
            f"g4 stellar cover: root mesh schema must be {ROOT_MESH_SCHEMA!r}"
        )
    if int(payload.get("group_bits", -1)) != GROUP_BITS or not bool(
        payload.get("complete_exact_regular_mesh", False)
    ):
        raise ValueError("g4 stellar cover: root mesh is not a complete exact g=4 mesh")
    anchors = list(verifier.load_anchors(payload))
    raw_roots = payload.get("root_cell_ledger")
    if not isinstance(raw_roots, list) or not raw_roots:
        raise ValueError("g4 stellar cover: root mesh lacks root_cell_ledger")
    roots = []
    identifiers = set()
    for row in raw_roots:
        identifier = str(row.get("cell_id", ""))
        if not identifier or identifier in identifiers:
            raise ValueError("g4 stellar cover: root cell identifiers must be unique")
        identifiers.add(identifier)
        indices = verifier.simplex_indices(row, len(anchors))
        roots.append(Cell(identifier, indices, 0))
    expected = int(payload.get("simplices", -1))
    if expected != len(roots):
        raise ValueError("g4 stellar cover: root mesh simplex count mismatch")
    geometry = verifier.verify_mesh(
        {"complete_cover": True, "simplices": len(roots)},
        tuple(anchors),
        [cell_geometry(cell) for cell in roots],
    )
    claimed = payload.get("exact_geometry")
    if not isinstance(claimed, dict):
        raise ValueError("g4 stellar cover: root mesh lacks exact_geometry")
    for key in ("simplices", "determinant_volume", "clipped_hull_determinant_volume"):
        if claimed.get(key) != geometry.get(key):
            raise ValueError(f"g4 stellar cover: root exact_geometry mismatch for {key}")
    audit = {
        "mode": "supplied_exact_regular_root_mesh",
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "schema": payload["schema"],
        "anchors": len(anchors),
        "root_cells": len(roots),
        "source_anchor_profiles_sha256": payload.get("source_anchor_profiles_sha256"),
        "lifting": payload.get("lifting"),
        "exact_geometry": claimed,
    }
    return anchors, tuple(roots), audit


def matrix_rank(rows: list[list[Fraction]]) -> int:
    if not rows:
        return 0
    matrix = [list(row) for row in rows]
    columns = len(matrix[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [value / scale for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def affine_dimension(indices: Iterable[int], anchors: tuple[Profile, ...]) -> int:
    selected = tuple(indices)
    if not selected:
        return -1
    origin = anchors[selected[0]][1:]
    rows = [
        [anchors[index][coordinate] - origin[coordinate - 1] for coordinate in range(1, 5)]
        for index in selected[1:]
    ]
    return matrix_rank(rows)


def hull_facets(anchors: tuple[Profile, ...]) -> tuple[frozenset[int], ...]:
    result = []
    for coordinate in range(CLASSES):
        result.append(
            frozenset(index for index, row in enumerate(anchors) if row[coordinate] == 0)
        )
    result.append(
        frozenset(
            index
            for index, row in enumerate(anchors)
            if sum(weight * row[weight] for weight in range(CLASSES))
            == MINIMUM_PHYSICAL_WEIGHT
        )
    )
    return tuple(result)


def pulling_triangulation(anchors: tuple[Profile, ...]) -> tuple[tuple[int, ...], ...]:
    """Triangulate the clipped hull using exact face incidences only."""

    facets = hull_facets(anchors)

    def triangulate(face: frozenset[int], dimension: int) -> list[tuple[int, ...]]:
        if dimension == 0:
            if len(face) != 1:
                raise RuntimeError("zero-dimensional pulling face is not one vertex")
            return [(next(iter(face)),)]
        pivot = min(face)
        opposite = set()
        for global_facet in facets:
            candidate = face & global_facet
            if pivot in candidate or not candidate:
                continue
            if affine_dimension(candidate, anchors) == dimension - 1:
                opposite.add(frozenset(candidate))
        if not opposite:
            raise RuntimeError("pulling triangulation found no opposite facets")
        result = []
        for facet in sorted(opposite, key=lambda row: tuple(sorted(row))):
            for simplex in triangulate(facet, dimension - 1):
                result.append((pivot, *simplex))
        return result

    roots = triangulate(frozenset(range(len(anchors))), GROUP_BITS)
    keys = [tuple(sorted(row)) for row in roots]
    if len(set(keys)) != len(keys):
        raise RuntimeError("pulling triangulation produced duplicate roots")
    return tuple(roots)


def profile_simplex_determinant(profiles: list[Profile]) -> Fraction:
    if len(profiles) != CLASSES:
        raise ValueError("g4 simplex determinant needs five profiles")
    origin = profiles[0][1:]
    return verifier.determinant(
        [
            [profiles[row][coordinate] - origin[coordinate - 1] for coordinate in range(1, 5)]
            for row in range(1, 5)
        ]
    )


def barycentric_coordinates(
    point: Profile, cell: Cell, anchors: list[Profile]
) -> tuple[Fraction, ...]:
    profiles = [anchors[index] for index in cell.anchor_indices]
    denominator = profile_simplex_determinant(profiles)
    if not denominator:
        raise RuntimeError(f"g4 stellar cover: degenerate cell {cell.identifier}")
    result = []
    for replaced in range(CLASSES):
        replacement = list(profiles)
        replacement[replaced] = point
        result.append(profile_simplex_determinant(replacement) / denominator)
    if sum(result, Fraction(0)) != 1:
        raise RuntimeError("g4 stellar cover: barycentric coordinates do not sum to one")
    return tuple(result)


def locate_anchor(
    point: Profile, leaves: dict[str, Cell], anchors: list[Profile]
) -> tuple[str, Cell | None, tuple[str, ...]]:
    strict = []
    boundary = []
    for identifier in sorted(leaves):
        cell = leaves[identifier]
        coordinates = barycentric_coordinates(point, cell, anchors)
        if all(value > 0 for value in coordinates):
            strict.append(cell)
        elif all(value >= 0 for value in coordinates):
            boundary.append(identifier)
    if len(strict) > 1:
        raise RuntimeError("g4 stellar cover: point is strictly inside multiple leaves")
    if strict:
        return "strict", strict[0], ()
    if boundary:
        return "boundary", None, tuple(boundary)
    return "outside", None, ()


def squared_distance(left: Profile, right: Profile) -> Fraction:
    return sum(
        ((left[index] - right[index]) / M) ** 2 for index in range(CLASSES)
    )


def order_atlas_anchors(
    rows: list[AtlasAnchor], mode: str, root_anchors: tuple[Profile, ...]
) -> list[AtlasAnchor]:
    if mode == "lexicographic":
        return sorted(rows, key=lambda row: row.profile)
    if mode == "worst-first":
        return sorted(rows, key=lambda row: (-row.worst_gap_bits, row.profile))
    if mode != "farthest-first":
        raise ValueError(f"unsupported atlas insertion order {mode!r}")
    remaining = list(sorted(rows, key=lambda row: row.profile))
    selected: list[AtlasAnchor] = []
    distances = {
        row.profile: min(squared_distance(row.profile, root) for root in root_anchors)
        for row in remaining
    }
    while remaining:
        chosen = max(
            remaining,
            key=lambda row: (
                distances[row.profile],
                row.worst_gap_bits,
                tuple(-value for value in row.profile),
            ),
        )
        selected.append(chosen)
        remaining.remove(chosen)
        for row in remaining:
            distances[row.profile] = min(
                distances[row.profile], squared_distance(row.profile, chosen.profile)
            )
    return selected


def compact_digest(value: Any) -> str:
    encoded = json.dumps(value, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def integer_profile_count_upper(
    cell: Cell, anchors: list[Profile], global_profile_count: int
) -> tuple[int, int, tuple[int, ...]]:
    """Bound lattice profiles by the smallest dropped-coordinate box.

    The simplex projection onto any four of its five profile coordinates is
    contained in the product of their inclusive integer ranges.  Dropping
    one coordinate is injective because all profiles have fixed total mass.
    """

    ranges = []
    for coordinate in range(CLASSES):
        values = [anchors[index][coordinate] for index in cell.anchor_indices]
        width = max(
            0,
            floor_fraction(max(values)) - ceil_fraction(min(values)) + 1,
        )
        ranges.append(width)
    candidates = []
    for dropped in range(CLASSES):
        product = math.prod(
            ranges[coordinate]
            for coordinate in range(CLASSES)
            if coordinate != dropped
        )
        candidates.append((product, dropped))
    raw, dropped = min(candidates)
    return min(raw, global_profile_count), dropped, tuple(ranges)


def log2sumexp_upperish(terms: list[float]) -> float:
    """Stable binary64 log-sum with explicit upward nudges (diagnostic only)."""

    if not terms:
        return -math.inf
    maximum = max(terms)
    scaled_sum = 0.0
    for value in terms:
        scaled = math.nextafter(2.0 ** (value - maximum), math.inf)
        scaled_sum = math.nextafter(scaled_sum + scaled, math.inf)
    return math.nextafter(maximum + math.log2(scaled_sum), math.inf)


class Evaluator:
    def __init__(self, witnesses: list[Witness]):
        self.witnesses = witnesses
        self.cache: dict[Profile, tuple[float, tuple[float, ...]]] = {}

    def evaluate(self, profile: Profile) -> tuple[float, tuple[float, ...]]:
        cached = self.cache.get(profile)
        if cached is not None:
            return cached
        matrix = np.asarray([[float(value) for value in profile]], dtype=np.float64)
        normalization = float(normalization_log2(GROUP_BITS, matrix)[0])
        values = tuple(row.value(profile, normalization) for row in self.witnesses)
        result = normalization, values
        self.cache[profile] = result
        return result


def mixture_candidates(values: np.ndarray, per_vertex: int) -> np.ndarray:
    finite = np.flatnonzero(np.all(np.isfinite(values), axis=0))
    if not len(finite):
        return finite
    if per_vertex > 0 and len(finite) > per_vertex:
        selected = set()
        for vertex in range(values.shape[0]):
            # Witness index is the explicit tie-break, so a capped diagnostic
            # does not depend on NumPy's sort stability for equal values.
            order = finite[np.lexsort((finite, values[vertex, finite]))]
            selected.update(map(int, order[:per_vertex]))
        finite = np.asarray(sorted(selected), dtype=np.int64)
    rows = values[:, finite]
    keep = np.ones(len(finite), dtype=bool)
    for local in range(len(finite)):
        dominated = np.all(rows <= rows[:, local, None], axis=0) & np.any(
            rows < rows[:, local, None], axis=0
        )
        dominated[local] = False
        if np.any(dominated):
            keep[local] = False
    return finite[keep]


def optimize_values(
    values: np.ndarray,
    witnesses: list[Witness],
    target: float,
    safety_bits: float,
    candidates_per_vertex: int,
    weight_tolerance: float,
    candidate_mode: str = "capped",
    reduced_cost_tolerance: float = 1e-9,
) -> dict[str, Any]:
    if candidate_mode not in {"capped", "all", "column-generation"}:
        raise ValueError(f"g4 stellar cover: unknown candidate mode {candidate_mode}")
    scores = np.max(values, axis=0)
    best = int(np.argmin(scores))
    if (
        candidate_mode == "capped"
        and math.isfinite(float(scores[best]))
        and scores[best] <= target - safety_bits
    ):
        return {
            "indices": (best,),
            "weights": (Fraction(1),),
            "score": float(scores[best]),
            "vertex_values": tuple(map(float, values[:, best])),
            "passes_uniform_target": True,
            "candidate_count": 1,
            "selection_mode": "single_witness_uniform_target_fast_path",
            "column_generation_rounds": 0,
            "full_column_scans": 0,
            "minimum_reduced_cost": None,
            "maximum_pricing_gap": None,
            "reduced_cost_tolerance": reduced_cost_tolerance,
        }

    eligible = np.flatnonzero(np.all(np.isfinite(values), axis=0))
    if not len(eligible):
        raise RuntimeError("g4 stellar cover: no finite witness on a cell")
    if candidate_mode == "all" or candidates_per_vertex == 0:
        candidates = mixture_candidates(values, 0)
    else:
        seed_count = max(1, candidates_per_vertex)
        candidates = mixture_candidates(values, seed_count)

    def solve(selected: np.ndarray):
        shifted_selected = values[:, selected] - target
        count = len(selected)
        solved = linprog(
            np.append(np.zeros(count), 1.0),
            A_ub=np.hstack((shifted_selected, -np.ones((values.shape[0], 1)))),
            b_ub=np.zeros(values.shape[0]),
            A_eq=np.asarray([[1.0] * count + [0.0]]),
            b_eq=np.asarray([1.0]),
            bounds=[(0.0, None)] * count + [(None, None)],
            method="highs-ds",
        )
        if solved.status == 2:
            raise RuntimeError("g4 stellar cover: finite mixture LP is infeasible")
        if not solved.success:
            raise RuntimeError(f"g4 stellar mixture LP failed: {solved.message}")
        return solved

    column_rounds = 0
    full_scans = 0
    minimum_reduced_cost: float | None = None
    maximum_pricing_gap: float | None = None
    result = solve(candidates)
    if candidate_mode == "column-generation":
        shifted_all = values[:, eligible] - target
        while True:
            # HiGHS reports nonpositive marginals for our <= rows.  The raw
            # game-dual weights are their negatives.  Validate and price with
            # these raw values; do not truncate/renormalize before the KKT test.
            dual = -np.asarray(result.ineqlin.marginals, dtype=np.float64)
            if (
                len(dual) != values.shape[0]
                or np.min(dual) < -1e-10
                or abs(float(np.sum(dual)) - 1.0) > 1e-7
            ):
                raise RuntimeError("g4 stellar mixture LP returned an invalid vertex dual")
            restricted_prices = dual @ (values[:, candidates] - target)
            restricted_gap = float(result.fun) - float(np.min(restricted_prices))
            if abs(restricted_gap) > 1e-7 * max(1.0, abs(float(result.fun))):
                raise RuntimeError(
                    "g4 stellar mixture LP failed restricted primal-dual KKT"
                )
            prices = dual @ shifted_all
            reduced = prices - float(result.fun)
            full_scans += 1
            minimum_reduced_cost = float(np.min(reduced))
            maximum_pricing_gap = float(result.fun) - float(np.min(prices))
            selected_set = set(map(int, candidates))
            violating = [
                (float(reduced[local]), int(index))
                for local, index in enumerate(eligible)
                if int(index) not in selected_set
                and float(reduced[local]) < -reduced_cost_tolerance
            ]
            if not violating:
                break
            violating.sort(key=lambda item: (item[0], item[1]))
            additions = [index for _cost, index in violating[: values.shape[0]]]
            candidates = np.asarray(
                sorted(selected_set.union(additions)), dtype=np.int64
            )
            result = solve(candidates)
            column_rounds += 1

    count = len(candidates)
    active = np.flatnonzero(result.x[:count] > weight_tolerance)
    if not len(active):
        raise RuntimeError("g4 stellar cover: mixture LP returned empty support")
    indices = tuple(int(candidates[index]) for index in active)
    raw = tuple(Fraction.from_float(float(result.x[index])) for index in active)
    total = sum(raw, Fraction(0))
    if total <= 0:
        raise RuntimeError("g4 stellar cover: mixture LP returned zero weight")
    weights = tuple(value / total for value in raw)
    mixed = tuple(
        math.fsum(float(weight) * float(values[vertex, index]) for index, weight in zip(indices, weights))
        for vertex in range(values.shape[0])
    )
    score = max(mixed)
    return {
        "indices": indices,
        "weights": weights,
        "score": score,
        "vertex_values": mixed,
        "passes_uniform_target": score <= target - safety_bits,
        "candidate_count": int(len(candidates)),
        "selection_mode": f"minimax_mixture_lp_{candidate_mode}",
        "column_generation_rounds": column_rounds,
        "full_column_scans": full_scans,
        "minimum_reduced_cost": minimum_reduced_cost,
        "maximum_pricing_gap": maximum_pricing_gap,
        "reduced_cost_tolerance": reduced_cost_tolerance,
        "sparse_basic_solution": len(indices) <= values.shape[0],
    }


def evaluate_cell(
    cell: Cell,
    anchors: list[Profile],
    evaluator: Evaluator,
    target: float,
    safety_bits: float,
    candidates_per_vertex: int,
    weight_tolerance: float,
) -> dict[str, Any]:
    values = np.asarray(
        [evaluator.evaluate(anchors[index])[1] for index in cell.anchor_indices],
        dtype=np.float64,
    )
    result = optimize_values(
        values,
        evaluator.witnesses,
        target,
        safety_bits,
        candidates_per_vertex,
        weight_tolerance,
    )
    result["best_single_witness_max_vertex_log2"] = float(
        np.min(np.max(values, axis=0))
    )
    return result


def mixture_record(certificate: dict[str, Any], witnesses: list[Witness]) -> list[dict[str, str]]:
    return [
        {
            "witness": witnesses[index].reference,
            "weight_exact": fraction_text(weight),
        }
        for index, weight in zip(certificate["indices"], certificate["weights"])
    ]


def cell_bound_record(
    cell: Cell,
    anchors: list[Profile],
    result: dict[str, Any],
    witnesses: list[Witness],
    global_profile_count: int,
) -> dict[str, Any]:
    count, dropped, ranges = integer_profile_count_upper(
        cell, anchors, global_profile_count
    )
    count_log2 = math.nextafter(math.log2(count), math.inf) if count else -math.inf
    contribution = (
        math.nextafter(result["score"] + count_log2, math.inf)
        if count
        else -math.inf
    )
    return {
        "mixture": mixture_record(result, witnesses),
        "maximum_vertex_log2": result["score"],
        "vertex_mixture_values_log2": list(result["vertex_values"]),
        "integer_profile_count_upper": count,
        "integer_profile_count_method": (
            "minimum_dropped_coordinate_inclusive_integer_bounding_box_"
            "capped_by_global_profile_count"
        ),
        "integer_profile_count_dropped_coordinate": dropped,
        "inclusive_integer_coordinate_range_widths": list(ranges),
        "cell_local_contribution_log2_upperish": contribution,
    }


def cell_geometry(cell: Cell) -> dict[str, Any]:
    return {
        "cell_id": cell.identifier,
        "anchor_indices": list(cell.anchor_indices),
        "depth": cell.depth,
    }


def centroid_profile(cell: Cell, anchors: list[Profile]) -> Profile:
    return tuple(
        sum((anchors[index][coordinate] for index in cell.anchor_indices), Fraction(0))
        / CLASSES
        for coordinate in range(CLASSES)
    )


def add_centroid(cell: Cell, anchors: list[Profile], lookup: dict[Profile, int]) -> int:
    centroid = centroid_profile(cell, anchors)
    existing = lookup.get(centroid)
    if existing is not None:
        return existing
    index = len(anchors)
    anchors.append(centroid)
    lookup[centroid] = index
    return index


def subdivide(cell: Cell, insertion_index: int) -> tuple[Cell, ...]:
    children = []
    for omitted in range(CLASSES):
        indices = (
            insertion_index,
            *cell.anchor_indices[:omitted],
            *cell.anchor_indices[omitted + 1 :],
        )
        children.append(Cell(f"{cell.identifier}.{omitted}", tuple(indices), cell.depth + 1))
    return tuple(children)


def insert_atlas_anchors(
    candidates: list[AtlasAnchor],
    order: str,
    limit: int | None,
    roots: tuple[Cell, ...],
    anchors: list[Profile],
    anchor_lookup: dict[Profile, int],
    subdivisions: list[dict[str, Any]],
    duplicate_source_rows: int,
) -> tuple[list[Cell], dict[str, Any]]:
    ordered = order_atlas_anchors(candidates, order, tuple(anchors))
    if limit is not None:
        ordered = ordered[:limit]
    leaves = {cell.identifier: cell for cell in roots}
    inserted_rows = []
    skipped_rows = []
    for candidate in ordered:
        existing = anchor_lookup.get(candidate.profile)
        if existing is not None:
            skipped_rows.append(
                {
                    "profile": [int(value) for value in candidate.profile],
                    "reason": "duplicate_existing_anchor",
                    "existing_anchor_index": existing,
                    "witness_references": list(candidate.references),
                }
            )
            continue
        disposition, parent, boundary_cells = locate_anchor(
            candidate.profile, leaves, anchors
        )
        if disposition != "strict" or parent is None:
            skipped_rows.append(
                {
                    "profile": [int(value) for value in candidate.profile],
                    "reason": (
                        "current_mesh_boundary" if disposition == "boundary" else "outside_clipped_hull"
                    ),
                    "boundary_cell_ids": list(boundary_cells),
                    "witness_references": list(candidate.references),
                }
            )
            continue

        insertion_index = len(anchors)
        anchors.append(candidate.profile)
        anchor_lookup[candidate.profile] = insertion_index
        children = subdivide(parent, insertion_index)
        del leaves[parent.identifier]
        leaves.update((child.identifier, child) for child in children)
        barycentric = barycentric_coordinates(candidate.profile, parent, anchors)
        subdivisions.append(
            {
                "parent_cell_id": parent.identifier,
                "insertion_anchor_index": insertion_index,
                "insertion_kind": "atlas_anchor",
                "insertion_barycentric_coordinates": [
                    fraction_text(value) for value in barycentric
                ],
                "children": [cell_geometry(child) for child in children],
            }
        )
        inserted_rows.append(
            {
                "profile": [int(value) for value in candidate.profile],
                "anchor_index": insertion_index,
                "parent_cell_id": parent.identifier,
                "worst_gap_bits": candidate.worst_gap_bits,
                "witness_references": list(candidate.references),
            }
        )

    inserted_payload = [row["profile"] for row in inserted_rows]
    skipped_payload = [
        {"profile": row["profile"], "reason": row["reason"]}
        for row in skipped_rows
    ]
    audit = {
        "enabled": True,
        "order": order,
        "unique_full_support_integer_candidates": len(candidates),
        "considered_candidates": len(ordered),
        "duplicate_source_rows": duplicate_source_rows,
        "inserted_anchors": len(inserted_rows),
        "skipped_anchors": len(skipped_rows),
        "inserted_profiles_sha256": compact_digest(inserted_payload),
        "skipped_profiles_and_reasons_sha256": compact_digest(skipped_payload),
        "inserted": inserted_rows,
        "skipped": skipped_rows,
    }
    return [leaves[key] for key in sorted(leaves)], audit


def render_report(
    *,
    anchors: list[Profile],
    roots: tuple[Cell, ...],
    subdivisions: list[dict[str, Any]],
    covered: dict[str, dict[str, Any]],
    failed: list[tuple[Cell, dict[str, Any]]],
    witnesses: list[Witness],
    sources: list[dict[str, str]],
    target: float,
    safety_bits: float,
    wave_rows: list[dict[str, Any]],
    atlas_insertion_audit: dict[str, Any],
    root_mesh_audit: dict[str, Any],
) -> dict[str, Any]:
    covered_rows = []
    for owner_rank, identifier in enumerate(sorted(covered)):
        row = dict(covered[identifier])
        row["owner_rank"] = owner_rank
        covered_rows.append(row)
    failed_rows = []
    residual_by_profile: dict[tuple[int, ...], dict[str, Any]] = {}
    for cell, result in sorted(
        failed, key=lambda item: item[1]["score"], reverse=True
    ):
        center = centroid_profile(cell, anchors)
        profile = rounded_integer_profile(center)
        best_score = result["best_single_witness_max_vertex_log2"]
        gap = result["score"] - target
        row = {
            **cell_geometry(cell),
            **cell_bound_record(
                cell,
                anchors,
                result,
                witnesses,
                profile_count(GROUP_BITS, N),
            ),
            "centroid_profile_exact": [fraction_text(value) for value in center],
            "residual_centroid_profile": profile,
            "best_single_witness_max_vertex_log2": best_score,
            "target_log2": target,
            "gap_bits": gap,
        }
        failed_rows.append(row)
        key = tuple(profile)
        current = residual_by_profile.get(key)
        if current is None or gap > float(current["gap_bits"]):
            residual_by_profile[key] = {
                "name": f"stellar_{cell.identifier}_centroid",
                "profile": profile,
                "value_log2": result["score"],
                "target_log2": target,
                "gap_bits": gap,
                "source_cell_id": cell.identifier,
                "depth": cell.depth,
            }
    residuals = sorted(
        residual_by_profile.values(), key=lambda row: float(row["gap_bits"]), reverse=True
    )
    complete = not failed_rows
    aggregation_rows = [*covered_rows, *failed_rows]
    contribution_rows = sorted(
        (
            {
                "cell_id": row["cell_id"],
                "maximum_vertex_log2": row["maximum_vertex_log2"],
                "integer_profile_count_upper": row["integer_profile_count_upper"],
                "contribution_log2_upperish": row[
                    "cell_local_contribution_log2_upperish"
                ],
                "uniform_target_passed": (
                    row["maximum_vertex_log2"] <= target - safety_bits
                ),
            }
            for row in aggregation_rows
        ),
        key=lambda row: row["contribution_log2_upperish"],
        reverse=True,
    )
    aggregate_union = log2sumexp_upperish(
        [row["contribution_log2_upperish"] for row in contribution_rows]
    )
    return {
        "schema": "packet-group-g4-recursive-stellar-cover-v1",
        "status": "DIAGNOSTIC_BINARY64_G4_RECURSIVE_EXACT_STELLAR_COVER",
        "group_bits": GROUP_BITS,
        "target_log2": target,
        "safety_bits": safety_bits,
        "profile_owner_rule": OWNER_RULE,
        "union_accounting_mode": "cell_local_profile_count_upper",
        "complete_cover": complete,
        "simplices": len(covered_rows) + len(failed_rows),
        "anchors": len(anchors),
        "witnesses": len(witnesses),
        "sources": sources,
        "anchor_profiles": [
            [fraction_text(value) for value in profile] for profile in anchors
        ],
        "root_cell_ledger": [cell_geometry(cell) for cell in roots],
        "stellar_subdivision_ledger": subdivisions,
        "covered_cell_ledger": covered_rows,
        "failed_cell_ledger": failed_rows,
        "uncovered_integer_residuals": residuals,
        "wave_checkpoint_ledger": wave_rows,
        "atlas_insertion_audit": atlas_insertion_audit,
        "root_mesh_audit": root_mesh_audit,
        "cell_local_aggregation_diagnostic": {
            "status": "DIAGNOSTIC_BINARY64_UPWARD_NUDGED_NOT_OUTWARD_CERTIFIED",
            "overlap_policy": "closed_cells_may_overlap_and_are_deliberately_overcounted",
            "terminal_leaves_aggregated": len(aggregation_rows),
            "global_profile_count": profile_count(GROUP_BITS, N),
            "union_log2_upperish": aggregate_union,
            "security_margin_bits_upperish": -40.0 - aggregate_union,
            "passed_40_bits_diagnostic": aggregate_union <= -40.0,
            "worst_contributors": contribution_rows[:20],
        },
        "summary": {
            "root_cells": len(roots),
            "stellar_subdivisions": len(subdivisions),
            "covered_leaves": len(covered_rows),
            "failed_leaves": len(failed_rows),
            "unique_residual_centroid_profiles": len(residuals),
        },
        "assumptions": [
            "The root pulling triangulation uses exact rational clipped-hull vertices.",
            "Every subdivision is the five exact barycentric cones of one failed parent.",
            "One unchanged fixed witness or fixed convex mixture passes all five leaf vertices.",
            "Binary64 witness evaluation and LP selection are diagnostic until outward verification.",
            "Incomplete ledgers retain failed leaves and do not assert a continuous cover.",
        ],
    }


def write_json(path: Path, report: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def run_self_test() -> None:
    root_anchors = clipped_hull_vertices()
    root_indices = pulling_triangulation(root_anchors)
    roots = tuple(Cell(f"r{index:04d}", row, 0) for index, row in enumerate(root_indices))
    root_rows = [cell_geometry(cell) for cell in roots]
    geometry = verifier.verify_mesh(
        {"complete_cover": True, "simplices": len(roots)}, root_anchors, root_rows
    )
    if geometry["determinant_volume"] != geometry["clipped_hull_determinant_volume"]:
        raise SystemExit("g4 stellar self-test: root volume mismatch")

    anchors = list(root_anchors)
    lookup = {profile: index for index, profile in enumerate(anchors)}
    # Use a deliberately non-centroid strict insertion to exercise the
    # verifier's generalized stellar schema.
    weights = (Fraction(1, 15), Fraction(2, 15), Fraction(3, 15), Fraction(4, 15), Fraction(5, 15))
    insertion_profile = tuple(
        sum(
            (weight * root_anchors[index][coordinate] for index, weight in zip(roots[0].anchor_indices, weights)),
            Fraction(0),
        )
        for coordinate in range(CLASSES)
    )
    location, parent, boundary = locate_anchor(
        insertion_profile, {cell.identifier: cell for cell in roots}, anchors
    )
    if location != "strict" or parent != roots[0] or boundary:
        raise SystemExit("g4 stellar self-test: exact strict locator mismatch")
    insertion = len(anchors)
    anchors.append(insertion_profile)
    lookup[insertion_profile] = insertion
    children = subdivide(roots[0], insertion)
    subdivision = {
        "parent_cell_id": roots[0].identifier,
        "insertion_anchor_index": insertion,
        "insertion_kind": "self_test_noncentroid",
        "children": [cell_geometry(cell) for cell in children],
    }
    leaves = [*roots[1:], *children]
    ledger = {
        "anchor_profiles": [
            [fraction_text(value) for value in profile] for profile in anchors
        ],
        "root_cell_ledger": root_rows,
        "stellar_subdivision_ledger": [subdivision],
        "covered_cell_ledger": [
            {
                **cell_geometry(cell),
                "owner_rank": rank,
                "mixture": [
                    {"witness": "full_bijection", "weight_exact": "1/1"}
                ],
            }
            for rank, cell in enumerate(leaves)
        ],
        "failed_cell_ledger": [],
    }
    parsed_anchors = verifier.load_anchors(ledger)
    stellar = verifier.stellar_leaf_rows(ledger, parsed_anchors)
    if stellar is None or len(stellar[0]) != len(leaves):
        raise SystemExit("g4 stellar self-test: verifier schema mismatch")
    for row in ledger["covered_cell_ledger"]:
        if verifier.parse_mixture(row) != (("full_bijection", Fraction(1)),):
            raise SystemExit("g4 stellar self-test: exact mixture schema mismatch")

    # Insert two non-centroid synthetic atlas anchors sequentially and ensure
    # both the producer locator and verifier accept the resulting leaf tree.
    atlas_anchors = list(root_anchors)
    atlas_lookup = {profile: index for index, profile in enumerate(atlas_anchors)}
    atlas_subdivisions: list[dict[str, Any]] = []
    first = AtlasAnchor(insertion_profile, 2.0, ("synthetic:0",))
    # Strictly inside child .0: positive mass on the first insertion and on
    # each of the four retained parent vertices.
    first_child_profile = tuple(
        (
            2 * insertion_profile[coordinate]
            + sum(
                root_anchors[index][coordinate]
                for index in roots[0].anchor_indices[1:]
            )
        )
        / 6
        for coordinate in range(CLASSES)
    )
    second = AtlasAnchor(first_child_profile, 1.0, ("synthetic:1",))
    atlas_leaves, audit = insert_atlas_anchors(
        [first, second],
        "worst-first",
        None,
        roots,
        atlas_anchors,
        atlas_lookup,
        atlas_subdivisions,
        0,
    )
    if audit["inserted_anchors"] != 2 or audit["skipped_anchors"]:
        raise SystemExit("g4 stellar self-test: atlas insertion audit mismatch")
    atlas_ledger = {
        "anchor_profiles": [
            [fraction_text(value) for value in profile] for profile in atlas_anchors
        ],
        "root_cell_ledger": root_rows,
        "stellar_subdivision_ledger": atlas_subdivisions,
        "covered_cell_ledger": [
            {
                **cell_geometry(cell),
                "owner_rank": rank,
                "mixture": [{"witness": "full_bijection", "weight_exact": "1/1"}],
            }
            for rank, cell in enumerate(atlas_leaves)
        ],
        "failed_cell_ledger": [],
    }
    atlas_stellar = verifier.stellar_leaf_rows(
        atlas_ledger, verifier.load_anchors(atlas_ledger)
    )
    if atlas_stellar is None or len(atlas_stellar[0]) != len(atlas_leaves):
        raise SystemExit("g4 stellar self-test: atlas-aligned verifier mismatch")

    # Exercise the exact-regular-root input contract without invoking Qhull.
    import tempfile

    root_mesh_payload = {
        "schema": ROOT_MESH_SCHEMA,
        "group_bits": GROUP_BITS,
        "complete_exact_regular_mesh": True,
        "anchors": len(root_anchors),
        "simplices": len(roots),
        "anchor_profiles": [
            [fraction_text(value) for value in profile] for profile in root_anchors
        ],
        "root_cell_ledger": root_rows,
        "lifting": {"self_test": True},
        "exact_geometry": {"mode": "self_test", **geometry},
        "source_anchor_profiles_sha256": compact_digest([]),
    }
    with tempfile.TemporaryDirectory() as directory:
        root_path = Path(directory) / "root.json"
        root_path.write_text(json.dumps(root_mesh_payload), encoding="utf-8")
        loaded_anchors, loaded_roots, root_audit = load_root_mesh(root_path)
    if (
        loaded_anchors != list(root_anchors)
        or loaded_roots != roots
        or root_audit["root_cells"] != len(roots)
    ):
        raise SystemExit("g4 stellar self-test: exact root-mesh contract mismatch")

    synthetic = np.full((5, 5), 1.0)
    np.fill_diagonal(synthetic, -10.0)
    fake = [
        Witness(str(index), 0.0, (0.0,) * 5, frozenset(range(5)), True)
        for index in range(5)
    ]
    certificate = optimize_values(synthetic, fake, 0.0, 0.1, 5, 1e-12)
    if certificate is None or sum(certificate["weights"], Fraction(0)) != 1:
        raise SystemExit("g4 stellar self-test: exact mixture did not close")
    above_target = optimize_values(synthetic, fake, -100.0, 0.1, 5, 1e-12)
    if above_target["passes_uniform_target"] or not above_target["weights"]:
        raise SystemExit("g4 stellar self-test: above-target mixture was not retained")
    global_count = profile_count(GROUP_BITS, N)
    count_upper, dropped, ranges = integer_profile_count_upper(
        roots[0], list(root_anchors), global_count
    )
    if not 0 < count_upper <= global_count or not 0 <= dropped < CLASSES or len(ranges) != CLASSES:
        raise SystemExit("g4 stellar self-test: invalid cell lattice-count upper")
    expected_logsum = 4.0 + math.log2(3.0)
    actual_logsum = log2sumexp_upperish([4.0, 4.0, 4.0])
    if not expected_logsum <= actual_logsum <= math.nextafter(expected_logsum, math.inf) * (1 + 1e-15):
        raise SystemExit("g4 stellar self-test: diagnostic logsumexp mismatch")
    print(f"self_test_root_cells={len(roots)}")
    print(f"self_test_stellar_leaves={len(leaves)}")
    print(f"self_test_root_volume={geometry['determinant_volume']}")
    print("status=PASS_PACKET_GROUP_G4_STELLAR_COVER_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, action="append", default=[])
    parser.add_argument(
        "--root-mesh",
        type=Path,
        help="exact regular root mesh from build_packet_group_g4_exact_regular_mesh.py",
    )
    parser.add_argument(
        "--atlas-insertion-order",
        choices=("farthest-first", "worst-first", "lexicographic"),
        default="farthest-first",
    )
    parser.add_argument(
        "--max-atlas-insertions",
        type=int,
        help="consider at most this many ordered unique full-support atlas anchors",
    )
    parser.add_argument("--disable-atlas-insertion", action="store_true")
    parser.add_argument("--max-waves", type=int, default=0)
    parser.add_argument("--max-leaves", type=int, default=1000000)
    parser.add_argument("--safety-bits", type=float, default=1e-5)
    parser.add_argument("--mixture-candidates-per-vertex", type=int, default=32)
    parser.add_argument("--mixture-weight-tolerance", type=float, default=1e-10)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if not args.atlas or args.output is None:
        raise SystemExit("g4 stellar cover: --atlas and --output are required")
    if (
        args.max_waves < 0
        or args.max_leaves <= 0
        or args.safety_bits < 0.0
        or args.mixture_candidates_per_vertex < 0
        or args.mixture_weight_tolerance < 0.0
        or (args.max_atlas_insertions is not None and args.max_atlas_insertions < 0)
    ):
        raise SystemExit("g4 stellar cover: invalid search parameter")

    witnesses, sources, atlas_candidates, duplicate_anchor_rows = load_witnesses(
        args.atlas
    )
    evaluator = Evaluator(witnesses)
    if args.root_mesh is None:
        anchors = list(clipped_hull_vertices())
        root_indices = pulling_triangulation(tuple(anchors))
        roots = tuple(
            Cell(f"r{index:04d}", indices, 0)
            for index, indices in enumerate(root_indices)
        )
        geometry = verifier.verify_mesh(
            {"complete_cover": True, "simplices": len(roots)},
            tuple(anchors),
            [cell_geometry(cell) for cell in roots],
        )
        root_mesh_audit = {
            "mode": "built_in_four_root_exact_pulling_mesh",
            "anchors": len(anchors),
            "root_cells": len(roots),
            "exact_geometry": geometry,
        }
    else:
        anchors, roots, root_mesh_audit = load_root_mesh(args.root_mesh)
    anchor_lookup = {profile: index for index, profile in enumerate(anchors)}
    subdivisions: list[dict[str, Any]] = []
    if args.root_mesh is not None or args.disable_atlas_insertion:
        active = list(roots)
        atlas_insertion_audit = {
            "enabled": False,
            "disabled_reason": (
                "supplied_root_mesh_already_contains_selected_anchors"
                if args.root_mesh is not None
                else "command_line"
            ),
            "order": args.atlas_insertion_order,
            "unique_full_support_integer_candidates": len(atlas_candidates),
            "considered_candidates": 0,
            "duplicate_source_rows": duplicate_anchor_rows,
            "inserted_anchors": 0,
            "skipped_anchors": 0,
            "inserted_profiles_sha256": compact_digest([]),
            "skipped_profiles_and_reasons_sha256": compact_digest([]),
            "inserted": [],
            "skipped": [],
        }
    else:
        active, atlas_insertion_audit = insert_atlas_anchors(
            atlas_candidates,
            args.atlas_insertion_order,
            args.max_atlas_insertions,
            roots,
            anchors,
            anchor_lookup,
            subdivisions,
            duplicate_anchor_rows,
        )
    target = -40.0 - math.log2(profile_count(GROUP_BITS, N))
    covered: dict[str, dict[str, Any]] = {}
    wave_rows = []
    print(
        f"g=4 roots={len(roots)} initial_leaves={len(active)} "
        f"atlas_insertions={atlas_insertion_audit['inserted_anchors']} "
        f"witnesses={len(witnesses)} "
        f"target={target:.12f} arithmetic=DIAGNOSTIC_BINARY64",
        flush=True,
    )

    for wave in range(args.max_waves + 1):
        failures: list[tuple[Cell, dict[str, Any]]] = []
        newly_covered = 0
        for cell in active:
            result = evaluate_cell(
                cell,
                anchors,
                evaluator,
                target,
                args.safety_bits,
                args.mixture_candidates_per_vertex,
                args.mixture_weight_tolerance,
            )
            if not result["passes_uniform_target"]:
                failures.append((cell, result))
                continue
            covered[cell.identifier] = {
                **cell_geometry(cell),
                **cell_bound_record(
                    cell,
                    anchors,
                    result,
                    witnesses,
                    profile_count(GROUP_BITS, N),
                ),
                "target_log2": target,
                "margin_bits": target - result["score"],
            }
            newly_covered += 1

        wave_rows.append(
            {
                "wave": wave,
                "evaluated_leaves": len(active),
                "newly_covered_leaves": newly_covered,
                "failed_leaves": len(failures),
                "total_covered_terminal_leaves": len(covered),
                "anchors": len(anchors),
            }
        )
        report = render_report(
            anchors=anchors,
            roots=roots,
            subdivisions=subdivisions,
            covered=covered,
            failed=failures,
            witnesses=witnesses,
            sources=sources,
            target=target,
            safety_bits=args.safety_bits,
            wave_rows=wave_rows,
            atlas_insertion_audit=atlas_insertion_audit,
            root_mesh_audit=root_mesh_audit,
        )
        write_json(args.output, report)
        if args.checkpoint_dir:
            write_json(args.checkpoint_dir / f"wave_{wave:04d}.json", report)
        print(
            f"wave={wave} evaluated={len(active)} newly_covered={newly_covered} "
            f"failed={len(failures)} anchors={len(anchors)}",
            flush=True,
        )
        if not failures or wave == args.max_waves:
            break
        projected_leaves = len(covered) + 5 * len(failures)
        if projected_leaves > args.max_leaves:
            print(
                f"stellar_stop=projected_leaves_{projected_leaves}_exceeds_{args.max_leaves}",
                flush=True,
            )
            break
        next_active = []
        for cell, _score in failures:
            centroid = add_centroid(cell, anchors, anchor_lookup)
            children = subdivide(cell, centroid)
            subdivisions.append(
                {
                    "parent_cell_id": cell.identifier,
                    "insertion_anchor_index": centroid,
                    "insertion_kind": "exact_barycenter",
                    "children": [cell_geometry(child) for child in children],
                }
            )
            next_active.extend(children)
        active = next_active

    print(f"complete_cover={report['complete_cover']}")
    print(f"covered_leaves={report['summary']['covered_leaves']}")
    print(f"failed_leaves={report['summary']['failed_leaves']}")
    print(f"residual_centroid_profiles={report['summary']['unique_residual_centroid_profiles']}")
    print("status=DIAGNOSTIC_BINARY64_G4_RECURSIVE_EXACT_STELLAR_COVER")


if __name__ == "__main__":
    main()
