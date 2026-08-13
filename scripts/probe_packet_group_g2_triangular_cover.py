#!/usr/bin/env python3
"""Adaptive triangular fixed-witness cover for the generalized ``g=2`` profile.

Write an integer packet profile as ``(a0,a1,a2)`` with

    a0 = M - a1 - a2,  M = N / 2.

For every frozen witness used here, its value on a profile is either affine,
or has the form

    constant - charge . profile - log2 Q(profile),

where ``Q`` is the exact packet-profile multinomial normalization.  The latter
function is convex on the nonnegative real simplex because its nonlinear part
is a sum of ``log Gamma(a_i+1)`` terms.  Consequently, one *fixed* witness that
is below the target at all three vertices of a triangle is below the target
throughout that triangle.  The same is true for one fixed convex mixture of
witnesses.  The script first tries the single-witness fast path, then a small
diagnostic LP whose one weight vector is held unchanged at all three vertices
and throughout the triangle.  Vertex-specific witness choices or mixture
weights are forbidden.

By default only profiles with ``a0,a1,a2 >= 1`` are covered.  The three exact
binary edges are emitted as explicit external-certificate hooks.  Set
``--interior-floor 0`` to include them in the triangular search instead.

The default initial mesh is a deterministic Delaunay triangulation of all
feasible integer witness anchors and the exact hull vertices.  Floating-point
Delaunay is used only for discovery: exact integer orientation, total area,
anchor-use, interior-edge multiplicity, and boundary-interval conformity are
checked before traversal.  ``--initial-mesh coarse`` retains the original
two- or three-triangle hull fan.

Arithmetic is binary64, matching ``probe_packet_group_outer_bijection_cover``.
The resulting ledger is therefore a diagnostic convex-cover artifact, not an
outward-rounded theorem certificate.  ``--safety-bits`` reserves a numerical
margin but is not a substitute for directed rounding.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import Delaunay, QhullError

from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import D, K, N, atom_count, normalization_log2


GROUP_BITS = 2
MINIMUM_WEIGHT = 21


Point = tuple[int, int]


def cross(left: Point, middle: Point, right: Point) -> int:
    return (middle[0] - left[0]) * (right[1] - left[1]) - (
        middle[1] - left[1]
    ) * (right[0] - left[0])


def canonical_triangle(vertices: Iterable[Point]) -> tuple[Point, Point, Point]:
    rows = tuple(vertices)
    if len(rows) != 3:
        raise ValueError("a triangle must have three vertices")
    orientation = cross(rows[0], rows[1], rows[2])
    if orientation == 0:
        raise ValueError("degenerate triangle")
    if orientation < 0:
        rows = (rows[0], rows[2], rows[1])
    first = min(range(3), key=lambda index: rows[index])
    return rows[first:] + rows[:first]


def area2(vertices: tuple[Point, Point, Point]) -> int:
    return cross(*vertices)


def edge_steps(left: Point, right: Point) -> int:
    return math.gcd(abs(right[0] - left[0]), abs(right[1] - left[1]))


def lattice_counts(vertices: tuple[Point, Point, Point]) -> tuple[int, int, int]:
    """Return ``(boundary, interior, total)`` using Pick's theorem."""

    boundary = sum(
        edge_steps(vertices[index], vertices[(index + 1) % 3])
        for index in range(3)
    )
    numerator = area2(vertices) - boundary + 2
    if numerator < 0 or numerator % 2:
        raise RuntimeError("triangle violates the lattice form of Pick's theorem")
    interior = numerator // 2
    return boundary, interior, boundary + interior


def point_in_triangle_inclusive(
    vertices: tuple[Point, Point, Point], point: Point
) -> bool:
    """Test membership in a canonical counter-clockwise triangle exactly."""

    return all(
        cross(vertices[index], vertices[(index + 1) % 3], point) >= 0
        for index in range(3)
    )


def enumerate_triangle_lattice_points(
    vertices: tuple[Point, Point, Point],
    expected_total: int,
    scan_limit: int,
    max_points: int,
) -> tuple[Point, ...] | None:
    """Enumerate a triangle's complete integer lattice intersection exactly.

    ``None`` means that a configured resource bound prevented enumeration; it
    never means that a sample was treated as exhaustive.  A Pick-count check
    makes an incomplete scan a hard error.
    """

    vertices = canonical_triangle(vertices)
    if expected_total > max_points:
        return None
    if expected_total == 3:
        # Pick's theorem proves that a unimodular triangle contains only its
        # three vertices, avoiding even a coordinate-span scan.
        return tuple(sorted(vertices))

    spans = (
        max(point[0] for point in vertices) - min(point[0] for point in vertices),
        max(point[1] for point in vertices) - min(point[1] for point in vertices),
    )
    axis = 0 if spans[0] <= spans[1] else 1
    if spans[axis] > scan_limit:
        return None

    first = min(point[axis] for point in vertices)
    last = max(point[axis] for point in vertices)
    points: set[Point] = set()
    for coordinate in range(first, last + 1):
        intersections: list[Fraction] = []
        for index in range(3):
            left = vertices[index]
            right = vertices[(index + 1) % 3]
            left_axis = left[axis]
            right_axis = right[axis]
            other_axis = 1 - axis
            if left_axis == right_axis:
                if coordinate == left_axis:
                    intersections.extend(
                        (Fraction(left[other_axis]), Fraction(right[other_axis]))
                    )
                continue
            low, high = sorted((left_axis, right_axis))
            if low <= coordinate <= high:
                intersections.append(
                    Fraction(left[other_axis])
                    + Fraction(
                        (coordinate - left_axis)
                        * (right[other_axis] - left[other_axis]),
                        right_axis - left_axis,
                    )
                )
        if len(intersections) < 2:
            continue
        lower = ceil_fraction(min(intersections))
        upper = floor_fraction(max(intersections))
        for other in range(lower, upper + 1):
            point = (
                (coordinate, other) if axis == 0 else (other, coordinate)
            )
            if point_in_triangle_inclusive(vertices, point):
                points.add(point)
        if len(points) > expected_total:
            raise RuntimeError("lattice enumeration exceeded the exact Pick count")

    result = tuple(sorted(points))
    if len(result) != expected_total:
        raise RuntimeError(
            "lattice enumeration disagrees with the exact Pick count: "
            f"enumerated={len(result)} expected={expected_total}"
        )
    return result


def compact_json_sha256(value) -> str:
    payload = json.dumps(value, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def strictly_inside(vertices: tuple[Point, Point, Point], point: Point) -> bool:
    return all(
        cross(vertices[index], vertices[(index + 1) % 3], point) > 0
        for index in range(3)
    )


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def scan_vertical(
    vertices: tuple[Point, Point, Point], coordinate: int
) -> Point | None:
    intersections: list[Fraction] = []
    for index in range(3):
        left = vertices[index]
        right = vertices[(index + 1) % 3]
        if left[0] == right[0]:
            if coordinate == left[0]:
                intersections.extend((Fraction(left[1]), Fraction(right[1])))
            continue
        low, high = sorted((left[0], right[0]))
        if low <= coordinate <= high:
            intersections.append(
                Fraction(left[1])
                + Fraction(
                    (coordinate - left[0]) * (right[1] - left[1]),
                    right[0] - left[0],
                )
            )
    if len(intersections) < 2:
        return None
    lower = min(intersections)
    upper = max(intersections)
    first = floor_fraction(lower) + 1
    last = ceil_fraction(upper) - 1
    if first > last:
        return None
    candidates = (first, last, (first + last) // 2)
    for second in candidates:
        point = (coordinate, second)
        if strictly_inside(vertices, point):
            return point
    return None


def find_interior_point(
    vertices: tuple[Point, Point, Point], scan_limit: int
) -> Point | None:
    """Find an integer point strictly inside a non-unimodular triangle."""

    sum_x = sum(point[0] for point in vertices)
    sum_y = sum(point[1] for point in vertices)
    center_x = sum_x // 3
    center_y = sum_y // 3
    for radius in range(9):
        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                if max(abs(dx), abs(dy)) != radius:
                    continue
                candidate = (center_x + dx, center_y + dy)
                if strictly_inside(vertices, candidate):
                    return candidate

    minimum_x = min(point[0] for point in vertices)
    maximum_x = max(point[0] for point in vertices)
    minimum_y = min(point[1] for point in vertices)
    maximum_y = max(point[1] for point in vertices)
    x_span = maximum_x - minimum_x
    y_span = maximum_y - minimum_y
    if min(x_span, y_span) > scan_limit:
        return None
    if x_span <= y_span:
        order = sorted(
            range(minimum_x + 1, maximum_x),
            key=lambda value: abs(3 * value - sum_x),
        )
        for value in order:
            candidate = scan_vertical(vertices, value)
            if candidate is not None:
                return candidate
        return None

    swapped = canonical_triangle((point[1], point[0]) for point in vertices)
    order = sorted(
        range(minimum_y + 1, maximum_y),
        key=lambda value: abs(3 * value - sum_y),
    )
    for value in order:
        candidate = scan_vertical(swapped, value)
        if candidate is not None:
            restored = (candidate[1], candidate[0])
            if strictly_inside(vertices, restored):
                return restored
    return None


def subdivide(
    vertices: tuple[Point, Point, Point], scan_limit: int
) -> tuple[tuple[Point, Point, Point], ...] | None:
    """Partition a lattice triangle into two or three smaller triangles."""

    edges = []
    for left_index, right_index, opposite_index in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        left = vertices[left_index]
        right = vertices[right_index]
        steps = edge_steps(left, right)
        if steps > 1:
            length2 = (right[0] - left[0]) ** 2 + (right[1] - left[1]) ** 2
            edges.append((length2, steps, left, right, vertices[opposite_index]))
    if edges:
        _length2, steps, left, right, opposite = max(edges)
        offset = steps // 2
        point = (
            left[0] + (right[0] - left[0]) * offset // steps,
            left[1] + (right[1] - left[1]) * offset // steps,
        )
        return (
            canonical_triangle((opposite, left, point)),
            canonical_triangle((opposite, point, right)),
        )

    _boundary, interior, _total = lattice_counts(vertices)
    if interior == 0:
        return None
    point = find_interior_point(vertices, scan_limit)
    if point is None:
        return None
    children = tuple(
        canonical_triangle(
            (point, vertices[index], vertices[(index + 1) % 3])
        )
        for index in range(3)
    )
    if sum(area2(child) for child in children) != area2(vertices):
        raise RuntimeError("interior subdivision does not preserve area")
    return children


def polygon_area2(vertices: tuple[Point, ...]) -> int:
    value = sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[index][1] * vertices[(index + 1) % len(vertices)][0]
        for index in range(len(vertices))
    )
    if value <= 0:
        raise ValueError("initial hull must be strictly counterclockwise")
    return value


def point_on_segment(point: Point, left: Point, right: Point) -> bool:
    if cross(left, right, point) != 0:
        return False
    return (
        min(left[0], right[0]) <= point[0] <= max(left[0], right[0])
        and min(left[1], right[1]) <= point[1] <= max(left[1], right[1])
    )


def point_in_convex_polygon(point: Point, hull: tuple[Point, ...]) -> bool:
    return all(
        cross(hull[index], hull[(index + 1) % len(hull)], point) >= 0
        for index in range(len(hull))
    )


def segment_parameter(point: Point, left: Point, right: Point) -> Fraction:
    if right[0] != left[0]:
        return Fraction(point[0] - left[0], right[0] - left[0])
    return Fraction(point[1] - left[1], right[1] - left[1])


def audit_initial_mesh(
    hull: tuple[Point, ...],
    anchors: tuple[Point, ...],
    triangles: tuple[tuple[Point, Point, Point], ...],
) -> dict:
    """Exact integer audit of a floating-point discovery triangulation."""

    if not triangles:
        raise ValueError("initial mesh contains no triangles")
    polygon_area = polygon_area2(hull)
    seen_triangles = set()
    used_points = set()
    edge_multiplicity: dict[tuple[Point, Point], int] = {}
    triangle_area = 0
    for triangle in triangles:
        if area2(triangle) <= 0:
            raise ValueError("initial mesh has a nonpositive triangle orientation")
        key = tuple(sorted(triangle))
        if key in seen_triangles:
            raise ValueError("initial mesh contains a duplicate triangle")
        seen_triangles.add(key)
        if not all(point_in_convex_polygon(point, hull) for point in triangle):
            raise ValueError("initial mesh has a vertex outside the exact hull")
        used_points.update(triangle)
        triangle_area += area2(triangle)
        for index in range(3):
            edge = tuple(sorted((triangle[index], triangle[(index + 1) % 3])))
            edge_multiplicity[edge] = edge_multiplicity.get(edge, 0) + 1
    if triangle_area != polygon_area:
        raise ValueError(
            f"initial mesh area mismatch: triangles={triangle_area} hull={polygon_area}"
        )
    missing = sorted(set(anchors) - used_points)
    if missing:
        raise ValueError(f"initial mesh omitted {len(missing)} anchor points")

    boundary_by_side: list[list[tuple[Fraction, Fraction]]] = [
        [] for _ in hull
    ]
    interior_edges = 0
    boundary_edges = 0
    for edge, multiplicity in edge_multiplicity.items():
        if multiplicity == 2:
            interior_edges += 1
            continue
        if multiplicity != 1:
            raise ValueError("initial mesh edge multiplicity is not one or two")
        matching = [
            index
            for index in range(len(hull))
            if point_on_segment(edge[0], hull[index], hull[(index + 1) % len(hull)])
            and point_on_segment(edge[1], hull[index], hull[(index + 1) % len(hull)])
        ]
        if len(matching) != 1:
            raise ValueError("initial mesh boundary edge is not on one exact hull side")
        side = matching[0]
        left = hull[side]
        right = hull[(side + 1) % len(hull)]
        parameters = sorted(
            (segment_parameter(edge[0], left, right), segment_parameter(edge[1], left, right))
        )
        boundary_by_side[side].append((parameters[0], parameters[1]))
        boundary_edges += 1

    for side, intervals in enumerate(boundary_by_side):
        cursor = Fraction(0)
        for low, high in sorted(intervals):
            if low != cursor or high <= low:
                raise ValueError(f"initial mesh boundary side {side} has a gap or overlap")
            cursor = high
        if cursor != 1:
            raise ValueError(f"initial mesh boundary side {side} has an uncovered tail")
    return {
        "exact_hull_area2": polygon_area,
        "exact_triangle_area2_sum": triangle_area,
        "anchor_points": len(anchors),
        "triangles": len(triangles),
        "interior_edges": interior_edges,
        "boundary_edges": boundary_edges,
        "all_anchors_used": True,
        "orientation_check": "PASS",
        "area_check": "PASS",
        "edge_conformity_check": "PASS",
    }


def feasible_anchor_points(
    witnesses: list["Witness"],
    hull: tuple[Point, ...],
    atoms: int,
    floor: int,
) -> tuple[tuple[Point, ...], dict]:
    exact_pentagon = ((21, 0), (atoms, 0), (0, atoms), (0, 11), (1, 10))
    raw = [*hull, *exact_pentagon]
    raw.extend(row.anchor for row in witnesses if row.anchor is not None)
    accepted = set()
    rejected = 0
    for point in raw:
        a1, a2 = point
        profile = (atoms - a1 - a2, a1, a2)
        if (
            min(profile) < floor
            or a1 + 2 * a2 < MINIMUM_WEIGHT
            or not point_in_convex_polygon(point, hull)
        ):
            rejected += 1
            continue
        accepted.add(point)
    accepted.update(hull)
    return tuple(sorted(accepted)), {
        "raw_anchor_candidates": len(raw),
        "rejected_anchor_candidates": rejected,
        "unique_feasible_anchor_points": len(accepted),
    }


def build_initial_mesh(
    mode: str,
    hull: tuple[Point, ...],
    witnesses: list["Witness"],
    atoms: int,
    floor: int,
) -> tuple[tuple[tuple[Point, Point, Point], ...], dict]:
    if mode == "coarse":
        anchors = tuple(hull)
        triangles = tuple(
            canonical_triangle((hull[0], hull[index], hull[index + 1]))
            for index in range(1, len(hull) - 1)
        )
        stats = audit_initial_mesh(hull, anchors, triangles)
        return triangles, {
            "mode": mode,
            "anchors_a1_a2": [list(point) for point in anchors],
            **stats,
        }

    anchors, anchor_stats = feasible_anchor_points(witnesses, hull, atoms, floor)
    coordinates = np.asarray(anchors, dtype=np.float64)
    try:
        discovery = Delaunay(coordinates)
    except QhullError as error:
        raise RuntimeError(f"anchor Delaunay discovery failed: {error}") from error
    triangles = tuple(
        sorted(
            canonical_triangle(tuple(anchors[int(index)] for index in simplex))
            for simplex in discovery.simplices
        )
    )
    stats = audit_initial_mesh(hull, anchors, triangles)
    return triangles, {
        "mode": mode,
        "discovery": "scipy.spatial.Delaunay on sorted integer anchors",
        "exact_postcheck": "orientation, area, anchor use, and edge conformity",
        "anchors_a1_a2": [list(point) for point in anchors],
        **anchor_stats,
        **stats,
    }


@dataclass(frozen=True)
class Witness:
    identifier: str
    reference: str
    name: str
    kind: str
    constant_log2: float
    charge: tuple[float, float, float]
    support: frozenset[int]
    subtract_normalization: bool
    source: str
    anchor: Point | None

    def eligible(self, profile: tuple[int, int, int]) -> bool:
        return all(count == 0 or index in self.support for index, count in enumerate(profile))

    def value(self, profile: tuple[int, int, int], normalization: float) -> float:
        if not self.eligible(profile):
            return math.inf
        linear = math.fsum(
            count * coefficient for count, coefficient in zip(profile, self.charge)
        )
        result = self.constant_log2 - linear
        if self.subtract_normalization:
            result -= normalization
        return result


@dataclass(frozen=True)
class Triangle:
    vertices: tuple[Point, Point, Point]
    depth: int
    serial: int


@dataclass(order=True)
class QueueItem:
    priority: int
    serial: int
    triangle: Triangle = field(compare=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows_from_json(path: Path) -> tuple[int, list[dict]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(data, list):
        rows = data
        group = int(rows[0].get("group_bits", GROUP_BITS)) if rows else GROUP_BITS
    else:
        rows = list(data.get("rows", []))
        default_group = rows[0].get("group_bits", GROUP_BITS) if rows else GROUP_BITS
        group = int(data.get("group_bits", default_group))
    return group, rows


def row_anchor(row: dict, atoms: int) -> Point | None:
    raw = row.get("profile")
    if not isinstance(raw, list) or len(raw) != 3:
        return None
    try:
        profile = tuple(int(value) for value in raw)
    except (TypeError, ValueError):
        return None
    if any(value < 0 for value in profile) or sum(profile) != atoms:
        return None
    return profile[1], profile[2]


def load_witnesses(
    outer_atlas: Path, combined_atlases: list[Path]
) -> tuple[list[Witness], list[dict]]:
    witnesses: list[Witness] = []
    sources = []
    atoms = atom_count(GROUP_BITS)
    group, rows = rows_from_json(outer_atlas)
    if group != GROUP_BITS:
        raise SystemExit("triangular cover: outer atlas must have group_bits=2")
    sources.append({"path": str(outer_atlas), "sha256": sha256(outer_atlas), "kind": "outer"})
    for index, row in enumerate(rows):
        charge = tuple(float(value) for value in row["charge"])
        if len(charge) != 3 or not all(math.isfinite(value) for value in charge):
            raise SystemExit("triangular cover: invalid outer witness charge")
        witnesses.append(
            Witness(
                f"outer:{index}",
                f"{outer_atlas.name}:{index}",
                str(row.get("name", f"outer_{index}")),
                "outer",
                float(row["constant_log2"]),
                charge,
                frozenset(range(3)),
                False,
                str(outer_atlas),
                row_anchor(row, atoms),
            )
        )

    for path in combined_atlases:
        group, rows = rows_from_json(path)
        if group != GROUP_BITS:
            raise SystemExit("triangular cover: combined atlas must have group_bits=2")
        sources.append({"path": str(path), "sha256": sha256(path), "kind": "combined"})
        for index, row in enumerate(rows):
            charge = tuple(float(value) for value in row["charge"])
            if len(charge) != 3 or not all(math.isfinite(value) for value in charge):
                raise SystemExit("triangular cover: invalid combined witness charge")
            support = frozenset(int(value) for value in row.get("support", range(3)))
            if not support.issubset(range(3)):
                raise SystemExit("triangular cover: invalid combined witness support")
            witnesses.append(
                Witness(
                    f"combined:{path.name}:{index}",
                    f"{path.name}:{index}",
                    str(row.get("name", f"combined_{path.stem}_{index}")),
                    "combined",
                    float(row["constant_log2"]),
                    charge,
                    support,
                    True,
                    str(path),
                    row_anchor(row, atoms),
                )
            )

    bijection_constant = K + N * math.log2(11) - (N - D) * math.log2(10)
    witnesses.append(
        Witness(
            "full_bijection",
            "full_bijection",
            "full_bijection",
            "bijection",
            bijection_constant,
            (0.0, 0.0, 0.0),
            frozenset(range(3)),
            True,
            "built-in formula from probe_packet_group_outer_bijection_cover.py",
            None,
        )
    )
    if not witnesses:
        raise SystemExit("triangular cover: no witnesses loaded")
    return witnesses, sources


class Evaluator:
    def __init__(self, atoms: int, witnesses: list[Witness]):
        self.atoms = atoms
        self.witnesses = witnesses
        self.cache: dict[Point, tuple[tuple[int, int, int], float, tuple[float, ...]]] = {}

    def evaluate(self, point: Point):
        if point in self.cache:
            return self.cache[point]
        a1, a2 = point
        profile = (self.atoms - a1 - a2, a1, a2)
        if min(profile) < 0:
            raise RuntimeError("triangular cover generated a profile outside the simplex")
        normalization = float(
            normalization_log2(GROUP_BITS, np.asarray(profile, dtype=np.float64))[0]
        )
        values = tuple(witness.value(profile, normalization) for witness in self.witnesses)
        result = (profile, normalization, values)
        self.cache[point] = result
        return result


def mixture_candidates(
    vertex_values: np.ndarray, candidates_per_vertex: int
) -> np.ndarray:
    """Return a small, nondominated candidate set for the three-row LP.

    Restricting the candidate set can only lose possible coverage; it cannot
    create a false certificate.  The default keeps the union of the best
    witnesses at each vertex and then removes componentwise dominated rows.
    """

    finite = np.flatnonzero(np.all(np.isfinite(vertex_values), axis=0))
    if not len(finite):
        return finite
    if candidates_per_vertex > 0 and len(finite) > candidates_per_vertex:
        chosen = set()
        for vertex in range(3):
            order = finite[np.argsort(vertex_values[vertex, finite])]
            chosen.update(map(int, order[:candidates_per_vertex]))
        finite = np.asarray(sorted(chosen), dtype=np.int64)

    rows = vertex_values[:, finite]
    keep = np.ones(len(finite), dtype=bool)
    for local in range(len(finite)):
        if not keep[local]:
            continue
        dominated = np.all(rows <= rows[:, local, None], axis=0) & np.any(
            rows < rows[:, local, None], axis=0
        )
        dominated[local] = False
        if np.any(dominated):
            keep[local] = False
    return finite[keep]


def solve_fixed_mixture(
    vertex_values: np.ndarray,
    candidate_indices: np.ndarray,
    target: float,
    safety_bits: float,
    max_support: int,
    weight_tolerance: float,
) -> dict | None:
    """Minimize the worst of three vertex values over one convex mixture."""

    if not len(candidate_indices):
        return None
    values = vertex_values[:, candidate_indices]
    # Subtract the target before solving to keep the small LP numerically
    # centered even when individual witness exponents are very large.
    shifted = values - target
    variables = len(candidate_indices)
    result = linprog(
        np.append(np.zeros(variables), 1.0),
        A_ub=np.hstack((shifted, -np.ones((3, 1)))),
        b_ub=np.zeros(3),
        A_eq=np.asarray([[1.0] * variables + [0.0]]),
        b_eq=np.asarray([1.0]),
        bounds=[(0.0, None)] * variables + [(None, None)],
        # Dual simplex returns a basic solution, normally using no more than
        # three or four component witnesses for this three-vertex problem.
        method="highs-ds",
    )
    if result.status == 2:
        return None
    if not result.success:
        raise RuntimeError(f"triangular mixture LP failed: {result.message}")

    active = np.flatnonzero(result.x[:variables] > weight_tolerance)
    if not len(active) or len(active) > max_support:
        return None
    selected = candidate_indices[active]
    raw_fractions = [Fraction.from_float(float(result.x[index])) for index in active]
    total = sum(raw_fractions, Fraction(0))
    if total <= 0:
        return None
    weights = tuple(value / total for value in raw_fractions)
    float_weights = tuple(float(value) for value in weights)
    component_values = vertex_values[:, selected]
    mixture_values = tuple(
        math.fsum(
            weight * float(component_values[vertex, local])
            for local, weight in enumerate(float_weights)
        )
        for vertex in range(3)
    )
    score = max(mixture_values)
    if score > target - safety_bits:
        return None
    return {
        "indices": tuple(map(int, selected)),
        "weights": weights,
        "weights_binary64": float_weights,
        "component_values": tuple(
            tuple(map(float, component_values[vertex])) for vertex in range(3)
        ),
        "mixture_values": mixture_values,
        "score": score,
        "lp_objective_gap": float(result.fun),
    }


def triangle_certificate(
    triangle: Triangle,
    evaluator: Evaluator,
    target: float,
    safety_bits: float,
    *,
    mixtures: bool,
    candidates_per_vertex: int,
    max_mixture_support: int,
    mixture_weight_tolerance: float,
):
    rows = [evaluator.evaluate(point) for point in triangle.vertices]
    physical_weights = [profile[1] + 2 * profile[2] for profile, _normalization, _values in rows]
    if max(physical_weights) < MINIMUM_WEIGHT:
        return "excluded", None, (), -math.inf, None
    best_index = None
    best_score = math.inf
    best_values: tuple[float, ...] = ()
    for index in range(len(evaluator.witnesses)):
        values = tuple(row[2][index] for row in rows)
        score = max(values)
        if score < best_score:
            best_index = index
            best_score = score
            best_values = values
    if best_index is not None and best_score <= target - safety_bits:
        return "covered", best_index, best_values, best_score, None
    if mixtures:
        vertex_values = np.asarray([row[2] for row in rows], dtype=np.float64)
        candidates = mixture_candidates(vertex_values, candidates_per_vertex)
        mixture = solve_fixed_mixture(
            vertex_values,
            candidates,
            target,
            safety_bits,
            max_mixture_support,
            mixture_weight_tolerance,
        )
        if mixture is not None:
            return "covered_mixture", None, (), mixture["score"], mixture
    return "unresolved", best_index, best_values, best_score, None


def point_residual(
    point: Point,
    evaluator: Evaluator,
    target: float,
    safety_bits: float,
) -> dict | None:
    profile, _normalization, values = evaluator.evaluate(point)
    physical_weight = profile[1] + 2 * profile[2]
    if physical_weight < MINIMUM_WEIGHT:
        return None
    best_index = int(np.argmin(np.asarray(values)))
    best_value = values[best_index]
    if best_value <= target - safety_bits:
        return None
    witness = evaluator.witnesses[best_index]
    return {
        "a1": point[0],
        "a2": point[1],
        "profile": list(profile),
        "physical_weight": physical_weight,
        "best_witness_id": witness.identifier,
        "best_witness_name": witness.name,
        "best_value_log2": best_value,
        "target_log2": target,
        "gap_bits": best_value - target,
    }


def best_point_assignment(
    point: Point,
    evaluator: Evaluator,
    target: float,
    safety_bits: float,
) -> dict | None:
    """Return the best eligible safe single-witness assignment, if one exists."""

    profile, _normalization, values = evaluator.evaluate(point)
    if profile[1] + 2 * profile[2] < MINIMUM_WEIGHT:
        return None
    eligible = [index for index, value in enumerate(values) if math.isfinite(value)]
    if not eligible:
        return None
    best_index = min(eligible, key=lambda index: (values[index], index))
    best_value = values[best_index]
    if best_value > target - safety_bits:
        return None
    witness = evaluator.witnesses[best_index]
    return {
        "profile": list(profile),
        "a1": point[0],
        "a2": point[1],
        "physical_weight": profile[1] + 2 * profile[2],
        "certificate_type": "single_witness",
        "witness_id": witness.identifier,
        "witness_reference": witness.reference,
        "witness": witness.reference,
        "witness_name": witness.name,
        "witness_kind": witness.kind,
        "value_log2": best_value,
        "target_log2": target,
        "margin_bits": target - best_value,
    }


TERMINAL_NONCERTIFICATE_FIELDS = {
    "certificate_type",
    "witness_id",
    "witness_reference",
    "witness_name",
    "witness_kind",
    "vertex_values_log2",
    "witness",
    "mixture_witness_ids",
    "mixture_witness_references",
    "mixture_witness_names",
    "mixture_witness_kinds",
    "mixture_weights",
    "mixture_weights_exact",
    "mixture_weights_binary64",
    "vertex_component_values_log2",
    "vertex_mixture_values_log2",
    "mixture_lp_objective_gap_log2",
}


def terminal_triangle_record(row: dict, profiles: list[list[int]]) -> dict:
    """Copy unresolved geometry without presenting its best witness as a proof."""

    record = {
        key: value for key, value in row.items() if key not in TERMINAL_NONCERTIFICATE_FIELDS
    }
    record["vertices"] = record["vertex_profiles"]
    record["terminal_integer_lattice_count"] = len(profiles)
    record["terminal_integer_profiles_sha256"] = compact_json_sha256(profiles)
    return record


def triangle_record(
    triangle: Triangle,
    evaluator: Evaluator,
    witness_index: int | None,
    values: tuple[float, ...],
    score: float,
    reason: str,
    mixture: dict | None = None,
) -> dict:
    boundary, interior, total = lattice_counts(triangle.vertices)
    record = {
        "triangle_id": triangle.serial,
        "depth": triangle.depth,
        "vertices_a1_a2": [list(point) for point in triangle.vertices],
        "vertex_profiles": [list(evaluator.evaluate(point)[0]) for point in triangle.vertices],
        "area2": area2(triangle.vertices),
        "boundary_lattice_points": boundary,
        "interior_lattice_points": interior,
        "total_lattice_points": total,
        "reason": reason,
        "best_max_vertex_value_log2": score,
    }
    if mixture is not None:
        indices = mixture["indices"]
        components = [evaluator.witnesses[index] for index in indices]
        rational_weights = [
            f"{weight.numerator}/{weight.denominator}"
            for weight in mixture["weights"]
        ]
        record.update(
            {
                "certificate_type": "fixed_convex_mixture",
                "vertices": [
                    list(evaluator.evaluate(point)[0]) for point in triangle.vertices
                ],
                "mixture_witness_ids": [row.identifier for row in components],
                "mixture_witness_references": [row.reference for row in components],
                "mixture_witness_names": [row.name for row in components],
                "mixture_witness_kinds": [row.kind for row in components],
                # The exact rational array is canonical for the outward
                # verifier.  The binary64 copy is diagnostic convenience.
                "mixture_weights": rational_weights,
                "mixture_weights_exact": rational_weights,
                "mixture_weights_binary64": list(mixture["weights_binary64"]),
                "vertex_component_values_log2": [
                    list(row) for row in mixture["component_values"]
                ],
                "vertex_mixture_values_log2": list(mixture["mixture_values"]),
                "mixture_lp_objective_gap_log2": mixture["lp_objective_gap"],
            }
        )
    elif witness_index is not None:
        witness = evaluator.witnesses[witness_index]
        record.update(
            {
                "certificate_type": "single_witness",
                "witness_id": witness.identifier,
                "witness_reference": witness.reference,
                "witness_name": witness.name,
                "witness_kind": witness.kind,
                "vertex_values_log2": list(values),
                "vertices": [
                    list(evaluator.evaluate(point)[0]) for point in triangle.vertices
                ],
                "witness": witness.reference,
            }
        )
    return record


def run_self_test() -> None:
    root = canonical_triangle(((0, 0), (4, 0), (0, 4)))
    pending = [root]
    leaves = []
    while pending:
        row = pending.pop()
        children = subdivide(row, 100)
        if children is None:
            leaves.append(row)
        else:
            if sum(area2(child) for child in children) != area2(row):
                raise SystemExit("triangular cover self-test: area mismatch")
            pending.extend(children)
    if sum(area2(row) for row in leaves) != area2(root):
        raise SystemExit("triangular cover self-test: leaf area mismatch")
    points = set()
    for row in leaves:
        boundary, interior, total = lattice_counts(row)
        if (boundary, interior, total) != (3, 0, 3):
            raise SystemExit("triangular cover self-test: non-unimodular leaf")
        enumerated = enumerate_triangle_lattice_points(row, total, 100, 100)
        if enumerated is None or set(enumerated) != set(row):
            raise SystemExit("triangular cover self-test: terminal enumeration mismatch")
        profiles = [[4 - a1 - a2, a1, a2] for a1, a2 in enumerated]
        if compact_json_sha256(profiles) != hashlib.sha256(
            json.dumps(profiles, separators=(",", ":")).encode("utf-8")
        ).hexdigest():
            raise SystemExit("triangular cover self-test: digest mismatch")
        points.update(enumerated)
    expected = {(a1, a2) for a1 in range(5) for a2 in range(5 - a1)}
    if points != expected:
        raise SystemExit("triangular cover self-test: lattice cover mismatch")
    _boundary, _interior, root_total = lattice_counts(root)
    root_enumerated = enumerate_triangle_lattice_points(root, root_total, 100, 100)
    if root_enumerated is None or set(root_enumerated) != expected:
        raise SystemExit("triangular cover self-test: full-cell enumeration mismatch")
    if enumerate_triangle_lattice_points(root, root_total, 100, 14) is not None:
        raise SystemExit("triangular cover self-test: enumeration cap ignored")
    smoke_hull = ((0, 0), (4, 0), (0, 4))
    smoke_anchors = ((0, 0), (0, 4), (1, 1), (2, 1), (4, 0))
    smoke_discovery = Delaunay(np.asarray(smoke_anchors, dtype=np.float64))
    smoke_triangles = tuple(
        canonical_triangle(tuple(smoke_anchors[int(index)] for index in simplex))
        for simplex in smoke_discovery.simplices
    )
    smoke_audit = audit_initial_mesh(smoke_hull, smoke_anchors, smoke_triangles)
    if smoke_audit["exact_hull_area2"] != smoke_audit["exact_triangle_area2_sum"]:
        raise SystemExit("triangular cover self-test: seeded mesh area mismatch")
    mixture = solve_fixed_mixture(
        np.asarray(
            [
                [-4.0, 1.0, 1.0],
                [1.0, -4.0, 1.0],
                [1.0, 1.0, -4.0],
            ]
        ),
        np.arange(3, dtype=np.int64),
        target=0.0,
        safety_bits=0.1,
        max_support=4,
        weight_tolerance=1e-12,
    )
    if mixture is None or len(mixture["indices"]) != 3:
        raise SystemExit("triangular cover self-test: mixture LP did not close")
    if sum(mixture["weights"], Fraction(0)) != 1:
        raise SystemExit("triangular cover self-test: mixture weights are not exact")
    print(f"self_test_terminal_triangles={len(leaves)}")
    print(f"self_test_lattice_points={len(points)}")
    print(f"self_test_seeded_triangles={len(smoke_triangles)}")
    print(
        "self_test_mixture_weights="
        + ",".join(
            f"{value.numerator}/{value.denominator}"
            for value in mixture["weights"]
        )
    )
    print("status=PASS_PACKET_GROUP_G2_TRIANGULAR_COVER_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path)
    parser.add_argument("--combined-atlas", type=Path, action="append", default=[])
    parser.add_argument("--interior-floor", type=int, default=1)
    parser.add_argument(
        "--initial-mesh", choices=("anchors", "coarse"), default="anchors"
    )
    parser.add_argument("--max-nodes", type=int, default=10000)
    parser.add_argument("--max-depth", type=int, default=128)
    parser.add_argument("--lattice-scan-limit", type=int, default=100000)
    parser.add_argument("--max-terminal-lattice-points", type=int, default=100000)
    parser.add_argument("--safety-bits", type=float, default=1e-5)
    parser.add_argument("--disable-mixtures", action="store_true")
    parser.add_argument("--mixture-candidates-per-vertex", type=int, default=32)
    parser.add_argument("--max-mixture-support", type=int, default=4)
    parser.add_argument("--mixture-weight-tolerance", type=float, default=1e-10)
    parser.add_argument("--progress-every", type=int, default=250)
    parser.add_argument("--show-residuals", type=int, default=20)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.atlas is None:
        raise SystemExit("triangular cover: --atlas is required")
    if (
        args.max_nodes <= 0
        or args.max_depth < 0
        or args.lattice_scan_limit <= 0
        or args.max_terminal_lattice_points <= 0
    ):
        raise SystemExit("triangular cover: invalid search limit")
    if args.interior_floor < 0 or args.safety_bits < 0.0:
        raise SystemExit("triangular cover: invalid floor or safety margin")
    if args.mixture_candidates_per_vertex < 0:
        raise SystemExit("triangular cover: invalid mixture candidate limit")
    if not 1 <= args.max_mixture_support <= 4:
        raise SystemExit("triangular cover: mixture support must lie in 1..4")
    if args.mixture_weight_tolerance < 0.0:
        raise SystemExit("triangular cover: invalid mixture weight tolerance")

    atoms = atom_count(GROUP_BITS)
    floor = args.interior_floor
    if 3 * floor > atoms:
        raise SystemExit("triangular cover: interior floor makes the simplex empty")
    witnesses, sources = load_witnesses(args.atlas, args.combined_atlas)
    evaluator = Evaluator(atoms, witnesses)
    target = -40.0 - math.log2(profile_count(GROUP_BITS, N))
    # Exact integer hulls of the physical-weight >=21 region.  With all
    # three classes allowed the hull is a pentagon.  After the exact edges
    # are deferred, its full-support part is a quadrilateral.  Using these
    # hulls avoids spending the adaptive budget on the 121 impossible
    # low-physical-weight profiles near the all-zero corner.
    if floor == 0:
        hull = ((21, 0), (atoms, 0), (0, atoms), (0, 11), (1, 10))
    elif floor == 1:
        hull = ((19, 1), (atoms - 2, 1), (1, atoms - 2), (1, 10))
    else:
        hull = (
            (floor, floor),
            (atoms - 2 * floor, floor),
            (floor, atoms - 2 * floor),
        )
    roots, initial_mesh_stats = build_initial_mesh(
        args.initial_mesh, hull, witnesses, atoms, floor
    )

    serial = -1
    queue = []
    for root in roots:
        serial += 1
        root_triangle = Triangle(root, 0, serial)
        heapq.heappush(queue, QueueItem(-area2(root), serial, root_triangle))
    covered = []
    excluded = []
    unresolved = []
    processed = 0
    print("g=2 adaptive triangular fixed-witness cover", flush=True)
    print(
        f"M={atoms} interior_floor={floor} witnesses={len(witnesses)} "
        f"target={target:.12f} safety_bits={args.safety_bits:.9g}",
        flush=True,
    )
    print(
        f"initial_mesh={args.initial_mesh} anchors={initial_mesh_stats['anchor_points']} "
        f"triangles={initial_mesh_stats['triangles']} exact_mesh_audit=PASS",
        flush=True,
    )
    print(
        "soundness=one_frozen_witness_or_one_fixed_convex_mixture_must_pass_all_vertices",
        flush=True,
    )
    print(
        "arithmetic=DIAGNOSTIC_BINARY64_NOT_OUTWARD_ROUNDED",
        flush=True,
    )
    if floor:
        print(
            "edge_policy=DEFERRED_PROFILES_WITH_ANY_CLASS_BELOW_INTERIOR_FLOOR",
            flush=True,
        )

    while queue and processed < args.max_nodes:
        item = heapq.heappop(queue)
        triangle = item.triangle
        processed += 1
        disposition, witness_index, values, score, mixture = triangle_certificate(
            triangle,
            evaluator,
            target,
            args.safety_bits,
            mixtures=not args.disable_mixtures,
            candidates_per_vertex=args.mixture_candidates_per_vertex,
            max_mixture_support=args.max_mixture_support,
            mixture_weight_tolerance=args.mixture_weight_tolerance,
        )
        if disposition in {"covered", "covered_mixture"}:
            covered.append(
                triangle_record(
                    triangle,
                    evaluator,
                    witness_index,
                    values,
                    score,
                    (
                        "same_witness_safe_at_all_vertices"
                        if mixture is None
                        else "fixed_convex_mixture_safe_at_all_vertices"
                    ),
                    mixture,
                )
            )
            continue
        if disposition == "excluded":
            excluded.append(
                triangle_record(
                    triangle,
                    evaluator,
                    None,
                    (),
                    score,
                    "all_profiles_have_physical_weight_below_21",
                )
            )
            continue
        boundary, interior, total = lattice_counts(triangle.vertices)
        reason = None
        children = None
        if triangle.depth >= args.max_depth:
            reason = "depth_limit"
        elif total == 3:
            reason = "unimodular_integer_leaf"
        else:
            children = subdivide(triangle.vertices, args.lattice_scan_limit)
            if children is None:
                reason = "no_lattice_subdivision_found"
        if reason is not None:
            unresolved.append(
                triangle_record(
                    triangle, evaluator, witness_index, values, score, reason
                )
            )
            continue
        for vertices in children or ():
            serial += 1
            child = Triangle(vertices, triangle.depth + 1, serial)
            heapq.heappush(queue, QueueItem(-area2(vertices), serial, child))
        if args.progress_every and processed % args.progress_every == 0:
            print(
                f"processed={processed} covered={len(covered)} excluded={len(excluded)} "
                f"open={len(queue)} unresolved={len(unresolved)}",
                flush=True,
            )

    for item in queue:
        triangle = item.triangle
        disposition, witness_index, values, score, mixture = triangle_certificate(
            triangle,
            evaluator,
            target,
            args.safety_bits,
            mixtures=not args.disable_mixtures,
            candidates_per_vertex=args.mixture_candidates_per_vertex,
            max_mixture_support=args.max_mixture_support,
            mixture_weight_tolerance=args.mixture_weight_tolerance,
        )
        unresolved.append(
            triangle_record(
                triangle,
                evaluator,
                witness_index,
                values,
                score,
                "node_budget",
            )
        )
    unresolved.sort(
        key=lambda row: (row["best_max_vertex_value_log2"], row["area2"]),
        reverse=True,
    )

    # Close small terminal cells on the integer lattice without claiming that
    # their real interiors are covered.  Boundary points shared by cells are
    # classified once and emitted once in the point-assignment ledger.
    terminal_triangles = []
    nonterminal_unresolved = []
    terminal_cell_profiles: list[tuple[dict, list[list[int]]]] = []
    point_assignments: dict[Point, dict] = {}
    residual_points: dict[Point, dict] = {}
    point_classification: dict[Point, tuple[dict | None, dict | None]] = {}
    terminal_lattice_occurrences = 0
    for row in unresolved:
        vertices = canonical_triangle(
            tuple(tuple(point) for point in row["vertices_a1_a2"])
        )
        points = enumerate_triangle_lattice_points(
            vertices,
            row["total_lattice_points"],
            args.lattice_scan_limit,
            args.max_terminal_lattice_points,
        )
        if points is None:
            nonterminal_unresolved.append(row)
            continue

        profiles = sorted(
            [list(evaluator.evaluate(point)[0]) for point in points]
        )
        terminal_lattice_occurrences += len(profiles)
        assigned_in_cell = 0
        unassigned_in_cell = 0
        for point in points:
            if point not in point_classification:
                assignment = best_point_assignment(
                    point, evaluator, target, args.safety_bits
                )
                residual = None
                if assignment is None:
                    residual = point_residual(
                        point, evaluator, target, args.safety_bits
                    )
                    if residual is None:
                        raise RuntimeError(
                            "terminal hull unexpectedly contains an unassignable "
                            "physical-weight-below-21 profile"
                        )
                point_classification[point] = (assignment, residual)
            assignment, residual = point_classification[point]
            if assignment is not None:
                point_assignments[point] = assignment
                assigned_in_cell += 1
            else:
                assert residual is not None
                residual_points[point] = residual
                unassigned_in_cell += 1

        terminal = terminal_triangle_record(row, profiles)
        terminal["terminal_assigned_integer_lattice_count"] = assigned_in_cell
        terminal["terminal_unassigned_integer_lattice_count"] = unassigned_in_cell
        terminal_triangles.append(terminal)
        terminal_cell_profiles.append((terminal, profiles))

    discrete_assignments = sorted(
        point_assignments.values(), key=lambda row: tuple(row["profile"])
    )
    residuals = sorted(
        residual_points.values(), key=lambda row: row["gap_bits"], reverse=True
    )
    residuals_exhaustive = not nonterminal_unresolved

    terminal_cell_profiles.sort(
        key=lambda item: tuple(tuple(profile) for profile in item[0]["vertices"])
    )
    terminal_cell_digest_rows = [
        {"vertices": row["vertices"], "profiles": profiles}
        for row, profiles in terminal_cell_profiles
    ]
    terminal_cell_lattice_sha256 = compact_json_sha256(terminal_cell_digest_rows)
    unique_terminal_profiles = sorted(
        {tuple(evaluator.evaluate(point)[0]) for point in point_classification}
    )
    terminal_unique_integer_profiles_sha256 = compact_json_sha256(
        [list(profile) for profile in unique_terminal_profiles]
    )

    deferred_regions = []
    if floor == 1:
        deferred_regions = [
            {
                "name": "edge_a0_zero",
                "constraint": "a0=0, a1+a2=M",
                "parameter": "0<=a1<=M, a2=M-a1",
            },
            {
                "name": "edge_a1_zero",
                "constraint": "a1=0, a0+a2=M",
                "parameter": "0<=a2<=M, a0=M-a2",
            },
            {
                "name": "edge_a2_zero",
                "constraint": "a2=0, a0+a1=M",
                "parameter": "0<=a1<=M, a0=M-a1",
            },
        ]
    elif floor > 1:
        deferred_regions = [
            {
                "name": "boundary_slabs",
                "constraint": f"min(a0,a1,a2)<{floor}",
                "note": "includes the exact edges and adjacent integer layers",
            }
        ]

    complete_continuous_interior = not unresolved
    complete_interior_integer = not nonterminal_unresolved and not residuals
    complete_global_continuous = complete_continuous_interior and floor == 0
    complete_global_integer = complete_interior_integer and floor == 0
    single_covered = sum(
        row.get("certificate_type") == "single_witness" for row in covered
    )
    mixture_covered = sum(
        row.get("certificate_type") == "fixed_convex_mixture" for row in covered
    )
    report = {
        "schema": "packet-group-g2-adaptive-triangular-cover-v3",
        "status": "DIAGNOSTIC_BINARY64_G2_ADAPTIVE_TRIANGULAR_FIXED_WITNESS_COVER",
        "parameters": {
            "group_bits": GROUP_BITS,
            "N": N,
            "K": K,
            "D": D,
            "M": atoms,
            "target_log2": target,
            "interior_floor": floor,
            "initial_mesh": args.initial_mesh,
            "safety_bits": args.safety_bits,
            "max_nodes": args.max_nodes,
            "max_depth": args.max_depth,
            "lattice_scan_limit": args.lattice_scan_limit,
            "max_terminal_lattice_points": args.max_terminal_lattice_points,
            "mixtures_enabled": not args.disable_mixtures,
            "mixture_candidates_per_vertex": args.mixture_candidates_per_vertex,
            "max_mixture_support": args.max_mixture_support,
            "mixture_weight_tolerance": args.mixture_weight_tolerance,
        },
        "assumptions": [
            "Each outer atlas row is a valid fixed affine upper bound on its stated g=2 domain.",
            (
                "Each combined atlas row is valid on profiles whose positive "
                "coordinates lie in its support."
            ),
            (
                "The built-in full-bijection formula is the one used by "
                "probe_packet_group_outer_bijection_cover.py."
            ),
            (
                "A triangle is covered only when one unchanged witness is at "
                "most target-safety_bits at all three vertices."
            ),
            (
                "If the single-witness test fails, one fixed nonnegative "
                "convex mixture whose exact rational weights sum to one may "
                "be used unchanged at all three vertices."
            ),
            (
                "Convexity of constant-charge.profile-log2(Q(profile)) then "
                "covers the full real triangle."
            ),
            (
                "Anchor-seeded Delaunay output is discovery data only; exact "
                "integer orientation, area, anchor-use, and boundary-edge "
                "checks must all pass before traversal."
            ),
            "Physical profile weights below 21 contain no nonzero committed EBCH outer word.",
            "Binary64 evaluations are diagnostic and require outward rounding before theorem use.",
            (
                "A discrete point assignment certifies only its named integer "
                "profile and does not cover any real neighborhood of that point."
            ),
        ],
        "sources": sources,
        "initial_mesh_audit": initial_mesh_stats,
        "witnesses": [
            {
                "witness_id": row.identifier,
                "witness_reference": row.reference,
                "name": row.name,
                "kind": row.kind,
                "constant_log2": row.constant_log2,
                "charge": list(row.charge),
                "support": sorted(row.support),
                "subtract_normalization": row.subtract_normalization,
                "source": row.source,
                "anchor_a1_a2": list(row.anchor) if row.anchor is not None else None,
            }
            for row in witnesses
        ],
        "summary": {
            "processed_nodes": processed,
            "covered_triangles": len(covered),
            "single_witness_triangles": single_covered,
            "fixed_mixture_triangles": mixture_covered,
            "excluded_triangles": len(excluded),
            "adaptive_uncovered_triangles": len(unresolved),
            "terminal_triangles": len(terminal_triangles),
            "nonterminal_unresolved_triangles": len(nonterminal_unresolved),
            "terminal_integer_lattice_occurrences": terminal_lattice_occurrences,
            "terminal_unique_integer_profiles": len(unique_terminal_profiles),
            "discrete_point_assignments": len(discrete_assignments),
            "uncovered_integer_residuals": len(residuals),
            "integer_residuals_exhaustive": residuals_exhaustive,
            "complete_interior_triangle_cover": complete_continuous_interior,
            "complete_interior_integer_cover": complete_interior_integer,
            "complete_global_continuous_cover": complete_global_continuous,
            "complete_global_integer_cover": complete_global_integer,
            "deferred_region_count": len(deferred_regions),
        },
        "deferred_edge_or_boundary_hooks": deferred_regions,
        "group_bits": GROUP_BITS,
        # Backward-compatible generic flag retains its continuous meaning.
        "complete_cover": complete_global_continuous,
        "complete_continuous_cover": complete_global_continuous,
        "complete_interior_continuous_cover": complete_continuous_interior,
        "complete_integer_cover": complete_global_integer,
        "complete_global_integer_cover": complete_global_integer,
        "complete_interior_integer_cover": complete_interior_integer,
        "integer_residuals_exhaustive": residuals_exhaustive,
        "terminal_cell_lattice_sha256": terminal_cell_lattice_sha256,
        "terminal_unique_integer_profiles_sha256": (
            terminal_unique_integer_profiles_sha256
        ),
        "terminal_lattice_digest_serialization": {
            "encoding": "UTF-8",
            "json_separators": [",", ":"],
            "trailing_newline": False,
            "profile_order": "lexicographic_[a0,a1,a2]",
            "cell_order": "lexicographic_canonical_CCW_vertex_profiles",
            "terminal_cell_lattice_sha256_scope": "cell-tagged_lattice_occurrences",
            "terminal_unique_integer_profiles_sha256_scope": "globally_deduplicated_profiles",
        },
        "triangles": covered,
        "triangle_to_witness_ledger": covered,
        "terminal_triangles": terminal_triangles,
        "discrete_point_assignments": discrete_assignments,
        "excluded_triangle_ledger": excluded,
        "unresolved_triangles": nonterminal_unresolved,
        "uncovered_integer_residuals": residuals,
    }

    print(f"processed_nodes={processed}")
    print(f"covered_triangles={len(covered)}")
    print(f"single_witness_triangles={single_covered}")
    print(f"fixed_mixture_triangles={mixture_covered}")
    print(f"excluded_triangles={len(excluded)}")
    print(f"adaptive_uncovered_triangles={len(unresolved)}")
    print(f"terminal_triangles={len(terminal_triangles)}")
    print(f"nonterminal_unresolved_triangles={len(nonterminal_unresolved)}")
    print(f"terminal_integer_lattice_occurrences={terminal_lattice_occurrences}")
    print(f"terminal_unique_integer_profiles={len(unique_terminal_profiles)}")
    print(f"discrete_point_assignments={len(discrete_assignments)}")
    print(f"uncovered_integer_residuals={len(residuals)}")
    print("integer_residuals_exhaustive=" + ("YES" if residuals_exhaustive else "NO"))
    print(
        "complete_interior_continuous_cover="
        + ("PASS" if complete_continuous_interior else "NO")
    )
    print(
        "complete_interior_integer_cover="
        + ("PASS" if complete_interior_integer else "NO")
    )
    print(
        "complete_global_continuous_cover="
        + ("PASS" if complete_global_continuous else "NO")
    )
    print(
        "complete_global_integer_cover="
        + ("PASS" if complete_global_integer else "NO")
    )
    for row in residuals[: args.show_residuals]:
        print(
            f"residual_profile={','.join(map(str,row['profile']))} "
            f"best={row['best_value_log2']:.9f} gap={row['gap_bits']:.9f} "
            f"witness={row['best_witness_name']}"
        )
    print("status=DIAGNOSTIC_BINARY64_G2_ADAPTIVE_TRIANGULAR_FIXED_WITNESS_COVER")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
