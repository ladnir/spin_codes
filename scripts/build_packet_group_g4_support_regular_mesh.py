#!/usr/bin/env python3
"""Build exact regular root meshes for every ``g=4`` profile support.

For exact support ``S``, write ``a_i=1+b_i`` on active classes and zero on
inactive classes.  The shifted counts form a simplex of residual mass
``R=M-|S|``, clipped by physical weight at least 21.  Exact hull vertices and
integer atlas anchors are lifted with deterministic rational heights.  Qhull
only proposes lower facets; exact rational support, topology, boundary, and
determinant-volume checks decide the emitted mesh.
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
from typing import Any, Iterable

import numpy as np
from scipy.spatial import ConvexHull, QhullError

import certify_packet_group_g4_anchor_mesh as g4
import probe_packet_group_g4_stellar_cover as cover
from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import N


SCHEMA = "packet-group-g4-support-stratified-exact-regular-mesh-v1"
PERTURB_PRIME = 1_000_003
PERTURB_SCALE = 1 << 20
OWNER_RULE = "minimum_cell_id_within_exact_support_stratum"
UNION_MODE = "cell_local_profile_count_upper"
MIXTURE_CANDIDATES_PER_VERTEX = 64

Profile = tuple[Fraction, ...]
Simplex = tuple[int, ...]


def compact_hash(value: Any) -> str:
    encoded = json.dumps(
        value, separators=(",", ":"), sort_keys=True, ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


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
        raise ValueError(f"support regular mesh: atlas {path} has no rows")
    return rows


def exact_integer_profile(row: dict[str, Any]) -> tuple[int, ...] | None:
    raw = row.get("profile")
    if not isinstance(raw, list) or len(raw) != 5:
        return None
    values = []
    for value in raw:
        if isinstance(value, (bool, float)):
            return None
        try:
            parsed = Fraction(str(value))
        except (ValueError, ZeroDivisionError):
            return None
        if parsed.denominator != 1:
            return None
        values.append(parsed.numerator)
    if any(value < 0 for value in values) or sum(values) != g4.M:
        return None
    if sum(index * value for index, value in enumerate(values)) < 21:
        return None
    return tuple(values)


def load_atlas_anchors(paths: list[Path]) -> tuple[dict[int, set[Profile]], list[dict[str, str]], int]:
    by_mask: dict[int, set[Profile]] = {mask: set() for mask in range(1, 32)}
    duplicate_rows = 0
    sources = []
    for path in paths:
        for row in artifact_rows(path):
            profile = exact_integer_profile(row)
            if profile is None:
                continue
            mask = sum((1 << index) for index, value in enumerate(profile) if value)
            exact = tuple(Fraction(value) for value in profile)
            if exact in by_mask[mask]:
                duplicate_rows += 1
            by_mask[mask].add(exact)
        sources.append({"path": str(path.resolve()), "sha256": file_sha256(path)})
    return by_mask, sources, duplicate_rows


def determinant(matrix: Iterable[Iterable[Fraction]]) -> Fraction:
    return g4.determinant(matrix)


def coordinate_chart(profile: Profile, support: tuple[int, ...]) -> tuple[Fraction, ...]:
    return tuple(profile[index] for index in support[1:])


def simplex_determinant(
    simplex: Simplex, anchors: tuple[Profile, ...], support: tuple[int, ...]
) -> Fraction:
    dimension = len(support) - 1
    if dimension == 0:
        return Fraction(1)
    points = [coordinate_chart(anchors[index], support) for index in simplex]
    return determinant(
        [
            [points[row][column] - points[0][column] for column in range(dimension)]
            for row in range(1, dimension + 1)
        ]
    )


def affine_dimension(
    indices: Iterable[int], anchors: tuple[Profile, ...], support: tuple[int, ...]
) -> int:
    selected = tuple(indices)
    if not selected:
        return -1
    dimension = len(support) - 1
    if len(selected) == 1:
        return 0
    points = [coordinate_chart(anchors[index], support) for index in selected]
    rows = [
        [points[row][column] - points[0][column] for column in range(dimension)]
        for row in range(1, len(points))
    ]
    # Tiny exact Gaussian rank.
    rank = 0
    for column in range(dimension):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [value / scale for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * other
                for value, other in zip(rows[row], rows[rank])
            ]
        rank += 1
    return rank


def integer_hull_vertices(support: tuple[int, ...]) -> tuple[Profile, ...]:
    """Vertices of the exact-support feasible integer-profile hull.

    Only strata containing class zero can be clipped: without class zero,
    every shifted-simplex pure vertex has enormous positive physical weight.
    At a clipped zero corner, a non-pure vertex with positive class-zero mass
    must have shifted physical weight below ``T+max(S)``.  Otherwise shifting
    one unit between class zero and a positive occupied class gives two
    feasible integer neighbors whose midpoint is the point.  The class-zero
    face itself has only the nonzero pure vertices.
    """

    residual = g4.M - len(support)
    base_weight = sum(support)
    threshold = 21 - base_weight
    pure = []
    for favored in support:
        shifted = [0] * 5
        shifted[favored] = residual
        profile = tuple(
            Fraction(shifted[index] + (1 if index in support else 0))
            for index in range(5)
        )
        if sum(index * profile[index] for index in range(5)) >= 21:
            pure.append(profile)
    if 0 not in support or threshold <= 0:
        return tuple(sorted(pure))

    positive = tuple(index for index in support if index)
    if not positive:
        return ()
    maximum_weight = max(positive)
    candidates = set(pure)

    def enumerate_positive(position: int, remaining_weight: int, values: list[int]) -> None:
        if position == len(positive):
            physical = sum(index * value for index, value in zip(positive, values))
            if not threshold <= physical < threshold + maximum_weight:
                return
            used = sum(values)
            if used > residual:
                return
            shifted = [0] * 5
            shifted[0] = residual - used
            for index, value in zip(positive, values):
                shifted[index] = value
            candidates.add(
                tuple(
                    Fraction(shifted[index] + (1 if index in support else 0))
                    for index in range(5)
                )
            )
            return
        weight = positive[position]
        maximum = (remaining_weight + maximum_weight - 1) // weight
        for value in range(maximum + 1):
            values.append(value)
            enumerate_positive(
                position + 1, remaining_weight - weight * value, values
            )
            values.pop()

    enumerate_positive(0, threshold + maximum_weight - 1, [])
    candidate_rows = tuple(sorted(candidates))
    dimension = len(support) - 1
    if dimension == 1:
        return (candidate_rows[0], candidate_rows[-1])
    coordinates = np.asarray(
        [[float(value) for value in coordinate_chart(row, support)] for row in candidate_rows]
    )
    hull = ConvexHull(coordinates, qhull_options="Qx Qt")
    vertices = tuple(candidate_rows[int(index)] for index in sorted(hull.vertices))
    return tuple(sorted(vertices))


def exact_hull_facet_vertex_sets(
    hull_vertices: tuple[Profile, ...], support: tuple[int, ...]
) -> tuple[frozenset[int], ...]:
    dimension = len(support) - 1
    if dimension == 0:
        return ()
    if dimension == 1:
        return (frozenset((0,)), frozenset((len(hull_vertices) - 1,)))
    import itertools

    facets = set()
    for seed in itertools.combinations(range(len(hull_vertices)), dimension):
        if affine_dimension(seed, hull_vertices, support) != dimension - 1:
            continue
        points = [coordinate_chart(hull_vertices[index], support) for index in seed]
        signs = []
        coplanar = set(seed)
        for index, profile in enumerate(hull_vertices):
            if index in coplanar:
                continue
            point = coordinate_chart(profile, support)
            matrix = [
                [points[row][column] - points[0][column] for column in range(dimension)]
                for row in range(1, dimension)
            ]
            matrix.append([point[column] - points[0][column] for column in range(dimension)])
            value = determinant(matrix)
            if value:
                signs.append(1 if value > 0 else -1)
            else:
                coplanar.add(index)
        if signs and (all(value > 0 for value in signs) or all(value < 0 for value in signs)):
            facets.add(frozenset(coplanar))
    return tuple(sorted(facets, key=lambda row: tuple(sorted(row))))


def hull_boundary_sets(
    anchors: tuple[Profile, ...],
    support: tuple[int, ...],
    hull_vertices: tuple[Profile, ...] | None = None,
) -> dict[str, frozenset[int]]:
    hull = anchors if hull_vertices is None else hull_vertices
    facet_vertices = exact_hull_facet_vertex_sets(hull, support)
    dimension = len(support) - 1
    result = {}
    for ordinal, facet in enumerate(facet_vertices):
        import itertools

        seed = next(
            candidate
            for candidate in itertools.combinations(sorted(facet), dimension)
            if affine_dimension(candidate, hull, support) == dimension - 1
        )
        seed_points = [coordinate_chart(hull[index], support) for index in seed]
        members = set()
        for index, profile in enumerate(anchors):
            point = coordinate_chart(profile, support)
            if dimension == 1:
                coplanar = point == seed_points[0]
            else:
                matrix = [
                    [seed_points[row][column] - seed_points[0][column] for column in range(dimension)]
                    for row in range(1, dimension)
                ]
                matrix.append([point[column] - seed_points[0][column] for column in range(dimension)])
                coplanar = determinant(matrix) == 0
            if coplanar:
                members.add(index)
        labels = [
            f"a{coordinate}=1"
            for coordinate in support
            if all(hull[index][coordinate] == 1 for index in facet)
        ]
        label = labels[0] if len(labels) == 1 else f"integer_hull_facet_{ordinal:03d}"
        result[label] = frozenset(members)
    return result


def pulling_triangulation(
    hull_vertices: tuple[Profile, ...], support: tuple[int, ...]
) -> tuple[Simplex, ...]:
    dimension = len(support) - 1
    if dimension == 0:
        return ((0,),)
    boundaries = tuple(hull_boundary_sets(hull_vertices, support, hull_vertices).values())

    def triangulate(face: frozenset[int], local_dimension: int) -> list[Simplex]:
        if local_dimension == 0:
            return [(min(face),)]
        pivot = min(face)
        opposite = {
            candidate
            for boundary in boundaries
            if (candidate := face & boundary)
            and pivot not in candidate
            and affine_dimension(candidate, hull_vertices, support)
            == local_dimension - 1
        }
        if not opposite:
            raise RuntimeError("support pulling triangulation found no opposite facet")
        result = []
        for facet in sorted(opposite, key=lambda row: tuple(sorted(row))):
            for simplex in triangulate(facet, local_dimension - 1):
                result.append((pivot, *simplex))
        return result

    rows = triangulate(frozenset(range(len(hull_vertices))), dimension)
    if len({tuple(sorted(row)) for row in rows}) != len(rows):
        raise RuntimeError("support pulling triangulation produced duplicates")
    return tuple(rows)


def expected_hull_volume(
    hull_vertices: tuple[Profile, ...], support: tuple[int, ...]
) -> Fraction:
    return sum(
        (abs(simplex_determinant(simplex, hull_vertices, support)) for simplex in pulling_triangulation(hull_vertices, support)),
        Fraction(0),
    )


def lifting_height(profile: Profile, support: tuple[int, ...], rank: int) -> Fraction:
    coordinates = coordinate_chart(profile, support)
    base = sum((value / g4.M) ** 2 for value in coordinates)
    residue = pow(rank + 1, 5, PERTURB_PRIME)
    return base + Fraction(residue, PERTURB_PRIME * PERTURB_SCALE)


def solve_linear(matrix, right) -> tuple[Fraction, ...]:
    rows = [list(map(Fraction, row)) + [Fraction(value)] for row, value in zip(matrix, right)]
    size = len(rows)
    for column in range(size):
        pivot = next((row for row in range(column, size) if rows[row][column]), None)
        if pivot is None:
            raise ValueError("singular support lower-facet interpolation")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [value / scale for value in rows[column]]
        for row in range(size):
            if row == column or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * other
                for value, other in zip(rows[row], rows[column])
            ]
    return tuple(row[-1] for row in rows)


def exact_lower_facet(
    simplex: Simplex,
    anchors: tuple[Profile, ...],
    heights: tuple[Fraction, ...],
    support: tuple[int, ...],
) -> Fraction:
    coordinates = [coordinate_chart(anchors[index], support) for index in simplex]
    coefficients = solve_linear(
        [[Fraction(1), *point] for point in coordinates],
        [heights[index] for index in simplex],
    )
    minimum = None
    selected = set(simplex)
    for index, (profile, height) in enumerate(zip(anchors, heights)):
        point = coordinate_chart(profile, support)
        plane = coefficients[0] + sum(
            coefficient * value
            for coefficient, value in zip(coefficients[1:], point)
        )
        slack = height - plane
        if index in selected:
            if slack:
                raise ValueError("support facet interpolation failed")
        elif slack <= 0:
            raise ValueError(
                f"support facet {simplex} is not strict below anchor {index}"
            )
        else:
            minimum = slack if minimum is None else min(minimum, slack)
    return Fraction(0) if minimum is None else minimum


def propose_facets(
    anchors: tuple[Profile, ...], heights: tuple[Fraction, ...], support: tuple[int, ...]
) -> list[Simplex]:
    dimension = len(support) - 1
    if dimension == 0:
        return [(0,)]
    if len(anchors) == dimension + 1:
        return [tuple(range(len(anchors)))]
    if dimension == 1:
        order = sorted(
            range(len(anchors)), key=lambda index: coordinate_chart(anchors[index], support)
        )
        return [tuple(sorted((left, right))) for left, right in zip(order, order[1:])]
    points = np.asarray(
        [
            [
                *(float(value / g4.M) for value in coordinate_chart(profile, support)),
                float(height),
            ]
            for profile, height in zip(anchors, heights)
        ]
    )
    try:
        hull = ConvexHull(points, qhull_options="Qx Qt Qc")
    except QhullError as error:
        raise ValueError(f"support lifted hull proposal failed: {error}") from error
    facets = {
        tuple(sorted(int(index) for index in simplex))
        for simplex, equation in zip(hull.simplices, hull.equations)
        if float(equation[-2]) < 0.0
    }
    return sorted(facets)


def permutation_sign(values: tuple[int, ...]) -> int:
    inversions = sum(
        values[left] > values[right]
        for left in range(len(values))
        for right in range(left + 1, len(values))
    )
    return -1 if inversions % 2 else 1


def verify_mesh(
    anchors: tuple[Profile, ...],
    simplices: list[Simplex],
    support: tuple[int, ...],
    hull_vertices: tuple[Profile, ...],
) -> dict[str, Any]:
    dimension = len(support) - 1
    if dimension == 0:
        if simplices != [(0,)] or anchors != hull_vertices:
            raise ValueError("zero-dimensional support mesh is not its unique profile")
        return {
            "dimension": 0,
            "simplices": 1,
            "determinant_volume": "1",
            "exact_hull_determinant_volume": "1",
            "boundary_facets": {},
        }
    facets: Counter[tuple[int, ...]] = Counter()
    orientations: Counter[tuple[int, ...]] = Counter()
    volume = Fraction(0)
    for simplex in simplices:
        determinant_value = simplex_determinant(simplex, anchors, support)
        if not determinant_value:
            raise ValueError(f"support mesh contains degenerate simplex {simplex}")
        oriented = simplex if determinant_value > 0 else (simplex[1], simplex[0], *simplex[2:])
        volume += abs(determinant_value)
        for omitted in range(dimension + 1):
            ordered = oriented[:omitted] + oriented[omitted + 1 :]
            facet = tuple(sorted(ordered))
            facets[facet] += 1
            orientations[facet] += ((-1) ** omitted) * permutation_sign(ordered)

    boundaries = hull_boundary_sets(anchors, support, hull_vertices)
    boundary_counts = Counter()
    for facet, multiplicity in facets.items():
        if multiplicity == 2:
            if orientations[facet] != 0:
                raise ValueError("support mesh interior orientation mismatch")
        elif multiplicity == 1:
            labels = [label for label, indices in boundaries.items() if set(facet) <= indices]
            if len(labels) != 1:
                raise ValueError(
                    f"support mesh unpaired facet {facet} has {len(labels)} "
                    f"boundary labels {labels}"
                )
            boundary_counts[labels[0]] += 1
        else:
            raise ValueError("support mesh facet multiplicity is not conforming")
    expected = expected_hull_volume(hull_vertices, support)
    if volume != expected:
        raise ValueError(f"support mesh volume mismatch {volume} != {expected}")
    if set(boundary_counts) != set(
        hull_boundary_sets(hull_vertices, support, hull_vertices)
    ):
        raise ValueError("support mesh does not expose every exact hull boundary")
    return {
        "dimension": dimension,
        "simplices": len(simplices),
        "determinant_volume": str(volume),
        "exact_hull_determinant_volume": str(expected),
        "boundary_facets": dict(sorted(boundary_counts.items())),
    }


def pure_residual_profile(support: tuple[int, ...], favored: int) -> Profile:
    residual = g4.M - len(support)
    return tuple(
        Fraction((1 if index in support else 0) + (residual if index == favored else 0))
        for index in range(5)
    )


def refine_pure_triangle_star(
    anchors: tuple[Profile, ...],
    rows: list[dict[str, Any]],
    support: tuple[int, ...],
    classes: tuple[int, int, int],
    mode: str,
) -> tuple[tuple[Profile, ...], list[dict[str, Any]], list[dict[str, Any]]]:
    if len(set(classes)) != 3 or any(value not in support for value in classes):
        raise ValueError("pure triangle classes must be three distinct active classes")
    mutable = list(anchors)
    lookup = {profile: index for index, profile in enumerate(mutable)}
    triangle = tuple(lookup[pure_residual_profile(support, value)] for value in classes)
    triangle_profiles = [mutable[index] for index in triangle]

    def add(profile: Profile) -> int:
        existing = lookup.get(profile)
        if existing is not None:
            return existing
        index = len(mutable)
        mutable.append(profile)
        lookup[profile] = index
        return index

    face_triangles: list[tuple[int, int, int]] = []
    refinement_fields: dict[str, Any]
    if mode == "m3":
        denominator = 3
        grid = {}
        grid_rows = []
        for first in range(denominator + 1):
            for second in range(denominator - first + 1):
                third = denominator - first - second
                numerators = (first, second, third)
                profile = tuple(
                    sum(
                        (
                            Fraction(numerators[local], denominator)
                            * triangle_profiles[local][coordinate]
                            for local in range(3)
                        ),
                        Fraction(0),
                    )
                    for coordinate in range(5)
                )
                index = add(profile)
                grid[numerators] = index
                grid_rows.append(
                    {
                        "barycentric_numerators": list(numerators),
                        "anchor_index": index,
                    }
                )
        # Use (j,k) as the last two barycentric coordinates.  Six upward and
        # three downward unit lattice triangles partition 3*Delta_2.
        for second in range(3):
            for third in range(3 - second):
                face_triangles.append(
                    (
                        grid[(3 - second - third, second, third)],
                        grid[(2 - second - third, second + 1, third)],
                        grid[(2 - second - third, second, third + 1)],
                    )
                )
        for second in range(2):
            for third in range(2 - second):
                face_triangles.append(
                    (
                        grid[(2 - second - third, second + 1, third)],
                        grid[(1 - second - third, second + 1, third + 1)],
                        grid[(2 - second - third, second, third + 1)],
                    )
                )
        refinement_fields = {
            "refinement_kind": "barycentric_face_grid",
            "grid_denominator": denominator,
            "grid_nodes": grid_rows,
        }
    elif mode == "six-chambers":
        midpoint_pairs = ((0, 1), (0, 2), (1, 2))
        midpoints = tuple(
            add(
                tuple(
                    (
                        triangle_profiles[left][coordinate]
                        + triangle_profiles[right][coordinate]
                    )
                    / 2
                    for coordinate in range(5)
                )
            )
            for left, right in midpoint_pairs
        )
        midpoint_by_pair = {
            frozenset(pair): index for pair, index in zip(midpoint_pairs, midpoints)
        }
        centroid = add(
            tuple(
                sum(
                    (profile[coordinate] for profile in triangle_profiles),
                    Fraction(0),
                )
                / 3
                for coordinate in range(5)
            )
        )
        permutations = tuple(itertools.permutations(range(3)))
        face_triangles = [
            (
                triangle[first],
                midpoint_by_pair[frozenset((first, second))],
                centroid,
            )
            for first, second, _third in permutations
        ]
        refinement_fields = {
            "refinement_kind": "ordered_barycentric_chambers",
            "edge_midpoint_anchor_indices": list(midpoints),
            "centroid_anchor_index": centroid,
        }
    else:
        raise ValueError(f"unsupported pure triangle refinement mode {mode!r}")
    triangle_set = set(triangle)
    parents = [row for row in rows if triangle_set <= set(row["anchor_indices"])]
    if not parents:
        raise ValueError("selected pure triangle has no incident 4-simplex")
    parent_ids = {str(row["cell_id"]) for row in parents}
    terminal = [row for row in rows if str(row["cell_id"]) not in parent_ids]
    children = []
    ledger_children = []
    for parent in parents:
        complement = tuple(
            index for index in parent["anchor_indices"] if index not in triangle_set
        )
        if len(complement) != 2:
            raise RuntimeError("pure triangle parent does not have complementary edge")
        for local, chamber in enumerate(face_triangles):
            identifier = f"{parent['cell_id']}.ptg{3 if mode == 'm3' else 2}t{local:02d}"
            child = {
                "cell_id": identifier,
                "parent_cell_id": str(parent["cell_id"]),
                "anchor_indices": [*chamber, *complement],
            }
            children.append(child)
            ledger_children.append(dict(child))
    terminal.extend(children)
    terminal.sort(key=lambda row: str(row["cell_id"]))
    refinement = {
        "refinement_id": "pure_triangle_classes_" + "_".join(map(str, classes)),
        "triangle_classes": list(classes),
        "triangle_anchor_indices": list(triangle),
        **refinement_fields,
        "parents": [dict(row) for row in parents],
        "children": ledger_children,
    }
    return tuple(mutable), terminal, [refinement]


def build_stratum(
    mask: int,
    atlas_anchors: set[Profile],
    evaluator: cover.Evaluator | None = None,
    witnesses: list[cover.Witness] | None = None,
    pure_triangle_classes: tuple[int, int, int] | None = None,
    pure_triangle_mode: str = "m3",
) -> dict[str, Any]:
    support = tuple(index for index in range(5) if mask & (1 << index))
    dimension = len(support) - 1
    hull_vertices = integer_hull_vertices(support)
    if not hull_vertices:
        return {
            "support_mask": mask,
            "support": list(support),
            "dimension": dimension,
            "feasible": False,
            "reason": "no_positive_support_profile_meets_physical_weight_21",
            "integer_profile_count": 0,
            "profile_owner_rule": OWNER_RULE,
            "union_accounting_mode": UNION_MODE,
        }
    anchors = tuple(sorted({*hull_vertices, *atlas_anchors}))
    for profile in anchors:
        if any((profile[index] > 0) != (index in support) for index in range(5)):
            raise ValueError(f"stratum {support} contains an anchor of different support")
    heights = tuple(
        lifting_height(profile, support, rank) for rank, profile in enumerate(anchors)
    )
    proposals = propose_facets(anchors, heights, support)
    simplices = [
        simplex
        for simplex in proposals
        if simplex_determinant(simplex, anchors, support) != 0
    ]
    if not simplices:
        raise ValueError(f"support {support} has no nondegenerate proposed cells")
    slacks = [exact_lower_facet(simplex, anchors, heights, support) for simplex in simplices]
    unrefined_rows = [
        {
            "cell_id": f"s{mask:02x}r{cell_id:06d}",
            "anchor_indices": list(simplex),
        }
        for cell_id, simplex in enumerate(simplices)
    ]
    star_ledger = []
    if pure_triangle_classes is not None and set(pure_triangle_classes) <= set(support):
        anchors, rows, star_ledger = refine_pure_triangle_star(
            anchors,
            unrefined_rows,
            support,
            pure_triangle_classes,
            pure_triangle_mode,
        )
        simplices = [tuple(row["anchor_indices"]) for row in rows]
    else:
        rows = unrefined_rows
    geometry = verify_mesh(anchors, simplices, support, hull_vertices)
    covered_rows = []
    target = -40.0 - math.log2(profile_count(g4.GROUP_BITS, N))
    global_count = profile_count(g4.GROUP_BITS, N)
    if evaluator is not None and witnesses is not None:
        for row, simplex in zip(rows, simplices):
            values = np.asarray(
                [evaluator.evaluate(anchors[index])[1] for index in simplex],
                dtype=np.float64,
            )
            result = cover.optimize_values(
                values,
                witnesses,
                target,
                1e-5,
                MIXTURE_CANDIDATES_PER_VERTEX,
                1e-10,
            )
            cell = cover.Cell(str(row["cell_id"]), simplex, 0)
            count, dropped, ranges = cover.integer_profile_count_upper(
                cell, list(anchors), global_count
            )
            covered_rows.append(
                {
                    **row,
                    "mixture": cover.mixture_record(result, witnesses),
                    "maximum_vertex_log2": result["score"],
                    "target_log2": target,
                    "uniform_target_passed": result["passes_uniform_target"],
                    "integer_profile_count_upper": count,
                    "integer_profile_count_method": (
                        "minimum_dropped_coordinate_inclusive_integer_bounding_box_"
                        "capped_by_global_profile_count"
                    ),
                    "integer_profile_count_dropped_coordinate": dropped,
                    "inclusive_integer_coordinate_range_widths": list(ranges),
                }
            )
    return {
        "support_mask": mask,
        "support": list(support),
        "dimension": dimension,
        "feasible": True,
        "shift_rule": "a_i=1+b_i on active support; a_i=0 otherwise",
        "residual_mass": g4.M - len(support),
        "profile_owner_rule": OWNER_RULE,
        "union_accounting_mode": UNION_MODE,
        "complete_cover": evaluator is not None,
        "anchors": len(anchors),
        "atlas_anchors": len(atlas_anchors),
        "simplices": len(rows),
        "anchor_profiles": [[str(value) for value in profile] for profile in anchors],
        "exact_hull_vertices": [
            [str(value) for value in profile] for profile in hull_vertices
        ],
        "exact_hull_facets": [
            {
                "facet_id": label,
                "hull_vertex_indices": sorted(indices),
            }
            for label, indices in hull_boundary_sets(
                hull_vertices, support, hull_vertices
            ).items()
        ],
        "root_cell_ledger": rows,
        "unrefined_root_cell_ledger": unrefined_rows,
        "pure_triangle_star_refinement_ledger": star_ledger,
        "covered_cell_ledger": covered_rows,
        "failed_cell_ledger": [],
        "exact_geometry": geometry,
        "lifting": {
            "base": "sum of squared support-chart coordinates divided by M^2",
            "perturbation": "((rank+1)^5 mod 1000003)/(1000003*2^20)",
            "rank_order": "lexicographic exact profile order within support stratum",
            "heights_sha256": compact_hash([str(value) for value in heights]),
            "minimum_exact_nonvertex_lower_slack": str(min(slacks)),
            "proposed_facets": len(proposals),
            "exact_projected_degenerate_proposals_discarded": len(proposals) - len(simplices),
        },
        "anchor_profiles_sha256": compact_hash(
            [[str(value) for value in profile] for profile in anchors]
        ),
        "root_cell_ledger_sha256": compact_hash(rows),
    }


def build(
    paths: list[Path],
    masks: set[int] | None = None,
    *,
    geometry_only: bool = False,
    pure_triangle_classes: tuple[int, int, int] | None = None,
    pure_triangle_mode: str = "m3",
) -> dict[str, Any]:
    by_mask, sources, duplicate_rows = load_atlas_anchors(paths)
    evaluator = None
    witnesses = None
    if paths and not geometry_only:
        witnesses, hardened_sources, _anchors, _duplicates = cover.load_witnesses(paths)
        if {row["sha256"] for row in sources} != {
            row["sha256"] for row in hardened_sources
        }:
            raise RuntimeError("support mesh source hashing disagrees with mixture loader")
        evaluator = cover.Evaluator(witnesses)
    selected = set(range(1, 32)) if masks is None else masks
    strata = [
        (
            build_stratum(
                mask,
                by_mask[mask],
                evaluator,
                witnesses,
                pure_triangle_classes if mask == 31 else None,
                pure_triangle_mode,
            )
            if mask in selected
            else {
                "support_mask": mask,
                "support": [index for index in range(5) if mask & (1 << index)],
                "dimension": mask.bit_count() - 1,
                "feasible": False,
                "reason": "not_selected_in_prototype_run",
                "integer_profile_count": 0,
                "profile_owner_rule": OWNER_RULE,
                "union_accounting_mode": UNION_MODE,
            }
        )
        for mask in range(1, 32)
    ]
    return {
        "schema": SCHEMA,
        "status": "EXACT_CERTIFIED_G4_SUPPORT_STRATIFIED_REGULAR_ROOT_MESH",
        "group_bits": g4.GROUP_BITS,
        "profile_owner_rule": OWNER_RULE,
        "union_accounting_mode": UNION_MODE,
        "complete_support_stratified_cover": (
            evaluator is not None and selected == set(range(1, 32))
        ),
        "support_strata": strata,
        "sources": sources,
        "duplicate_atlas_anchor_rows": duplicate_rows,
        "support_strata_sha256": compact_hash(strata),
        "assumptions": [
            "Every support mask is a disjoint exact-support integer-profile stratum.",
            "Qhull proposes lifted facets; exact lower support and exact mesh audits decide output.",
            "Stratum-local anchor indices never refer across support masks.",
        ],
    }


def run_self_test() -> None:
    # Exercise dimensions 0 through 4 using hull-only strata, including the
    # clipped full-support stratum requested as the first prototype target.
    masks = {1 << 4, (1 << 0) | (1 << 4), 0b10101, 0b10111, 0b11111}
    rows = [build_stratum(mask, set()) for mask in sorted(masks)]
    dimensions = {row["dimension"] for row in rows if row["feasible"]}
    if dimensions != set(range(5)):
        raise SystemExit(f"support regular mesh self-test missed dimensions: {dimensions}")
    full = next(row for row in rows if row["support_mask"] == 31)
    if (
        not full["feasible"]
        or full["dimension"] != 4
        or len(full["exact_hull_vertices"]) != 15
        or not full["root_cell_ledger"]
    ):
        raise SystemExit("support regular mesh self-test: full support did not close")
    synthetic_witnesses = [
        cover.Witness(
            reference="synthetic",
            constant_log2=-1.0e9,
            charge=(0.0,) * 5,
            support=frozenset(range(5)),
            subtract_normalization=True,
        )
    ]
    synthetic_evaluator = cover.Evaluator(synthetic_witnesses)
    covered_full = build_stratum(31, set(), synthetic_evaluator, synthetic_witnesses)
    if (
        not covered_full["complete_cover"]
        or len(covered_full["covered_cell_ledger"])
        != len(covered_full["root_cell_ledger"])
        or any(not row["mixture"] for row in covered_full["covered_cell_ledger"])
    ):
        raise SystemExit("support regular mesh self-test: mixture integration failed")
    # A local face-grid replacement is intentionally rejected: its edge
    # nodes hang against cells incident to a triangle edge but not the full
    # triangle.  Such nodes must instead enter the global regular
    # retriangulation.  Keep this as a regression test for mesh soundness.
    try:
        build_stratum(31, set(), None, None, (1, 2, 3), "m3")
    except ValueError as error:
        if "unpaired facet" not in str(error):
            raise
    else:
        raise SystemExit("support regular mesh self-test: unsafe local star was accepted")
    print(f"self_test_feasible_strata={len(rows)}")
    print(f"self_test_full_support_simplices={full['simplices']}")
    print("status=PASS_G4_SUPPORT_STRATIFIED_REGULAR_MESH_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, action="append", default=[])
    parser.add_argument("--support-mask", type=lambda value: int(value, 0), action="append")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--geometry-only", action="store_true")
    parser.add_argument(
        "--pure-triangle-classes",
        default="",
        help="disabled: face-grid nodes require global regular retriangulation",
    )
    parser.add_argument(
        "--pure-triangle-mode",
        choices=("m3", "six-chambers"),
        default="m3",
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.output is None:
        raise SystemExit("support regular mesh: --output is required")
    masks = None if args.support_mask is None else set(args.support_mask)
    if masks is not None and (not masks or not masks <= set(range(1, 32))):
        raise SystemExit("support regular mesh: masks must lie in 1..31")
    if args.pure_triangle_classes.strip():
        raise SystemExit(
            "support regular mesh: local pure-triangle refinement is nonconforming; "
            "add face-grid nodes as global atlas anchors and retriangulate"
        )
    else:
        pure_classes = None
    report = build(
        args.atlas,
        masks,
        geometry_only=args.geometry_only,
        pure_triangle_classes=pure_classes,
        pure_triangle_mode=args.pure_triangle_mode,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    feasible = sum(row["feasible"] for row in report["support_strata"])
    simplices = sum(row.get("simplices", 0) for row in report["support_strata"])
    print(f"support_strata=31 feasible={feasible} simplices={simplices}")
    print("exact_support_mesh_audit=PASS")


if __name__ == "__main__":
    main()
