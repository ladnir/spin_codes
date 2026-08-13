#!/usr/bin/env python3
"""Exact/outward verifier for an enriched g=4 anchor-mesh certificate.

The Qhull/LP producer is discovery code.  A certifiable copy of its ledger must
replace binary64 mixture weights by exact rationals, attach SHA-256 digests to
all witness sources, and state the deterministic integer-profile owner rule.
This verifier reconstructs every 4-simplex over ``Fraction``, audits its
oriented facets and the exact clipped-hull volume, outward-hardens only used
witnesses, and performs one profile-count union bound (never one per cell).
Its preferred recursive form replaces a cell by the five cones from any exact
rational point with strictly positive barycentric coordinates; exact
barycenters are the symmetric special case.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import certify_packet_group_triangle_ledger as g2
from outward_log2 import Interval, LN2, ln_factorial, log2_int, self_check
from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import BLOCK_BITS, D, K, N, atom_count, profile_classes


GROUP_BITS = 4
DIMENSION = 4
M = atom_count(GROUP_BITS)
MINIMUM_PHYSICAL_WEIGHT = 21
FLAT_OWNER_RULE = "minimum_cell_id_among_closed_simplices"
STELLAR_OWNER_RULE = "minimum_leaf_owner_rank_among_closed_simplices"
SUPPORT_OWNER_RULE = "minimum_cell_id_within_exact_support_stratum"
SUPPORT_SCHEMA = "packet-group-g4-support-stratified-exact-regular-mesh-v1"
CELL_COUNT_METHOD = (
    "minimum_dropped_coordinate_inclusive_integer_bounding_box_"
    "capped_by_global_profile_count"
)


def parse_fraction(value: Any) -> Fraction:
    """Parse exact certificate data; binary64 values are deliberately refused."""

    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError("certificate rationals must not be JSON binary64 numbers")
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        try:
            return Fraction(value.strip())
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError(f"invalid exact rational {value!r}") from error
    if isinstance(value, dict) and set(value) >= {"numerator", "denominator"}:
        return Fraction(int(value["numerator"]), int(value["denominator"]))
    raise ValueError(f"unsupported exact rational encoding {value!r}")


def determinant(matrix: Iterable[Iterable[Fraction]]) -> Fraction:
    """Exact fraction-free determinant for the tiny fixed 4x4 matrices."""

    rows = [[Fraction(value) for value in row] for row in matrix]
    size = len(rows)
    if any(len(row) != size for row in rows):
        raise ValueError("determinant matrix must be square")
    if not size:
        return Fraction(1)
    sign = 1
    previous = Fraction(1)
    for column in range(size - 1):
        pivot = next((row for row in range(column, size) if rows[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            sign = -sign
        pivot_value = rows[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                rows[row][target] = (
                    rows[row][target] * pivot_value
                    - rows[row][column] * rows[column][target]
                ) / previous
            rows[row][column] = Fraction(0)
        previous = pivot_value
    return sign * rows[-1][-1]


def permutation_sign(values: tuple[int, ...]) -> int:
    inversions = sum(
        values[left] > values[right]
        for left in range(len(values))
        for right in range(left + 1, len(values))
    )
    return -1 if inversions % 2 else 1


def parse_anchor_profiles(
    raw: Any, label: str = "anchor mesh"
) -> tuple[tuple[Fraction, ...], ...]:
    if not isinstance(raw, list) or not raw:
        raise ValueError(f"{label} lacks anchor_profiles")
    anchors = []
    for index, profile in enumerate(raw):
        if not isinstance(profile, list) or len(profile) != GROUP_BITS + 1:
            raise ValueError(f"anchor {index} is not a five-class profile")
        point = tuple(parse_fraction(value) for value in profile)
        if any(value < 0 for value in point) or sum(point) != M:
            raise ValueError(f"invalid anchor profile {index}: {point}")
        if sum(weight * point[weight] for weight in range(5)) < MINIMUM_PHYSICAL_WEIGHT:
            raise ValueError(f"anchor {index} lies below the clipped physical-weight hull")
        anchors.append(point)
    if len(set(anchors)) != len(anchors):
        raise ValueError(f"{label} anchor_profiles contains duplicates")
    return tuple(anchors)


def load_anchors(ledger: dict[str, Any]) -> tuple[tuple[Fraction, ...], ...]:
    return parse_anchor_profiles(ledger.get("anchor_profiles"))


def cell_rows(ledger: dict[str, Any]) -> list[dict[str, Any]]:
    covered = ledger.get("covered_cell_ledger", [])
    failed = ledger.get("failed_cell_ledger", [])
    if not isinstance(covered, list) or not isinstance(failed, list):
        raise ValueError("cell ledgers must be arrays")
    if failed or not bool(ledger.get("complete_cover", False)):
        raise ValueError("anchor mesh producer did not assert a complete real cover")
    rows = list(covered)
    expected = int(ledger.get("simplices", -1))
    if expected != len(rows):
        raise ValueError(f"simplex ledger length mismatch: {len(rows)} != {expected}")
    identifiers = [int(row.get("cell_id", -1)) for row in rows]
    if len(set(identifiers)) != len(rows) or set(identifiers) != set(range(expected)):
        raise ValueError("cell_id values must be exactly 0..simplices-1")
    return sorted(rows, key=lambda row: int(row["cell_id"]))


def cell_key(row: dict[str, Any]) -> str:
    return str(row.get("cell_id"))


def stellar_leaf_rows(
    ledger: dict[str, Any], anchors: tuple[tuple[Fraction, ...], ...]
) -> tuple[list[dict[str, Any]], dict[str, Any]] | None:
    """Verify a recursive exact-interior-point subdivision ledger, if present."""

    roots = ledger.get("root_cell_ledger")
    subdivisions = ledger.get("stellar_subdivision_ledger")
    if roots is None and subdivisions is None:
        return None
    if not isinstance(roots, list) or not roots:
        raise ValueError("stellar certificate requires a nonempty root_cell_ledger")
    if not isinstance(subdivisions, list):
        raise ValueError("stellar_subdivision_ledger must be an array")
    leaves = ledger.get("covered_cell_ledger", [])
    if not isinstance(leaves, list) or not leaves:
        raise ValueError("stellar certificate requires covered leaf cells")

    nodes: dict[str, tuple[int, ...]] = {}
    for row in roots:
        key = cell_key(row)
        if not key or key in nodes:
            raise ValueError("root cell identifiers must be unique strings")
        nodes[key] = simplex_indices(row, len(anchors))
    root_geometry = verify_mesh(
        {"complete_cover": True, "simplices": len(roots)}, anchors, roots
    )

    children_by_parent = {}
    for row in subdivisions:
        parent = str(row.get("parent_cell_id"))
        if parent in children_by_parent:
            raise ValueError(f"stellar parent {parent!r} is subdivided twice")
        raw_children = row.get("children")
        if not isinstance(raw_children, list) or len(raw_children) != 5:
            raise ValueError("every stellar subdivision must list exactly five children")
        children_by_parent[parent] = (row, raw_children)

    visiting = set()
    terminal_nodes: dict[str, tuple[int, ...]] = {}
    subdivision_count = 0
    insertion_kinds = Counter()

    def descend(identifier: str, vertices: tuple[int, ...]) -> None:
        nonlocal subdivision_count
        if identifier in visiting:
            raise ValueError("stellar subdivision ledger contains a cycle")
        item = children_by_parent.get(identifier)
        if item is None:
            terminal_nodes[identifier] = vertices
            return
        visiting.add(identifier)
        row, raw_children = item
        raw_insertion = row.get("insertion_anchor_index", row.get("centroid_anchor_index", -1))
        insertion_index = int(raw_insertion)
        if not 0 <= insertion_index < len(anchors):
            raise ValueError(f"stellar parent {identifier!r} has invalid insertion index")
        expected_centroid = tuple(
            sum((anchors[index][coordinate] for index in vertices), Fraction(0)) / 5
            for coordinate in range(5)
        )
        barycentric = barycentric_coordinates(anchors[insertion_index], vertices, anchors)
        if any(coordinate <= 0 for coordinate in barycentric):
            raise ValueError(
                f"stellar parent {identifier!r} insertion is not strictly interior; "
                f"barycentric={tuple(map(str, barycentric))}"
            )
        reported = row.get("insertion_barycentric_coordinates")
        if reported is not None:
            parsed = tuple(parse_fraction(value) for value in reported)
            if parsed != barycentric:
                raise ValueError(
                    f"stellar parent {identifier!r} reported barycentric coordinates mismatch"
                )
        insertion_kind = "exact_barycenter" if anchors[insertion_index] == expected_centroid else "exact_strict_interior"
        expected_children = {
            tuple(sorted((insertion_index,) + vertices[:omitted] + vertices[omitted + 1 :]))
            for omitted in range(5)
        }
        actual_children = {}
        for child in raw_children:
            if not isinstance(child, dict):
                raise ValueError("stellar child records must be objects")
            child_id = cell_key(child)
            child_vertices = simplex_indices(child, len(anchors))
            child_key = tuple(sorted(child_vertices))
            if child_id in nodes or child_id in actual_children:
                raise ValueError(f"duplicate stellar cell identifier {child_id!r}")
            actual_children[child_id] = child_vertices
            nodes[child_id] = child_vertices
        if {tuple(sorted(value)) for value in actual_children.values()} != expected_children:
            raise ValueError(f"stellar parent {identifier!r} children are not its five insertion cones")
        insertion_kinds[insertion_kind] += 1
        subdivision_count += 1
        for child_id, child_vertices in actual_children.items():
            descend(child_id, child_vertices)
        visiting.remove(identifier)

    for identifier, vertices in list(nodes.items()):
        if identifier in {str(row.get("parent_cell_id")) for row in subdivisions} or identifier in {
            cell_key(row) for row in roots
        }:
            # Only original roots start recursion; children are reached by their parent.
            if any(cell_key(root) == identifier for root in roots):
                descend(identifier, vertices)

    unused_parents = set(children_by_parent) - set(nodes)
    if unused_parents:
        raise ValueError(f"stellar subdivision references unreachable parents: {sorted(unused_parents)}")
    leaf_map = {cell_key(row): row for row in leaves}
    if len(leaf_map) != len(leaves) or set(leaf_map) != set(terminal_nodes):
        raise ValueError("covered leaf ledger is not exactly the recursive terminal-node set")
    for identifier, row in leaf_map.items():
        if tuple(sorted(simplex_indices(row, len(anchors)))) != tuple(
            sorted(terminal_nodes[identifier])
        ):
            raise ValueError(f"covered leaf {identifier!r} vertices mismatch recursive geometry")
    owner_ranks = [int(row.get("owner_rank", -1)) for row in leaves]
    if set(owner_ranks) != set(range(len(leaves))):
        raise ValueError("stellar leaf owner_rank values must be exactly 0..leaf_count-1")
    if ledger.get("failed_cell_ledger"):
        raise ValueError("stellar certificate contains failed leaf cells")
    return list(leaf_map.values()), {
        "mode": "recursive_exact_centroid_stellar_subdivision",
        "root_geometry": root_geometry,
        "root_cells": len(roots),
        "stellar_subdivisions": subdivision_count,
        "insertion_kinds": dict(sorted(insertion_kinds.items())),
        "covered_leaves": len(leaves),
    }


def simplex_indices(row: dict[str, Any], anchor_count: int) -> tuple[int, ...]:
    raw = row.get("anchor_indices")
    if not isinstance(raw, list) or len(raw) != DIMENSION + 1:
        raise ValueError("every g=4 cell must contain five anchor indices")
    indices = tuple(int(value) for value in raw)
    if len(set(indices)) != len(indices) or min(indices) < 0 or max(indices) >= anchor_count:
        raise ValueError(f"invalid simplex anchor indices {indices}")
    return indices


def simplex_determinant(
    indices: tuple[int, ...], anchors: tuple[tuple[Fraction, ...], ...]
) -> Fraction:
    points = [anchors[index][1:] for index in indices]
    return determinant(
        [[points[row][column] - points[0][column] for column in range(4)] for row in range(1, 5)]
    )


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def simplex_integer_profile_count_upper(
    indices: tuple[int, ...],
    anchors: tuple[tuple[Fraction, ...], ...],
) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    """Cheap rigorous integer-profile count via the best coordinate box.

    Any four coordinates uniquely determine the fifth because their sum is
    M.  For each omitted coordinate, bound the other four by their exact
    simplex coordinate ranges and take the smallest product.
    """

    best = None
    for omitted in range(5):
        ranges = []
        count = 1
        for coordinate in range(5):
            if coordinate == omitted:
                continue
            values = [anchors[index][coordinate] for index in indices]
            low = ceil_fraction(min(values))
            high = floor_fraction(max(values))
            width = max(0, high - low + 1)
            ranges.append((low, high))
            count *= width
        count = min(count, profile_count(GROUP_BITS, N))
        candidate = (count, omitted, tuple(ranges))
        if best is None or candidate < best:
            best = candidate
    assert best is not None
    return best


def support_hull_vertices(support: tuple[int, ...]) -> tuple[tuple[Fraction, ...], ...]:
    """Exact positive-coordinate support hull clipped at physical weight 21."""

    residual = M - len(support)
    base = [Fraction(1) if index in support else Fraction(0) for index in range(5)]
    pure = {}
    for favored in support:
        profile = list(base)
        profile[favored] += residual
        pure[favored] = tuple(profile)
    physical = {
        favored: sum(index * value for index, value in enumerate(profile))
        for favored, profile in pure.items()
    }
    vertices = {
        profile for favored, profile in pure.items() if physical[favored] >= MINIMUM_PHYSICAL_WEIGHT
    }
    for invalid in support:
        if physical[invalid] >= MINIMUM_PHYSICAL_WEIGHT:
            continue
        for valid in support:
            if physical[valid] < MINIMUM_PHYSICAL_WEIGHT:
                continue
            denominator = physical[valid] - physical[invalid]
            toward_valid = Fraction(MINIMUM_PHYSICAL_WEIGHT - physical[invalid], denominator)
            profile = tuple(
                (1 - toward_valid) * left + toward_valid * right
                for left, right in zip(pure[invalid], pure[valid])
            )
            vertices.add(profile)
    return tuple(sorted(vertices))


def integer_hull_candidates(support: tuple[int, ...]) -> tuple[tuple[Fraction, ...], ...]:
    """Finite exact candidate set whose convex hull is the stratum integer hull.

    For strata containing zero, only the shifted corner band
    ``T <= w.b < T+max(S)`` can contain a non-pure vertex.  Above it, an
    occupied positive coordinate permits opposite unit transfers to/from
    class zero, expressing the point as a midpoint.  The b0=0 face has only
    the nonzero shifted-simplex pure vertices as possible vertices.
    """

    residual = M - len(support)
    base_weight = sum(support)
    threshold = MINIMUM_PHYSICAL_WEIGHT - base_weight
    candidates = set()
    for favored in support:
        shifted = [0] * 5
        shifted[favored] = residual
        profile = tuple(
            Fraction(shifted[index] + (1 if index in support else 0))
            for index in range(5)
        )
        if sum(index * profile[index] for index in range(5)) >= MINIMUM_PHYSICAL_WEIGHT:
            candidates.add(profile)
    if 0 not in support or threshold <= 0:
        return tuple(sorted(candidates))
    positive = tuple(index for index in support if index)
    if not positive:
        return ()
    maximum_weight = max(positive)

    def descend(position: int, used: int, physical: int, values: list[int]) -> None:
        if position == len(positive):
            if threshold <= physical < threshold + maximum_weight:
                shifted = [0] * 5
                shifted[0] = residual - used
                for coordinate, value in zip(positive, values):
                    shifted[coordinate] = value
                candidates.add(
                    tuple(
                        Fraction(shifted[index] + (1 if index in support else 0))
                        for index in range(5)
                    )
                )
            return
        coordinate = positive[position]
        maximum = min(
            residual - used,
            max(0, (threshold + maximum_weight - 1 - physical) // coordinate),
        )
        for value in range(maximum + 1):
            values.append(value)
            descend(position + 1, used + value, physical + coordinate * value, values)
            values.pop()

    descend(0, 0, 0, [])
    return tuple(sorted(candidates))


def support_affine_dimension(
    indices: Iterable[int],
    profiles: tuple[tuple[Fraction, ...], ...],
    support: tuple[int, ...],
) -> int:
    selected = tuple(indices)
    if not selected:
        return -1
    dimension = len(support) - 1
    if len(selected) == 1:
        return 0
    rows = [
        [profiles[index][coordinate] - profiles[selected[0]][coordinate] for coordinate in support[1:]]
        for index in selected[1:]
    ]
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
            rows[row] = [value - factor * other for value, other in zip(rows[row], rows[rank])]
        rank += 1
    return rank


def exact_integer_hull_facets(
    vertices: tuple[tuple[Fraction, ...], ...], support: tuple[int, ...]
) -> tuple[frozenset[int], ...]:
    """Enumerate every maximal exact supporting facet of a tiny d<=4 hull."""

    import itertools

    dimension = len(support) - 1
    if dimension == 0:
        return ()
    if dimension == 1:
        order = sorted(range(len(vertices)), key=lambda index: vertices[index][support[1]])
        return (frozenset((order[0],)), frozenset((order[-1],)))
    facets = set()
    for seed in itertools.combinations(range(len(vertices)), dimension):
        if support_affine_dimension(seed, vertices, support) != dimension - 1:
            continue
        base = vertices[seed[0]]
        signs = []
        coplanar = set(seed)
        for index, profile in enumerate(vertices):
            if index in coplanar:
                continue
            matrix = [
                [vertices[row][coordinate] - base[coordinate] for coordinate in support[1:]]
                for row in seed[1:]
            ]
            matrix.append([profile[coordinate] - base[coordinate] for coordinate in support[1:]])
            value = determinant(matrix)
            if value:
                signs.append(1 if value > 0 else -1)
            else:
                coplanar.add(index)
        if signs and (all(sign > 0 for sign in signs) or all(sign < 0 for sign in signs)):
            facets.add(frozenset(coplanar))
    return tuple(sorted(facets, key=lambda facet: tuple(sorted(facet))))


def point_is_inside_exact_facets(
    point: tuple[Fraction, ...],
    vertices: tuple[tuple[Fraction, ...], ...],
    support: tuple[int, ...],
    facets: tuple[frozenset[int], ...],
) -> bool:
    dimension = len(support) - 1
    if dimension == 0:
        return point == vertices[0]
    if dimension == 1:
        coordinate = support[1]
        return min(row[coordinate] for row in vertices) <= point[coordinate] <= max(
            row[coordinate] for row in vertices
        )
    for facet in facets:
        seed = next(
            candidate
            for candidate in __import__("itertools").combinations(sorted(facet), dimension)
            if support_affine_dimension(candidate, vertices, support) == dimension - 1
        )
        base = vertices[seed[0]]

        def side(profile: tuple[Fraction, ...]) -> Fraction:
            matrix = [
                [vertices[row][coordinate] - base[coordinate] for coordinate in support[1:]]
                for row in seed[1:]
            ]
            matrix.append([profile[coordinate] - base[coordinate] for coordinate in support[1:]])
            return determinant(matrix)

        reference = next(side(profile) for profile in vertices if side(profile))
        value = side(point)
        if value and (value > 0) != (reference > 0):
            return False
    return True


def integer_hull_vertices_from_certificate(
    support: tuple[int, ...], raw_vertices: Any, raw_facets: Any
) -> tuple[tuple[tuple[Fraction, ...], ...], tuple[frozenset[int], ...]]:
    """Independently establish that claimed vertices/facets equal the integer hull."""

    if not isinstance(raw_vertices, list):
        raise ValueError(f"support {support} lacks exact_hull_vertices")
    vertices = tuple(
        tuple(parse_fraction(value) for value in profile) for profile in raw_vertices
    )
    if any(len(profile) != 5 for profile in vertices) or len(set(vertices)) != len(vertices):
        raise ValueError(f"support {support} has malformed exact_hull_vertices")
    candidates = integer_hull_candidates(support)
    if not candidates:
        if vertices:
            raise ValueError(f"infeasible support {support} claims integer hull vertices")
        return (), ()
    if not set(vertices) <= set(candidates):
        raise ValueError(f"support {support} claims a vertex outside the exact candidate lemma")
    facets = exact_integer_hull_facets(vertices, support)
    if not facets and len(support) > 1:
        raise ValueError(f"support {support} claimed vertices do not span a hull")
    if len(support) > 1:
        for index in range(len(vertices)):
            incident = [set(facet) for facet in facets if index in facet]
            if not incident or set.intersection(*incident) != {index}:
                raise ValueError(f"support {support} exact_hull_vertices contains a nonvertex")
    if any(not point_is_inside_exact_facets(point, vertices, support, facets) for point in candidates):
        raise ValueError(f"support {support} integer-hull candidate lies outside claimed facets")
    if not isinstance(raw_facets, list):
        raise ValueError(f"support {support} lacks exact_hull_facets")
    reported = []
    identifiers = set()
    for row in raw_facets:
        if not isinstance(row, dict):
            raise ValueError("exact_hull_facets rows must be objects")
        identifier = str(row.get("facet_id", ""))
        if not identifier or identifier in identifiers:
            raise ValueError(f"support {support} has invalid facet identifiers")
        identifiers.add(identifier)
        indices = frozenset(int(value) for value in row.get("hull_vertex_indices", []))
        if not indices or min(indices) < 0 or max(indices) >= len(vertices):
            raise ValueError(f"support {support} has invalid hull facet indices")
        reported.append(indices)
    if set(reported) != set(facets) or len(reported) != len(facets):
        raise ValueError(f"support {support} exact_hull_facets mismatch")
    return vertices, facets


def integer_hull_pulling_cells(
    vertices: tuple[tuple[Fraction, ...], ...],
    support: tuple[int, ...],
    facets: tuple[frozenset[int], ...],
) -> tuple[tuple[int, ...], ...]:
    dimension = len(support) - 1
    if dimension == 0:
        return ((0,),)

    def triangulate(face: frozenset[int], local_dimension: int) -> list[tuple[int, ...]]:
        if local_dimension == 0:
            return [(min(face),)]
        pivot = min(face)
        opposite = {
            candidate
            for boundary in facets
            if (candidate := face & boundary)
            and pivot not in candidate
            and support_affine_dimension(candidate, vertices, support) == local_dimension - 1
        }
        if not opposite:
            raise ValueError(f"support {support} exact facets do not define a pulling triangulation")
        result = []
        for facet in sorted(opposite, key=lambda row: tuple(sorted(row))):
            for simplex in triangulate(facet, local_dimension - 1):
                result.append((pivot, *simplex))
        return result

    cells = triangulate(frozenset(range(len(vertices))), dimension)
    if len({tuple(sorted(cell)) for cell in cells}) != len(cells):
        raise ValueError(f"support {support} pulling triangulation contains duplicates")
    return tuple(cells)


def support_simplex_determinant(
    indices: tuple[int, ...],
    anchors: tuple[tuple[Fraction, ...], ...],
    support: tuple[int, ...],
) -> Fraction:
    dimension = len(support) - 1
    points = [anchors[index] for index in indices]
    return determinant(
        [
            [points[row][coordinate] - points[0][coordinate] for coordinate in support[1:]]
            for row in range(1, dimension + 1)
        ]
    )


def support_hull_determinant_volume(support: tuple[int, ...]) -> Fraction:
    dimension = len(support) - 1
    vertices = support_hull_vertices(support)
    if dimension == 0:
        return Fraction(1) if vertices else Fraction(0)
    residual = M - len(support)
    threshold = MINIMUM_PHYSICAL_WEIGHT - sum(support)
    below = Fraction(0)
    for weight in support:
        numerator = max(0, threshold - residual * weight) ** dimension
        denominator = math.prod(other - weight for other in support if other != weight)
        below += Fraction(numerator, denominator)
    return Fraction(residual**dimension) - below


def stratum_cell_indices(
    row: dict[str, Any], anchor_count: int, dimension: int
) -> tuple[int, ...]:
    raw = row.get("anchor_indices")
    if not isinstance(raw, list) or len(raw) != dimension + 1:
        raise ValueError(
            f"support dimension {dimension} cell must contain {dimension + 1} anchors"
        )
    indices = tuple(int(value) for value in raw)
    if len(set(indices)) != len(indices) or min(indices) < 0 or max(indices) >= anchor_count:
        raise ValueError(f"invalid support-stratum cell anchor indices {indices}")
    return indices


def verify_pure_triangle_star_refinements(
    stratum: dict[str, Any],
    anchors: tuple[tuple[Fraction, ...], ...],
    terminal_rows: list[dict[str, Any]],
    support: tuple[int, ...],
) -> dict[str, Any] | None:
    """Check a claimed local face-grid replacement before the global audit.

    This is only a necessary ledger-consistency check.  Face grids with edge
    nodes generally create hanging facets outside the full-triangle star; the
    mandatory ``verify_support_stratum_mesh`` call on the resulting terminal
    mesh independently rejects those unless a global conforming
    retriangulation has also closed every adjacent face.
    """

    ledger = stratum.get("pure_triangle_star_refinement_ledger")
    if ledger is None:
        return None
    if not isinstance(ledger, list) or not ledger:
        raise ValueError("pure_triangle_star_refinement_ledger must be a nonempty array")
    dimension = len(support) - 1
    if dimension != 4:
        raise ValueError("pure-triangle complementary-edge refinement currently requires dimension 4")
    base_rows = stratum.get("unrefined_root_cell_ledger")
    if not isinstance(base_rows, list) or not base_rows:
        raise ValueError("pure-triangle refinement requires unrefined_root_cell_ledger")
    current: dict[str, dict[str, Any]] = {}
    for row in base_rows:
        identifier = cell_key(row)
        stratum_cell_indices(row, len(anchors), dimension)
        if identifier in current:
            raise ValueError("unrefined root cell identifiers are not unique")
        current[identifier] = row
    total_parents = 0
    total_children = 0
    for refinement in ledger:
        if not isinstance(refinement, dict):
            raise ValueError("pure-triangle refinement rows must be objects")
        triangle = tuple(int(value) for value in refinement.get("triangle_anchor_indices", []))
        if len(triangle) != 3 or len(set(triangle)) != 3 or min(triangle) < 0 or max(triangle) >= len(anchors):
            raise ValueError("pure-triangle refinement has invalid triangle_anchor_indices")
        favored = []
        for index in triangle:
            profile = anchors[index]
            excess = [coordinate for coordinate in support if profile[coordinate] > 1]
            if len(excess) != 1 or any(
                profile[coordinate] != 1 for coordinate in support if coordinate != excess[0]
            ):
                raise ValueError("selected pure triangle contains a non-pure shifted-simplex vertex")
            favored.append(excess[0])
        if len(set(favored)) != 3:
            raise ValueError("selected pure triangle repeats a favored class")
        refinement_kind = str(refinement.get("refinement_kind", "ordered_barycentric_chambers"))
        face_cells: list[tuple[tuple[int, int, int], tuple[int, ...]]] = []
        if refinement_kind in {"ordered_barycentric_chambers", "six_ordered_barycentric_chambers"}:
            midpoint_indices = tuple(
                int(value) for value in refinement.get("edge_midpoint_anchor_indices", [])
            )
            if len(midpoint_indices) != 3 or len(set(midpoint_indices)) != 3:
                raise ValueError("pure-triangle refinement must list three distinct edge midpoints")
            pairs = ((0, 1), (0, 2), (1, 2))
            midpoint_by_pair = {}
            for pair, midpoint_index in zip(pairs, midpoint_indices):
                if not 0 <= midpoint_index < len(anchors):
                    raise ValueError("pure-triangle midpoint anchor index is invalid")
                expected = tuple(
                    (anchors[triangle[pair[0]]][coordinate] + anchors[triangle[pair[1]]][coordinate]) / 2
                    for coordinate in range(5)
                )
                if anchors[midpoint_index] != expected:
                    raise ValueError(f"pure-triangle edge midpoint {pair} is not exact")
                midpoint_by_pair[frozenset(pair)] = midpoint_index
            centroid_index = int(refinement.get("centroid_anchor_index", -1))
            if not 0 <= centroid_index < len(anchors):
                raise ValueError("pure-triangle centroid anchor index is invalid")
            centroid = tuple(
                sum((anchors[index][coordinate] for index in triangle), Fraction(0)) / 3
                for coordinate in range(5)
            )
            if anchors[centroid_index] != centroid:
                raise ValueError("pure-triangle centroid is not exact")
            for permutation in itertools.permutations(range(3)):
                face_cells.append(
                    (
                        permutation,
                        (
                            triangle[permutation[0]],
                            midpoint_by_pair[frozenset((permutation[0], permutation[1]))],
                            centroid_index,
                        ),
                    )
                )
        elif refinement_kind == "barycentric_face_grid":
            if int(refinement.get("grid_denominator", -1)) != 3:
                raise ValueError("pure-triangle grid_denominator must equal 3")
            raw_nodes = refinement.get("grid_nodes")
            if not isinstance(raw_nodes, list) or len(raw_nodes) != 10:
                raise ValueError("m=3 pure-triangle grid must list exactly ten nodes")
            nodes: dict[tuple[int, int, int], int] = {}
            for row in raw_nodes:
                numerators = tuple(int(value) for value in row.get("barycentric_numerators", []))
                anchor_index = int(row.get("anchor_index", -1))
                if (
                    len(numerators) != 3
                    or min(numerators) < 0
                    or sum(numerators) != 3
                    or numerators in nodes
                    or not 0 <= anchor_index < len(anchors)
                ):
                    raise ValueError("m=3 pure-triangle grid has an invalid node")
                expected = tuple(
                    sum(
                        (Fraction(numerators[position], 3) * anchors[triangle[position]][coordinate]
                         for position in range(3)),
                        Fraction(0),
                    )
                    for coordinate in range(5)
                )
                if anchors[anchor_index] != expected:
                    raise ValueError("m=3 pure-triangle grid node is not its exact barycentric point")
                nodes[numerators] = anchor_index
            expected_nodes = {
                (left, middle, 3 - left - middle)
                for left in range(4)
                for middle in range(4 - left)
            }
            if set(nodes) != expected_nodes:
                raise ValueError("m=3 pure-triangle grid does not contain all i+j+k=3 nodes")
            for numerators in itertools.combinations(sorted(nodes), 3):
                if all(
                    sum(abs(left[position] - right[position]) for position in range(3)) == 2
                    for left, right in itertools.combinations(numerators, 2)
                ):
                    face_cells.append((numerators[0], tuple(nodes[item] for item in numerators)))
            if len(face_cells) != 9:
                raise ValueError("m=3 pure-triangle grid does not induce exactly nine elementary triangles")
        else:
            raise ValueError(f"unsupported pure-triangle refinement_kind {refinement_kind!r}")

        incident = {
            identifier: row
            for identifier, row in current.items()
            if set(triangle) <= set(stratum_cell_indices(row, len(anchors), dimension))
        }
        if not incident:
            raise ValueError("pure-triangle selected triangle has an empty current star")
        raw_parents = refinement.get("parents")
        if not isinstance(raw_parents, list):
            raise ValueError("pure-triangle refinement lacks parents")
        parents = {cell_key(row): row for row in raw_parents}
        if len(parents) != len(raw_parents) or set(parents) != set(incident):
            raise ValueError("pure-triangle parents are not all and only the complete incident star")
        for identifier, row in parents.items():
            if tuple(sorted(stratum_cell_indices(row, len(anchors), dimension))) != tuple(
                sorted(stratum_cell_indices(incident[identifier], len(anchors), dimension))
            ):
                raise ValueError("pure-triangle reported parent geometry mismatch")

        raw_children = refinement.get("children")
        if not isinstance(raw_children, list):
            raise ValueError("pure-triangle refinement lacks children")
        children_by_parent: dict[str, list[dict[str, Any]]] = {identifier: [] for identifier in parents}
        for child in raw_children:
            parent_id = str(child.get("parent_cell_id", ""))
            if parent_id not in children_by_parent:
                # Canonical IDs also encode the parent and are accepted when
                # the redundant parent_cell_id field is omitted.
                matches = [identifier for identifier in parents if cell_key(child).startswith(identifier + ".pt")]
                if len(matches) != 1:
                    raise ValueError("pure-triangle child has no unique parent")
                parent_id = matches[0]
            children_by_parent[parent_id].append(child)
        new_rows = {}
        for parent_id, parent in parents.items():
            parent_indices = stratum_cell_indices(parent, len(anchors), dimension)
            complement = tuple(index for index in parent_indices if index not in triangle)
            if len(complement) != 2:
                raise ValueError("pure-triangle parent does not have an exact complementary edge")
            expected = {
                tuple(sorted((*face_cell, *complement))): label
                for label, face_cell in face_cells
            }
            children = children_by_parent[parent_id]
            if len(children) != len(face_cells):
                raise ValueError("pure-triangle parent child count does not match its face partition")
            actual = {}
            child_volume = Fraction(0)
            for child in children:
                indices = stratum_cell_indices(child, len(anchors), dimension)
                key = tuple(sorted(indices))
                if key in actual:
                    raise ValueError("pure-triangle children contain duplicate geometry")
                actual[key] = child
                child_det = support_simplex_determinant(indices, anchors, support)
                if not child_det:
                    raise ValueError("pure-triangle refinement produced a degenerate child")
                child_volume += abs(child_det)
            if set(actual) != set(expected):
                raise ValueError("pure-triangle children are not the exact face partition joined to the edge")
            parent_volume = abs(support_simplex_determinant(parent_indices, anchors, support))
            if child_volume != parent_volume:
                raise ValueError("pure-triangle child volumes do not sum exactly to parent volume")
            for key, child in actual.items():
                if not cell_key(child).startswith(parent_id + "."):
                    raise ValueError("pure-triangle child ID is not namespaced by its parent")
                if refinement_kind in {"ordered_barycentric_chambers", "six_ordered_barycentric_chambers"}:
                    expected_suffix = ".pts" + "".join(map(str, expected[key]))
                    if not cell_key(child).endswith(expected_suffix):
                        raise ValueError("pure-triangle child ID does not match chamber permutation")
                elif not any(
                    cell_key(child).endswith(f".ptg3t{ordinal:02d}")
                    for ordinal in range(9)
                ):
                    raise ValueError("m=3 pure-triangle child ID has no canonical triangle ordinal")
                if cell_key(child) in current or cell_key(child) in new_rows:
                    raise ValueError("pure-triangle child identifier is not unique")
                new_rows[cell_key(child)] = child
            if refinement_kind == "barycentric_face_grid" and {
                cell_key(child).rsplit(".ptg3t", 1)[-1] for child in actual.values()
            } != {f"{ordinal:02d}" for ordinal in range(9)}:
                raise ValueError("m=3 pure-triangle child IDs do not use ordinals 00..08 exactly")
        for parent_id in parents:
            del current[parent_id]
        current.update(new_rows)
        total_parents += len(parents)
        total_children += len(new_rows)

    terminal = {cell_key(row): row for row in terminal_rows}
    if len(terminal) != len(terminal_rows) or set(terminal) != set(current):
        raise ValueError("pure-triangle refinements do not reconstruct the exact terminal mesh")
    for identifier, row in terminal.items():
        if tuple(sorted(stratum_cell_indices(row, len(anchors), dimension))) != tuple(
            sorted(stratum_cell_indices(current[identifier], len(anchors), dimension))
        ):
            raise ValueError("pure-triangle terminal cell geometry mismatch")
    return {
        "refinements": len(ledger),
        "replaced_parents": total_parents,
        "created_children": total_children,
        "terminal_cells": len(current),
    }


def verify_support_stratum_mesh(
    support: tuple[int, ...],
    anchors: tuple[tuple[Fraction, ...], ...],
    rows: list[dict[str, Any]],
    hull_vertices: tuple[tuple[Fraction, ...], ...] | None = None,
    hull_facets: tuple[frozenset[int], ...] | None = None,
) -> dict[str, Any]:
    """Exact dimension-0..4 simplicial audit for one positive-support hull."""

    dimension = len(support) - 1
    expected_vertices = support_hull_vertices(support) if hull_vertices is None else hull_vertices
    if hull_facets is None:
        hull_facets = exact_integer_hull_facets(expected_vertices, support)
    pulling_cells = integer_hull_pulling_cells(expected_vertices, support, hull_facets) if expected_vertices else ()
    expected_volume = sum(
        (abs(support_simplex_determinant(cell, expected_vertices, support)) for cell in pulling_cells),
        Fraction(0),
    )
    if not expected_vertices or expected_volume == 0:
        if anchors or rows:
            raise ValueError(f"infeasible support {support} contains geometry")
        return {"support": list(support), "dimension": dimension, "feasible": False}
    expected_set = set(expected_vertices)
    for profile in anchors:
        if any(profile[index] <= 0 for index in support) or any(
            profile[index] != 0 for index in range(5) if index not in support
        ):
            raise ValueError(f"anchor {profile} does not have exact support {support}")
        if sum(index * profile[index] for index in support) < MINIMUM_PHYSICAL_WEIGHT:
            raise ValueError(f"support anchor {profile} lies below physical weight 21")
    if not expected_set <= set(anchors):
        raise ValueError(f"support {support} anchor list omits exact hull vertices")
    if dimension == 0:
        if len(rows) != 1 or len(anchors) != 1:
            raise ValueError(f"zero-dimensional support {support} needs one point cell")
        stratum_cell_indices(rows[0], len(anchors), 0)
        return {
            "support": list(support),
            "dimension": 0,
            "feasible": True,
            "cells": 1,
            "determinant_volume": "1",
        }

    facets: Counter[tuple[int, ...]] = Counter()
    orientations: Counter[tuple[int, ...]] = Counter()
    seen = set()
    volume = Fraction(0)
    for row in rows:
        indices = stratum_cell_indices(row, len(anchors), dimension)
        key = tuple(sorted(indices))
        if key in seen:
            raise ValueError(f"duplicate support-stratum simplex {key}")
        seen.add(key)
        points = [anchors[index] for index in indices]
        det = support_simplex_determinant(indices, anchors, support)
        if not det:
            raise ValueError(f"degenerate support-stratum simplex {key}")
        oriented = indices if det > 0 else (indices[1], indices[0], *indices[2:])
        volume += abs(det)
        for omitted in range(dimension + 1):
            order = oriented[:omitted] + oriented[omitted + 1 :]
            facet = tuple(sorted(order))
            facets[facet] += 1
            orientations[facet] += ((-1) ** omitted) * permutation_sign(order)

    boundary = Counter()
    interior = 0
    for facet, multiplicity in facets.items():
        if multiplicity == 2:
            if orientations[facet] != 0:
                raise ValueError(f"support {support} interior facet orientation mismatch")
            interior += 1
            continue
        if multiplicity != 1:
            raise ValueError(f"support {support} facet multiplicity is {multiplicity}")
        profiles = [anchors[index] for index in facet]
        labels = []
        for ordinal, hull_facet in enumerate(hull_facets):
            seed = next(
                candidate
                for candidate in __import__("itertools").combinations(sorted(hull_facet), dimension)
                if support_affine_dimension(candidate, expected_vertices, support) == dimension - 1
            )
            base = expected_vertices[seed[0]]
            if all(
                determinant(
                    [
                        [expected_vertices[row][coordinate] - base[coordinate] for coordinate in support[1:]]
                        for row in seed[1:]
                    ]
                    + [[profile[coordinate] - base[coordinate] for coordinate in support[1:]]]
                ) == 0
                for profile in profiles
            ):
                labels.append(f"integer_hull_facet_{ordinal:03d}")
        if len(labels) != 1:
            raise ValueError(f"support {support} has an unpaired nonboundary facet {facet}")
        boundary[labels[0]] += 1
    if volume != expected_volume:
        raise ValueError(
            f"support {support} exact volume mismatch: {volume} != {expected_volume}"
        )
    return {
        "support": list(support),
        "dimension": dimension,
        "feasible": True,
        "cells": len(rows),
        "interior_facets": interior,
        "boundary_facets": dict(sorted(boundary.items())),
        "determinant_volume": str(volume),
        "expected_determinant_volume": str(expected_volume),
        "integer_hull_facets": len(hull_facets),
        "exact_hull_vertices": [[str(value) for value in profile] for profile in expected_vertices],
    }


def support_stratum_contexts(
    ledger: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]] | None:
    """Load and exactly audit all 31 disjoint positive-support strata."""

    raw_strata = ledger.get("support_strata")
    if raw_strata is None:
        return None
    if ledger.get("schema") != SUPPORT_SCHEMA:
        raise ValueError(f"support-stratified ledger schema must be {SUPPORT_SCHEMA!r}")
    if not bool(ledger.get("complete_support_stratified_cover", False)):
        raise ValueError("producer did not assert a complete support-stratified cover")
    if not isinstance(raw_strata, list) or len(raw_strata) != 31:
        raise ValueError("support_strata must contain exactly the 31 nonempty supports")

    by_mask: dict[int, dict[str, Any]] = {}
    for stratum in raw_strata:
        if not isinstance(stratum, dict):
            raise ValueError("support-stratum records must be objects")
        mask = int(stratum.get("support_mask", -1))
        if mask in by_mask or not 1 <= mask < 32:
            raise ValueError(f"invalid or duplicate support mask {mask}")
        by_mask[mask] = stratum
    if set(by_mask) != set(range(1, 32)):
        raise ValueError("support masks must be exactly 1..31")

    contexts: list[dict[str, Any]] = []
    geometry_rows = []
    for mask in range(1, 32):
        stratum = by_mask[mask]
        support = tuple(index for index in range(5) if mask & (1 << index))
        if stratum.get("support") != list(support):
            raise ValueError(f"support mask {mask} has a mismatched support list")
        dimension = len(support) - 1
        if int(stratum.get("dimension", -1)) != dimension:
            raise ValueError(f"support mask {mask} has a mismatched dimension")
        if stratum.get("profile_owner_rule") != SUPPORT_OWNER_RULE:
            raise ValueError(f"support mask {mask} has a mismatched profile_owner_rule")
        if stratum.get("union_accounting_mode") != "cell_local_profile_count_upper":
            raise ValueError(f"support mask {mask} has a mismatched union_accounting_mode")
        candidates = integer_hull_candidates(support)
        feasible = bool(candidates)
        if bool(stratum.get("feasible", False)) != feasible:
            raise ValueError(f"support mask {mask} has a mismatched feasibility flag")
        if "residual_mass" in stratum and parse_fraction(stratum["residual_mass"]) != M - len(support):
            raise ValueError(f"support mask {mask} has a mismatched residual_mass")

        raw_anchors = stratum.get("anchor_profiles", [])
        roots = stratum.get("root_cell_ledger", [])
        covered = stratum.get("covered_cell_ledger")
        failed = stratum.get("failed_cell_ledger", [])
        if not isinstance(roots, list) or not isinstance(failed, list):
            raise ValueError(f"support mask {mask} has malformed cell ledgers")
        if failed:
            raise ValueError(f"support mask {mask} contains failed cells")
        if not feasible:
            if raw_anchors or roots or covered:
                raise ValueError(f"infeasible support mask {mask} contains geometry")
            if int(stratum.get("integer_profile_count", 0)) != 0:
                raise ValueError(f"infeasible support mask {mask} has nonzero profile count")
            geometry_rows.append(
                {"support_mask": mask, "support": list(support), "dimension": dimension, "feasible": False}
            )
            continue
        if not bool(stratum.get("complete_cover", False)):
            raise ValueError(f"support mask {mask} does not assert complete_cover")
        anchors = parse_anchor_profiles(raw_anchors, f"support mask {mask}")
        expected_vertices, expected_facets = integer_hull_vertices_from_certificate(
            support,
            stratum.get("exact_hull_vertices"),
            stratum.get("exact_hull_facets"),
        )

        rows = covered if covered is not None else roots
        if not isinstance(rows, list) or not rows:
            raise ValueError(f"feasible support mask {mask} has no covered cells")
        if covered is not None and roots:
            root_geometry = {
                (cell_key(row), tuple(int(value) for value in row.get("anchor_indices", [])))
                for row in roots
            }
            covered_geometry = {
                (cell_key(row), tuple(int(value) for value in row.get("anchor_indices", [])))
                for row in covered
            }
            if root_geometry != covered_geometry:
                raise ValueError(f"support mask {mask} root/covered cell geometry mismatch")
        identifiers = [cell_key(row) for row in rows]
        if any(identifier in {"", "None"} for identifier in identifiers) or len(set(identifiers)) != len(rows):
            raise ValueError(f"support mask {mask} cell identifiers are not unique")
        refinement_geometry = verify_pure_triangle_star_refinements(
            stratum, anchors, rows, support
        )
        geometry = verify_support_stratum_mesh(
            support, anchors, rows, expected_vertices, expected_facets
        )
        if refinement_geometry is not None:
            geometry["pure_triangle_star_refinement"] = refinement_geometry
        geometry["support_mask"] = mask
        geometry_rows.append(geometry)
        for row in sorted(rows, key=cell_key):
            contexts.append(
                {
                    "identifier": f"s{mask:02x}:{cell_key(row)}",
                    "row": row,
                    "anchors": anchors,
                    "support": support,
                    "dimension": dimension,
                }
            )
    return contexts, {
        "mode": "31_disjoint_exact_positive_support_strata",
        "strata": geometry_rows,
        "feasible_strata": sum(row["feasible"] for row in geometry_rows),
        "covered_cells": len(contexts),
    }


def barycentric_coordinates(
    point: tuple[Fraction, ...],
    vertices: tuple[int, ...],
    anchors: tuple[tuple[Fraction, ...], ...],
) -> tuple[Fraction, ...]:
    """Return exact barycentric coordinates of ``point`` in a 4-simplex."""

    base = anchors[vertices[0]][1:]
    columns = [
        tuple(anchors[vertices[index]][coordinate] - anchors[vertices[0]][coordinate] for coordinate in range(1, 5))
        for index in range(1, 5)
    ]
    matrix = [[columns[column][row] for column in range(4)] for row in range(4)]
    denominator = determinant(matrix)
    if not denominator:
        raise ValueError("cannot compute barycentric coordinates in a degenerate simplex")
    target = [point[coordinate] - anchors[vertices[0]][coordinate] for coordinate in range(1, 5)]
    tail = []
    for column in range(4):
        replaced = [list(row) for row in matrix]
        for row in range(4):
            replaced[row][column] = target[row]
        tail.append(determinant(replaced) / denominator)
    return (Fraction(1) - sum(tail, Fraction(0)), *tail)


def boundary_type(facet: tuple[int, ...], anchors) -> str | None:
    profiles = [anchors[index] for index in facet]
    equalities = []
    for coordinate in range(5):
        if all(profile[coordinate] == 0 for profile in profiles):
            equalities.append(f"a{coordinate}=0")
    if all(
        sum(weight * profile[weight] for weight in range(5)) == MINIMUM_PHYSICAL_WEIGHT
        for profile in profiles
    ):
        equalities.append("physical_weight=21")
    if len(equalities) != 1:
        return None
    return equalities[0]


def verify_mesh(ledger: dict[str, Any], anchors, rows) -> dict[str, Any]:
    facets: Counter[tuple[int, ...]] = Counter()
    orientations: Counter[tuple[int, ...]] = Counter()
    seen_simplices = set()
    volume = Fraction(0)
    for row in rows:
        indices = simplex_indices(row, len(anchors))
        key = tuple(sorted(indices))
        if key in seen_simplices:
            raise ValueError(f"duplicate simplex {key}")
        seen_simplices.add(key)
        det = simplex_determinant(indices, anchors)
        if not det:
            raise ValueError(f"degenerate simplex {key}")
        oriented = indices if det > 0 else (indices[1], indices[0], *indices[2:])
        volume += abs(det)
        for omitted in range(5):
            face_order = oriented[:omitted] + oriented[omitted + 1 :]
            face = tuple(sorted(face_order))
            facets[face] += 1
            orientations[face] += ((-1) ** omitted) * permutation_sign(face_order)

    boundary_counts = Counter()
    interior = 0
    for facet, multiplicity in facets.items():
        if multiplicity == 2:
            if orientations[facet] != 0:
                raise ValueError(f"interior facet orientation mismatch {facet}")
            interior += 1
        elif multiplicity == 1:
            kind = boundary_type(facet, anchors)
            if kind is None:
                raise ValueError(f"unpaired facet is not on one exact clipped-hull boundary: {facet}")
            boundary_counts[kind] += 1
        else:
            raise ValueError(f"facet multiplicity {multiplicity} is not conforming")

    expected_volume = Fraction(M**4) - Fraction(MINIMUM_PHYSICAL_WEIGHT**4, math.factorial(4))
    if volume != expected_volume:
        raise ValueError(f"exact mesh volume mismatch: {volume} != {expected_volume}")
    required_boundaries = {f"a{index}=0" for index in range(5)} | {"physical_weight=21"}
    if set(boundary_counts) != required_boundaries:
        raise ValueError("mesh does not expose all six exact clipped-hull boundary facets")
    return {
        "simplices": len(rows),
        "interior_facets": interior,
        "boundary_facets": dict(sorted(boundary_counts.items())),
        "determinant_volume": str(volume),
        "clipped_hull_determinant_volume": str(expected_volume),
    }


def source_paths(ledger_path: Path, ledger: dict[str, Any]) -> list[Path]:
    rows = ledger.get("sources")
    if not isinstance(rows, list) or not rows or not all(isinstance(row, dict) for row in rows):
        raise ValueError("certifiable sources must be {path,sha256} objects, not path strings")
    result = []
    for row in rows:
        raw = Path(str(row.get("path", "")).replace("\\", "/"))
        # Remote proof runs record absolute source paths.  A downloaded ledger
        # remains independently replayable when the digest-matching source is
        # placed beside it under the same basename.
        candidates = (
            raw,
            ledger_path.parent / raw,
            Path.cwd() / raw,
            ledger_path.parent / raw.name,
            Path.cwd() / raw.name,
        )
        path = next((candidate for candidate in candidates if candidate.exists()), None)
        if path is None:
            raise ValueError(f"cannot resolve witness source {raw}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != str(row.get("sha256", "")):
            raise ValueError(f"witness source digest mismatch for {path}")
        result.append(path)
    return result


def parse_mixture(row: dict[str, Any]) -> tuple[tuple[str, Fraction], ...]:
    raw = row.get("mixture")
    if not isinstance(raw, list) or not raw:
        raise ValueError(f"cell {row.get('cell_id')} lacks a fixed mixture")
    combined: dict[str, Fraction] = {}
    for component in raw:
        if not isinstance(component, dict):
            raise ValueError("mixture components must be objects")
        name = str(component.get("witness", component.get("witness_reference", "")))
        if not name:
            raise ValueError("mixture component lacks a witness reference")
        exact = component.get("weight_exact", component.get("weight"))
        weight = parse_fraction(exact)
        if weight < 0:
            raise ValueError("mixture weights must be nonnegative")
        combined[name] = combined.get(name, Fraction(0)) + weight
    if sum(combined.values(), Fraction(0)) != 1:
        raise ValueError("fixed mixture weights must sum exactly to one")
    return tuple(sorted((name, weight) for name, weight in combined.items() if weight))


def configure_generic_hardener() -> None:
    """Configure the shared outward implementation for g=4 before any use."""

    g2.GROUP_BITS = GROUP_BITS
    g2.BLOCK_ATOMS = BLOCK_BITS // GROUP_BITS
    g2._OUTWARD_NORMALIZATION_CACHE.clear()


def harden_full_bijection() -> dict[str, Any]:
    constant = Interval.exact(K) + log2_int(11).times_int(N) - log2_int(10).times_int(N - D)
    zero = Interval.exact(0)
    return {
        "name": "full_bijection",
        "constant": constant,
        "charges": tuple(zero for _ in range(5)),
        "subtract_normalization": True,
        "report": {"name": "full_bijection", "kind": "bijection"},
    }


def scale_interval(value: Interval, factor: Fraction) -> Interval:
    return value * (
        Interval.exact(factor.numerator) / Interval.exact(factor.denominator)
    )


def outward_fractional_normalization(profile: tuple[Fraction, ...]) -> Interval:
    """Rigorous coarse multinomial continuation at rational mesh anchors.

    For a nonintegral count, log-convexity of Gamma bounds log Gamma(x+1)
    above by the chord between its adjacent integer factorial values.  The
    universal lower bound ``log Gamma(x+1) > -1`` for x>=0 is deliberately
    loose but rigorous and only widens the unused lower branch endpoint.
    """

    natural = ln_factorial(M)
    for count in profile:
        if count.denominator == 1:
            gamma = ln_factorial(count.numerator)
        else:
            floor = count.numerator // count.denominator
            fraction = count - floor
            left = ln_factorial(floor)
            right = ln_factorial(floor + 1)
            chord = scale_interval(left, 1 - fraction) + scale_interval(right, fraction)
            gamma = Interval(Decimal(-1), chord.hi)
        natural = natural - gamma
    result = natural / LN2
    for count, classes in zip(profile, profile_classes(GROUP_BITS)):
        if count:
            result = result + scale_interval(log2_int(classes), count)
    return result


def evaluate_profile(profile, mixture, hardened) -> tuple[Interval, list[Any]]:
    total = Interval.exact(0)
    components = []
    for name, weight in mixture:
        value = hardened[name]["constant"]
        for count, charge in zip(profile, hardened[name]["charges"]):
            if not count:
                continue
            if charge is None:
                raise ValueError(f"witness {name!r} is support-ineligible at {profile}")
            value = value - charge * (
                Interval.exact(count.numerator) / Interval.exact(count.denominator)
            )
        if hardened[name].get("subtract_normalization", True):
            value = value - outward_fractional_normalization(profile)
        weight_interval = Interval.exact(weight.numerator) / Interval.exact(weight.denominator)
        total = total + value * weight_interval
        components.append((name, weight, value))
    return total, components


def require_support_eligible(
    identifier: str,
    support: tuple[int, ...],
    mixture: tuple[tuple[str, Fraction], ...],
    hardened: dict[str, Any],
) -> None:
    for name, _weight in mixture:
        ineligible = [
            coordinate for coordinate in support
            if hardened[name]["charges"][coordinate] is None
        ]
        if ineligible:
            raise ValueError(
                f"cell {identifier!r} mixture witness {name!r} is "
                f"ineligible on active coordinates {ineligible}"
            )


def support_self_test() -> None:
    """Small exact smokes for dimensions 0, 1, and the full-support 4-hull."""

    full = (0, 1, 2, 3, 4)
    candidates = integer_hull_candidates(full)
    # Use discovery code only to supply a claimed vertex list to the smoke;
    # the verifier below independently rebuilds every exact supporting facet
    # and checks containment of its separately enumerated candidate set.
    from build_packet_group_g4_support_regular_mesh import integer_hull_vertices

    anchors = integer_hull_vertices(full)
    facets = exact_integer_hull_facets(anchors, full)
    assert all(point_is_inside_exact_facets(point, anchors, full, facets) for point in candidates)
    cells = integer_hull_pulling_cells(anchors, full, facets)
    report = verify_support_stratum_mesh(
        full,
        anchors,
        [{"cell_id": f"t{index}", "anchor_indices": list(cell)} for index, cell in enumerate(cells)],
        anchors,
        facets,
    )
    assert len(anchors) == 15 and report["cells"] == 18
    segment_support = (0, 1)
    segment = integer_hull_candidates(segment_support)
    segment_facets = exact_integer_hull_facets(segment, segment_support)
    verify_support_stratum_mesh(
        segment_support, segment, [{"cell_id": "segment", "anchor_indices": [0, 1]}],
        segment, segment_facets,
    )
    point_support = (1,)
    point = integer_hull_candidates(point_support)
    verify_support_stratum_mesh(
        point_support, point, [{"cell_id": "point", "anchor_indices": [0]}], point, ()
    )
    assert not integer_hull_candidates((0,))
    for mask in range(1, 32):
        support = tuple(index for index in range(5) if mask & (1 << index))
        assert support_hull_determinant_volume(support) >= 0
    try:
        require_support_eligible(
            "reject", (0, 2), (("w", Fraction(1)),),
            {"w": {"charges": (Interval.exact(0), None, None, None, None)}},
        )
    except ValueError:
        pass
    else:
        raise AssertionError("support-ineligible mixture was accepted")

    # Local m=3 ledger smoke on a one-simplex domain.  This checks exact node,
    # chamber, complement, and volume reconstruction only; production ledgers
    # additionally pass the mandatory global terminal-mesh conformity audit.
    residual = M - 5
    pure = []
    for favored in range(5):
        profile = [Fraction(1) for _ in range(5)]
        profile[favored] += residual
        pure.append(tuple(profile))
    triangle = (1, 2, 3)
    grid_anchors = list(pure)
    grid_nodes = []
    node_indices = {}
    for numerators in sorted(
        (left, middle, 3 - left - middle)
        for left in range(4)
        for middle in range(4 - left)
    ):
        profile = tuple(
            sum(
                (Fraction(numerators[position], 3) * pure[triangle[position]][coordinate]
                 for position in range(3)),
                Fraction(0),
            )
            for coordinate in range(5)
        )
        if profile in grid_anchors:
            index = grid_anchors.index(profile)
        else:
            index = len(grid_anchors)
            grid_anchors.append(profile)
        node_indices[numerators] = index
        grid_nodes.append({"barycentric_numerators": list(numerators), "anchor_index": index})
    face_triangles = [
        numerators
        for numerators in itertools.combinations(sorted(node_indices), 3)
        if all(
            sum(abs(left[position] - right[position]) for position in range(3)) == 2
            for left, right in itertools.combinations(numerators, 2)
        )
    ]
    parent = {"cell_id": "parent", "anchor_indices": [0, 1, 2, 3, 4]}
    children = [
        {
            "cell_id": f"parent.ptg3t{ordinal:02d}",
            "parent_cell_id": "parent",
            "anchor_indices": [0, 4, *(node_indices[item] for item in face)],
        }
        for ordinal, face in enumerate(face_triangles)
    ]
    refinement = {
        "unrefined_root_cell_ledger": [parent],
        "pure_triangle_star_refinement_ledger": [
            {
                "refinement_id": "smoke_grid",
                "refinement_kind": "barycentric_face_grid",
                "triangle_anchor_indices": list(triangle),
                "grid_denominator": 3,
                "grid_nodes": grid_nodes,
                "parents": [parent],
                "children": children,
            }
        ],
    }
    star = verify_pure_triangle_star_refinements(
        refinement, tuple(grid_anchors), children, full
    )
    assert star == {
        "refinements": 1,
        "replaced_parents": 1,
        "created_children": 9,
        "terminal_cells": 9,
    }
    print("support-stratified exact geometry self-test: PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    self_check()
    if args.self_test:
        support_self_test()
        return
    if args.ledger is None:
        parser.error("--ledger is required unless --self-test is used")
    ledger = json.loads(args.ledger.read_text(encoding="utf-8"))
    if int(ledger.get("group_bits", -1)) != GROUP_BITS:
        raise SystemExit("g4 anchor verifier: ledger group_bits must equal 4")

    support_mode = support_stratum_contexts(ledger)
    if support_mode is not None:
        contexts, geometry = support_mode
        owner_rule = SUPPORT_OWNER_RULE
        if ledger.get("profile_owner_rule") != owner_rule:
            raise SystemExit(f"g4 anchor verifier: profile_owner_rule must be {owner_rule!r}")
        if ledger.get("union_accounting_mode") != "cell_local_profile_count_upper":
            raise ValueError("support-stratified covers require cell-local union accounting")
    else:
        anchors = load_anchors(ledger)
        stellar = stellar_leaf_rows(ledger, anchors)
        if stellar is None:
            owner_rule = FLAT_OWNER_RULE
            if ledger.get("profile_owner_rule") != owner_rule:
                raise SystemExit(f"g4 anchor verifier: profile_owner_rule must be {owner_rule!r}")
            rows = cell_rows(ledger)
            geometry = {"mode": "flat_exact_simplicial_mesh", **verify_mesh(ledger, anchors, rows)}
        else:
            owner_rule = STELLAR_OWNER_RULE
            if ledger.get("profile_owner_rule") != owner_rule:
                raise SystemExit(f"g4 anchor verifier: profile_owner_rule must be {owner_rule!r}")
            rows, geometry = stellar
        contexts = [
            {"identifier": cell_key(row), "row": row, "anchors": anchors,
             "support": None, "dimension": DIMENSION}
            for row in rows
        ]
    paths = source_paths(args.ledger, ledger)
    witnesses = g2.load_witnesses(paths)
    mixtures = [(context, parse_mixture(context["row"])) for context in contexts]
    used = sorted({name for _context, mixture in mixtures for name, _weight in mixture})
    missing = [name for name in used if name != "full_bijection" and name not in witnesses]
    if missing:
        raise ValueError("missing witness references: " + ", ".join(missing))

    configure_generic_hardener()
    split_caps = g2.split_cap_table()
    hardened = {}
    for name in used:
        if name == "full_bijection":
            hardened[name] = harden_full_bijection()
        else:
            row, path = witnesses[name]
            hardened[name] = g2.harden_witness(name, row, path, split_caps, args.iterations)

    # Eligibility is a property of a whole positive-support stratum, not just
    # of its sampled vertices.  A positive mixture component must therefore
    # have a finite charge on every active class.
    for context, mixture in mixtures:
        support = context["support"]
        if support is None:
            continue
        require_support_eligible(context["identifier"], support, mixture, hardened)

    maximum = None
    worst = None
    digest = hashlib.sha256()
    component_evaluations = 0
    cell_records = []
    for context, mixture in mixtures:
        row = context["row"]
        anchors = context["anchors"]
        identifier = context["identifier"]
        indices = stratum_cell_indices(row, len(anchors), context["dimension"])
        cell_maximum = None
        for vertex in indices:
            value, components = evaluate_profile(anchors[vertex], mixture, hardened)
            component_evaluations += len(components)
            record = [
                identifier,
                vertex,
                [[name, str(weight), str(part.lo), str(part.hi)] for name, weight, part in components],
                str(value.lo),
                str(value.hi),
            ]
            digest.update(json.dumps(record, separators=(",", ":")).encode())
            if cell_maximum is None:
                cell_maximum = value
            else:
                cell_maximum = Interval(
                    max(cell_maximum.lo, value.lo),
                    max(cell_maximum.hi, value.hi),
                )
            if maximum is None or value.hi > maximum.hi:
                maximum = value
                worst = (identifier, vertex, mixture)
        assert cell_maximum is not None
        local_count, omitted, ranges = simplex_integer_profile_count_upper(
            indices, anchors
        )
        reported_count = row.get("integer_profile_count_upper")
        if reported_count is not None and int(reported_count) != local_count:
            raise ValueError(
                f"cell {identifier!r} integer_profile_count_upper mismatch: "
                f"{reported_count} != {local_count}"
            )
        reported_method = row.get("integer_profile_count_method")
        if reported_method is not None and str(reported_method) != CELL_COUNT_METHOD:
            raise ValueError(
                f"cell {identifier!r} has unsupported integer_profile_count_method "
                f"{reported_method!r}"
            )
        if local_count <= 0:
            # A zero count proves that this real cell contains no integer
            # profile and contributes no term to the union ledger.
            cell_term = None
        else:
            cell_term = cell_maximum + log2_int(local_count)
        count_record = [
            identifier,
            str(local_count),
            omitted,
            [list(pair) for pair in ranges],
            str(cell_maximum.lo),
            str(cell_maximum.hi),
            None if cell_term is None else [str(cell_term.lo), str(cell_term.hi)],
        ]
        digest.update(json.dumps(count_record, separators=(",", ":")).encode())
        cell_records.append(
            {
                "cell_id": identifier,
                "integer_profile_count_upper": local_count,
                "omitted_coordinate": omitted,
                "coordinate_ranges": [list(pair) for pair in ranges],
                "maximum_branch_log2_interval": [
                    str(cell_maximum.lo),
                    str(cell_maximum.hi),
                ],
                "union_term_log2_interval": (
                    None if cell_term is None else [str(cell_term.lo), str(cell_term.hi)]
                ),
                "_term": cell_term,
            }
        )
    assert maximum is not None and worst is not None

    count = profile_count(GROUP_BITS, N)
    count_log = log2_int(count)
    accounting_mode = str(ledger.get("union_accounting_mode", "global_profile_count"))
    if accounting_mode == "global_profile_count":
        union = maximum + count_log
        worst_contribution = None
        union_description = "one owner per integer profile; global profile-count upper"
    elif accounting_mode == "cell_local_profile_count_upper":
        terms = [record["_term"] for record in cell_records if record["_term"] is not None]
        if not terms:
            raise ValueError("cell-local union ledger has no nonempty integer-profile cells")
        union = g2._log2_sum_exp(terms)
        worst_contribution = max(
            (record for record in cell_records if record["_term"] is not None),
            key=lambda record: record["_term"].hi,
        )
        union_description = (
            "outward log2 sum_c n_c*2^U_c; n_c is the verifier's best "
            "four-coordinate exact bounding-box upper"
        )
    else:
        raise ValueError(f"unsupported union_accounting_mode {accounting_mode!r}")
    margin = Decimal(-40) - union.hi
    passed = margin >= 0
    report = {
        "status": "OUTWARD_CERTIFIED_G4_ANCHOR_MESH_2^-40" if passed else "OUTWARD_G4_ANCHOR_MESH_DID_NOT_CLOSE",
        "group_bits": GROUP_BITS,
        "passed": passed,
        "exact_geometry": geometry,
        "profile_owner_rule": owner_rule,
        "integer_profile_completeness": (
            "every feasible integer profile has one unique positive-coordinate support; "
            "each exact support hull is partitioned and its minimum closed-cell id owns boundaries"
            if support_mode is not None
            else "exact real partition contains every clipped-hull integer profile; "
            "minimum closed-cell id assigns shared-boundary profiles uniquely"
        ),
        "used_witnesses": len(used),
        "used_fixed_mixtures": len(set(mixture for _context, mixture in mixtures)),
        "used_cell_vertex_inequalities": sum(
            context["dimension"] + 1 for context, _mixture in mixtures
        ),
        "used_component_vertex_evaluations": component_evaluations,
        "inequality_sha256": digest.hexdigest(),
        "worst_cell_id": worst[0],
        "worst_anchor_index": worst[1],
        "worst_mixture": [{"witness": name, "weight": str(weight)} for name, weight in worst[2]],
        "maximum_branch_log2_interval": [str(maximum.lo), str(maximum.hi)],
        "profile_count_upper": count,
        "profile_count_log2_interval": [str(count_log.lo), str(count_log.hi)],
        "union_accounting_mode": accounting_mode,
        "union_accounting": union_description,
        "cell_local_count_method": CELL_COUNT_METHOD,
        "sum_of_cell_integer_profile_count_uppers": (
            sum(record["integer_profile_count_upper"] for record in cell_records)
            if accounting_mode == "cell_local_profile_count_upper"
            else None
        ),
        "cell_local_union_terms": [
            {key: value for key, value in record.items() if key != "_term"}
            for record in cell_records
        ] if accounting_mode == "cell_local_profile_count_upper" else None,
        "worst_cell_union_contribution": (
            None
            if worst_contribution is None
            else {key: value for key, value in worst_contribution.items() if key != "_term"}
        ),
        "union_log2_interval": [str(union.lo), str(union.hi)],
        "certified_margin_bits": str(margin),
        "witness_reports": [hardened[name]["report"] for name in used],
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    if not passed:
        raise SystemExit("g4 anchor verifier: outward 40-bit union did not close")


if __name__ == "__main__":
    main()
