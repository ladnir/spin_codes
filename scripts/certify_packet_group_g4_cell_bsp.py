#!/usr/bin/env python3
"""Independently replay one exact per-cell g=4 dominance BSP.

The source support-stratified mesh is fully audited first.  The BSP then
partitions one already-certified root simplex by exact rational half-spaces.
Every leaf owner is outward-hardened and evaluated at every exact leaf vertex.
This script certifies the selected cell bound only; it does not claim the
complete g=4 union until such records are integrated for every required cell.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import certify_packet_group_g4_anchor_mesh as g4
import certify_packet_group_triangle_ledger as g2
from outward_log2 import Interval, log2_int, self_check
from packet_group_drive_stratified import profile_count


DIMENSION = 4
CLASSES = 5
SUPPORT = (0, 1, 2, 3, 4)
CHART_INDICES = (0, 1, 2, 3)
DEPENDENT_INDEX = 4
BSP_SCHEMA = "packet-group-g4-cell-dominance-bsp-v1"
Profile = tuple[Fraction, ...]
Constraint = tuple[tuple[Fraction, ...], Fraction]
Affine = tuple[Fraction, tuple[Fraction, ...]]


def configure_support(support_mask: int) -> None:
    """Select the affine chart for one exact positive-support stratum."""

    global DIMENSION, SUPPORT, CHART_INDICES, DEPENDENT_INDEX
    support = tuple(index for index in range(CLASSES) if support_mask & (1 << index))
    if len(support) < 2:
        raise ValueError("g4 cell BSP verifier requires positive dimension")
    SUPPORT = support
    DIMENSION = len(support) - 1
    CHART_INDICES = support[:-1]
    DEPENDENT_INDEX = support[-1]


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def solve_linear(
    matrix: Iterable[Iterable[Fraction]], right: Iterable[Fraction]
) -> tuple[Fraction, ...] | None:
    rows = [list(map(Fraction, row)) + [Fraction(value)] for row, value in zip(matrix, right)]
    n = len(rows)
    if n == 0 or any(len(row) != n + 1 for row in rows):
        raise ValueError("g4 cell BSP verifier: malformed linear system")
    for column in range(n):
        pivot = next((row for row in range(column, n) if rows[row][column]), None)
        if pivot is None:
            return None
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [value / scale for value in rows[column]]
        for row in range(n):
            if row == column or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                left - factor * right_value
                for left, right_value in zip(rows[row], rows[column])
            ]
    return tuple(row[-1] for row in rows)


def simplex_constraints(vertices: tuple[Profile, ...]) -> tuple[Constraint, ...]:
    augmented = [
        [vertex[index] for index in CHART_INDICES] + [Fraction(1)]
        for vertex in vertices
    ]
    if len(augmented) != DIMENSION + 1:
        raise ValueError("g4 cell BSP verifier: root has the wrong dimension")
    constraints = []
    for owner in range(DIMENSION + 1):
        target = [Fraction(int(index == owner)) for index in range(DIMENSION + 1)]
        barycentric = solve_linear(augmented, target)
        if barycentric is None:
            raise ValueError("g4 cell BSP verifier: degenerate root simplex")
        constraints.append(
            (tuple(-value for value in barycentric[:DIMENSION]), barycentric[-1])
        )
    return tuple(constraints)


def enumerate_vertices(constraints: tuple[Constraint, ...]) -> tuple[Profile, ...]:
    vertices: set[Profile] = set()
    for active in itertools.combinations(range(len(constraints)), DIMENSION):
        solution = solve_linear(
            [constraints[index][0] for index in active],
            [constraints[index][1] for index in active],
        )
        if solution is None:
            continue
        if all(
            sum(coefficient * value for coefficient, value in zip(left, solution)) <= right
            for left, right in constraints
        ):
            profile = [Fraction(0)] * CLASSES
            for index, value in zip(CHART_INDICES, solution):
                profile[index] = value
            profile[DEPENDENT_INDEX] = Fraction(g4.M) - sum(solution, Fraction(0))
            vertices.add(tuple(profile))
    return tuple(sorted(vertices))


def parse_affine(row: Any) -> Affine:
    if not isinstance(row, dict):
        raise ValueError("g4 cell BSP verifier: malformed hyperplane")
    coefficients = tuple(g4.parse_fraction(value) for value in row.get("coefficients", []))
    constant = g4.parse_fraction(row.get("constant"))
    if len(coefficients) != CLASSES or not any((constant, *coefficients)):
        raise ValueError("g4 cell BSP verifier: invalid hyperplane")
    expected = "constant+dot(coefficients,profile)<=0"
    if row.get("nonpositive_child_semantics") != expected:
        raise ValueError("g4 cell BSP verifier: hyperplane semantics mismatch")
    return constant, coefficients


def affine_value(affine: Affine, profile: Profile) -> Fraction:
    return affine[0] + sum(
        coefficient * value for coefficient, value in zip(affine[1], profile)
    )


def affine_constraint(affine: Affine, nonpositive: bool) -> Constraint:
    constant, coefficients = affine
    chart = tuple(
        coefficients[index] - coefficients[DEPENDENT_INDEX]
        for index in CHART_INDICES
    )
    chart_constant = constant + coefficients[DEPENDENT_INDEX] * g4.M
    if nonpositive:
        return chart, -chart_constant
    return tuple(-value for value in chart), chart_constant


def serialized_vertices(vertices: tuple[Profile, ...]) -> list[list[str]]:
    return [[fraction_text(value) for value in profile] for profile in vertices]


def parse_serialized_vertices(raw: Any) -> tuple[Profile, ...]:
    if not isinstance(raw, list):
        raise ValueError("g4 cell BSP verifier: malformed serialized vertices")
    result = []
    for profile in raw:
        if not isinstance(profile, list) or len(profile) != CLASSES:
            raise ValueError("g4 cell BSP verifier: malformed serialized profile")
        result.append(tuple(g4.parse_fraction(value) for value in profile))
    return tuple(result)


def resolve_source(artifact_path: Path, row: Any, label: str) -> Path:
    if not isinstance(row, dict):
        raise ValueError(f"g4 cell BSP verifier: missing {label} source record")
    raw = Path(str(row.get("path", "")).replace("\\", "/"))
    candidates = (
        raw,
        artifact_path.parent / raw,
        artifact_path.parent / raw.name,
        Path.cwd() / raw.name,
    )
    path = next((candidate for candidate in candidates if candidate.exists()), None)
    if path is None:
        raise ValueError(f"g4 cell BSP verifier: cannot resolve {label} {raw}")
    if file_sha256(path) != str(row.get("sha256", "")):
        raise ValueError(f"g4 cell BSP verifier: {label} digest mismatch")
    return path


def audit_source_stratum(
    ledger: dict[str, Any], support_mask: int, cell_id: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Audit one exact support stratum without claiming the other 30 close."""

    raw = ledger.get("support_strata")
    if not isinstance(raw, list):
        raise ValueError("g4 cell BSP verifier: source lacks support strata")
    stratum = next(
        (row for row in raw if int(row.get("support_mask", -1)) == support_mask), None
    )
    if stratum is None:
        raise ValueError("g4 cell BSP verifier: requested support stratum is absent")
    support = tuple(index for index in range(CLASSES) if support_mask & (1 << index))
    dimension = len(support) - 1
    if (
        stratum.get("support") != list(support)
        or int(stratum.get("dimension", -1)) != dimension
        or not bool(stratum.get("feasible", False))
        or not bool(stratum.get("complete_cover", False))
        or stratum.get("failed_cell_ledger", [])
        or stratum.get("profile_owner_rule") != g4.SUPPORT_OWNER_RULE
        or stratum.get("union_accounting_mode") != "cell_local_profile_count_upper"
    ):
        raise ValueError("g4 cell BSP verifier: malformed or incomplete source stratum")
    anchors = g4.parse_anchor_profiles(
        stratum.get("anchor_profiles", []), f"support mask {support_mask}"
    )
    expected_vertices, expected_facets = g4.integer_hull_vertices_from_certificate(
        support,
        stratum.get("exact_hull_vertices"),
        stratum.get("exact_hull_facets"),
    )
    rows = stratum.get("covered_cell_ledger")
    if not isinstance(rows, list) or not rows:
        raise ValueError("g4 cell BSP verifier: source stratum has no covered cells")
    geometry = g4.verify_support_stratum_mesh(
        support, anchors, rows, expected_vertices, expected_facets
    )
    matches = [row for row in rows if str(row.get("cell_id")) == cell_id]
    if len(matches) != 1:
        raise ValueError("g4 cell BSP verifier: selected source cell is not unique")
    return {
        "identifier": f"s{support_mask:02x}:{cell_id}",
        "row": matches[0],
        "anchors": anchors,
        "support": support,
        "dimension": dimension,
    }, {"mode": "single_exact_positive_support_stratum", **geometry}


def replay_tree(
    artifact: dict[str, Any], root_constraints: tuple[Constraint, ...]
) -> tuple[list[dict[str, Any]], dict[str, tuple[Profile, ...]]]:
    raw_nodes = artifact.get("nodes")
    if not isinstance(raw_nodes, list) or not raw_nodes:
        raise ValueError("g4 cell BSP verifier: empty node ledger")
    nodes: dict[str, dict[str, Any]] = {}
    for row in raw_nodes:
        if not isinstance(row, dict):
            raise ValueError("g4 cell BSP verifier: node rows must be objects")
        identifier = str(row.get("node_id", ""))
        if not identifier or identifier in nodes:
            raise ValueError("g4 cell BSP verifier: duplicate or empty node id")
        nodes[identifier] = row

    used: set[str] = set()
    leaves: list[dict[str, Any]] = []
    leaf_vertices: dict[str, tuple[Profile, ...]] = {}

    def visit(identifier: str, constraints: tuple[Constraint, ...], depth: int) -> None:
        if identifier in used or identifier not in nodes:
            raise ValueError("g4 cell BSP verifier: cyclic, repeated, or missing child")
        used.add(identifier)
        row = nodes[identifier]
        if int(row.get("depth", -1)) != depth:
            raise ValueError(f"g4 cell BSP verifier: node {identifier} depth mismatch")
        vertices = enumerate_vertices(constraints)
        if len(vertices) < DIMENSION + 1:
            raise ValueError(f"g4 cell BSP verifier: node {identifier} is not full-dimensional")
        kind = row.get("kind")
        if kind == "leaf":
            owner = str(row.get("owner_witness", ""))
            if not owner:
                raise ValueError(f"g4 cell BSP verifier: leaf {identifier} lacks an owner")
            if int(row.get("vertex_count", -1)) != len(vertices):
                raise ValueError(f"g4 cell BSP verifier: leaf {identifier} vertex count mismatch")
            if str(row.get("vertices_sha256", "")) != compact_hash(serialized_vertices(vertices)):
                raise ValueError(f"g4 cell BSP verifier: leaf {identifier} vertex digest mismatch")
            leaves.append(row)
            leaf_vertices[identifier] = vertices
            return
        if kind != "split":
            raise ValueError(f"g4 cell BSP verifier: node {identifier} has invalid kind")
        affine = parse_affine(row.get("hyperplane"))
        signs = [affine_value(affine, profile) for profile in vertices]
        if not any(value < 0 for value in signs) or not any(value > 0 for value in signs):
            raise ValueError(f"g4 cell BSP verifier: split {identifier} does not cut the interior")
        negative = str(row.get("nonpositive_child", ""))
        positive = str(row.get("nonnegative_child", ""))
        visit(negative, constraints + (affine_constraint(affine, True),), depth + 1)
        visit(positive, constraints + (affine_constraint(affine, False),), depth + 1)

    visit(str(artifact.get("root_node", "")), root_constraints, 0)
    if used != set(nodes):
        raise ValueError("g4 cell BSP verifier: unreachable nodes are present")
    return leaves, leaf_vertices


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bsp", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    self_check()

    artifact = json.loads(args.bsp.read_text(encoding="utf-8"))
    if artifact.get("schema") != BSP_SCHEMA or int(artifact.get("group_bits", -1)) != 4:
        raise ValueError("g4 cell BSP verifier: unsupported artifact schema")
    if not bool(artifact.get("complete_exact_partition_claimed", False)):
        raise ValueError("g4 cell BSP verifier: producer did not claim a complete partition")

    ledger_path = resolve_source(args.bsp, artifact.get("source_ledger"), "source ledger")
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    expected_identifier = f"s{int(artifact['support_mask']):02x}:{artifact['cell_id']}"
    context, geometry = audit_source_stratum(
        ledger, int(artifact["support_mask"]), str(artifact["cell_id"])
    )
    configure_support(int(artifact["support_mask"]))
    if context["identifier"] != expected_identifier:
        raise ValueError("g4 cell BSP verifier: selected cell identifier mismatch")
    row = context["row"]
    anchors = context["anchors"]
    indices = g4.stratum_cell_indices(row, len(anchors), context["dimension"])
    if list(indices) != list(map(int, artifact.get("root_anchor_indices", []))):
        raise ValueError("g4 cell BSP verifier: root anchor indices mismatch")
    root_vertices = tuple(anchors[index] for index in indices)
    if root_vertices != parse_serialized_vertices(artifact.get("root_anchor_profiles")):
        raise ValueError("g4 cell BSP verifier: root anchor profiles mismatch")
    root_constraints = simplex_constraints(root_vertices)
    if enumerate_vertices(root_constraints) != tuple(sorted(root_vertices)):
        raise ValueError("g4 cell BSP verifier: root reconstruction mismatch")

    local_count, omitted, ranges = g4.simplex_integer_profile_count_upper(indices, anchors)
    if local_count != int(row.get("integer_profile_count_upper", -1)):
        raise ValueError("g4 cell BSP verifier: source cell count mismatch")
    if local_count != int(artifact.get("root_integer_profile_count_upper", -1)):
        raise ValueError("g4 cell BSP verifier: BSP cell count mismatch")

    leaves, vertices_by_leaf = replay_tree(artifact, root_constraints)
    paths = g4.source_paths(ledger_path, ledger)
    witnesses = g2.load_witnesses(paths)
    used = sorted({str(leaf["owner_witness"]) for leaf in leaves})
    missing = [name for name in used if name != "full_bijection" and name not in witnesses]
    if missing:
        raise ValueError("g4 cell BSP verifier: missing owners " + ", ".join(missing))

    g4.configure_generic_hardener()
    split_caps = g2.split_cap_table()
    hardened: dict[str, Any] = {}
    for name in used:
        if name == "full_bijection":
            hardened[name] = g4.harden_full_bijection()
        else:
            witness_row, witness_path = witnesses[name]
            hardened[name] = g2.harden_witness(
                name, witness_row, witness_path, split_caps, args.iterations
            )
        g4.require_support_eligible(
            expected_identifier, context["support"], ((name, Fraction(1)),), hardened
        )

    digest = hashlib.sha256()
    maximum: Interval | None = None
    worst: tuple[str, int, str] | None = None
    evaluations = 0
    leaf_reports = []
    for leaf in sorted(leaves, key=lambda value: str(value["node_id"])):
        identifier = str(leaf["node_id"])
        owner = str(leaf["owner_witness"])
        leaf_maximum: Interval | None = None
        for vertex_index, profile in enumerate(vertices_by_leaf[identifier]):
            value, components = g4.evaluate_profile(
                profile, ((owner, Fraction(1)),), hardened
            )
            evaluations += 1
            digest.update(
                json.dumps(
                    [
                        identifier,
                        vertex_index,
                        serialized_vertices((profile,))[0],
                        owner,
                        str(value.lo),
                        str(value.hi),
                    ],
                    separators=(",", ":"),
                ).encode()
            )
            leaf_maximum = value if leaf_maximum is None else Interval(
                max(leaf_maximum.lo, value.lo), max(leaf_maximum.hi, value.hi)
            )
            if maximum is None or value.hi > maximum.hi:
                maximum = value
                worst = (identifier, vertex_index, owner)
        assert leaf_maximum is not None
        leaf_reports.append(
            {
                "node_id": identifier,
                "owner_witness": owner,
                "vertex_count": len(vertices_by_leaf[identifier]),
                "maximum_log2_interval": [str(leaf_maximum.lo), str(leaf_maximum.hi)],
            }
        )
    assert maximum is not None and worst is not None

    global_count = profile_count(4, g4.N)
    uniform_target = Interval.exact(-40) - log2_int(global_count)
    cell_term = maximum + log2_int(local_count)
    report = {
        "status": "OUTWARD_CERTIFIED_G4_SINGLE_CELL_BSP",
        "scope": "one source-mesh cell only; not a complete g4 union certificate",
        "group_bits": 4,
        "source_ledger_sha256": file_sha256(ledger_path),
        "source_geometry": geometry,
        "cell_id": expected_identifier,
        "root_anchor_indices": list(indices),
        "integer_profile_count_upper": local_count,
        "omitted_coordinate": omitted,
        "coordinate_ranges": [list(pair) for pair in ranges],
        "exact_partition_replayed": True,
        "bsp_nodes": len(artifact["nodes"]),
        "bsp_leaves": len(leaves),
        "used_owner_witnesses": len(used),
        "used_leaf_vertex_inequalities": evaluations,
        "inequality_sha256": digest.hexdigest(),
        "maximum_branch_log2_interval": [str(maximum.lo), str(maximum.hi)],
        "cell_union_term_log2_interval": [str(cell_term.lo), str(cell_term.hi)],
        "uniform_per_profile_target_log2_interval": [
            str(uniform_target.lo), str(uniform_target.hi)
        ],
        "all_leaves_pass_uniform_target": maximum.hi <= uniform_target.lo,
        "worst_leaf": worst[0],
        "worst_leaf_vertex": worst[1],
        "worst_owner_witness": worst[2],
        "leaf_reports": leaf_reports,
        "witness_reports": [hardened[name]["report"] for name in used],
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    if not report["all_leaves_pass_uniform_target"]:
        raise SystemExit("g4 cell BSP verifier: selected cell does not close")


if __name__ == "__main__":
    main()
