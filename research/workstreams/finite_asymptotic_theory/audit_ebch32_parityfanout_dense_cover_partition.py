#!/usr/bin/env python3
"""Audit the recursive partition represented by the dense-cover receipt."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

from diagnose_ebch32_parityfanout_ba_three_band_cover_all_q import (
    L,
    MIN_OCCUPATION,
    bounding_box_lattice_count,
    initial_tetrahedra,
    subdivide,
)


WORKSTREAM = Path(__file__).resolve().parent
INPUT = WORKSTREAM / os.environ.get(
    "SPIN_EBCH_DENSE_COVER_OUTPUT",
    "ebch32_parityfanout31x33_ba3_B256_three_band_cover_all_q_d11.json",
)
OUTPUT = WORKSTREAM / os.environ.get(
    "SPIN_EBCH_PARTITION_AUDIT_OUTPUT",
    "ebch32_parityfanout31x33_ba3_B256_dense_cover_partition_audit.json",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def signature(vertices: np.ndarray) -> tuple[tuple[str, str, str], ...]:
    rows = [tuple(float(value).hex() for value in row) for row in vertices]
    return tuple(sorted(rows))


def determinant3(matrix: list[list[int]]) -> int:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def six_volume(vertices: np.ndarray) -> int:
    integer_vertices = [[int(value) for value in row] for row in vertices]
    base = integer_vertices[0]
    columns = [
        [integer_vertices[column][row] - base[row] for column in range(1, 4)]
        for row in range(3)
    ]
    return abs(determinant3(columns))


def main() -> None:
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    leaves: dict[tuple[tuple[str, str, str], ...], dict[str, object]] = {}

    def add_leaf(kind: str, vertices, depth: int, index: int) -> None:
        array = np.asarray(vertices, dtype=np.float64)
        key = signature(array)
        if key in leaves:
            raise ArithmeticError(
                f"duplicate leaf geometry: {kind}[{index}] duplicates "
                f"{leaves[key]['kind']}[{leaves[key]['index']}]"
            )
        leaves[key] = {
            "kind": kind,
            "index": index,
            "recorded_depth": depth,
            "vertices": array,
        }

    for index, row in enumerate(payload["tetrahedra"]):
        add_leaf(
            "accepted",
            row["vertices_active_counts"],
            int(row["depth"]),
            index,
        )
    for index, row in enumerate(payload["pending_worklist"]):
        add_leaf("pending", row["vertices"], int(row["depth"]), index)
    for index, row in enumerate(payload["rejected_tetrahedra"]):
        add_leaf("rejected", row["vertices"], int(row["depth"]), index)

    root_cells = initial_tetrahedra()
    expected_root_vertices = [
        {(0, 0, MIN_OCCUPATION), (0, 0, L), (0, L, 0), (L, 0, 0)},
        {(0, 0, MIN_OCCUPATION), (0, L, 0), (L, 0, 0), (0, MIN_OCCUPATION, 0)},
        {(0, 0, MIN_OCCUPATION), (MIN_OCCUPATION, 0, 0), (L, 0, 0), (0, MIN_OCCUPATION, 0)},
    ]
    actual_root_vertices = [
        {tuple(int(value) for value in row) for row in root}
        for root in root_cells
    ]
    if actual_root_vertices != expected_root_vertices:
        raise ArithmeticError("initial tetrahedra do not match the frozen pulling triangulation")
    root_six_volumes = [six_volume(root) for root in root_cells]
    shell_six_volume = L**3 - MIN_OCCUPATION**3
    if sum(root_six_volumes) != shell_six_volume:
        raise ArithmeticError("initial tetrahedron volumes do not sum to the shell volume")
    maximum_depth = max(
        int(row["recorded_depth"]) for row in leaves.values()
    )
    expected_nodes = 2 * len(leaves) - len(root_cells)
    consumed: set[tuple[tuple[str, str, str], ...]] = set()
    stack = [(root, 0) for root in root_cells]
    visited_nodes = 0
    internal_nodes = 0
    depth_histogram: dict[str, int] = {}
    kind_histogram = {"accepted": 0, "pending": 0, "rejected": 0}

    while stack:
        vertices, depth = stack.pop()
        visited_nodes += 1
        if visited_nodes > expected_nodes:
            raise ArithmeticError(
                "receipt leaves do not form the deterministic subdivision tree"
            )
        key = signature(vertices)
        leaf = leaves.get(key)
        if leaf is not None:
            if int(leaf["recorded_depth"]) != depth:
                raise ArithmeticError("recorded leaf depth disagrees with tree depth")
            if key in consumed:
                raise ArithmeticError("subdivision tree reached a leaf twice")
            consumed.add(key)
            kind = str(leaf["kind"])
            kind_histogram[kind] += 1
            depth_histogram[str(depth)] = depth_histogram.get(str(depth), 0) + 1
            continue
        if depth >= maximum_depth:
            raise ArithmeticError("subdivision tree contains an uncovered branch")
        internal_nodes += 1
        stack.extend((child, depth + 1) for child in subdivide(vertices))

    if consumed != set(leaves):
        raise ArithmeticError("receipt contains leaves outside the subdivision tree")
    if visited_nodes != expected_nodes:
        raise ArithmeticError("full binary forest node-count identity failed")
    if internal_nodes != len(leaves) - len(root_cells):
        raise ArithmeticError("full binary forest internal-node identity failed")

    accepted_box_count = 0
    for row in payload["tetrahedra"]:
        vertices = np.asarray(row["vertices_active_counts"], dtype=np.float64)
        geometric_count = bounding_box_lattice_count(vertices)
        mode = row.get("optimizer_mode")
        recorded = int(row["bounding_box_lattice_count"])
        if mode == "enumerated-lattice-box":
            if not 0 < recorded <= geometric_count:
                raise ArithmeticError("enumerated-box count is inconsistent")
        elif mode in ("empty-lattice-box", "empty-admissible-lattice-box"):
            if recorded != 0:
                raise ArithmeticError("empty box has a nonzero recorded count")
        elif recorded != geometric_count:
            raise ArithmeticError("continuous-cell box count mismatch")
        accepted_box_count += recorded

    result = {
        "schema": "ebch32-parityfanout31x33-b256-dense-cover-partition-audit-v1",
        "status": "EXACT_SUBDIVISION_TREE_AUDIT",
        "claim": {
            "partition_valid": True,
            "all_receipt_leaves_consumed_once": True,
            "recorded_depths_valid": True,
            "accepted_box_counts_valid": True,
            "leaf_count": len(leaves),
            "root_count": len(root_cells),
            "root_six_volumes": root_six_volumes,
            "shell_six_volume": shell_six_volume,
            "root_volume_identity_valid": True,
            "explicit_pulling_triangulation_valid": True,
            "internal_node_count": internal_nodes,
            "visited_node_count": visited_nodes,
            "maximum_leaf_depth": maximum_depth,
            "kind_histogram": kind_histogram,
            "depth_histogram": depth_histogram,
            "sum_recorded_accepted_box_counts_with_overlap": accepted_box_count,
        },
        "scope": (
            "The audit proves that the receipt leaves are exactly the leaves "
            "of the deterministic longest-edge subdivision forest generated "
            "from the explicit three-cell pulling triangulation. The exact "
            "integer volume identity is checked here. The audit does not "
            "prove the numerical witness bounds."
        ),
        "dependencies": [
            {"path": INPUT.name, "sha256": sha256(INPUT)},
            {
                "path": Path(__file__).name,
                "sha256": sha256(Path(__file__).resolve()),
            },
            {
                "path": "diagnose_ebch32_parityfanout_ba_three_band_cover_all_q.py",
                "sha256": sha256(
                    WORKSTREAM
                    / "diagnose_ebch32_parityfanout_ba_three_band_cover_all_q.py"
                ),
            },
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
