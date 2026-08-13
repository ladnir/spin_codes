#!/usr/bin/env python3
"""Extract exact pointwise BSP failures and integer profiles for witness tuning.

The BSP partition is replayed exactly.  A residual is emitted only for a leaf
vertex whose current frozen-atlas lower envelope is above the cell's uniform
target (including its discovery safety allowance).  Exact rational geometry is
kept for audit; the integer profile is a separate, explicitly labelled tuning
anchor and is never used as proof evidence for the rational point.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import certify_packet_group_g4_cell_bsp as certify
import probe_packet_group_g4_cell_bsp as bsp
import probe_packet_group_g4_stellar_cover as cover


SCHEMA = "packet-group-g4-bsp-pointwise-residuals-v1"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def round_tuning_profile(profile: bsp.Profile) -> tuple[int, ...]:
    """Round inside the exact full-support, fixed-sum, weight>=21 domain."""

    floors = [value.numerator // value.denominator for value in profile]
    remainder = cover.M - sum(floors)
    order = sorted(
        range(len(profile)),
        key=lambda index: (
            -(profile[index] - floors[index]),
            -index,
        ),
    )
    rounded = floors[:]
    for index in order[:remainder]:
        rounded[index] += 1
    if sum(rounded) != cover.M or any(value < 1 for value in rounded):
        raise RuntimeError("g4 BSP residual extractor: full-support rounding failed")

    weight = sum(index * value for index, value in enumerate(rounded))
    while weight < 21:
        donor = next(
            (
                index
                for index in range(len(rounded))
                if rounded[index] > 1 and index < len(rounded) - 1
            ),
            None,
        )
        if donor is None:
            raise RuntimeError("g4 BSP residual extractor: weight repair failed")
        rounded[donor] -= 1
        rounded[-1] += 1
        weight += len(rounded) - 1 - donor
    return tuple(rounded)


def replay_cell(cell: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, tuple[bsp.Profile, ...]]]:
    root_vertices = tuple(bsp.parse_profile(row) for row in cell["root_anchor_profiles"])
    constraints = bsp.simplex_constraints(root_vertices)
    if bsp.enumerate_vertices(constraints) != tuple(sorted(root_vertices)):
        raise RuntimeError("g4 BSP residual extractor: root reconstruction failed")
    return certify.replay_tree(cell, constraints)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bsp-batch", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    batch = json.loads(args.bsp_batch.read_text(encoding="utf-8"))
    schema = batch.get("schema")
    if schema not in (
        "packet-group-g4-cell-dominance-bsp-batch-v1",
        bsp.BSP_SCHEMA,
    ):
        raise ValueError("g4 BSP residual extractor: unsupported batch schema")
    witnesses, sources, _anchors, duplicates = cover.load_witnesses(args.atlas)
    evaluator = cover.Evaluator(witnesses)
    if schema == bsp.BSP_SCHEMA:
        cells = [batch]
        safety = float(batch["safety_bits"])
    else:
        cells = batch["cells"]
        safety = float(batch["discovery"]["safety_bits"])

    by_profile: dict[tuple[int, ...], dict[str, Any]] = {}
    exact_failures = 0
    failure_leaves = 0
    for cell in cells:
        leaves, vertices_by_leaf = replay_cell(cell)
        target = float(cell["uniform_target_log2"])
        threshold = target - safety
        cell_contribution = float(
            cell.get("diagnostic", {}).get(
                "bsp_cell_union_contribution_log2", float("-inf")
            )
        )
        for leaf in leaves:
            if leaf.get("status") != "pointwise_atlas_failure":
                continue
            failure_leaves += 1
            leaf_id = str(leaf["node_id"])
            for vertex_index, profile in enumerate(vertices_by_leaf[leaf_id]):
                values = evaluator.evaluate(profile)[1]
                owner = min(
                    range(len(witnesses)),
                    key=lambda index: (float(values[index]), witnesses[index].reference),
                )
                envelope = float(values[owner])
                if not math.isfinite(envelope) or envelope <= threshold:
                    continue
                exact_failures += 1
                rounded = round_tuning_profile(profile)
                row = {
                    "name": "bsp_" + "_".join(map(str, rounded)),
                    "profile": list(rounded),
                    "source_cell_id": str(cell["cell_id"]),
                    "source_leaf_id": leaf_id,
                    "source_leaf_vertex_index": vertex_index,
                    "exact_rational_profile": [fraction_text(value) for value in profile],
                    "diagnostic_current_envelope_log2": envelope,
                    "diagnostic_uniform_threshold_log2": threshold,
                    "diagnostic_gap_bits": envelope - threshold,
                    "diagnostic_best_witness": witnesses[owner].reference,
                    "tuning_profile_semantics": (
                        "largest-remainder full-support integer approximation; discovery only"
                    ),
                }
                origin = {
                    "source_cell_id": str(cell["cell_id"]),
                    "source_leaf_id": leaf_id,
                    "source_leaf_vertex_index": vertex_index,
                    "source_cell_union_contribution_log2": cell_contribution,
                    "exact_rational_profile": row["exact_rational_profile"],
                    "diagnostic_gap_bits": row["diagnostic_gap_bits"],
                }
                prior = by_profile.get(rounded)
                if prior is None:
                    row["source_occurrences"] = [origin]
                    by_profile[rounded] = row
                    continue

                occurrence_key = (
                    origin["source_cell_id"],
                    origin["source_leaf_id"],
                    origin["source_leaf_vertex_index"],
                    tuple(origin["exact_rational_profile"]),
                )
                prior_keys = {
                    (
                        item["source_cell_id"],
                        item["source_leaf_id"],
                        item["source_leaf_vertex_index"],
                        tuple(item["exact_rational_profile"]),
                    )
                    for item in prior["source_occurrences"]
                }
                occurrences = prior["source_occurrences"]
                if occurrence_key not in prior_keys:
                    occurrences.append(origin)
                if row["diagnostic_gap_bits"] > prior["diagnostic_gap_bits"]:
                    row["source_occurrences"] = occurrences
                    by_profile[rounded] = row

    residuals = sorted(
        by_profile.values(),
        key=lambda row: (-float(row["diagnostic_gap_bits"]), row["profile"]),
    )
    for row in residuals:
        row["source_occurrences"].sort(
            key=lambda item: (
                -float(item["source_cell_union_contribution_log2"]),
                item["source_cell_id"],
                item["source_leaf_id"],
                item["source_leaf_vertex_index"],
            )
        )
        row["source_cell_ids"] = sorted(
            {item["source_cell_id"] for item in row["source_occurrences"]}
        )
        row["diagnostic_max_source_cell_union_contribution_log2"] = max(
            float(item["source_cell_union_contribution_log2"])
            for item in row["source_occurrences"]
        )
    report = {
        "schema": SCHEMA,
        "status": "DIAGNOSTIC_G4_BSP_POINTWISE_RESIDUAL_LEDGER",
        "group_bits": 4,
        "source_bsp_batch": {
            "path": str(args.bsp_batch.resolve()),
            "sha256": file_sha256(args.bsp_batch),
        },
        "sources": sources,
        "duplicate_atlas_anchor_rows": duplicates,
        "exact_pointwise_failure_vertex_occurrences": exact_failures,
        "pointwise_failure_leaves": failure_leaves,
        "unique_integer_tuning_profiles": len(residuals),
        "uncovered_integer_residuals": residuals,
        "proof_scope": (
            "The exact rational profiles identify atlas failures. Rounded integer profiles "
            "are witness-tuning inputs only and do not certify those rational vertices."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "exact_pointwise_failure_vertex_occurrences": exact_failures,
                "pointwise_failure_leaves": failure_leaves,
                "unique_integer_tuning_profiles": len(residuals),
                "largest_gap_bits": residuals[0]["diagnostic_gap_bits"] if residuals else None,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
