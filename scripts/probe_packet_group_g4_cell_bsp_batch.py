#!/usr/bin/env python3
"""Run the exact-dyadic g=4 per-cell BSP discovery on a contribution queue."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
from pathlib import Path
from typing import Any

import probe_packet_group_g4_cell_bsp as bsp
import probe_packet_group_g4_stellar_cover as cover


BATCH_SCHEMA = "packet-group-g4-cell-dominance-bsp-batch-v1"
_WITNESSES: list[cover.Witness] | None = None


def initialize_worker(atlas_paths: list[str]) -> None:
    global _WITNESSES
    _WITNESSES = cover.load_witnesses([Path(path) for path in atlas_paths])[0]


def build_cell(task: dict[str, Any]) -> dict[str, Any]:
    if _WITNESSES is None:
        raise RuntimeError("g4 BSP batch worker was not initialized")
    bsp.configure_support(int(task["support_mask"]))
    root_vertices = tuple(bsp.parse_profile(row) for row in task["root_anchor_profiles"])
    constraints = bsp.simplex_constraints(root_vertices)
    if bsp.enumerate_vertices(constraints) != tuple(sorted(root_vertices)):
        raise RuntimeError("g4 BSP batch: root reconstruction failed")
    builder = bsp.Builder(
        _WITNESSES,
        float(task["target_log2"]),
        float(task["safety_bits"]),
        int(task["max_depth"]),
        int(task["top_witnesses_per_vertex"]),
    )
    root_node = builder.build_node(constraints, 0)
    maximum = max(builder.leaf_scores)
    count = int(task["integer_profile_count_upper"])
    contribution = math.nextafter(maximum + math.log2(count), math.inf)
    statuses = {
        status: builder.leaf_statuses.count(status)
        for status in sorted(set(builder.leaf_statuses))
    }
    return {
        "support_mask": int(task["support_mask"]),
        "cell_id": task["cell_id"],
        "root_anchor_indices": task["root_anchor_indices"],
        "root_anchor_profiles": task["root_anchor_profiles"],
        "root_integer_profile_count_upper": count,
        "uniform_target_log2": float(task["target_log2"]),
        "root_node": root_node,
        "nodes": builder.nodes,
        "complete_exact_partition_claimed": True,
        "diagnostic": {
            "root_maximum_vertex_log2": float(task["root_maximum_vertex_log2"]),
            "root_union_contribution_log2": float(task["root_union_contribution_log2"]),
            "bsp_cell_maximum_log2": maximum,
            "bsp_cell_union_contribution_log2": contribution,
            "improvement_bits": float(task["root_maximum_vertex_log2"]) - maximum,
            "leaf_count": len(builder.leaf_scores),
            "leaf_status_counts": statuses,
            "all_leaves_pass_uniform_target": all(
                status == "covered" for status in builder.leaf_statuses
            ),
        },
    }


def logsumexp2(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log2(math.fsum(2.0 ** (value - maximum) for value in values))


def write_report(path: Path, report: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument("--support-mask", type=lambda value: int(value, 0), action="append")
    parser.add_argument("--all-lower-supports", action="store_true")
    parser.add_argument("--top-cells", type=int, default=32)
    parser.add_argument("--max-depth", type=int, default=5)
    parser.add_argument("--top-witnesses-per-vertex", type=int, default=2)
    parser.add_argument("--safety-bits", type=float, default=1e-5)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        args.top_cells < 1
        or args.max_depth < 0
        or args.top_witnesses_per_vertex < 1
        or args.safety_bits < 0
        or args.workers < 1
        or args.workers > 8
    ):
        parser.error("invalid BSP batch limits")
    if args.all_lower_supports and args.support_mask:
        parser.error("choose either --all-lower-supports or --support-mask")

    ledger = json.loads(args.ledger.read_text(encoding="utf-8"))
    masks = (
        list(range(1, 31))
        if args.all_lower_supports
        else (args.support_mask if args.support_mask else [31])
    )
    strata = {
        int(row["support_mask"]): row
        for row in ledger["support_strata"]
        if int(row["support_mask"]) in masks and row.get("feasible")
    }
    if set(strata) != {mask for mask in masks if mask != 1}:
        # Support {0} is the unique infeasible nonempty support.
        raise ValueError("g4 BSP batch: requested support strata are missing")
    ranked = sorted(
        (
            (float(row["cell_local_contribution_log2_upperish"]), mask, row)
            for mask, stratum in strata.items()
            for row in stratum["covered_cell_ledger"]
        ),
        reverse=True,
        key=lambda item: (item[0], item[1], str(item[2].get("cell_id"))),
    )
    selected = ranked[: args.top_cells]
    tasks = []
    for _term, mask, row in selected:
        anchors = strata[mask]["anchor_profiles"]
        indices = list(map(int, row["anchor_indices"]))
        tasks.append(
            {
                "support_mask": mask,
                "cell_id": str(row["cell_id"]),
                "root_anchor_indices": indices,
                "root_anchor_profiles": [anchors[index] for index in indices],
                "integer_profile_count_upper": int(row["integer_profile_count_upper"]),
                "target_log2": float(row["target_log2"]),
                "root_maximum_vertex_log2": float(row["maximum_vertex_log2"]),
                "root_union_contribution_log2": float(
                    row["cell_local_contribution_log2_upperish"]
                ),
                "safety_bits": args.safety_bits,
                "max_depth": args.max_depth,
                "top_witnesses_per_vertex": args.top_witnesses_per_vertex,
            }
        )

    witnesses, sources, _atlas_anchors, duplicates = cover.load_witnesses(args.atlas)
    if args.workers == 1:
        global _WITNESSES
        _WITNESSES = witnesses
        results = [build_cell(task) for task in tasks]
    else:
        results = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=initialize_worker,
            initargs=([str(path.resolve()) for path in args.atlas],),
        ) as pool:
            future_to_cell = {pool.submit(build_cell, task): task["cell_id"] for task in tasks}
            for future in concurrent.futures.as_completed(future_to_cell):
                result = future.result()
                results.append(result)
                print(
                    f"completed={result['cell_id']} "
                    f"term={result['diagnostic']['bsp_cell_union_contribution_log2']:.6f} "
                    f"leaves={result['diagnostic']['leaf_count']}",
                    flush=True,
                )
    results.sort(key=lambda row: (int(row["support_mask"]), row["cell_id"]))

    replacements = {
        (int(row["support_mask"]), row["cell_id"]): float(
            row["diagnostic"]["bsp_cell_union_contribution_log2"]
        )
        for row in results
    }
    updated_terms = [
        replacements.get((mask, str(row["cell_id"])), term)
        for term, mask, row in ranked
    ]
    report = {
        "schema": BATCH_SCHEMA,
        "status": "DIAGNOSTIC_G4_CELL_DOMINANCE_BSP_BATCH",
        "group_bits": 4,
        "source_ledger": {
            "path": str(args.ledger.resolve()),
            "sha256": bsp.file_sha256(args.ledger),
        },
        "sources": sources,
        "duplicate_atlas_anchor_rows": duplicates,
        "support_mask": masks[0] if len(masks) == 1 else masks,
        "discovery": {
            "top_cells": args.top_cells,
            "max_depth": args.max_depth,
            "top_witnesses_per_vertex": args.top_witnesses_per_vertex,
            "safety_bits": args.safety_bits,
            "workers": args.workers,
            "host": os.environ.get("COMPUTERNAME", os.uname().nodename if hasattr(os, "uname") else "unknown"),
        },
        "cells": results,
        "diagnostic": {
            "source_cell_local_union_log2": logsumexp2([term for term, _mask, _row in ranked]),
            "updated_cell_local_union_log2": logsumexp2(updated_terms),
            "processed_cells": len(results),
            "all_processed_cells_pass_uniform_target": all(
                row["diagnostic"]["all_leaves_pass_uniform_target"] for row in results
            ),
            "largest_updated_term_log2": max(updated_terms),
            "largest_unprocessed_term_log2": max(
                (
                    term
                    for term, mask, row in ranked
                    if (mask, str(row["cell_id"])) not in replacements
                ),
                default=float("-inf"),
            ),
        },
    }
    write_report(args.output, report)
    print(json.dumps(report["diagnostic"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
