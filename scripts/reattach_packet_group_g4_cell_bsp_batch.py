#!/usr/bin/env python3
"""Re-optimize leaf owners on a frozen exact g=4 BSP batch.

The split geometry is copied unchanged.  Only leaf witness ownership and the
associated binary64 discovery diagnostics are recomputed against an expanded
frozen atlas.  This makes atlas enrichment monotone for a fixed partition and
avoids losing a useful old partition to a different greedy BSP rebuild.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

import extract_packet_group_g4_bsp_failures as extract
import probe_packet_group_g4_cell_bsp_batch as batch
import probe_packet_group_g4_stellar_cover as cover


def reattach_cell(
    cell: dict[str, Any],
    evaluator: cover.Evaluator,
    witnesses: list[cover.Witness],
    safety_bits: float,
) -> dict[str, Any]:
    result = copy.deepcopy(cell)
    leaves, vertices_by_leaf = extract.replay_cell(result)
    node_by_id = {str(row["node_id"]): row for row in result["nodes"]}
    threshold = float(result["uniform_target_log2"]) - safety_bits
    scores: list[float] = []
    statuses: list[str] = []

    for old_leaf in leaves:
        identifier = str(old_leaf["node_id"])
        leaf = node_by_id[identifier]
        values = np.asarray(
            [evaluator.evaluate(profile)[1] for profile in vertices_by_leaf[identifier]],
            dtype=np.float64,
        )
        owner_scores = np.max(values, axis=0)
        owner = int(np.argmin(owner_scores))
        owner_score = float(owner_scores[owner])
        envelope_maximum = float(np.max(np.min(values, axis=1)))
        if owner_score <= threshold:
            status = "covered"
        elif envelope_maximum > threshold:
            status = "pointwise_atlas_failure"
        else:
            status = "reattached_interpolation_gap"
        leaf["owner_witness"] = witnesses[owner].reference
        leaf["status"] = status
        leaf["diagnostic_owner_maximum_log2"] = owner_score
        leaf["diagnostic_envelope_maximum_at_vertices_log2"] = envelope_maximum
        scores.append(owner_score)
        statuses.append(status)

    maximum = max(scores)
    count = int(result["root_integer_profile_count_upper"])
    contribution = math.nextafter(maximum + math.log2(count), math.inf)
    status_counts = Counter(statuses)
    diagnostic = result["diagnostic"]
    diagnostic.update(
        {
            "bsp_cell_maximum_log2": maximum,
            "bsp_cell_union_contribution_log2": contribution,
            "improvement_bits": float(diagnostic["root_maximum_vertex_log2"])
            - maximum,
            "leaf_count": len(scores),
            "leaf_status_counts": dict(sorted(status_counts.items())),
            "all_leaves_pass_uniform_target": all(
                status == "covered" for status in statuses
            ),
        }
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bsp-batch", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.bsp_batch.read_text(encoding="utf-8"))
    if source.get("schema") != batch.BATCH_SCHEMA:
        raise SystemExit("g4 BSP reattachment: unsupported batch schema")
    witnesses, sources, _anchors, duplicates = cover.load_witnesses(args.atlas)
    evaluator = cover.Evaluator(witnesses)
    safety_bits = float(source["discovery"]["safety_bits"])
    cells = [
        reattach_cell(cell, evaluator, witnesses, safety_bits)
        for cell in source["cells"]
    ]

    report = copy.deepcopy(source)
    report["cells"] = cells
    report["sources"] = sources
    report["duplicate_atlas_anchor_rows"] = duplicates
    report["status"] = "DIAGNOSTIC_G4_FROZEN_BSP_LEAF_REATTACHMENT"
    report["discovery"]["reattached_from"] = {
        "path": str(args.bsp_batch.resolve()),
        "sha256": batch.bsp.file_sha256(args.bsp_batch),
        "semantics": "exact split geometry frozen; leaf owners re-optimized",
    }
    terms = [
        float(cell["diagnostic"]["bsp_cell_union_contribution_log2"])
        for cell in cells
    ]
    report["diagnostic"].update(
        {
            "updated_cell_local_union_log2": batch.logsumexp2(terms),
            "largest_updated_term_log2": max(terms),
            "all_processed_cells_pass_uniform_target": all(
                cell["diagnostic"]["all_leaves_pass_uniform_target"]
                for cell in cells
            ),
        }
    )
    batch.write_report(args.output, report)
    print(
        json.dumps(
            {
                "cells": len(cells),
                "witnesses": len(witnesses),
                "largest_updated_term_log2": max(terms),
                "all_processed_cells_pass_uniform_target": report["diagnostic"][
                    "all_processed_cells_pass_uniform_target"
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
