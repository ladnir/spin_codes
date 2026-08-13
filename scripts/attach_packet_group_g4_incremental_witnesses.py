#!/usr/bin/env python3
"""Attach new singleton witnesses only where an exact g=4 mesh currently fails."""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import numpy as np

import build_packet_group_g4_support_regular_mesh as mesh_builder
import probe_packet_group_g4_stellar_cover as cover
from reattach_packet_group_g4_support_mixtures import (
    geometry_digest,
    logsumexp2,
    parse_profile,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.mesh.read_text(encoding="utf-8"))
    result = copy.deepcopy(source)
    before_geometry = geometry_digest(source["support_strata"])
    witnesses, new_sources, _anchors, duplicates = cover.load_witnesses(args.atlas)
    evaluator = cover.Evaluator(witnesses)
    improved = 0
    examined = 0
    contributions = []
    for stratum in result["support_strata"]:
        if not stratum.get("feasible", False) or not stratum.get("complete_cover", False):
            continue
        anchors = [parse_profile(profile) for profile in stratum["anchor_profiles"]]
        target = float(stratum["covered_cell_ledger"][0]["target_log2"])
        for cell in stratum["covered_cell_ledger"]:
            if not bool(cell.get("uniform_target_passed", False)):
                examined += 1
                indices = tuple(int(value) for value in cell["anchor_indices"])
                values = np.asarray(
                    [evaluator.evaluate(anchors[index])[1] for index in indices],
                    dtype=np.float64,
                )
                scores = np.max(values, axis=0)
                owner = int(np.argmin(scores))
                score = float(scores[owner])
                if score < float(cell["maximum_vertex_log2"]):
                    improved += 1
                    cell["mixture"] = [
                        {"witness": witnesses[owner].reference, "weight_exact": "1/1"}
                    ]
                    cell["maximum_vertex_log2"] = score
                    cell["vertex_mixture_values_log2"] = values[:, owner].tolist()
                    cell["uniform_target_passed"] = score <= target - 1e-5
                    count = int(cell["integer_profile_count_upper"])
                    cell["cell_local_contribution_log2_upperish"] = math.nextafter(
                        score + math.log2(count), math.inf
                    )
                    cell["mixture_selection_mode"] = "incremental_best_singleton"
                    cell["mixture_candidate_mode"] = "all_new_singletons"
                    cell["eligible_witness_count"] = int(
                        np.sum(np.all(np.isfinite(values), axis=0))
                    )
                    cell["lp_candidate_count_after_cap_and_dominance"] = len(witnesses)
                    cell["mixture_candidate_selection_complete"] = True
                    cell["sparse_basic_solution"] = True
            contributions.append(float(cell["cell_local_contribution_log2_upperish"]))
        covered = stratum["covered_cell_ledger"]
        stratum["uniform_target_passed_cells"] = sum(
            bool(cell["uniform_target_passed"]) for cell in covered
        )
        stratum["cell_local_security_diagnostic_log2"] = logsumexp2(
            [float(cell["cell_local_contribution_log2_upperish"]) for cell in covered]
        )

    after_geometry = geometry_digest(result["support_strata"])
    if after_geometry != before_geometry:
        raise RuntimeError("incremental witness attachment changed exact geometry")
    existing = {row["sha256"]: row for row in result.get("sources", [])}
    for row in new_sources:
        existing[row["sha256"]] = row
    result["sources"] = list(existing.values())
    result["duplicate_atlas_anchor_rows"] = int(
        result.get("duplicate_atlas_anchor_rows", 0)
    ) + int(duplicates)
    result["support_strata_sha256"] = mesh_builder.compact_hash(
        result["support_strata"]
    )
    result["cell_local_security_diagnostic_log2"] = logsumexp2(contributions)
    result["incremental_witness_attachment"] = {
        "examined_currently_nonuniform_cells": examined,
        "improved_cells": improved,
        "new_witnesses": len(witnesses),
        "new_sources": new_sources,
        "geometry_digest": before_geometry,
        "geometry_preserved": True,
        "selection": "all new singleton witnesses evaluated at every failed-cell vertex",
    }
    result["status"] = "G4_SUPPORT_MESH_INCREMENTAL_SINGLETONS_ATTACHED"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    remaining = sum(
        not bool(cell.get("uniform_target_passed", False))
        for stratum in result["support_strata"]
        for cell in stratum.get("covered_cell_ledger", [])
    )
    print(
        json.dumps(
            {
                "examined": examined,
                "improved": improved,
                "remaining_nonuniform_cells": remaining,
                "cell_local_security_diagnostic_log2": result[
                    "cell_local_security_diagnostic_log2"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
