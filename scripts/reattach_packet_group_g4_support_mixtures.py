#!/usr/bin/env python3
"""Re-evaluate mixtures on an existing exact g=4 support mesh.

This is deliberately geometry-free: anchor profiles and terminal root cells are
copied from a previously audited support-mesh artifact.  Only witness values,
mixtures, cell-local count bounds, and aggregate diagnostics are recomputed.
The default ``--candidate-mode column-generation`` repeatedly prices every
eligible finite atlas column while keeping each LP small.  Use mode ``all``
for the monolithic fallback, or ``--candidates-per-vertex -1`` as a shorthand
for column generation with a one-candidate-per-vertex seed.  Mode ``capped``
is only a discovery heuristic and is explicitly labelled as such.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import build_packet_group_g4_support_regular_mesh as mesh_builder
import probe_packet_group_g4_stellar_cover as cover
from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import N


GROUP_BITS = 4
DEFAULT_CANDIDATES_PER_VERTEX = 64

Profile = tuple[Fraction, ...]


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_profile(raw: Any) -> Profile:
    if not isinstance(raw, list) or len(raw) != 5:
        raise ValueError("mixture reattachment: malformed anchor profile")
    if any(isinstance(value, (bool, float)) for value in raw):
        raise ValueError("mixture reattachment: anchor profiles must be exact")
    return tuple(Fraction(str(value)) for value in raw)


def logsumexp2(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return -math.inf
    maximum = max(finite)
    return maximum + math.log2(math.fsum(2.0 ** (value - maximum) for value in finite))


def geometry_digest(strata: list[dict[str, Any]]) -> str:
    """Hash every geometry-bearing field, excluding mixtures and diagnostics."""
    keys = (
        "support_mask",
        "support",
        "dimension",
        "feasible",
        "anchor_profiles",
        "exact_hull_vertices",
        "exact_hull_facets",
        "root_cell_ledger",
        "unrefined_root_cell_ledger",
        "pure_triangle_star_refinement_ledger",
        "exact_geometry",
        "lifting",
        "anchor_profiles_sha256",
        "root_cell_ledger_sha256",
    )
    projection = [
        {key: row[key] for key in keys if key in row}
        for row in strata
    ]
    return mesh_builder.compact_hash(projection)


def reattach_stratum(
    source: dict[str, Any],
    evaluator: cover.Evaluator,
    witnesses: list[cover.Witness],
    candidates_per_vertex: int,
    safety_bits: float,
    candidate_mode: str,
    reduced_cost_tolerance: float,
) -> tuple[dict[str, Any], list[float]]:
    row = copy.deepcopy(source)
    if not row.get("feasible", False):
        row["covered_cell_ledger"] = []
        row["failed_cell_ledger"] = []
        row["complete_cover"] = False
        return row, []
    anchors = [parse_profile(profile) for profile in row["anchor_profiles"]]
    roots = row.get("root_cell_ledger")
    if not isinstance(roots, list) or not roots:
        raise ValueError(
            f"mixture reattachment: support {row.get('support_mask')} has no root cells"
        )
    target = -40.0 - math.log2(profile_count(GROUP_BITS, N))
    global_count = profile_count(GROUP_BITS, N)
    covered = []
    contributions = []
    for geometry in roots:
        cell = cover.Cell(
            str(geometry["cell_id"]),
            tuple(int(value) for value in geometry["anchor_indices"]),
            int(geometry.get("depth", 0)),
        )
        values = np.asarray(
            [evaluator.evaluate(anchors[index])[1] for index in cell.anchor_indices],
            dtype=np.float64,
        )
        eligible = np.flatnonzero(np.all(np.isfinite(values), axis=0))
        result = cover.optimize_values(
            values,
            witnesses,
            target,
            safety_bits,
            candidates_per_vertex,
            1e-10,
            candidate_mode,
            reduced_cost_tolerance,
        )
        bound = cover.cell_bound_record(cell, anchors, result, witnesses, global_count)
        contribution = float(bound["cell_local_contribution_log2_upperish"])
        contributions.append(contribution)
        covered.append(
            {
                **geometry,
                **bound,
                "target_log2": target,
                "uniform_target_passed": result["passes_uniform_target"],
                "eligible_witness_count": int(len(eligible)),
                "lp_candidate_count_after_cap_and_dominance": result["candidate_count"],
                "mixture_selection_mode": result["selection_mode"],
                "mixture_candidate_cap_per_vertex": candidates_per_vertex,
                "mixture_candidate_selection_complete": (
                    candidate_mode in {"all", "column-generation"}
                    or candidates_per_vertex == 0
                ),
                "mixture_candidate_mode": candidate_mode,
                "column_generation_rounds": result["column_generation_rounds"],
                "full_column_scans": result["full_column_scans"],
                "minimum_reduced_cost": result["minimum_reduced_cost"],
                "maximum_pricing_gap": result["maximum_pricing_gap"],
                "reduced_cost_tolerance": result["reduced_cost_tolerance"],
                "sparse_basic_solution": result.get("sparse_basic_solution", True),
            }
        )
    row["covered_cell_ledger"] = covered
    row["failed_cell_ledger"] = []
    row["complete_cover"] = True
    row["uniform_target_passed_cells"] = sum(
        bool(cell["uniform_target_passed"]) for cell in covered
    )
    row["cell_local_security_diagnostic_log2"] = logsumexp2(contributions)
    row["mixture_candidate_cap_per_vertex"] = candidates_per_vertex
    row["mixture_candidate_selection_complete"] = (
        candidate_mode in {"all", "column-generation"}
        or candidates_per_vertex == 0
    )
    row["column_generation_rounds_total"] = sum(
        int(cell["column_generation_rounds"]) for cell in covered
    )
    row["full_column_scans_total"] = sum(
        int(cell["full_column_scans"]) for cell in covered
    )
    reduced_costs = [
        float(cell["minimum_reduced_cost"])
        for cell in covered
        if cell["minimum_reduced_cost"] is not None
    ]
    row["minimum_final_reduced_cost"] = min(reduced_costs, default=None)
    pricing_gaps = [
        float(cell["maximum_pricing_gap"])
        for cell in covered
        if cell["maximum_pricing_gap"] is not None
    ]
    row["maximum_final_pricing_gap"] = max(pricing_gaps, default=None)
    row["all_mixtures_sparse_basic"] = all(
        bool(cell["sparse_basic_solution"]) for cell in covered
    )
    return row, contributions


def reattach(
    payload: dict[str, Any],
    witnesses: list[cover.Witness],
    witness_sources: list[dict[str, str]],
    candidates_per_vertex: int,
    safety_bits: float,
    candidate_mode: str = "capped",
    reduced_cost_tolerance: float = 1e-9,
) -> dict[str, Any]:
    if int(payload.get("group_bits", GROUP_BITS)) != GROUP_BITS:
        raise ValueError("mixture reattachment: input is not g=4")
    strata = payload.get("support_strata")
    if not isinstance(strata, list) or not strata:
        raise ValueError("mixture reattachment: input lacks support_strata")
    source_geometry_sha256 = geometry_digest(strata)
    evaluator = cover.Evaluator(witnesses)
    attached = []
    contributions = []
    for source in strata:
        row, local = reattach_stratum(
            source,
            evaluator,
            witnesses,
            candidates_per_vertex,
            safety_bits,
            candidate_mode,
            reduced_cost_tolerance,
        )
        attached.append(row)
        contributions.extend(local)
    result = copy.deepcopy(payload)
    result["support_strata"] = attached
    attached_geometry_sha256 = geometry_digest(attached)
    if attached_geometry_sha256 != source_geometry_sha256:
        raise RuntimeError("mixture reattachment: geometry changed during reattachment")
    result["geometry_sources"] = payload.get("sources", [])
    result["sources"] = witness_sources
    result["support_strata_sha256"] = mesh_builder.compact_hash(attached)
    masks = [int(row["support_mask"]) for row in attached]
    has_prototype_placeholders = any(
        row.get("reason") == "not_selected_in_prototype_run" for row in attached
    )
    result["complete_support_stratified_cover"] = (
        sorted(masks) == list(range(1, 32))
        and not has_prototype_placeholders
        and all(not row.get("feasible", False) or row.get("complete_cover", False)
                for row in attached)
    )
    result["mixture_reattachment"] = {
        "candidate_policy": (
            "dual_column_generation_with_full_eligible_column_scans"
            if candidate_mode == "column-generation"
            else "all_eligible_then_pareto_prune"
            if candidate_mode == "all" or candidates_per_vertex == 0
            else "union_of_per_vertex_top_k_then_pareto_prune"
        ),
        "candidate_mode": candidate_mode,
        "candidates_per_vertex": candidates_per_vertex,
        "candidate_selection_complete": (
            candidate_mode in {"all", "column-generation"}
            or candidates_per_vertex == 0
        ),
        "safety_bits": safety_bits,
        "geometry_recomputed": False,
        "source_geometry_sha256": source_geometry_sha256,
        "output_geometry_sha256": attached_geometry_sha256,
        "geometry_digest_preserved": True,
        "reduced_cost_tolerance": reduced_cost_tolerance,
        "warning": (
            None
            if candidate_mode in {"all", "column-generation"}
            or candidates_per_vertex == 0
            else "positive per-vertex cap is diagnostic and may exclude an LP-useful witness"
        ),
    }
    result["cell_local_security_diagnostic_log2"] = logsumexp2(contributions)
    result["column_generation_rounds_total"] = sum(
        int(row.get("column_generation_rounds_total", 0)) for row in attached
    )
    result["full_column_scans_total"] = sum(
        int(row.get("full_column_scans_total", 0)) for row in attached
    )
    final_reduced_costs = [
        float(row["minimum_final_reduced_cost"])
        for row in attached
        if row.get("minimum_final_reduced_cost") is not None
    ]
    result["minimum_final_reduced_cost"] = min(final_reduced_costs, default=None)
    final_pricing_gaps = [
        float(row["maximum_final_pricing_gap"])
        for row in attached
        if row.get("maximum_final_pricing_gap") is not None
    ]
    result["maximum_final_pricing_gap"] = max(final_pricing_gaps, default=None)
    result["all_mixtures_sparse_basic"] = all(
        bool(row.get("all_mixtures_sparse_basic", True)) for row in attached
    )
    result["status"] = "G4_SUPPORT_MESH_MIXTURES_REATTACHED"
    return result


def run_self_test() -> None:
    geometry = mesh_builder.build([], {31}, geometry_only=True)
    witnesses = [
        cover.Witness(
            "synthetic", -1.0e9, (0.0,) * 5, frozenset(range(5)), True
        )
    ]
    result = reattach(geometry, witnesses, [], 0, 1e-5, "all")
    full = next(row for row in result["support_strata"] if row["support_mask"] == 31)
    if (
        not full["complete_cover"]
        or len(full["covered_cell_ledger"]) != len(full["root_cell_ledger"])
        or not result["mixture_reattachment"]["candidate_selection_complete"]
        or any(not row["mixture"] for row in full["covered_cell_ledger"])
    ):
        raise SystemExit("mixture reattachment self-test failed")
    # Per-vertex top-k is not an LP-safe reduction.  The balanced witness is
    # second at both vertices but beats every mixture of the two vertex winners.
    cap_values = np.asarray([[0.0, 100.0, 49.0], [100.0, 0.0, 49.0]])
    cap_witnesses = [
        cover.Witness(f"cap_{index}", 0.0, (0.0,) * 5, frozenset(range(5)), True)
        for index in range(3)
    ]
    capped = cover.optimize_values(cap_values, cap_witnesses, -100.0, 0.0, 1, 1e-12)
    complete = cover.optimize_values(
        cap_values, cap_witnesses, -100.0, 0.0, 0, 1e-12, "all"
    )
    generated = cover.optimize_values(
        cap_values,
        cap_witnesses,
        -100.0,
        0.0,
        1,
        1e-12,
        "column-generation",
        1e-10,
    )
    if not (
        abs(capped["score"] - 50.0) < 1e-8
        and abs(complete["score"] - 49.0) < 1e-8
        and abs(generated["score"] - complete["score"]) < 1e-8
        and generated["full_column_scans"] >= 1
        and generated["minimum_reduced_cost"] >= -1e-10
        and generated["maximum_pricing_gap"] <= 1e-10
    ):
        raise SystemExit("mixture reattachment self-test did not expose unsafe cap")
    rng = np.random.default_rng(20260812)
    random_cases = 24
    for case in range(random_cases):
        vertices = 2 + case % 4
        columns = 6 + case
        values = rng.normal(size=(vertices, columns)) * 20.0
        fake = [
            cover.Witness(
                f"random_{case}_{index}",
                0.0,
                (0.0,) * 5,
                frozenset(range(5)),
                True,
            )
            for index in range(columns)
        ]
        all_columns = cover.optimize_values(
            values, fake, -100.0, 0.0, 0, 1e-12, "all"
        )
        column_generated = cover.optimize_values(
            values,
            fake,
            -100.0,
            0.0,
            1,
            1e-12,
            "column-generation",
            1e-10,
        )
        if (
            abs(all_columns["score"] - column_generated["score"]) > 1e-7
            or column_generated["minimum_reduced_cost"] < -1e-10
            or column_generated["maximum_pricing_gap"] > 1e-10
            or not column_generated["sparse_basic_solution"]
        ):
            raise SystemExit(
                "mixture reattachment self-test: column generation disagrees "
                f"on random case {case}: all={all_columns['score']} "
                f"cg={column_generated['score']} rc={column_generated['minimum_reduced_cost']}"
            )
    print(f"self_test_cells={len(full['covered_cell_ledger'])}")
    print("self_test_cap_counterexample=capped_50_complete_49")
    print(f"self_test_random_all_column_equivalence={random_cases}")
    print("status=PASS_G4_SUPPORT_MIXTURE_REATTACHMENT_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path)
    parser.add_argument("--atlas", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--candidates-per-vertex", type=int, default=DEFAULT_CANDIDATES_PER_VERTEX
    )
    parser.add_argument(
        "--candidate-mode",
        choices=("capped", "all", "column-generation"),
        default="column-generation",
    )
    parser.add_argument("--reduced-cost-tolerance", type=float, default=1e-9)
    parser.add_argument("--safety-bits", type=float, default=1e-5)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.mesh is None or not args.atlas or args.output is None:
        raise SystemExit("mixture reattachment: --mesh, --atlas, and --output are required")
    candidate_mode = args.candidate_mode
    candidates_per_vertex = args.candidates_per_vertex
    if candidates_per_vertex == -1:
        candidate_mode = "column-generation"
        candidates_per_vertex = 1
    if candidates_per_vertex < 0 or args.safety_bits < 0 or args.reduced_cost_tolerance <= 0:
        raise SystemExit("mixture reattachment: invalid candidate cap or safety margin")
    payload = json.loads(args.mesh.read_text(encoding="utf-8"))
    witnesses, sources, _anchors, duplicate_rows = cover.load_witnesses(args.atlas)
    result = reattach(
        payload,
        witnesses,
        sources,
        candidates_per_vertex,
        args.safety_bits,
        candidate_mode,
        args.reduced_cost_tolerance,
    )
    result["mixture_reattachment"]["source_mesh"] = {
        "path": str(args.mesh.resolve()),
        "sha256": file_sha256(args.mesh),
    }
    result["duplicate_atlas_anchor_rows"] = duplicate_rows
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    cells = sum(len(row.get("covered_cell_ledger", [])) for row in result["support_strata"])
    print(f"reattached_cells={cells}")
    print(f"candidate_selection_complete={candidate_mode in {'all', 'column-generation'} or candidates_per_vertex == 0}")
    print("status=G4_SUPPORT_MESH_MIXTURES_REATTACHED")


if __name__ == "__main__":
    main()
