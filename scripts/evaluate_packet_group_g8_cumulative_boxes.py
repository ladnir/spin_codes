#!/usr/bin/env python3
"""Evaluate fixed witnesses and mixtures on exact full-support cumulative cells.

The cell geometry and counts are exact.  Witness scores, minimax discovery,
and rational-mixture evaluation use binary64 arithmetic.  Consequently, this
script emits a diagnostic artifact and never claims an outward proof.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

# Keep BLAS and HiGHS use serialized.  This lane must not create a concurrent
# benchmark accidentally through a threaded numerical runtime.
for _variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_variable] = "1"

import numpy as np

import plan_packet_group_g8_lowdim_root_cover as minimax
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import produce_packet_group_g8_highdim_dominance_shards as dominance
import prototype_packet_group_g8_cumulative_boxes as geometry


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "G8_SUPPORT_MANIFEST.json"
DEFAULT_ATLAS = ROOT / "out" / "g8_support_seed_atlas.json"
DEFAULT_SUPPLEMENTARY = ROOT / "out" / "g8_highdim_supplementary_final32.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_full_support_cumulative_h2_diagnostic.json"
SCHEMA = "permute-conv.packet-group-g8-cumulative-box-diagnostic.v1"
EXPECTED_MANIFEST_SHA256 = (
    "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
)
EXPECTED_ATLAS_SHA256 = (
    "c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e"
)
EXPECTED_SUPPLEMENTARY_SHA256 = (
    "d50af6e55d5bc1d0aa0a9de65025684fdfe9989cce0014c808a059c2ca356d6c"
)
FULL_SUPPORT = tuple(range(9))
LEVELS = 2


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_vertex(row: list[str]) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in row)


def selector(result: dict[str, Any], bank: coordinate.AtlasBank) -> dict[str, Any]:
    components = [
        {
            "weight": fraction_text(weight),
            "witness": bank.references[index],
        }
        for index, weight in zip(result["indices"], result["weights"])
    ]
    if len(components) == 1 and components[0]["weight"] == "1/1":
        return {"kind": "witness", "witness": components[0]["witness"]}
    return {"kind": "mixture", "components": components}


def singleton_selector(index: int, bank: coordinate.AtlasBank) -> dict[str, Any]:
    return {"kind": "witness", "witness": bank.references[index]}


def log2_sum(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log2(math.fsum(math.exp2(value - maximum) for value in values))


def worst_record(cells: list[dict[str, Any]], method: str, field: str) -> dict[str, Any]:
    index = max(range(len(cells)), key=lambda item: cells[item][method][field])
    cell = cells[index]
    return {
        "cell": index,
        "bin_indices": cell["bin_indices"],
        "exact_count": cell["exact_count"],
        "candidate_upper_log2": cell[method]["candidate_upper_log2"],
        "contribution_log2": cell[method]["contribution_log2"],
        "uniform_target_gap_bits": cell[method]["uniform_target_gap_bits"],
    }


def build_bank(
    manifest_path: Path, atlas_path: Path, supplementary_path: Path
) -> tuple[coordinate.AtlasBank, list[dict[str, Any]], dict[str, str]]:
    (
        _manifest,
        atlas,
        source,
        bindings,
        manifest_digest,
        atlas_digest,
    ) = coordinate.load_frozen_inputs(manifest_path, atlas_path)
    supplementary_digest = coordinate.sha256_path(supplementary_path)
    if manifest_digest != EXPECTED_MANIFEST_SHA256:
        raise ValueError("manifest digest changed")
    if atlas_digest != EXPECTED_ATLAS_SHA256:
        raise ValueError("base atlas digest changed")
    if supplementary_digest != EXPECTED_SUPPLEMENTARY_SHA256:
        raise ValueError("supplementary atlas digest changed")
    base = coordinate.build_atlas_bank(atlas, source, bindings)
    extra = dominance.load_supplementary_atlas(
        supplementary_path, manifest_digest, atlas_digest
    )
    combined = dominance.combined_bank(base, [extra])
    base_source = {
        "source_id": source["source_id"],
        "path": str(atlas_path),
        "sha256": atlas_digest,
        "schema": source["schema"],
        "row_count": len(bindings),
        "manifest_bound": True,
    }
    supplementary_source = dict(extra[3])
    supplementary_source["manifest_bound"] = False
    return combined, [base_source, supplementary_source], {
        "manifest_sha256": manifest_digest,
        "base_atlas_sha256": atlas_digest,
        "supplementary_atlas_sha256": supplementary_digest,
    }


def evaluate(
    manifest_path: Path,
    atlas_path: Path,
    supplementary_path: Path,
    *,
    denominator: int,
    seed_per_vertex: int,
    reduced_cost_tolerance: float,
    weight_tolerance: float,
) -> dict[str, Any]:
    started = time.perf_counter()
    bank, sources, bindings = build_bank(
        manifest_path, atlas_path, supplementary_path
    )
    structural = geometry.build(FULL_SUPPORT, LEVELS, 10000)
    cells = structural["cells"]
    if (
        len(cells) != 165
        or structural["summary"]["maximum_vertices_per_cell"] != 81
        or structural["summary"]["total_vertex_incidences"] != 6435
    ):
        raise RuntimeError("full-support h=2 cumulative geometry changed")

    # Score all vertex incidences in one dense pass.  The same matrix then
    # serves each small cell LP without repeated affine evaluation.
    vertices = [parse_vertex(row) for cell in cells for row in cell["vertices"]]
    profiles = np.asarray(vertices, dtype=np.float64)
    normalizations = coordinate.diagnostic_normalization(profiles, geometry.MASS)
    all_values = (
        bank.constants[:, None] - bank.charges @ profiles.T - normalizations[None, :]
    ).T
    if not np.all(np.isfinite(all_values)):
        raise RuntimeError("combined witness bank produced a non-finite score")

    target = minimax.UNIFORM_TARGET_LOG2
    evaluated = []
    offset = 0
    total_rounds = 0
    total_scans = 0
    for cell_index, cell in enumerate(cells):
        vertex_count = len(cell["vertices"])
        values = all_values[offset : offset + vertex_count]
        offset += vertex_count
        single = minimax.best_singleton(values)
        mixed = minimax.optimize_minimax_mixture(
            values,
            target,
            seed_per_vertex=seed_per_vertex,
            reduced_cost_tolerance=reduced_cost_tolerance,
            weight_tolerance=weight_tolerance,
            denominator=denominator,
        )
        count = int(cell["exact_count"])
        count_log2 = math.log2(count)
        total_rounds += int(mixed["column_generation_rounds"])
        total_scans += int(mixed["full_column_scans"])
        single_score = float(single["score"])
        mixed_score = float(mixed["score"])
        evaluated.append(
            {
                **cell,
                "cell": cell_index,
                "singleton": {
                    "candidate_status": (
                        "PASS_BINARY64_REQUIRES_OUTWARD_REPLAY"
                        if single_score <= target
                        else "FAIL_BINARY64"
                    ),
                    "candidate_upper_log2": single_score,
                    "uniform_target_gap_bits": single_score - target,
                    "contribution_log2": count_log2 + single_score,
                    "selector": singleton_selector(int(single["index"]), bank),
                },
                "mixture": {
                    "candidate_status": (
                        "PASS_BINARY64_REQUIRES_OUTWARD_REPLAY"
                        if mixed_score <= target
                        else "FAIL_BINARY64"
                    ),
                    "candidate_upper_log2": mixed_score,
                    "binary64_lp_upper_log2": float(mixed["binary64_lp_score"]),
                    "rationalization_loss_bits": float(
                        mixed["rationalization_loss_bits"]
                    ),
                    "uniform_target_gap_bits": mixed_score - target,
                    "contribution_log2": count_log2 + mixed_score,
                    "component_count": len(mixed["indices"]),
                    "selector": selector(mixed, bank),
                    "column_generation_rounds": int(
                        mixed["column_generation_rounds"]
                    ),
                    "full_column_scans": int(mixed["full_column_scans"]),
                    "minimum_reduced_cost": mixed["minimum_reduced_cost"],
                    "maximum_pricing_gap": mixed["maximum_pricing_gap"],
                },
            }
        )
    if offset != len(vertices):
        raise RuntimeError("cell slices did not consume all vertex incidences")

    singleton_contributions = [row["singleton"]["contribution_log2"] for row in evaluated]
    mixture_contributions = [row["mixture"]["contribution_log2"] for row in evaluated]
    singleton_aggregate = log2_sum(singleton_contributions)
    mixture_aggregate = log2_sum(mixture_contributions)
    exact_root_count = sum(int(row["exact_count"]) for row in evaluated)
    if exact_root_count != math.comb(geometry.MASS - 1, 8):
        raise RuntimeError("diagnostic cells lost exact full-support profiles")

    summary = {
        **structural["summary"],
        "eligible_witnesses": len(bank.references),
        "base_witnesses": len(bank.references) - 32,
        "supplementary_witnesses": 32,
        "singleton_pass_cells": sum(
            row["singleton"]["candidate_upper_log2"] <= target for row in evaluated
        ),
        "mixture_pass_cells": sum(
            row["mixture"]["candidate_upper_log2"] <= target for row in evaluated
        ),
        "mixture_strictly_improves_cells": sum(
            row["mixture"]["candidate_upper_log2"]
            < row["singleton"]["candidate_upper_log2"] - 1e-12
            for row in evaluated
        ),
        "singleton_worst_bound": worst_record(evaluated, "singleton", "candidate_upper_log2"),
        "mixture_worst_bound": worst_record(evaluated, "mixture", "candidate_upper_log2"),
        "singleton_worst_contribution": worst_record(evaluated, "singleton", "contribution_log2"),
        "mixture_worst_contribution": worst_record(evaluated, "mixture", "contribution_log2"),
        "singleton_aggregate_log2_union_diagnostic": singleton_aggregate,
        "mixture_aggregate_log2_union_diagnostic": mixture_aggregate,
        "singleton_aggregate_gap_from_2^-40_bits": singleton_aggregate + 40.0,
        "mixture_aggregate_gap_from_2^-40_bits": mixture_aggregate + 40.0,
        "column_generation_rounds": total_rounds,
        "full_column_scans": total_scans,
        "runtime_seconds": time.perf_counter() - started,
    }
    return {
        "schema": SCHEMA,
        "status": "DIAGNOSTIC_BINARY64_EXACT_GEOMETRY_REQUIRES_OUTWARD_REPLAY",
        **bindings,
        "witness_sources": sources,
        "support_mask": structural["support_mask"],
        "active_classes": list(FULL_SUPPORT),
        "dimension": 8,
        "levels": LEVELS,
        "bins": 4,
        "exact_root_count": str(exact_root_count),
        "uniform_profile_target_log2": target,
        "mixture_discovery": {
            "method": "binary64-safe-full-column-scan-minimax-v1",
            "rational_denominator": denominator,
            "seed_per_vertex": seed_per_vertex,
            "reduced_cost_tolerance": reduced_cost_tolerance,
            "weight_tolerance": weight_tolerance,
        },
        "geometry": {
            "method": "balanced-shifted-cumulative-boxes-h2-v1",
            "integer_bin_intervals": structural["integer_bin_intervals"],
            "ownership": structural["ownership"],
            "canonical_prefix_splits": structural["canonical_prefix_splits"],
            "count_method": "run-product-multiset-count-v1",
        },
        "summary": summary,
        "cells": evaluated,
        "scope_limit": (
            "Exact ownership, vertices, and counts only. Binary64 affine scores and "
            "mixture discovery require independent outward replay. This artifact is "
            "not a shard, certificate, completion claim, or probability proof."
        ),
    }


def run_self_test() -> None:
    values = np.asarray([[1.0, -3.0], [-3.0, 1.0]], dtype=np.float64)
    single = minimax.best_singleton(values)
    mixed = minimax.optimize_minimax_mixture(
        values,
        0.0,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-10,
        weight_tolerance=1e-12,
        denominator=1024,
    )
    fake_bank = coordinate.AtlasBank(
        np.zeros(2),
        np.zeros((2, 9)),
        (
            {"source_id": "x", "source_sha256": "0" * 64, "row": 0, "row_sha256": "1" * 64},
            {"source_id": "x", "source_sha256": "0" * 64, "row": 1, "row_sha256": "2" * 64},
        ),
    )
    selected = selector(mixed, fake_bank)
    if not (
        single["score"] == 1.0
        and mixed["score"] == -1.0
        and mixed["weights"] == (Fraction(1, 2), Fraction(1, 2))
        and selected["kind"] == "mixture"
        and [row["weight"] for row in selected["components"]] == ["1/2", "1/2"]
    ):
        raise SystemExit("cumulative diagnostic minimax self-test failed")
    structural = geometry.build(FULL_SUPPORT, LEVELS, 10000)
    if (
        len(structural["cells"]) != 165
        or sum(int(row["exact_count"]) for row in structural["cells"])
        != math.comb(geometry.MASS - 1, 8)
        or structural["summary"]["total_vertex_incidences"] != 6435
    ):
        raise SystemExit("cumulative diagnostic geometry self-test failed")
    print("synthetic_singleton_score=1.0")
    print("synthetic_rational_mixture=1/2,1/2")
    print("full_support_h2_cells=165")
    print("full_support_h2_vertex_incidences=6435")
    print("status=PASS_G8_CUMULATIVE_BOX_DIAGNOSTIC_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--supplementary", type=Path, default=DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mixture-denominator", type=int, default=1 << 30)
    parser.add_argument("--seed-per-vertex", type=int, default=1)
    parser.add_argument("--reduced-cost-tolerance", type=float, default=1e-9)
    parser.add_argument("--weight-tolerance", type=float, default=1e-12)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if (
        args.mixture_denominator < 1
        or args.seed_per_vertex < 1
        or args.reduced_cost_tolerance <= 0.0
        or args.weight_tolerance <= 0.0
    ):
        parser.error("invalid minimax discovery parameter")
    report = evaluate(
        args.manifest,
        args.atlas,
        args.supplementary,
        denominator=args.mixture_denominator,
        seed_per_vertex=args.seed_per_vertex,
        reduced_cost_tolerance=args.reduced_cost_tolerance,
        weight_tolerance=args.weight_tolerance,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    digest = coordinate.write_canonical(args.output, report)
    summary = report["summary"]
    print(
        f"cells={summary['nonempty_cells']} witnesses={summary['eligible_witnesses']} "
        f"singleton_pass={summary['singleton_pass_cells']} "
        f"mixture_pass={summary['mixture_pass_cells']}"
    )
    print(
        "mixture_worst_bound_log2="
        f"{summary['mixture_worst_bound']['candidate_upper_log2']:.12f} "
        "mixture_aggregate_log2="
        f"{summary['mixture_aggregate_log2_union_diagnostic']:.12f}"
    )
    print(f"runtime_seconds={summary['runtime_seconds']:.3f}")
    print(f"output={args.output}")
    print(f"sha256={digest}")
    print("status=DIAGNOSTIC_BINARY64_EXACT_GEOMETRY_REQUIRES_OUTWARD_REPLAY")


if __name__ == "__main__":
    main()
