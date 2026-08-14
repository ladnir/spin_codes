#!/usr/bin/env python3
"""Replay fixed adaptive leaves with additional supplementary witnesses.

This tool never reconstructs or changes split geometry.  It verifies the
stored active-leaf vertex payloads and exact counts, then recomputes rational
minimax mixtures over a larger bound witness bank.  Binary64 scores still
require independent outward replay.
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

for _variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_variable] = "1"

import numpy as np

import evaluate_packet_group_g8_cumulative_boxes as h2
import plan_packet_group_g8_lowdim_root_cover as minimax
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import produce_packet_group_g8_highdim_dominance_shards as dominance
import run_packet_group_g8_adaptive_cumulative_engine as engine


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CHECKPOINT = ROOT / "out" / "g8_full_support_adaptive_cumulative_checkpoint.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_full_support_adaptive_atlas_replay.json"
CHECKPOINT_SHA256 = "9d3407b2d4a5e1b6463aca831b2bdec79da7d0a14341ad0a83f67ff9cbf8da0c"
SCHEMA = "permute-conv.packet-group-g8-adaptive-atlas-replay.v1"


def log2_sum(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log2(
        math.fsum(math.exp2(value - maximum) for value in values)
    )


def load_checkpoint(path: Path) -> tuple[dict[str, Any], str]:
    digest = coordinate.sha256_path(path)
    if digest != CHECKPOINT_SHA256:
        raise ValueError(
            f"adaptive checkpoint digest changed: {digest}; expected {CHECKPOINT_SHA256}"
        )
    checkpoint = json.loads(path.read_bytes())
    if checkpoint.get("schema") != engine.SCHEMA:
        raise ValueError("unexpected adaptive checkpoint schema")
    if checkpoint.get("source_bindings", {}).get("manifest_sha256") != h2.EXPECTED_MANIFEST_SHA256:
        raise ValueError("adaptive checkpoint binds another manifest")
    return checkpoint, digest


def load_extended_bank(
    manifest_path: Path,
    atlas_path: Path,
    original_supplementary_path: Path,
    added_paths: list[Path],
) -> tuple[coordinate.AtlasBank, list[dict[str, Any]], dict[str, Any]]:
    base_bank, original_sources, bindings = h2.build_bank(
        manifest_path, atlas_path, original_supplementary_path
    )
    additions = []
    seen_digests = {source["sha256"] for source in original_sources}
    for path in added_paths:
        addition = dominance.load_supplementary_atlas(
            path,
            bindings["manifest_sha256"],
            bindings["base_atlas_sha256"],
        )
        source = dict(addition[3])
        if source["sha256"] in seen_digests:
            raise ValueError(f"added atlas duplicates an existing source: {path}")
        seen_digests.add(source["sha256"])
        source["manifest_bound"] = False
        additions.append(addition)
    if not additions:
        raise ValueError("at least one added supplementary atlas is required")
    bank = dominance.combined_bank(base_bank, additions)
    sources = [*original_sources, *(addition[3] for addition in additions)]
    return bank, sources, {
        **bindings,
        "original_supplementary_sha256": bindings["supplementary_atlas_sha256"],
        "added_supplementary_sha256": [addition[3]["sha256"] for addition in additions],
    }


def active_leaves(checkpoint: dict[str, Any]) -> list[dict[str, Any]]:
    leaves = [
        node
        for node in checkpoint["nodes"].values()
        if node.get("state") == "ACTIVE_LEAF"
    ]
    leaves.sort(key=lambda node: node["node_id"])
    if not leaves:
        raise ValueError("adaptive checkpoint has no active leaves")
    return leaves


def validate_leaf_payload(node: dict[str, Any]) -> tuple[tuple[int, ...], ...]:
    profiles = tuple(tuple(map(int, row)) for row in node["exact_vertices"])
    if (
        not profiles
        or len(profiles) != int(node["exact_vertex_count"])
        or len(set(profiles)) != len(profiles)
        or any(len(profile) != 9 for profile in profiles)
        or any(any(value < 1 for value in profile) for profile in profiles)
        or any(sum(profile) != 262144 for profile in profiles)
    ):
        raise ValueError(f"malformed stored leaf vertices: {node['node_id']}")
    digest = coordinate.sha256_bytes(
        coordinate.canonical_bytes([list(profile) for profile in profiles])
    )
    if digest != node["exact_vertex_sha256"]:
        raise ValueError(f"stored leaf vertex digest changed: {node['node_id']}")
    if int(node["exact_count"]) <= 0:
        raise ValueError(f"stored leaf count is not positive: {node['node_id']}")
    return profiles


def replay_leaf(
    node: dict[str, Any], bank: coordinate.AtlasBank, denominator: int
) -> dict[str, Any]:
    profiles = validate_leaf_payload(node)
    matrix = np.asarray(profiles, dtype=np.float64)
    normalizations = coordinate.diagnostic_normalization(matrix, 262144)
    values = (
        bank.constants[:, None] - bank.charges @ matrix.T - normalizations[None, :]
    ).T
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"extended bank produced a non-finite score: {node['node_id']}")
    mixed = minimax.optimize_minimax_mixture(
        values,
        minimax.UNIFORM_TARGET_LOG2,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-9,
        weight_tolerance=1e-12,
        denominator=denominator,
    )
    upper = float(mixed["score"])
    count = int(node["exact_count"])
    old_upper = float(node["candidate_upper_log2"])
    old_contribution = float(node["contribution_log2"])
    contribution = math.log2(count) + upper
    worst_index = int(np.argmax(np.asarray(mixed["vertex_values"])))
    return {
        "node_id": node["node_id"],
        "h2_cell": int(node["h2_cell"]),
        "depth": int(node["depth"]),
        "exact_count": node["exact_count"],
        "exact_vertex_count": len(profiles),
        "exact_vertex_sha256": node["exact_vertex_sha256"],
        "old_candidate_upper_log2": old_upper,
        "candidate_upper_log2": upper,
        "upper_improvement_bits": old_upper - upper,
        "old_contribution_log2": old_contribution,
        "contribution_log2": contribution,
        "contribution_improvement_bits": old_contribution - contribution,
        "worst_vertex": list(profiles[worst_index]),
        "selector": h2.selector(mixed, bank),
        "component_count": len(mixed["indices"]),
        "candidate_status": (
            "PASS_BINARY64_REQUIRES_OUTWARD_REPLAY"
            if upper <= minimax.UNIFORM_TARGET_LOG2
            else "FAIL_BINARY64"
        ),
        "diagnostic": {
            "binary64_lp_upper_log2": float(mixed["binary64_lp_score"]),
            "rationalization_loss_bits": float(mixed["rationalization_loss_bits"]),
            "column_generation_rounds": int(mixed["column_generation_rounds"]),
            "full_column_scans": int(mixed["full_column_scans"]),
        },
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    checkpoint, checkpoint_digest = load_checkpoint(args.checkpoint)
    bank, sources, bindings = load_extended_bank(
        args.manifest,
        args.atlas,
        args.original_supplementary,
        args.added_supplementary,
    )
    checkpoint_sources = {
        source["sha256"]
        for source in checkpoint["source_bindings"]["witness_sources"]
    }
    original_sources = {source["sha256"] for source in sources[:2]}
    if checkpoint_sources != original_sources:
        raise ValueError("checkpoint witness-source bindings changed")
    leaves = active_leaves(checkpoint)
    exact_count = sum(int(node["exact_count"]) for node in leaves)
    if exact_count != math.comb(262143, 8):
        raise ValueError("stored active-leaf counts do not recover the exact root census")
    old_contributions = [float(node["contribution_log2"]) for node in leaves]
    replayed = [replay_leaf(node, bank, args.mixture_denominator) for node in leaves]
    contributions = [float(node["contribution_log2"]) for node in replayed]
    old_aggregate = log2_sum(old_contributions)
    aggregate = log2_sum(contributions)
    worst = max(replayed, key=lambda node: (node["contribution_log2"], node["node_id"]))
    return {
        "schema": SCHEMA,
        "status": "DIAGNOSTIC_GEOMETRY_FREE_ATLAS_REPLAY_REQUIRES_OUTWARD_REPLAY",
        "checkpoint_binding": {
            "path": str(args.checkpoint),
            "sha256": checkpoint_digest,
            "schema": checkpoint["schema"],
            "completed_waves": int(checkpoint["completed_waves"]),
        },
        **bindings,
        "witness_sources": sources,
        "mixture_discovery": {
            "method": "binary64-safe-full-column-scan-minimax-v1",
            "rational_denominator": args.mixture_denominator,
            "seed_per_vertex": 1,
            "reduced_cost_tolerance": 1e-9,
            "weight_tolerance": 1e-12,
        },
        "geometry_policy": (
            "geometry-free: reload stored active-leaf vertices, digests, and exact "
            "counts; do not reconstruct hulls or create splits"
        ),
        "leaves": replayed,
        "summary": {
            "active_leaf_count": len(leaves),
            "exact_root_count": str(exact_count),
            "eligible_witness_count": len(bank.references),
            "added_witness_count": len(bank.references) - 542,
            "old_aggregate_log2_union_diagnostic": old_aggregate,
            "aggregate_log2_union_diagnostic": aggregate,
            "aggregate_improvement_bits": old_aggregate - aggregate,
            "worst_leaf_id": worst["node_id"],
            "worst_contribution_log2": worst["contribution_log2"],
            "binary64_pass_leaf_count": sum(
                node["candidate_upper_log2"] <= minimax.UNIFORM_TARGET_LOG2
                for node in replayed
            ),
            "strictly_improved_leaf_count": sum(
                node["contribution_improvement_bits"] > 0.0 for node in replayed
            ),
            "runtime_seconds": time.perf_counter() - started,
        },
        "scope_limit": (
            "Stored geometry and counts are replay inputs, not reconstructed claims. "
            "Binary64 mixtures require independent outward replay. This artifact is "
            "not a certificate, completion claim, or probability proof."
        ),
    }


def run_self_test() -> None:
    old = np.asarray([[1.0, -3.0], [-3.0, 1.0]], dtype=np.float64)
    old_result = minimax.optimize_minimax_mixture(
        old,
        0.0,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-10,
        weight_tolerance=1e-12,
        denominator=1024,
    )
    extended = np.column_stack((old, np.asarray([-2.0, -2.0])))
    new_result = minimax.optimize_minimax_mixture(
        extended,
        0.0,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-10,
        weight_tolerance=1e-12,
        denominator=1024,
    )
    payload = [[1] * 8 + [262136]]
    digest = coordinate.sha256_bytes(coordinate.canonical_bytes(payload))
    node = {
        "node_id": "synthetic",
        "exact_vertices": payload,
        "exact_vertex_count": 1,
        "exact_vertex_sha256": digest,
        "exact_count": "1",
    }
    if not (
        validate_leaf_payload(node) == (tuple(payload[0]),)
        and old_result["score"] == -1.0
        and new_result["score"] == -2.0
        and new_result["indices"] == (2,)
        and abs(log2_sum([3.0, 3.0]) - 4.0) < 1e-12
    ):
        raise SystemExit("geometry-free atlas replay self-test failed")
    print("stored_vertex_digest_validation=PASS")
    print("synthetic_old_mixture_score=-1.0")
    print("synthetic_added_witness_score=-2.0")
    print("expanded_aggregate=PASS")
    print("new_split_count=0")
    print("status=PASS_G8_GEOMETRY_FREE_ATLAS_REPLAY_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=h2.DEFAULT_ATLAS)
    parser.add_argument(
        "--original-supplementary", type=Path, default=h2.DEFAULT_SUPPLEMENTARY
    )
    parser.add_argument(
        "--added-supplementary", type=Path, action="append", default=[]
    )
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mixture-denominator", type=int, default=1 << 30)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if not args.added_supplementary or args.mixture_denominator < 1:
        parser.error("at least one --added-supplementary and positive denominator are required")
    report = build(args)
    digest = coordinate.write_canonical(args.output, report)
    summary = report["summary"]
    print(f"active_leaves={summary['active_leaf_count']}")
    print(f"eligible_witnesses={summary['eligible_witness_count']}")
    print(f"aggregate_log2={summary['aggregate_log2_union_diagnostic']:.12f}")
    print(f"aggregate_improvement_bits={summary['aggregate_improvement_bits']:.12f}")
    print(f"runtime_seconds={summary['runtime_seconds']:.3f}")
    print(f"output={args.output}")
    print(f"sha256={digest}")
    print("status=DIAGNOSTIC_GEOMETRY_FREE_ATLAS_REPLAY_REQUIRES_OUTWARD_REPLAY")


if __name__ == "__main__":
    main()
