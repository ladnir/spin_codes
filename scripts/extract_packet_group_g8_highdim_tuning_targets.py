#!/usr/bin/env python3
"""Extract unique tuning targets from unresolved high-dimensional g=8 shards.

This is a read-only bridge from adaptive coverage discovery to witness tuning.
It validates frozen input and batch digests, deduplicates exact worst-vertex
profiles, and preserves every batch/shard/parent/node contribution provenance.
It does not optimize, evaluate, or certify a witness.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "G8_SUPPORT_MANIFEST.json"
DEFAULT_ATLAS = ROOT / "out" / "g8_support_seed_atlas.json"
EXPECTED_MANIFEST_SHA256 = (
    "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
)
EXPECTED_ATLAS_SHA256 = (
    "c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e"
)
MANIFEST_SCHEMA = "packet-group-g8-support-manifest-v1"
ATLAS_SCHEMA = "permute-conv.packet-group-g8-support-seed-atlas.v1"
BATCH_SCHEMA = "permute-conv.packet-group-g8-highdim-adaptive-batch.v1"
SHARD_SCHEMA = "packet-group-g8-support-shard-v1"
OUTPUT_SCHEMA = "permute-conv.packet-group-g8-tuning-target-catalogue.v1"
CLASSES = 9


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def write_canonical(path: Path, value: Any) -> str:
    payload = canonical_bytes(value) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return sha256_bytes(payload)


def parse_mask(value: Any) -> int:
    if not isinstance(value, str) or value != value.lower() or not value.startswith("0x"):
        raise ValueError("support mask is not canonical hexadecimal")
    parsed = int(value, 16)
    if value != f"0x{parsed:03x}" or not 0 < parsed < (1 << CLASSES):
        raise ValueError("support mask is outside nine classes")
    return parsed


def validate_frozen_inputs(manifest_path: Path, atlas_path: Path):
    manifest_digest = sha256_path(manifest_path)
    atlas_digest = sha256_path(atlas_path)
    if manifest_digest != EXPECTED_MANIFEST_SHA256:
        raise ValueError(
            f"manifest digest changed: {manifest_digest}; expected {EXPECTED_MANIFEST_SHA256}"
        )
    if atlas_digest != EXPECTED_ATLAS_SHA256:
        raise ValueError(
            f"atlas digest changed: {atlas_digest}; expected {EXPECTED_ATLAS_SHA256}"
        )
    manifest = json.loads(manifest_path.read_bytes())
    atlas = json.loads(atlas_path.read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA or atlas.get("schema") != ATLAS_SCHEMA:
        raise ValueError("frozen manifest or atlas schema changed")
    sources = manifest.get("witness_sources", [])
    if len(sources) != 1 or sources[0].get("sha256") != atlas_digest:
        raise ValueError("manifest does not bind the supplied witness atlas")
    if int(manifest["parameters"]["M"]) != 262144:
        raise ValueError("manifest profile mass changed")
    return manifest, atlas, sources[0], manifest_digest, atlas_digest


def parent_map(shard: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for node in shard.get("nodes", []):
        if node.get("state") != "SPLIT":
            continue
        parent = str(node.get("node_id"))
        split = node.get("split")
        if not isinstance(split, dict):
            raise ValueError(f"split node {parent!r} lacks split data")
        coefficients = split.get("coefficients")
        threshold = split.get("threshold")
        if not isinstance(coefficients, list) or len(coefficients) != CLASSES:
            raise ValueError(f"split node {parent!r} has malformed coefficients")
        for side in ("left", "right"):
            child = node.get(side)
            if not isinstance(child, str) or child in result:
                raise ValueError("shard has a missing child or duplicate parent")
            result[child] = {
                "parent_node_id": parent,
                "child_side": side,
                "split_coefficients": [int(value) for value in coefficients],
                "split_threshold": str(threshold),
            }
    return result


def parse_worst_profile(value: Any, total: int) -> tuple[int, ...]:
    if not isinstance(value, list) or len(value) != CLASSES:
        raise ValueError("unresolved worst_vertex must contain nine integers")
    if any(isinstance(item, bool) for item in value):
        raise ValueError("worst_vertex contains a Boolean")
    profile = tuple(int(item) for item in value)
    if any(count < 0 for count in profile) or sum(profile) != total:
        raise ValueError("worst_vertex is not a nonnegative profile of mass M")
    return profile


def validate_terminal_ownership(
    profile: tuple[int, ...], shard: dict[str, Any], expected_node: str
) -> None:
    nodes = {str(node.get("node_id")): node for node in shard.get("nodes", [])}
    node_id = str(shard.get("root_node", "r"))
    visited: set[str] = set()
    while True:
        if node_id in visited or node_id not in nodes:
            raise ValueError("shard ownership path has a cycle or missing node")
        visited.add(node_id)
        node = nodes[node_id]
        if node.get("state") != "SPLIT":
            if node_id != expected_node:
                raise ValueError(
                    f"worst_vertex belongs to {node_id!r}, not claimed node {expected_node!r}"
                )
            return
        split = node.get("split", {})
        coefficients = split.get("coefficients")
        if not isinstance(coefficients, list) or len(coefficients) != CLASSES:
            raise ValueError("ownership split has malformed coefficients")
        value = sum(int(left) * right for left, right in zip(coefficients, profile))
        threshold = int(split.get("threshold"))
        node_id = str(node.get("left") if value <= threshold else node.get("right"))


def occurrence_from_node(
    batch_path: Path,
    batch_digest: str,
    shard_path: Path,
    shard_digest: str,
    shard: dict[str, Any],
    node: dict[str, Any],
    parents: dict[str, dict[str, Any]],
    total: int,
) -> tuple[tuple[int, ...], dict[str, Any]] | None:
    if node.get("state") != "UNRESOLVED":
        return None
    diagnostic = node.get("diagnostic")
    if not isinstance(diagnostic, dict) or diagnostic.get("worst_vertex") is None:
        return None
    profile = parse_worst_profile(diagnostic["worst_vertex"], total)
    support_mask = parse_mask(shard.get("support_mask"))
    exact_mask = sum((1 << index) for index, count in enumerate(profile) if count)
    if exact_mask != support_mask:
        raise ValueError("worst_vertex escaped its exact support shard")
    contribution = diagnostic.get("collapsed_candidate_contribution_log2")
    upper = diagnostic.get("candidate_upper_log2")
    residual = diagnostic.get("positive_residual_bits")
    if not all(
        isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
        for value in (contribution, upper, residual)
    ):
        raise ValueError("unresolved node lacks finite diagnostic contribution fields")
    node_id = str(node.get("node_id"))
    validate_terminal_ownership(profile, shard, node_id)
    count = node.get("count")
    if not isinstance(count, dict) or "value" not in count:
        raise ValueError("unresolved node lacks its owned-count record")
    selector = diagnostic.get("candidate_selector")
    provenance = {
        "batch_index": {
            "path": str(batch_path),
            "sha256": batch_digest,
        },
        "shard": {
            "path": str(shard_path),
            "sha256": shard_digest,
            "support_mask": shard["support_mask"],
            "active_classes": shard.get("active_classes"),
            "root_count": str(shard.get("root_count")),
        },
        "node_id": node_id,
        "parent": parents.get(node_id),
        "owned_count": {
            "kind": count.get("kind"),
            "method": count.get("method"),
            "value": str(count["value"]),
        },
        "candidate_selector": selector,
        "candidate_upper_log2": float(upper),
        "positive_residual_bits": float(residual),
        "collapsed_candidate_contribution_log2": float(contribution),
        "unresolved_reason": diagnostic.get("unresolved_reason"),
        "diagnostic_arithmetic": diagnostic.get("arithmetic"),
    }
    return profile, provenance


def collect_batches(
    batch_paths: list[Path], manifest_digest: str, atlas_digest: str, total: int
) -> tuple[dict[tuple[int, ...], list[dict[str, Any]]], list[dict[str, Any]], int]:
    collected: dict[tuple[int, ...], list[dict[str, Any]]] = defaultdict(list)
    input_records = []
    unresolved_without_vertex = 0
    seen_shards: set[str] = set()
    for batch_path in batch_paths:
        batch_digest = sha256_path(batch_path)
        batch = json.loads(batch_path.read_bytes())
        if batch.get("schema") != BATCH_SCHEMA:
            raise ValueError(f"unexpected batch schema in {batch_path}")
        if batch.get("manifest_sha256") != manifest_digest:
            raise ValueError(f"batch {batch_path} uses another manifest")
        if batch.get("atlas_sha256") != atlas_digest:
            raise ValueError(f"batch {batch_path} uses another witness atlas")
        shard_records = batch.get("shards")
        if not isinstance(shard_records, list):
            raise ValueError(f"batch {batch_path} has no shard list")
        input_record = {
            "path": str(batch_path),
            "sha256": batch_digest,
            "shards": [],
        }
        for record in shard_records:
            locator = Path(str(record.get("path", "")))
            shard_path = locator if locator.is_absolute() else batch_path.parent / locator
            shard_digest = sha256_path(shard_path)
            if shard_digest != record.get("sha256"):
                raise ValueError(f"shard digest mismatch: {shard_path}")
            if shard_digest in seen_shards:
                raise ValueError(f"duplicate shard input: {shard_path}")
            seen_shards.add(shard_digest)
            shard = json.loads(shard_path.read_bytes())
            if shard.get("schema") != SHARD_SCHEMA:
                raise ValueError(f"unexpected shard schema in {shard_path}")
            if shard.get("manifest_sha256") != manifest_digest:
                raise ValueError(f"shard {shard_path} uses another manifest")
            if shard.get("state") != "INCOMPLETE":
                raise ValueError("target extractor accepts only explicitly incomplete shards")
            parents = parent_map(shard)
            input_record["shards"].append(
                {
                    "path": str(shard_path),
                    "sha256": shard_digest,
                    "support_mask": shard.get("support_mask"),
                }
            )
            for node in shard.get("nodes", []):
                if node.get("state") == "UNRESOLVED" and not isinstance(
                    node.get("diagnostic", {}).get("worst_vertex"), list
                ):
                    unresolved_without_vertex += 1
                occurrence = occurrence_from_node(
                    batch_path,
                    batch_digest,
                    shard_path,
                    shard_digest,
                    shard,
                    node,
                    parents,
                    total,
                )
                if occurrence is not None:
                    profile, provenance = occurrence
                    collected[profile].append(provenance)
        input_records.append(input_record)
    return collected, input_records, unresolved_without_vertex


def build_targets(
    collected: dict[tuple[int, ...], list[dict[str, Any]]], maximum: int
) -> tuple[list[dict[str, Any]], int]:
    prepared = []
    for profile, occurrences in collected.items():
        occurrences.sort(
            key=lambda row: (
                -row["collapsed_candidate_contribution_log2"],
                row["shard"]["support_mask"],
                row["node_id"],
                row["shard"]["sha256"],
            )
        )
        maximum_contribution = max(
            row["collapsed_candidate_contribution_log2"] for row in occurrences
        )
        maximum_residual = max(row["positive_residual_bits"] for row in occurrences)
        profile_list = list(profile)
        profile_digest = sha256_bytes(canonical_bytes(profile_list))
        prepared.append(
            {
                "profile": profile_list,
                "profile_sha256": profile_digest,
                "support_mask": f"0x{sum(1 << index for index, count in enumerate(profile) if count):03x}",
                "physical_weight": str(
                    sum(index * count for index, count in enumerate(profile))
                ),
                "occurrence_count": len(occurrences),
                "maximum_collapsed_candidate_contribution_log2": maximum_contribution,
                "maximum_positive_residual_bits": maximum_residual,
                "provenance": occurrences,
            }
        )
    prepared.sort(
        key=lambda row: (
            -row["maximum_collapsed_candidate_contribution_log2"],
            -row["maximum_positive_residual_bits"],
            row["profile"],
        )
    )
    unique_before_limit = len(prepared)
    if maximum:
        prepared = prepared[:maximum]
    for rank, row in enumerate(prepared, 1):
        row["rank"] = rank
        row["target_id"] = f"target-{rank:06d}-{row['profile_sha256'][:12]}"
    return prepared, unique_before_limit


def synthetic_self_test() -> None:
    profile = [262136] + [1] * 8
    shard = {
        "support_mask": "0x1ff",
        "active_classes": list(range(9)),
        "root_count": "99",
        "nodes": [
            {
                "node_id": "r",
                "state": "SPLIT",
                "left": "r0",
                "right": "r1",
                "split": {"coefficients": [1] + [0] * 8, "threshold": "100"},
            },
            {
                "node_id": "r0",
                "state": "EMPTY",
            },
            {
                "node_id": "r1",
                "state": "UNRESOLVED",
                "count": {"kind": "upper", "method": "coordinate-box-ie-v1", "value": "59"},
                "diagnostic": {
                    "worst_vertex": profile,
                    "candidate_selector": {"kind": "witness", "witness": {"row": 2}},
                    "candidate_upper_log2": 110.0,
                    "positive_residual_bits": 210.0,
                    "collapsed_candidate_contribution_log2": 115.0,
                    "unresolved_reason": "MAX_DEPTH_REACHED",
                    "arithmetic": "binary64-discovery-only",
                },
            },
        ],
    }
    repeated_shard = {
        "support_mask": "0x1ff",
        "active_classes": list(range(9)),
        "root_count": "99",
        "root_node": "r",
        "nodes": [
            {
                "node_id": "r",
                "state": "UNRESOLVED",
                "count": {
                    "kind": "upper",
                    "method": "coordinate-box-ie-v1",
                    "value": "99",
                },
                "diagnostic": {
                    "worst_vertex": profile,
                    "candidate_selector": {"kind": "witness", "witness": {"row": 1}},
                    "candidate_upper_log2": 100.0,
                    "positive_residual_bits": 200.0,
                    "collapsed_candidate_contribution_log2": 105.0,
                    "unresolved_reason": "NODE_BUDGET_EXHAUSTED",
                    "arithmetic": "binary64-discovery-only",
                },
            }
        ],
    }
    parents = parent_map(shard)
    collected: dict[tuple[int, ...], list[dict[str, Any]]] = defaultdict(list)
    for node in shard["nodes"]:
        occurrence = occurrence_from_node(
            Path("synthetic_batch.json"),
            "a" * 64,
            Path("synthetic_shard.json"),
            "b" * 64,
            shard,
            node,
            parents,
            262144,
        )
        if occurrence:
            target_profile, provenance = occurrence
            collected[target_profile].append(provenance)
    repeated = occurrence_from_node(
        Path("synthetic_batch_2.json"),
        "c" * 64,
        Path("synthetic_shard_2.json"),
        "d" * 64,
        repeated_shard,
        repeated_shard["nodes"][0],
        {},
        262144,
    )
    if repeated:
        target_profile, provenance = repeated
        collected[target_profile].append(provenance)
    targets, unique = build_targets(collected, 0)
    if (
        unique != 1
        or len(targets) != 1
        or targets[0]["occurrence_count"] != 2
        or targets[0]["maximum_collapsed_candidate_contribution_log2"] != 115.0
        or targets[0]["provenance"][0]["parent"]["child_side"] != "right"
    ):
        raise AssertionError("synthetic tuning-target deduplication failed")
    print("synthetic_unique_targets=1")
    print("synthetic_occurrences=2")
    print("synthetic_max_contribution_log2=115.0")
    print("status=SYNTHETIC_TUNING_TARGET_SMOKE_PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-index", type=Path, action="append", default=[])
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--max-targets",
        type=int,
        default=0,
        help="retain the highest-ranked prefix; zero retains every unique target",
    )
    parser.add_argument("--synthetic-self-test", action="store_true")
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if not args.batch_index or args.output is None:
        parser.error("--batch-index and --output are required")
    if args.max_targets < 0:
        parser.error("--max-targets must be nonnegative")
    try:
        manifest, _atlas, source, manifest_digest, atlas_digest = validate_frozen_inputs(
            args.manifest, args.atlas
        )
        total = int(manifest["parameters"]["M"])
        collected, inputs, missing = collect_batches(
            args.batch_index, manifest_digest, atlas_digest, total
        )
        targets, unique_before_limit = build_targets(collected, args.max_targets)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise SystemExit(f"g=8 tuning-target extractor: {error}") from error
    output = {
        "schema": OUTPUT_SCHEMA,
        "status": "DIAGNOSTIC_TUNING_TARGETS_ONLY",
        "scope": (
            "exact profiles deduplicated from unresolved binary64 worst vertices; "
            "no witness optimization, evaluation, or coverage claim"
        ),
        "manifest_sha256": manifest_digest,
        "atlas_source_id": source["source_id"],
        "atlas_sha256": atlas_digest,
        "canonical_json": "sorted-compact-json-v1",
        "inputs": inputs,
        "unresolved_nodes_with_worst_vertex": sum(
            row["occurrence_count"]
            for row in build_targets(collected, 0)[0]
        ),
        "unresolved_nodes_without_worst_vertex": missing,
        "unique_targets_before_limit": unique_before_limit,
        "targets_emitted": len(targets),
        "max_targets": args.max_targets,
        "ranking": (
            "descending maximum collapsed candidate contribution, then descending "
            "positive residual, then lexicographic exact profile"
        ),
        "targets": targets,
    }
    digest = write_canonical(args.output, output)
    print(f"unique_targets={unique_before_limit}")
    print(f"targets_emitted={len(targets)}")
    print(f"output={args.output}")
    print(f"output_sha256={digest}")
    print("status=DIAGNOSTIC_TUNING_TARGETS_ONLY")


if __name__ == "__main__":
    main()
