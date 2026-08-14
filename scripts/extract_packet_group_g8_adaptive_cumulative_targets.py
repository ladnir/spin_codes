#!/usr/bin/env python3
"""Extract canonical tuning targets from an adaptive cumulative checkpoint.

The extractor is read-only with respect to its checkpoint and witness inputs.
It ranks active leaves by their stored contribution, deduplicates exact stored
worst vertices, and emits the tuning-target catalogue consumed by
``tune_packet_group_g8_supplementary_atlas.py``.  It performs no optimization.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import produce_packet_group_g8_highdim_dominance_shards as dominance


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CHECKPOINT = ROOT / "out" / "g8_full_support_adaptive_cumulative_checkpoint.json"
DEFAULT_SUPPLEMENTARY = ROOT / "out" / "g8_highdim_supplementary_final32.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_adaptive_cumulative_tuning_targets_top16.json"
CHECKPOINT_SCHEMA = "permute-conv.packet-group-g8-adaptive-cumulative-checkpoint.v1"
SUPPLEMENTARY_SCHEMA = "permute-conv.packet-group-g8-supplementary-witness-atlas.v1"
OUTPUT_SCHEMA = "permute-conv.packet-group-g8-tuning-target-catalogue.v1"
EXPECTED_MANIFEST_SHA256 = (
    "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
)
EXPECTED_BASE_ATLAS_SHA256 = (
    "c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e"
)
EXPECTED_SUPPLEMENTARY_SHA256 = (
    "d50af6e55d5bc1d0aa0a9de65025684fdfe9989cce0014c808a059c2ca356d6c"
)
CLASSES = 9
NEIGHBORHOOD_RADIUS = 2048
CORE_TARGET_COUNT = 12
TRANSITION_ENDPOINTS = (
    {
        "transition_id": "h2:162/L:y5:left",
        "parent_id": "h2:162/L",
        "coordinate": 5,
        "regime": "left",
        "profile": [196602, 1, 2, 15031, 1, 1, 50504, 1, 1],
    },
    {
        "transition_id": "h2:162/L:y5:right",
        "parent_id": "h2:162/L",
        "coordinate": 5,
        "regime": "right",
        "profile": [196602, 1, 2, 15031, 1, 50504, 1, 1, 1],
    },
    {
        "transition_id": "h2:161/R:y3:left",
        "parent_id": "h2:161/R",
        "coordinate": 3,
        "regime": "left",
        "profile": [196602, 1, 1, 2, 65534, 1, 1, 1, 1],
    },
    {
        "transition_id": "h2:161/R:y3:right",
        "parent_id": "h2:161/R",
        "coordinate": 3,
        "regime": "right",
        "profile": [196602, 1, 1, 65535, 1, 1, 1, 1, 1],
    },
)


def validate_profile(value: Any, total: int) -> tuple[int, ...]:
    if (
        not isinstance(value, list)
        or len(value) != CLASSES
        or any(not isinstance(item, int) or isinstance(item, bool) or item <= 0 for item in value)
    ):
        raise ValueError("adaptive worst vertex must contain nine positive integers")
    profile = tuple(value)
    if sum(profile) != total:
        raise ValueError("adaptive worst vertex has the wrong total mass")
    return profile


def validate_selector(value: Any) -> dict[str, Any]:
    if not isinstance(value, dict) or value.get("kind") not in ("witness", "mixture"):
        raise ValueError("active leaf has a malformed current selector")
    # Canonical serialization catches non-JSON values and non-finite floats.
    coordinate.canonical_bytes(value)
    return value


def transported_mass(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    distance_sum = sum(abs(a - b) for a, b in zip(left, right))
    if distance_sum & 1:
        raise ValueError("equal-mass profiles have a nonintegral transported distance")
    return distance_sum // 2


def logsumexp2(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log2(sum(2.0 ** (value - maximum) for value in values))


def witness_lookup(bank: coordinate.AtlasBank) -> dict[tuple[str, int], tuple[float, np.ndarray]]:
    result = {}
    for index, reference in enumerate(bank.references):
        key = (str(reference["source_id"]), int(reference["row"]))
        if key in result:
            raise ValueError("combined bank contains a duplicate witness reference")
        result[key] = (float(bank.constants[index]), bank.charges[index])
    return result


def reconstruct_selector_maximum(
    leaf: dict[str, Any],
    lookup: dict[tuple[str, int], tuple[float, np.ndarray]],
    total: int,
) -> tuple[tuple[int, ...], float, dict[str, Any]]:
    selector = validate_selector(leaf.get("selector"))
    vertices = tuple(validate_profile(row, total) for row in leaf.get("exact_vertices", []))
    if not vertices or len(vertices) != int(leaf.get("exact_vertex_count", -1)):
        raise ValueError("active leaf exact vertex catalogue is malformed")
    if selector["kind"] == "witness":
        components = [{"weight": "1/1", "witness": selector["witness"]}]
    else:
        components = selector.get("components")
        if not isinstance(components, list) or not components:
            raise ValueError("active leaf mixture is empty")
    constant = 0.0
    charges = np.zeros(CLASSES, dtype=np.float64)
    weight_sum = Fraction(0)
    for component in components:
        weight = Fraction(str(component["weight"]))
        if weight < 0:
            raise ValueError("selector weight is negative")
        reference = component["witness"]
        key = (str(reference["source_id"]), int(reference["row"]))
        affine = lookup.get(key)
        if affine is None:
            raise ValueError(f"selector references unavailable witness {key!r}")
        scalar = float(weight)
        constant += scalar * affine[0]
        charges += scalar * affine[1]
        weight_sum += weight
    if weight_sum != 1:
        raise ValueError("selector weights do not sum exactly to one")
    matrix = np.asarray(vertices, dtype=np.float64)
    scores = constant - charges @ matrix.T - coordinate.diagnostic_normalization(matrix, total)
    maximum = float(np.max(scores))
    maximizers = [vertices[index] for index, value in enumerate(scores) if float(value) == maximum]
    profile = min(maximizers)
    stored = float(leaf.get("candidate_upper_log2"))
    tolerance = max(0.01, float(leaf.get("diagnostic", {}).get("rationalization_loss_bits", 0.0)) + 1e-6)
    if abs(maximum - stored) > tolerance:
        raise ValueError(
            f"selector replay changed at {leaf.get('node_id')}: {maximum} != {stored}"
        )
    return profile, maximum, {
        "stored_candidate_upper_log2": stored,
        "replayed_candidate_upper_log2": maximum,
        "absolute_replay_difference_bits": abs(maximum - stored),
        "maximum_tie_count": len(maximizers),
        "maximum_tie_rule": "lexicographically least exact profile",
    }


def ancestor_path(nodes: dict[str, dict[str, Any]], leaf: dict[str, Any]) -> list[dict[str, Any]]:
    reverse = []
    current = leaf
    visited: set[str] = set()
    while True:
        identifier = str(current.get("node_id", ""))
        if not identifier or identifier in visited:
            raise ValueError("adaptive checkpoint parent chain is cyclic or malformed")
        visited.add(identifier)
        step = {
            "node_id": identifier,
            "depth": int(current.get("depth", 0)),
            "branch": current.get("branch"),
        }
        if current.get("state") == "SPLIT":
            step["split"] = current.get("split")
        reverse.append(step)
        parent = current.get("parent_id")
        if parent is None:
            break
        current = nodes.get(str(parent))
        if current is None:
            raise ValueError(f"adaptive checkpoint omits parent {parent!r}")
    reverse.reverse()
    expected_depths = list(range(len(reverse)))
    if [int(row["depth"]) for row in reverse] != expected_depths:
        raise ValueError("adaptive checkpoint path depths are not contiguous")
    return reverse


def occurrence_from_leaf(
    checkpoint_digest: str,
    nodes: dict[str, dict[str, Any]],
    leaf: dict[str, Any],
    total: int,
    uniform_target: float,
    lookup: dict[tuple[str, int], tuple[float, np.ndarray]],
) -> tuple[tuple[int, ...], dict[str, Any]] | None:
    if leaf.get("state") != "ACTIVE_LEAF":
        return None
    profile, upper, replay = reconstruct_selector_maximum(leaf, lookup, total)
    selector = validate_selector(leaf.get("selector"))
    contribution = math.log2(int(leaf["exact_count"])) + upper
    if not math.isfinite(upper) or not math.isfinite(contribution):
        raise ValueError("active leaf has a non-finite diagnostic")
    selector_digest = coordinate.sha256_bytes(coordinate.canonical_bytes(selector))
    provenance = {
        "checkpoint": {
            "sha256": checkpoint_digest,
            "schema": CHECKPOINT_SCHEMA,
        },
        "node_id": str(leaf["node_id"]),
        "parent_id": leaf.get("parent_id"),
        "branch": leaf.get("branch"),
        "depth": int(leaf.get("depth", 0)),
        "origin": leaf.get("origin"),
        "h2_cell": int(leaf["h2_cell"]),
        "h2_bin_indices": list(leaf["h2_bin_indices"]),
        "path": ancestor_path(nodes, leaf),
        "cell": {
            "lower": list(leaf["lower"]),
            "upper": list(leaf["upper"]),
            "exact_count": str(leaf["exact_count"]),
            "exact_vertex_count": int(leaf["exact_vertex_count"]),
            "exact_vertex_sha256": str(leaf["exact_vertex_sha256"]),
        },
        "current_selector": selector,
        "current_selector_sha256": selector_digest,
        "current_selector_diagnostic": {
            "candidate_status": leaf.get("candidate_status"),
            "candidate_upper_log2": upper,
            "contribution_log2": contribution,
            "positive_residual_bits": max(0.0, upper - uniform_target),
            "component_count": int(leaf.get("component_count", 0)),
            "diagnostic": leaf.get("diagnostic"),
            "selector_replay": replay,
        },
    }
    return profile, provenance


def collect_targets(
    checkpoint: dict[str, Any],
    checkpoint_digest: str,
    total: int,
    uniform_target: float,
    maximum: int,
    lookup: dict[tuple[str, int], tuple[float, np.ndarray]],
    prior_rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    if checkpoint.get("schema") != CHECKPOINT_SCHEMA:
        raise ValueError("unexpected adaptive cumulative checkpoint schema")
    nodes = checkpoint.get("nodes")
    if not isinstance(nodes, dict) or not nodes:
        raise ValueError("adaptive checkpoint node map is empty")
    for identifier, node in nodes.items():
        if not isinstance(node, dict) or node.get("node_id") != identifier:
            raise ValueError("adaptive checkpoint node map is not canonical")

    collected: dict[tuple[int, ...], list[dict[str, Any]]] = defaultdict(list)
    active_count = 0
    missing_count = 0
    for identifier in sorted(nodes):
        node = nodes[identifier]
        if node.get("state") != "ACTIVE_LEAF":
            continue
        active_count += 1
        if "worst_vertex" not in node:
            missing_count += 1
        occurrence = occurrence_from_leaf(
            checkpoint_digest, nodes, node, total, uniform_target, lookup
        )
        if occurrence is None:
            continue
        profile, provenance = occurrence
        collected[profile].append(provenance)

    prepared = []
    for profile, occurrences in collected.items():
        occurrences.sort(
            key=lambda row: (
                -row["current_selector_diagnostic"]["contribution_log2"],
                row["node_id"],
            )
        )
        contributions = [
            row["current_selector_diagnostic"]["contribution_log2"]
            for row in occurrences
        ]
        maximum_contribution = max(contributions)
        acquisition_contribution = logsumexp2(contributions)
        maximum_residual = max(
            row["current_selector_diagnostic"]["positive_residual_bits"]
            for row in occurrences
        )
        profile_list = list(profile)
        profile_digest = coordinate.sha256_bytes(coordinate.canonical_bytes(profile_list))
        prepared.append(
            {
                "profile": profile_list,
                "profile_sha256": profile_digest,
                "support_mask": "0x1ff",
                "physical_weight": str(
                    sum(index * count for index, count in enumerate(profile))
                ),
                "occurrence_count": len(occurrences),
                "maximum_collapsed_candidate_contribution_log2": maximum_contribution,
                "acquisition_contribution_log2": acquisition_contribution,
                "maximum_positive_residual_bits": maximum_residual,
                "provenance": occurrences,
                "original_h2_cells": sorted(set(row["h2_cell"] for row in occurrences)),
                "target_role": "active-worst",
            }
        )
    prior_profiles = [tuple(int(value) for value in row["profile"]) for row in prior_rows]
    prior_set = set(prior_profiles)
    excluded_prior = sum(tuple(row["profile"]) in prior_set for row in prepared)
    prepared = [row for row in prepared if tuple(row["profile"]) not in prior_set]
    for row in prepared:
        profile = tuple(row["profile"])
        distances = [transported_mass(profile, prior) for prior in prior_profiles]
        nearest_distance = min(distances)
        nearest_index = min(
            index for index, distance in enumerate(distances) if distance == nearest_distance
        )
        prior = prior_rows[nearest_index]
        row["nearest_prior_target"] = {
            "supplementary_row": nearest_index,
            "target_id": prior.get("source_target", {}).get("target_id"),
            "profile_sha256": coordinate.sha256_bytes(
                coordinate.canonical_bytes(prior["profile"])
            ),
            "transported_mass_distance": nearest_distance,
        }
    prepared.sort(
        key=lambda row: (
            -row["acquisition_contribution_log2"],
            -row["maximum_positive_residual_bits"],
            row["profile"],
        )
    )
    unique_before_limit = len(prepared)

    # Higher-ranked representatives own transported-mass neighborhoods.  The
    # documented core has no top-16 coverage exception, but record every
    # suppressed neighbor so a future audit can detect a policy change.
    neighborhood_survivors = []
    suppressed_neighbors = []
    for row in prepared:
        profile = tuple(row["profile"])
        near = next(
            (
                kept
                for kept in neighborhood_survivors
                if transported_mass(profile, tuple(kept["profile"])) < NEIGHBORHOOD_RADIUS
            ),
            None,
        )
        if near is None:
            neighborhood_survivors.append(row)
        else:
            suppressed_neighbors.append(
                {
                    "profile_sha256": row["profile_sha256"],
                    "retained_profile_sha256": near["profile_sha256"],
                    "transported_mass_distance": transported_mass(
                        profile, tuple(near["profile"])
                    ),
                }
            )

    # Apply the original-cell cap in E-rank order.  The farthest-first audit
    # records each retained row's separation from the previously retained
    # core, while E-rank remains the primary acquisition order specified by
    # the policy and recommended table.
    core = []
    per_cell: defaultdict[int, int] = defaultdict(int)
    for row in neighborhood_survivors:
        cell = int(row["provenance"][0]["h2_cell"])
        if per_cell[cell] >= 2:
            continue
        prior_core = [tuple(item["profile"]) for item in core]
        row["selection_audit"] = {
            "e_rank_before_diversity": prepared.index(row) + 1,
            "original_h2_cell": cell,
            "cell_ordinal": per_cell[cell] + 1,
            "minimum_transported_distance_to_earlier_core": (
                min(transported_mass(tuple(row["profile"]), profile) for profile in prior_core)
                if prior_core
                else None
            ),
        }
        core.append(row)
        per_cell[cell] += 1
        if len(core) == CORE_TARGET_COUNT:
            break
    farthest_pool = list(core)
    farthest_order = []
    if farthest_pool:
        farthest_order.append(farthest_pool.pop(0))
    while farthest_pool:
        chosen = max(
            farthest_pool,
            key=lambda row: (
                min(
                    transported_mass(tuple(row["profile"]), tuple(prior["profile"]))
                    for prior in farthest_order
                ),
                row["acquisition_contribution_log2"],
                tuple(-value for value in row["profile"]),
            ),
        )
        farthest_pool.remove(chosen)
        farthest_order.append(chosen)
    for diversity_rank, row in enumerate(farthest_order, 1):
        row["selection_audit"]["farthest_first_diversity_rank"] = diversity_rank
    prepared = core
    if maximum > CORE_TARGET_COUNT:
        for transition in TRANSITION_ENDPOINTS[: maximum - CORE_TARGET_COUNT]:
            profile = validate_profile(transition["profile"], total)
            if profile in prior_set or any(tuple(row["profile"]) == profile for row in prepared):
                raise ValueError("transition endpoint is not exact-new")
            parent = nodes.get(transition["parent_id"])
            if parent is None or parent.get("state") != "SPLIT":
                raise ValueError("transition endpoint parent is absent")
            distances = [transported_mass(profile, prior) for prior in prior_profiles]
            nearest_distance = min(distances)
            nearest_index = min(
                index for index, distance in enumerate(distances) if distance == nearest_distance
            )
            profile_list = list(profile)
            profile_digest = coordinate.sha256_bytes(coordinate.canonical_bytes(profile_list))
            prepared.append(
                {
                    "profile": profile_list,
                    "profile_sha256": profile_digest,
                    "support_mask": "0x1ff",
                    "physical_weight": str(
                        sum(index * count for index, count in enumerate(profile))
                    ),
                    "occurrence_count": 1,
                    "maximum_collapsed_candidate_contribution_log2": float(
                        parent["contribution_log2"]
                    ),
                    "acquisition_contribution_log2": float(parent["contribution_log2"]),
                    "maximum_positive_residual_bits": max(
                        0.0, float(parent["candidate_upper_log2"]) - uniform_target
                    ),
                    "original_h2_cells": [int(parent["h2_cell"])],
                    "target_role": "transition-endpoint",
                    "transition": {
                        "transition_id": transition["transition_id"],
                        "parent_id": transition["parent_id"],
                        "coordinate": transition["coordinate"],
                        "regime": transition["regime"],
                        "split": parent["split"],
                        "split_diagnostic": parent.get("split_diagnostic"),
                    },
                    "nearest_prior_target": {
                        "supplementary_row": nearest_index,
                        "target_id": prior_rows[nearest_index]
                        .get("source_target", {})
                        .get("target_id"),
                        "profile_sha256": coordinate.sha256_bytes(
                            coordinate.canonical_bytes(prior_rows[nearest_index]["profile"])
                        ),
                        "transported_mass_distance": nearest_distance,
                    },
                    "provenance": [
                        {
                            "checkpoint": {
                                "sha256": checkpoint_digest,
                                "schema": CHECKPOINT_SCHEMA,
                            },
                            "transition_id": transition["transition_id"],
                            "parent_id": transition["parent_id"],
                            "path": ancestor_path(nodes, parent),
                            "current_selector": parent["selector"],
                            "current_selector_sha256": coordinate.sha256_bytes(
                                coordinate.canonical_bytes(parent["selector"])
                            ),
                        }
                    ],
                }
            )
    prepared = prepared[:maximum]
    for rank, row in enumerate(prepared, 1):
        row["rank"] = rank
        row["target_id"] = f"target-{rank:06d}-{row['profile_sha256'][:12]}"
    return prepared, {
        "active_leaf_count": active_count,
        "active_leaves_without_stored_worst_vertex": missing_count,
        "unique_targets_before_limit": unique_before_limit,
        "excluded_prior_exact_profiles": excluded_prior,
        "suppressed_neighborhood_count": len(suppressed_neighbors),
    }


def validate_bindings(
    checkpoint: dict[str, Any],
    manifest_digest: str,
    base_digest: str,
    supplementary_digest: str,
) -> None:
    bindings = checkpoint.get("source_bindings")
    if not isinstance(bindings, dict):
        raise ValueError("checkpoint source bindings are missing")
    expected = {
        "manifest_sha256": manifest_digest,
        "base_atlas_sha256": base_digest,
        "supplementary_atlas_sha256": supplementary_digest,
    }
    for key, digest in expected.items():
        if bindings.get(key) != digest:
            raise ValueError(f"checkpoint {key} binding changed")


def build_catalogue(
    checkpoint: dict[str, Any],
    checkpoint_digest: str,
    checkpoint_path: Path,
    manifest_digest: str,
    base_digest: str,
    supplementary_digest: str,
    supplementary_path: Path,
    total: int,
    uniform_target: float,
    maximum: int,
    lookup: dict[tuple[str, int], tuple[float, np.ndarray]],
    prior_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    targets, summary = collect_targets(
        checkpoint,
        checkpoint_digest,
        total,
        uniform_target,
        maximum,
        lookup,
        prior_rows,
    )
    return {
        "schema": OUTPUT_SCHEMA,
        "status": "DIAGNOSTIC_TUNING_TARGETS_ONLY",
        "scope": (
            "exact active-leaf worst vertices extracted from one adaptive cumulative "
            "checkpoint; no witness optimization, evaluation, or coverage claim"
        ),
        "manifest_sha256": manifest_digest,
        "atlas_source_id": "atlas-3934dae",
        "atlas_sha256": base_digest,
        "supplementary_atlas_sha256": supplementary_digest,
        "canonical_json": "sorted-compact-json-v1",
        "inputs": [
            {
                "kind": "adaptive-cumulative-checkpoint",
                "path": str(checkpoint_path),
                "schema": CHECKPOINT_SCHEMA,
                "sha256": checkpoint_digest,
            },
            {
                "kind": "supplementary-witness-atlas",
                "path": str(supplementary_path),
                "schema": SUPPLEMENTARY_SCHEMA,
                "sha256": supplementary_digest,
            },
        ],
        "active_leaf_count": summary["active_leaf_count"],
        "active_leaves_without_stored_worst_vertex": summary[
            "active_leaves_without_stored_worst_vertex"
        ],
        "excluded_prior_exact_profiles": summary["excluded_prior_exact_profiles"],
        "suppressed_neighborhood_count": summary["suppressed_neighborhood_count"],
        "unique_targets_before_limit": summary["unique_targets_before_limit"],
        "targets_emitted": len(targets),
        "max_targets": maximum,
        "ranking": "audit-policy-g8-witness-enrichment-targets-from-181-leaves-v1",
        "selection_policy": {
            "acquisition": (
                "E=log2(exact_count)+replayed selector maximum; duplicate profiles "
                "receive log-sum-exp of leaf E values"
            ),
            "maximum_tie": "lexicographically least exact profile",
            "exclude_prior_exact_profiles": True,
            "prior_target_source_sha256": supplementary_digest,
            "transported_mass": "half-l1-v1",
            "neighborhood_radius_exclusive": NEIGHBORHOOD_RADIUS,
            "maximum_core_per_original_h2_cell": 2,
            "core_target_count": CORE_TARGET_COUNT,
            "core_order": "descending acquisition contribution with transported-distance audit",
            "farthest_first": (
                "highest-E seed, then maximum minimum transported distance; recorded as "
                "farthest_first_diversity_rank on the capped recommended core"
            ),
            "variant": "diverse16-with-four-transition-endpoints" if maximum > 12 else "core12",
            "transition_endpoints": [row["transition_id"] for row in TRANSITION_ENDPOINTS]
            if maximum > 12
            else [],
        },
        "targets": targets,
    }


def synthetic_self_test() -> None:
    profile_a = [131069, 1, 80565, 1, 1, 1, 1, 50504, 1]
    profile_b = [196602, 1, 15033, 1, 50503, 1, 1, 1, 1]
    selector = {
        "kind": "mixture",
        "components": [
            {
                "weight": "1/1",
                "witness": {
                    "source_id": "atlas-3934dae",
                    "source_sha256": EXPECTED_BASE_ATLAS_SHA256,
                    "row": 1,
                    "row_sha256": "a" * 64,
                },
            }
        ],
    }

    def leaf(identifier: str, parent: str, profile: list[int], count: int):
        upper = -float(
            coordinate.diagnostic_normalization(
                np.asarray([profile], dtype=np.float64), 262144
            )[0]
        )
        return {
            "node_id": identifier,
            "parent_id": parent,
            "branch": identifier[-1],
            "depth": 1,
            "state": "ACTIVE_LEAF",
            "origin": "synthetic",
            "h2_cell": 1,
            "h2_bin_indices": [0] * 8,
            "lower": [0] * 8,
            "upper": [262135] * 8,
            "exact_count": str(count),
            "exact_vertex_count": 1,
            "exact_vertex_sha256": "b" * 64,
            "exact_vertices": [profile],
            "candidate_status": "FAIL_BINARY64",
            "candidate_upper_log2": upper,
            "contribution_log2": math.log2(count) + upper,
            "component_count": 1,
            "selector": selector,
            "worst_vertex": profile,
            "diagnostic": {"synthetic": True},
        }

    checkpoint = {
        "schema": CHECKPOINT_SCHEMA,
        "nodes": {
            "h2:001": {
                "node_id": "h2:001",
                "parent_id": None,
                "branch": None,
                "depth": 0,
                "state": "SPLIT",
                "split": {"k": 0, "threshold": "5"},
            },
            "h2:001/L": leaf("h2:001/L", "h2:001", profile_a, 5),
            "h2:001/R": leaf("h2:001/R", "h2:001", profile_b, 1000),
            "h2:001/D": leaf("h2:001/D", "h2:001", profile_b, 10),
        },
    }
    lookup = {("atlas-3934dae", 1): (0.0, np.zeros(CLASSES, dtype=np.float64))}
    prior_rows = [
        {
            "profile": [262136] + [1] * 8,
            "source_target": {"target_id": "synthetic-prior"},
        }
    ]
    targets, summary = collect_targets(
        checkpoint,
        "c" * 64,
        262144,
        -168.0,
        maximum=2,
        lookup=lookup,
        prior_rows=prior_rows,
    )
    if (
        summary["active_leaf_count"] != 3
        or summary["unique_targets_before_limit"] != 2
        or len(targets) != 2
        or targets[0]["profile"] != profile_b
        or targets[0]["occurrence_count"] != 2
        or targets[0]["provenance"][0]["node_id"] != "h2:001/R"
        or [step["node_id"] for step in targets[0]["provenance"][0]["path"]]
        != ["h2:001", "h2:001/R"]
    ):
        raise AssertionError("synthetic adaptive target extraction failed")
    # Exercise the tuner's public target validation without invoking tuning.
    from tune_packet_group_g8_supplementary_atlas import select_targets

    if select_targets({"targets": targets}, 2) != targets:
        raise AssertionError("tuner rejected synthetic adaptive targets")
    print("synthetic_unique_targets=2")
    print("synthetic_occurrences=3")
    print("status=SYNTHETIC_ADAPTIVE_TARGET_EXTRACTION_PASS")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--manifest", type=Path, default=coordinate.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=coordinate.DEFAULT_ATLAS)
    parser.add_argument("--supplementary-atlas", type=Path, default=DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-targets", type=int, default=16)
    parser.add_argument("--synthetic-self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if not 1 <= args.max_targets <= 16:
        raise ValueError("max-targets must lie in [1,16]")

    manifest_digest = coordinate.sha256_path(args.manifest)
    base_digest = coordinate.sha256_path(args.atlas)
    supplementary_digest = coordinate.sha256_path(args.supplementary_atlas)
    checkpoint_digest = coordinate.sha256_path(args.checkpoint)
    if manifest_digest != EXPECTED_MANIFEST_SHA256:
        raise ValueError("manifest digest changed")
    if base_digest != EXPECTED_BASE_ATLAS_SHA256:
        raise ValueError("base atlas digest changed")
    if supplementary_digest != EXPECTED_SUPPLEMENTARY_SHA256:
        raise ValueError("supplementary atlas digest changed")
    manifest = json.loads(args.manifest.read_bytes())
    atlas = json.loads(args.atlas.read_bytes())
    supplementary = json.loads(args.supplementary_atlas.read_bytes())
    checkpoint = json.loads(args.checkpoint.read_bytes())
    if supplementary.get("schema") != SUPPLEMENTARY_SCHEMA:
        raise ValueError("supplementary atlas schema changed")
    if supplementary.get("manifest_sha256") != manifest_digest:
        raise ValueError("supplementary atlas binds another manifest")
    validate_bindings(checkpoint, manifest_digest, base_digest, supplementary_digest)
    manifest_source = manifest["witness_sources"][0]
    base_bindings = manifest_source["rows"]
    base_bank = coordinate.build_atlas_bank(atlas, manifest_source, base_bindings)
    supplementary_bank = dominance.load_supplementary_atlas(
        args.supplementary_atlas, manifest_digest, base_digest
    )
    bank = dominance.combined_bank(base_bank, [supplementary_bank])
    lookup = witness_lookup(bank)
    prior_rows = supplementary["rows"]

    total = int(manifest["parameters"]["M"])
    probability_target = float(
        manifest["parameters"]["probability_target_log2"].split("/")[0]
    )
    uniform_target = probability_target - math.log2(
        int(manifest["census"]["feasible_profile_count"])
    )
    catalogue = build_catalogue(
        checkpoint,
        checkpoint_digest,
        args.checkpoint,
        manifest_digest,
        base_digest,
        supplementary_digest,
        args.supplementary_atlas,
        total,
        uniform_target,
        args.max_targets,
        lookup,
        prior_rows,
    )
    digest = coordinate.write_canonical(args.output, catalogue)
    print(f"active_leaf_count={catalogue['active_leaf_count']}")
    print(f"unique_targets={catalogue['unique_targets_before_limit']}")
    print(f"targets_emitted={catalogue['targets_emitted']}")
    print(f"output={args.output}")
    print(f"output_sha256={digest}")
    print("status=DIAGNOSTIC_TUNING_TARGETS_ONLY")


if __name__ == "__main__":
    main()
