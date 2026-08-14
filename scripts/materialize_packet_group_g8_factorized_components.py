#!/usr/bin/env python3
"""Freeze the positive-weight g=8 factorized components and their selectors.

This is a deterministic materialization step, not an outward proof replay.  It
reconstructs each selected inner Collatz vector from its frozen binary64
parameters, separates inner and outer rows, and writes content-addressed
selector bindings plus a factorized v2 manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

from packet_group_drive_stratified import block_histograms, point_caps
from packet_group_profile_bound import split_cap_table
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import best_witness
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "out" / "g8_factorized_geometry_persistence_wave.json"
GEOMETRY_SHA256 = "636103bd034315c9c4381f3319d3d908f98be93e607a8efe3d1fd4fcff6bfa4e"
BASE_MANIFEST = ROOT / "G8_SUPPORT_MANIFEST.json"
BASE_MANIFEST_SHA256 = "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
VERIFIER = ROOT / "scripts" / "certify_packet_group_g8_factorized_selector.py"

SOURCE_PATHS = {
    "c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e":
        ROOT / "out" / "g8_support_seed_atlas.json",
    "d50af6e55d5bc1d0aa0a9de65025684fdfe9989cce0014c808a059c2ca356d6c":
        ROOT / "out" / "g8_highdim_supplementary_final32.json",
    "530f042f0f4be3fd64ad166c68380302dbb5abb745ffb3745a4fd4ce88b1236c":
        ROOT / "out" / "g8_adaptive_cumulative_supplementary_diverse16.json",
    "edc9c8b1714024d99af4c8d97352a8d20ec5f880202962667030d1ea2d9ae6f3":
        ROOT / "out" / "g8_cell_aware_factorized_columns" / "candidate_catalogue.json",
    "b6d324477840ea5708776c777b45ca0f3272ea3c3e9a217729a1984301d37250":
        ROOT / "out" / "g8_iterative_factorized_columns" / "round_001" / "candidate_catalogue.json",
}

INNER_SOURCE_ID = "g8-factorized-inner-components-v1"
OUTER_SOURCE_ID = "g8-factorized-outer-components-v1"
INNER_SCHEMA = "permute-conv.packet-group-g8-factorized-inner-components.v1"
OUTER_SCHEMA = "permute-conv.packet-group-g8-factorized-outer-components.v1"
SELECTOR_SCHEMA = "permute-conv.packet-group-g8-factorized-selector-bindings.v1"
MANIFEST_SCHEMA = "packet-group-g8-factorized-manifest-v2"
SELECTOR_KIND = "independent-sum-mixture-v1"
NORMALIZATION_RULE = "packet-profile-orbit-v1"
INNER_INTERFACE = "uniform-conditional-on-complete-outer-word-v1"
COMBINATION_RULE = "sum-marginals-subtract-normalization-once-v1"


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest_path(path: Path) -> str:
    return digest_bytes(path.read_bytes())


def digest_row(row: Any) -> str:
    return digest_bytes(canonical_bytes(row))


def write_canonical(path: Path, value: Any) -> str:
    payload = canonical_bytes(value) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(payload)
    os.replace(temporary, path)
    return digest_bytes(payload)


def reference_key(reference: dict[str, Any]) -> tuple[str, str, int, str]:
    required = ("source_id", "source_sha256", "row", "row_sha256")
    if set(reference) != set(required):
        raise ValueError("component reference fields are not canonical")
    source_id = str(reference["source_id"])
    source_digest = str(reference["source_sha256"])
    row = int(reference["row"])
    row_digest = str(reference["row_sha256"])
    if not source_id or source_digest not in SOURCE_PATHS or row < 0:
        raise ValueError("component reference is outside the frozen source closure")
    return source_id, source_digest, row, row_digest


def parse_weight(text: Any) -> Fraction:
    if not isinstance(text, str) or "/" not in text:
        raise ValueError("selector weight must be an explicit rational")
    numerator, denominator = text.split("/", 1)
    value = Fraction(int(numerator), int(denominator))
    if value <= 0 or f"{value.numerator}/{value.denominator}" != text:
        raise ValueError("selector weight must be positive and canonical")
    return value


def selected_references(geometry: dict[str, Any]) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, Any]]]:
    selected: dict[str, dict[tuple[str, str, int, str], dict[str, Any]]] = {
        "inner": {}, "outer": {},
    }
    leaves = []
    for node in sorted(geometry["nodes"].values(), key=lambda row: row["node_id"]):
        if node["state"] != "ACTIVE_LEAF":
            continue
        selector = node.get("selector", {})
        if selector.get("kind") != "independent-components":
            raise ValueError("active leaf lacks the frozen factorized selector")
        leaf = {
            "node_id": node["node_id"],
            "h2_cell": int(node["h2_cell"]),
            "exact_count": str(node["exact_count"]),
            "exact_vertex_sha256": node["exact_vertex_sha256"],
            "marginals": {},
        }
        for role in ("inner", "outer"):
            components = selector.get(role, {}).get("components")
            if not isinstance(components, list) or not components:
                raise ValueError(f"{role} marginal is empty")
            total = Fraction(0)
            parsed = []
            seen = set()
            for component in components:
                weight = parse_weight(component.get("weight"))
                reference = component.get("component")
                if not isinstance(reference, dict):
                    raise ValueError("selector component reference is malformed")
                key = reference_key(reference)
                if key in seen:
                    raise ValueError("selector marginal repeats one component")
                seen.add(key)
                total += weight
                selected[role][key] = reference
                parsed.append((key, weight))
            if total != 1:
                raise ValueError(f"{role} marginal does not sum exactly to one")
            leaf["marginals"][role] = parsed
        leaves.append(leaf)
    if len(leaves) != 197:
        raise ValueError("frozen factorized geometry no longer has 197 active leaves")
    return {
        role: [selected[role][key] for key in sorted(selected[role])]
        for role in ("inner", "outer")
    }, leaves


def load_sources() -> dict[str, dict[str, Any]]:
    result = {}
    for expected, path in SOURCE_PATHS.items():
        actual = digest_path(path)
        if actual != expected:
            raise ValueError(f"frozen component source digest changed: {path}")
        result[expected] = json.loads(path.read_bytes())
    return result


def resolve_row(reference: dict[str, Any], sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    source_id, source_digest, index, expected_row_digest = reference_key(reference)
    artifact = sources[source_digest]
    rows = artifact.get("rows")
    if not isinstance(rows, list) or not 0 <= index < len(rows):
        raise ValueError("component row is outside its source")
    row = rows[index]
    if source_digest in tuple(SOURCE_PATHS)[:3]:
        actual = digest_row(row)
    else:
        actual = str(row.get("parameter_digest", ""))
        if actual != digest_row(row.get("parameters")):
            raise ValueError("candidate parameter digest is inconsistent")
        lane = source_id.rsplit("-", 1)[-1]
        if row.get("lane") != lane or row.get("parameters", {}).get("lane") != lane:
            raise ValueError("candidate role does not match its source identifier")
    if actual != expected_row_digest:
        raise ValueError("component row digest does not match its reference")
    return row


def role_parameters(row: dict[str, Any], role: str) -> dict[str, Any]:
    if "parameters" in row:
        parameters = row["parameters"]
    else:
        parameters = row[role]
    if not isinstance(parameters, dict):
        raise ValueError("component parameters are malformed")
    return parameters


def freeze_inner(reference: dict[str, Any], row: dict[str, Any], caps: Any) -> dict[str, Any]:
    parameters = role_parameters(row, "inner")
    fugacities = np.asarray(parameters.get("fugacities"), dtype=np.float64)
    pole = float(parameters.get("pole"))
    iterations = int(parameters.get("final_iterations", parameters.get("collatz_best_iteration", -1)))
    if (
        fugacities.shape != (9,) or not np.all(np.isfinite(fugacities))
        or np.any(fugacities <= 0.0) or not 0.0 < pole < 1.0
        or not 0 <= iterations <= 64
    ):
        raise ValueError("selected inner parameters are outside the replay contract")
    kernel = SharedDriveStratifiedKernel(
        block_histograms(8, fugacities), point_caps(8, fugacities), caps, pole
    )
    _eigenvalue, _domination, values, _worst, selected, _score = best_witness(
        kernel, iterations, INNER_BLOCKS
    )
    recorded = parameters.get("collatz_best_iteration")
    if recorded is not None and selected != int(recorded):
        raise ValueError("deterministic Collatz replay selected another iteration")
    vector = [float(value) for value in values]
    if len(vector) != 65 or any(not math.isfinite(value) or value <= 0.0 for value in vector):
        raise ValueError("deterministic Collatz replay produced an invalid vector")
    return {
        "component_role": "inner",
        "origin_reference": reference,
        "parameters": {
            "pole": pole,
            "fugacities": [float(value) for value in fugacities],
            "collatz_vector": vector,
        },
        "materialization": {
            "method": "binary64-shared-drive-best-power-trajectory-v1",
            "maximum_iteration": iterations,
            "selected_iteration": selected,
            "inner_blocks": int(INNER_BLOCKS),
        },
    }


def freeze_outer(reference: dict[str, Any], row: dict[str, Any]) -> dict[str, Any]:
    parameters = role_parameters(row, "outer")
    logs = [float(value) for value in parameters.get("log_variables", [])]
    coefficients = [float(value) for value in parameters.get("band_coefficients", [])]
    theta = float(parameters.get("pair_cauchy_theta"))
    if (
        len(logs) != 9 or logs[0] != 0.0 or len(coefficients) != 3
        or not all(math.isfinite(value) for value in (*logs, *coefficients, theta))
    ):
        raise ValueError("selected outer parameters are outside the replay contract")
    return {
        "component_role": "outer",
        "origin_reference": reference,
        "parameters": {
            "log_variables": logs,
            "band_coefficients": coefficients,
            "pair_cauchy_theta": theta,
        },
    }


def source_artifact(role: str, rows: list[dict[str, Any]], geometry_digest: str) -> dict[str, Any]:
    return {
        "schema": INNER_SCHEMA if role == "inner" else OUTER_SCHEMA,
        "status": "FROZEN_PARAMETERS_PENDING_INDEPENDENT_OUTWARD_REPLAY",
        "component_role": role,
        "source_geometry_sha256": geometry_digest,
        "row_order": "lexicographic-origin-reference-v1",
        "rows": rows,
    }


def declarations(artifact: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"row": index, "row_sha256": digest_row(row)}
        for index, row in enumerate(artifact["rows"])
    ]


def materialized_reference(source_id: str, source_digest: str, index: int, row: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_id": source_id,
        "source_sha256": source_digest,
        "row": index,
        "row_sha256": digest_row(row),
    }


def build_selector_artifact(
    leaves: list[dict[str, Any]], maps: dict[str, dict[tuple[str, str, int, str], dict[str, Any]]],
    geometry_digest: str,
) -> dict[str, Any]:
    rows = []
    for leaf in leaves:
        selector = {
            "kind": SELECTOR_KIND,
            "normalization": NORMALIZATION_RULE,
        }
        for role in ("inner", "outer"):
            selector[f"{role}_marginal"] = [
                {"component": maps[role][key], "weight": f"{weight.numerator}/{weight.denominator}"}
                for key, weight in leaf["marginals"][role]
            ]
        rows.append({
            "node_id": leaf["node_id"], "h2_cell": leaf["h2_cell"],
            "exact_count": leaf["exact_count"],
            "exact_vertex_sha256": leaf["exact_vertex_sha256"],
            "selector": selector,
        })
    return {
        "schema": SELECTOR_SCHEMA,
        "status": "FROZEN_RATIONAL_SELECTORS_PENDING_OUTWARD_REPLAY",
        "source_geometry_sha256": geometry_digest,
        "active_leaf_count": len(rows),
        "rows": rows,
    }


def source_binding(source_id: str, path: Path, digest: str, artifact: dict[str, Any], role: str) -> dict[str, Any]:
    return {
        "source_id": source_id,
        "path": path.relative_to(ROOT).as_posix(),
        "sha256": digest,
        "schema": artifact["schema"],
        "component_role": role,
        "rows": declarations(artifact),
    }


def validate_emitted(manifest_path: Path) -> dict[str, int]:
    manifest = json.loads(manifest_path.read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("emitted manifest schema changed")
    if manifest.get("recombination_contract") != {
        "selector_kind": SELECTOR_KIND,
        "inner_interface": INNER_INTERFACE,
        "combination_rule": COMBINATION_RULE,
        "normalization": NORMALIZATION_RULE,
    }:
        raise ValueError("emitted recombination contract changed")
    for binding_name in ("supersedes", "geometry_source"):
        binding = manifest[binding_name]
        if digest_path(ROOT / binding["path"]) != binding["sha256"]:
            raise ValueError(f"emitted {binding_name} digest mismatch")
    if manifest["geometry_source"]["active_leaf_count"] != 197:
        raise ValueError("emitted geometry leaf count changed")
    for source in manifest["fixed_sources"]:
        if digest_path(ROOT / source["path"]) != source["sha256"]:
            raise ValueError(f"emitted fixed-source digest mismatch: {source['source_id']}")
    catalog: dict[tuple[str, int], tuple[str, str]] = {}
    counts = {}
    for source in manifest["component_sources"]:
        path = ROOT / source["path"]
        if digest_path(path) != source["sha256"]:
            raise ValueError("emitted component source digest mismatch")
        artifact = json.loads(path.read_bytes())
        if artifact["schema"] != source["schema"] or artifact["component_role"] != source["component_role"]:
            raise ValueError("emitted component source metadata mismatch")
        if declarations(artifact) != source["rows"]:
            raise ValueError("emitted component row catalogue mismatch")
        counts[source["component_role"]] = len(artifact["rows"])
        for index, row in enumerate(artifact["rows"]):
            catalog[(source["source_id"], index)] = (digest_row(row), source["component_role"])
            if source["component_role"] == "inner":
                values = row["parameters"]["collatz_vector"]
                if len(values) != 65 or any(not math.isfinite(value) or value <= 0.0 for value in values):
                    raise ValueError("emitted Collatz vector is invalid")
    selector_binding = manifest["selector_source"]
    selector_path = ROOT / selector_binding["path"]
    if digest_path(selector_path) != selector_binding["sha256"]:
        raise ValueError("emitted selector source digest mismatch")
    selectors = json.loads(selector_path.read_bytes())
    if selectors.get("schema") != SELECTOR_SCHEMA:
        raise ValueError("emitted selector source schema changed")
    if selectors["active_leaf_count"] != len(selectors["rows"]) or len(selectors["rows"]) != 197:
        raise ValueError("emitted selector leaf catalogue is incomplete")
    used = {"inner": set(), "outer": set()}
    for leaf in selectors["rows"]:
        for role in ("inner", "outer"):
            marginal = leaf["selector"][f"{role}_marginal"]
            if sum((parse_weight(row["weight"]) for row in marginal), Fraction(0)) != 1:
                raise ValueError("emitted marginal mass changed")
            for component in marginal:
                reference = component["component"]
                key = (reference["source_id"], int(reference["row"]))
                if key not in catalog or catalog[key] != (reference["row_sha256"], role):
                    raise ValueError("emitted selector reference is missing or cross-role")
                used[role].add(key)
    if any(len(used[role]) != counts[role] for role in used):
        raise ValueError("emitted role source contains an unused row")
    return {"leaves": 197, **counts}


def run(args: argparse.Namespace) -> dict[str, Any]:
    if digest_path(args.geometry) != GEOMETRY_SHA256:
        raise ValueError("frozen factorized geometry digest changed")
    if digest_path(args.base_manifest) != BASE_MANIFEST_SHA256:
        raise ValueError("base support manifest digest changed")
    geometry = json.loads(args.geometry.read_bytes())
    references, leaves = selected_references(geometry)
    sources = load_sources()
    caps = split_cap_table()
    inner_rows = [freeze_inner(reference, resolve_row(reference, sources), caps) for reference in references["inner"]]
    outer_rows = [freeze_outer(reference, resolve_row(reference, sources)) for reference in references["outer"]]
    inner_artifact = source_artifact("inner", inner_rows, GEOMETRY_SHA256)
    outer_artifact = source_artifact("outer", outer_rows, GEOMETRY_SHA256)
    inner_digest = write_canonical(args.inner_output, inner_artifact)
    outer_digest = write_canonical(args.outer_output, outer_artifact)
    maps = {"inner": {}, "outer": {}}
    for role, rows, source_id, source_digest in (
        ("inner", inner_rows, INNER_SOURCE_ID, inner_digest),
        ("outer", outer_rows, OUTER_SOURCE_ID, outer_digest),
    ):
        for index, row in enumerate(rows):
            maps[role][reference_key(row["origin_reference"])] = materialized_reference(
                source_id, source_digest, index, row
            )
    selector_artifact = build_selector_artifact(leaves, maps, GEOMETRY_SHA256)
    selector_digest = write_canonical(args.selector_output, selector_artifact)
    base = json.loads(args.base_manifest.read_bytes())
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "run_id": "g8-factorized-197-leaf-materialization-v2",
        "status": "MATERIALIZED_PARAMETERS_PENDING_INDEPENDENT_OUTWARD_REPLAY",
        "supersedes": {
            "path": args.base_manifest.relative_to(ROOT).as_posix(),
            "sha256": BASE_MANIFEST_SHA256,
            "schema": base["schema"],
        },
        "parameters": base["parameters"],
        "arithmetic": {**base["arithmetic"], "canonical_json": "sorted-compact-json-v1"},
        "recombination_contract": {
            "selector_kind": SELECTOR_KIND,
            "inner_interface": INNER_INTERFACE,
            "combination_rule": COMBINATION_RULE,
            "normalization": NORMALIZATION_RULE,
        },
        "geometry_source": {
            "path": args.geometry.relative_to(ROOT).as_posix(),
            "sha256": GEOMETRY_SHA256,
            "schema": geometry["schema"],
            "active_leaf_count": 197,
        },
        "selector_source": {
            "path": args.selector_output.relative_to(ROOT).as_posix(),
            "sha256": selector_digest,
            "schema": SELECTOR_SCHEMA,
        },
        "component_sources": [
            source_binding(INNER_SOURCE_ID, args.inner_output, inner_digest, inner_artifact, "inner"),
            source_binding(OUTER_SOURCE_ID, args.outer_output, outer_digest, outer_artifact, "outer"),
        ],
        "fixed_sources": [
            {"source_id": "materializer", "path": Path(__file__).resolve().relative_to(ROOT).as_posix(), "sha256": digest_path(Path(__file__).resolve())},
            {"source_id": "factorized-verifier", "path": VERIFIER.relative_to(ROOT).as_posix(), "sha256": digest_path(VERIFIER)},
            *[
                {"source_id": f"origin-{digest[:16]}", "path": path.relative_to(ROOT).as_posix(), "sha256": digest}
                for digest, path in SOURCE_PATHS.items()
            ],
        ],
        "scope_limit": "No component was outward-evaluated. No vertex or union endpoint is certified by this materialization.",
    }
    manifest_digest = write_canonical(args.manifest_output, manifest)
    counts = validate_emitted(args.manifest_output)
    return {"manifest_sha256": manifest_digest, **counts}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=GEOMETRY)
    parser.add_argument("--base-manifest", type=Path, default=BASE_MANIFEST)
    parser.add_argument("--inner-output", type=Path, default=ROOT / "out" / "g8_factorized_inner_components_v1.json")
    parser.add_argument("--outer-output", type=Path, default=ROOT / "out" / "g8_factorized_outer_components_v1.json")
    parser.add_argument("--selector-output", type=Path, default=ROOT / "out" / "g8_factorized_selector_bindings_v1.json")
    parser.add_argument("--manifest-output", type=Path, default=ROOT / "G8_FACTORIZED_MANIFEST_V2.json")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    result = validate_emitted(args.manifest_output) if args.validate_only else run(args)
    for key, value in result.items():
        print(f"{key}={value}")
    print("status=PASS_G8_FACTORIZED_COMPONENT_MATERIALIZATION")


if __name__ == "__main__":
    main()
