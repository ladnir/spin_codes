#!/usr/bin/env python3
"""Outward-replay the frozen 197-leaf g=8 factorized selector ledger."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any

import certify_packet_group_g8_factorized_selector as factorized
import certify_packet_group_triangle_ledger as outward
from outward_log2 import Interval, log2_int


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "G8_FACTORIZED_MANIFEST_V2.json"
SELECTORS = ROOT / "out" / "g8_factorized_selector_bindings_v1.json"
GEOMETRY = ROOT / "out" / "g8_factorized_geometry_persistence_wave.json"
DEFAULT_CERTIFICATE = ROOT / "out" / "g8_factorized_full_support_certificate_v1.json"
DEFAULT_RECEIPT = ROOT / "out" / "g8_factorized_full_support_outward_receipt.json"
RECEIPT_SCHEMA = "permute-conv.packet-group-g8-factorized-outward-receipt.v1"

_WORKER_CONTEXT: dict[str, Any] | None = None


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_canonical(path: Path, value: Any) -> str:
    payload = canonical_bytes(value) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(payload)
    os.replace(temporary, path)
    return hashlib.sha256(payload).hexdigest()


def _initialize_worker(manifest_text: str) -> None:
    global _WORKER_CONTEXT
    path = Path(manifest_text)
    manifest = json.loads(path.read_bytes())
    _WORKER_CONTEXT = factorized.validate_manifest(manifest, path)


def _interval_wire(value: Interval) -> list[str]:
    return [str(value.lo), str(value.hi)]


def _harden_worker(reference: tuple[str, int]) -> dict[str, Any]:
    if _WORKER_CONTEXT is None:
        raise RuntimeError("factorized hardening worker was not initialized")
    hardened = factorized.harden_reference(reference, _WORKER_CONTEXT)
    return {
        "source_id": reference[0],
        "row": reference[1],
        "role": hardened["role"],
        "constant_log2_interval": _interval_wire(hardened["constant"]),
        "charge_log2_intervals": [
            None if value is None else _interval_wire(value)
            for value in hardened["charges"]
        ],
    }


def _interval_from_wire(value: list[str]) -> Interval:
    return Interval(Decimal(value[0]), Decimal(value[1]))


def _restore_hardened(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": f"{row['source_id']}:{row['row']}",
        "role": row["role"],
        "constant": _interval_from_wire(row["constant_log2_interval"]),
        "charges": tuple(
            None if value is None else _interval_from_wire(value)
            for value in row["charge_log2_intervals"]
        ),
    }


def build_certificate(geometry: dict[str, Any], selectors: dict[str, Any]) -> dict[str, Any]:
    by_leaf = {row["node_id"]: row for row in selectors["rows"]}
    if len(by_leaf) != 197 or set(by_leaf) != {
        row["node_id"] for row in geometry["nodes"].values()
        if row["state"] == "ACTIVE_LEAF"
    }:
        raise ValueError("selector binding and active geometry leaves disagree")
    nodes = {}
    for identifier, source in geometry["nodes"].items():
        node = {
            key: source[key]
            for key in (
                "node_id", "parent_id", "h2_cell", "h2_bin_indices", "lower",
                "upper", "exact_count",
            )
        }
        if source["state"] == "ACTIVE_LEAF":
            binding = by_leaf[identifier]
            if (
                binding["exact_count"] != source["exact_count"]
                or binding["exact_vertex_sha256"] != source["exact_vertex_sha256"]
            ):
                raise ValueError("selector binding changed its exact leaf geometry")
            node.update(state="CERTIFIED_LEAF", selector=binding["selector"])
        elif source["state"] == "SPLIT":
            split = source["split"]
            node.update(
                state="SPLIT", left=source["left_id"], right=source["right_id"],
                split={
                    "coordinate": int(split["k"]),
                    "threshold": int(split["t"]),
                    "profile_coefficients": split["coefficients"],
                    "profile_threshold": int(split["threshold"]),
                },
            )
        else:
            raise ValueError("factorized geometry contains an unresolved node")
        nodes[identifier] = node
    return {
        "schema": factorized.CERTIFICATE_SCHEMA,
        "status": "FROZEN_PARAMETERS_PENDING_OUTWARD_REPLAY",
        "source_geometry_sha256": sha256_path(GEOMETRY),
        "source_selector_sha256": sha256_path(SELECTORS),
        "nodes": nodes,
    }


def maximum_interval(values: list[Interval]) -> Interval:
    if not values:
        raise ValueError("cannot maximize an empty interval list")
    return Interval(max(value.lo for value in values), max(value.hi for value in values))


def run(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    manifest = json.loads(args.manifest.read_bytes())
    context = factorized.validate_manifest(manifest, args.manifest)
    selectors = json.loads(args.selectors.read_bytes())
    geometry = json.loads(args.geometry.read_bytes())
    if manifest["selector_source"]["sha256"] != sha256_path(args.selectors):
        raise ValueError("selector artifact does not match the factorized manifest")
    if manifest["geometry_source"]["sha256"] != sha256_path(args.geometry):
        raise ValueError("geometry artifact does not match the factorized manifest")
    certificate = build_certificate(geometry, selectors)
    certificate_digest = write_canonical(args.certificate_output, certificate)
    replay = factorized.replay_cumulative_geometry(
        certificate, expected_leaves=197
    )
    references = sorted(context["component_rows"])
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=_initialize_worker,
        initargs=(str(args.manifest.resolve()),),
    ) as executor:
        component_rows = list(executor.map(_harden_worker, references, chunksize=1))
    context["hardened"] = {
        (row["source_id"], int(row["row"])): _restore_hardened(row)
        for row in component_rows
    }

    leaf_results = []
    inequality_digest = hashlib.sha256()
    for leaf in sorted(replay["leaves"], key=lambda row: row["node"]["node_id"]):
        node = leaf["node"]
        selector = factorized.parse_selector(node["selector"], context)
        values = []
        for profile in leaf["vertices"]:
            value = factorized.evaluate_selector(profile, selector, context)
            values.append(value)
            inequality_digest.update(canonical_bytes({
                "node_id": node["node_id"],
                "profile": [str(value) for value in profile],
                "selector": node["selector"],
                "value_log2_interval": _interval_wire(value),
            }))
        maximum = maximum_interval(values)
        leaf_results.append({
            "node_id": node["node_id"],
            "exact_count": node["exact_count"],
            "vertex_count": len(values),
            "maximum_log2_interval": _interval_wire(maximum),
        })

    maxima = [_interval_from_wire(row["maximum_log2_interval"]) for row in leaf_results]
    global_maximum = maximum_interval(maxima)
    collapsed = global_maximum + log2_int(replay["root_count"])
    expanded_terms = [
        value + log2_int(int(row["exact_count"]))
        for value, row in zip(maxima, leaf_results)
    ]
    expanded_upper = outward._log2_sum_exp(
        Interval.exact(value.hi) for value in expanded_terms
    )
    expanded = Interval(max(value.lo for value in expanded_terms), expanded_upper.hi)
    aggregation = expanded if expanded.hi <= collapsed.hi else collapsed
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "status": "VERIFIED_COMPLETE" if aggregation.hi <= -40 else "VERIFIED_COMPLETE_DID_NOT_CLOSE",
        "passed": aggregation.hi <= -40,
        "manifest_sha256": sha256_path(args.manifest),
        "verifier_source_sha256": sha256_path(Path(__file__).resolve()),
        "factorized_core_source_sha256": sha256_path(Path(factorized.__file__).resolve()),
        "certificate_sha256": certificate_digest,
        "geometry_sha256": sha256_path(args.geometry),
        "selector_sha256": sha256_path(args.selectors),
        "root_count": str(replay["root_count"]),
        "root_count_count": replay["root_count_count"],
        "leaf_count": replay["leaf_count"],
        "component_count": len(component_rows),
        "inner_component_count": sum(row["role"] == "inner" for row in component_rows),
        "outer_component_count": sum(row["role"] == "outer" for row in component_rows),
        "vertex_inequality_count": sum(row["vertex_count"] for row in leaf_results),
        "inequality_sha256": inequality_digest.hexdigest(),
        "collapsed_log2_interval": _interval_wire(collapsed),
        "expanded_log2_interval": _interval_wire(expanded),
        "aggregation_used": "expanded" if aggregation is expanded else "collapsed",
        "selected_log2_interval": _interval_wire(aggregation),
        "pure_python_outward_point_caps": True,
        "components": component_rows,
        "leaves": leaf_results,
        "runtime_seconds": time.perf_counter() - started,
    }
    receipt_digest = write_canonical(args.receipt_output, receipt)
    return {"receipt_sha256": receipt_digest, **receipt}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--selectors", type=Path, default=SELECTORS)
    parser.add_argument("--geometry", type=Path, default=GEOMETRY)
    parser.add_argument("--certificate-output", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--receipt-output", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        raise SystemExit("workers must lie in [1,4]")
    result = run(args)
    for key in (
        "receipt_sha256", "status", "passed", "leaf_count", "component_count",
        "vertex_inequality_count", "aggregation_used", "selected_log2_interval",
        "runtime_seconds",
    ):
        print(f"{key}={result[key]}")


if __name__ == "__main__":
    main()
