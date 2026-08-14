#!/usr/bin/env python3
"""Freeze corrected BL2 outer columns and replay the frozen 197-leaf geometry."""

from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import replay_packet_group_g8_independent_component_mixtures as mixtures
import run_packet_group_g8_adaptive_cumulative_engine as engine
import run_packet_group_g8_factorized_geometry_persistence_wave as geometry_wave
import probe_packet_group_g8_outer_family_prices as probe
import probe_packet8_three_band_linear_bl2 as bl2


ROOT = Path(__file__).resolve().parents[1]
PROBE = ROOT / "out" / "g8_outer_family_price_probe_top8.json"
PROBE_SHA256 = "797e097f39fdfea4b30bd62a7c8c7efb11450beb0446776b5d3ccc07da85c6e5"
GEOMETRY = probe.GEOMETRY
V2_INNER = probe.V2_INNER
V2_OUTER = probe.V2_OUTER
DEFAULT_SOURCE = ROOT / "out" / "g8_corrected_bl2_outer_components.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_corrected_bl2_factorized_replay_197.json"
SOURCE_SCHEMA = "permute-conv.packet-group-g8-corrected-bl2-outer-components.v1"
SCHEMA = "permute-conv.packet-group-g8-corrected-bl2-factorized-replay.v1"
DENOMINATOR = 1 << 30


def key(reference: dict[str, Any]):
    return str(reference["source_sha256"]), int(reference["row"]), str(reference["row_sha256"])


def select_v2(full: geometry_wave.FactorizedBank, artifact: dict[str, Any], role: str):
    references = full.inner_references if role == "inner" else full.outer_references
    constants = full.inner_constants if role == "inner" else full.outer_constants
    charges = full.inner_charges if role == "inner" else full.outer_charges
    lookup = {key(reference): index for index, reference in enumerate(references)}
    indices = [lookup[key(row["origin_reference"])] for row in artifact["rows"]]
    if len(indices) != len(set(indices)):
        raise ValueError(f"v2 {role} source repeats a component")
    take = np.asarray(indices, dtype=np.int64)
    return constants[take], charges[take], tuple(references[index] for index in indices)


def graph_correction(q: np.ndarray) -> float:
    # Same orientation as probe_packet_group_fixed_atlas.py: new minus old.
    return 128.0 * max(
        float(q[new] - q[old])
        for old in range(9) for new in range(9) if abs(new - old) <= 1
    )


def reconstruct_row(row: dict[str, Any]) -> dict[str, Any]:
    profile = np.asarray([float(Fraction(value)) for value in row["barycenter_exact"]])
    parameters = row["three_band_linear_bl"]["parameters"]
    log_t = np.asarray(parameters["log_variables"], dtype=np.float64)
    q = log_t / math.log(2.0)
    coefficients = np.asarray(parameters["band_coefficients"], dtype=np.float64)
    group_logs = bl2.group_log_moments(log_t)
    p0, p1 = float(coefficients[0]), float(coefficients[1])
    norm0 = p0 * bl2.logsumexp(bl2.BINOMIAL_LOG_PROBABILITIES + group_logs / p0)
    norm1 = p1 * bl2.logsumexp(bl2.BINOMIAL_LOG_PROBABILITIES + group_logs / p1)
    unpunctured_constant = (
        bl2.K * math.log(2.0) + 256.0 * (42.0 * norm0 + 86.0 * norm1)
    ) / math.log(2.0)
    correction = graph_correction(q)
    constant = unpunctured_constant + correction
    affine = constant - float(profile @ q)
    scalar = float(row["three_band_linear_bl"]["price_log2"]) + correction
    error = abs(affine - scalar)
    if error > 1e-7:
        raise ValueError("corrected BL2 scalar-affine identity failed")
    return {
        "component_role": "outer", "family": "three-band-linear-bl2-graph-corrected-v1",
        "source_node_id": row["node_id"], "source_h2_cell": int(row["h2_cell"]),
        "source_barycenter_exact": row["barycenter_exact"],
        "parameters": {
            "log_variables": log_t.tolist(), "band_coefficients": coefficients.tolist(),
            "unpunctured_constant_log2": unpunctured_constant,
            "graph_replacement_log2": correction,
            "constant_log2": constant, "charge_log2": q.tolist(),
        },
        "source_identity": {
            "corrected_scalar_log2": scalar, "corrected_affine_log2": affine,
            "absolute_error_bits": error,
        },
    }


def source_artifact(probe_artifact: dict[str, Any]) -> dict[str, Any]:
    rows = [reconstruct_row(row) for row in probe_artifact["rows"]]
    return {
        "schema": SOURCE_SCHEMA, "status": "FROZEN_BINARY64_PARAMETERS_REQUIRES_OUTWARD_REPLAY",
        "source_probe_sha256": PROBE_SHA256,
        "graph_replacement_rule": "128*max_{abs(new-old)<=1}(q_new-q_old)-v1",
        "rows": rows,
    }


def selector(rows, references):
    return {"components": [{"weight": f"{weight.numerator}/{weight.denominator}", "component": references[index]} for index, weight in rows]}


def replay_leaf(node, inner, outer, inner_refs, outer_refs):
    matrix = np.asarray(node["exact_vertices"], dtype=np.float64)
    normalization = coordinate.diagnostic_normalization(matrix, 262144)
    iv = inner[0][None, :] - matrix @ inner[1].T
    ov = outer[0][None, :] - matrix @ outer[1].T
    mixed = mixtures.independent_minimax(iv, ov, normalization, DENOMINATOR)
    return {
        "score": float(mixed["score"]),
        "inner_selector": selector(mixed["inner"], inner_refs),
        "outer_selector": selector(mixed["outer"], outer_refs),
        "new_outer_positive": any(index >= len(outer_refs) - 8 for index, _weight in mixed["outer"]),
        "new_outer_rows": [index - (len(outer_refs) - 8) for index, _weight in mixed["outer"] if index >= len(outer_refs) - 8],
    }


def build(args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, Any]]:
    started = time.perf_counter()
    if coordinate.sha256_path(args.probe) != PROBE_SHA256:
        raise ValueError("frozen BL2 probe changed")
    geometry = json.loads(args.geometry.read_bytes())
    probe_artifact = json.loads(args.probe.read_bytes())
    frozen = source_artifact(probe_artifact)
    source_digest = engine.atomic_write(args.source, frozen)
    checkpoint = json.loads(geometry_wave.FACTORIZED.read_bytes())
    base, bank_bindings = mixtures.load_component_bank(args.manifest, args.atlas, args.old_supplementary, args.new_supplementary)
    full = geometry_wave.build_bank(base, checkpoint["candidates"])
    v2_inner = json.loads(args.v2_inner.read_bytes())
    v2_outer = json.loads(args.v2_outer.read_bytes())
    inner_c, inner_q, inner_refs = select_v2(full, v2_inner, "inner")
    outer_c, outer_q, outer_refs = select_v2(full, v2_outer, "outer")
    new_c = np.asarray([row["parameters"]["constant_log2"] for row in frozen["rows"]])
    new_q = np.asarray([row["parameters"]["charge_log2"] for row in frozen["rows"]])
    new_refs = tuple({
        "source_id": "g8-corrected-bl2-outer-v1", "source_sha256": source_digest,
        "row": index, "row_sha256": coordinate.sha256_bytes(coordinate.canonical_bytes(row)),
    } for index, row in enumerate(frozen["rows"]))
    active = sorted((node for node in geometry["nodes"].values() if node["state"] == "ACTIVE_LEAF"), key=lambda node: node["node_id"])
    rows = []
    for node in active:
        baseline = replay_leaf(node, (inner_c, inner_q), (outer_c, outer_q), inner_refs, outer_refs)
        enlarged_refs = (*outer_refs, *new_refs)
        enlarged = replay_leaf(node, (inner_c, inner_q), (np.concatenate((outer_c, new_c)), np.vstack((outer_q, new_q))), inner_refs, enlarged_refs)
        chosen = enlarged if enlarged["score"] <= baseline["score"] else baseline
        count = int(node["exact_count"])
        rows.append({
            "node_id": node["node_id"], "h2_cell": int(node["h2_cell"]),
            "exact_count": node["exact_count"], "exact_vertex_sha256": node["exact_vertex_sha256"],
            "old_upper_log2": baseline["score"], "candidate_upper_log2": chosen["score"],
            "old_contribution_log2": math.log2(count) + baseline["score"],
            "improvement_bits": baseline["score"] - chosen["score"],
            "contribution_log2": math.log2(count) + chosen["score"],
            "inner_selector": chosen["inner_selector"], "outer_selector": chosen["outer_selector"],
            "new_outer_positive": chosen["new_outer_positive"], "new_outer_rows": chosen["new_outer_rows"],
            "fallback_used": chosen is baseline,
        })
    old_aggregate = engine.log2_sum([row["old_contribution_log2"] for row in rows])
    new_aggregate = engine.log2_sum([row["contribution_log2"] for row in rows])
    old_worst = max(rows, key=lambda row: (row["old_contribution_log2"], row["node_id"]))
    new_worst = max(rows, key=lambda row: (row["contribution_log2"], row["node_id"]))
    used = sorted({index for row in rows for index in row["new_outer_rows"]})
    report = {
        "schema": SCHEMA, "status": "DIAGNOSTIC_RATIONAL_REPLAY_REQUIRES_OUTWARD_REPLAY",
        "source_bindings": {**bank_bindings, "geometry_sha256": coordinate.sha256_path(args.geometry), "probe_sha256": PROBE_SHA256, "v2_inner_sha256": coordinate.sha256_path(args.v2_inner), "v2_outer_sha256": coordinate.sha256_path(args.v2_outer), "corrected_bl2_source_sha256": source_digest},
        "configuration": {"inner_components": len(inner_refs), "old_outer_components": len(outer_refs), "new_outer_components": 8, "mixture_denominator": DENOMINATOR, "old_selector_fallback": True},
        "leaves": rows,
        "summary": {
            "active_leaf_count": len(rows), "old_aggregate_log2": old_aggregate,
            "aggregate_log2": new_aggregate, "aggregate_improvement_bits": old_aggregate - new_aggregate,
            "old_worst_leaf_id": old_worst["node_id"], "old_worst_log2": old_worst["old_contribution_log2"],
            "worst_leaf_id": new_worst["node_id"], "worst_log2": new_worst["contribution_log2"],
            "worst_improvement_bits": old_worst["old_contribution_log2"] - new_worst["contribution_log2"],
            "new_rows_used": used, "new_rows_used_count": len(used),
            "leaves_using_new_rows": sum(row["new_outer_positive"] for row in rows),
            "fallback_leaf_count": sum(row["fallback_used"] for row in rows),
            "no_regression": all(row["improvement_bits"] >= -1e-8 for row in rows),
            "residual_closes_minus40": new_aggregate <= -40.0,
        },
        "runtime_seconds": time.perf_counter() - started,
        "scope_limit": "Binary64 component materialization and rational replay; requires independent outward replay.",
    }
    return frozen, report


def self_test(args):
    frozen = source_artifact(json.loads(args.probe.read_bytes()))
    errors = [row["source_identity"]["absolute_error_bits"] for row in frozen["rows"]]
    corrections = [row["parameters"]["graph_replacement_log2"] for row in frozen["rows"]]
    if len(errors) != 8 or max(errors) > 1e-7 or min(corrections) < 0:
        raise SystemExit("corrected BL2 materialization self-test failed")
    print(f"maximum_source_identity_error_bits={max(errors):.12g}")
    print("graph_orientation_new_minus_old=PASS")
    print("status=PASS_G8_CORRECTED_BL2_COMPONENT_SELF_TEST")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, default=PROBE); parser.add_argument("--geometry", type=Path, default=GEOMETRY)
    parser.add_argument("--v2-inner", type=Path, default=V2_INNER); parser.add_argument("--v2-outer", type=Path, default=V2_OUTER)
    parser.add_argument("--manifest", type=Path, default=mixtures.h2.DEFAULT_MANIFEST); parser.add_argument("--atlas", type=Path, default=mixtures.h2.DEFAULT_ATLAS)
    parser.add_argument("--old-supplementary", type=Path, default=mixtures.h2.DEFAULT_SUPPLEMENTARY); parser.add_argument("--new-supplementary", type=Path, default=mixtures.NEW_ATLAS)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE); parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT); parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test: self_test(args); return
    _source, report = build(args); digest = engine.atomic_write(args.output, report)
    print(f"aggregate_improvement_bits={report['summary']['aggregate_improvement_bits']:.12f}")
    print(f"aggregate_log2={report['summary']['aggregate_log2']:.12f}")
    print(f"worst_log2={report['summary']['worst_log2']:.12f}")
    print(f"residual_closes={report['summary']['residual_closes_minus40']}")
    print(f"runtime_seconds={report['runtime_seconds']:.3f}")
    print(f"output={args.output}"); print(f"sha256={digest}")


if __name__ == "__main__": main()
