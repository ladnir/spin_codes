#!/usr/bin/env python3
"""Replay corrected BL2 and exact-total-spectrum outer columns on 197 leaves."""

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
import replay_packet_group_g8_corrected_bl2_components as bl2_replay
import replay_packet_group_g8_independent_component_mixtures as mixtures
import run_packet_group_g8_adaptive_cumulative_engine as engine
import run_packet_group_g8_factorized_geometry_persistence_wave as geometry_wave


ROOT = Path(__file__).resolve().parents[1]
BL2_SOURCE = bl2_replay.DEFAULT_SOURCE
BL2_SOURCE_SHA256 = "cbc91782e8b137459e888003a75bc16217af5250951693c3e3ed1e67000064a5"
DEFAULT_SOURCE = ROOT / "out" / "g8_corrected_total_spectrum_outer_components.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_corrected_bl2_total_spectrum_replay_197.json"
SOURCE_SCHEMA = "permute-conv.packet-group-g8-corrected-total-spectrum-outer-components.v1"
SCHEMA = "permute-conv.packet-group-g8-corrected-bl2-total-spectrum-replay.v1"


def reconstruct_total_row(row: dict[str, Any]) -> dict[str, Any]:
    profile = np.asarray([float(Fraction(value)) for value in row["barycenter_exact"]])
    family = row["total_spectrum"]
    q = np.asarray(family["parameters"]["log_variables"], dtype=np.float64) / math.log(2.0)
    unpunctured_price = float(family["price_log2"])
    unpunctured_constant = unpunctured_price + float(profile @ q)
    correction = bl2_replay.graph_correction(q)
    constant = unpunctured_constant + correction
    scalar = unpunctured_price + correction
    affine = constant - float(profile @ q)
    error = abs(scalar - affine)
    if error > 1e-7:
        raise ValueError("corrected total-spectrum scalar-affine identity failed")
    return {
        "component_role": "outer",
        "family": "exact-total-spectrum-envelope-graph-corrected-v1",
        "source_node_id": row["node_id"], "source_h2_cell": int(row["h2_cell"]),
        "source_barycenter_exact": row["barycenter_exact"],
        "parameters": {
            "log_variables": family["parameters"]["log_variables"],
            "log_beta": family["parameters"]["log_beta"],
            "log_alpha": family["parameters"]["log_alpha"],
            "log_local_enumerator": family["parameters"]["log_local_enumerator"],
            "active_envelope_weights": family["parameters"]["active_envelope_weights"],
            "unpunctured_constant_log2": unpunctured_constant,
            "graph_replacement_log2": correction,
            "constant_log2": constant, "charge_log2": q.tolist(),
        },
        "source_identity": {
            "corrected_scalar_log2": scalar, "corrected_affine_log2": affine,
            "absolute_error_bits": error,
        },
    }


def total_source(probe: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": SOURCE_SCHEMA,
        "status": "FROZEN_BINARY64_PARAMETERS_REQUIRES_OUTWARD_REPLAY",
        "source_probe_sha256": bl2_replay.PROBE_SHA256,
        "graph_replacement_rule": "128*max_{abs(new-old)<=1}(q_new-q_old)-v1",
        "rows": [reconstruct_total_row(row) for row in probe["rows"]],
    }


def source_arrays(source: dict[str, Any], source_id: str, digest: str):
    constants = np.asarray([row["parameters"]["constant_log2"] for row in source["rows"]])
    charges = np.asarray([row["parameters"]["charge_log2"] for row in source["rows"]])
    references = tuple({
        "source_id": source_id, "source_sha256": digest, "row": index,
        "row_sha256": coordinate.sha256_bytes(coordinate.canonical_bytes(row)),
    } for index, row in enumerate(source["rows"]))
    return constants, charges, references


def build(args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, Any]]:
    started = time.perf_counter()
    if coordinate.sha256_path(args.probe) != bl2_replay.PROBE_SHA256:
        raise ValueError("frozen outer-family probe changed")
    if coordinate.sha256_path(args.bl2_source) != BL2_SOURCE_SHA256:
        raise ValueError("corrected BL2 source changed")
    probe = json.loads(args.probe.read_bytes())
    total = total_source(probe)
    total_digest = engine.atomic_write(args.total_source, total)
    bl2_source = json.loads(args.bl2_source.read_bytes())
    geometry = json.loads(args.geometry.read_bytes())
    checkpoint = json.loads(geometry_wave.FACTORIZED.read_bytes())
    base, bank_bindings = mixtures.load_component_bank(
        args.manifest, args.atlas, args.old_supplementary, args.new_supplementary
    )
    full = geometry_wave.build_bank(base, checkpoint["candidates"])
    inner = bl2_replay.select_v2(full, json.loads(args.v2_inner.read_bytes()), "inner")
    old_outer = bl2_replay.select_v2(full, json.loads(args.v2_outer.read_bytes()), "outer")
    bc, bq, br = source_arrays(bl2_source, "g8-corrected-bl2-outer-v1", BL2_SOURCE_SHA256)
    tc, tq, tr = source_arrays(total, "g8-corrected-total-spectrum-outer-v1", total_digest)
    enlarged = (
        np.concatenate((old_outer[0], bc, tc)),
        np.vstack((old_outer[1], bq, tq)),
        (*old_outer[2], *br, *tr),
    )
    active = sorted(
        (node for node in geometry["nodes"].values() if node["state"] == "ACTIVE_LEAF"),
        key=lambda node: node["node_id"],
    )
    rows = []
    for node in active:
        baseline = bl2_replay.replay_leaf(node, (inner[0], inner[1]), (old_outer[0], old_outer[1]), inner[2], old_outer[2])
        candidate = bl2_replay.replay_leaf(node, (inner[0], inner[1]), (enlarged[0], enlarged[1]), inner[2], enlarged[2])
        chosen = candidate if candidate["score"] <= baseline["score"] else baseline
        count_log = math.log2(int(node["exact_count"]))
        used = []
        if chosen is candidate:
            for component in chosen["outer_selector"]["components"]:
                ref = component["component"]
                if ref["source_sha256"] in (BL2_SOURCE_SHA256, total_digest):
                    used.append({"family": "bl2" if ref["source_sha256"] == BL2_SOURCE_SHA256 else "total_spectrum", "row": int(ref["row"])})
        rows.append({
            "node_id": node["node_id"], "h2_cell": int(node["h2_cell"]),
            "exact_count": node["exact_count"], "exact_vertex_sha256": node["exact_vertex_sha256"],
            "old_upper_log2": baseline["score"], "candidate_upper_log2": chosen["score"],
            "old_contribution_log2": count_log + baseline["score"],
            "contribution_log2": count_log + chosen["score"],
            "improvement_bits": baseline["score"] - chosen["score"],
            "inner_selector": chosen["inner_selector"], "outer_selector": chosen["outer_selector"],
            "new_components": used, "fallback_used": chosen is baseline,
        })
    old_aggregate = engine.log2_sum([row["old_contribution_log2"] for row in rows])
    aggregate = engine.log2_sum([row["contribution_log2"] for row in rows])
    old_worst = max(rows, key=lambda row: (row["old_contribution_log2"], row["node_id"]))
    worst = max(rows, key=lambda row: (row["contribution_log2"], row["node_id"]))
    bl2_used = sorted({item["row"] for row in rows for item in row["new_components"] if item["family"] == "bl2"})
    total_used = sorted({item["row"] for row in rows for item in row["new_components"] if item["family"] == "total_spectrum"})
    report = {
        "schema": SCHEMA, "status": "DIAGNOSTIC_RATIONAL_REPLAY_REQUIRES_OUTWARD_REPLAY",
        "source_bindings": {**bank_bindings, "geometry_sha256": coordinate.sha256_path(args.geometry), "v2_inner_sha256": coordinate.sha256_path(args.v2_inner), "v2_outer_sha256": coordinate.sha256_path(args.v2_outer), "bl2_source_sha256": BL2_SOURCE_SHA256, "total_spectrum_source_sha256": total_digest},
        "configuration": {"inner_components": len(inner[2]), "v2_outer_components": len(old_outer[2]), "corrected_bl2_components": 8, "corrected_total_spectrum_components": 8, "mixture_denominator": bl2_replay.DENOMINATOR, "old_selector_fallback": True},
        "leaves": rows,
        "summary": {
            "active_leaf_count": len(rows), "old_aggregate_log2": old_aggregate,
            "aggregate_log2": aggregate, "aggregate_improvement_bits": old_aggregate - aggregate,
            "old_worst_leaf_id": old_worst["node_id"], "old_worst_log2": old_worst["old_contribution_log2"],
            "worst_leaf_id": worst["node_id"], "worst_log2": worst["contribution_log2"],
            "worst_improvement_bits": old_worst["old_contribution_log2"] - worst["contribution_log2"],
            "bl2_rows_used": bl2_used, "total_spectrum_rows_used": total_used,
            "leaves_using_bl2": sum(any(item["family"] == "bl2" for item in row["new_components"]) for row in rows),
            "leaves_using_total_spectrum": sum(any(item["family"] == "total_spectrum" for item in row["new_components"]) for row in rows),
            "fallback_leaf_count": sum(row["fallback_used"] for row in rows),
            "no_regression": all(row["improvement_bits"] >= -1e-8 for row in rows),
            "residual_closes_minus40": aggregate <= -40.0,
        },
        "runtime_seconds": time.perf_counter() - started,
        "scope_limit": "Binary64 materialization and rational replay; requires independent outward replay.",
    }
    return total, report


def self_test(args):
    total = total_source(json.loads(args.probe.read_bytes()))
    errors = [row["source_identity"]["absolute_error_bits"] for row in total["rows"]]
    corrections = [row["parameters"]["graph_replacement_log2"] for row in total["rows"]]
    if len(errors) != 8 or max(errors) > 1e-7 or min(corrections) < 0:
        raise SystemExit("corrected total-spectrum source self-test failed")
    print(f"maximum_total_spectrum_identity_error_bits={max(errors):.12g}")
    print("ordered_adjacent_graph_orientation=PASS")
    print("status=PASS_G8_CORRECTED_TOTAL_SPECTRUM_SELF_TEST")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, default=bl2_replay.PROBE); parser.add_argument("--geometry", type=Path, default=bl2_replay.GEOMETRY)
    parser.add_argument("--v2-inner", type=Path, default=bl2_replay.V2_INNER); parser.add_argument("--v2-outer", type=Path, default=bl2_replay.V2_OUTER)
    parser.add_argument("--bl2-source", type=Path, default=BL2_SOURCE); parser.add_argument("--total-source", type=Path, default=DEFAULT_SOURCE); parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--manifest", type=Path, default=mixtures.h2.DEFAULT_MANIFEST); parser.add_argument("--atlas", type=Path, default=mixtures.h2.DEFAULT_ATLAS)
    parser.add_argument("--old-supplementary", type=Path, default=mixtures.h2.DEFAULT_SUPPLEMENTARY); parser.add_argument("--new-supplementary", type=Path, default=mixtures.NEW_ATLAS)
    parser.add_argument("--self-test", action="store_true"); args = parser.parse_args()
    if args.self_test: self_test(args); return
    _source, report = build(args); digest = engine.atomic_write(args.output, report)
    print(f"aggregate_improvement_bits={report['summary']['aggregate_improvement_bits']:.12f}")
    print(f"aggregate_log2={report['summary']['aggregate_log2']:.12f}")
    print(f"worst_log2={report['summary']['worst_log2']:.12f}")
    print(f"residual_closes={report['summary']['residual_closes_minus40']}")
    print(f"runtime_seconds={report['runtime_seconds']:.3f}")
    print(f"output={args.output}"); print(f"sha256={digest}")


if __name__ == "__main__": main()
