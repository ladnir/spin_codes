#!/usr/bin/env python3
"""Run one fresh corrected BL2 plus total-spectrum outer column round."""

from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import extract_packet_group_g8_component_pricing_target as pricing
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import probe_packet_group_g8_outer_family_prices as family_probe
import replay_packet_group_g8_corrected_bl2_components as corrected
import replay_packet_group_g8_corrected_bl2_total_spectrum as combined
import replay_packet_group_g8_independent_component_mixtures as mixtures
import run_packet_group_g8_adaptive_cumulative_engine as engine
import run_packet_group_g8_factorized_geometry_persistence_wave as geometry_wave


ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "out" / "g8_corrected_bl2_total_spectrum_replay_197.json"
BASELINE_SHA256 = "bc7235bf483218ea645a8eb128f837fa0bc2417a5c938fe001f5267a16332233"
DEFAULT_BL2_SOURCE = ROOT / "out" / "g8_fresh_round_bl2_outer_components.json"
DEFAULT_TOTAL_SOURCE = ROOT / "out" / "g8_fresh_round_total_spectrum_outer_components.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_fresh_outer_family_round_replay_197.json"
SOURCE_SCHEMA = "permute-conv.packet-group-g8-fresh-outer-family-components.v1"
SCHEMA = "permute-conv.packet-group-g8-fresh-outer-family-round.v1"
TOP = 8
GATE = 65536.0


def source_arrays(source, source_id, digest):
    return combined.source_arrays(source, source_id, digest)


def current_bank(args):
    checkpoint = json.loads(geometry_wave.FACTORIZED.read_bytes())
    base, bindings = mixtures.load_component_bank(
        args.manifest, args.atlas, args.old_supplementary, args.new_supplementary
    )
    full = geometry_wave.build_bank(base, checkpoint["candidates"])
    inner = corrected.select_v2(full, json.loads(args.v2_inner.read_bytes()), "inner")
    outer = corrected.select_v2(full, json.loads(args.v2_outer.read_bytes()), "outer")
    bl2_source = json.loads(args.corrected_bl2.read_bytes())
    total_source = json.loads(args.corrected_total.read_bytes())
    bc, bq, br = source_arrays(bl2_source, "g8-corrected-bl2-outer-v1", combined.BL2_SOURCE_SHA256)
    total_digest = coordinate.sha256_path(args.corrected_total)
    tc, tq, tr = source_arrays(total_source, "g8-corrected-total-spectrum-outer-v1", total_digest)
    outer = (np.concatenate((outer[0], bc, tc)), np.vstack((outer[1], bq, tq)), (*outer[2], *br, *tr))
    return inner, outer, bindings


def rational_dual(node, inner, outer):
    matrix = np.asarray(node["exact_vertices"], dtype=np.float64)
    normalization = coordinate.diagnostic_normalization(matrix, 262144)
    iv = inner[0][None, :] - matrix @ inner[1].T
    ov = outer[0][None, :] - matrix @ outer[1].T
    lp = pricing.solve_component_lp(iv, ov, normalization)
    return pricing.rationalize_distribution(
        -np.asarray(lp.ineqlin.marginals, dtype=np.float64),
        geometry_wave.factorized.one_shot.DUAL_DENOMINATOR,
    )


def barycenter(node, dual):
    profiles = tuple(tuple(map(int, row)) for row in node["exact_vertices"])
    return pricing.exact_barycenter(profiles, dual)


def optimized_row(node_id, h2_cell, exact_barycenter, family):
    profile = np.asarray([float(value) for value in exact_barycenter])
    if family == "bl2":
        value, result = family_probe.three_band.optimize_outer(
            profile.tolist(), optimize_band_coefficients=True, log_fugacity_bound=40.0
        )
        q = np.concatenate(([0.0], result.x[:8])) / math.log(2.0)
        p1 = float(result.x[8])
        raw_price = float(value / math.log(2.0))
        parameters = {"log_variables": (q * math.log(2.0)).tolist(), "band_coefficients": [1.0 - p1, p1, p1]}
        optimizer = result
        name = "three-band-linear-bl2-graph-corrected-v1"
    else:
        value, result, parts = family_probe.total_spectrum.optimize_outer(profile.tolist())
        q = np.concatenate(([0.0], result.x[:8])) / math.log(2.0)
        raw_price = float(value / math.log(2.0))
        parameters = {
            "log_variables": (q * math.log(2.0)).tolist(), "log_beta": float(result.x[8]),
            "log_alpha": float(parts[1]), "log_local_enumerator": float(parts[2]),
        }
        optimizer = result
        name = "exact-total-spectrum-envelope-graph-corrected-v1"
    graph = corrected.graph_correction(q)
    constant = raw_price + float(profile @ q) + graph
    affine = constant - float(profile @ q)
    scalar = raw_price + graph
    return {
        "component_role": "outer", "family": name, "source_node_id": node_id,
        "source_h2_cell": h2_cell,
        "source_barycenter_exact": [pricing.fraction_text(value) for value in exact_barycenter],
        "parameters": {**parameters, "graph_replacement_log2": graph, "constant_log2": constant, "charge_log2": q.tolist()},
        "source_identity": {"corrected_scalar_log2": scalar, "corrected_affine_log2": affine, "absolute_error_bits": abs(scalar - affine)},
        "optimizer": {"success": bool(optimizer.success), "iterations": int(optimizer.nit), "evaluations": int(optimizer.nfev), "message": str(optimizer.message)},
    }


def source_artifact(family, rows):
    return {
        "schema": SOURCE_SCHEMA, "status": "FROZEN_BINARY64_PARAMETERS_REQUIRES_OUTWARD_REPLAY",
        "family": family, "baseline_sha256": BASELINE_SHA256,
        "graph_replacement_rule": "128*max_{abs(new-old)<=1}(q_new-q_old)-v1",
        "rows": rows,
    }


def build(args):
    started = time.perf_counter()
    if coordinate.sha256_path(args.baseline) != BASELINE_SHA256:
        raise ValueError("combined baseline changed")
    baseline_artifact = json.loads(args.baseline.read_bytes())
    geometry = json.loads(args.geometry.read_bytes())
    nodes = {node["node_id"]: node for node in geometry["nodes"].values() if node["state"] == "ACTIVE_LEAF"}
    inner, outer, bindings = current_bank(args)
    baseline_rows = []
    for stored in baseline_artifact["leaves"]:
        node = nodes[stored["node_id"]]
        replay = corrected.replay_leaf(node, (inner[0], inner[1]), (outer[0], outer[1]), inner[2], outer[2])
        dual = rational_dual(node, inner, outer)
        contribution = math.log2(int(node["exact_count"])) + replay["score"]
        baseline_rows.append({"node": node, "replay": replay, "dual": dual, "contribution": contribution})
        if abs(contribution - float(stored["contribution_log2"])) > 1e-7:
            raise ValueError("current selector replay differs from combined baseline")
    ranked = sorted(baseline_rows, key=lambda row: (-row["contribution"], row["node"]["node_id"]))
    targets, seen = [], set()
    for row in ranked:
        cell = int(row["node"]["h2_cell"])
        if cell in seen: continue
        seen.add(cell); targets.append(row)
        if len(targets) == TOP: break
    candidates = {"bl2": [], "total_spectrum": []}
    pricing_rows = []
    for target in targets:
        node = target["node"]
        center = barycenter(node, target["dual"])
        profile = np.asarray([float(value) for value in center])
        mu = float(np.min(outer[0] - outer[1] @ profile))
        for family in ("bl2", "total_spectrum"):
            row = optimized_row(node["node_id"], int(node["h2_cell"]), center, family)
            price = row["source_identity"]["corrected_affine_log2"]
            reduced = price - mu
            row["pricing"] = {"current_outer_mu_log2": mu, "candidate_price_log2": price, "reduced_cost_bits": reduced}
            row["accepted"] = row["optimizer"]["success"] and row["source_identity"]["absolute_error_bits"] <= 1e-7 and reduced < 0.0
            pricing_rows.append({"node_id": node["node_id"], "h2_cell": int(node["h2_cell"]), "family": family, **row["pricing"], "accepted": row["accepted"]})
            if row["accepted"]: candidates[family].append(row)
    bl2_source = source_artifact("bl2", candidates["bl2"])
    total_source = source_artifact("total_spectrum", candidates["total_spectrum"])
    bl2_digest = engine.atomic_write(args.bl2_source, bl2_source)
    total_digest = engine.atomic_write(args.total_source, total_source)
    bc, bq, br = source_arrays(bl2_source, "g8-fresh-round-bl2-v1", bl2_digest)
    tc, tq, tr = source_arrays(total_source, "g8-fresh-round-total-spectrum-v1", total_digest)
    enlarged = (np.concatenate((outer[0], bc, tc)), np.vstack((outer[1], bq, tq)), (*outer[2], *br, *tr))
    leaves = []
    for base_row in baseline_rows:
        node, old = base_row["node"], base_row["replay"]
        candidate = corrected.replay_leaf(node, (inner[0], inner[1]), (enlarged[0], enlarged[1]), inner[2], enlarged[2])
        chosen = candidate if candidate["score"] <= old["score"] else old
        count_log = math.log2(int(node["exact_count"]))
        used = [component["component"] for component in chosen["outer_selector"]["components"] if component["component"]["source_sha256"] in (bl2_digest, total_digest)] if chosen is candidate else []
        leaves.append({"node_id": node["node_id"], "h2_cell": int(node["h2_cell"]), "exact_count": node["exact_count"], "exact_vertex_sha256": node["exact_vertex_sha256"], "old_upper_log2": old["score"], "candidate_upper_log2": chosen["score"], "old_contribution_log2": count_log + old["score"], "contribution_log2": count_log + chosen["score"], "improvement_bits": old["score"] - chosen["score"], "inner_selector": chosen["inner_selector"], "outer_selector": chosen["outer_selector"], "new_components": used, "fallback_used": chosen is old})
    old_aggregate = engine.log2_sum([row["old_contribution_log2"] for row in leaves]); aggregate = engine.log2_sum([row["contribution_log2"] for row in leaves])
    old_worst = max(leaves, key=lambda row: (row["old_contribution_log2"], row["node_id"])); worst = max(leaves, key=lambda row: (row["contribution_log2"], row["node_id"]))
    report = {
        "schema": SCHEMA, "status": "DIAGNOSTIC_RATIONAL_REPLAY_REQUIRES_OUTWARD_REPLAY",
        "source_bindings": {**bindings, "baseline_sha256": BASELINE_SHA256, "fresh_bl2_source_sha256": bl2_digest, "fresh_total_spectrum_source_sha256": total_digest},
        "configuration": {"targets": TOP, "jobs": 16, "serialized": True, "old_outer_components": len(outer[2]), "accepted_bl2": len(bc), "accepted_total_spectrum": len(tc), "mixture_denominator": corrected.DENOMINATOR},
        "pricing": pricing_rows, "leaves": leaves,
        "summary": {"active_leaf_count": len(leaves), "old_aggregate_log2": old_aggregate, "aggregate_log2": aggregate, "aggregate_improvement_bits": old_aggregate - aggregate, "old_worst_leaf_id": old_worst["node_id"], "old_worst_log2": old_worst["old_contribution_log2"], "worst_leaf_id": worst["node_id"], "worst_log2": worst["contribution_log2"], "worst_improvement_bits": old_worst["old_contribution_log2"] - worst["contribution_log2"], "leaves_using_new": sum(bool(row["new_components"]) for row in leaves), "new_rows_used": sorted({(ref["source_sha256"], ref["row"]) for row in leaves for ref in row["new_components"]}), "fallback_leaf_count": sum(row["fallback_used"] for row in leaves), "no_regression": all(row["improvement_bits"] >= -1e-8 for row in leaves), "improvement_gate_bits": GATE, "improvement_gate_passes": old_aggregate - aggregate >= GATE, "residual_closes_minus40": aggregate <= -40.0},
        "runtime_seconds": time.perf_counter() - started,
        "scope_limit": "Binary64 discovery and rational replay; requires independent outward replay.",
    }
    return report


def self_test(args):
    if coordinate.sha256_path(args.baseline) != BASELINE_SHA256: raise SystemExit("baseline binding failed")
    q = np.asarray([0., 1., 0., 0., 0., 0., 0., 0., 0.]); expected = 128.0
    if corrected.graph_correction(q) != expected: raise SystemExit("graph correction synthetic failed")
    print("combined_baseline_binding=PASS"); print("ordered_adjacent_graph_correction=PASS"); print("status=PASS_G8_FRESH_OUTER_FAMILY_ROUND_SELF_TEST")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=BASELINE); parser.add_argument("--geometry", type=Path, default=corrected.GEOMETRY)
    parser.add_argument("--v2-inner", type=Path, default=corrected.V2_INNER); parser.add_argument("--v2-outer", type=Path, default=corrected.V2_OUTER)
    parser.add_argument("--corrected-bl2", type=Path, default=combined.BL2_SOURCE); parser.add_argument("--corrected-total", type=Path, default=combined.DEFAULT_SOURCE)
    parser.add_argument("--manifest", type=Path, default=mixtures.h2.DEFAULT_MANIFEST); parser.add_argument("--atlas", type=Path, default=mixtures.h2.DEFAULT_ATLAS); parser.add_argument("--old-supplementary", type=Path, default=mixtures.h2.DEFAULT_SUPPLEMENTARY); parser.add_argument("--new-supplementary", type=Path, default=mixtures.NEW_ATLAS)
    parser.add_argument("--bl2-source", type=Path, default=DEFAULT_BL2_SOURCE); parser.add_argument("--total-source", type=Path, default=DEFAULT_TOTAL_SOURCE); parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT); parser.add_argument("--self-test", action="store_true"); args=parser.parse_args()
    if args.self_test: self_test(args); return
    report=build(args); digest=engine.atomic_write(args.output,report); s=report["summary"]
    print(f"aggregate_improvement_bits={s['aggregate_improvement_bits']:.12f}"); print(f"aggregate_log2={s['aggregate_log2']:.12f}"); print(f"gate_passes={s['improvement_gate_passes']}"); print(f"residual_closes={s['residual_closes_minus40']}"); print(f"runtime_seconds={report['runtime_seconds']:.3f}"); print(f"output={args.output}"); print(f"sha256={digest}")


if __name__ == "__main__": main()
