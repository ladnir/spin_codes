#!/usr/bin/env python3
"""Run one exact geometry wave with the frozen enlarged factorized bank."""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

for _variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_variable] = "1"

import numpy as np

import extract_packet_group_g8_component_pricing_target as pricing
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import replay_packet_group_g8_adaptive_checkpoint_with_atlas as stored
import replay_packet_group_g8_independent_component_mixtures as component
import run_packet_group_g8_adaptive_cumulative_engine as engine
import run_packet_group_g8_adaptive_independent_persistence_wave as persistence
import run_packet_group_g8_iterative_factorized_columns as factorized


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "out" / "g8_full_support_adaptive_independent_persistence_wave.json"
GEOMETRY_SHA256 = "32a7eb88d09f20725a8fd4fe5a727c6f3bc925ee9e90febcc753c2ed9cfa0739"
FACTORIZED = ROOT / "out" / "g8_iterative_factorized_columns" / "checkpoint.json"
FACTORIZED_SHA256 = "359dcd1cb4a048816d360872593df23cd957e9063198dcfc304e7e8e37bb8d27"
ROUND1_SHA256 = "edc9c8b1714024d99af4c8d97352a8d20ec5f880202962667030d1ea2d9ae6f3"
ROUND2_SHA256 = "b6d324477840ea5708776c777b45ca0f3272ea3c3e9a217729a1984301d37250"
DEFAULT_OUTPUT = ROOT / "out" / "g8_factorized_geometry_persistence_wave.json"
SCHEMA = "permute-conv.packet-group-g8-factorized-geometry-persistence-wave.v1"
DENOMINATOR = 1 << 30
TOP_LEAVES = 8
MINIMUM_AGGREGATE_IMPROVEMENT = 65536.0


@dataclass(frozen=True)
class FactorizedBank:
    inner_constants: np.ndarray
    inner_charges: np.ndarray
    inner_references: tuple[dict[str, Any], ...]
    outer_constants: np.ndarray
    outer_charges: np.ndarray
    outer_references: tuple[dict[str, Any], ...]


def build_bank(base: component.ComponentBank, candidates: list[dict[str, Any]]) -> FactorizedBank:
    ic, iz, ir = factorized.pool(base, candidates, "inner")
    oc, oz, outer_references = factorized.pool(base, candidates, "outer")
    return FactorizedBank(ic, iz, ir, oc, oz, outer_references)


def score_bounds(
    lower: tuple[int, ...], upper: tuple[int, ...], bank: FactorizedBank,
    node_id: str, h2_cell: int, depth: int,
) -> dict[str, Any] | None:
    vertices = persistence.geometry.chain_vertices(lower, upper)
    if not vertices or len(vertices) > persistence.geometry.VERTEX_CAP_PER_CHILD:
        return None
    count = persistence.geometry.exact_chain_count(lower, upper)
    if count <= 0:
        return None
    profiles = [list(persistence.geometry.profile_from_chain(vertex)) for vertex in vertices]
    matrix = np.asarray(profiles, dtype=np.float64)
    normalization = coordinate.diagnostic_normalization(matrix, 262144)
    inner_values = bank.inner_constants[None, :] - matrix @ bank.inner_charges.T
    outer_values = bank.outer_constants[None, :] - matrix @ bank.outer_charges.T
    mixed = component.independent_minimax(
        inner_values, outer_values, normalization, DENOMINATOR
    )
    lp = pricing.solve_component_lp(inner_values, outer_values, normalization)
    dual = pricing.rationalize_distribution(
        -np.asarray(lp.ineqlin.marginals, dtype=np.float64),
        factorized.one_shot.DUAL_DENOMINATOR,
    )
    score = float(mixed["score"])
    digest = coordinate.sha256_bytes(coordinate.canonical_bytes(profiles))
    worst = int(np.argmax(np.asarray(mixed["vertex_values"])))
    return {
        "lower": list(lower), "upper": list(upper), "exact_count": str(count),
        "count": {
            "kind": "exact", "method": "bounded-monotone-chain-prefix-dp-v1",
            "parameters": {"lower": list(lower), "upper": list(upper), "free_mass": str(262144 - 9)},
            "value": str(count),
        },
        "exact_vertex_count": len(profiles), "exact_vertices": profiles,
        "exact_vertex_sha256": digest, "candidate_upper_log2": score,
        "contribution_log2": math.log2(count) + score, "worst_vertex": profiles[worst],
        "selector": {
            "kind": "independent-components",
            "inner": factorized.selector(mixed["inner"], bank.inner_references),
            "outer": factorized.selector(mixed["outer"], bank.outer_references),
        },
        "inner_support_size": len(mixed["inner"]),
        "outer_support_size": len(mixed["outer"]),
        "component_count": len(mixed["inner"]) + len(mixed["outer"]),
        "vertex_dual": [
            {"vertex": index, "weight": pricing.fraction_text(weight)}
            for index, weight in dual
        ],
        "candidate_status": "DIAGNOSTIC_BINARY64_REQUIRES_OUTWARD_REPLAY",
        "diagnostic": {
            "binary64_lp_upper_log2": float(lp.fun),
            "rationalization_loss_bits": float(mixed["rationalization_loss_bits"]),
        },
    }


def leading_pairs(vertices: tuple[tuple[int, ...], ...], bank: FactorizedBank):
    matrix = np.asarray(
        [persistence.geometry.profile_from_chain(vertex) for vertex in vertices],
        dtype=np.float64,
    )
    inner = bank.inner_constants[None, :] - matrix @ bank.inner_charges.T
    outer = bank.outer_constants[None, :] - matrix @ bank.outer_charges.T
    return tuple(
        (int(np.argmin(inner[row])), int(np.argmin(outer[row])))
        for row in range(len(vertices))
    )


def transition_candidates(vertices: tuple[tuple[int, ...], ...], bank: FactorizedBank):
    pairs = leading_pairs(vertices, bank)
    profiles = tuple(persistence.geometry.profile_from_chain(vertex) for vertex in vertices)
    proposals = {}
    for k in range(8):
        groups = {}
        for index, vertex in enumerate(vertices):
            groups.setdefault(vertex[:k] + vertex[k + 1:], []).append(index)
        for indices in groups.values():
            indices.sort(key=lambda index: (vertices[index][k], index))
            for left_index, right_index in zip(indices, indices[1:]):
                left_value, right_value = vertices[left_index][k], vertices[right_index][k]
                left_pair, right_pair = pairs[left_index], pairs[right_index]
                if left_value >= right_value or left_pair == right_pair:
                    continue
                lc = bank.inner_constants[left_pair[0]] + bank.outer_constants[left_pair[1]]
                rc = bank.inner_constants[right_pair[0]] + bank.outer_constants[right_pair[1]]
                lz = bank.inner_charges[left_pair[0]] + bank.outer_charges[left_pair[1]]
                rz = bank.inner_charges[right_pair[0]] + bank.outer_charges[right_pair[1]]
                dc, dz = float(lc - rc), lz - rz
                dl = dc - float(dz @ np.asarray(profiles[left_index]))
                dr = dc - float(dz @ np.asarray(profiles[right_index]))
                thresholds, crossing = set(), None
                if math.isfinite(dl) and math.isfinite(dr) and dl != dr and (
                    dl <= 0.0 <= dr or dr <= 0.0 <= dl
                ):
                    crossing = left_value + (-dl / (dr - dl)) * (right_value - left_value)
                    if math.isfinite(crossing):
                        thresholds.update((math.floor(crossing), math.ceil(crossing) - 1))
                if not thresholds:
                    thresholds.add((left_value + right_value) // 2)
                for threshold in thresholds:
                    if not left_value <= threshold < right_value:
                        continue
                    record = {
                        "coordinate": k, "threshold": str(threshold),
                        "left_leading_inner": bank.inner_references[left_pair[0]],
                        "left_leading_outer": bank.outer_references[left_pair[1]],
                        "right_leading_inner": bank.inner_references[right_pair[0]],
                        "right_leading_outer": bank.outer_references[right_pair[1]],
                        "crossing_coordinate_binary64": crossing,
                    }
                    key = (k, threshold)
                    if key not in proposals or coordinate.canonical_bytes(record) < coordinate.canonical_bytes(proposals[key]):
                        proposals[key] = record
    return proposals


def refine_leaf(node: dict[str, Any], bank: FactorizedBank) -> dict[str, Any]:
    lower, upper = tuple(map(int, node["lower"])), tuple(map(int, node["upper"]))
    parent = score_bounds(lower, upper, bank, node["node_id"], int(node["h2_cell"]), int(node["depth"]))
    if parent is None or parent["exact_count"] != node["exact_count"] or parent["exact_vertex_sha256"] != node["exact_vertex_sha256"]:
        raise RuntimeError("selected parent exact geometry changed")
    vertices = persistence.geometry.chain_vertices(lower, upper)
    proposals, candidates, cap_rejections = transition_candidates(vertices, bank), [], 0
    for (k, threshold), transition in sorted(proposals.items()):
        left_bounds = persistence.geometry.child_bounds(lower, upper, k, threshold, True)
        right_bounds = persistence.geometry.child_bounds(lower, upper, k, threshold, False)
        if left_bounds is None or right_bounds is None:
            continue
        if max(len(persistence.geometry.chain_vertices(*left_bounds)), len(persistence.geometry.chain_vertices(*right_bounds))) > persistence.geometry.VERTEX_CAP_PER_CHILD:
            cap_rejections += 1
            continue
        left = score_bounds(*left_bounds, bank, f"{node['node_id']}/L", int(node["h2_cell"]), int(node["depth"]) + 1)
        right = score_bounds(*right_bounds, bank, f"{node['node_id']}/R", int(node["h2_cell"]), int(node["depth"]) + 1)
        if left is None or right is None:
            continue
        if int(left["exact_count"]) + int(right["exact_count"]) != int(parent["exact_count"]):
            raise RuntimeError("child counts do not sum to parent")
        expanded = engine.log2_sum([left["contribution_log2"], right["contribution_log2"]])
        maximum = max(left["candidate_upper_log2"], right["candidate_upper_log2"])
        candidates.append({
            "key": (expanded, maximum, -min(int(left["exact_count"]), int(right["exact_count"])), left["exact_vertex_count"] + right["exact_vertex_count"], k, threshold),
            "improvement_bits": parent["contribution_log2"] - expanded,
            "expanded_child_contribution_log2": expanded,
            "collapsed_child_contribution_log2": math.log2(int(parent["exact_count"])) + maximum,
            "k": k, "t": threshold, "transition": transition, "left": left, "right": right,
        })
    positive = [row for row in candidates if row["improvement_bits"] > 0.0]
    base = {"parent_replay": parent, "transition_candidate_count": len(proposals), "evaluated_candidate_count": len(candidates), "child_cap_rejection_count": cap_rejections}
    if not positive:
        return {**base, "state": "NO_VALID_SPLIT"}
    chosen = min(positive, key=lambda row: row["key"])
    return {
        **base, "state": "SPLIT_DIAGNOSTIC_ONLY",
        "improvement_bits": chosen["improvement_bits"],
        "expanded_child_contribution_log2": chosen["expanded_child_contribution_log2"],
        "collapsed_child_contribution_log2": chosen["collapsed_child_contribution_log2"],
        "split": {"k": chosen["k"], "t": str(chosen["t"]), "coefficients": [1 if index <= chosen["k"] else 0 for index in range(9)], "threshold": str(chosen["t"] + chosen["k"] + 1)},
        "diagnostic": {"selected_transition": chosen["transition"]},
        "left": chosen["left"], "right": chosen["right"],
    }


def install_replay(state: dict[str, Any], replay: dict[str, Any]) -> None:
    rows = {row["node_id"]: row for row in replay["leaves"]}
    active = engine.active_leaf_ids(state)
    if set(rows) != set(active):
        raise ValueError("factorized replay leaf set changed")
    for node_id in active:
        node, row = state["nodes"][node_id], rows[node_id]
        if row["exact_count"] != node["exact_count"] or row["exact_vertex_sha256"] != node["exact_vertex_sha256"]:
            raise ValueError("factorized replay geometry binding changed")
        node["candidate_upper_log2"] = row["candidate_upper_log2"]
        node["contribution_log2"] = row["contribution_log2"]
        node["selector"] = {"kind": "independent-components", "inner": row["inner_selector"], "outer": row["outer_selector"]}
        node["vertex_dual"] = row["vertex_dual"]


def load_inputs(args: argparse.Namespace):
    digests = {"geometry_sha256": coordinate.sha256_path(args.geometry), "factorized_checkpoint_sha256": coordinate.sha256_path(args.factorized_checkpoint)}
    if digests != {"geometry_sha256": GEOMETRY_SHA256, "factorized_checkpoint_sha256": FACTORIZED_SHA256}:
        raise ValueError("frozen geometry or factorized checkpoint changed")
    geometry_state = json.loads(args.geometry.read_bytes())
    factorized_state = json.loads(args.factorized_checkpoint.read_bytes())
    origins = {row["origin_catalogue_sha256"] for row in factorized_state["candidates"]}
    if origins != {ROUND1_SHA256, ROUND2_SHA256} or len(factorized_state["candidates"]) != 24:
        raise ValueError("factorized candidate catalogue closure changed")
    base, bindings = component.load_component_bank(args.manifest, args.atlas, args.old_supplementary, args.new_supplementary)
    bank = build_bank(base, factorized_state["candidates"])
    replay = factorized.replay_all(geometry_state, base, factorized_state["candidates"])
    if coordinate.canonical_bytes(replay) != coordinate.canonical_bytes(factorized_state["current_replay"]):
        raise ValueError("full factorized bank replay differs from frozen checkpoint")
    return geometry_state, factorized_state, bank, replay, {**digests, **bindings, "round1_catalogue_sha256": ROUND1_SHA256, "round2_catalogue_sha256": ROUND2_SHA256}


def build(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    geometry_state, factorized_state, bank, replay, bindings = load_inputs(args)
    state = copy.deepcopy(geometry_state)
    install_replay(state, replay)
    before = engine.global_diagnostic(state)
    ranked = factorized.distinct_h2_ranking(replay)
    selected = [row["node_id"] for row in ranked[:TOP_LEAVES]]
    results = {node_id: refine_leaf(state["nodes"][node_id], bank) for node_id in selected}
    outcome = persistence.apply_one_wave(state, selected, results)
    if not outcome["committed"]:
        raise RuntimeError(f"geometry wave did not commit: {outcome['stop_reason']}")
    after = engine.global_diagnostic(state)
    aggregate_improvement = before["aggregate_log2_union_diagnostic"] - after["aggregate_log2_union_diagnostic"]
    no_cap = (
        outcome["actual_child_vertex_incidences"] <= persistence.geometry.VERTEX_CAP_WAVE
        and all(results[node_id][side]["exact_vertex_count"] <= persistence.geometry.VERTEX_CAP_PER_CHILD for node_id in outcome["accepted_leaf_ids"] for side in ("left", "right"))
    )
    no_regression = outcome["expanded_selected_frontier_strictly_decreases"] and aggregate_improvement >= -1e-8
    outcome["gate"] = {
        "minimum_aggregate_improvement_bits": MINIMUM_AGGREGATE_IMPROVEMENT,
        "aggregate_improvement_bits": aggregate_improvement,
        "aggregate_improvement_passes": aggregate_improvement >= MINIMUM_AGGREGATE_IMPROVEMENT,
        "no_vertex_cap_violation": no_cap, "no_expanded_regression": no_regression,
        "overall_passes": aggregate_improvement >= MINIMUM_AGGREGATE_IMPROVEMENT and no_cap and no_regression,
        "outward_replay_performed": False,
    }
    outcome["splits"] = [{
        "parent_id": node_id, "h2_cell": int(state["nodes"][node_id]["h2_cell"]),
        "split": results[node_id]["split"], "parent_improvement_bits": results[node_id]["improvement_bits"],
        "expanded_child_contribution_log2": results[node_id]["expanded_child_contribution_log2"],
        "collapsed_child_contribution_log2": results[node_id]["collapsed_child_contribution_log2"],
        "child_vertex_counts": [results[node_id]["left"]["exact_vertex_count"], results[node_id]["right"]["exact_vertex_count"]],
    } for node_id in outcome["accepted_leaf_ids"]]
    outcome["runtime_seconds"] = time.perf_counter() - started
    state.update(
        schema=SCHEMA, status="STOPPED_DIAGNOSTIC_ONE_FACTORIZED_GEOMETRY_WAVE",
        policy="top8-distinct-h2-factorized-transition-prefix-split-v1",
        source_bindings={**bindings, "inner_component_count": len(bank.inner_references), "outer_component_count": len(bank.outer_references)},
        factorized_geometry_wave=outcome, global_diagnostic=after,
        stop_reason="MAX_FACTORIZED_GEOMETRY_WAVES",
        scope_limit="Exact geometry/counting; binary64 rational mixtures require outward replay.",
    )
    return state


def self_test(args: argparse.Namespace) -> None:
    _geometry, _state, bank, replay, _bindings = load_inputs(args)
    if len(bank.inner_references) != 566 or len(bank.outer_references) != 574:
        raise SystemExit("full factorized bank cardinality failed")
    cells = [row["h2_cell"] for row in factorized.distinct_h2_ranking(replay)[:8]]
    if len(cells) != 8 or len(set(cells)) != 8:
        raise SystemExit("distinct-h2 selection failed")
    print("full_bank_inner_566_outer_574=PASS")
    print("frozen_189_leaf_replay_equivalence=PASS")
    print("top8_distinct_h2=" + ",".join(f"{cell:03d}" for cell in cells))
    print("status=PASS_G8_FACTORIZED_GEOMETRY_PERSISTENCE_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=component.h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=component.h2.DEFAULT_ATLAS)
    parser.add_argument("--old-supplementary", type=Path, default=component.h2.DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--new-supplementary", type=Path, default=component.NEW_ATLAS)
    parser.add_argument("--geometry", type=Path, default=GEOMETRY)
    parser.add_argument("--factorized-checkpoint", type=Path, default=FACTORIZED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(args)
        return
    report = build(args)
    digest = engine.atomic_write(args.output, report)
    wave = report["factorized_geometry_wave"]
    print(f"accepted_splits={wave['accepted_split_count']}")
    print(f"aggregate_improvement_bits={wave['aggregate_improvement_bits']:.12f}")
    print(f"worst_improvement_bits={wave['global_worst_improvement_bits']:.12f}")
    print(f"gate_passes={wave['gate']['overall_passes']}")
    print(f"runtime_seconds={wave['runtime_seconds']:.3f}")
    print(f"output={args.output}")
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
