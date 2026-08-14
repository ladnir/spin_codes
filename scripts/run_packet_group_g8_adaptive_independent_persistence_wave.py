#!/usr/bin/env python3
"""Run one exact cumulative persistence wave with independent components."""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import time
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

import evaluate_packet_group_g8_adaptive_cumulative_wave as geometry
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import replay_packet_group_g8_adaptive_checkpoint_with_atlas as stored
import replay_packet_group_g8_independent_component_mixtures as component
import run_packet_group_g8_adaptive_cumulative_engine as engine


ROOT = Path(__file__).resolve().parents[1]
FIXED_REPLAY = ROOT / "out" / "g8_full_support_adaptive_independent_component_replay.json"
FIXED_REPLAY_SHA256 = "0c734b79e1085e47c0fc2ea90d524d24cb6a67e2b2980ed143c962b37123cf6c"
DEFAULT_OUTPUT = ROOT / "out" / "g8_full_support_adaptive_independent_persistence_wave.json"
SCHEMA = "permute-conv.packet-group-g8-independent-persistence-checkpoint.v1"
TOP_LEAVES = 8


def score_bounds(
    lower: tuple[int, ...],
    upper: tuple[int, ...],
    bank: component.ComponentBank,
    denominator: int,
    node_id: str,
    h2_cell: int,
    depth: int,
) -> dict[str, Any] | None:
    vertices = geometry.chain_vertices(lower, upper)
    if not vertices or len(vertices) > geometry.VERTEX_CAP_PER_CHILD:
        return None
    count = geometry.exact_chain_count(lower, upper)
    if count <= 0:
        return None
    profiles = [list(geometry.profile_from_chain(vertex)) for vertex in vertices]
    digest = coordinate.sha256_bytes(coordinate.canonical_bytes(profiles))
    temporary = {
        "node_id": node_id,
        "h2_cell": h2_cell,
        "depth": depth,
        "exact_count": str(count),
        "exact_vertex_count": len(profiles),
        "exact_vertices": profiles,
        "exact_vertex_sha256": digest,
        "candidate_upper_log2": 0.0,
        "contribution_log2": 0.0,
    }
    replayed = component.replay_leaf(temporary, bank, denominator)
    return {
        "lower": list(lower),
        "upper": list(upper),
        "exact_count": str(count),
        "count": {
            "kind": "exact",
            "method": "bounded-monotone-chain-prefix-dp-v1",
            "parameters": {
                "lower": list(lower),
                "upper": list(upper),
                "free_mass": str(262144 - 9),
            },
            "value": str(count),
        },
        "exact_vertex_count": len(profiles),
        "exact_vertices": profiles,
        "exact_vertex_sha256": digest,
        "candidate_upper_log2": replayed["candidate_upper_log2"],
        "contribution_log2": replayed["contribution_log2"],
        "worst_vertex": replayed["worst_vertex"],
        "selector": {
            "kind": "independent-components",
            "inner": replayed["inner_selector"],
            "outer": replayed["outer_selector"],
        },
        "component_count": replayed["inner_support_size"] + replayed["outer_support_size"],
        "inner_support_size": replayed["inner_support_size"],
        "outer_support_size": replayed["outer_support_size"],
        "candidate_status": "DIAGNOSTIC_BINARY64_REQUIRES_OUTWARD_REPLAY",
        "diagnostic": replayed["diagnostic"],
    }


def leading_pairs(
    vertices: tuple[tuple[int, ...], ...], bank: component.ComponentBank
) -> tuple[tuple[int, int], ...]:
    profiles = np.asarray(
        [geometry.profile_from_chain(vertex) for vertex in vertices], dtype=np.float64
    )
    inner = bank.inner_constants[None, :] - profiles @ bank.inner_charges.T
    outer = bank.outer_constants[None, :] - profiles @ bank.outer_charges.T
    return tuple(
        (int(np.argmin(inner[row])), int(np.argmin(outer[row])))
        for row in range(len(vertices))
    )


def transition_candidates(
    vertices: tuple[tuple[int, ...], ...], bank: component.ComponentBank
) -> dict[tuple[int, int], dict[str, Any]]:
    pairs = leading_pairs(vertices, bank)
    profiles = tuple(geometry.profile_from_chain(vertex) for vertex in vertices)
    proposals: dict[tuple[int, int], dict[str, Any]] = {}
    for k in range(8):
        groups: dict[tuple[int, ...], list[int]] = {}
        for index, vertex in enumerate(vertices):
            groups.setdefault(vertex[:k] + vertex[k + 1 :], []).append(index)
        for indices in groups.values():
            indices.sort(key=lambda index: (vertices[index][k], index))
            for left_index, right_index in zip(indices, indices[1:]):
                left_value = vertices[left_index][k]
                right_value = vertices[right_index][k]
                left_pair = pairs[left_index]
                right_pair = pairs[right_index]
                if left_value >= right_value or left_pair == right_pair:
                    continue
                left_constant = (
                    bank.inner_constants[left_pair[0]] + bank.outer_constants[left_pair[1]]
                )
                right_constant = (
                    bank.inner_constants[right_pair[0]] + bank.outer_constants[right_pair[1]]
                )
                left_charge = bank.inner_charges[left_pair[0]] + bank.outer_charges[left_pair[1]]
                right_charge = bank.inner_charges[right_pair[0]] + bank.outer_charges[right_pair[1]]
                delta_constant = float(left_constant - right_constant)
                delta_charge = left_charge - right_charge
                d_left = delta_constant - float(delta_charge @ np.asarray(profiles[left_index]))
                d_right = delta_constant - float(delta_charge @ np.asarray(profiles[right_index]))
                thresholds = set()
                crossing = None
                if (
                    math.isfinite(d_left)
                    and math.isfinite(d_right)
                    and d_left != d_right
                    and ((d_left <= 0.0 <= d_right) or (d_right <= 0.0 <= d_left))
                ):
                    crossing = left_value + (-d_left / (d_right - d_left)) * (
                        right_value - left_value
                    )
                    if math.isfinite(crossing):
                        thresholds.update((math.floor(crossing), math.ceil(crossing) - 1))
                if not thresholds:
                    thresholds.add((left_value + right_value) // 2)
                for threshold in thresholds:
                    if not left_value <= threshold < right_value:
                        continue
                    record = {
                        "coordinate": k,
                        "threshold": str(threshold),
                        "left_leading_inner": bank.references[left_pair[0]],
                        "left_leading_outer": bank.references[left_pair[1]],
                        "right_leading_inner": bank.references[right_pair[0]],
                        "right_leading_outer": bank.references[right_pair[1]],
                        "crossing_coordinate_binary64": crossing,
                    }
                    key = (k, threshold)
                    if key not in proposals or coordinate.canonical_bytes(record) < coordinate.canonical_bytes(proposals[key]):
                        proposals[key] = record
    return proposals


def replay_parent(
    node: dict[str, Any], bank: component.ComponentBank, denominator: int
) -> tuple[dict[str, Any], tuple[tuple[int, ...], ...]]:
    lower = tuple(map(int, node["lower"]))
    upper = tuple(map(int, node["upper"]))
    exact = score_bounds(
        lower, upper, bank, denominator, node["node_id"], int(node["h2_cell"]), int(node["depth"])
    )
    if exact is None:
        raise RuntimeError("selected parent failed exact component replay")
    if (
        exact["exact_count"] != node["exact_count"]
        or exact["exact_vertex_sha256"] != node["exact_vertex_sha256"]
        or exact["exact_vertex_count"] != int(node["exact_vertex_count"])
    ):
        raise RuntimeError("selected parent geometry changed during component replay")
    return exact, geometry.chain_vertices(lower, upper)


def refine_leaf(
    node: dict[str, Any], bank: component.ComponentBank, denominator: int
) -> dict[str, Any]:
    parent, vertices = replay_parent(node, bank, denominator)
    lower = tuple(map(int, node["lower"]))
    upper = tuple(map(int, node["upper"]))
    proposals = transition_candidates(vertices, bank)
    candidates = []
    cap_rejections = 0
    for (k, threshold), transition in sorted(proposals.items()):
        left_bounds = geometry.child_bounds(lower, upper, k, threshold, True)
        right_bounds = geometry.child_bounds(lower, upper, k, threshold, False)
        if left_bounds is None or right_bounds is None:
            continue
        if (
            len(geometry.chain_vertices(*left_bounds)) > geometry.VERTEX_CAP_PER_CHILD
            or len(geometry.chain_vertices(*right_bounds)) > geometry.VERTEX_CAP_PER_CHILD
        ):
            cap_rejections += 1
            continue
        left = score_bounds(
            *left_bounds, bank, denominator, f"{node['node_id']}/L", int(node["h2_cell"]), int(node["depth"]) + 1
        )
        right = score_bounds(
            *right_bounds, bank, denominator, f"{node['node_id']}/R", int(node["h2_cell"]), int(node["depth"]) + 1
        )
        if left is None or right is None:
            continue
        if int(left["exact_count"]) + int(right["exact_count"]) != int(parent["exact_count"]):
            raise RuntimeError("component child counts do not sum to parent")
        expanded = engine.log2_sum([left["contribution_log2"], right["contribution_log2"]])
        maximum_bound = max(left["candidate_upper_log2"], right["candidate_upper_log2"])
        collapsed = math.log2(int(parent["exact_count"])) + maximum_bound
        candidates.append(
            {
                "key": (
                    expanded,
                    maximum_bound,
                    -min(int(left["exact_count"]), int(right["exact_count"])),
                    left["exact_vertex_count"] + right["exact_vertex_count"],
                    k,
                    threshold,
                ),
                "improvement_bits": parent["contribution_log2"] - expanded,
                "expanded_child_contribution_log2": expanded,
                "collapsed_child_contribution_log2": collapsed,
                "k": k,
                "t": threshold,
                "transition": transition,
                "left": left,
                "right": right,
            }
        )
    positive = [row for row in candidates if row["improvement_bits"] > 0.0]
    base = {
        "parent_replay": parent,
        "transition_candidate_count": len(proposals),
        "evaluated_candidate_count": len(candidates),
        "child_cap_rejection_count": cap_rejections,
    }
    if not positive:
        return {**base, "state": "NO_VALID_SPLIT"}
    chosen = min(positive, key=lambda row: row["key"])
    return {
        **base,
        "state": "SPLIT_DIAGNOSTIC_ONLY",
        "improvement_bits": chosen["improvement_bits"],
        "expanded_child_contribution_log2": chosen["expanded_child_contribution_log2"],
        "collapsed_child_contribution_log2": chosen["collapsed_child_contribution_log2"],
        "split": {
            "k": chosen["k"],
            "t": str(chosen["t"]),
            "coefficients": [1 if index <= chosen["k"] else 0 for index in range(9)],
            "threshold": str(chosen["t"] + chosen["k"] + 1),
        },
        "diagnostic": {"selected_transition": chosen["transition"]},
        "left": chosen["left"],
        "right": chosen["right"],
    }


def install_fixed_replay(state: dict[str, Any], artifact: dict[str, Any]) -> None:
    rows = {row["node_id"]: row for row in artifact["leaves"]}
    active = engine.active_leaf_ids(state)
    if set(rows) != set(active):
        raise ValueError("fixed component replay leaf set differs from checkpoint")
    for node_id in active:
        node = state["nodes"][node_id]
        row = rows[node_id]
        if row["exact_count"] != node["exact_count"] or row["exact_vertex_sha256"] != node["exact_vertex_sha256"]:
            raise ValueError("fixed component replay geometry binding changed")
        node["candidate_upper_log2"] = row["candidate_upper_log2"]
        node["contribution_log2"] = row["contribution_log2"]
        node["selector"] = {
            "kind": "independent-components",
            "inner": row["inner_selector"],
            "outer": row["outer_selector"],
        }
        node["component_count"] = row["inner_support_size"] + row["outer_support_size"]


def child_node(parent: dict[str, Any], side: str, diagnostic: dict[str, Any]) -> dict[str, Any]:
    node = engine.child_node(parent, side, diagnostic)
    node["inner_support_size"] = diagnostic["inner_support_size"]
    node["outer_support_size"] = diagnostic["outer_support_size"]
    return node


def apply_one_wave(
    state: dict[str, Any], selected: list[str], results: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    before = engine.global_diagnostic(state)
    accepted = [node_id for node_id in selected if results[node_id]["state"] == "SPLIT_DIAGNOSTIC_ONLY"]
    if not accepted:
        return {"committed": False, "stop_reason": "NO_VALID_SPLIT"}
    incidences = sum(results[node_id][side]["exact_vertex_count"] for node_id in accepted for side in ("left", "right"))
    if incidences > geometry.VERTEX_CAP_WAVE:
        return {"committed": False, "stop_reason": "WAVE_VERTEX_INCIDENCE_CAP"}
    selected_before = engine.log2_sum([results[node_id]["parent_replay"]["contribution_log2"] for node_id in selected])
    for node_id in selected:
        replayed = results[node_id]["parent_replay"]
        state["nodes"][node_id]["candidate_upper_log2"] = replayed["candidate_upper_log2"]
        state["nodes"][node_id]["contribution_log2"] = replayed["contribution_log2"]
        state["nodes"][node_id]["selector"] = replayed["selector"]
    for node_id in accepted:
        parent = state["nodes"][node_id]
        result = results[node_id]
        left = child_node(parent, "L", result["left"])
        right = child_node(parent, "R", result["right"])
        state["nodes"][left["node_id"]] = left
        state["nodes"][right["node_id"]] = right
        parent.update(
            state="SPLIT",
            split=result["split"],
            left_id=left["node_id"],
            right_id=right["node_id"],
            split_diagnostic={
                "improvement_bits": result["improvement_bits"],
                "expanded_child_contribution_log2": result["expanded_child_contribution_log2"],
                "collapsed_child_contribution_log2": result["collapsed_child_contribution_log2"],
                **result["diagnostic"],
            },
        )
    after = engine.global_diagnostic(state)
    selected_after = engine.log2_sum(
        [
            value
            for node_id in selected
            for value in (
                [results[node_id]["left"]["contribution_log2"], results[node_id]["right"]["contribution_log2"]]
                if node_id in accepted
                else [results[node_id]["parent_replay"]["contribution_log2"]]
            )
        ]
    )
    return {
        "committed": True,
        "before": before,
        "after": after,
        "selected_leaf_ids": selected,
        "accepted_leaf_ids": accepted,
        "accepted_split_count": len(accepted),
        "global_worst_improvement_bits": before["worst_contribution_log2"] - after["worst_contribution_log2"],
        "aggregate_improvement_bits": before["aggregate_log2_union_diagnostic"] - after["aggregate_log2_union_diagnostic"],
        "expanded_selected_frontier_before_log2": selected_before,
        "expanded_selected_frontier_after_log2": selected_after,
        "expanded_selected_frontier_strictly_decreases": selected_after < selected_before,
        "actual_child_vertex_incidences": incidences,
        "stop_reason": "MAX_COMPONENT_PERSISTENCE_WAVES",
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    checkpoint, checkpoint_digest = stored.load_checkpoint(args.checkpoint)
    if coordinate.sha256_path(args.fixed_replay) != FIXED_REPLAY_SHA256:
        raise ValueError("fixed independent replay digest changed")
    fixed = json.loads(args.fixed_replay.read_bytes())
    bank, bindings = component.load_component_bank(
        args.manifest, args.atlas, args.old_supplementary, args.new_supplementary
    )
    state = copy.deepcopy(checkpoint)
    install_fixed_replay(state, fixed)
    state["global_diagnostic"] = engine.global_diagnostic(state)
    selected = engine.top_leaf_ids(state, TOP_LEAVES)
    results = {node_id: refine_leaf(state["nodes"][node_id], bank, args.mixture_denominator) for node_id in selected}
    outcome = apply_one_wave(state, selected, results)
    outcome["runtime_seconds"] = time.perf_counter() - started
    state.update(
        schema=SCHEMA,
        status="STOPPED_DIAGNOSTIC_MAX_COMPONENT_PERSISTENCE_WAVES",
        policy="top8-independent-component-transition-prefix-split-v1",
        source_bindings={
            **bindings,
            "checkpoint_sha256": checkpoint_digest,
            "fixed_component_replay_sha256": FIXED_REPLAY_SHA256,
            "witness_sources": list(bank.sources),
        },
        component_persistence_wave=outcome,
        global_diagnostic=engine.global_diagnostic(state),
        stop_reason=outcome["stop_reason"],
        scope_limit="Exact geometry/counting; binary64 independent mixtures require outward replay.",
    )
    return state


def run_self_test() -> None:
    inner = np.asarray([[0.0, 2.0], [2.0, 0.0]])
    outer = np.asarray([[3.0, 0.0], [0.0, 3.0]])
    normalization = np.zeros(2)
    first = component.independent_minimax(inner, outer, normalization, 1 << 20)
    second = component.independent_minimax(inner, outer, normalization, 1 << 20)
    if first != second or abs(first["binary64_lp_score"] - 2.5) > 1e-12:
        raise SystemExit("component persistence scorer differs from fixed replay")
    print(f"synthetic_fixed_replay_rational_score={first['score']:.12f}")
    print(f"synthetic_persistence_rational_score={second['score']:.12f}")
    print("independent_component_scorer_equivalence=PASS")
    print("status=PASS_G8_INDEPENDENT_PERSISTENCE_WAVE_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=component.h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=component.h2.DEFAULT_ATLAS)
    parser.add_argument("--old-supplementary", type=Path, default=component.h2.DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--new-supplementary", type=Path, default=component.NEW_ATLAS)
    parser.add_argument("--checkpoint", type=Path, default=stored.DEFAULT_CHECKPOINT)
    parser.add_argument("--fixed-replay", type=Path, default=FIXED_REPLAY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mixture-denominator", type=int, default=1 << 30)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    report = build(args)
    digest = engine.atomic_write(args.output, report)
    row = report["component_persistence_wave"]
    print(f"accepted_splits={row['accepted_split_count']}")
    print(f"aggregate_improvement_bits={row['aggregate_improvement_bits']:.12f}")
    print(f"worst_improvement_bits={row['global_worst_improvement_bits']:.12f}")
    print(f"active_leaves={report['global_diagnostic']['active_leaf_count']}")
    print(f"runtime_seconds={row['runtime_seconds']:.3f}")
    print(f"stop_reason={report['stop_reason']}")
    print(f"output={args.output}")
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
