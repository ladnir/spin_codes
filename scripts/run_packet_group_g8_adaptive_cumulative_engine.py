#!/usr/bin/env python3
"""Resume bounded adaptive cumulative-prefix discovery for full-support g=8."""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import tempfile
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

import evaluate_packet_group_g8_adaptive_cumulative_wave as wave
import evaluate_packet_group_g8_cumulative_boxes as h2
import plan_packet_group_g8_lowdim_root_cover as minimax
import produce_packet_group_g8_highdim_adaptive_shards as coordinate


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CHECKPOINT = ROOT / "out" / "g8_full_support_adaptive_cumulative_checkpoint.json"
SCHEMA = "permute-conv.packet-group-g8-adaptive-cumulative-checkpoint.v1"
POLICY = "top8-witness-transition-prefix-split-v1"
TOP_LEAVES = 8
MINIMUM_WORST_IMPROVEMENT = 16.0
LOW_IMPROVEMENT_LIMIT = 2
TWO_WAVE_GATE_LENGTH = 2
MINIMUM_ACCEPTED_SPLITS = 6
MAXIMUM_PROJECTED_LEVELS = 8
MAXIMUM_PROJECTED_LEAVES = 2048


def log2_sum(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log2(
        math.fsum(math.exp2(value - maximum) for value in values)
    )


def active_leaf_ids(state: dict[str, Any]) -> list[str]:
    return sorted(
        node_id
        for node_id, node in state["nodes"].items()
        if node["state"] == "ACTIVE_LEAF"
    )


def top_leaf_ids(state: dict[str, Any], limit: int = TOP_LEAVES) -> list[str]:
    leaves = active_leaf_ids(state)
    leaves.sort(
        key=lambda node_id: (
            -float(state["nodes"][node_id]["contribution_log2"]),
            node_id,
        )
    )
    return leaves[:limit]


def global_diagnostic(state: dict[str, Any]) -> dict[str, Any]:
    leaves = active_leaf_ids(state)
    if not leaves:
        raise RuntimeError("adaptive checkpoint has no active leaves")
    contributions = [float(state["nodes"][node]["contribution_log2"]) for node in leaves]
    worst = max(
        leaves,
        key=lambda node: (float(state["nodes"][node]["contribution_log2"]), node),
    )
    return {
        "active_leaf_count": len(leaves),
        "aggregate_log2_union_diagnostic": log2_sum(contributions),
        "worst_leaf_id": worst,
        "worst_contribution_log2": float(state["nodes"][worst]["contribution_log2"]),
        "worst_candidate_upper_log2": float(
            state["nodes"][worst]["candidate_upper_log2"]
        ),
    }


def atomic_write(path: Path, value: dict[str, Any]) -> str:
    payload = coordinate.canonical_bytes(value) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return coordinate.sha256_bytes(payload)


def source_bindings(args: argparse.Namespace) -> tuple[coordinate.AtlasBank, dict[str, Any]]:
    bank, sources, digests = h2.build_bank(
        args.manifest, args.atlas, args.supplementary
    )
    return bank, {
        **digests,
        "h2_artifact_sha256": wave.PARENT_SHA256,
        "h2_artifact_schema": h2.SCHEMA,
        "witness_sources": sources,
    }


def initial_node(cell: dict[str, Any]) -> dict[str, Any]:
    lower, upper = wave.parent_bounds(cell)
    profiles = [list(wave.parse_profile(row)) for row in cell["vertices"]]
    node_id = f"h2:{int(cell['cell']):03d}"
    return {
        "node_id": node_id,
        "parent_id": None,
        "branch": None,
        "depth": 0,
        "state": "ACTIVE_LEAF",
        "h2_cell": int(cell["cell"]),
        "h2_bin_indices": cell["bin_indices"],
        "lower": list(lower),
        "upper": list(upper),
        "exact_count": cell["exact_count"],
        "exact_vertex_count": len(profiles),
        "exact_vertices": profiles,
        "exact_vertex_sha256": coordinate.sha256_bytes(
            coordinate.canonical_bytes(profiles)
        ),
        "candidate_upper_log2": float(cell["mixture"]["candidate_upper_log2"]),
        "contribution_log2": float(cell["mixture"]["contribution_log2"]),
        "selector": cell["mixture"]["selector"],
        "component_count": int(cell["mixture"]["component_count"]),
        "candidate_status": cell["mixture"]["candidate_status"],
        "origin": "frozen-h2-diagnostic",
    }


def initialize_state(
    args: argparse.Namespace, bindings: dict[str, Any]
) -> dict[str, Any]:
    if coordinate.sha256_path(args.h2_artifact) != wave.PARENT_SHA256:
        raise ValueError("h2 artifact digest changed")
    parent = json.loads(args.h2_artifact.read_bytes())
    if parent.get("schema") != h2.SCHEMA or len(parent.get("cells", [])) != 165:
        raise ValueError("unexpected h2 artifact")
    discovery = parent.get("mixture_discovery", {})
    if (
        int(discovery.get("rational_denominator", -1)) != args.mixture_denominator
        or int(discovery.get("seed_per_vertex", -1)) != 1
        or float(discovery.get("reduced_cost_tolerance", -1.0)) != 1e-9
        or float(discovery.get("weight_tolerance", -1.0)) != 1e-12
    ):
        raise ValueError("engine settings differ from the h2 mixture settings")
    nodes = {}
    for cell in parent["cells"]:
        node = initial_node(cell)
        nodes[node["node_id"]] = node
    state = {
        "schema": SCHEMA,
        "status": "RUNNABLE_DIAGNOSTIC_REQUIRES_OUTWARD_REPLAY",
        "policy": POLICY,
        "source_bindings": bindings,
        "configuration": {
            "top_leaves_per_wave": TOP_LEAVES,
            "minimum_global_worst_improvement_bits": MINIMUM_WORST_IMPROVEMENT,
            "low_improvement_consecutive_wave_limit": LOW_IMPROVEMENT_LIMIT,
            "vertex_cap_per_child": wave.VERTEX_CAP_PER_CHILD,
            "vertex_cap_per_wave": wave.VERTEX_CAP_WAVE,
            "mixture_denominator": args.mixture_denominator,
            "two_wave_gate_length": TWO_WAVE_GATE_LENGTH,
            "minimum_accepted_splits_per_wave": MINIMUM_ACCEPTED_SPLITS,
            "maximum_projected_levels": MAXIMUM_PROJECTED_LEVELS,
            "maximum_projected_leaves": MAXIMUM_PROJECTED_LEAVES,
            "parent_replay": "exact-count-hull-and-identical-minimax-v1",
            "expanded_objective": "log2(2^(log2(nL)+UL)+2^(log2(nR)+UR))",
            "collapsed_objective": "log2(nParent)+max(UL,UR)",
        },
        "completed_waves": 0,
        "low_improvement_streak": 0,
        "stop_reason": None,
        "nodes": nodes,
        "waves": [],
        "global_diagnostic": None,
        "scope_limit": (
            "Exact integer ownership, counts, and vertices only. Binary64 split "
            "selection and rational-mixture scoring require independent outward "
            "replay. This checkpoint is not a certificate or completion claim."
        ),
    }
    state["global_diagnostic"] = global_diagnostic(state)
    return state


def validate_resumed_state(
    state: dict[str, Any], args: argparse.Namespace, bindings: dict[str, Any]
) -> None:
    if state.get("schema") != SCHEMA or state.get("policy") != POLICY:
        raise ValueError("unexpected adaptive checkpoint schema or policy")
    if state.get("source_bindings") != bindings:
        raise ValueError("adaptive checkpoint source bindings changed")
    expected_configuration = {
        "top_leaves_per_wave": TOP_LEAVES,
        "minimum_global_worst_improvement_bits": MINIMUM_WORST_IMPROVEMENT,
        "low_improvement_consecutive_wave_limit": LOW_IMPROVEMENT_LIMIT,
        "vertex_cap_per_child": wave.VERTEX_CAP_PER_CHILD,
        "vertex_cap_per_wave": wave.VERTEX_CAP_WAVE,
        "mixture_denominator": args.mixture_denominator,
        "two_wave_gate_length": TWO_WAVE_GATE_LENGTH,
        "minimum_accepted_splits_per_wave": MINIMUM_ACCEPTED_SPLITS,
        "maximum_projected_levels": MAXIMUM_PROJECTED_LEVELS,
        "maximum_projected_leaves": MAXIMUM_PROJECTED_LEAVES,
        "parent_replay": "exact-count-hull-and-identical-minimax-v1",
        "expanded_objective": "log2(2^(log2(nL)+UL)+2^(log2(nR)+UR))",
        "collapsed_objective": "log2(nParent)+max(UL,UR)",
    }
    if state.get("configuration") != expected_configuration:
        raise ValueError("adaptive checkpoint configuration changed")
    if int(state.get("completed_waves", -1)) != len(state.get("waves", [])):
        raise ValueError("adaptive checkpoint wave history is malformed")
    if state.get("global_diagnostic") != global_diagnostic(state):
        raise ValueError("adaptive checkpoint global diagnostic is stale")


def load_or_initialize(
    args: argparse.Namespace, bindings: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    if args.resume_from is None:
        return initialize_state(args, bindings), None
    resume_digest = coordinate.sha256_path(args.resume_from)
    state = json.loads(args.resume_from.read_bytes())
    validate_resumed_state(state, args, bindings)
    if (
        state.get("stop_reason") == "MAX_WAVES"
        and int(state["completed_waves"]) < args.max_waves
    ):
        state["stop_reason"] = None
        state["status"] = "RUNNABLE_DIAGNOSTIC_REQUIRES_OUTWARD_REPLAY"
    return state, {
        "path": str(args.resume_from),
        "sha256": resume_digest,
        "completed_waves": int(state["completed_waves"]),
    }


def replay_parent(
    node: dict[str, Any], bank: coordinate.AtlasBank, denominator: int
) -> tuple[dict[str, Any], tuple[tuple[int, ...], ...], np.ndarray]:
    lower = tuple(map(int, node["lower"]))
    upper = tuple(map(int, node["upper"]))
    closed = wave.close_bounds(lower, upper)
    if closed != (lower, upper):
        raise RuntimeError("stored cumulative bounds are not monotone-closed")
    exact_count = wave.exact_chain_count(lower, upper)
    if exact_count != int(node["exact_count"]):
        raise RuntimeError("stored parent count failed exact replay")
    vertices = wave.chain_vertices(lower, upper)
    profiles = tuple(sorted(wave.profile_from_chain(vertex) for vertex in vertices))
    stored_profiles = {tuple(map(int, row)) for row in node["exact_vertices"]}
    if set(profiles) != stored_profiles:
        raise RuntimeError("stored parent hull failed exact replay")
    digest = coordinate.sha256_bytes(
        coordinate.canonical_bytes([list(profile) for profile in profiles])
    )
    if digest != node["exact_vertex_sha256"]:
        raise RuntimeError("stored parent vertex digest failed exact replay")
    values = wave.score_matrix(vertices, bank)
    mixed = minimax.optimize_minimax_mixture(
        values,
        minimax.UNIFORM_TARGET_LOG2,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-9,
        weight_tolerance=1e-12,
        denominator=denominator,
    )
    upper_score = float(mixed["score"])
    replay = {
        "exact_count": str(exact_count),
        "exact_vertex_count": len(vertices),
        "exact_vertex_sha256": digest,
        "candidate_upper_log2": upper_score,
        "contribution_log2": math.log2(exact_count) + upper_score,
        "selector": h2.selector(mixed, bank),
        "component_count": len(mixed["indices"]),
        "binary64_lp_upper_log2": float(mixed["binary64_lp_score"]),
        "rationalization_loss_bits": float(mixed["rationalization_loss_bits"]),
    }
    return replay, vertices, values


def refine_leaf(
    node: dict[str, Any], bank: coordinate.AtlasBank, denominator: int
) -> dict[str, Any]:
    lower = tuple(map(int, node["lower"]))
    upper = tuple(map(int, node["upper"]))
    parent_replay, vertices, values = replay_parent(node, bank, denominator)
    leaders = wave.leading_witnesses(values, bank)
    proposals = wave.transition_candidates(vertices, values, leaders, bank)
    parent_count = int(parent_replay["exact_count"])
    parent_contribution = float(parent_replay["contribution_log2"])
    candidates = []
    cap_rejections = 0
    for (k, threshold), transition in sorted(proposals.items()):
        left_bounds = wave.child_bounds(lower, upper, k, threshold, True)
        right_bounds = wave.child_bounds(lower, upper, k, threshold, False)
        if left_bounds is None or right_bounds is None:
            continue
        left_vertices = wave.chain_vertices(*left_bounds)
        right_vertices = wave.chain_vertices(*right_bounds)
        if (
            len(left_vertices) > wave.VERTEX_CAP_PER_CHILD
            or len(right_vertices) > wave.VERTEX_CAP_PER_CHILD
        ):
            cap_rejections += 1
            continue
        left = wave.child_diagnostic(left_bounds, bank, denominator)
        right = wave.child_diagnostic(right_bounds, bank, denominator)
        if left is None or right is None:
            continue
        if int(left["exact_count"]) + int(right["exact_count"]) != parent_count:
            raise RuntimeError("resumable child counts do not sum to the parent")
        expanded_contribution = log2_sum(
            [left["contribution_log2"], right["contribution_log2"]]
        )
        maximum_bound = max(left["candidate_upper_log2"], right["candidate_upper_log2"])
        collapsed_contribution = math.log2(parent_count) + maximum_bound
        minimum_count = min(int(left["exact_count"]), int(right["exact_count"]))
        total_vertices = left["exact_vertex_count"] + right["exact_vertex_count"]
        candidates.append(
            {
                "key": (
                    expanded_contribution,
                    maximum_bound,
                    -minimum_count,
                    total_vertices,
                    k,
                    threshold,
                ),
                "improvement_bits": parent_contribution - expanded_contribution,
                "expanded_child_contribution_log2": expanded_contribution,
                "collapsed_child_contribution_log2": collapsed_contribution,
                "k": k,
                "t": threshold,
                "transition": transition,
                "left": left,
                "right": right,
            }
        )
    positive = [row for row in candidates if row["improvement_bits"] > 0.0]
    if not positive:
        return {
            "state": "NO_VALID_SPLIT",
            "transition_candidate_count": len(proposals),
            "evaluated_candidate_count": len(candidates),
            "child_cap_rejection_count": cap_rejections,
            "parent_replay": parent_replay,
        }
    chosen = min(positive, key=lambda row: row["key"])
    coefficients = [1 if index <= chosen["k"] else 0 for index in range(9)]
    return {
        "state": "SPLIT_DIAGNOSTIC_ONLY",
        "transition_candidate_count": len(proposals),
        "evaluated_candidate_count": len(candidates),
        "child_cap_rejection_count": cap_rejections,
        "parent_replay": parent_replay,
        "improvement_bits": chosen["improvement_bits"],
        "expanded_child_contribution_log2": chosen["expanded_child_contribution_log2"],
        "collapsed_child_contribution_log2": chosen["collapsed_child_contribution_log2"],
        "split": {
            "k": chosen["k"],
            "t": str(chosen["t"]),
            "coefficients": coefficients,
            "threshold": str(chosen["t"] + chosen["k"] + 1),
        },
        "diagnostic": {"selected_transition": chosen["transition"]},
        "left": chosen["left"],
        "right": chosen["right"],
    }


def child_node(
    parent: dict[str, Any], side: str, diagnostic: dict[str, Any]
) -> dict[str, Any]:
    node_id = f"{parent['node_id']}/{side}"
    return {
        "node_id": node_id,
        "parent_id": parent["node_id"],
        "branch": side,
        "depth": int(parent["depth"]) + 1,
        "state": "ACTIVE_LEAF",
        "h2_cell": parent["h2_cell"],
        "h2_bin_indices": parent["h2_bin_indices"],
        "lower": diagnostic["lower"],
        "upper": diagnostic["upper"],
        "exact_count": diagnostic["exact_count"],
        "count": diagnostic["count"],
        "exact_vertex_count": diagnostic["exact_vertex_count"],
        "exact_vertices": diagnostic["exact_vertices"],
        "exact_vertex_sha256": diagnostic["exact_vertex_sha256"],
        "candidate_upper_log2": diagnostic["candidate_upper_log2"],
        "contribution_log2": diagnostic["contribution_log2"],
        "worst_vertex": diagnostic["worst_vertex"],
        "selector": diagnostic["selector"],
        "component_count": diagnostic["component_count"],
        "candidate_status": diagnostic["candidate_status"],
        "diagnostic": diagnostic["diagnostic"],
        "origin": "adaptive-prefix-split",
    }


def apply_wave(
    state: dict[str, Any], selected: list[str], results: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    for node_id in selected:
        node = state["nodes"][node_id]
        replay = results[node_id]["parent_replay"]
        node["candidate_upper_log2"] = replay["candidate_upper_log2"]
        node["contribution_log2"] = replay["contribution_log2"]
        node["selector"] = replay["selector"]
        node["component_count"] = replay["component_count"]
        node["last_parent_replay"] = replay
    before = global_diagnostic(state)
    accepted = [
        node_id
        for node_id in selected
        if results[node_id]["state"] == "SPLIT_DIAGNOSTIC_ONLY"
    ]
    if not accepted:
        state["stop_reason"] = "NO_VALID_SPLIT"
        state["status"] = "STOPPED_DIAGNOSTIC_NO_VALID_SPLIT"
        return {"committed": False, "stop_reason": state["stop_reason"]}
    incidences = sum(
        int(results[node_id][side]["exact_vertex_count"])
        for node_id in accepted
        for side in ("left", "right")
    )
    if incidences > wave.VERTEX_CAP_WAVE:
        state["stop_reason"] = "WAVE_VERTEX_INCIDENCE_CAP"
        state["status"] = "STOPPED_DIAGNOSTIC_VERTEX_CAP"
        return {"committed": False, "stop_reason": state["stop_reason"]}
    selected_before = log2_sum(
        [float(results[node_id]["parent_replay"]["contribution_log2"]) for node_id in selected]
    )
    for node_id in accepted:
        parent = state["nodes"][node_id]
        result = results[node_id]
        children = {}
        for side, label in (("left", "L"), ("right", "R")):
            child = child_node(parent, label, result[side])
            if int(child["exact_vertex_count"]) > wave.VERTEX_CAP_PER_CHILD:
                state["stop_reason"] = "CHILD_VERTEX_CAP"
                state["status"] = "STOPPED_DIAGNOSTIC_VERTEX_CAP"
                return {"committed": False, "stop_reason": state["stop_reason"]}
            state["nodes"][child["node_id"]] = child
            children[side] = child["node_id"]
        parent["state"] = "SPLIT"
        parent["split"] = result["split"]
        parent["left_id"] = children["left"]
        parent["right_id"] = children["right"]
        parent["split_diagnostic"] = {
            "improvement_bits": result["improvement_bits"],
            "expanded_child_contribution_log2": result[
                "expanded_child_contribution_log2"
            ],
            "collapsed_child_contribution_log2": result[
                "collapsed_child_contribution_log2"
            ],
            "transition_candidate_count": result["transition_candidate_count"],
            "evaluated_candidate_count": result["evaluated_candidate_count"],
            "child_cap_rejection_count": result["child_cap_rejection_count"],
            **result["diagnostic"],
        }
    after = global_diagnostic(state)
    selected_after_terms = []
    for node_id in selected:
        result = results[node_id]
        if node_id in accepted:
            selected_after_terms.extend(
                [result["left"]["contribution_log2"], result["right"]["contribution_log2"]]
            )
        else:
            selected_after_terms.append(result["parent_replay"]["contribution_log2"])
    selected_after = log2_sum(selected_after_terms)
    selected_frontier_decreases = selected_after < selected_before
    improvement = before["worst_contribution_log2"] - after["worst_contribution_log2"]
    if improvement < MINIMUM_WORST_IMPROVEMENT:
        state["low_improvement_streak"] = int(state["low_improvement_streak"]) + 1
    else:
        state["low_improvement_streak"] = 0
    record = {
        "wave": int(state["completed_waves"]) + 1,
        "selected_leaf_ids": selected,
        "accepted_leaf_ids": accepted,
        "accepted_split_count": len(accepted),
        "before": before,
        "after": after,
        "global_worst_improvement_bits": improvement,
        "expanded_selected_frontier_before_log2": selected_before,
        "expanded_selected_frontier_after_log2": selected_after,
        "expanded_selected_frontier_strictly_decreases": selected_frontier_decreases,
        "actual_child_vertex_incidences": incidences,
        "splits": [
            {
                "parent_id": node_id,
                "split": results[node_id]["split"],
                "parent_improvement_bits": results[node_id]["improvement_bits"],
                "expanded_child_contribution_log2": results[node_id][
                    "expanded_child_contribution_log2"
                ],
                "collapsed_child_contribution_log2": results[node_id][
                    "collapsed_child_contribution_log2"
                ],
            }
            for node_id in accepted
        ],
        "unsplit_attempts": [
            {
                "parent_id": node_id,
                "reason": results[node_id]["state"],
                "transition_candidate_count": results[node_id]["transition_candidate_count"],
                "child_cap_rejection_count": results[node_id]["child_cap_rejection_count"],
            }
            for node_id in selected
            if node_id not in accepted
        ],
    }
    state["waves"].append(record)
    state["completed_waves"] = int(state["completed_waves"]) + 1
    state["global_diagnostic"] = after
    if int(state["low_improvement_streak"]) >= LOW_IMPROVEMENT_LIMIT:
        state["stop_reason"] = "LOW_GLOBAL_WORST_IMPROVEMENT_TWO_CONSECUTIVE_WAVES"
        state["status"] = "STOPPED_DIAGNOSTIC_LOW_IMPROVEMENT"
    return {"committed": True, "record": record, "stop_reason": state["stop_reason"]}


def lower_quartile(values: list[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    rank = max(1, math.ceil(len(ordered) / 4))
    return ordered[rank - 1]


def assess_two_wave_gate(state: dict[str, Any]) -> dict[str, Any] | None:
    if int(state["completed_waves"]) < TWO_WAVE_GATE_LENGTH:
        return None
    recent = state["waves"][-TWO_WAVE_GATE_LENGTH:]
    accepted_pass = all(
        int(record["accepted_split_count"]) >= MINIMUM_ACCEPTED_SPLITS
        for record in recent
    )
    expanded_pass = all(
        bool(record["expanded_selected_frontier_strictly_decreases"])
        for record in recent
    )
    improvements = [
        float(split["parent_improvement_bits"])
        for record in recent
        for split in record["splits"]
    ]
    quartile = lower_quartile(improvements)
    remaining_gap = max(
        0.0,
        float(state["global_diagnostic"]["aggregate_log2_union_diagnostic"]) + 40.0,
    )
    projected_levels = (
        math.inf if quartile <= 0.0 else math.ceil(remaining_gap / quartile)
    )
    finite_projection = min(
        int(projected_levels) if math.isfinite(projected_levels) else MAXIMUM_PROJECTED_LEVELS + 1,
        MAXIMUM_PROJECTED_LEVELS + 1,
    )
    projected_leaves = (
        int(state["global_diagnostic"]["active_leaf_count"])
        + TOP_LEAVES * ((1 << finite_projection) - 1)
    )
    levels_pass = projected_levels <= MAXIMUM_PROJECTED_LEVELS
    leaves_pass = projected_leaves <= MAXIMUM_PROJECTED_LEAVES
    gate = {
        "waves": [int(record["wave"]) for record in recent],
        "minimum_accepted_splits_per_wave": MINIMUM_ACCEPTED_SPLITS,
        "accepted_splits_per_wave": [
            int(record["accepted_split_count"]) for record in recent
        ],
        "accepted_split_gate_passes": accepted_pass,
        "expanded_selected_frontier_strict_decrease_passes": expanded_pass,
        "q25_parent_improvement_bits": quartile,
        "q25_method": "nearest-rank-ceil-one-quarter-v1",
        "remaining_aggregate_gap_to_2^-40_bits": remaining_gap,
        "projected_further_levels": (
            projected_levels if math.isfinite(projected_levels) else None
        ),
        "maximum_projected_further_levels": MAXIMUM_PROJECTED_LEVELS,
        "projected_levels_gate_passes": levels_pass,
        "projected_leaf_count": projected_leaves,
        "maximum_projected_leaf_count": MAXIMUM_PROJECTED_LEAVES,
        "projected_leaf_gate_passes": leaves_pass,
        "passes": accepted_pass and expanded_pass and levels_pass and leaves_pass,
    }
    state["two_wave_gate"] = gate
    return gate


def run_engine(
    state: dict[str, Any], bank: coordinate.AtlasBank, args: argparse.Namespace
) -> str:
    while int(state["completed_waves"]) < args.max_waves and state["stop_reason"] is None:
        selected = top_leaf_ids(state)
        if not selected:
            state["stop_reason"] = "NO_ACTIVE_LEAF"
            state["status"] = "STOPPED_DIAGNOSTIC_NO_VALID_SPLIT"
            break
        started = time.perf_counter()
        results = {
            node_id: refine_leaf(
                state["nodes"][node_id], bank, args.mixture_denominator
            )
            for node_id in selected
        }
        outcome = apply_wave(state, selected, results)
        if outcome.get("committed"):
            outcome["record"]["runtime_seconds"] = time.perf_counter() - started
            gate = assess_two_wave_gate(state)
            if gate is not None and not gate["passes"]:
                state["stop_reason"] = "TWO_WAVE_GO_GATE_FAILED"
                state["status"] = "STOPPED_DIAGNOSTIC_TWO_WAVE_GATE"
        atomic_write(args.checkpoint, state)
        if not outcome.get("committed") or state["stop_reason"] is not None:
            break
    if state["stop_reason"] is None and int(state["completed_waves"]) >= args.max_waves:
        state["stop_reason"] = "MAX_WAVES"
        state["status"] = "STOPPED_DIAGNOSTIC_MAX_WAVES"
    return atomic_write(args.checkpoint, state)


def synthetic_state() -> dict[str, Any]:
    nodes = {
        f"n{index}": {
            "node_id": f"n{index}",
            "state": "ACTIVE_LEAF",
            "contribution_log2": float(100 - index),
            "candidate_upper_log2": float(90 - index),
        }
        for index in range(12)
    }
    state = {
        "nodes": nodes,
        "completed_waves": 0,
        "low_improvement_streak": 0,
        "stop_reason": None,
        "waves": [],
    }
    state["global_diagnostic"] = global_diagnostic(state)
    return state


def synthetic_step(state: dict[str, Any]) -> None:
    selected = top_leaf_ids(state)
    before = global_diagnostic(state)
    for node_id in selected:
        parent = state["nodes"][node_id]
        parent["state"] = "SPLIT"
        for side, loss in (("L", 20.0), ("R", 21.0)):
            child_id = f"{node_id}/{side}"
            state["nodes"][child_id] = {
                "node_id": child_id,
                "state": "ACTIVE_LEAF",
                "contribution_log2": parent["contribution_log2"] - loss,
                "candidate_upper_log2": parent["candidate_upper_log2"] - loss,
            }
    after = global_diagnostic(state)
    state["completed_waves"] += 1
    state["waves"].append(
        {
            "wave": state["completed_waves"],
            "selected_leaf_ids": selected,
            "accepted_leaf_ids": selected,
            "accepted_split_count": len(selected),
            "before": before,
            "after": after,
            "expanded_selected_frontier_before_log2": log2_sum(
                [before["worst_contribution_log2"] - index for index in range(len(selected))]
            ),
            "expanded_selected_frontier_after_log2": log2_sum(
                [
                    state["nodes"][f"{node_id}/{side}"]["contribution_log2"]
                    for node_id in selected
                    for side in ("L", "R")
                ]
            ),
            "expanded_selected_frontier_strictly_decreases": True,
            "splits": [
                {"parent_id": node_id, "parent_improvement_bits": 20.0}
                for node_id in selected
            ],
        }
    )
    state["global_diagnostic"] = after


def run_self_test() -> None:
    uninterrupted = synthetic_state()
    for _ in range(4):
        synthetic_step(uninterrupted)
    resumed = synthetic_state()
    for _ in range(2):
        synthetic_step(resumed)
    with tempfile.TemporaryDirectory() as directory:
        checkpoint = Path(directory) / "checkpoint.json"
        atomic_write(checkpoint, resumed)
        resumed = json.loads(checkpoint.read_bytes())
        for _ in range(2):
            synthetic_step(resumed)
        atomic_write(checkpoint, resumed)
        if json.loads(checkpoint.read_bytes()) != uninterrupted:
            raise SystemExit("adaptive checkpoint resume-equivalence self-test failed")
    if top_leaf_ids(uninterrupted) != top_leaf_ids(resumed):
        raise SystemExit("adaptive checkpoint top-leaf self-test failed")
    gate_state = synthetic_state()
    synthetic_step(gate_state)
    synthetic_step(gate_state)
    gate = assess_two_wave_gate(gate_state)
    if gate is None or not gate["passes"]:
        raise SystemExit("adaptive checkpoint two-wave gate self-test failed")
    expanded = log2_sum([10.0, 10.0])
    if abs(expanded - 11.0) > 1e-12:
        raise SystemExit("adaptive expanded-objective self-test failed")
    print("uninterrupted_waves=4")
    print("resumed_waves=4")
    print("atomic_checkpoint_roundtrip=PASS")
    print("top8_selection_resume_equivalence=PASS")
    print("expanded_two_child_objective=PASS")
    print("two_wave_stop_go_gate=PASS")
    print("status=PASS_G8_ADAPTIVE_CUMULATIVE_ENGINE_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=h2.DEFAULT_ATLAS)
    parser.add_argument("--supplementary", type=Path, default=h2.DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--h2-artifact", type=Path, default=wave.DEFAULT_PARENT)
    parser.add_argument("--resume-from", type=Path)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--max-waves", type=int, required=False)
    parser.add_argument("--mixture-denominator", type=int, default=1 << 30)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.max_waves is None:
        parser.error("--max-waves is required outside --self-test")
    if args.max_waves < 0 or args.mixture_denominator < 1:
        parser.error("invalid engine limit")
    bank, bindings = source_bindings(args)
    state, resume_binding = load_or_initialize(args, bindings)
    if resume_binding is not None:
        state["last_resume_source"] = resume_binding
    digest = run_engine(state, bank, args)
    print(f"completed_waves={state['completed_waves']}")
    print(f"active_leaves={state['global_diagnostic']['active_leaf_count']}")
    print(f"worst_contribution_log2={state['global_diagnostic']['worst_contribution_log2']:.12f}")
    print(f"stop_reason={state['stop_reason']}")
    print(f"checkpoint={args.checkpoint}")
    print(f"sha256={digest}")
    print(f"status={state['status']}")


if __name__ == "__main__":
    main()
