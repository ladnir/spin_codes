#!/usr/bin/env python3
"""Run resumable cell-aware factorized column-generation rounds for g=8."""

from __future__ import annotations

import argparse
import copy
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import extract_packet_group_g8_component_pricing_target as pricing
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import replay_packet_group_g8_adaptive_checkpoint_with_atlas as stored
import replay_packet_group_g8_independent_component_mixtures as component
import run_packet_group_g8_adaptive_cumulative_engine as cumulative
import run_packet_group_g8_cell_aware_factorized_columns as one_shot


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "out" / "g8_full_support_adaptive_independent_persistence_wave.json"
GEOMETRY_SHA256 = one_shot.INPUT_SHA256
PRIOR_MAIN = ROOT / "out" / "g8_cell_aware_factorized_columns.json"
PRIOR_MAIN_SHA256 = "ccf82e0316db84fac37f7a4fdd78097aa60268eff4be494ed335eede698c5bb2"
PRIOR_CATALOGUE = ROOT / "out" / "g8_cell_aware_factorized_columns" / "candidate_catalogue.json"
PRIOR_CATALOGUE_SHA256 = "edc9c8b1714024d99af4c8d97352a8d20ec5f880202962667030d1ea2d9ae6f3"
DEFAULT_WORK_DIR = ROOT / "out" / "g8_iterative_factorized_columns"
DEFAULT_CHECKPOINT = DEFAULT_WORK_DIR / "checkpoint.json"
SCHEMA = "permute-conv.packet-group-g8-iterative-factorized-columns.v1"
ROUND_SCHEMA = "permute-conv.packet-group-g8-factorized-round.v1"
LANE_POLICIES = {
    "both": ("inner", "outer"),
    "inner": ("inner",),
    "outer": ("outer",),
}


def candidate_reference(row: dict[str, Any], catalogue_digest: str, ordinal: int) -> dict[str, Any]:
    return {
        "source_id": f"factorized-{catalogue_digest[:16]}-{row['lane']}",
        "source_sha256": catalogue_digest,
        "row": ordinal,
        "row_sha256": row["parameter_digest"],
    }


def initial_candidates(catalogue: dict[str, Any], digest: str) -> list[dict[str, Any]]:
    result = []
    for ordinal, row in enumerate(catalogue["rows"]):
        if not row.get("accepted"):
            continue
        copied = copy.deepcopy(row)
        copied["reference"] = candidate_reference(copied, digest, ordinal)
        copied["origin_catalogue_sha256"] = digest
        result.append(copied)
    return result


def pool(
    bank: component.ComponentBank, candidates: list[dict[str, Any]], lane: str
) -> tuple[np.ndarray, np.ndarray, tuple[dict[str, Any], ...]]:
    selected = [row for row in candidates if row["lane"] == lane and row.get("accepted")]
    base_constants = bank.inner_constants if lane == "inner" else bank.outer_constants
    base_charges = bank.inner_charges if lane == "inner" else bank.outer_charges
    if not selected:
        return base_constants, base_charges, bank.references
    constants = np.concatenate(
        (base_constants, np.asarray([row["parameters"]["constant_log2"] for row in selected]))
    )
    charges = np.vstack(
        (base_charges, np.asarray([row["parameters"]["charge_log2"] for row in selected]))
    )
    references = (*bank.references, *(row["reference"] for row in selected))
    return constants, charges, references


def selector(
    rows: tuple[tuple[int, Fraction], ...], references: tuple[dict[str, Any], ...]
) -> dict[str, Any]:
    return {
        "components": [
            {
                "weight": f"{weight.numerator}/{weight.denominator}",
                "component": references[index],
            }
            for index, weight in rows
        ]
    }


def replay_all(
    geometry: dict[str, Any], bank: component.ComponentBank, candidates: list[dict[str, Any]]
) -> dict[str, Any]:
    inner_constants, inner_charges, inner_refs = pool(bank, candidates, "inner")
    outer_constants, outer_charges, outer_refs = pool(bank, candidates, "outer")
    leaves = []
    for node in sorted(
        (row for row in geometry["nodes"].values() if row["state"] == "ACTIVE_LEAF"),
        key=lambda row: row["node_id"],
    ):
        profiles = stored.validate_leaf_payload(node)
        matrix = np.asarray(profiles, dtype=np.float64)
        normalization = coordinate.diagnostic_normalization(matrix, 262144)
        inner_values = inner_constants[None, :] - matrix @ inner_charges.T
        outer_values = outer_constants[None, :] - matrix @ outer_charges.T
        mixed = component.independent_minimax(
            inner_values, outer_values, normalization, one_shot.MIXTURE_DENOMINATOR
        )
        lp = pricing.solve_component_lp(inner_values, outer_values, normalization)
        dual_raw = -np.asarray(lp.ineqlin.marginals, dtype=np.float64)
        dual = pricing.rationalize_distribution(dual_raw, one_shot.DUAL_DENOMINATOR)
        upper = float(mixed["score"])
        leaves.append(
            {
                "node_id": node["node_id"],
                "h2_cell": int(node["h2_cell"]),
                "exact_count": node["exact_count"],
                "exact_vertex_count": node["exact_vertex_count"],
                "exact_vertex_sha256": node["exact_vertex_sha256"],
                "candidate_upper_log2": upper,
                "contribution_log2": math.log2(int(node["exact_count"])) + upper,
                "inner_selector": selector(mixed["inner"], inner_refs),
                "outer_selector": selector(mixed["outer"], outer_refs),
                "vertex_dual": [
                    {"vertex": index, "weight": pricing.fraction_text(weight)}
                    for index, weight in dual
                ],
                "binary64_primal_score": float(lp.fun),
                "rationalization_loss_bits": float(mixed["rationalization_loss_bits"]),
            }
        )
    aggregate = stored.log2_sum([row["contribution_log2"] for row in leaves])
    worst = max(leaves, key=lambda row: (row["contribution_log2"], row["node_id"]))
    return {
        "leaves": leaves,
        "summary": {
            "active_leaf_count": len(leaves),
            "inner_component_count": len(inner_refs),
            "outer_component_count": len(outer_refs),
            "aggregate_log2_union_diagnostic": aggregate,
            "worst_leaf_id": worst["node_id"],
            "worst_contribution_log2": worst["contribution_log2"],
        },
    }


def diversified_targets(
    replay: dict[str, Any], geometry: dict[str, Any], inner_pool, outer_pool
) -> list[dict[str, Any]]:
    nodes = geometry["nodes"]
    ranked = distinct_h2_ranking(replay)
    selected = []
    inner_constants, inner_charges, _inner_refs = inner_pool
    outer_constants, outer_charges, _outer_refs = outer_pool
    for row in ranked:
        node = nodes[row["node_id"]]
        profiles = stored.validate_leaf_payload(node)
        dual = [(item["vertex"], Fraction(item["weight"])) for item in row["vertex_dual"]]
        barycenter = pricing.exact_barycenter(profiles, dual)
        dual_float = np.zeros(len(profiles))
        for index, weight in dual:
            dual_float[index] = float(weight)
        matrix = np.asarray(profiles, dtype=np.float64)
        inner_values = inner_constants[None, :] - matrix @ inner_charges.T
        outer_values = outer_constants[None, :] - matrix @ outer_charges.T
        selected.append(
            {
                "node_id": row["node_id"],
                "h2_cell": row["h2_cell"],
                "contribution_log2": row["contribution_log2"],
                "exact_vertex_sha256": row["exact_vertex_sha256"],
                "profiles": [list(profile) for profile in profiles],
                "dual": [(index, pricing.fraction_text(weight)) for index, weight in dual],
                "barycenter_exact": [pricing.fraction_text(value) for value in barycenter],
                "barycenter_float": one_shot.fractional_profile(
                    [pricing.fraction_text(value) for value in barycenter]
                ),
                "inner_mu": float(np.min(dual_float @ inner_values)),
                "outer_mu": float(np.min(dual_float @ outer_values)),
            }
        )
        if len(selected) == one_shot.TARGETS:
            return selected
    raise RuntimeError("could not select eight diversified iterative targets")


def distinct_h2_ranking(replay: dict[str, Any]) -> list[dict[str, Any]]:
    """Rank leaves by contribution, retaining only the first leaf per h2 parent."""
    ranked = sorted(
        replay["leaves"], key=lambda row: (-row["contribution_log2"], row["node_id"])
    )
    selected = []
    seen = set()
    for row in ranked:
        if row["h2_cell"] in seen:
            continue
        seen.add(row["h2_cell"])
        selected.append(row)
    return selected


def round_configuration(
    state: dict[str, Any], round_number: int, targets: list[dict[str, Any]],
    lanes: tuple[str, ...],
) -> str:
    return coordinate.sha256_bytes(
        coordinate.canonical_bytes(
            {
                "state_seed_sha256": state["seed_bindings"]["prior_catalogue_sha256"],
                "round": round_number,
                "lanes": list(lanes),
                "prior_parameter_digests": sorted(
                    row["parameter_digest"] for row in state["candidates"]
                ),
                "targets": [(row["node_id"], row["barycenter_exact"]) for row in targets],
                "tuning": {
                    "inner_poles": [0.2, 0.5, 0.8],
                    "inner_exponents": [0.25, 0.5, 0.75],
                    "inner_screen_iterations": 1,
                    "inner_final_iterations": 3,
                    "outer_mode": "asymmetric",
                    "outer_max_iterations": 60,
                    "outer_max_evaluations": 800,
                    "tau_price_bits": one_shot.TAU_PRICE,
                },
            }
        )
    )


def round_checkpoint_path(work_dir: Path, round_number: int, task: dict[str, Any]) -> Path:
    return work_dir / f"round_{round_number:03d}" / f"{task['lane']}_{task['h2_cell']:03d}.json"


def load_matching_lane_checkpoint(
    path: Path, configuration: str, task: dict[str, Any]
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    row = json.loads(path.read_bytes())
    if (
        row.get("configuration_digest") != configuration
        or row.get("lane") != task["lane"]
        or row.get("node_id") != task["node_id"]
        or row.get("barycenter_exact") != task["barycenter_exact"]
    ):
        return None
    return row


def run_one_round(
    state: dict[str, Any], geometry: dict[str, Any], bank: component.ComponentBank,
    work_dir: Path, lanes: tuple[str, ...],
) -> dict[str, Any]:
    number = int(state["completed_rounds"]) + 1
    pre = replay_all(geometry, bank, state["candidates"])
    inner_pool = pool(bank, state["candidates"], "inner")
    outer_pool = pool(bank, state["candidates"], "outer")
    targets = diversified_targets(pre, geometry, inner_pool, outer_pool)
    configuration = round_configuration(state, number, targets, lanes)
    tasks = [{**target, "lane": lane} for target in targets for lane in lanes]
    completed = []
    missing = []
    for task in tasks:
        path = round_checkpoint_path(work_dir, number, task)
        row = load_matching_lane_checkpoint(path, configuration, task)
        if row is None:
            missing.append(task)
        else:
            completed.append(row)
    # Deliberately reuse the bounded worker implementation and its eight-worker cap.
    if missing:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=one_shot.WORKERS) as executor:
            futures = {executor.submit(one_shot.tune_job, task): task for task in missing}
            for future in as_completed(futures):
                row = future.result()
                row["configuration_digest"] = configuration
                completed.append(row)
                one_shot.atomic_json(
                    round_checkpoint_path(work_dir, number, row), row
                )
    completed.sort(key=lambda row: (row["h2_cell"], row["lane"]))
    target_by_node = {row["node_id"]: row for row in targets}
    seen_parameters = {row["parameter_digest"] for row in state["candidates"]}
    new_rows = []
    for row in completed:
        row["affine_barycenter_identity_error_bits"] = one_shot.identity_error(
            target_by_node[row["node_id"]], row
        )
        row["support_eligible"] = all(
            math.isfinite(float(value)) for value in row["parameters"]["charge_log2"]
        )
        row["accepted"] = (
            row["support_eligible"]
            and row["affine_barycenter_identity_error_bits"] <= 1e-6
            and row["reduced_cost"] < -one_shot.TAU_PRICE
            and row["parameter_digest"] not in seen_parameters
        )
        if row["accepted"]:
            seen_parameters.add(row["parameter_digest"])
            new_rows.append(row)
    round_catalogue = {
        "schema": "permute-conv.packet-group-g8-factorized-round-catalogue.v1",
        "round": number,
        "configuration_digest": configuration,
        "lanes": list(lanes),
        "rows": completed,
    }
    catalogue_path = work_dir / f"round_{number:03d}" / "candidate_catalogue.json"
    catalogue_digest = one_shot.atomic_json(catalogue_path, round_catalogue)
    for ordinal, row in enumerate(new_rows):
        row["reference"] = candidate_reference(row, catalogue_digest, ordinal)
        row["origin_catalogue_sha256"] = catalogue_digest
    state["candidates"].extend(copy.deepcopy(new_rows))
    post = replay_all(geometry, bank, state["candidates"])
    selected_ids = {row["node_id"] for row in targets}
    before_by_id = {row["node_id"]: row for row in pre["leaves"]}
    after_by_id = {row["node_id"]: row for row in post["leaves"]}
    selected_16k = sum(
        before_by_id[node_id]["contribution_log2"]
        - after_by_id[node_id]["contribution_log2"] >= 16384.0
        for node_id in selected_ids
    )
    aggregate_improvement = (
        pre["summary"]["aggregate_log2_union_diagnostic"]
        - post["summary"]["aggregate_log2_union_diagnostic"]
    )
    worst_improvement = (
        pre["summary"]["worst_contribution_log2"]
        - post["summary"]["worst_contribution_log2"]
    )
    new_inner_positive = any(
        component["component"]["source_sha256"] == catalogue_digest
        for leaf in post["leaves"] for component in leaf["inner_selector"]["components"]
    )
    new_outer_positive = any(
        component["component"]["source_sha256"] == catalogue_digest
        for leaf in post["leaves"] for component in leaf["outer_selector"]["components"]
    )
    no_regression = all(
        after_by_id[node_id]["candidate_upper_log2"]
        <= before_by_id[node_id]["candidate_upper_log2"] + 1e-8
        for node_id in before_by_id
    )
    required_lane_positive = (
        ("inner" not in lanes or new_inner_positive)
        and ("outer" not in lanes or new_outer_positive)
    )
    correctness = no_regression and required_lane_positive
    payoff = selected_16k >= 4 and (
        aggregate_improvement >= 32768.0 or worst_improvement >= 32768.0
    )
    record = {
        "schema": ROUND_SCHEMA,
        "round": number,
        "configuration_digest": configuration,
        "resumed_job_count": len(tasks) - len(missing),
        "new_job_count": len(missing),
        "targets": [
            {key: row[key] for key in ("node_id", "h2_cell", "contribution_log2", "barycenter_exact", "inner_mu", "outer_mu")}
            for row in targets
        ],
        "lane_results": completed,
        "candidate_catalogue": {
            "path": str(catalogue_path), "sha256": catalogue_digest
        },
        "accepted_inner_count": sum(row["accepted"] and row["lane"] == "inner" for row in completed),
        "accepted_outer_count": sum(row["accepted"] and row["lane"] == "outer" for row in completed),
        "pre_replay": pre,
        "post_replay": post,
        "aggregate_improvement_bits": aggregate_improvement,
        "worst_improvement_bits": worst_improvement,
        "gate": {
            "required_lanes": list(lanes),
            "rational_replay_no_regression": no_regression,
            "new_inner_positive": new_inner_positive,
            "new_outer_positive": new_outer_positive,
            "selected_leaves_improved_16384": selected_16k,
            "correctness_passes": correctness,
            "payoff_passes": payoff,
            "overall_passes": correctness and payoff,
            "outward_replay_performed": False,
        },
    }
    return record


def initialize(args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, Any], component.ComponentBank]:
    bindings = {
        "geometry_sha256": coordinate.sha256_path(args.geometry),
        "prior_main_sha256": coordinate.sha256_path(args.prior_main),
        "prior_catalogue_sha256": coordinate.sha256_path(args.prior_catalogue),
    }
    if bindings != {
        "geometry_sha256": GEOMETRY_SHA256,
        "prior_main_sha256": PRIOR_MAIN_SHA256,
        "prior_catalogue_sha256": PRIOR_CATALOGUE_SHA256,
    }:
        raise ValueError("iterative factorized seed binding changed")
    geometry = json.loads(args.geometry.read_bytes())
    prior = json.loads(args.prior_main.read_bytes())
    catalogue = json.loads(args.prior_catalogue.read_bytes())
    bank, bank_bindings = component.load_component_bank(
        args.manifest, args.atlas, args.old_supplementary, args.new_supplementary
    )
    state = {
        "schema": SCHEMA,
        "status": "RUNNABLE_DIAGNOSTIC_FACTORIZED_ROUNDS",
        "seed_bindings": {**bindings, **bank_bindings},
        "configuration": {
            "targets_per_round": one_shot.TARGETS,
            "workers": one_shot.WORKERS,
            "tau_price_bits": one_shot.TAU_PRICE,
            "fixed_active_leaf_count": 189,
            "lanes": list(args.lanes),
            "jobs_per_round": one_shot.TARGETS * len(args.lanes),
        },
        "completed_rounds": 0,
        "candidates": initial_candidates(catalogue, bindings["prior_catalogue_sha256"]),
        "rounds": [],
        "current_replay": prior["replay"],
        "stop_reason": None,
        "scope_limit": "Binary64 discovery and rational replay only; no outward claim.",
    }
    return state, geometry, bank


def validate_resume(state: dict[str, Any], fresh: dict[str, Any]) -> None:
    if state.get("schema") != SCHEMA:
        raise ValueError("unexpected iterative checkpoint schema")
    for key in ("seed_bindings", "configuration"):
        if state.get(key) != fresh.get(key):
            raise ValueError(f"iterative checkpoint {key} changed")
    if int(state.get("completed_rounds", -1)) != len(state.get("rounds", [])):
        raise ValueError("iterative checkpoint round history is malformed")


def run(args: argparse.Namespace) -> dict[str, Any]:
    fresh, geometry, bank = initialize(args)
    if args.resume_from:
        state = json.loads(args.resume_from.read_bytes())
        validate_resume(state, fresh)
        if state.get("stop_reason") == "MAX_ROUNDS" and state["completed_rounds"] < args.max_rounds:
            state["stop_reason"] = None
            state["status"] = "RUNNABLE_DIAGNOSTIC_FACTORIZED_ROUNDS"
    else:
        state = fresh
    while state["completed_rounds"] < args.max_rounds and state["stop_reason"] is None:
        started = time.perf_counter()
        record = run_one_round(state, geometry, bank, args.work_dir, args.lanes)
        record["runtime_seconds"] = time.perf_counter() - started
        state["rounds"].append(record)
        state["completed_rounds"] += 1
        state["current_replay"] = record["post_replay"]
        if record["accepted_inner_count"] + record["accepted_outer_count"] == 0:
            state["stop_reason"] = "NO_ACCEPTED_COLUMN"
            state["status"] = "STOPPED_DIAGNOSTIC_NO_ACCEPTED_COLUMN"
        elif not record["gate"]["overall_passes"]:
            state["stop_reason"] = "ROUND_GATE_FAILED"
            state["status"] = "STOPPED_DIAGNOSTIC_ROUND_GATE"
        one_shot.atomic_json(args.checkpoint, state)
    if state["stop_reason"] is None and state["completed_rounds"] >= args.max_rounds:
        state["stop_reason"] = "MAX_ROUNDS"
        state["status"] = "STOPPED_DIAGNOSTIC_MAX_ROUNDS"
    one_shot.atomic_json(args.checkpoint, state)
    return state


def synthetic_state() -> dict[str, Any]:
    return {
        "completed_rounds": 0,
        "candidates": [{"parameter_digest": "seed", "lane": "inner"}],
        "rounds": [],
        "current_replay": {"summary": {"aggregate_log2_union_diagnostic": 100.0}},
    }


def synthetic_round(state: dict[str, Any], number: int) -> None:
    row = {
        "parameter_digest": f"row-{number}", "lane": "inner", "accepted": True,
        "parameters": {"constant_log2": float(number), "charge_log2": [0.0] * 9},
    }
    state["candidates"].append(row)
    record = {"round": number, "accepted_inner_count": 1, "accepted_outer_count": 0}
    state["rounds"].append(record)
    state["completed_rounds"] += 1


def self_test() -> None:
    if LANE_POLICIES != {
        "both": ("inner", "outer"),
        "inner": ("inner",),
        "outer": ("outer",),
    }:
        raise SystemExit("iterative CLI lane policy mapping failed")
    uninterrupted = synthetic_state()
    synthetic_round(uninterrupted, 1); synthetic_round(uninterrupted, 2)
    resumed = synthetic_state(); synthetic_round(resumed, 1)
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "checkpoint.json"
        one_shot.atomic_json(path, resumed)
        resumed = json.loads(path.read_bytes())
        synthetic_round(resumed, 2)
        one_shot.atomic_json(path, resumed)
        if json.loads(path.read_bytes()) != uninterrupted:
            raise SystemExit("iterative factorized resume equivalence failed")
    task = {"lane": "outer", "node_id": "n", "h2_cell": 1, "barycenter_exact": ["1/1"]}
    row = {**task, "configuration_digest": "x"}
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "lane.json"
        one_shot.atomic_json(path, row)
        if load_matching_lane_checkpoint(path, "x", task) != row:
            raise SystemExit("iterative lane checkpoint matching failed")
        wrong_lane = {**task, "lane": "inner"}
        if load_matching_lane_checkpoint(path, "x", wrong_lane) is not None:
            raise SystemExit("iterative lane checkpoint crossed lane policy")
        if load_matching_lane_checkpoint(path, "different", task) is not None:
            raise SystemExit("iterative lane checkpoint crossed configuration")
    prior = json.loads(PRIOR_MAIN.read_bytes())
    target_cells = [
        row["h2_cell"] for row in distinct_h2_ranking(prior["replay"])[
            :one_shot.TARGETS
        ]
    ]
    if target_cells != [146, 102, 72, 63, 71, 101, 108, 22]:
        raise SystemExit(f"iterative audited target ranking changed: {target_cells}")
    targets = [{"node_id": f"h2:{cell:03d}", "h2_cell": cell} for cell in target_cells]
    lanes = ("outer",)
    tasks = [{**target, "lane": lane} for target in targets for lane in lanes]
    if (
        len(tasks) != 8
        or any(task["lane"] != "outer" for task in tasks)
        or [task["h2_cell"] for task in tasks] != [146,102,72,63,71,101,108,22]
    ):
        raise SystemExit("iterative outer-only task construction failed")
    catalogue = [{"parameter_digest": "seed-inner", "lane": "inner"}]
    prefix = copy.deepcopy(catalogue)
    catalogue.extend(
        {"parameter_digest": f"outer-{task['h2_cell']:03d}", "lane": task["lane"]}
        for task in tasks
    )
    if catalogue[:len(prefix)] != prefix or any(
        row["lane"] != "outer" for row in catalogue[len(prefix):]
    ):
        raise SystemExit("iterative outer-only catalogue append failed")
    print("two_round_uninterrupted_equals_resume=PASS")
    print("atomic_outer_lane_checkpoint_match_and_policy_rejection=PASS")
    print("outer_only_append_only_candidate_catalogue=PASS")
    print("outer_only_eight_job_policy=PASS")
    print("cli_lane_policies_both_inner_outer=PASS")
    print("audited_target_order_146_102_072_063_071_101_108_022=PASS")
    print("status=PASS_G8_ITERATIVE_FACTORIZED_COLUMNS_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=component.h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=component.h2.DEFAULT_ATLAS)
    parser.add_argument("--old-supplementary", type=Path, default=component.h2.DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--new-supplementary", type=Path, default=component.NEW_ATLAS)
    parser.add_argument("--geometry", type=Path, default=GEOMETRY)
    parser.add_argument("--prior-main", type=Path, default=PRIOR_MAIN)
    parser.add_argument("--prior-catalogue", type=Path, default=PRIOR_CATALOGUE)
    parser.add_argument("--resume-from", type=Path)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--max-rounds", type=int)
    parser.add_argument(
        "--lanes", choices=tuple(LANE_POLICIES), default="both",
        help="component-tuning lanes for each round (default: both)",
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(); return
    if args.max_rounds is None or args.max_rounds < 0:
        parser.error("--max-rounds must be a nonnegative integer")
    args.lanes = LANE_POLICIES[args.lanes]
    state = run(args)
    print(f"completed_rounds={state['completed_rounds']}")
    print(f"candidate_count={len(state['candidates'])}")
    print(f"aggregate_log2={state['current_replay']['summary']['aggregate_log2_union_diagnostic']:.12f}")
    print(f"stop_reason={state['stop_reason']}")
    print(f"checkpoint={args.checkpoint}")


if __name__ == "__main__":
    main()
