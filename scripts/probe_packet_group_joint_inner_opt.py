#!/usr/bin/env python3
"""Jointly optimize a fixed-profile inner witness from frozen seed tilts."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
from pathlib import Path

for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

import numpy as np
from scipy.optimize import minimize

from packet_group_outer_profile import D
from packet_group_profile_bound import inner_probability, split_cap_table


_SPLIT_CAPS = None


def artifact_rows(path: Path) -> list[dict]:
    value = json.loads(path.read_text(encoding="utf-8"))
    rows = value if isinstance(value, list) else value.get("rows")
    if not isinstance(rows, list):
        raise ValueError(f"joint inner optimizer: no rows in {path}")
    return rows


def initialize_worker() -> None:
    global _SPLIT_CAPS
    _SPLIT_CAPS = split_cap_table()


def optimize_one(task: tuple[int, dict, dict, int, int, int]) -> tuple[int, dict]:
    index, target, seed, maxiter, maxfev, witness_iterations = task
    if _SPLIT_CAPS is None:
        raise RuntimeError("joint inner optimizer worker was not initialized")
    group_bits = int(target["group_bits"])
    profile = [int(value) for value in target["profile"]]
    profile_array = np.asarray(profile, dtype=np.int64)
    if np.any(profile_array <= 0):
        raise ValueError("joint inner optimizer currently requires full support")
    anchor = int(np.argmax(profile_array))
    movable = [coordinate for coordinate in range(group_bits + 1) if coordinate != anchor]
    seed_fugacities = np.maximum(
        np.asarray(seed["fugacities"], dtype=np.float64), 1e-300
    )
    seed_fugacities /= seed_fugacities[anchor]
    initial = np.asarray(
        [math.log(float(seed_fugacities[i])) for i in movable]
        + [math.log(float(seed["pole"]))],
        dtype=np.float64,
    )
    bounds = [(-50.0, 50.0)] * len(movable) + [
        (math.log(0.01), math.log(0.95))
    ]
    evaluations = 0

    def unpack(point: np.ndarray) -> tuple[float, np.ndarray]:
        fugacities = np.ones(group_bits + 1, dtype=np.float64)
        for coordinate, value in zip(movable, point[:-1]):
            fugacities[coordinate] = math.exp(float(value))
        return math.exp(float(point[-1])), fugacities

    def objective(point: np.ndarray) -> float:
        nonlocal evaluations
        evaluations += 1
        pole, fugacities = unpack(point)
        probability, _details = inner_probability(
            group_bits,
            profile,
            pole,
            fugacities,
            _SPLIT_CAPS,
            iterations=witness_iterations,
        )
        return float(probability)

    result = minimize(
        objective,
        initial,
        method="Powell",
        bounds=bounds,
        options={
            "maxiter": maxiter,
            "maxfev": maxfev,
            "xtol": 2e-5,
            "ftol": 1e-9,
        },
    )
    candidates = [(float(result.fun), np.asarray(result.x, dtype=np.float64))]
    candidates.append((objective(initial), initial))
    _, best_point = min(candidates, key=lambda row: row[0])
    pole, fugacities = unpack(best_point)
    inner_value, inner_details = inner_probability(
        group_bits,
        profile,
        pole,
        fugacities,
        _SPLIT_CAPS,
        iterations=max(160, witness_iterations),
    )
    inner_details = {
        **inner_details,
        "joint_optimizer": "scipy_powell",
        "joint_optimizer_success": bool(result.success),
        "joint_optimizer_message": str(result.message),
        "joint_optimizer_iterations": int(result.nit),
        "joint_optimizer_evaluations": evaluations,
        "joint_optimizer_seed": seed["name"],
    }
    outer_charge = np.asarray(target["outer_charge"], dtype=np.float64)
    vector = np.asarray(profile, dtype=np.float64)
    inner_constant = float(inner_details["inner_mgf_log2"]) - D * math.log2(pole)
    charge = outer_charge + np.log2(fugacities)
    constant = float(target["outer_constant_log2"]) + inner_constant
    combined = constant - float(vector @ charge) - float(target["normalization_log2"])
    row = {
        **target,
        "name": f"{target['name']}__powell_{seed['name']}",
        "pole": pole,
        "fugacities": fugacities.tolist(),
        "inner_constant_log2": inner_constant,
        "constant_log2": constant,
        "charge": charge.tolist(),
        "combined_log2": combined,
        "margin_bits": float(target["target_log2"]) - combined,
        "inner_probability_log2": inner_value,
        "inner_details": inner_details,
        "tuning_method": "joint_powell_log_fugacities",
        "status": "DIAGNOSTIC_BINARY64_JOINT_FIXED_PROFILE_WITNESS",
    }
    return index, row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-artifact", type=Path, required=True)
    parser.add_argument("--target-name", action="append", default=[])
    parser.add_argument(
        "--max-targets",
        type=int,
        default=0,
        help="when target names are omitted, use at most this many leading rows; 0 means all",
    )
    parser.add_argument("--seed-artifact", type=Path, action="append", required=True)
    parser.add_argument("--seed-name", action="append", default=[])
    parser.add_argument(
        "--match-seed-profile",
        action="store_true",
        help="pair each target only with seeds having the identical profile",
    )
    parser.add_argument("--maxiter", type=int, default=16)
    parser.add_argument(
        "--maxfev",
        type=int,
        default=400,
        help="hard cap on objective evaluations per target",
    )
    parser.add_argument("--witness-iterations", type=int, default=96)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        args.maxiter < 1
        or args.maxfev < 1
        or args.witness_iterations < 1
        or args.workers < 1
        or args.max_targets < 0
    ):
        parser.error(
            "maxiter, maxfev, witness iterations, and workers must be positive; "
            "max targets must be nonnegative"
        )

    requested = set(args.target_name)
    target_rows = artifact_rows(args.target_artifact)
    if requested:
        targets = [row for row in target_rows if row.get("name") in requested]
        if {row.get("name") for row in targets} != requested:
            raise ValueError("joint inner optimizer: target name is absent or ambiguous")
    else:
        targets = target_rows[: args.max_targets or None]
    if not targets:
        raise ValueError("joint inner optimizer: no target rows selected")
    seeds = [row for path in args.seed_artifact for row in artifact_rows(path)]
    if args.seed_name:
        selected = set(args.seed_name)
        seeds = [row for row in seeds if row.get("name") in selected]
        if {row.get("name") for row in seeds} != selected:
            raise ValueError("joint inner optimizer: requested seed name is missing")
    if not seeds:
        raise ValueError("joint inner optimizer: no seeds")

    pairs = [
        (target, seed)
        for target in targets
        for seed in seeds
        if not args.match_seed_profile
        or tuple(map(int, target["profile"])) == tuple(map(int, seed["profile"]))
    ]
    if args.match_seed_profile:
        matched = {tuple(map(int, target["profile"])) for target, _seed in pairs}
        missing = [
            target["name"]
            for target in targets
            if tuple(map(int, target["profile"])) not in matched
        ]
        if missing:
            raise ValueError(
                "joint inner optimizer: no same-profile seed for " + ", ".join(missing)
            )
    tasks = [
        (index, target, seed, args.maxiter, args.maxfev, args.witness_iterations)
        for index, (target, seed) in enumerate(pairs)
    ]
    rows: list[dict | None] = [None] * len(tasks)
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=min(args.workers, len(tasks)), initializer=initialize_worker
    ) as pool:
        futures = [pool.submit(optimize_one, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            index, row = future.result()
            rows[index] = row
            print(
                f"target={row['profile']} seed={row['inner_details']['joint_optimizer_seed']} "
                f"combined={row['combined_log2']:.9f} margin={row['margin_bits']:.9f}",
                flush=True,
            )
    completed = [row for row in rows if row is not None]
    report = {
        "status": "DIAGNOSTIC_BINARY64_JOINT_FIXED_PROFILE_ATLAS",
        "group_bits": int(targets[0]["group_bits"]),
        "complete": len(completed) == len(tasks),
        "target_artifact": str(args.target_artifact),
        "target_names": args.target_name,
        "maxiter": args.maxiter,
        "witness_iterations": args.witness_iterations,
        "match_seed_profile": args.match_seed_profile,
        "rows": completed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
