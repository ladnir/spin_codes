#!/usr/bin/env python3
"""Retune each residual profile from its named strongest inner witness."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
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

from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import D, N, normalization_log2
from packet_group_profile_bound import inner_probability, split_cap_table


_SPLIT_CAPS = None


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact_rows(path: Path) -> list[dict]:
    value = json.loads(path.read_text(encoding="utf-8"))
    rows = value if isinstance(value, list) else value.get("rows")
    if not isinstance(rows, list):
        raise ValueError(f"seeded residual optimizer: no rows in {path}")
    return rows


def initialize_worker() -> None:
    global _SPLIT_CAPS
    _SPLIT_CAPS = split_cap_table()


def retarget_seed(seed: dict, residual: dict) -> dict:
    group_bits = int(seed["group_bits"])
    profile = [int(value) for value in residual["profile"]]
    if len(profile) != group_bits + 1 or any(value <= 0 for value in profile):
        raise ValueError("seeded residual optimizer requires a full-support profile")
    vector = np.asarray(profile, dtype=np.float64)
    normalization = float(normalization_log2(group_bits, vector)[0])
    target = -40.0 - math.log2(profile_count(group_bits, N))
    return {
        **seed,
        "name": str(residual.get("name", "residual")),
        "profile": profile,
        "physical_weight": sum(i * count for i, count in enumerate(profile)),
        "support": list(range(group_bits + 1)),
        "normalization_log2": normalization,
        "target_log2": target,
        "residual_seed_reference": str(residual["diagnostic_best_witness"]),
        "residual_initial_envelope_log2": float(
            residual["diagnostic_current_envelope_log2"]
        ),
        "residual_initial_gap_bits": float(residual["diagnostic_gap_bits"]),
        "source_cell_ids": residual.get("source_cell_ids", []),
    }


def optimize_one(
    task: tuple[int, dict, int, int, int],
) -> tuple[int, dict]:
    index, target, maxiter, maxfev, witness_iterations = task
    if _SPLIT_CAPS is None:
        raise RuntimeError("seeded residual optimizer worker was not initialized")
    group_bits = int(target["group_bits"])
    profile = [int(value) for value in target["profile"]]
    vector = np.asarray(profile, dtype=np.float64)
    anchor = int(np.argmax(vector))
    movable = [coordinate for coordinate in range(group_bits + 1) if coordinate != anchor]
    seed_fugacities = np.maximum(
        np.asarray(target["fugacities"], dtype=np.float64), 1e-300
    )
    seed_fugacities /= seed_fugacities[anchor]
    initial = np.asarray(
        [math.log(float(seed_fugacities[i])) for i in movable]
        + [math.log(float(target["pole"]))],
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

    initial_value = objective(initial)
    result = minimize(
        objective,
        initial,
        method="Powell",
        bounds=bounds,
        options={
            "maxiter": maxiter,
            "maxfev": maxfev,
            "xtol": 1e-5,
            "ftol": 2e-10,
        },
    )
    candidates = [
        (initial_value, initial),
        (float(result.fun), np.asarray(result.x, dtype=np.float64)),
    ]
    _screen_value, best_point = min(candidates, key=lambda item: item[0])
    pole, fugacities = unpack(best_point)
    final_iterations = max(320, witness_iterations)
    inner_value, inner_details = inner_probability(
        group_bits,
        profile,
        pole,
        fugacities,
        _SPLIT_CAPS,
        iterations=final_iterations,
    )
    inner_details = {
        **inner_details,
        "seeded_optimizer": "scipy_powell_log_fugacities",
        "seeded_optimizer_success": bool(result.success),
        "seeded_optimizer_message": str(result.message),
        "seeded_optimizer_iterations": int(result.nit),
        "seeded_optimizer_evaluations": evaluations,
        "seeded_optimizer_objective_iterations": witness_iterations,
        "seeded_optimizer_seed_reference": target["residual_seed_reference"],
        "seeded_optimizer_initial_inner_log2": initial_value,
    }
    inner_constant = float(inner_details["inner_mgf_log2"]) - D * math.log2(pole)
    outer_charge = np.asarray(target["outer_charge"], dtype=np.float64)
    outer_constant = float(target["outer_constant_log2"])
    charge = outer_charge + np.log2(fugacities)
    constant = outer_constant + inner_constant
    combined = (
        constant
        - float(vector @ charge)
        - float(target["normalization_log2"])
    )
    return index, {
        **target,
        "name": str(target["name"]) + "__seeded_powell",
        "pole": pole,
        "fugacities": fugacities.tolist(),
        "inner_constant_log2": inner_constant,
        "inner_probability_log2": inner_value,
        "inner_details": inner_details,
        "constant_log2": constant,
        "charge": charge.tolist(),
        "combined_log2": combined,
        "margin_bits": float(target["target_log2"]) - combined,
        "tuning_method": "residual_named_seed_powell",
        "status": "DIAGNOSTIC_BINARY64_SEEDED_RESIDUAL_INNER_WITNESS",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--residual-ledger", type=Path, required=True)
    parser.add_argument("--seed-artifact", type=Path, action="append", required=True)
    parser.add_argument("--max-residuals", type=int, default=8)
    parser.add_argument("--maxiter", type=int, default=20)
    parser.add_argument("--maxfev", type=int, default=600)
    parser.add_argument("--witness-iterations", type=int, default=160)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        args.max_residuals < 1
        or args.maxiter < 1
        or args.maxfev < 1
        or args.witness_iterations < 1
        or args.workers < 1
        or args.workers > 8
    ):
        parser.error("positive parameters required; workers must be in [1,8]")

    seed_index: dict[str, dict] = {}
    sources = []
    for path in args.seed_artifact:
        rows = artifact_rows(path)
        for index, row in enumerate(rows):
            reference = f"{path.name}:{index}"
            if reference in seed_index:
                raise ValueError(f"duplicate seed reference: {reference}")
            seed_index[reference] = row
        sources.append({"path": str(path), "sha256": digest(path)})

    ledger = json.loads(args.residual_ledger.read_text(encoding="utf-8"))
    residuals = ledger.get("uncovered_integer_residuals", [])[: args.max_residuals]
    if not residuals:
        raise ValueError("seeded residual optimizer: no residuals selected")
    targets = []
    for residual in residuals:
        reference = str(residual["diagnostic_best_witness"])
        if reference not in seed_index:
            raise ValueError(f"missing named residual seed: {reference}")
        targets.append(retarget_seed(seed_index[reference], residual))

    tasks = [
        (index, target, args.maxiter, args.maxfev, args.witness_iterations)
        for index, target in enumerate(targets)
    ]
    rows: list[dict | None] = [None] * len(tasks)
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=min(args.workers, len(tasks)), initializer=initialize_worker
    ) as pool:
        for future in concurrent.futures.as_completed(
            pool.submit(optimize_one, task) for task in tasks
        ):
            index, row = future.result()
            rows[index] = row
            print(
                f"residual={index + 1}/{len(tasks)} profile={row['profile']} "
                f"inner={row['inner_probability_log2']:.9f} "
                f"seeded_combined={row['combined_log2']:.9f}",
                flush=True,
            )

    report = {
        "schema": "packet-group-g4-seeded-residual-inner-witnesses-v1",
        "status": "DIAGNOSTIC_BINARY64_SEEDED_RESIDUAL_INNER_ATLAS",
        "group_bits": 4,
        "complete": all(row is not None for row in rows),
        "parameters": {
            "maxiter": args.maxiter,
            "maxfev": args.maxfev,
            "witness_iterations": args.witness_iterations,
            "workers": args.workers,
        },
        "residual_ledger": {
            "path": str(args.residual_ledger),
            "sha256": digest(args.residual_ledger),
        },
        "sources": sources,
        "rows": [row for row in rows if row is not None],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
