#!/usr/bin/env python3
"""Retune one fixed-profile witness from frozen neighboring witness tilts."""

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

from packet_group_outer_profile import D
from packet_group_profile_bound import (
    inner_probability,
    refine_inner_full_fugacities,
    split_cap_table,
)


_SPLIT_CAPS = None


def artifact_rows(path: Path) -> list[dict]:
    value = json.loads(path.read_text(encoding="utf-8"))
    rows = value if isinstance(value, list) else value.get("rows")
    if not isinstance(rows, list):
        raise ValueError(f"reseed witness: no rows in {path}")
    return rows


def initialize_worker() -> None:
    global _SPLIT_CAPS
    _SPLIT_CAPS = split_cap_table()


def tune_one(task: tuple[int, dict, dict, int]) -> tuple[int, dict]:
    index, target, seed, passes = task
    if _SPLIT_CAPS is None:
        raise RuntimeError("reseed witness worker was not initialized")
    profile = [int(value) for value in target["profile"]]
    refined = refine_inner_full_fugacities(
        int(target["group_bits"]),
        profile,
        _SPLIT_CAPS,
        initial_pole=float(seed["pole"]),
        initial_fugacities=np.asarray(seed["fugacities"], dtype=np.float64),
        log_steps=(
            math.log(16.0),
            math.log(4.0),
            math.log(2.0),
            math.log(2.0**0.5),
            math.log(2.0**0.25),
            math.log(2.0**0.125),
            math.log(2.0**0.0625),
        ),
        passes=passes,
        witness_iterations=16,
    )
    inner_value, inner_details = inner_probability(
        int(target["group_bits"]),
        profile,
        refined[1],
        refined[2],
        _SPLIT_CAPS,
        iterations=40,
    )
    inner_details = {
        **inner_details,
        "full_coordinate_trace": refined[4],
        "full_coordinate_evaluations": refined[5],
        "full_coordinate_anchor": refined[6],
        "reseeded_from": seed["name"],
    }
    outer_charge = np.asarray(target["outer_charge"], dtype=np.float64)
    vector = np.asarray(profile, dtype=np.float64)
    fugacities = np.asarray(refined[2], dtype=np.float64)
    inner_constant = float(inner_details["inner_mgf_log2"]) - D * math.log2(
        float(refined[1])
    )
    charge = outer_charge + np.log2(fugacities)
    constant = float(target["outer_constant_log2"]) + inner_constant
    combined = constant - float(vector @ charge) - float(target["normalization_log2"])
    row = {
        **target,
        "name": f"{target['name']}__seed_{seed['name']}",
        "pole": float(refined[1]),
        "fugacities": fugacities.tolist(),
        "inner_constant_log2": inner_constant,
        "constant_log2": constant,
        "charge": charge.tolist(),
        "combined_log2": combined,
        "margin_bits": float(target["target_log2"]) - combined,
        "inner_probability_log2": inner_value,
        "inner_details": inner_details,
        "tuning_method": "reseeded_full_coordinate",
        "status": "DIAGNOSTIC_BINARY64_RESEEDED_FIXED_PROFILE_WITNESS",
    }
    return index, row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-artifact", type=Path, required=True)
    parser.add_argument("--target-name", action="append", required=True)
    parser.add_argument("--seed-artifact", type=Path, action="append", required=True)
    parser.add_argument("--seed-name", action="append", default=[])
    parser.add_argument("--passes", type=int, default=12)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.passes < 1 or args.workers < 1:
        parser.error("passes and workers must be positive")

    requested_targets = set(args.target_name)
    targets = [
        row
        for row in artifact_rows(args.target_artifact)
        if row.get("name") in requested_targets
    ]
    if {row.get("name") for row in targets} != requested_targets:
        raise ValueError("reseed witness: a target name is absent or ambiguous")
    seeds = [row for path in args.seed_artifact for row in artifact_rows(path)]
    if args.seed_name:
        selected = set(args.seed_name)
        seeds = [row for row in seeds if row.get("name") in selected]
        if {row.get("name") for row in seeds} != selected:
            raise ValueError("reseed witness: a requested seed name is missing")
    if not seeds:
        raise ValueError("reseed witness: no seed rows")

    tasks = [
        (index, target, seed, args.passes)
        for index, (target, seed) in enumerate(
            (target, seed) for target in targets for seed in seeds
        )
    ]
    rows: list[dict | None] = [None] * len(tasks)
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=min(args.workers, len(tasks)), initializer=initialize_worker
    ) as pool:
        for future in concurrent.futures.as_completed(
            [pool.submit(tune_one, task) for task in tasks]
        ):
            index, row = future.result()
            rows[index] = row
            print(
                f"seed={row['inner_details']['reseeded_from']} "
                f"combined={row['combined_log2']:.9f} margin={row['margin_bits']:.9f}",
                flush=True,
            )
    completed = [row for row in rows if row is not None]
    report = {
        "status": "DIAGNOSTIC_BINARY64_RESEEDED_FIXED_PROFILE_ATLAS",
        "group_bits": int(targets[0]["group_bits"]),
        "complete": len(completed) == len(tasks),
        "target_artifact": str(args.target_artifact),
        "target_names": args.target_name,
        "passes": args.passes,
        "rows": completed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
