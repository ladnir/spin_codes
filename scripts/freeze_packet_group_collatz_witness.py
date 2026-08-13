#!/usr/bin/env python3
"""Freeze the best finite-block Collatz vector in existing witness rows."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np

from packet_group_outer_profile import D
from packet_group_profile_bound import inner_probability, split_cap_table


_SPLIT_CAPS = None


def rows_from_artifact(value):
    rows = value if isinstance(value, list) else value.get("rows")
    if not isinstance(rows, list):
        raise ValueError("Collatz freezer requires an artifact with witness rows")
    return rows


def freeze_row(row: dict, split_caps, iterations: int) -> dict:
    profile = np.asarray(row["profile"], dtype=np.int64)
    fugacities = np.asarray(row["fugacities"], dtype=np.float64)
    pole = float(row["pole"])
    probability, details = inner_probability(
        int(row["group_bits"]),
        profile.tolist(),
        pole,
        fugacities,
        split_caps,
        iterations=iterations,
    )
    normalization = float(row["normalization_log2"])
    if abs(float(details["normalization_log2"]) - normalization) > 1e-7:
        raise ValueError("stored and recomputed profile normalizations disagree")
    inner_constant = float(details["inner_mgf_log2"]) - D * math.log2(pole)
    outer_charge = np.asarray(row["outer_charge"], dtype=np.float64)
    inner_charge = np.log2(fugacities)
    constant = float(row["outer_constant_log2"]) + inner_constant
    charge = outer_charge + inner_charge
    combined = constant - float(profile @ charge) - normalization
    direct_combined = (
        float(row["outer_constant_log2"])
        - float(profile @ outer_charge)
        + probability
    )
    if abs(combined - direct_combined) > 2e-7:
        raise ValueError("combined witness reconstruction is inconsistent")
    return {
        **row,
        "inner_constant_log2": inner_constant,
        "inner_probability_log2": probability,
        "inner_details": details,
        "collatz_vector": details["collatz_vector"],
        "constant_log2": constant,
        "charge": charge.tolist(),
        "combined_log2": combined,
        "margin_bits": float(row["target_log2"]) - combined,
        "collatz_freeze_iterations": iterations,
        "status": "DIAGNOSTIC_BINARY64_FROZEN_BEST_COLLATZ_VECTOR",
    }


def initialize_worker() -> None:
    global _SPLIT_CAPS
    _SPLIT_CAPS = split_cap_table()


def freeze_task(task: tuple[int, dict, int]) -> tuple[int, dict]:
    index, row, iterations = task
    if _SPLIT_CAPS is None:
        raise RuntimeError("Collatz freezer worker was not initialized")
    return index, freeze_row(row, _SPLIT_CAPS, iterations)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--name", action="append", default=[])
    parser.add_argument("--iterations", type=int, default=320)
    parser.add_argument(
        "--workers", type=int, default=max(1, min(16, os.cpu_count() or 1))
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.iterations < 1 or args.workers < 1:
        parser.error("iterations and workers must be positive")

    rows = []
    sources = []
    for path in args.input:
        artifact = json.loads(path.read_text(encoding="utf-8"))
        rows.extend(rows_from_artifact(artifact))
        sources.append(
            {
                "path": str(path.resolve()),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    requested = set(args.name)
    selected = [row for row in rows if not requested or row.get("name") in requested]
    if requested and requested != {row.get("name") for row in selected}:
        raise ValueError("one or more requested witness names are absent")
    tasks = [(index, row, args.iterations) for index, row in enumerate(selected)]
    if args.workers == 1:
        initialize_worker()
        results = [freeze_task(task) for task in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers, initializer=initialize_worker
        ) as executor:
            results = list(executor.map(freeze_task, tasks, chunksize=1))
    frozen = [row for _index, row in sorted(results)]
    report = {
        "schema": "packet-group-frozen-best-collatz-witnesses-v1",
        "sources": sources,
        "iterations": args.iterations,
        "workers": args.workers,
        "rows": frozen,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summaries = [
        {
            "name": row["name"],
            "best_iteration": row["inner_details"]["collatz_best_iteration"],
            "combined_log2": row["combined_log2"],
            "margin_bits": row["margin_bits"],
        }
        for row in frozen
    ]
    if len(summaries) <= 20:
        printed = summaries
    else:
        selected_iterations = [row["best_iteration"] for row in summaries]
        printed = {
            "rows": len(summaries),
            "minimum_best_iteration": min(selected_iterations),
            "maximum_best_iteration": max(selected_iterations),
            "mean_best_iteration": sum(selected_iterations) / len(selected_iterations),
            "minimum_margin_bits": min(row["margin_bits"] for row in summaries),
            "maximum_margin_bits": max(row["margin_bits"] for row in summaries),
        }
    print(json.dumps(printed, indent=2))


if __name__ == "__main__":
    main()
