#!/usr/bin/env python3
"""Tune a batch of generalized packet-profile witnesses in worker processes."""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# One proof-search worker owns one numerical stream.  Do not let NumPy/SciPy
# silently multiply that by a full BLAS thread team.
for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

from packet_group_profile_bound import split_cap_table
from packet_group_outer_profile import atom_count
from probe_packet_group_fixed_atlas import parse_named_profile, tune_fixed_witness


_SPLIT_CAPS = None


def initialize_worker() -> None:
    global _SPLIT_CAPS
    _SPLIT_CAPS = split_cap_table()


def tune_one(task):
    (
        index,
        group_bits,
        name,
        profile,
        screen_iterations,
        coordinate_iterations,
        exact_graph_puncture,
        full_coordinate_passes,
        full_coordinate_fine_steps,
    ) = task
    if _SPLIT_CAPS is None:
        raise RuntimeError("parallel fixed-atlas worker was not initialized")
    row = tune_fixed_witness(
        group_bits,
        name,
        profile,
        _SPLIT_CAPS,
        screen_iterations=screen_iterations,
        coordinate_iterations=coordinate_iterations,
        exact_graph_puncture=exact_graph_puncture,
        full_coordinate_passes=full_coordinate_passes,
        full_coordinate_fine_steps=full_coordinate_fine_steps,
    )
    return index, row


def write_checkpoint(path: Path, group_bits: int, rows) -> None:
    completed = [row for row in rows if row is not None]
    report = {
        "status": "DIAGNOSTIC_BINARY64_PARALLEL_FIXED_PROFILE_ATLAS",
        "group_bits": group_bits,
        "complete": len(completed) == len(rows),
        "completed_witnesses": len(completed),
        "requested_witnesses": len(rows),
        "rows": completed,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--profile", action="append", default=[])
    parser.add_argument(
        "--all-support-uniform",
        action="store_true",
        help="add one nearly uniform exact profile for every feasible nonzero class support",
    )
    parser.add_argument("--residual-ledger", type=Path)
    parser.add_argument("--max-residuals", type=int, default=16)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument("--screen-iterations", type=int, default=8)
    parser.add_argument("--coordinate-iterations", type=int, default=3)
    parser.add_argument("--full-coordinate-passes", type=int, default=2)
    parser.add_argument("--full-coordinate-fine-steps", action="store_true")
    parser.add_argument(
        "--exact-graph-puncture",
        action="store_true",
        help="enable the exact graph-spectrum/puncture-averaged g=4 outer branch",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers <= 0:
        raise SystemExit("parallel fixed atlas: workers must be positive")
    if args.max_residuals <= 0:
        raise SystemExit("parallel fixed atlas: max residuals must be positive")
    if (
        args.screen_iterations <= 0
        or args.coordinate_iterations <= 0
        or args.full_coordinate_passes <= 0
    ):
        raise SystemExit("parallel fixed atlas: iteration counts must be positive")

    profile_specs = list(args.profile)
    if args.all_support_uniform:
        atoms = atom_count(args.group)
        for mask in range(1, 1 << (args.group + 1)):
            support = [index for index in range(args.group + 1) if mask & (1 << index)]
            if support == [0]:
                continue
            quotient, remainder = divmod(atoms, len(support))
            profile = [0] * (args.group + 1)
            # Give the remainder to the heaviest active classes.  This is
            # deterministic and keeps the nonzero-weight feasibility check
            # maximally conservative near the clipped zero corner.
            for index in support:
                profile[index] = quotient
            for index in reversed(support[-remainder:] if remainder else []):
                profile[index] += 1
            if sum(index * count for index, count in enumerate(profile)) < 21:
                continue
            label = "_".join(map(str, support))
            profile_specs.append(f"support_{label}_uniform:{','.join(map(str, profile))}")
    if args.residual_ledger:
        ledger = json.loads(args.residual_ledger.read_text(encoding="utf-8"))
        residuals = ledger.get("uncovered_integer_residuals", [])
        for index, residual in enumerate(residuals[: args.max_residuals]):
            values = [int(value) for value in residual["profile"]]
            name = residual.get("name", f"residual_{index:04d}")
            profile_specs.append(f"{name}:{','.join(map(str, values))}")
    if not profile_specs:
        raise SystemExit("parallel fixed atlas: provide --profile or --residual-ledger")
    profiles = [parse_named_profile(text, args.group) for text in profile_specs]
    unique_profiles = []
    seen = set()
    for name, profile in profiles:
        key = tuple(profile)
        if key not in seen:
            seen.add(key)
            unique_profiles.append((name, profile))
    profiles = unique_profiles
    tasks = [
        (
            index,
            args.group,
            name,
            profile,
            args.screen_iterations,
            args.coordinate_iterations,
            args.exact_graph_puncture,
            args.full_coordinate_passes,
            args.full_coordinate_fine_steps,
        )
        for index, (name, profile) in enumerate(profiles)
    ]
    rows = [None] * len(tasks)
    with ProcessPoolExecutor(
        max_workers=min(args.workers, len(tasks)), initializer=initialize_worker
    ) as executor:
        futures = [executor.submit(tune_one, task) for task in tasks]
        for future in as_completed(futures):
            index, row = future.result()
            rows[index] = row
            write_checkpoint(args.output, args.group, rows)
            print(
                f"witness={sum(item is not None for item in rows)}/{len(rows)} "
                f"name={row['name']} combined={row['combined_log2']:.9f} "
                f"margin={row['margin_bits']:.9f}",
                flush=True,
            )
    write_checkpoint(args.output, args.group, rows)
    print("status=DIAGNOSTIC_BINARY64_PARALLEL_FIXED_PROFILE_ATLAS")


if __name__ == "__main__":
    main()
