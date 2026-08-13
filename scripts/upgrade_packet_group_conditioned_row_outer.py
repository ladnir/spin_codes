#!/usr/bin/env python3
"""Attach optimized one-conditioned-row outers to frozen packet-group witnesses."""

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

from probe_packet_group_conditioned_row_outer import (
    load_split_spectrum,
    optimize_conditioned_outer,
    optimize_conditioned_outer_asymmetric,
)


_SPECTRA = None


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact_rows(path: Path) -> list[dict]:
    value = json.loads(path.read_text(encoding="utf-8"))
    rows = value if isinstance(value, list) else value.get("rows", [])
    if not isinstance(rows, list):
        raise ValueError(f"conditioned-row upgrader: no rows in {path}")
    return rows


def initialize_worker(spectrum01: str, punctured01: str, spectrum12: str) -> None:
    global _SPECTRA
    _SPECTRA = (
        load_split_spectrum(Path(spectrum01)),
        load_split_spectrum(Path(punctured01), count_field="pair_count", divisor=42),
        load_split_spectrum(Path(spectrum12)),
    )


def upgrade(task: tuple[int, dict, bool]) -> tuple[int, dict]:
    index, row, asymmetric = task
    if _SPECTRA is None:
        raise RuntimeError("conditioned-row upgrader worker is not initialized")
    group_bits = int(row.get("group_bits", -1))
    if group_bits not in (1, 2, 4, 8, 16, 32, 64) or "fugacities" not in row:
        raise ValueError("conditioned-row upgrader requires a supported inner witness")
    profile = [int(value) for value in row["profile"]]
    details = row.get("outer_details", {})
    initial_logs = details.get("log_variables")
    if initial_logs is None:
        initial_logs = [0.0] * (group_bits + 1)
    initial_logs = np.asarray(initial_logs, dtype=np.float64)
    if initial_logs.shape != (group_bits + 1,):
        raise ValueError("conditioned-row upgrader: outer variable dimension mismatch")
    initial_band1 = float(details.get("band1_coefficient", 0.6))
    if asymmetric:
        initial_coefficients = details.get("band_coefficients")
        initial_band2 = (
            float(initial_coefficients[2])
            if isinstance(initial_coefficients, list) and len(initial_coefficients) == 3
            else initial_band1
        )
        value, variables, coefficients, theta, result, conditioned = (
            optimize_conditioned_outer_asymmetric(
                profile,
                initial_logs,
                initial_band1,
                initial_band2,
                0.5,
                *_SPECTRA,
                group_bits,
            )
        )
        band1 = coefficients[1]
        outer_type = "conditioned_row_exact_graph_asymmetric"
    else:
        value, variables, band1, theta, result, conditioned = optimize_conditioned_outer(
            profile,
            initial_logs,
            initial_band1,
            0.5,
            *_SPECTRA,
            1,
            group_bits,
        )
        coefficients = (1.0 - band1, band1, band1)
        outer_type = "conditioned_row_exact_graph"
    outer_charge = variables / math.log(2.0)
    profile_vector = np.asarray(profile, dtype=np.float64)
    outer_constant = value + float(profile_vector @ outer_charge)
    fugacities = np.asarray(row["fugacities"], dtype=np.float64)
    if fugacities.shape != (group_bits + 1,) or np.any(fugacities < 0.0):
        raise ValueError("conditioned-row upgrader: malformed fugacity vector")
    if any(count and fugacity == 0.0 for count, fugacity in zip(profile, fugacities)):
        raise ValueError("conditioned-row upgrader: occupied class has zero fugacity")
    with np.errstate(divide="ignore"):
        inner_charge = np.log2(fugacities)
    charge = outer_charge + inner_charge
    constant = outer_constant + float(row["inner_constant_log2"])
    combined = (
        constant
        - sum(count * float(charge[index]) for index, count in enumerate(profile) if count)
        - float(row["normalization_log2"])
    )
    target = float(row["target_log2"])
    return index, {
        **row,
        "name": str(row.get("name", "witness")) + "__conditioned_row",
        "outer_type": outer_type,
        "outer_details": {
            "log_variables": variables.tolist(),
            "band1_coefficient": band1,
            "band_coefficients": list(coefficients),
            "pair_cauchy_theta": theta,
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_iterations": int(result.nit),
            "optimizer_evaluations": int(result.nfev),
            **conditioned,
        },
        "outer_log2": value,
        "outer_constant_log2": outer_constant,
        "outer_charge": outer_charge.tolist(),
        "constant_log2": constant,
        "charge": charge.tolist(),
        "combined_log2": combined,
        "margin_bits": target - combined,
        "status": "DIAGNOSTIC_BINARY64_PACKET_GROUP_CONDITIONED_ROW_UPGRADED_WITNESS",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument("--artifact", type=Path, action="append", required=True)
    parser.add_argument("--name", action="append", default=[])
    parser.add_argument(
        "--reference",
        action="append",
        default=[],
        help="select BASENAME.json:INDEX from the supplied artifacts",
    )
    parser.add_argument(
        "--residual-ledger",
        type=Path,
        help="also select every diagnostic_best_witness named by this ledger",
    )
    parser.add_argument("--max-rows", type=int, default=0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--asymmetric",
        action="store_true",
        help="optimize over the valid p0+p1=1, p2>=p1 BL boundary",
    )
    parser.add_argument(
        "--spectrum01", type=Path, default=root / "out" / "ebch85_band01_split_spectrum.csv"
    )
    parser.add_argument(
        "--punctured01",
        type=Path,
        default=root / "out" / "ebch84_punctured_band01_split_spectrum.csv",
    )
    parser.add_argument(
        "--spectrum12", type=Path, default=root / "out" / "ebch86_band12_split_spectrum.csv"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 8 or args.max_rows < 0:
        parser.error("workers must be in [1,8] and max-rows must be nonnegative")

    indexed = [
        (f"{artifact.name}:{index}", row)
        for artifact in args.artifact
        for index, row in enumerate(artifact_rows(artifact))
    ]
    requested_references = set(args.reference)
    if args.residual_ledger is not None:
        residual = json.loads(args.residual_ledger.read_text(encoding="utf-8"))
        requested_references.update(
            str(row["diagnostic_best_witness"])
            for row in residual.get("uncovered_integer_residuals", [])
        )
    if requested_references:
        indexed = [item for item in indexed if item[0] in requested_references]
        if {reference for reference, _row in indexed} != requested_references:
            raise ValueError("conditioned-row upgrader: requested reference is missing")
    rows = [row for _reference, row in indexed]
    if args.name:
        requested = set(args.name)
        rows = [row for row in rows if row.get("name") in requested]
        if {row.get("name") for row in rows} != requested:
            raise ValueError("conditioned-row upgrader: requested witness is missing")
    rows = rows[: args.max_rows or None]
    if not rows:
        raise ValueError("conditioned-row upgrader: no rows selected")
    group_widths = {int(row.get("group_bits", -1)) for row in rows}
    if len(group_widths) != 1:
        raise ValueError("conditioned-row upgrader: selected rows mix packet widths")
    group_bits = next(iter(group_widths))
    tasks = [(index, row, args.asymmetric) for index, row in enumerate(rows)]
    completed: list[dict | None] = [None] * len(tasks)
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=min(args.workers, len(tasks)),
        initializer=initialize_worker,
        initargs=(str(args.spectrum01), str(args.punctured01), str(args.spectrum12)),
    ) as pool:
        for future in concurrent.futures.as_completed(
            pool.submit(upgrade, task) for task in tasks
        ):
            index, row = future.result()
            completed[index] = row
            print(
                f"profile={row['profile']} combined={row['combined_log2']:.9f} "
                f"margin={row['margin_bits']:.9f}",
                flush=True,
            )
    report = {
        "status": "DIAGNOSTIC_BINARY64_PACKET_GROUP_CONDITIONED_ROW_ATLAS",
        "group_bits": group_bits,
        "asymmetric": args.asymmetric,
        "complete": all(row is not None for row in completed),
        "source_artifacts": [
            {"path": str(path), "sha256": digest(path)} for path in args.artifact
        ],
        "sources": {
            str(args.spectrum01): digest(args.spectrum01),
            str(args.punctured01): digest(args.punctured01),
            str(args.spectrum12): digest(args.spectrum12),
        },
        "rows": [row for row in completed if row is not None],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
