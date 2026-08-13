#!/usr/bin/env python3
"""Build a diagnostic support-local witness atlas for packet width g=8.

The job tunes one inner witness and one conditioned-row outer witness at a
deterministic representative of each feasible exact support.  Independent
workers process supports in parallel.  Each completed row is checkpointed
atomically.  The result seeds later exact shard refinement; it is not a
coverage certificate.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import itertools
import json
import math
import os
import sys
import time
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

from packet_group_drive_stratified import feasible_profile_count
from packet_group_native import load_shared_drive_apply
from packet_group_outer_profile import D, N, normalization_log2
from packet_group_profile_bound import split_cap_table
from probe_packet_group_conditioned_row_outer import load_split_spectrum
from probe_packet_group_g8_low_support_witness_bank import (
    GROUP_BITS,
    MINIMUM_PHYSICAL_WEIGHT,
    representative_profile,
    tune_inner_seed,
    tune_outer_seed,
)


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_SUPPORTS = 510
_SPLIT_CAPS = None
_SPECTRA = None


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def support_catalogue() -> list[dict]:
    """Return one deterministic representative for every feasible support."""

    rows = []
    ordinal = 0
    for size in range(1, GROUP_BITS + 2):
        for support in itertools.combinations(range(GROUP_BITS + 1), size):
            if support == (0,):
                continue
            profile = representative_profile(support)
            physical_weight = sum(index * count for index, count in enumerate(profile))
            if physical_weight < MINIMUM_PHYSICAL_WEIGHT:
                raise AssertionError(f"representative violates distance cutoff: {support}")
            rows.append(
                {
                    "ordinal": ordinal,
                    "support": list(support),
                    "support_mask": sum(1 << index for index in support),
                    "dimension": size - 1,
                    "profile": profile,
                    "physical_weight": physical_weight,
                }
            )
            ordinal += 1
    if len(rows) != EXPECTED_SUPPORTS:
        raise AssertionError(f"support census changed: {len(rows)}")
    return rows


def initialize_worker(spectrum01: str, punctured01: str, spectrum12: str) -> None:
    global _SPLIT_CAPS, _SPECTRA
    _SPLIT_CAPS = split_cap_table()
    _SPECTRA = (
        load_split_spectrum(Path(spectrum01)),
        load_split_spectrum(Path(punctured01), count_field="pair_count", divisor=42),
        load_split_spectrum(Path(spectrum12)),
    )
    native = load_shared_drive_apply()
    if native is None:
        raise RuntimeError("support-atlas workers require the native inner kernel")


def build_witness(task: tuple[dict, dict]) -> dict:
    row, configuration = task
    if _SPLIT_CAPS is None or _SPECTRA is None:
        raise RuntimeError("worker was not initialized")
    started = time.perf_counter()
    inner = tune_inner_seed(
        row["ordinal"],
        row,
        _SPLIT_CAPS,
        configuration["inner_screen_iterations"],
        configuration["inner_final_iterations"],
    )
    outer = tune_outer_seed(
        row["ordinal"],
        row,
        _SPECTRA,
        configuration["conditioned_row_mode"],
        configuration["outer_max_iterations"],
        configuration["outer_max_evaluations"],
    )
    profile = np.asarray(row["profile"], dtype=np.float64)
    combined_constant = float(inner["constant_log2"]) + float(outer["constant_log2"])
    combined_charge = np.asarray(inner["charge"]) + np.asarray(outer["charge"])
    combined = combined_constant - float(profile @ combined_charge)
    combined -= float(normalization_log2(GROUP_BITS, profile)[0])
    local = float(inner["source_inner_log2"]) + float(outer["anchor_outer_log2"])
    if not math.isclose(combined, local, rel_tol=0.0, abs_tol=2e-8):
        raise AssertionError("combined affine witness failed local replay")
    return {
        **row,
        "local_combined_log2": local,
        "inner": inner,
        "outer": outer,
        "affine": {
            "constant_log2": combined_constant,
            "charge_log2": combined_charge.tolist(),
            "normalization": "multinomial packet-profile normalization",
        },
        "elapsed_seconds": time.perf_counter() - started,
    }


def checkpoint_path(directory: Path, support_mask: int) -> Path:
    return directory / f"support_{support_mask:03x}.json"


def write_json_atomic(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def load_completed(directory: Path, configuration_digest: str) -> dict[int, dict]:
    completed = {}
    if not directory.exists():
        return completed
    for path in sorted(directory.glob("support_*.json")):
        artifact = json.loads(path.read_text(encoding="utf-8"))
        if artifact.get("configuration_digest") != configuration_digest:
            raise ValueError(f"checkpoint configuration mismatch: {path}")
        row = artifact["row"]
        mask = int(row["support_mask"])
        if mask in completed:
            raise ValueError(f"duplicate support checkpoint: {mask}")
        completed[mask] = row
    return completed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--inner-screen-iterations", type=int, default=1)
    parser.add_argument("--inner-final-iterations", type=int, default=3)
    parser.add_argument(
        "--conditioned-row-mode", choices=("symmetric", "asymmetric"), default="asymmetric"
    )
    parser.add_argument("--outer-max-iterations", type=int, default=60)
    parser.add_argument("--outer-max-evaluations", type=int, default=800)
    parser.add_argument(
        "--spectrum01", type=Path, default=ROOT / "out" / "ebch85_band01_split_spectrum.csv"
    )
    parser.add_argument(
        "--punctured01",
        type=Path,
        default=ROOT / "out" / "ebch84_punctured_band01_split_spectrum.csv",
    )
    parser.add_argument(
        "--spectrum12", type=Path, default=ROOT / "out" / "ebch86_band12_split_spectrum.csv"
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")

    configuration = {
        "group_bits": GROUP_BITS,
        "inner_screen_iterations": args.inner_screen_iterations,
        "inner_final_iterations": args.inner_final_iterations,
        "conditioned_row_mode": args.conditioned_row_mode,
        "outer_max_iterations": args.outer_max_iterations,
        "outer_max_evaluations": args.outer_max_evaluations,
        "spectrum01_sha256": sha256(args.spectrum01),
        "punctured01_sha256": sha256(args.punctured01),
        "spectrum12_sha256": sha256(args.spectrum12),
    }
    encoded_configuration = json.dumps(
        configuration, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    configuration_digest = hashlib.sha256(encoded_configuration).hexdigest()
    catalogue = support_catalogue()
    completed = load_completed(args.checkpoint_dir, configuration_digest)
    pending = [row for row in catalogue if row["support_mask"] not in completed]
    print(
        f"supports={len(catalogue)} resumed={len(completed)} pending={len(pending)} "
        f"workers={args.workers}",
        flush=True,
    )
    started = time.perf_counter()
    initializer_args = (
        str(args.spectrum01.resolve()),
        str(args.punctured01.resolve()),
        str(args.spectrum12.resolve()),
    )
    tasks = [(row, configuration) for row in pending]
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=initialize_worker,
        initargs=initializer_args,
    ) as executor:
        futures = {executor.submit(build_witness, task): task[0] for task in tasks}
        for finished, future in enumerate(concurrent.futures.as_completed(futures), 1):
            source = futures[future]
            row = future.result()
            completed[row["support_mask"]] = row
            write_json_atomic(
                checkpoint_path(args.checkpoint_dir, row["support_mask"]),
                {"configuration_digest": configuration_digest, "row": row},
            )
            if finished == 1 or finished % 10 == 0 or finished == len(pending):
                print(
                    f"completed={len(completed)}/{len(catalogue)} "
                    f"last_support={source['support']} local={row['local_combined_log2']:.6f}",
                    flush=True,
                )

    rows = [completed[row["support_mask"]] for row in catalogue]
    target = -40.0 - math.log2(
        feasible_profile_count(GROUP_BITS, N, MINIMUM_PHYSICAL_WEIGHT)
    )
    local_margins = [target - float(row["local_combined_log2"]) for row in rows]
    report = {
        "schema": "permute-conv.packet-group-g8-support-seed-atlas.v1",
        "status": "DIAGNOSTIC_BINARY64_SUPPORT_LOCAL_ATLAS",
        "scope": "one locally tuned witness at one representative of each exact support",
        "configuration": configuration,
        "configuration_digest": configuration_digest,
        "workers": args.workers,
        "support_count": len(rows),
        "uniform_profile_target_log2": target,
        "local_representatives_passing": sum(margin >= 0.0 for margin in local_margins),
        "worst_local_margin_bits": min(local_margins),
        "worst_local_support": rows[int(np.argmin(local_margins))]["support"],
        "rows": rows,
        "elapsed_seconds_this_invocation": time.perf_counter() - started,
        "argv": sys.argv,
        "blockers": [
            "binary64 witnesses require outward hardening before certification",
            "one representative does not cover its exact support stratum",
            "no exact slab or BSP ownership ledger is attached",
            "the atlas does not perform the global union aggregation",
        ],
    }
    write_json_atomic(args.output, report)
    print(
        f"passing={report['local_representatives_passing']}/{len(rows)} "
        f"worst_margin={report['worst_local_margin_bits']:.6f}",
        flush=True,
    )
    print(f"output={args.output}", flush=True)


if __name__ == "__main__":
    main()
