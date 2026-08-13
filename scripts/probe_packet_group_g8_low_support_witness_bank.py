#!/usr/bin/env python3
"""Probe one representative profile in every low-support g=8 stratum.

This is a binary64 reconnaissance tool, not a certificate.  It constructs one
deterministic representative for each feasible support of size at most three,
tunes one explicitly selected bounded seed bank, freezes the resulting inner
witnesses, and tests their affine extensions with conditioned-row outer
witnesses.  The paired-local mode tunes one outer witness at each selected
seed.  Passing a representative does not establish coverage of its support
stratum.
"""

from __future__ import annotations

import argparse
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

from packet_group_drive_stratified import feasible_profile_count, profile_classes
from packet_group_native import load_shared_drive_apply
from packet_group_outer_profile import D, N, atom_count, normalization_log2
from packet_group_profile_bound import inner_probability, split_cap_table, tune_inner_density
from probe_packet_group_conditioned_row_outer import (
    evaluate_conditioned_outer,
    load_split_spectrum,
)
from survey_packet_group_g8_full_support_witness_reuse import optimize_outer_seed


GROUP_BITS = 8
CLASSES = GROUP_BITS + 1
ATOMS = atom_count(GROUP_BITS)
MINIMUM_PHYSICAL_WEIGHT = 21
EXPECTED_STRATA_BY_DIMENSION = {0: 8, 1: 36, 2: 84}
BASELINE_SEED_SUPPORTS = ((1,), (0, 1), (0, 4, 8))
HIGH_CLASS_SEED_SUPPORTS = (
    (7,),
    (8,),
    (0, 7),
    (0, 8),
    (7, 8),
    (0, 7, 8),
)
SEED_BANKS = {
    "baseline": BASELINE_SEED_SUPPORTS,
    "enriched-high-class": BASELINE_SEED_SUPPORTS + HIGH_CLASS_SEED_SUPPORTS,
}
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTER_BANK = ROOT / "scripts" / "packet_group_g8_conditioned_row_regression.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def allocate_positive(weights: list[int], total: int) -> list[int]:
    """Allocate ``total`` proportionally, with positive stable coordinates."""

    if not weights or total < len(weights) or any(weight <= 0 for weight in weights):
        raise ValueError("positive allocation has invalid weights or mass")
    result = [1] * len(weights)
    residual = total - len(weights)
    mass = sum(weights)
    raw = [residual * weight / mass for weight in weights]
    base = [math.floor(value) for value in raw]
    for index, value in enumerate(base):
        result[index] += value
    missing = total - sum(result)
    order = sorted(range(len(weights)), key=lambda index: (-(raw[index] - base[index]), index))
    for index in order[:missing]:
        result[index] += 1
    return result


def representative_profile(support: tuple[int, ...]) -> list[int]:
    """Return the concrete-configuration barycenter of one exact support."""

    allocated = allocate_positive(
        [math.comb(GROUP_BITS, weight) for weight in support], ATOMS
    )
    profile = [0] * CLASSES
    for weight, count in zip(support, allocated):
        profile[weight] = count
    return profile


def low_support_catalogue() -> list[dict]:
    rows = []
    ordinal = 0
    for size in (1, 2, 3):
        for support in itertools.combinations(range(CLASSES), size):
            profile = representative_profile(support)
            physical_weight = sum(index * count for index, count in enumerate(profile))
            if physical_weight < MINIMUM_PHYSICAL_WEIGHT:
                continue
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
    counts = {
        dimension: sum(row["dimension"] == dimension for row in rows)
        for dimension in EXPECTED_STRATA_BY_DIMENSION
    }
    if counts != EXPECTED_STRATA_BY_DIMENSION or len(rows) != 128:
        raise RuntimeError(f"low-support census changed: {counts}")
    return rows


def load_outer_bank(path: Path, spectra: tuple[np.ndarray, np.ndarray, np.ndarray]) -> list[dict]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("schema") != "permute-conv.packet-group-g8-conditioned-row-regression.v1":
        raise ValueError("unexpected frozen conditioned-row bank schema")
    spectrum01, punctured01, spectrum12 = spectra
    rows = []
    for source in artifact["rows"]:
        profile = [int(value) for value in source["profile"]]
        variables = np.asarray(source["log_variables"], dtype=np.float64)
        value, _details = evaluate_conditioned_outer(
            profile,
            variables,
            float(source["band1_coefficient"]),
            float(source["pair_cauchy_theta"]),
            spectrum01,
            punctured01,
            spectrum12,
            group_bits=GROUP_BITS,
        )
        if not math.isclose(value, float(source["outer_log2"]), rel_tol=0.0, abs_tol=1e-9):
            raise RuntimeError(f"frozen outer witness changed: {source['name']}")
        charge = variables / math.log(2.0)
        rows.append(
            {
                "name": str(source["name"]),
                "anchor_profile": profile,
                "anchor_outer_log2": value,
                "constant_log2": value + float(np.asarray(profile) @ charge),
                "charge": charge.tolist(),
                "log_variables": variables.tolist(),
                "band1_coefficient": float(source["band1_coefficient"]),
                "pair_cauchy_theta": float(source["pair_cauchy_theta"]),
            }
        )
    return rows


def tune_outer_seed(
    seed_index: int,
    row: dict,
    spectra: tuple[np.ndarray, np.ndarray, np.ndarray],
    mode: str,
    maximum_iterations: int,
    maximum_evaluations: int,
) -> dict:
    """Tune and freeze one support-local conditioned-row outer witness."""

    started = time.perf_counter()
    profile = row["profile"]
    result = optimize_outer_seed(
        profile,
        spectra,
        mode,
        maximum_iterations,
        maximum_evaluations,
    )
    variables = np.asarray(result["log_variables"], dtype=np.float64)
    charge = variables / math.log(2.0)
    return {
        "name": f"outer_seed_{seed_index}_support_" + "_".join(map(str, row["support"])),
        "source_ordinal": row["ordinal"],
        "source_support": row["support"],
        "anchor_profile": profile,
        "anchor_outer_log2": float(result["outer_log2"]),
        "constant_log2": float(result["outer_log2"])
        + float(np.asarray(profile, dtype=np.float64) @ charge),
        "charge": charge.tolist(),
        "log_variables": variables.tolist(),
        "band_coefficients": result["band_coefficients"],
        "pair_cauchy_theta": float(result["pair_cauchy_theta"]),
        "optimizer_success": bool(result["optimizer_success"]),
        "optimizer_iterations": int(result["optimizer_iterations"]),
        "optimizer_evaluations": int(result["optimizer_evaluations"]),
        "optimizer_message": str(result["optimizer_message"]),
        "elapsed_seconds": time.perf_counter() - started,
    }


def tune_inner_seed(
    seed_index: int,
    row: dict,
    split_caps,
    screen_iterations: int,
    final_iterations: int,
) -> dict:
    started = time.perf_counter()
    profile = row["profile"]
    screened = tune_inner_density(
        GROUP_BITS,
        profile,
        split_caps,
        poles=(0.2, 0.5, 0.8),
        exponents=(0.25, 0.5, 0.75),
        iterations=screen_iterations,
    )
    pole = float(screened[1])
    fugacities = np.maximum(np.asarray(screened[3], dtype=np.float64), 2.0**-80)
    fugacities /= float(np.max(fugacities))
    probability, details = inner_probability(
        GROUP_BITS,
        profile,
        pole,
        fugacities,
        split_caps,
        iterations=final_iterations,
    )
    constant = float(details["inner_mgf_log2"]) - D * math.log2(pole)
    reconstructed = (
        constant
        - float(np.asarray(profile) @ np.log2(fugacities))
        - float(normalization_log2(GROUP_BITS, np.asarray(profile))[0])
    )
    if not math.isclose(probability, reconstructed, rel_tol=0.0, abs_tol=1e-8):
        raise RuntimeError("inner affine reconstruction failed")
    return {
        "name": f"inner_seed_{seed_index}_support_" + "_".join(map(str, row["support"])),
        "source_ordinal": row["ordinal"],
        "source_support": row["support"],
        "source_profile": profile,
        "source_inner_log2": float(probability),
        "pole": pole,
        "fugacities": fugacities.tolist(),
        "charge": np.log2(fugacities).tolist(),
        "constant_log2": constant,
        "screen_evaluations": 9,
        "screen_iterations": screen_iterations,
        "final_iterations": final_iterations,
        "collatz_best_iteration": int(details["collatz_best_iteration"]),
        "elapsed_seconds": time.perf_counter() - started,
    }


def evaluate_catalogue(
    catalogue: list[dict], inner_bank: list[dict], outer_bank: list[dict], target: float
) -> tuple[list[dict], list[dict]]:
    profiles = np.asarray([row["profile"] for row in catalogue], dtype=np.float64)
    normalizations = normalization_log2(GROUP_BITS, profiles)
    columns = []
    witnesses = []
    for inner in inner_bank:
        for outer in outer_bank:
            charge = np.asarray(inner["charge"]) + np.asarray(outer["charge"])
            constant = float(inner["constant_log2"]) + float(outer["constant_log2"])
            columns.append(constant - profiles @ charge - normalizations)
            witnesses.append(
                {
                    "name": f"{inner['name']}__outer_{outer['name']}",
                    "inner": inner["name"],
                    "outer": outer["name"],
                    "constant_log2": constant,
                    "charge": charge.tolist(),
                }
            )
    matrix = np.asarray(columns, dtype=np.float64)
    leaders = np.argmin(matrix, axis=0)
    values = matrix[leaders, np.arange(len(catalogue))]
    rows = []
    for index, source in enumerate(catalogue):
        value = float(values[index])
        rows.append(
            {
                **source,
                "best_witness": witnesses[int(leaders[index])]["name"],
                "combined_log2": value,
                "uniform_profile_margin_bits": target - value,
                "representative_passes_uniform_target": value <= target,
            }
        )
    return rows, witnesses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seed-bank",
        choices=tuple(SEED_BANKS),
        required=True,
        help=(
            "baseline tunes the original three supports; enriched-high-class "
            "retains them and appends six deterministic supports near classes 7 and 8"
        ),
    )
    parser.add_argument("--screen-iterations", type=int, default=1)
    parser.add_argument("--final-iterations", type=int, default=4)
    parser.add_argument(
        "--outer-bank-mode",
        choices=("frozen-regression", "paired-local"),
        default="paired-local",
        help="use the five regression witnesses or tune one outer witness per selected seed",
    )
    parser.add_argument("--outer-bank", type=Path, default=DEFAULT_OUTER_BANK)
    parser.add_argument(
        "--conditioned-row-mode", choices=("symmetric", "asymmetric"), default="asymmetric"
    )
    parser.add_argument("--outer-max-iterations", type=int, default=80)
    parser.add_argument("--outer-max-evaluations", type=int, default=1200)
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
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.screen_iterations < 1 or args.final_iterations < 1:
        parser.error("iteration counts must be positive")

    started = time.perf_counter()
    catalogue = low_support_catalogue()
    by_support = {tuple(row["support"]): row for row in catalogue}
    seed_supports = SEED_BANKS[args.seed_bank]
    if seed_supports[: len(BASELINE_SEED_SUPPORTS)] != BASELINE_SEED_SUPPORTS:
        raise RuntimeError("selected seed bank does not retain the baseline prefix")
    seed_rows = [by_support[support] for support in seed_supports]
    spectra = (
        load_split_spectrum(args.spectrum01),
        load_split_spectrum(args.punctured01, count_field="pair_count", divisor=42),
        load_split_spectrum(args.spectrum12),
    )
    native = load_shared_drive_apply()
    native_path = getattr(native, "_packet_group_library_path", None) if native else None
    split_caps = split_cap_table()
    inner_bank = []
    for seed_index, row in enumerate(seed_rows):
        witness = tune_inner_seed(
            seed_index, row, split_caps, args.screen_iterations, args.final_iterations
        )
        inner_bank.append(witness)
        print(
            f"seed={seed_index + 1}/{len(seed_rows)} support={row['support']} "
            f"inner={witness['source_inner_log2']:.9f} "
            f"elapsed={witness['elapsed_seconds']:.3f}s",
            flush=True,
        )

    if args.outer_bank_mode == "frozen-regression":
        outer_bank = load_outer_bank(args.outer_bank, spectra)
    else:
        outer_bank = []
        for seed_index, row in enumerate(seed_rows):
            witness = tune_outer_seed(
                seed_index,
                row,
                spectra,
                args.conditioned_row_mode,
                args.outer_max_iterations,
                args.outer_max_evaluations,
            )
            outer_bank.append(witness)
            print(
                f"outer_seed={seed_index + 1}/{len(seed_rows)} support={row['support']} "
                f"outer={witness['anchor_outer_log2']:.9f} "
                f"elapsed={witness['elapsed_seconds']:.3f}s",
                flush=True,
            )

    feasible_profiles = feasible_profile_count(
        GROUP_BITS, N, MINIMUM_PHYSICAL_WEIGHT
    )
    target = -40.0 - math.log2(feasible_profiles)
    rows, combined_bank = evaluate_catalogue(catalogue, inner_bank, outer_bank, target)
    by_dimension = []
    for dimension in sorted(EXPECTED_STRATA_BY_DIMENSION):
        local = [row for row in rows if row["dimension"] == dimension]
        worst = max(local, key=lambda row: row["combined_log2"])
        by_dimension.append(
            {
                "dimension": dimension,
                "representatives_tested": len(local),
                "representatives_passing_uniform_target": sum(
                    row["representative_passes_uniform_target"] for row in local
                ),
                "worst_representative_support": worst["support"],
                "worst_representative_profile": worst["profile"],
                "worst_representative_combined_log2": worst["combined_log2"],
                "worst_representative_margin_bits": worst["uniform_profile_margin_bits"],
                "worst_representative_witness": worst["best_witness"],
            }
        )
    worst_rows = sorted(rows, key=lambda row: row["combined_log2"], reverse=True)
    report = {
        "schema": "permute-conv.packet-group-g8-low-support-recon.v3",
        "status": "DIAGNOSTIC_BINARY64_G8_LOW_SUPPORT_REPRESENTATIVES",
        "scope": (
            "one deterministic representative in each feasible support stratum of "
            "dimension 0, 1, or 2; no stratum coverage or completeness claim"
        ),
        "group_bits": GROUP_BITS,
        "atoms": ATOMS,
        "minimum_physical_weight": MINIMUM_PHYSICAL_WEIGHT,
        "feasible_profile_count": feasible_profiles,
        "uniform_profile_target_log2": target,
        "experiment_question": (
            "Does tuning both inner and conditioned-row outer witnesses at nine "
            "support-local seeds repair affine transfer on the 128 representatives?"
        ),
        "configuration": {
            "seed_bank": args.seed_bank,
            "seed_count": len(seed_supports),
            "baseline_seed_supports": [
                list(value) for value in BASELINE_SEED_SUPPORTS
            ],
            "high_class_seed_supports": [
                list(value) for value in HIGH_CLASS_SEED_SUPPORTS
            ],
            "selected_seed_supports": [list(value) for value in seed_supports],
            "screen_iterations": args.screen_iterations,
            "final_iterations": args.final_iterations,
            "outer_bank_mode": args.outer_bank_mode,
            "conditioned_row_mode": args.conditioned_row_mode,
            "outer_max_iterations": args.outer_max_iterations,
            "outer_max_evaluations": args.outer_max_evaluations,
            "serialized_compute": True,
            "thread_environment": {
                variable: os.environ.get(variable)
                for variable in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "BLIS_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            },
            "argv": sys.argv,
        },
        "sources": {
            "outer_bank": (
                {"path": str(args.outer_bank), "sha256": digest(args.outer_bank)}
                if args.outer_bank_mode == "frozen-regression"
                else None
            ),
            "spectrum01": {"path": str(args.spectrum01), "sha256": digest(args.spectrum01)},
            "punctured01": {"path": str(args.punctured01), "sha256": digest(args.punctured01)},
            "spectrum12": {"path": str(args.spectrum12), "sha256": digest(args.spectrum12)},
            "conditioned_row_evaluator": {
                "path": str(ROOT / "scripts" / "probe_packet_group_conditioned_row_outer.py"),
                "sha256": digest(ROOT / "scripts" / "probe_packet_group_conditioned_row_outer.py"),
            },
            "native_inner_library": native_path,
        },
        "representative_rule": (
            "positive largest-remainder allocation proportional to binom(8,w) "
            "within the exact support, with stable class-index ties"
        ),
        "representatives_tested": len(rows),
        "representatives_passing_uniform_target": sum(
            row["representative_passes_uniform_target"] for row in rows
        ),
        "strata_by_dimension": by_dimension,
        "inner_witness_bank": inner_bank,
        "outer_witness_bank": outer_bank,
        "combined_affine_witness_count": len(combined_bank),
        "combined_affine_witness_bank": combined_bank,
        "worst_representatives": worst_rows[:10],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
        "blockers": [
            "binary64 inner witnesses lack independent outward replay",
            "one representative cannot bound the other integer profiles in its stratum",
            "no exact slab or BSP ownership ledger exists for these 128 strata",
            "the selected tuned inner/outer bank is a transfer experiment, not a complete witness basis",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        f"representatives={len(rows)} passing={report['representatives_passing_uniform_target']} "
        f"worst={worst_rows[0]['combined_log2']:.9f} "
        f"margin={worst_rows[0]['uniform_profile_margin_bits']:.9f}",
        flush=True,
    )
    print(f"output={args.output}", flush=True)
    print("status=DIAGNOSTIC_BINARY64_G8_LOW_SUPPORT_REPRESENTATIVES", flush=True)


if __name__ == "__main__":
    main()
