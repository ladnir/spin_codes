#!/usr/bin/env python3
"""Run a finite, reproducible g=8 profile reconnaissance.

This is a binary64 discovery tool, not a domain certificate.  It constructs a
fixed list of profiles, tunes one inexpensive shared-drive inner witness per
profile, and compares several valid profile-local outer branches.  The
conditioned-row branch uses the generalized evaluator from
``probe_packet_group_conditioned_row_outer.py`` without modifying it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from packet_group_drive_stratified import feasible_profile_count, profile_count
from packet_group_outer_profile import D, N, atom_count, fixed_total_weight_outer
from packet_group_profile_bound import (
    inner_probability,
    outer_probability,
    refine_inner_sparse_fugacities,
    split_cap_table,
    tune_inner_density,
)
from probe_packet_group_conditioned_row_outer import (
    evaluate_conditioned_outer,
    load_split_spectrum,
    optimize_conditioned_outer,
)


GROUP_BITS = 8
CLASSES = GROUP_BITS + 1
ATOMS = atom_count(GROUP_BITS)
MINIMUM_NONZERO_OUTER_WEIGHT = 21
FEASIBLE_PROFILES = feasible_profile_count(
    GROUP_BITS, N, MINIMUM_NONZERO_OUTER_WEIGHT
)


def allocate_rational(weights: list[Fraction], total: int = ATOMS) -> list[int]:
    """Round a rational distribution by largest remainder and stable index."""

    mass = sum(weights)
    if mass <= 0:
        raise ValueError("profile distribution has zero mass")
    scaled = [weight * total / mass for weight in weights]
    result = [value.numerator // value.denominator for value in scaled]
    missing = total - sum(result)
    order = sorted(
        range(len(weights)),
        key=lambda index: (-(scaled[index] - result[index]), index),
    )
    for index in order[:missing]:
        result[index] += 1
    if sum(result) != total:
        raise AssertionError("largest-remainder allocation lost mass")
    return result


def pure_profile(weight: int) -> list[int]:
    profile = [0] * CLASSES
    profile[weight] = ATOMS
    return profile


def binomial_profile(numerator: int, denominator: int) -> list[int]:
    p = Fraction(numerator, denominator)
    weights = [
        Fraction(math.comb(GROUP_BITS, weight))
        * p**weight
        * (1 - p) ** (GROUP_BITS - weight)
        for weight in range(CLASSES)
    ]
    return allocate_rational(weights)


def profile_catalogue() -> list[dict]:
    """Return the ordered deterministic profile catalogue."""

    distance_edge = [ATOMS - D, D] + [0] * (CLASSES - 2)
    smoke = [
        {
            "name": "pure_1",
            "family": "pure_class",
            "construction": "all 2^18 packets have weight 1",
            "profile": pure_profile(1),
        },
        {
            "name": "distance_edge_0_1",
            "family": "two_class_support",
            "construction": "c1=D and c0=2^18-D; physical weight is exactly D",
            "profile": distance_edge,
        },
        {
            "name": "adjacent_3_4_half",
            "family": "two_class_support",
            "construction": "equal mass on packet weights 3 and 4",
            "profile": allocate_rational(
                [Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2),
                 Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0),
                 Fraction(0)]
            ),
        },
        {
            "name": "binomial_3_32",
            "family": "full_support_binomial",
            "construction": "stable largest-remainder rounding of Binomial(8,3/32)",
            "profile": binomial_profile(3, 32),
        },
        {
            "name": "binomial_1_2",
            "family": "full_support_binomial",
            "construction": "exact Binomial(8,1/2) class distribution",
            "profile": binomial_profile(1, 2),
        },
    ]
    extras = [
        *(
            {
                "name": f"pure_{weight}",
                "family": "pure_class",
                "construction": f"all 2^18 packets have weight {weight}",
                "profile": pure_profile(weight),
            }
            for weight in range(2, CLASSES)
        ),
        {
            "name": "symmetric_1_7_half",
            "family": "two_class_support",
            "construction": "equal mass on packet weights 1 and 7",
            "profile": allocate_rational(
                [Fraction(0), Fraction(1, 2), Fraction(0), Fraction(0),
                 Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2),
                 Fraction(0)]
            ),
        },
        {
            "name": "three_class_0_4_8",
            "family": "three_class_support",
            "construction": "mass 1/4,1/2,1/4 on packet weights 0,4,8",
            "profile": allocate_rational(
                [Fraction(1, 4), Fraction(0), Fraction(0), Fraction(0),
                 Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0),
                 Fraction(1, 4)]
            ),
        },
        {
            "name": "class_uniform",
            "family": "full_support_uniform",
            "construction": "stable largest-remainder allocation across nine classes",
            "profile": allocate_rational([Fraction(1, CLASSES)] * CLASSES),
        },
        *(
            {
                "name": f"binomial_{numerator}_{denominator}",
                "family": "full_support_binomial",
                "construction": (
                    "stable largest-remainder rounding of "
                    f"Binomial(8,{numerator}/{denominator})"
                ),
                "profile": binomial_profile(numerator, denominator),
            }
            for numerator, denominator in (
                (1, 16),
                (1, 8),
                (1, 4),
                (3, 4),
                (15, 16),
            )
        ),
    ]
    # Preserve the smoke ordering and remove duplicate names introduced by
    # the extended pure-class family.
    seen = set()
    result = []
    for row in smoke + extras:
        if row["name"] not in seen:
            seen.add(row["name"])
            result.append(row)
    return result


def select_profiles(suite: str, names: list[str], maximum: int) -> list[dict]:
    rows = profile_catalogue()
    if names:
        by_name = {row["name"]: row for row in rows}
        unknown = sorted(set(names) - set(by_name))
        if unknown:
            raise ValueError(f"unknown profile names: {', '.join(unknown)}")
        rows = [by_name[name] for name in names]
    elif suite == "smoke":
        rows = rows[:5]
    if maximum:
        rows = rows[:maximum]
    return rows


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tune_inner(
    profile: list[int], split_caps, screen_iterations: int, final_iterations: int
) -> tuple[float, float, np.ndarray, dict, str, int]:
    initial = tune_inner_density(
        GROUP_BITS,
        profile,
        split_caps,
        poles=(0.2, 0.5, 0.8),
        exponents=(0.25, 0.5, 0.75),
        iterations=screen_iterations,
    )
    pole = float(initial[1])
    fugacities = np.asarray(initial[3], dtype=np.float64)
    method = "density_grid"
    search_evaluations = 9
    if any(count == 0 for count in profile):
        sparse = refine_inner_sparse_fugacities(
            GROUP_BITS,
            profile,
            split_caps,
            initial_pole=pole,
            initial_fugacities=fugacities,
            log_steps=(math.log(4.0), math.log(2.0)),
            witness_iterations=screen_iterations,
            zero_floor=2.0**-16,
        )
        pole = float(sparse[1])
        fugacities = np.asarray(sparse[2], dtype=np.float64)
        fugacities[np.asarray(profile, dtype=np.int64) == 0] = 0.0
        method = "support_sparse_exact_face"
        search_evaluations += int(sparse[5])
    value, details = inner_probability(
        GROUP_BITS,
        profile,
        pole,
        fugacities,
        split_caps,
        iterations=final_iterations,
    )
    return value, pole, fugacities, details, method, search_evaluations


def evaluate_profile(
    descriptor: dict,
    split_caps,
    conditioned_mode: str,
    spectra: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
    screen_iterations: int,
    final_iterations: int,
) -> dict:
    started = time.perf_counter()
    profile = descriptor["profile"]
    inner_value, pole, fugacities, inner_details, inner_method, evaluations = (
        tune_inner(profile, split_caps, screen_iterations, final_iterations)
    )
    linear_value, linear_details = outer_probability(
        GROUP_BITS, profile, optimize_band_coefficients=True, fast=True
    )
    physical_weight = sum(index * count for index, count in enumerate(profile))
    total_value, total_details = fixed_total_weight_outer(physical_weight)
    branches = {
        "linear_bl_fast": {
            "outer_log2": float(linear_value),
            "band1_coefficient": float(linear_details["band1_coefficient"]),
            "log_variables": linear_details["log_variables"],
            "optimizer_success": bool(linear_details["optimizer_success"]),
        },
        "total_weight": {
            "outer_log2": float(total_value),
            "charge_per_bit": float(total_details["charge_per_bit"]),
        },
    }
    if conditioned_mode != "off":
        if spectra is None:
            raise AssertionError("conditioned spectra were not loaded")
        spectrum01, punctured01, spectrum12 = spectra
        initial_variables = np.asarray(
            linear_details["log_variables"], dtype=np.float64
        )
        initial_band1 = float(linear_details["band1_coefficient"])
        if conditioned_mode == "optimize":
            (
                conditioned_value,
                conditioned_variables,
                conditioned_band1,
                conditioned_theta,
                result,
                _details,
            ) = optimize_conditioned_outer(
                profile,
                initial_variables,
                initial_band1,
                0.5,
                spectrum01,
                punctured01,
                spectrum12,
                1,
                GROUP_BITS,
            )
            optimizer = {
                "success": bool(result.success),
                "iterations": int(result.nit),
                "evaluations": int(result.nfev),
                "message": str(result.message),
            }
        else:
            conditioned_variables = initial_variables
            conditioned_band1 = initial_band1
            conditioned_theta = 0.5
            conditioned_value, _details = evaluate_conditioned_outer(
                profile,
                conditioned_variables,
                conditioned_band1,
                conditioned_theta,
                spectrum01,
                punctured01,
                spectrum12,
                1,
                None,
                GROUP_BITS,
            )
            optimizer = None
        branches[f"conditioned_row_{conditioned_mode}"] = {
            "outer_log2": float(conditioned_value),
            "band1_coefficient": float(conditioned_band1),
            "pair_cauchy_theta": float(conditioned_theta),
            "log_variables": conditioned_variables.tolist(),
            "optimizer": optimizer,
        }
    best_name, best = min(
        branches.items(), key=lambda item: float(item[1]["outer_log2"])
    )
    combined = float(best["outer_log2"]) + float(inner_value)
    target = -40.0 - math.log2(FEASIBLE_PROFILES)
    compact_inner = {
        "probability_log2": float(inner_value),
        "method": inner_method,
        "pole": pole,
        "fugacities": fugacities.tolist(),
        "search_evaluations": evaluations,
        "collatz_best_iteration": int(inner_details["collatz_best_iteration"]),
        "collatz_max_iterations": int(inner_details["collatz_max_iterations"]),
        "worst_state": int(inner_details["worst_state"]),
        "inner_mgf_log2": float(inner_details["inner_mgf_log2"]),
        "normalization_log2": float(inner_details["normalization_log2"]),
    }
    return {
        **descriptor,
        "support": [index for index, count in enumerate(profile) if count],
        "physical_weight": physical_weight,
        "relative_weight": physical_weight / N,
        "inner": compact_inner,
        "outer_branches": branches,
        "selected_outer": best_name,
        "combined_log2": combined,
        "uniform_profile_target_log2": target,
        "uniform_profile_margin_bits": target - combined,
        "single_profile_40bit_margin_bits": -40.0 - combined,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument("--suite", choices=("smoke", "recon"), default="smoke")
    parser.add_argument("--profile", action="append", default=[])
    parser.add_argument("--max-profiles", type=int, default=0)
    parser.add_argument("--screen-iterations", type=int, default=1)
    parser.add_argument("--final-iterations", type=int, default=4)
    parser.add_argument(
        "--conditioned", choices=("off", "evaluate", "optimize"), default="evaluate"
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
    if args.max_profiles < 0 or args.screen_iterations <= 0 or args.final_iterations <= 0:
        raise SystemExit("g=8 bounded recon: invalid positive iteration/profile limit")
    try:
        selected = select_profiles(args.suite, args.profile, args.max_profiles)
    except ValueError as error:
        raise SystemExit(f"g=8 bounded recon: {error}") from error
    if not selected:
        raise SystemExit("g=8 bounded recon: no profiles selected")

    spectra = None
    sources = {}
    if args.conditioned != "off":
        spectra = (
            load_split_spectrum(args.spectrum01),
            load_split_spectrum(args.punctured01, count_field="pair_count", divisor=42),
            load_split_spectrum(args.spectrum12),
        )
        sources = {
            "spectrum01": {"path": str(args.spectrum01), "sha256": sha256(args.spectrum01)},
            "punctured01": {"path": str(args.punctured01), "sha256": sha256(args.punctured01)},
            "spectrum12": {"path": str(args.spectrum12), "sha256": sha256(args.spectrum12)},
        }
    evaluator_path = Path(__file__).with_name(
        "probe_packet_group_conditioned_row_outer.py"
    )
    sources["conditioned_row_evaluator"] = {
        "path": str(evaluator_path),
        "sha256": sha256(evaluator_path),
    }

    split_caps = split_cap_table()
    rows = []
    run_started = time.perf_counter()
    for index, descriptor in enumerate(selected, 1):
        row = evaluate_profile(
            descriptor,
            split_caps,
            args.conditioned,
            spectra,
            args.screen_iterations,
            args.final_iterations,
        )
        rows.append(row)
        print(
            f"profile={index}/{len(selected)} name={row['name']} "
            f"combined={row['combined_log2']:.9f} "
            f"uniform_margin={row['uniform_profile_margin_bits']:.9f} "
            f"outer={row['selected_outer']} elapsed={row['elapsed_seconds']:.3f}s",
            flush=True,
        )
    worst = sorted(rows, key=lambda row: row["combined_log2"], reverse=True)
    tested_union = float(
        logsumexp(
            np.asarray([row["combined_log2"] for row in rows]) * math.log(2.0)
        )
        / math.log(2.0)
    )
    report = {
        "status": "DIAGNOSTIC_BINARY64_G8_BOUNDED_PROFILE_RECON",
        "scope": (
            "finite deterministic profile set only; no interpolation, domain coverage, "
            "or completeness claim"
        ),
        "group_bits": GROUP_BITS,
        "atoms": ATOMS,
        "N": N,
        "D": D,
        "ambient_profile_count": profile_count(GROUP_BITS, N),
        "minimum_nonzero_outer_weight": MINIMUM_NONZERO_OUTER_WEIGHT,
        "feasible_profile_count": FEASIBLE_PROFILES,
        "uniform_profile_target_log2": -40.0 - math.log2(FEASIBLE_PROFILES),
        "configuration": {
            "suite": args.suite,
            "requested_profiles": args.profile,
            "max_profiles": args.max_profiles,
            "screen_iterations": args.screen_iterations,
            "final_iterations": args.final_iterations,
            "conditioned": args.conditioned,
            "argv": sys.argv,
        },
        "sources": sources,
        "profiles_evaluated": len(rows),
        "tested_set_union_log2": tested_union,
        "tested_set_40bit_margin_bits": -40.0 - tested_union,
        "worst_profiles": [
            {
                "name": row["name"],
                "profile": row["profile"],
                "combined_log2": row["combined_log2"],
                "uniform_profile_margin_bits": row["uniform_profile_margin_bits"],
                "selected_outer": row["selected_outer"],
            }
            for row in worst[: min(5, len(worst))]
        ],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - run_started,
        "blockers_to_end_to_end_survey": [
            "binary64 witnesses and optimized parameters lack independent outward replay",
            "the finite profile set gives no coverage of the eight-dimensional profile simplex",
            "no exact BSP, support-slab ownership, or complete profile-count aggregation exists for g=8",
            "conditioned-row generalization still needs frozen g=8 regression vectors and an outward implementation",
            "profile-local inner tuning is too expensive to repeat over a dense survey without witness reuse",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_BINARY64_G8_BOUNDED_PROFILE_RECON")


if __name__ == "__main__":
    main()
