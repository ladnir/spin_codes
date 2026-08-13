#!/usr/bin/env python3
"""Bounded diagnostic survey of frozen-witness reuse on full-support g=8 profiles.

The script optimizes a small, deterministic bank of combined outer/inner
witnesses at seed profiles.  It then freezes each witness as

    constant - <profile, charge> - normalization(profile)

and evaluates the bank on a larger deterministic profile catalogue.  This is
a binary64 reconnaissance tool.  It does not interpolate between profiles and
does not certify any part of the eight-dimensional profile simplex.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from packet_group_drive_stratified import feasible_profile_count, profile_classes
from packet_group_native import load_point_caps, load_shared_drive_apply
from packet_group_outer_profile import D, N, atom_count, normalization_log2
from packet_group_profile_bound import inner_probability, split_cap_table, tune_inner_density
from probe_packet_group_conditioned_row_outer import (
    evaluate_conditioned_outer,
    load_split_spectrum,
)


GROUP_BITS = 8
CLASSES = GROUP_BITS + 1
ATOMS = atom_count(GROUP_BITS)
MINIMUM_NONZERO_OUTER_WEIGHT = 21
FEASIBLE_PROFILES = feasible_profile_count(
    GROUP_BITS, N, MINIMUM_NONZERO_OUTER_WEIGHT
)
UNIFORM_TARGET_LOG2 = -40.0 - math.log2(FEASIBLE_PROFILES)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def allocate_positive(weights: list[Fraction] | np.ndarray) -> list[int]:
    """Allocate ATOMS with stable largest remainder and one item per class."""

    if len(weights) != CLASSES:
        raise ValueError("full-support allocation requires nine weights")
    values = [Fraction(value) for value in weights]
    mass = sum(values)
    if mass <= 0 or any(value < 0 for value in values):
        raise ValueError("profile weights must be nonnegative with positive mass")
    remaining = ATOMS - CLASSES
    scaled = [value * remaining / mass for value in values]
    floors = [value.numerator // value.denominator for value in scaled]
    missing = remaining - sum(floors)
    order = sorted(
        range(CLASSES),
        key=lambda index: (-(scaled[index] - floors[index]), index),
    )
    result = [value + 1 for value in floors]
    for index in order[:missing]:
        result[index] += 1
    if sum(result) != ATOMS or any(value <= 0 for value in result):
        raise AssertionError("positive largest-remainder allocation failed")
    return result


def binomial_weights(numerator: int, denominator: int) -> list[Fraction]:
    p = Fraction(numerator, denominator)
    if not 0 < p < 1:
        raise ValueError("binomial parameter must lie strictly between zero and one")
    return [
        Fraction(math.comb(GROUP_BITS, weight))
        * p**weight
        * (1 - p) ** (GROUP_BITS - weight)
        for weight in range(CLASSES)
    ]


def add_profile(
    rows: list[dict], seen: set[tuple[int, ...]], name: str, family: str, weights
) -> None:
    profile = allocate_positive(weights)
    key = tuple(profile)
    if key in seen:
        return
    seen.add(key)
    rows.append({"name": name, "family": family, "profile": profile})


def seed_catalogue() -> list[dict]:
    """Return the original deterministic density-oriented seed bank."""

    rows: list[dict] = []
    seen: set[tuple[int, ...]] = set()
    for numerator, denominator in (
        (1, 16),
        (15, 16),
        (1, 4),
        (3, 4),
        (1, 2),
    ):
        add_profile(
            rows,
            seen,
            f"seed_binomial_{numerator}_{denominator}",
            "seed_binomial",
            binomial_weights(numerator, denominator),
        )
    add_profile(
        rows,
        seen,
        "seed_class_uniform",
        "seed_class_uniform",
        [Fraction(1, CLASSES)] * CLASSES,
    )
    add_profile(
        rows,
        seen,
        "seed_alternating_skew",
        "seed_skew",
        [Fraction(16 if index % 2 else 1, 1) for index in range(CLASSES)],
    )
    return rows


def stable_skew_vectors(
    random_seed: int, concentration: float, count: int
) -> np.ndarray:
    """Return a stream stable across changes to other skew-family sizes."""

    tags = {0.03: 3, 0.1: 10, 0.3: 30, 1.0: 100}
    tag = tags[concentration]
    sequence = np.random.SeedSequence([random_seed, GROUP_BITS, tag])
    rng = np.random.default_rng(sequence)
    return rng.dirichlet(np.full(CLASSES, concentration), count)


def targeted_seed_catalogue(random_seed: int) -> list[dict]:
    """Return stable non-binomial seeds for support-family enrichment."""

    rows: list[dict] = []
    seen: set[tuple[int, ...]] = set()
    add_profile(
        rows,
        seen,
        "class_uniform",
        "class_uniform",
        [Fraction(1, CLASSES)] * CLASSES,
    )
    ridge = [Fraction(0)] * CLASSES
    ridge[0] = ridge[8] = Fraction(1, 2)
    add_profile(
        rows,
        seen,
        "near_ridge_0_8_50_100",
        "near_boundary_ridge",
        ridge,
    )
    for dominant in (0, 4, 8):
        weights = [Fraction(0)] * CLASSES
        weights[dominant] = Fraction(1)
        add_profile(
            rows,
            seen,
            f"near_vertex_{dominant}",
            "near_boundary_vertex",
            weights,
        )
    for concentration in (0.03, 0.3, 1.0):
        vector = stable_skew_vectors(random_seed, concentration, 1)[0]
        add_profile(
            rows,
            seen,
            f"skew_{concentration:g}_000",
            f"seeded_skew_{concentration:g}",
            [Fraction(float(value)) for value in vector],
        )
    return rows


def survey_catalogue(
    random_seed: int, random_per_concentration: int, ridge_ratios: tuple[int, ...]
) -> list[dict]:
    """Build deterministic binomial, uniform, skew, and near-face families."""

    rows: list[dict] = []
    seen: set[tuple[int, ...]] = set()

    # A dense one-dimensional curve through full support.
    for numerator in range(1, 64):
        add_profile(
            rows,
            seen,
            f"binomial_{numerator}_64",
            "binomial_density",
            binomial_weights(numerator, 64),
        )

    add_profile(
        rows,
        seen,
        "class_uniform",
        "class_uniform",
        [Fraction(1, CLASSES)] * CLASSES,
    )
    for base in (2, 4, 8, 16, 32, 64):
        add_profile(
            rows,
            seen,
            f"class_uniform_tilt_up_{base}",
            "class_uniform_tilt",
            [Fraction(base**index, 1) for index in range(CLASSES)],
        )
        add_profile(
            rows,
            seen,
            f"class_uniform_tilt_down_{base}",
            "class_uniform_tilt",
            [Fraction(base ** (GROUP_BITS - index), 1) for index in range(CLASSES)],
        )

    # These rows are one count away from a proper face in every non-dominant
    # coordinate.  They remain full-support profiles.
    for dominant in range(CLASSES):
        weights = [Fraction(0)] * CLASSES
        weights[dominant] = Fraction(1)
        add_profile(
            rows,
            seen,
            f"near_vertex_{dominant}",
            "near_boundary_vertex",
            weights,
        )

    # Pair ridges probe many directions near codimension-seven faces.
    for left in range(CLASSES):
        for right in range(left + 1, CLASSES):
            for ratio in ridge_ratios:
                weights = [Fraction(0)] * CLASSES
                weights[left] = Fraction(ratio, 100)
                weights[right] = Fraction(100 - ratio, 100)
                add_profile(
                    rows,
                    seen,
                    f"near_ridge_{left}_{right}_{ratio}_100",
                    "near_boundary_ridge",
                    weights,
                )

    # A fixed PCG64 stream makes the skewed interior family reproducible.
    for concentration in (0.03, 0.1, 0.3, 1.0):
        for index, vector in enumerate(
            stable_skew_vectors(random_seed, concentration, random_per_concentration)
        ):
            # Convert binary64 draws to their exact dyadic values before the
            # stable integer allocation.
            add_profile(
                rows,
                seen,
                f"skew_{concentration:g}_{index:03d}",
                f"seeded_skew_{concentration:g}",
                [Fraction(float(value)) for value in vector],
            )

    return rows


def select_seed_profiles(
    profiles: list[dict],
    random_seed: int,
    preset: str,
    requested_names: list[str],
    requested_families: list[str],
    limit: int,
) -> tuple[list[dict], list[dict]]:
    """Select a deterministic union of preset, exact-name, and family seeds."""

    density = seed_catalogue()
    targeted = targeted_seed_catalogue(random_seed)
    available: list[dict] = []
    seen_names: set[str] = set()
    for descriptor in density + targeted + profiles:
        if descriptor["name"] not in seen_names:
            seen_names.add(descriptor["name"])
            available.append(descriptor)
    by_name = {row["name"]: row for row in available}
    available_families = {row["family"] for row in available}
    unknown_names = sorted(set(requested_names) - set(by_name))
    unknown_families = sorted(set(requested_families) - available_families)
    if unknown_names:
        raise ValueError(f"unknown seed profile names: {', '.join(unknown_names)}")
    if unknown_families:
        raise ValueError(f"unknown seed families: {', '.join(unknown_families)}")

    selected_names: list[str] = []
    if preset in ("density", "enriched"):
        selected_names.extend(row["name"] for row in density)
    if preset == "enriched":
        selected_names.extend(row["name"] for row in targeted)
    selected_names.extend(requested_names)
    requested_family_set = set(requested_families)
    selected_names.extend(
        row["name"] for row in available if row["family"] in requested_family_set
    )

    selected: list[dict] = []
    seen_profiles: set[tuple[int, ...]] = set()
    for name in selected_names:
        descriptor = by_name[name]
        key = tuple(descriptor["profile"])
        if key not in seen_profiles:
            seen_profiles.add(key)
            selected.append(descriptor)
    if limit:
        selected = selected[:limit]
    if not selected:
        raise ValueError("seed selection is empty")
    return selected, available


def empirical_start(profile: list[int]) -> np.ndarray:
    classes = np.asarray(profile_classes(GROUP_BITS), dtype=np.float64)
    counts = np.asarray(profile, dtype=np.float64)
    positive = counts > 0.0
    if not np.any(positive):
        raise ValueError("outer seed profile has empty support")
    log_values = np.full(CLASSES, -30.0, dtype=np.float64)
    log_values[positive] = np.log(counts[positive]) - np.log(classes[positive])
    reference = log_values[0] if positive[0] else float(np.max(log_values[positive]))
    log_values[positive] -= reference
    return np.clip(log_values, -30.0, 30.0)


def optimize_outer_seed(
    profile: list[int],
    spectra: tuple[np.ndarray, np.ndarray, np.ndarray],
    mode: str,
    maximum_iterations: int,
    maximum_evaluations: int,
) -> dict:
    """Bounded local optimization of one conditioned-row outer witness."""

    spectrum01, punctured01, spectrum12 = spectra

    def unpack(point: np.ndarray):
        variables = np.concatenate(([0.0], point[:GROUP_BITS]))
        band1 = float(point[GROUP_BITS])
        if mode == "symmetric":
            coefficients = (1.0 - band1, band1, band1)
            theta = float(point[GROUP_BITS + 1])
        else:
            fraction = float(point[GROUP_BITS + 1])
            band2 = band1 + (0.999 - band1) * fraction
            coefficients = (1.0 - band1, band1, band2)
            theta = float(point[GROUP_BITS + 2])
        return variables, coefficients, theta

    def objective(point: np.ndarray) -> float:
        variables, coefficients, theta = unpack(point)
        value, _details = evaluate_conditioned_outer(
            profile,
            variables,
            coefficients[1],
            theta,
            spectrum01,
            punctured01,
            spectrum12,
            1,
            coefficients,
            GROUP_BITS,
        )
        return value

    base = empirical_start(profile)[1:]
    if mode == "symmetric":
        starts = (
            np.concatenate((base, [0.6, 0.5])),
            np.asarray([0.0] * GROUP_BITS + [0.6, 0.5]),
        )
        bounds = [(-40.0, 40.0)] * GROUP_BITS + [(0.500001, 0.999), (0.0, 1.0)]
    else:
        starts = (
            np.concatenate((base, [0.6, 0.95, 0.5])),
            np.asarray([0.0] * GROUP_BITS + [0.6, 0.95, 0.5]),
        )
        bounds = (
            [(-40.0, 40.0)] * GROUP_BITS
            + [(0.500001, 0.998), (0.0, 1.0), (0.0, 1.0)]
        )
    candidates = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "maxiter": maximum_iterations,
                "maxfun": maximum_evaluations,
                "ftol": 1e-13,
                "gtol": 1e-8,
            },
        )
        candidates.append((objective(result.x), result))
    value, result = min(candidates, key=lambda row: row[0])
    variables, coefficients, theta = unpack(result.x)
    checked, details = evaluate_conditioned_outer(
        profile,
        variables,
        coefficients[1],
        theta,
        spectrum01,
        punctured01,
        spectrum12,
        1,
        coefficients,
        GROUP_BITS,
    )
    if checked != value:
        raise AssertionError("conditioned-row outer replay changed its value")
    return {
        "outer_log2": float(checked),
        "log_variables": variables.tolist(),
        "band_coefficients": list(coefficients),
        "pair_cauchy_theta": theta,
        "optimizer_success": bool(result.success),
        "optimizer_iterations": int(result.nit),
        "optimizer_evaluations": int(result.nfev),
        "optimizer_message": str(result.message),
        "normal_tile_log": float(details["normal_tile_log"]),
        "graph_average_log": float(details["graph_average_log"]),
    }


def optimize_seed_witness(
    descriptor: dict,
    split_caps,
    spectra: tuple[np.ndarray, np.ndarray, np.ndarray],
    args: argparse.Namespace,
) -> dict:
    started = time.perf_counter()
    profile = descriptor["profile"]
    inner = tune_inner_density(
        GROUP_BITS,
        profile,
        split_caps,
        poles=tuple(args.inner_poles),
        exponents=tuple(args.inner_exponents),
        iterations=args.inner_screen_iterations,
    )
    pole = float(inner[1])
    fugacities = np.asarray(inner[3], dtype=np.float64)
    inner_value, inner_details = inner_probability(
        GROUP_BITS,
        profile,
        pole,
        fugacities,
        split_caps,
        iterations=args.inner_final_iterations,
    )
    outer = optimize_outer_seed(
        profile,
        spectra,
        args.outer_mode,
        args.outer_max_iterations,
        args.outer_max_evaluations,
    )

    profile_vector = np.asarray(profile, dtype=np.float64)
    outer_charge = np.asarray(outer["log_variables"], dtype=np.float64) / math.log(2.0)
    inner_charge = np.log2(fugacities)
    outer_constant = float(outer["outer_log2"]) + float(profile_vector @ outer_charge)
    inner_constant = float(inner_details["inner_mgf_log2"]) - D * math.log2(pole)
    combined_constant = outer_constant + inner_constant
    combined_charge = outer_charge + inner_charge
    normalization = float(normalization_log2(GROUP_BITS, profile_vector)[0])
    replay = combined_constant - float(profile_vector @ combined_charge) - normalization
    local_combined = float(outer["outer_log2"]) + float(inner_value)
    if not math.isclose(replay, local_combined, rel_tol=0.0, abs_tol=2e-8):
        raise AssertionError("frozen affine witness does not replay its seed value")

    return {
        "name": descriptor["name"],
        "seed_family": descriptor["family"],
        "seed_profile": profile,
        "seed_physical_weight": sum(index * count for index, count in enumerate(profile)),
        "seed_combined_log2": local_combined,
        "seed_uniform_target_margin_bits": UNIFORM_TARGET_LOG2 - local_combined,
        "outer": outer,
        "inner": {
            "probability_log2": float(inner_value),
            "pole": pole,
            "fugacities": fugacities.tolist(),
            "grid_evaluations": len(args.inner_poles) * len(args.inner_exponents),
            "screen_iterations": args.inner_screen_iterations,
            "final_iterations": args.inner_final_iterations,
            "inner_mgf_log2": float(inner_details["inner_mgf_log2"]),
            "collatz_best_iteration": int(inner_details["collatz_best_iteration"]),
            "worst_state": int(inner_details["worst_state"]),
        },
        "affine": {
            "constant_log2": combined_constant,
            "charge_log2": combined_charge.tolist(),
            "normalization": "multinomial packet-profile normalization",
        },
        "elapsed_seconds": time.perf_counter() - started,
    }


def evaluate_bank(profiles: list[dict], witnesses: list[dict], show_worst: int) -> dict:
    profile_matrix = np.asarray([row["profile"] for row in profiles], dtype=np.float64)
    normalizations = normalization_log2(GROUP_BITS, profile_matrix)
    columns = []
    for witness in witnesses:
        affine = witness["affine"]
        columns.append(
            float(affine["constant_log2"])
            - profile_matrix @ np.asarray(affine["charge_log2"], dtype=np.float64)
            - normalizations
        )
    matrix = np.asarray(columns)
    leaders = np.argmin(matrix, axis=0)
    values = matrix[leaders, np.arange(len(profiles))]
    margins = UNIFORM_TARGET_LOG2 - values
    order = np.argsort(values)[::-1]
    leader_counts = Counter(int(index) for index in leaders)
    active_shares = np.asarray(
        [count / len(profiles) for count in leader_counts.values()], dtype=np.float64
    )
    entropy = -float(np.sum(active_shares * np.log2(active_shares)))

    rows = []
    for index, descriptor in enumerate(profiles):
        leader = int(leaders[index])
        rows.append(
            {
                **descriptor,
                "support_size": sum(count > 0 for count in descriptor["profile"]),
                "physical_weight": sum(
                    weight * count for weight, count in enumerate(descriptor["profile"])
                ),
                "best_witness": witnesses[leader]["name"],
                "combined_log2": float(values[index]),
                "uniform_target_margin_bits": float(margins[index]),
            }
        )
    reuse = []
    for index, witness in enumerate(witnesses):
        selected = np.flatnonzero(leaders == index)
        reuse.append(
            {
                "witness": witness["name"],
                "leader_count": int(len(selected)),
                "leader_fraction": float(len(selected) / len(profiles)),
                "worst_led_combined_log2": (
                    None if not len(selected) else float(np.max(values[selected]))
                ),
                "worst_led_margin_bits": (
                    None if not len(selected) else float(np.min(margins[selected]))
                ),
            }
        )
    family_counts = Counter(row["family"] for row in rows)
    family_worst = {}
    for family in sorted(family_counts):
        selected = [row for row in rows if row["family"] == family]
        worst = max(selected, key=lambda row: row["combined_log2"])
        family_worst[family] = {
            "profiles": len(selected),
            "worst_name": worst["name"],
            "worst_combined_log2": worst["combined_log2"],
            "worst_margin_bits": worst["uniform_target_margin_bits"],
            "best_witness": worst["best_witness"],
        }
    return {
        "profiles_evaluated": len(profiles),
        "covered_at_uniform_profile_target": int(np.sum(margins >= 0.0)),
        "uncovered_at_uniform_profile_target": int(np.sum(margins < 0.0)),
        "active_witnesses": len(leader_counts),
        "leader_entropy_bits": entropy,
        "maximum_leader_fraction": float(np.max(active_shares)),
        "reuse": reuse,
        "family_worst": family_worst,
        "worst_profiles": [rows[int(index)] for index in order[:show_worst]],
        "rows": rows,
    }


def compare_with_density_baseline(
    profiles: list[dict], witnesses: list[dict], enriched: dict, show_worst: int
) -> dict:
    """Measure coverage and exponent improvements from non-density witnesses."""

    density_names = {row["name"] for row in seed_catalogue()}
    baseline_witnesses = [row for row in witnesses if row["name"] in density_names]
    added_witnesses = [row for row in witnesses if row["name"] not in density_names]
    if not baseline_witnesses or not added_witnesses:
        return {
            "available": False,
            "reason": "comparison requires both density and non-density witnesses",
            "density_witnesses": [row["name"] for row in baseline_witnesses],
            "added_witnesses": [row["name"] for row in added_witnesses],
        }

    baseline = evaluate_bank(profiles, baseline_witnesses, show_worst)
    baseline_rows = baseline["rows"]
    enriched_rows = enriched["rows"]
    improvements = np.asarray(
        [
            before["combined_log2"] - after["combined_log2"]
            for before, after in zip(baseline_rows, enriched_rows)
        ],
        dtype=np.float64,
    )
    newly_covered = [
        after
        for before, after in zip(baseline_rows, enriched_rows)
        if before["uniform_target_margin_bits"] < 0.0
        and after["uniform_target_margin_bits"] >= 0.0
    ]
    family_gains = Counter(row["family"] for row in newly_covered)
    improved_order = np.argsort(improvements)[::-1]
    return {
        "available": True,
        "density_witnesses": [row["name"] for row in baseline_witnesses],
        "added_witnesses": [row["name"] for row in added_witnesses],
        "baseline_covered": baseline["covered_at_uniform_profile_target"],
        "enriched_covered": enriched["covered_at_uniform_profile_target"],
        "newly_covered": len(newly_covered),
        "newly_covered_by_family": dict(sorted(family_gains.items())),
        "profiles_strictly_improved": int(np.sum(improvements > 1e-8)),
        "median_improvement_bits": float(np.median(improvements)),
        "maximum_improvement_bits": float(np.max(improvements)),
        "baseline_worst_profiles": baseline["worst_profiles"],
        "largest_improvements": [
            {
                "name": enriched_rows[int(index)]["name"],
                "family": enriched_rows[int(index)]["family"],
                "improvement_bits": float(improvements[int(index)]),
                "enriched_combined_log2": enriched_rows[int(index)]["combined_log2"],
                "enriched_best_witness": enriched_rows[int(index)]["best_witness"],
            }
            for index in improved_order[:show_worst]
        ],
    }


def parse_float_list(text: str) -> list[float]:
    values = [float(item) for item in text.split(",") if item]
    if not values or any(not math.isfinite(value) for value in values):
        raise argparse.ArgumentTypeError("expected a nonempty comma-separated float list")
    return values


def parse_int_tuple(text: str) -> tuple[int, ...]:
    values = tuple(int(item) for item in text.split(",") if item)
    if not values or any(not 0 < value < 100 for value in values):
        raise argparse.ArgumentTypeError("ridge ratios must lie strictly between 0 and 100")
    return values


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--seed-preset", choices=("none", "density", "enriched"), default="density"
    )
    parser.add_argument(
        "--seed-name",
        action="append",
        default=[],
        help="add one exact profile name from the deterministic catalogue",
    )
    parser.add_argument(
        "--seed-family",
        action="append",
        default=[],
        help="add every profile in one deterministic catalogue family",
    )
    parser.add_argument(
        "--seed-limit",
        type=int,
        default=0,
        help="stable prefix cap after selection; zero keeps every selected seed",
    )
    parser.add_argument("--random-seed", type=int, default=20260813)
    parser.add_argument("--random-per-concentration", type=int, default=48)
    parser.add_argument("--ridge-ratios", type=parse_int_tuple, default=(1, 10, 50, 90, 99))
    parser.add_argument("--inner-poles", type=parse_float_list, default=[0.2, 0.5, 0.8])
    parser.add_argument("--inner-exponents", type=parse_float_list, default=[0.25, 0.5, 0.75])
    parser.add_argument("--inner-screen-iterations", type=int, default=2)
    parser.add_argument("--inner-final-iterations", type=int, default=6)
    parser.add_argument("--outer-mode", choices=("symmetric", "asymmetric"), default="asymmetric")
    parser.add_argument("--outer-max-iterations", type=int, default=160)
    parser.add_argument("--outer-max-evaluations", type=int, default=2400)
    parser.add_argument("--show-worst", type=int, default=12)
    parser.add_argument("--catalogue-only", action="store_true")
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
    args = parser.parse_args()
    if (
        args.seed_limit < 0
        or args.random_per_concentration < 0
        or args.inner_screen_iterations <= 0
        or args.inner_final_iterations <= 0
        or args.outer_max_iterations <= 0
        or args.outer_max_evaluations <= 0
        or args.show_worst <= 0
    ):
        raise SystemExit("g=8 full-support reuse survey: invalid positive budget")
    if any(not 0.0 < value < 1.0 for value in args.inner_poles):
        raise SystemExit("g=8 full-support reuse survey: inner poles must lie in (0,1)")

    profiles = survey_catalogue(
        args.random_seed, args.random_per_concentration, args.ridge_ratios
    )
    try:
        seeds, available_seeds = select_seed_profiles(
            profiles,
            args.random_seed,
            args.seed_preset,
            args.seed_name,
            args.seed_family,
            args.seed_limit,
        )
    except ValueError as error:
        raise SystemExit(f"g=8 full-support reuse survey: {error}") from error
    if not profiles or any(
        len(row["profile"]) != CLASSES
        or sum(row["profile"]) != ATOMS
        or any(count <= 0 for count in row["profile"])
        for row in seeds + profiles
    ):
        raise AssertionError("catalogue contains a malformed full-support profile")
    if args.catalogue_only:
        counts = Counter(row["family"] for row in profiles)
        available_counts = Counter(row["family"] for row in available_seeds)
        print(f"seed_profiles={len(seeds)}")
        print(f"seed_names={json.dumps([row['name'] for row in seeds])}")
        print(
            "available_seed_families="
            f"{json.dumps(dict(sorted(available_counts.items())), sort_keys=True)}"
        )
        print(f"survey_profiles={len(profiles)}")
        print(f"families={json.dumps(dict(sorted(counts.items())), sort_keys=True)}")
        print("status=CATALOGUE_ONLY_NO_OPTIMIZATION")
        return

    native_apply = load_shared_drive_apply()
    native_caps = load_point_caps()
    if native_apply is None or native_caps is None:
        raise SystemExit(
            "g=8 full-support reuse survey: native inner apply and point-cap kernels are required"
        )
    spectra = (
        load_split_spectrum(args.spectrum01),
        load_split_spectrum(args.punctured01, count_field="pair_count", divisor=42),
        load_split_spectrum(args.spectrum12),
    )
    split_caps = split_cap_table()
    witnesses = []
    started = time.perf_counter()
    for index, seed in enumerate(seeds, 1):
        witness = optimize_seed_witness(seed, split_caps, spectra, args)
        witnesses.append(witness)
        print(
            f"seed={index}/{len(seeds)} name={witness['name']} "
            f"combined={witness['seed_combined_log2']:.9f} "
            f"margin={witness['seed_uniform_target_margin_bits']:.9f} "
            f"elapsed={witness['elapsed_seconds']:.3f}s",
            flush=True,
        )

    survey = evaluate_bank(profiles, witnesses, args.show_worst)
    survey["marginal_vs_density_preset"] = compare_with_density_baseline(
        profiles, witnesses, survey, args.show_worst
    )
    source_paths = {
        "script": Path(__file__),
        "conditioned_row_evaluator": Path(__file__).with_name(
            "probe_packet_group_conditioned_row_outer.py"
        ),
        "inner_evaluator": Path(__file__).with_name("packet_group_profile_bound.py"),
        "native_library": Path(native_apply._packet_group_library_path),
        "spectrum01": args.spectrum01,
        "punctured01": args.punctured01,
        "spectrum12": args.spectrum12,
    }
    report = {
        "status": "DIAGNOSTIC_BINARY64_G8_FULL_SUPPORT_WITNESS_REUSE",
        "scope": (
            "finite deterministic full-support sample only; frozen affine binary64 "
            "witnesses; no interpolation, outward replay, or domain coverage claim"
        ),
        "group_bits": GROUP_BITS,
        "atoms": ATOMS,
        "N": N,
        "D": D,
        "minimum_nonzero_outer_weight": MINIMUM_NONZERO_OUTER_WEIGHT,
        "feasible_profile_count": FEASIBLE_PROFILES,
        "uniform_profile_target_log2": UNIFORM_TARGET_LOG2,
        "configuration": {
            "argv": sys.argv,
            "seed_preset": args.seed_preset,
            "requested_seed_names": args.seed_name,
            "requested_seed_families": args.seed_family,
            "seed_limit": args.seed_limit,
            "random_seed": args.random_seed,
            "random_per_concentration": args.random_per_concentration,
            "ridge_ratios": list(args.ridge_ratios),
            "inner_poles": args.inner_poles,
            "inner_exponents": args.inner_exponents,
            "inner_screen_iterations": args.inner_screen_iterations,
            "inner_final_iterations": args.inner_final_iterations,
            "outer_mode": args.outer_mode,
            "outer_max_iterations": args.outer_max_iterations,
            "outer_max_evaluations": args.outer_max_evaluations,
        },
        "sources": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in source_paths.items()
        },
        "native_backend_required_and_loaded": True,
        "witnesses": witnesses,
        "survey": survey,
        "elapsed_seconds": time.perf_counter() - started,
        "blockers_to_coverage_claim": [
            "the sample is finite and does not cover the full-support simplex",
            "the optimized and frozen parameters use binary64 arithmetic",
            "the witness bank has no independent outward replay",
            "no exact slab or BSP ownership ledger aggregates all full-support profiles",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    worst = survey["worst_profiles"][0]
    print(
        f"profiles={survey['profiles_evaluated']} active_witnesses={survey['active_witnesses']} "
        f"max_share={survey['maximum_leader_fraction']:.6f}"
    )
    print(
        f"worst={worst['name']} combined={worst['combined_log2']:.9f} "
        f"margin={worst['uniform_target_margin_bits']:.9f} "
        f"witness={worst['best_witness']}"
    )
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_BINARY64_G8_FULL_SUPPORT_WITNESS_REUSE")


if __name__ == "__main__":
    main()
