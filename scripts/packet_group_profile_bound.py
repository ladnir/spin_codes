#!/usr/bin/env python3
"""Fixed-profile shared-drive bound and inexpensive cross-g tuning."""

from __future__ import annotations

import math

import numpy as np

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from packet_group_drive_stratified import (
    block_histograms,
    point_caps,
    profile_classes,
    validate_group,
)
from packet_group_outer_profile import (
    D,
    optimize_outer,
    normalization_log2,
)
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import best_witness
from probe_packet8_weight_transfer import build_split_caps, load_exact
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


GRAPH_REPLACEMENTS = 128


def split_cap_table():
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    return build_split_caps(spectrum, load_exact(EXACT_SLICES))


def inner_probability(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    pole: float,
    fugacities: np.ndarray,
    split_caps,
    *,
    iterations: int = 48,
) -> tuple[float, dict[str, float | int]]:
    validate_group(group_bits)
    histograms = block_histograms(group_bits, fugacities)
    caps = point_caps(group_bits, fugacities)
    kernel = SharedDriveStratifiedKernel(
        histograms, caps, split_caps, pole
    )
    (
        eigenvalue,
        domination,
        values,
        worst_state,
        best_iteration,
        mgf,
    ) = best_witness(kernel, iterations, INNER_BLOCKS)
    charge = sum(
        count * math.log2(float(fugacity))
        for count, fugacity in zip(profile, fugacities)
        if count
    )
    normalization = float(
        normalization_log2(group_bits, np.asarray(profile, dtype=np.float64))[0]
    )
    probability = mgf - charge - normalization - D * math.log2(pole)
    return probability, {
        "inner_mgf_log2": mgf,
        "charge_log2": charge,
        "normalization_log2": normalization,
        "lambda_log2": math.log2(eigenvalue),
        "domination_log2": math.log2(domination),
        "worst_state": worst_state,
        "collatz_selection": "best_finite_block_bound_on_power_trajectory",
        "collatz_best_iteration": best_iteration,
        "collatz_max_iterations": iterations,
        "collatz_vector": values.tolist(),
    }


def density_fugacities(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    exponent: float,
    *,
    pseudocount: float = 0.5,
) -> np.ndarray:
    classes = np.asarray(profile_classes(group_bits), dtype=np.float64)
    profile_array = np.asarray(profile, dtype=np.float64)
    log_density = np.log(profile_array + pseudocount) - np.log(classes)
    log_density -= float(np.max(log_density))
    values = np.exp(exponent * log_density)
    return np.maximum(values, 1e-300)


def shaped_density_fugacities(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    parameters: np.ndarray,
    *,
    pseudocount: float = 0.5,
) -> np.ndarray:
    """Density tilt plus linear/quadratic normalized-weight corrections."""

    exponent, linear, quadratic = map(float, parameters)
    classes = np.asarray(profile_classes(group_bits), dtype=np.float64)
    profile_array = np.asarray(profile, dtype=np.float64)
    log_density = np.log(profile_array + pseudocount) - np.log(classes)
    coordinate = np.arange(group_bits + 1, dtype=np.float64) / group_bits - 0.5
    log_values = (
        exponent * log_density
        + linear * coordinate
        + quadratic * coordinate * coordinate
    )
    log_values -= float(np.max(log_values))
    return np.maximum(np.exp(log_values), 1e-300)


def refine_inner_shaped_density(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    split_caps,
    *,
    initial_pole: float,
    initial_exponent: float,
    initial_parameters: np.ndarray | None = None,
    coordinate_iterations: int = 4,
    passes: int = 2,
    witness_iterations: int = 48,
    linear_bound: float = 12.0,
    quadratic_bound: float = 24.0,
):
    """Small deterministic coordinate refinement used at a margin crossing."""

    point = (
        np.asarray(initial_parameters, dtype=np.float64).copy()
        if initial_parameters is not None
        else np.asarray(
            [initial_exponent, 0.0, 0.0, math.log(initial_pole)],
            dtype=np.float64,
        )
    )
    if point.shape != (4,):
        raise ValueError("shaped-density initial point must have four entries")
    bounds = (
        (0.0, 1.5),
        (-linear_bound, linear_bound),
        (-quadratic_bound, quadratic_bound),
        (math.log(0.03), math.log(0.7)),
    )
    evaluations = 0

    def evaluate(candidate: np.ndarray):
        nonlocal evaluations
        evaluations += 1
        fugacities = shaped_density_fugacities(group_bits, profile, candidate[:3])
        pole = math.exp(float(candidate[3]))
        probability, details = inner_probability(
            group_bits,
            profile,
            pole,
            fugacities,
            split_caps,
            iterations=witness_iterations,
        )
        return probability, pole, fugacities, details, candidate.copy()

    best = evaluate(point)
    for _pass in range(passes):
        start = point.copy()
        for coordinate_index, (low_bound, high_bound) in enumerate(bounds):
            low = low_bound
            high = high_bound
            for _ in range(coordinate_iterations):
                left = (2.0 * low + high) / 3.0
                right = (low + 2.0 * high) / 3.0
                left_point = point.copy()
                right_point = point.copy()
                left_point[coordinate_index] = left
                right_point[coordinate_index] = right
                left_value = evaluate(left_point)
                right_value = evaluate(right_point)
                if left_value[0] <= right_value[0]:
                    high = right
                    if left_value[0] < best[0]:
                        best = left_value
                else:
                    low = left
                    if right_value[0] < best[0]:
                        best = right_value
            point[coordinate_index] = (low + high) / 2.0
            current = evaluate(point)
            if current[0] < best[0]:
                best = current
        if float(np.max(np.abs(point - start))) < 1e-4:
            break
    final = evaluate(point)
    if final[0] < best[0]:
        best = final
    return (*best, evaluations)


def refine_inner_sparse_fugacities(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    split_caps,
    *,
    initial_pole: float,
    initial_fugacities: np.ndarray,
    log_steps: tuple[float, ...] = (
        math.log(16.0),
        math.log(4.0),
        math.log(2.0),
    ),
    witness_iterations: int = 12,
    zero_floor: float = 1e-64,
    verbose: bool = False,
):
    """Coordinate-refine the occupied classes of a sparse profile.

    A zero-count class has no coefficient-extraction charge, so its optimal
    fugacity may lie on the zero boundary.  Pinning such classes to a tiny
    positive floor both represents that limit and avoids wasting coordinates
    on pseudocount artifacts.  Common rescaling of every fugacity cancels from
    the homogeneous coefficient bound; the most populous occupied class is
    therefore held fixed as an anchor.
    """

    profile_array = np.asarray(profile, dtype=np.int64)
    if profile_array.shape != (group_bits + 1,):
        raise ValueError("profile has the wrong number of classes")
    occupied = np.flatnonzero(profile_array)
    if not len(occupied):
        raise ValueError("the all-zero profile is excluded")
    anchor = int(occupied[np.argmax(profile_array[occupied])])
    movable = [int(index) for index in occupied if int(index) != anchor]
    values = np.asarray(initial_fugacities, dtype=np.float64).copy()
    if values.shape != (group_bits + 1,):
        raise ValueError("initial fugacities have the wrong number of classes")
    values[profile_array == 0] = zero_floor
    values = np.maximum(values, 1e-300)
    values /= values[anchor]
    evaluations = 0

    def evaluate(candidate: np.ndarray):
        nonlocal evaluations
        evaluations += 1
        probability, details = inner_probability(
            group_bits,
            profile,
            initial_pole,
            candidate,
            split_caps,
            iterations=witness_iterations,
        )
        return probability, details

    best_probability, best_details = evaluate(values)
    trace = []
    for step in log_steps:
        accepted = 0
        for index in movable:
            candidates = []
            for direction in (-1.0, 1.0):
                candidate = values.copy()
                candidate[index] *= math.exp(direction * step)
                candidate = np.maximum(candidate, 1e-300)
                probability, details = evaluate(candidate)
                candidates.append((probability, candidate, details, direction))
            choice = min(candidates, key=lambda row: row[0])
            if choice[0] < best_probability:
                best_probability = choice[0]
                values = choice[1]
                best_details = choice[2]
                accepted += 1
                trace.append(
                    {
                        "step": float(step),
                        "class": index,
                        "direction": float(choice[3]),
                        "probability_log2": float(best_probability),
                    }
                )
                if verbose:
                    print(
                        "sparse-refine",
                        float(step),
                        index,
                        float(best_probability),
                        flush=True,
                    )
        if verbose:
            print(
                "sparse-step",
                float(step),
                accepted,
                float(best_probability),
                flush=True,
            )
    return (
        best_probability,
        initial_pole,
        values,
        best_details,
        trace,
        evaluations,
        anchor,
    )


def refine_inner_full_fugacities(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    split_caps,
    *,
    initial_pole: float,
    initial_fugacities: np.ndarray,
    log_steps: tuple[float, ...] = (
        math.log(16.0),
        math.log(4.0),
        math.log(2.0),
        math.log(math.sqrt(2.0)),
    ),
    passes: int = 2,
    witness_iterations: int = 24,
):
    """Coordinate-refine every class fugacity plus the pole.

    The shaped-density family is intentionally cheap, but a full-support
    profile can need a tilt outside its quadratic three-parameter manifold.
    Common fugacity scaling cancels from coefficient extraction, so the most
    populous class is fixed to one and only the remaining class ratios move.
    """

    profile_array = np.asarray(profile, dtype=np.int64)
    if profile_array.shape != (group_bits + 1,):
        raise ValueError("profile has the wrong number of classes")
    if np.any(profile_array <= 0):
        raise ValueError("full-fugacity refinement requires full support")
    anchor = int(np.argmax(profile_array))
    movable = [index for index in range(group_bits + 1) if index != anchor]
    values = np.maximum(np.asarray(initial_fugacities, dtype=np.float64), 1e-300)
    if values.shape != (group_bits + 1,):
        raise ValueError("initial fugacities have the wrong number of classes")
    values = values / values[anchor]
    pole = float(initial_pole)
    evaluations = 0

    def evaluate(candidate_pole: float, candidate_values: np.ndarray):
        nonlocal evaluations
        evaluations += 1
        probability, details = inner_probability(
            group_bits,
            profile,
            candidate_pole,
            candidate_values,
            split_caps,
            iterations=witness_iterations,
        )
        return probability, details

    best_probability, best_details = evaluate(pole, values)
    trace = []
    for _pass in range(passes):
        accepted_in_pass = 0
        for step in log_steps:
            for index in movable:
                candidates = []
                for direction in (-1.0, 1.0):
                    candidate = values.copy()
                    candidate[index] *= math.exp(direction * step)
                    probability, details = evaluate(pole, candidate)
                    candidates.append((probability, candidate, details, direction))
                choice = min(candidates, key=lambda row: row[0])
                if choice[0] < best_probability:
                    best_probability = choice[0]
                    values = choice[1]
                    best_details = choice[2]
                    accepted_in_pass += 1
                    trace.append(
                        {
                            "kind": "fugacity",
                            "step": float(step),
                            "class": index,
                            "direction": float(choice[3]),
                            "probability_log2": float(best_probability),
                        }
                    )
            pole_candidates = []
            for direction in (-1.0, 1.0):
                candidate_pole = min(
                    0.95, max(0.01, pole * math.exp(direction * step))
                )
                probability, details = evaluate(candidate_pole, values)
                pole_candidates.append(
                    (probability, candidate_pole, details, direction)
                )
            choice = min(pole_candidates, key=lambda row: row[0])
            if choice[0] < best_probability:
                best_probability = choice[0]
                pole = choice[1]
                best_details = choice[2]
                accepted_in_pass += 1
                trace.append(
                    {
                        "kind": "pole",
                        "step": float(step),
                        "direction": float(choice[3]),
                        "probability_log2": float(best_probability),
                    }
                )
        if accepted_in_pass == 0:
            break
    return (
        best_probability,
        pole,
        values,
        best_details,
        trace,
        evaluations,
        anchor,
    )


def tune_inner_density(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    split_caps,
    *,
    poles=(0.1, 0.2, 0.3, 0.4, 0.5),
    exponents=(0.0, 0.25, 0.5, 0.75, 1.0),
    iterations: int = 48,
):
    candidates = []
    for exponent in exponents:
        fugacities = density_fugacities(group_bits, profile, exponent)
        for pole in poles:
            probability, details = inner_probability(
                group_bits,
                profile,
                pole,
                fugacities,
                split_caps,
                iterations=iterations,
            )
            candidates.append(
                (probability, pole, exponent, fugacities.copy(), details)
            )
    return min(candidates, key=lambda row: row[0])


def outer_probability(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    *,
    optimize_band_coefficients: bool = True,
    fast: bool = False,
):
    coefficient_log, result = optimize_outer(
        group_bits,
        profile,
        optimize_band_coefficients=optimize_band_coefficients,
        fast=fast,
    )
    log_variables = np.concatenate(([0.0], result.x[:group_bits]))
    replacement_log = GRAPH_REPLACEMENTS * max(
        log_variables[new] - log_variables[old]
        for old in range(group_bits + 1)
        for new in range(group_bits + 1)
        if abs(new - old) <= 1
    )
    value = (coefficient_log + replacement_log) / math.log(2.0)
    return value, {
        "unpunctured_outer_log2": coefficient_log / math.log(2.0),
        "graph_replacement_log2": replacement_log / math.log(2.0),
        "log_variables": log_variables.tolist(),
        "band1_coefficient": (
            float(result.x[group_bits]) if optimize_band_coefficients else 0.5
        ),
        "optimizer_success": bool(result.success),
    }
