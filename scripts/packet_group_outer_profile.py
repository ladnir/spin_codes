#!/usr/bin/env python3
"""Outer packet-profile formulas parameterized by permutation atom width g."""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln, logsumexp

from packet_group_drive_stratified import BLOCK_BITS, profile_classes, validate_group
from punctured_ebch_outer import full_spectrum


N = 1 << 21
K = 1 << 20
D = 9 * N // 100
GROUPS = N // BLOCK_BITS
BINOMIAL_LOG_PROBABILITIES = np.asarray(
    [math.log(math.comb(BLOCK_BITS, weight)) - BLOCK_BITS * math.log(2.0)
     for weight in range(BLOCK_BITS + 1)]
)


def atom_count(group_bits: int) -> int:
    validate_group(group_bits)
    return N // group_bits


def group_log_moments(group_bits: int, log_variables: np.ndarray) -> np.ndarray:
    atoms = validate_group(group_bits)
    if log_variables.shape != (group_bits + 1,):
        raise ValueError("outer packet variable vector has wrong shape")
    log_polynomial = np.asarray(
        [
            math.log(classes) + value
            for classes, value in zip(profile_classes(group_bits), log_variables)
        ],
        dtype=np.float64,
    )
    log_power = np.array([0.0])
    for _ in range(atoms):
        next_power = np.full(len(log_power) + group_bits, -np.inf)
        for left, left_value in enumerate(log_power):
            next_power[left : left + group_bits + 1] = np.logaddexp(
                next_power[left : left + group_bits + 1],
                left_value + log_polynomial,
            )
        log_power = next_power
    return np.asarray(
        [
            log_power[weight] - math.log(math.comb(BLOCK_BITS, weight))
            for weight in range(BLOCK_BITS + 1)
        ]
    )


def optimize_outer(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    *,
    optimize_band_coefficients: bool = True,
    log_variable_bound: float = 40.0,
    fast: bool = False,
):
    validate_group(group_bits)
    if len(profile) != group_bits + 1 or sum(profile) != atom_count(group_bits):
        raise ValueError("outer profile has wrong dimension or mass")
    profile_vector = np.asarray(profile, dtype=np.float64)

    def objective(point: np.ndarray) -> float:
        log_variables = np.concatenate(([0.0], point[:group_bits]))
        moments = group_log_moments(group_bits, log_variables)
        band1 = float(point[group_bits]) if optimize_band_coefficients else 0.5
        band0 = 1.0 - band1
        norm0 = band0 * logsumexp(
            BINOMIAL_LOG_PROBABILITIES + moments / band0
        )
        norm1 = band1 * logsumexp(
            BINOMIAL_LOG_PROBABILITIES + moments / band1
        )
        return (
            K * math.log(2.0)
            + 256.0 * (42.0 * norm0 + 86.0 * norm1)
            - float(profile_vector @ log_variables)
        )

    classes = profile_classes(group_bits)
    empirical = np.full(group_bits + 1, -24.0)
    for weight, count in enumerate(profile):
        if count:
            empirical[weight] = math.log(count / classes[weight])
    empirical -= empirical[0]
    base_starts = [empirical[1:]] if fast else [np.zeros(group_bits), empirical[1:]]
    starts = []
    if optimize_band_coefficients:
        for coefficient in ((0.9,) if fast else (0.5, 0.6, 0.75, 0.9)):
            starts.extend(
                np.concatenate((start, [coefficient])) for start in base_starts
            )
        bounds = [(-log_variable_bound, log_variable_bound)] * group_bits + [
            (0.5, 0.999)
        ]
    else:
        starts = base_starts
        bounds = [(-log_variable_bound, log_variable_bound)] * group_bits
    candidates = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-9},
        )
        candidates.append((objective(result.x), result))
    return min(candidates, key=lambda row: row[0])


def evaluate_linear_outer(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    log_variables: np.ndarray,
    band1_coefficient: float,
) -> float:
    """Evaluate one frozen valid linear-BL outer witness in natural logs."""

    validate_group(group_bits)
    if len(profile) != group_bits + 1 or sum(profile) != atom_count(group_bits):
        raise ValueError("linear outer profile has wrong dimension or mass")
    if log_variables.shape != (group_bits + 1,) or log_variables[0] != 0.0:
        raise ValueError("linear outer variables must have shape g+1 and t0=1")
    if not 0.5 <= band1_coefficient < 1.0:
        raise ValueError("invalid linear-BL band coefficient")
    moments = group_log_moments(group_bits, log_variables)
    band0 = 1.0 - band1_coefficient
    norm0 = band0 * logsumexp(
        BINOMIAL_LOG_PROBABILITIES + moments / band0
    )
    norm1 = band1_coefficient * logsumexp(
        BINOMIAL_LOG_PROBABILITIES + moments / band1_coefficient
    )
    return (
        K * math.log(2.0)
        + 256.0 * (42.0 * norm0 + 86.0 * norm1)
        - float(np.asarray(profile, dtype=np.float64) @ log_variables)
    )


def optimize_spectrum_outer(
    group_bits: int,
    profile: list[int] | tuple[int, ...],
    *,
    smoothing: float = 1e-4,
):
    """Total-EBCH-spectrum envelope, valid independently of tile grouping."""

    validate_group(group_bits)
    if len(profile) != group_bits + 1 or sum(profile) != atom_count(group_bits):
        raise ValueError("spectrum outer profile has wrong dimension or mass")
    spectrum = full_spectrum()
    spectrum_weights = np.asarray(
        [weight for weight, count in enumerate(spectrum) if count], dtype=np.float64
    )
    spectrum_logs = np.asarray(
        [math.log(spectrum[int(weight)]) for weight in spectrum_weights]
    )
    profile_vector = np.asarray(profile, dtype=np.float64)

    def components(point: np.ndarray, *, smooth: bool):
        log_variables = np.concatenate(([0.0], point[:group_bits]))
        log_beta = float(point[group_bits])
        moments = group_log_moments(group_bits, log_variables)
        rows = moments - np.arange(BLOCK_BITS + 1) * log_beta
        log_alpha = (
            smoothing * logsumexp(rows / smoothing)
            if smooth
            else float(np.max(rows))
        )
        log_enumerator = logsumexp(
            spectrum_logs + spectrum_weights * log_beta
        )
        outer_log = GROUPS * log_alpha + (K // 64) * log_enumerator
        return (
            outer_log - float(profile_vector @ log_variables),
            log_alpha,
            log_enumerator,
            rows,
        )

    classes = profile_classes(group_bits)
    empirical = np.zeros(group_bits + 1)
    if profile[0]:
        for weight in range(1, group_bits + 1):
            empirical[weight] = (
                math.log(profile[weight] / classes[weight] / profile[0])
                if profile[weight]
                else -16.0
            )
    density = sum(weight * count for weight, count in enumerate(profile)) / N
    beta_start = (
        math.log(density / (1.0 - density)) if 0.0 < density < 1.0 else 0.0
    )
    starts = [
        np.zeros(group_bits + 1),
        np.concatenate((empirical[1:], [beta_start])),
    ]
    bounds = [(-40.0, 40.0)] * (group_bits + 1)
    candidates = []
    for start in starts:
        result = minimize(
            lambda point: components(point, smooth=True)[0],
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 3000, "ftol": 1e-13, "gtol": 1e-9},
        )
        exact = components(result.x, smooth=False)
        candidates.append((exact[0], result, exact))
    return min(candidates, key=lambda row: row[0])


def pure_profile_total_weight_outer(
    group_bits: int,
    class_weight: int,
    *,
    replacements: int = 128,
):
    """Exact-spectrum Chernoff bound for one pure physical atom profile.

    Replacing at most ``replacements`` data bits by graph bits changes total
    weight by at most that amount.  We sum the common-pole bounds over the
    resulting pre-replacement weight interval.
    """

    validate_group(group_bits)
    if not 0 <= class_weight <= group_bits:
        raise ValueError("pure class weight outside atom")
    physical_weight = atom_count(group_bits) * class_weight
    return total_weight_outer(physical_weight, replacements=replacements)


def total_weight_outer(physical_weight: int, *, replacements: int = 128):
    """Exact-spectrum common-pole bound from physical total weight alone."""

    if not 0 <= physical_weight <= N:
        raise ValueError("physical total weight outside code length")
    low = max(0, physical_weight - replacements)
    high = min(N, physical_weight + replacements)
    spectrum = full_spectrum()
    weights = np.asarray(
        [weight for weight, count in enumerate(spectrum) if count], dtype=np.float64
    )
    logs = np.asarray([math.log(spectrum[int(weight)]) for weight in weights])

    def objective(log_pole: float) -> float:
        local = logsumexp(logs + weights * log_pole)
        interval_weights = np.arange(low, high + 1, dtype=np.float64)
        interval = logsumexp(-interval_weights * log_pole)
        return (K // 64) * local + interval

    result = minimize_scalar(
        objective,
        bounds=(-40.0, 40.0),
        method="bounded",
        options={"xatol": 1e-12, "maxiter": 1000},
    )
    return objective(float(result.x)) / math.log(2.0), result, (low, high)


def fixed_total_weight_outer(physical_weight: int, *, replacements: int = 128):
    """Freeze the optimized total-weight pole as a reusable affine witness.

    The graph replacement changes physical weight by a delta in
    ``[-replacements, replacements]``. Summing over that entire range gives a
    globally valid affine majorant, including near weights 0 and N where some
    of the overcounted shifted weights are outside the physical interval.
    """

    optimized, result, interval = total_weight_outer(
        physical_weight, replacements=replacements
    )
    log_pole = float(result.x)
    spectrum = full_spectrum()
    weights = np.asarray(
        [weight for weight, count in enumerate(spectrum) if count], dtype=np.float64
    )
    logs = np.asarray([math.log(spectrum[int(weight)]) for weight in weights])
    deltas = np.arange(-replacements, replacements + 1, dtype=np.float64)
    constant_log2 = (
        (K // 64) * logsumexp(logs + weights * log_pole)
        + logsumexp(-deltas * log_pole)
    ) / math.log(2.0)
    charge_per_bit = log_pole / math.log(2.0)
    anchor_log2 = constant_log2 - physical_weight * charge_per_bit
    return anchor_log2, {
        "constant_log2": constant_log2,
        "charge_per_bit": charge_per_bit,
        "log_pole": log_pole,
        "optimized_clamped_log2": optimized,
        "optimized_interval": interval,
        "replacements": replacements,
    }


def normalization_log2(group_bits: int, profiles: np.ndarray) -> np.ndarray:
    validate_group(group_bits)
    rows = np.atleast_2d(np.asarray(profiles, dtype=np.float64))
    if rows.shape[1] != group_bits + 1:
        raise ValueError("normalization profile has wrong dimension")
    classes = np.asarray(profile_classes(group_bits), dtype=np.float64)
    return (
        gammaln(atom_count(group_bits) + 1)
        - np.sum(gammaln(rows + 1), axis=1)
        + rows @ np.log(classes)
    ) / math.log(2.0)


def profile_count_value(group_bits: int) -> int:
    return math.comb(atom_count(group_bits) + group_bits, group_bits)


def graph_replacement_ratio(variables: tuple[int, ...]) -> Fraction:
    """Worst monomial ratio when one physical bit replaces one data bit."""

    group_bits = len(variables) - 1
    validate_group(group_bits)
    if any(value <= 0 for value in variables):
        raise ValueError("outer variables must be positive")
    return max(
        Fraction(variables[new], variables[old])
        for old in range(group_bits + 1)
        for new in range(group_bits + 1)
        if abs(new - old) <= 1
    )


def symmetric_s2(group_bits: int, variables: tuple[int, ...]) -> Fraction:
    atoms = validate_group(group_bits)
    if len(variables) != group_bits + 1:
        raise ValueError("outer variables have wrong dimension")
    polynomial = [
        math.comb(group_bits, weight) * variables[weight]
        for weight in range(group_bits + 1)
    ]
    power = [1]
    for _ in range(atoms):
        result = [0] * (len(power) + len(polynomial) - 1)
        for left, left_value in enumerate(power):
            for right, right_value in enumerate(polynomial):
                result[left + right] += left_value * right_value
        power = result
    return sum(
        Fraction(value * value, math.comb(BLOCK_BITS, weight) * (1 << BLOCK_BITS))
        for weight, value in enumerate(power)
    )
