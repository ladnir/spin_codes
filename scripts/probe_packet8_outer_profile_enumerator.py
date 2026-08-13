#!/usr/bin/env python3
"""Diagnostic outer packet-weight-profile enumerator bound.

This probe combines the exact EBCH [128,64,22] total spectrum with the exact
random-lane packet profile of one physical 64-bit group.  It deliberately
uses only the total binary weight of an outer word, so it is valid for every
grouping/layout before the graph/puncture correction (and may be loose).

For positive profile variables ``t_j``, a group containing ``w`` one-bits has
exact lane-averaged moment

    R_w(t) = [y^w] (sum_j C(8,j) t_j y^j)^8 / C(64,w).

If ``R_w(t) <= alpha beta^w`` for every ``w``, then the complete direct sum of
16384 EBCH blocks obeys

    A_outer(t) <= alpha^32768 E_EBCH(beta)^16384.

Positive coefficient extraction consequently bounds the expected number of
outer words with any supplied packet-weight profile.  The optimization is
binary64 and the current model omits the exact 128 puncture/graph replacements,
so this is a diagnostic gate, not a theorem certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from probe_packet8_weight_profile_scalar import CLASSES, parse_profile
from punctured_ebch_outer import full_spectrum


N = 1 << 21
K = 1 << 20
M = N // 8
D = 9 * N // 100
GROUPS = N // 64
DATA_BLOCKS = K // 64
SMOOTH_TAU = 1e-4


def group_log_moments(log_t: np.ndarray) -> np.ndarray:
    scale = float(np.max(log_t))
    polynomial = np.array(
        [classes * math.exp(value - scale) for classes, value in zip(CLASSES, log_t)],
        dtype=np.float64,
    )
    power = np.array([1.0])
    for _ in range(8):
        power = np.convolve(power, polynomial)
    result = np.empty(65, dtype=np.float64)
    for weight in range(65):
        result[weight] = (
            8.0 * scale
            + math.log(float(power[weight]))
            - math.log(math.comb(64, weight))
        )
    return result


def spectrum_logs() -> tuple[np.ndarray, np.ndarray]:
    spectrum = full_spectrum()
    weights = np.array([weight for weight, count in enumerate(spectrum) if count])
    logs = np.array([math.log(spectrum[weight]) for weight in weights])
    return weights, logs


def optimize_outer(profile: list[int]):
    spectrum_weights, spectrum_count_logs = spectrum_logs()
    profile_vector = np.array(profile, dtype=np.float64)

    def components(point: np.ndarray, *, smooth: bool):
        log_t = np.concatenate(([0.0], point[:8]))
        log_beta = float(point[8])
        group_logs = group_log_moments(log_t)
        envelope_rows = group_logs - np.arange(65) * log_beta
        if smooth:
            log_alpha = SMOOTH_TAU * logsumexp(envelope_rows / SMOOTH_TAU)
        else:
            log_alpha = float(np.max(envelope_rows))
        log_enumerator = logsumexp(
            spectrum_count_logs + spectrum_weights * log_beta
        )
        outer_log = GROUPS * log_alpha + DATA_BLOCKS * log_enumerator
        coefficient_log = outer_log - float(np.dot(profile_vector, log_t))
        return coefficient_log, log_alpha, log_enumerator, envelope_rows

    # Empirical per-concrete-byte frequencies give the unconstrained profile
    # saddle.  The total-weight tilt initializes beta at the corresponding bit
    # density.  The all-zero start is retained as an independent basin check.
    empirical = np.zeros(9, dtype=np.float64)
    if profile[0]:
        for weight in range(1, 9):
            if profile[weight]:
                empirical[weight] = math.log(
                    profile[weight] / CLASSES[weight] / profile[0]
                )
            else:
                empirical[weight] = -16.0
    density = sum(weight * count for weight, count in enumerate(profile)) / N
    beta_start = math.log(density / (1.0 - density)) if 0.0 < density < 1.0 else 0.0
    starts = [np.zeros(9), np.concatenate((empirical[1:], [beta_start]))]
    candidates = []
    bounds = [(-20.0, 20.0)] * 9
    for start in starts:
        result = minimize(
            lambda point: components(point, smooth=True)[0],
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "ftol": 1e-13, "gtol": 1e-9},
        )
        exact_components = components(result.x, smooth=False)
        candidates.append((exact_components[0], result, exact_components))
    return min(candidates, key=lambda row: row[0])


def log2_profile_normalization(profile: list[int]) -> float:
    result = math.lgamma(M + 1) - sum(math.lgamma(count + 1) for count in profile)
    result += sum(
        count * math.log(classes)
        for count, classes in zip(profile, CLASSES)
        if count
    )
    return result / math.log(2.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"packet8 outer profile: {error}") from error

    coefficient_log, result, components = optimize_outer(profile)
    _value, log_alpha, log_enumerator, envelope_rows = components
    coefficient_log2 = coefficient_log / math.log(2.0)
    normalization_log2 = log2_profile_normalization(profile)
    # Exact pole z=1/10, as in certify_packet8_profile_high_branches.py.
    inner_probability_log2 = (
        N * math.log2(11)
        - (N - D) * math.log2(10)
        - normalization_log2
    )
    combined_log2 = coefficient_log2 + inner_probability_log2
    active_rows = np.flatnonzero(
        envelope_rows >= np.max(envelope_rows) - 1e-7
    )
    print("packet-8 outer profile enumerator probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"outer_profile_count_log2_upper={coefficient_log2:.9f}")
    print(f"trivial_family_bits={K}")
    print(f"outer_saving_vs_2^K={K-coefficient_log2:.9f}")
    print(f"profile_normalization_log2={normalization_log2:.9f}")
    print(f"bijection_inner_probability_log2={inner_probability_log2:.9f}")
    print(f"combined_first_moment_log2={combined_log2:.9f}")
    print(f"log_alpha={log_alpha:.12f} log_local_enumerator={log_enumerator:.12f}")
    print(f"active_envelope_weights={','.join(map(str, active_rows))}")
    print(
        "log_t="
        + ",".join(f"{value:.9f}" for value in np.concatenate(([0.0], result.x[:8])))
    )
    print(f"log_beta={result.x[8]:.9f}")
    print("model=UNPUNCTURED_16384_BLOCK_TOTAL_SPECTRUM")
    print("status=DIAGNOSTIC_BINARY64_OUTER_PROFILE_ENUMERATOR")


if __name__ == "__main__":
    main()
