#!/usr/bin/env python3
"""Linear Brascamp--Lieb packet-profile outer bound.

The exact kernel certificate proves that the three fixed-band kernels of one
EBCH message block form a direct sum.  For every subspace ``V`` of the global
message space this implies

    sum_(band,tile) dim L_(band,tile)(V) >= 2 dim(V),

where ``L_(band,tile)`` exposes the projected EBCH rows used by one tile.
The finite-field linear Brascamp--Lieb inequality therefore applies with
coefficient 1/2 to every tile factor.  More generally, the direct-sum proof
allows band coefficients ``p_b`` whenever every pair sums to at least one.
Since each fixed-band projection is surjective, each tile norm is a product
of iid column moments.  At the symmetric point:

    S_2(t) = 2^-64 sum_w C(64,w) R_w(t)^2,
    A_outer(t) <= 2^K S_2(t)^(B/2).

This script performs binary64 coefficient optimization for a supplied exact
packet-weight profile.  The inequality is rigorous; floating optimization and
the graph/puncture replacement remain diagnostic until outward hardened.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from probe_packet8_three_band_profile_enumerator import (
    BINOMIAL_LOG_PROBABILITIES,
    D,
    GROUPS,
    K,
    N,
    group_log_moments,
    log2_profile_normalization,
)
from probe_packet8_weight_profile_scalar import CLASSES, parse_profile


EXPONENT = 2.0


def optimize_outer(
    profile: list[int],
    *,
    optimize_band_coefficients: bool = False,
    log_fugacity_bound: float = 40.0,
):
    if not 0.0 < log_fugacity_bound <= 40.0:
        raise ValueError("outer log-fugacity bound must lie in (0,40]")
    profile_vector = np.array(profile, dtype=np.float64)

    def objective(point: np.ndarray):
        log_t = np.concatenate(([0.0], point[:8]))
        group_logs = group_log_moments(log_t)
        band1_coefficient = float(point[8]) if optimize_band_coefficients else 0.5
        band0_coefficient = 1.0 - band1_coefficient
        band0_norm = band0_coefficient * logsumexp(
            BINOMIAL_LOG_PROBABILITIES + group_logs / band0_coefficient
        )
        band1_norm = band1_coefficient * logsumexp(
            BINOMIAL_LOG_PROBABILITIES + group_logs / band1_coefficient
        )
        return (
            K * math.log(2.0)
            + 256.0 * (42.0 * band0_norm + 86.0 * band1_norm)
            - float(np.dot(profile_vector, log_t))
        )

    empirical = np.full(9, -24.0)
    for weight, count in enumerate(profile):
        if count:
            empirical[weight] = math.log(count / CLASSES[weight])
    empirical -= empirical[0]
    base_starts = [np.zeros(8), empirical[1:]]
    starts = []
    if optimize_band_coefficients:
        for coefficient in (0.5, 0.6, 0.7, 0.85):
            starts.extend(
                np.concatenate((start, [coefficient])) for start in base_starts
            )
    else:
        starts = base_starts
    candidates = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=(
                [(-log_fugacity_bound, log_fugacity_bound)] * 8
                + [(0.5, 0.999)]
                if optimize_band_coefficients
                else [(-log_fugacity_bound, log_fugacity_bound)] * 8
            ),
            options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-9},
        )
        candidates.append((objective(result.x), result))
    return min(candidates, key=lambda row: row[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--optimize-band-coefficients", action="store_true")
    parser.add_argument("--log-fugacity-bound", type=float, default=40.0)
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"linear BL2 profile: {error}") from error

    coefficient_log, result = optimize_outer(
        profile,
        optimize_band_coefficients=args.optimize_band_coefficients,
        log_fugacity_bound=args.log_fugacity_bound,
    )
    coefficient_log2 = coefficient_log / math.log(2.0)
    normalization_log2 = log2_profile_normalization(profile)
    inner_log2 = (
        N * math.log2(11) - (N - D) * math.log2(10) - normalization_log2
    )
    combined = coefficient_log2 + inner_log2
    print("packet-8 three-band linear Brascamp--Lieb profile probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"outer_profile_count_log2_upper={coefficient_log2:.12f}")
    print(f"profile_normalization_log2={normalization_log2:.12f}")
    print(f"bijection_inner_probability_log2={inner_log2:.12f}")
    print(f"combined_first_moment_log2={combined:.12f}")
    print(
        "log_t="
        + ",".join(
            f"{value:.9f}" for value in np.concatenate(([0.0], result.x[:8]))
        )
    )
    band1_coefficient = (
        float(result.x[8]) if args.optimize_band_coefficients else 0.5
    )
    band0_coefficient = 1.0 - band1_coefficient
    print(
        "band_coefficients="
        f"{band0_coefficient:.12f},{band1_coefficient:.12f},"
        f"{band1_coefficient:.12f}"
    )
    print(
        "holder_exponents="
        f"{1/band0_coefficient:.12f},{1/band1_coefficient:.12f},"
        f"{1/band1_coefficient:.12f}"
    )
    print("model=UNPUNCTURED_THREE_BAND_LINEAR_BL")
    print("status=DIAGNOSTIC_BINARY64_RIGOROUS_INEQUALITY")


if __name__ == "__main__":
    main()
