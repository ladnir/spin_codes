#!/usr/bin/env python3
"""Optimize the valid fractional Finner cover for the three band tilings.

Assign fractional cover weights ``x0,x1,x2`` to the three tile-factor bands,
with ``x0+x1+x2=1``.  Generalized Finner gives an ``L_(1/xb)`` norm for band
``b``.  Bands one and two are symmetric, so set ``x1=x2=(1-x0)/2`` and jointly
optimize ``x0`` with the packet-weight fugacities.

The equal-cover value ``x0=x1=x2=1/3`` reproduces the existing exponent-3
three-band bound exactly.  This remains an unpunctured binary64 diagnostic.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logsumexp

from probe_packet8_three_band_profile_enumerator import (
    BINOMIAL_LOG_PROBABILITIES,
    K,
    group_log_moments,
)
from probe_packet8_weight_profile_scalar import parse_profile


TILES = 256
BAND0 = 42
BAND1 = 43


def optimize(profile: list[int]):
    profile_vector = np.array(profile, dtype=np.float64)

    def unpack(point):
        log_t = np.concatenate(([0.0], point[:8]))
        # Keep every exponent comfortably finite while allowing the optimizer
        # to approach a boundary if the profile calls for it.
        x0 = 0.02 + 0.94 * float(expit(point[8]))
        x1 = (1.0 - x0) / 2.0
        return log_t, x0, x1

    def objective(point):
        log_t, x0, x1 = unpack(point)
        group_logs = group_log_moments(log_t)

        def norm_term(x):
            return x * logsumexp(BINOMIAL_LOG_PROBABILITIES + group_logs / x)

        outer = K * math.log(2.0) + TILES * (
            BAND0 * norm_term(x0) + 2 * BAND1 * norm_term(x1)
        )
        return outer - float(np.dot(profile_vector, log_t))

    empirical = np.full(9, -24.0)
    classes = [math.comb(8, weight) for weight in range(9)]
    for weight, count in enumerate(profile):
        if count:
            empirical[weight] = math.log(count / classes[weight])
    empirical -= empirical[0]
    starts = [
        np.concatenate((np.zeros(8), [math.log((1 / 3 - 0.02) / (0.96 - 1 / 3))])),
        np.concatenate((empirical[1:], [math.log((1 / 3 - 0.02) / (0.96 - 1 / 3))])),
    ]
    candidates = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=[(-40.0, 40.0)] * 8 + [(-8.0, 8.0)],
            options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-9},
        )
        candidates.append((objective(result.x), result, unpack(result.x)))
    return min(candidates, key=lambda row: row[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"fractional Finner: {error}") from error
    coefficient_log, result, (log_t, x0, x1) = optimize(profile)
    print("packet-8 three-band fractional-Finner profile probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"outer_profile_count_log2_upper={coefficient_log/math.log(2.0):.12f}")
    print(f"fractional_cover={x0:.12f},{x1:.12f},{x1:.12f}")
    print(f"holder_exponents={1/x0:.9f},{1/x1:.9f},{1/x1:.9f}")
    print(f"log_t={','.join(f'{value:.9f}' for value in log_t)}")
    print("status=DIAGNOSTIC_BINARY64_VALID_FRACTIONAL_FINNER_COVER")


if __name__ == "__main__":
    main()
