#!/usr/bin/env python3
"""Three-band Finner bound for the outer packet-weight enumerator.

For profile variables ``t_j`` let ``R_w(t)`` be the exact expected product of
the eight packet variables after a uniform lane permutation of a 64-bit group
of weight ``w``:

    R_w(t) = [y^w] (sum_j C(8,j)t_j y^j)^8 / C(64,w).

Fix the within-band coordinate permutations.  A tile factor is the product of
``R_w`` over its 42 or 43 physical columns.  Each data-message block occurs in
exactly one tile factor in each of the three bands.  Finner's inequality with
exponent three therefore gives, under independent uniform block messages,

    E product_tiles F_tile
      <= product_tiles (E F_tile^3)^(1/3).

Every fixed-band EBCH projection is surjective.  Hence the 64 projected words
inside one tile are independent uniform binary vectors and its columns are
iid.  Defining

    S_3(t) = 2^-64 sum_w C(64,w) R_w(t)^3,

the complete unpunctured 16384-block data enumerator satisfies

    A_outer(t) <= 2^K S_3(t)^(B/3),   B=N/64.

Positive coefficient extraction then bounds any exact packet-weight profile.
The inequality is layout-global and does not enumerate collision components.
This probe uses binary64 optimization and has not yet inserted the exact graph
word/puncture replacement, so its output is diagnostic.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from probe_packet8_weight_profile_scalar import CLASSES, parse_profile


N = 1 << 21
K = 1 << 20
M = N // 8
D = 9 * N // 100
GROUPS = N // 64


def group_log_moments(log_t: np.ndarray) -> np.ndarray:
    scale = float(np.max(log_t))
    polynomial = np.array(
        [classes * math.exp(value - scale) for classes, value in zip(CLASSES, log_t)],
        dtype=np.float64,
    )
    power = np.array([1.0])
    for _ in range(8):
        power = np.convolve(power, polynomial)
    return np.array(
        [
            8.0 * scale
            + math.log(float(power[weight]))
            - math.log(math.comb(64, weight))
            for weight in range(65)
        ],
        dtype=np.float64,
    )


BINOMIAL_LOG_PROBABILITIES = np.array(
    [math.log(math.comb(64, weight)) - 64 * math.log(2.0) for weight in range(65)]
)


def optimize_outer(profile: list[int]):
    profile_vector = np.array(profile, dtype=np.float64)

    def objective(point: np.ndarray):
        log_t = np.concatenate(([0.0], point))
        group_logs = group_log_moments(log_t)
        log_s3 = logsumexp(BINOMIAL_LOG_PROBABILITIES + 3.0 * group_logs)
        coefficient_log = (
            K * math.log(2.0)
            + (GROUPS / 3.0) * log_s3
            - float(np.dot(profile_vector, log_t))
        )
        return coefficient_log, log_s3, group_logs

    empirical = np.zeros(8, dtype=np.float64)
    if profile[0]:
        for weight in range(1, 9):
            empirical[weight - 1] = (
                math.log(profile[weight] / CLASSES[weight] / profile[0])
                if profile[weight]
                else -24.0
            )
    starts = [np.zeros(8), empirical]
    candidates = []
    for start in starts:
        result = minimize(
            lambda point: objective(point)[0],
            start,
            method="L-BFGS-B",
            bounds=[(-30.0, 30.0)] * 8,
            options={"maxiter": 3000, "ftol": 1e-14, "gtol": 1e-9},
        )
        candidates.append((objective(result.x)[0], result, objective(result.x)))
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
    parser.add_argument(
        "--inner-probability-log2",
        type=float,
        help="optional stronger conditional OFF/LIVE probability exponent",
    )
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"packet8 three-band profile: {error}") from error

    coefficient_log, result, components = optimize_outer(profile)
    _value, log_s3, group_logs = components
    coefficient_log2 = coefficient_log / math.log(2.0)
    normalization_log2 = log2_profile_normalization(profile)
    bijection_probability_log2 = (
        N * math.log2(11)
        - (N - D) * math.log2(10)
        - normalization_log2
    )
    chosen_inner = (
        bijection_probability_log2
        if args.inner_probability_log2 is None
        else min(bijection_probability_log2, args.inner_probability_log2)
    )
    combined_log2 = coefficient_log2 + chosen_inner
    print("packet-8 three-band Finner profile enumerator probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"outer_profile_count_log2_upper={coefficient_log2:.9f}")
    print(f"trivial_family_bits={K}")
    print(f"outer_saving_vs_2^K={K-coefficient_log2:.9f}")
    print(f"profile_normalization_log2={normalization_log2:.9f}")
    print(f"bijection_inner_probability_log2={bijection_probability_log2:.9f}")
    if args.inner_probability_log2 is not None:
        print(f"supplied_inner_probability_log2={args.inner_probability_log2:.9f}")
    print(f"chosen_inner_probability_log2={chosen_inner:.9f}")
    print(f"combined_first_moment_log2={combined_log2:.9f}")
    print(f"log_s3={log_s3:.12f}")
    print(
        "log_t="
        + ",".join(f"{value:.9f}" for value in np.concatenate(([0.0], result.x)))
    )
    dominant_group_weight = int(
        np.argmax(BINOMIAL_LOG_PROBABILITIES + 3.0 * group_logs)
    )
    print(f"dominant_s3_group_weight={dominant_group_weight}")
    print("model=UNPUNCTURED_THREE_BAND_FINNER")
    print("status=DIAGNOSTIC_BINARY64_THREE_BAND_PROFILE_ENUMERATOR")


if __name__ == "__main__":
    main()
