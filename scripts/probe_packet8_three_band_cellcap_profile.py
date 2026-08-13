#!/usr/bin/env python3
"""Three-band degree-8 profile probe from rigorous BCH cell caps.

For each fixed-band weight triple ``(a,b,c)``, every exact aggregate row that
contains it is an upper bound on its unknown exact count.  Convert the smallest
such cap into a membership-probability cap ``p_cap(a,b,c)``.  The solver-free
symmetric factors

    log f_0(a) = (1/3) max_{b,c} log p_cap(a,b,c)

and analogously for the other bands satisfy
``p_cap(a,b,c) <= f_0(a) f_1(b) f_2(c)``.

With packet-weight fugacities, apply degree-8 Finner independently to all
three real packet-variable bands.  This is a deliberately coarse diagnostic
for whether a fully optimized triple rank-one domination can beat the current
exponent-3 message Finner profile bound.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from probe_bch_three_band_cell_cap_motifs import cell_caps
from probe_bch_three_band_transport_lp import BAND_SIZES, exact_rows, triples
from probe_packet8_weight_profile_scalar import parse_profile


TILES = 256
FINNER_EXPONENT = 8
ROWS_PER_BAND = 64 * TILES
ROW_FACTOR_WEIGHT = ROWS_PER_BAND / FINNER_EXPONENT


def symmetric_factors():
    states = triples()
    rows = exact_rows(states)
    caps, _diagonal, _by_total = cell_caps(states, rows)
    logs = np.full(len(states), -math.inf, dtype=np.float64)
    for index, ((a, b, c), cap) in enumerate(zip(states, caps)):
        if cap:
            logs[index] = (
                math.log(cap)
                - math.log(math.comb(BAND_SIZES[0], a))
                - math.log(math.comb(BAND_SIZES[1], b))
                - math.log(math.comb(BAND_SIZES[2], c))
            )
    factors = [np.full(size + 1, -math.inf) for size in BAND_SIZES]
    for (a, b, c), value in zip(states, logs):
        if not math.isfinite(value):
            continue
        factors[0][a] = max(factors[0][a], value / 3.0)
        factors[1][b] = max(factors[1][b], value / 3.0)
        factors[2][c] = max(factors[2][c], value / 3.0)
    if any(np.any(~np.isfinite(factor)) for factor in factors):
        raise SystemExit("three-band cellcap profile: uncovered band weight")
    minimum_slack = math.inf
    for (a, b, c), value in zip(states, logs):
        if math.isfinite(value):
            minimum_slack = min(
                minimum_slack,
                float(factors[0][a] + factors[1][b] + factors[2][c] - value),
            )
    if minimum_slack < -1e-12:
        raise SystemExit("three-band cellcap profile: factor domination failed")
    return tuple(factors), minimum_slack, sum(math.isfinite(value) for value in logs)


def optimize_profile(profile: list[int], factors):
    log_one_coeff = np.full(9, -math.inf)
    log_zero_coeff = np.full(9, -math.inf)
    for weight in range(9):
        if weight:
            log_one_coeff[weight] = math.log(math.comb(7, weight - 1))
        if weight <= 7:
            log_zero_coeff[weight] = math.log(math.comb(7, weight))
    profile_vector = np.array(profile, dtype=np.float64)
    log_binomials = [
        np.array([math.log(math.comb(size, weight)) for weight in range(size + 1)])
        for size in BAND_SIZES
    ]
    weight_vectors = [np.arange(size + 1, dtype=np.float64) for size in BAND_SIZES]

    def unpack(point):
        return np.concatenate(([0.0], point))

    def components(point):
        log_t = unpack(point)
        log_a = logsumexp(log_one_coeff + log_t)
        log_c = logsumexp(log_zero_coeff + log_t)
        alpha = np.exp(log_one_coeff + log_t - log_a)
        beta = np.exp(log_zero_coeff + log_t - log_c)
        scores = [
            log_binomial
            + weights * log_a
            + (size - weights) * log_c
            + 8.0 * factor
            for size, weights, log_binomial, factor in zip(
                BAND_SIZES, weight_vectors, log_binomials, factors
            )
        ]
        return log_t, float(log_a), float(log_c), alpha, beta, scores

    def objective(point):
        log_t, _a, _c, _alpha, _beta, scores = components(point)
        return ROW_FACTOR_WEIGHT * sum(logsumexp(score) for score in scores) - float(
            np.dot(profile_vector, log_t)
        )

    def gradient(point):
        _log_t, _a, _c, alpha, beta, scores = components(point)
        result = np.zeros(9, dtype=np.float64)
        for size, weights, score in zip(BAND_SIZES, weight_vectors, scores):
            probability = np.exp(score - logsumexp(score))
            expected = float(np.dot(weights, probability))
            result += ROW_FACTOR_WEIGHT * (
                expected * alpha + (size - expected) * beta
            )
        result -= profile_vector
        return result[1:]

    initial = np.full(9, -24.0)
    for weight, count in enumerate(profile):
        if count:
            initial[weight] = math.log(count / math.comb(8, weight))
    initial -= initial[0]
    result = minimize(
        objective,
        initial[1:],
        jac=gradient,
        method="L-BFGS-B",
        bounds=[(-40.0, 40.0)] * 8,
        options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-9},
    )
    return result, unpack(result.x), objective(result.x)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"three-band cellcap profile: {error}") from error
    factors, slack, active_cells = symmetric_factors()
    result, log_t, coefficient_log = optimize_profile(profile, factors)
    print("packet-8 three-band cell-cap structured profile probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"active_triple_cells={active_cells}")
    print(f"factor_minimum_slack={slack:.12e}")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"outer_profile_count_log2_upper={coefficient_log/math.log(2.0):.12f}")
    print(f"log_t={','.join(f'{value:.9f}' for value in log_t)}")
    for band, factor in enumerate(factors):
        print(
            f"factor{band}_range={float(np.min(factor)):.9f},"
            f"{float(np.max(factor)):.9f}"
        )
    print("status=DIAGNOSTIC_BINARY64_THREE_BAND_CELLCAP_PROFILE")


if __name__ == "__main__":
    main()
