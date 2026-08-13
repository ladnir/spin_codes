#!/usr/bin/env python3
"""Exact-length degree-8 outer bound for a packet-weight profile.

Give every eight-bit value of Hamming weight ``j`` the same fugacity ``t_j``.
For one fixed internal bit position, the unnormalized symbol masses are

    A(t) = sum_j C(7,j-1) t_j   (bit one),
    C(t) = sum_j C(7,j)   t_j   (bit zero).

After degree-8 Finner, a band row of length ``n`` and factor ``f`` contributes

    sum_w C(n,w) A(t)^w C(t)^(n-w) f(w)^8.

Use separate rank-one factors for the ordinary 42/43 projected EBCH rows and
the exact 128 averaged 41/43 punctured rows.  The first 85 data coordinates
determine every full data word.  If the requested global profile is ``a_j``,
its first-two-band counts ``b_j`` obey ``0 <= b_j <= a_j`` and sum to the
known 174080 first-two-band packet slots.  The exact coefficient charge is
the minimum of ``sum_j b_j log(t_j)`` over those capped counts; a greedy fill
in increasing ``log(t_j)`` computes it exactly.  This avoids enumerating the
third band without discarding its known size.

For every nonzero fixed data message, the random graph input is uniform in 24
bits and the final packet assignment specifies all replacement bits, giving
an additional exact ``-24`` graph match charge.

The spectra and reduction are exact.  The optimized feasible witness is a
binary64 diagnostic until exported and outward checked.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from probe_packet8_constant_punctured_structured_upper import (
    BLOCKS,
    FINNER_DEGREE,
    FULL_FACTOR_WEIGHT,
    FULL_N0,
    GRAPH_BITS,
    N1,
    PUNCTURED_FACTOR_WEIGHT,
    PUNCTURED_N0,
    SCALE,
    load_full,
    load_punctured,
)
from probe_packet8_weight_profile_scalar import CLASSES, parse_profile


FIRST_TWO_PACKETS = 8 * 256 * (FULL_N0 + N1)


def optimize(full_cells, punctured_cells, profile: list[int]):
    # point = f0[1..42], f1[0..43], g0[1..41], g1[0..43], log(t_1..t_8),
    # with log(t_0)=0 fixing the common scaling gauge.
    off_f1 = FULL_N0
    off_g0 = off_f1 + N1 + 1
    off_g1 = off_g0 + PUNCTURED_N0
    off_t = off_g1 + N1 + 1
    variables = off_t + 8
    profile_vector = np.array(profile, dtype=np.float64)

    log_binomials = {
        n: np.array([math.log(math.comb(n, weight)) for weight in range(n + 1)])
        for n in (FULL_N0, PUNCTURED_N0, N1)
    }
    row_weights = {
        n: np.arange(n + 1, dtype=np.float64)
        for n in (FULL_N0, PUNCTURED_N0, N1)
    }
    log_one_coeff = np.full(9, -math.inf)
    log_zero_coeff = np.full(9, -math.inf)
    for weight in range(9):
        if weight:
            log_one_coeff[weight] = math.log(math.comb(7, weight - 1))
        if weight <= 7:
            log_zero_coeff[weight] = math.log(math.comb(7, weight))

    def unpack(point):
        f0 = np.concatenate(([0.0], point[:off_f1]))
        f1 = point[off_f1:off_g0]
        g0 = np.concatenate(([0.0], point[off_g0:off_g1]))
        g1 = point[off_g1:off_t]
        log_t = np.concatenate(([0.0], point[off_t:variables]))
        return f0, f1, g0, g1, log_t

    def minimum_subprofile(log_t):
        remaining = FIRST_TWO_PACKETS
        selected = np.zeros(9, dtype=np.float64)
        for weight in sorted(range(9), key=lambda index: (log_t[index], index)):
            take = min(remaining, profile[weight])
            selected[weight] = take
            remaining -= take
            if not remaining:
                break
        if remaining:
            raise SystemExit("structured profile outer: profile mass is too small")
        return selected

    def symbol_masses(log_t):
        log_a = logsumexp(log_one_coeff + log_t)
        log_c = logsumexp(log_zero_coeff + log_t)
        alpha = np.exp(log_one_coeff + log_t - log_a)
        beta = np.exp(log_zero_coeff + log_t - log_c)
        return float(log_a), float(log_c), alpha, beta

    def row_score(n, factor, log_a, log_c):
        weights = row_weights[n]
        return (
            log_binomials[n]
            + weights * log_a
            + (n - weights) * log_c
            + 8.0 * factor
        )

    full_fraction = FULL_FACTOR_WEIGHT / SCALE
    punctured_fraction = PUNCTURED_FACTOR_WEIGHT / SCALE

    def components(point):
        f0, f1, g0, g1, log_t = unpack(point)
        log_a, log_c, alpha, beta = symbol_masses(log_t)
        scores = (
            row_score(FULL_N0, f0, log_a, log_c),
            row_score(N1, f1, log_a, log_c),
            row_score(PUNCTURED_N0, g0, log_a, log_c),
            row_score(N1, g1, log_a, log_c),
        )
        return scores, alpha, beta

    def objective(point):
        (sf0, sf1, sg0, sg1), _alpha, _beta = components(point)
        log_t = unpack(point)[-1]
        selected = minimum_subprofile(log_t)
        return (
            full_fraction * (logsumexp(sf0) + logsumexp(sf1))
            + punctured_fraction * (logsumexp(sg0) + logsumexp(sg1))
            - float(np.dot(selected, log_t)) / SCALE
        )

    def gradient(point):
        (sf0, sf1, sg0, sg1), alpha, beta = components(point)
        log_t = unpack(point)[-1]
        selected = minimum_subprofile(log_t)
        pf0 = np.exp(sf0 - logsumexp(sf0))
        pf1 = np.exp(sf1 - logsumexp(sf1))
        pg0 = np.exp(sg0 - logsumexp(sg0))
        pg1 = np.exp(sg1 - logsumexp(sg1))

        def symbol_gradient(n, probability):
            expected_weight = float(np.dot(row_weights[n], probability))
            return expected_weight * alpha + (n - expected_weight) * beta

        t_gradient = (
            full_fraction
            * (
                symbol_gradient(FULL_N0, pf0)
                + symbol_gradient(N1, pf1)
            )
            + punctured_fraction
            * (
                symbol_gradient(PUNCTURED_N0, pg0)
                + symbol_gradient(N1, pg1)
            )
            - selected / SCALE
        )
        return np.concatenate(
            (
                8.0 * full_fraction * pf0[1:],
                8.0 * full_fraction * pf1,
                8.0 * punctured_fraction * pg0[1:],
                8.0 * punctured_fraction * pg1,
                t_gradient[1:],
            )
        )

    rows = len(full_cells) + len(punctured_cells)
    jacobian = np.zeros((rows, variables), dtype=np.float64)
    rhs = np.zeros(rows, dtype=np.float64)
    for row, (a, b, value) in enumerate(full_cells):
        if a:
            jacobian[row, a - 1] = 1.0
        jacobian[row, off_f1 + b] = 1.0
        rhs[row] = value
    base = len(full_cells)
    for local, (a, b, value) in enumerate(punctured_cells):
        row = base + local
        if a:
            jacobian[row, off_g0 + a - 1] = 1.0
        jacobian[row, off_g1 + b] = 1.0
        rhs[row] = value

    def constraints(point):
        return jacobian @ point - rhs

    initial_f1 = np.full(N1 + 1, -100.0)
    for _a, b, value in full_cells:
        initial_f1[b] = max(initial_f1[b], value)
    initial_g1 = np.full(N1 + 1, -100.0)
    for _a, b, value in punctured_cells:
        initial_g1[b] = max(initial_g1[b], value)
    empirical = np.full(9, -24.0)
    for weight, count in enumerate(profile):
        if count:
            empirical[weight] = math.log(count / CLASSES[weight])
    empirical -= empirical[0]
    initial = np.concatenate(
        (
            np.zeros(FULL_N0),
            initial_f1,
            np.zeros(PUNCTURED_N0),
            initial_g1,
            empirical[1:],
        )
    )
    bounds = [(-40.0, 40.0)] * variables
    for index in range(off_t, variables):
        bounds[index] = (-40.0, 40.0)
    result = minimize(
        objective,
        initial,
        jac=gradient,
        method="SLSQP",
        bounds=bounds,
        constraints={"type": "ineq", "fun": constraints, "jac": lambda _p: jacobian},
        options={"maxiter": 10000, "ftol": 1e-12, "disp": False},
    )

    point = np.array(result.x, copy=True)
    raw_full = float(np.min(constraints(point)[: len(full_cells)]))
    raw_punctured = float(np.min(constraints(point)[len(full_cells) :]))
    repair_full = max(0.0, -raw_full) + 1e-12
    repair_punctured = max(0.0, -raw_punctured) + 1e-12
    point[off_f1:off_g0] += repair_full
    point[off_g1:off_t] += repair_punctured
    minimum_slack = float(np.min(constraints(point)))
    return (
        result,
        unpack(point),
        raw_full,
        raw_punctured,
        repair_full,
        repair_punctured,
        minimum_slack,
        SCALE * objective(point),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument("--profile", required=True)
    parser.add_argument(
        "--full-spectrum",
        type=Path,
        default=root / "out" / "ebch85_band01_split_spectrum.csv",
    )
    parser.add_argument(
        "--punctured-spectrum",
        type=Path,
        default=root / "out" / "ebch84_punctured_band01_split_spectrum.csv",
    )
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"structured profile outer: {error}") from error
    result = optimize(
        load_full(args.full_spectrum),
        load_punctured(args.punctured_spectrum),
        profile,
    )
    (
        optimizer,
        factors,
        raw_full,
        raw_punctured,
        repair_full,
        repair_punctured,
        slack,
        coefficient_log,
    ) = result
    f0, f1, g0, g1, log_t = factors
    data_log2 = coefficient_log / math.log(2.0)
    print("packet-8 exact-length structured weight-profile outer probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"optimizer_success={optimizer.success} message={optimizer.message}")
    print(f"raw_full_minimum_slack={raw_full:.12e}")
    print(f"raw_punctured_minimum_slack={raw_punctured:.12e}")
    print(f"full_uniform_repair={repair_full:.12e}")
    print(f"punctured_uniform_repair={repair_punctured:.12e}")
    print(f"minimum_constraint_slack={slack:.12e}")
    print(f"outer_profile_data_count_log2_upper={data_log2:.12f}")
    print(f"graph_match_log2_upper={-GRAPH_BITS}")
    print(f"outer_profile_with_graph_log2_upper={data_log2-GRAPH_BITS:.12f}")
    print(f"log_t={','.join(f'{value:.9f}' for value in log_t)}")
    print(f"f0_range={float(np.min(f0)):.9f},{float(np.max(f0)):.9f}")
    print(f"f1_range={float(np.min(f1)):.9f},{float(np.max(f1)):.9f}")
    print(f"g0_range={float(np.min(g0)):.9f},{float(np.max(g0)):.9f}")
    print(f"g1_range={float(np.min(g1)):.9f},{float(np.max(g1)):.9f}")
    print("model=ACTUAL_EXACT_LENGTH_DEGREE8_PACKET_WEIGHT_PROFILE_FINNER")
    print("status=DIAGNOSTIC_BINARY64_STRUCTURED_WEIGHT_PROFILE_OUTER")


if __name__ == "__main__":
    main()
