#!/usr/bin/env python3
"""Hybrid three-band degree-8 Finner bound for constant packets.

For band weights ``a,b,c`` a full EBCH word must satisfy both projected-code
conditions on bands (0,1) and (1,2).  Hence its membership probability under
the independent within-band coordinate permutations is at most

    min(p01(a,b), p12(b,c)) <= sqrt(p01(a,b) p12(b,c)).

Dominate the two exact pair probabilities by positive rank-one factors,

    p01(a,b) <= exp(u0[a] + u1[b]),
    p12(b,c) <= exp(v1[b] + v2[c]).

The geometric mean then separates into one factor per band.  Every physical
constant-packet variable occurs in exactly eight block factors in its band,
so degree-8 Finner bounds each real tile.  This retains the actual balanced
lane partitions and all three fixed band tilings.  The reduction is rigorous;
the optimized witness remains binary64 diagnostic until exported and checked
outward.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


N0 = 42
N1 = 43
N2 = 43
TILES = 256


def load_cells(path: Path, left: int, right: int):
    cells = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            a = int(row["band0_weight"])
            b = int(row["band1_weight"])
            count = int(row["count"])
            if not (0 <= a <= left and 0 <= b <= right):
                raise SystemExit("three-band structured: spectrum shape mismatch")
            cells.append(
                (
                    a,
                    b,
                    math.log(count)
                    - math.log(math.comb(left, a))
                    - math.log(math.comb(right, b)),
                )
            )
    if not cells:
        raise SystemExit("three-band structured: empty spectrum")
    return cells


def optimize(cells01, cells12, cumulative_packets: int | None):
    # point = u0[1..42], u1[0..43], v1[1..43], v2[0..43].
    off_u0 = 0
    off_u1 = off_u0 + N0
    off_v1 = off_u1 + N1 + 1
    off_v2 = off_v1 + N1
    factor_variables = off_v2 + N2 + 1
    variables = factor_variables + (1 if cumulative_packets is not None else 0)

    def unpack(point):
        u0 = np.concatenate(([0.0], point[off_u0:off_u1]))
        u1 = point[off_u1:off_v1]
        v1 = np.concatenate(([0.0], point[off_v1:off_v2]))
        v2 = point[off_v2:factor_variables]
        log_q = float(point[-1]) if cumulative_packets is not None else 0.0
        return u0, u1, v1, v2, log_q

    log_binomial = [
        np.array([math.log(math.comb(n, a)) for a in range(n + 1)])
        for n in (N0, N1, N2)
    ]

    def objective(point):
        u0, u1, v1, v2, log_q = unpack(point)
        scores0 = log_binomial[0] + np.arange(N0 + 1) * log_q + 4.0 * u0
        scores1 = (
            log_binomial[1] + np.arange(N1 + 1) * log_q + 4.0 * (u1 + v1)
        )
        scores2 = log_binomial[2] + np.arange(N2 + 1) * log_q + 4.0 * v2
        value = 8.0 * (
            logsumexp(scores0) + logsumexp(scores1) + logsumexp(scores2)
        )
        if cumulative_packets is not None:
            value -= (cumulative_packets / TILES) * log_q
        return value

    def gradient(point):
        u0, u1, v1, v2, log_q = unpack(point)
        scores0 = log_binomial[0] + np.arange(N0 + 1) * log_q + 4.0 * u0
        scores1 = (
            log_binomial[1] + np.arange(N1 + 1) * log_q + 4.0 * (u1 + v1)
        )
        scores2 = log_binomial[2] + np.arange(N2 + 1) * log_q + 4.0 * v2
        probability0 = np.exp(scores0 - logsumexp(scores0))
        probability1 = np.exp(scores1 - logsumexp(scores1))
        probability2 = np.exp(scores2 - logsumexp(scores2))
        result = np.concatenate(
            (
                32.0 * probability0[1:],
                32.0 * probability1,
                32.0 * probability1[1:],
                32.0 * probability2,
            )
        )
        if cumulative_packets is not None:
            q_gradient = 8.0 * (
                float(np.dot(np.arange(N0 + 1), probability0))
                + float(np.dot(np.arange(N1 + 1), probability1))
                + float(np.dot(np.arange(N2 + 1), probability2))
            ) - cumulative_packets / TILES
            result = np.concatenate((result, [q_gradient]))
        return result

    rows = len(cells01) + len(cells12)
    jacobian = np.zeros((rows, variables), dtype=np.float64)
    rhs = np.zeros(rows, dtype=np.float64)
    for row, (a, b, value) in enumerate(cells01):
        if a:
            jacobian[row, off_u0 + a - 1] = 1.0
        jacobian[row, off_u1 + b] = 1.0
        rhs[row] = value
    base = len(cells01)
    for local, (b, c, value) in enumerate(cells12):
        row = base + local
        if b:
            jacobian[row, off_v1 + b - 1] = 1.0
        jacobian[row, off_v2 + c] = 1.0
        rhs[row] = value

    def constraints(point):
        return jacobian @ point - rhs

    # Feasible rank-one starts with the left factor fixed to one.
    initial_u1 = np.full(N1 + 1, -100.0)
    for _a, b, value in cells01:
        initial_u1[b] = max(initial_u1[b], value)
    initial_v2 = np.full(N2 + 1, -100.0)
    for _b, c, value in cells12:
        initial_v2[c] = max(initial_v2[c], value)
    initial = np.concatenate(
        (np.zeros(N0), initial_u1, np.zeros(N1), initial_v2)
    )
    if cumulative_packets is not None:
        initial_q = math.log(max(cumulative_packets / (TILES * 8 * 128), 1e-12))
        initial = np.concatenate((initial, [initial_q]))
    bounds = [(-40.0, 40.0)] * variables
    if cumulative_packets is not None:
        bounds[-1] = (-40.0, 0.0)
    result = minimize(
        objective,
        initial,
        jac=gradient,
        method="SLSQP",
        bounds=bounds,
        constraints={"type": "ineq", "fun": constraints, "jac": lambda _p: jacobian},
        options={"maxiter": 8000, "ftol": 1e-12, "disp": False},
    )

    point = np.array(result.x, copy=True)
    raw_slack = float(np.min(constraints(point)))
    # Raising both right factors uniformly repairs both constraint families.
    repair = max(0.0, -raw_slack) + 1e-12
    point[off_u1:off_v1] += repair
    point[off_v2:factor_variables] += repair
    slack = float(np.min(constraints(point)))
    reported = objective(point)
    if cumulative_packets is not None:
        u0, u1, v1, v2, log_q = unpack(point)
        scores0 = log_binomial[0] + np.arange(N0 + 1) * log_q + 4.0 * u0
        scores1 = (
            log_binomial[1] + np.arange(N1 + 1) * log_q + 4.0 * (u1 + v1)
        )
        scores2 = log_binomial[2] + np.arange(N2 + 1) * log_q + 4.0 * v2
        full_log = TILES * 8.0 * (
            logsumexp(scores0) + logsumexp(scores1) + logsumexp(scores2)
        )
        nonzero_log = (
            full_log + math.log1p(-math.exp(-full_log))
            if full_log > 50.0
            else math.log(math.expm1(full_log))
        )
        reported = nonzero_log - cumulative_packets * log_q
    return result, unpack(point), raw_slack, repair, slack, reported


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument(
        "--spectrum01",
        type=Path,
        default=root / "out" / "ebch85_band01_split_spectrum.csv",
    )
    parser.add_argument(
        "--spectrum12",
        type=Path,
        default=root / "out" / "ebch86_band12_split_spectrum.csv",
    )
    parser.add_argument("--cumulative-packets", type=int)
    args = parser.parse_args()
    if args.cumulative_packets is not None and args.cumulative_packets <= 0:
        raise SystemExit("three-band structured: cumulative packets must be positive")
    cells01 = load_cells(args.spectrum01, N0, N1)
    cells12 = load_cells(args.spectrum12, N1, N2)
    result, factors, raw_slack, repair, slack, tile_log = optimize(
        cells01, cells12, args.cumulative_packets
    )
    u0, u1, v1, v2, log_q = factors
    full_log2 = (
        tile_log / math.log(2.0)
        if args.cumulative_packets is not None
        else TILES * tile_log / math.log(2.0)
    )
    print("constant-packet hybrid three-band structured upper probe")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"raw_minimum_constraint_slack={raw_slack:.12e}")
    print(f"uniform_right_factor_repair={repair:.12e}")
    print(f"minimum_constraint_slack={slack:.12e}")
    print(f"one_three_band_tile_log2={tile_log/math.log(2.0):.12f}")
    print(f"full_expected_solution_count_log2_upper={full_log2:.12f}")
    if args.cumulative_packets is not None:
        print(f"cumulative_total_packet_support={args.cumulative_packets}")
        print(f"log_q={log_q:.12f}")
    print(f"u0_range={float(np.min(u0)):.9f},{float(np.max(u0)):.9f}")
    print(f"u1_range={float(np.min(u1)):.9f},{float(np.max(u1)):.9f}")
    print(f"v1_range={float(np.min(v1)):.9f},{float(np.max(v1)):.9f}")
    print(f"v2_range={float(np.min(v2)):.9f},{float(np.max(v2)):.9f}")
    print("model=ACTUAL_THREE_BAND_BALANCED_PARTITIONS_DEGREE8_FINNER")
    print("status=DIAGNOSTIC_BINARY64_HYBRID_THREE_BAND_CONSTANT_PACKET_UPPER")


if __name__ == "__main__":
    main()
