#!/usr/bin/env python3
"""Structured-lift upper bound for constant packet assignments.

For a block whose selected physical words in bands zero and one have weights
``a,b``, independent within-band coordinate permutations make its probability
of satisfying the exact projected EBCH constraint

    p(a,b) = A[a,b] / (C(42,a) C(43,b)).

Choose positive factors ``f0,f1`` with ``p(a,b)<=f0(a)f1(b)``.  This separates
the two band tilings.  In one real band tile, summing all packet assignments
gives ``2^(8n)`` possibilities.  For fixed balanced lane partitions, every
independent packet bit occurs in exactly eight block-row factors.  Finner with
exponent eight therefore gives the rigorous tile bound

    2^(8n) E[f(Bin(n,1/2))^8]^8.

The script numerically minimizes the product of the two tile bounds subject to
the exact local domination constraints.  The combinatorial reduction is valid
for the structured per-column balanced lane partitions; binary64 constrained
optimization makes the resulting number diagnostic until an outward witness
is exported and checked.
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
TILES = 256


def load_cells(path: Path):
    cells = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            a = int(row["band0_weight"])
            b = int(row["band1_weight"])
            count = int(row["count"])
            log_probability = (
                math.log(count)
                - math.log(math.comb(N0, a))
                - math.log(math.comb(N1, b))
            )
            cells.append((a, b, log_probability))
    if not cells:
        raise SystemExit("constant structured: empty local spectrum")
    return cells


def optimize_factors(
    cells,
    density0: float | None,
    density1: float | None,
    cumulative_packets: int | None,
):
    # point = u[1..42], v[0..43], and optionally log(q0),log(q1), with
    # u[0]=0 fixing the factorization gauge.
    log_binomial0 = np.array([math.log(math.comb(N0, a)) for a in range(N0 + 1)])
    log_binomial1 = np.array([math.log(math.comb(N1, b)) for b in range(N1 + 1)])
    weighted = density0 is not None or density1 is not None
    if weighted and (density0 is None or density1 is None):
        raise ValueError("both densities must be supplied")
    if weighted and cumulative_packets is not None:
        raise ValueError("exact band densities and cumulative support are exclusive")
    cumulative = cumulative_packets is not None

    def unpack(point):
        u = np.concatenate(([0.0], point[:N0]))
        v = point[N0 : N0 + N1 + 1]
        log_q0 = float(point[-2]) if weighted else (
            float(point[-1]) if cumulative else 0.0
        )
        log_q1 = float(point[-1]) if weighted or cumulative else 0.0
        return u, v, log_q0, log_q1

    def objective(point):
        u, v, log_q0, log_q1 = unpack(point)
        tile_value = 8.0 * logsumexp(
            log_binomial0 + np.arange(N0 + 1) * log_q0 + 8.0 * u
        ) + 8.0 * logsumexp(
            log_binomial1 + np.arange(N1 + 1) * log_q1 + 8.0 * v
        )
        if weighted:
            tile_value -= 8.0 * N0 * density0 * log_q0
            tile_value -= 8.0 * N1 * density1 * log_q1
            return tile_value
        if cumulative:
            # Optimize the inclusive bound.  SLSQP may visit infeasible points
            # where full_log < 0, so subtract the exact zero term only after a
            # feasible witness has been found.  Work per tile to keep the
            # optimizer and its gradient well scaled.
            return tile_value - (cumulative_packets / TILES) * log_q0
        return tile_value

    def gradient(point):
        u, v, log_q0, log_q1 = unpack(point)
        score0 = log_binomial0 + np.arange(N0 + 1) * log_q0 + 8.0 * u
        score1 = log_binomial1 + np.arange(N1 + 1) * log_q1 + 8.0 * v
        probability0 = np.exp(score0 - logsumexp(score0))
        probability1 = np.exp(score1 - logsumexp(score1))
        tile_gradient = np.concatenate(
            (64.0 * probability0[1:], 64.0 * probability1)
        )
        if weighted:
            return np.concatenate(
                (
                    tile_gradient,
                    [
                        8.0 * float(np.dot(np.arange(N0 + 1), probability0))
                        - 8.0 * N0 * density0,
                        8.0 * float(np.dot(np.arange(N1 + 1), probability1))
                        - 8.0 * N1 * density1,
                    ],
                )
            )
        if cumulative:
            common_q_gradient = 8.0 * (
                float(np.dot(np.arange(N0 + 1), probability0))
                + float(np.dot(np.arange(N1 + 1), probability1))
            )
            return np.concatenate(
                (
                    tile_gradient,
                    [
                        common_q_gradient - cumulative_packets / TILES
                    ],
                )
            )
        return tile_gradient

    constraint_a = np.array([a for a, _b, _p in cells], dtype=np.int64)
    constraint_b = np.array([b for _a, b, _p in cells], dtype=np.int64)
    constraint_p = np.array([p for _a, _b, p in cells], dtype=np.float64)

    def constraints(point):
        u, v, _log_q0, _log_q1 = unpack(point)
        return u[constraint_a] + v[constraint_b] - constraint_p

    constraint_jacobian = np.zeros((len(cells), N0 + N1 + 1), dtype=np.float64)
    for row, (a, b, _value) in enumerate(cells):
        if a:
            constraint_jacobian[row, a - 1] = 1.0
        constraint_jacobian[row, N0 + b] = 1.0

    def constraint_jacobian_for(point):
        if point.size == constraint_jacobian.shape[1]:
            return constraint_jacobian
        padding = point.size - constraint_jacobian.shape[1]
        return np.pad(constraint_jacobian, ((0, 0), (0, padding)))

    # u=0 and v_b=max_a log p(a,b) is feasible.
    initial_v = np.full(N1 + 1, -100.0)
    for _a, b, value in cells:
        initial_v[b] = max(initial_v[b], value)
    initial = np.concatenate((np.zeros(N0), initial_v))
    if weighted:
        def logit(value):
            clipped = min(max(value, 1e-8), 1.0 - 1e-8)
            return math.log(clipped / (1.0 - clipped))

        initial = np.concatenate((initial, [logit(density0), logit(density1)]))
    elif cumulative:
        projected_variables = TILES * 8 * (N0 + N1)
        initial_density = min(cumulative_packets / projected_variables, 0.5)
        initial_log_q = math.log(max(initial_density, 1e-12))
        initial = np.concatenate((initial, [initial_log_q]))
    bounds = [(-30.0, 30.0)] * len(initial)
    if cumulative:
        bounds[-1] = (-30.0, 0.0)
    result = minimize(
        objective,
        initial,
        jac=gradient,
        method="SLSQP",
        bounds=bounds,
        constraints={
            "type": "ineq",
            "fun": constraints,
            "jac": constraint_jacobian_for,
        },
        options={"maxiter": 4000, "ftol": 1e-12, "disp": False},
    )
    point = np.array(result.x, copy=True)
    raw_minimum_slack = float(np.min(constraints(point)))
    feasibility_repair = max(0.0, -raw_minimum_slack) + 1e-12
    if feasibility_repair:
        point[N0 : N0 + N1 + 1] += feasibility_repair
    u, v, log_q0, log_q1 = unpack(point)
    minimum_slack = float(np.min(constraints(point)))
    reported_objective = objective(point)
    if cumulative:
        score0 = log_binomial0 + np.arange(N0 + 1) * log_q0 + 8.0 * u
        score1 = log_binomial1 + np.arange(N1 + 1) * log_q1 + 8.0 * v
        full_log = TILES * (
            8.0 * logsumexp(score0) + 8.0 * logsumexp(score1)
        )
        nonzero_log = (
            full_log + math.log1p(-math.exp(-full_log))
            if full_log > 50.0
            else math.log(math.expm1(full_log))
        )
        reported_objective = nonzero_log - cumulative_packets * log_q0
    return (
        result,
        u,
        v,
        log_q0,
        log_q1,
        raw_minimum_slack,
        feasibility_repair,
        minimum_slack,
        reported_objective,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spectrum",
        type=Path,
        default=Path(__file__).resolve().parent.parent
        / "out"
        / "ebch85_band01_split_spectrum.csv",
    )
    parser.add_argument("--density0", type=float)
    parser.add_argument("--density1", type=float)
    parser.add_argument(
        "--cumulative-packets",
        type=int,
        help="bound nonzero solutions whose first-two-band support is at most this",
    )
    args = parser.parse_args()
    if (args.density0 is None) != (args.density1 is None):
        raise SystemExit("constant structured: supply both densities or neither")
    if args.cumulative_packets is not None and args.density0 is not None:
        raise SystemExit("constant structured: cumulative and exact-density modes are exclusive")
    if args.cumulative_packets is not None and args.cumulative_packets <= 0:
        raise SystemExit("constant structured: cumulative packet count must be positive")
    if args.density0 is not None and not (
        0.0 < args.density0 < 1.0 and 0.0 < args.density1 < 1.0
    ):
        raise SystemExit("constant structured: densities must lie in (0,1)")
    cells = load_cells(args.spectrum)
    (
        result,
        u,
        v,
        log_q0,
        log_q1,
        raw_minimum_slack,
        feasibility_repair,
        minimum_slack,
        tile_log,
    ) = optimize_factors(cells, args.density0, args.density1, args.cumulative_packets)
    full_log2 = (
        tile_log / math.log(2.0)
        if args.cumulative_packets is not None
        else TILES * tile_log / math.log(2.0)
    )
    print("constant-packet structured two-band upper probe")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"raw_minimum_constraint_slack={raw_minimum_slack:.12e}")
    print(f"uniform_v_feasibility_repair={feasibility_repair:.12e}")
    print(f"minimum_constraint_slack={minimum_slack:.12e}")
    print(f"one_tile_pair_log2={tile_log/math.log(2.0):.12f}")
    print(f"full_expected_solution_count_log2_upper={full_log2:.12f}")
    if args.density0 is not None:
        print(f"density0={args.density0:.12f} density1={args.density1:.12f}")
        print(f"log_q0={log_q0:.12f} log_q1={log_q1:.12f}")
    if args.cumulative_packets is not None:
        print(f"cumulative_projected_support={args.cumulative_packets}")
        print(f"log_q={log_q0:.12f}")
    print(f"u_range={float(np.min(u)):.9f},{float(np.max(u)):.9f}")
    print(f"v_range={float(np.min(v)):.9f},{float(np.max(v)):.9f}")
    print("model=ACTUAL_BALANCED_COLUMN_PARTITIONS_WITH_DEGREE8_FINNER")
    print("status=DIAGNOSTIC_BINARY64_STRUCTURED_CONSTANT_PACKET_UPPER")


if __name__ == "__main__":
    main()
