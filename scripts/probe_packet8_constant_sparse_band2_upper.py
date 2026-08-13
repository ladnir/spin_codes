#!/usr/bin/env python3
"""Sparse constant-packet bound with a determined-band-two weight tilt.

The first 85 EBCH coordinates determine the full codeword.  For first-two
band weights ``a,b``, let ``p(a,b)`` be the exact projected-code membership
probability and let ``p0(a,b)`` count the subcode whose band-two word is zero.
For ``0 < z <= 1`` the exact band-two weight moment satisfies

    sum_c p(a,b,c) z^c
      <= p0(a,b) + z^max(1,22-a-b) (p(a,b) - p0(a,b)).

The exponent uses the exact EBCH minimum distance 22: any continuation not
in the band-two-zero subcode must have positive band-two weight, and its total
weight is at least 22.

Factor this upper moment as ``f0(a) f1(b)`` and apply degree-8 Finner to the
actual balanced packet variables in bands zero and one.  A full constant word
with at most ``h`` active packets obeys the joint weight budget
``8 h_01 + c_2 <= 8h``.  A common bit fugacity ``t`` therefore uses
``q=t^8`` for first-two-band packets and ``z=t`` for determined band-two
bits.  Removing the exact zero assignment gives the cumulative nonzero
first-moment bound.

The combinatorial reduction and input spectra are exact.  SLSQP and the
reported witness evaluation use binary64; export and outward-check the final
factor witness before theorem use.
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


def load_counts(path: Path) -> dict[tuple[int, int], int]:
    result = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            result[(int(row["band0_weight"]), int(row["band1_weight"]))] = int(
                row["count"]
            )
    return result


def optimize(full_counts, zero_counts, support_cap: int):
    cells = []
    for (a, b), count in full_counts.items():
        zero = zero_counts.get((a, b), 0)
        if not 0 <= zero <= count:
            raise SystemExit("sparse band2: zero-subcode count exceeds full count")
        denominator_log = math.log(math.comb(N0, a)) + math.log(math.comb(N1, b))
        cells.append((a, b, count, zero, denominator_log))

    # point = u[1..42], v[0..43], log(q), log(z), with u[0]=0.
    off_v = N0
    off_q = off_v + N1 + 1
    off_z = off_q + 1
    variables = off_z + 1
    log_binomial0 = np.array([math.log(math.comb(N0, a)) for a in range(N0 + 1)])
    log_binomial1 = np.array([math.log(math.comb(N1, b)) for b in range(N1 + 1)])
    weights0 = np.arange(N0 + 1, dtype=np.float64)
    weights1 = np.arange(N1 + 1, dtype=np.float64)

    def unpack(point):
        u = np.concatenate(([0.0], point[:N0]))
        v = point[off_v:off_q]
        return u, v, float(point[off_q]), float(point[off_z])

    def objective(point):
        u, v, log_q, log_z = unpack(point)
        score0 = log_binomial0 + weights0 * log_q + 8.0 * u
        score1 = log_binomial1 + weights1 * log_q + 8.0 * v
        return (
            8.0 * (logsumexp(score0) + logsumexp(score1))
            - (8.0 * support_cap / TILES) * log_z
        )

    def gradient(point):
        u, v, log_q, _log_z = unpack(point)
        score0 = log_binomial0 + weights0 * log_q + 8.0 * u
        score1 = log_binomial1 + weights1 * log_q + 8.0 * v
        probability0 = np.exp(score0 - logsumexp(score0))
        probability1 = np.exp(score1 - logsumexp(score1))
        q_gradient = 8.0 * (
            float(np.dot(weights0, probability0))
            + float(np.dot(weights1, probability1))
        )
        return np.concatenate(
            (
                64.0 * probability0[1:],
                64.0 * probability1,
                [q_gradient, -8.0 * support_cap / TILES],
            )
        )

    def mixed_log_and_derivative(
        a: int,
        b: int,
        count: int,
        zero: int,
        denominator_log: float,
        log_z: float,
    ):
        other = count - zero
        minimum_c = max(1, 22 - a - b)
        if zero == 0:
            return (
                math.log(other) + minimum_c * log_z - denominator_log,
                float(minimum_c),
            )
        if other == 0:
            return math.log(zero) - denominator_log, 0.0
        left = math.log(zero)
        right = math.log(other) + minimum_c * log_z
        value = np.logaddexp(left, right)
        derivative = minimum_c * math.exp(right - value)
        return float(value - denominator_log), derivative

    def constraints(point):
        u, v, _log_q, log_z = unpack(point)
        return np.array(
            [
                u[a]
                + v[b]
                - mixed_log_and_derivative(a, b, count, zero, denominator, log_z)[0]
                for a, b, count, zero, denominator in cells
            ],
            dtype=np.float64,
        )

    def constraint_jacobian(point):
        _u, _v, _log_q, log_z = unpack(point)
        result = np.zeros((len(cells), variables), dtype=np.float64)
        for row, (a, b, count, zero, denominator) in enumerate(cells):
            if a:
                result[row, a - 1] = 1.0
            result[row, off_v + b] = 1.0
            result[row, off_z] = -mixed_log_and_derivative(
                a, b, count, zero, denominator, log_z
            )[1]
        return result

    initial_log_z = -0.5
    initial_log_q = 8.0 * initial_log_z
    initial_v = np.full(N1 + 1, -100.0)
    for a, b, count, zero, denominator in cells:
        value, _derivative = mixed_log_and_derivative(
            a, b, count, zero, denominator, initial_log_z
        )
        initial_v[b] = max(initial_v[b], value)
    initial = np.concatenate(
        (np.zeros(N0), initial_v, [initial_log_q, initial_log_z])
    )
    bounds = [(-50.0, 50.0)] * variables
    bounds[off_q] = (-40.0, 0.0)
    bounds[off_z] = (-5.0, 0.0)

    def joint_fugacity(point):
        return np.array([point[off_q] - 8.0 * point[off_z]])

    joint_jacobian = np.zeros((1, variables), dtype=np.float64)
    joint_jacobian[0, off_q] = 1.0
    joint_jacobian[0, off_z] = -8.0
    result = minimize(
        objective,
        initial,
        jac=gradient,
        method="SLSQP",
        bounds=bounds,
        constraints=[
            {
                "type": "ineq",
                "fun": constraints,
                "jac": constraint_jacobian,
            },
            {
                "type": "eq",
                "fun": joint_fugacity,
                "jac": lambda _point: joint_jacobian,
            },
        ],
        options={"maxiter": 8000, "ftol": 1e-12, "disp": False},
    )

    point = np.array(result.x, copy=True)
    raw_slack = float(np.min(constraints(point)))
    repair = max(0.0, -raw_slack) + 1e-12
    point[off_v:off_q] += repair
    slack = float(np.min(constraints(point)))
    u, v, log_q, log_z = unpack(point)
    score0 = log_binomial0 + weights0 * log_q + 8.0 * u
    score1 = log_binomial1 + weights1 * log_q + 8.0 * v
    full_log = TILES * 8.0 * (logsumexp(score0) + logsumexp(score1))
    nonzero_log = (
        full_log + math.log1p(-math.exp(-full_log))
        if full_log > 50.0
        else math.log(math.expm1(full_log))
    )
    cumulative_log = nonzero_log - 8.0 * support_cap * log_z
    return result, u, v, log_q, log_z, raw_slack, repair, slack, cumulative_log


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument("--support-cap", type=int, required=True)
    parser.add_argument(
        "--full-spectrum",
        type=Path,
        default=root / "out" / "ebch85_band01_split_spectrum.csv",
    )
    parser.add_argument(
        "--band2-zero-spectrum",
        type=Path,
        default=root / "out" / "ebch_band2_zero_band01_split_spectrum.csv",
    )
    args = parser.parse_args()
    if args.support_cap <= 0:
        raise SystemExit("sparse band2: support cap must be positive")
    result = optimize(
        load_counts(args.full_spectrum),
        load_counts(args.band2_zero_spectrum),
        args.support_cap,
    )
    optimizer, u, v, log_q, log_z, raw_slack, repair, slack, cumulative_log = result
    print("constant-packet sparse determined-band-two upper probe")
    print(f"support_cap={args.support_cap}")
    print(f"optimizer_success={optimizer.success} message={optimizer.message}")
    print(f"raw_minimum_constraint_slack={raw_slack:.12e}")
    print(f"uniform_v_feasibility_repair={repair:.12e}")
    print(f"minimum_constraint_slack={slack:.12e}")
    print(f"cumulative_nonzero_expected_count_log2_upper={cumulative_log/math.log(2.0):.12f}")
    print(f"log_q={log_q:.12f} log_z={log_z:.12f}")
    print(f"u_range={float(np.min(u)):.9f},{float(np.max(u)):.9f}")
    print(f"v_range={float(np.min(v)):.9f},{float(np.max(v)):.9f}")
    print("model=ACTUAL_TWO_BAND_BALANCED_PARTITIONS_WITH_DETERMINED_BAND2_TILT")
    print("status=DIAGNOSTIC_BINARY64_SPARSE_CONSTANT_PACKET_UPPER")


if __name__ == "__main__":
    main()
