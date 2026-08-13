#!/usr/bin/env python3
"""Degree-8 constant-packet bound with the exact 128 punctured rows.

The exact-length layout punctures one band-zero coordinate from 128 distinct
data blocks.  Use separate rank-one factors for the ordinary 42/43 projected
EBCH constraint and the exact averaged 41/43 punctured constraint.  Among the
16384 data-block factors, 16256 are ordinary and 128 are punctured.  Every
final packet variable occurs in at most eight row factors, so degree-8 Finner
remains valid even for the variables containing a graph replacement.

For a nonzero fixed data message, the sampled 24 x K graph matrix makes its
24-bit graph input uniform.  A fixed final constant-packet assignment uniquely
specifies all 128 required graph replacement bits; at most one of the 2^24
graph inputs encodes to that vector.  The final output therefore includes an
exact ``-24`` graph charge after the data-family count.

The spectra and combinatorial reduction are exact.  The factor witness and
objective evaluation are binary64 diagnostics until exported and verified
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


FULL_N0 = 42
PUNCTURED_N0 = 41
N1 = 43
BLOCKS = 16384
PUNCTURED_BLOCKS = 128
FINNER_DEGREE = 8
FULL_FACTOR_WEIGHT = (BLOCKS - PUNCTURED_BLOCKS) / FINNER_DEGREE
PUNCTURED_FACTOR_WEIGHT = PUNCTURED_BLOCKS / FINNER_DEGREE
SCALE = BLOCKS / FINNER_DEGREE
GRAPH_BITS = 24


def load_full(path: Path):
    cells = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            a = int(row["band0_weight"])
            b = int(row["band1_weight"])
            count = int(row["count"])
            cells.append(
                (
                    a,
                    b,
                    math.log(count)
                    - math.log(math.comb(FULL_N0, a))
                    - math.log(math.comb(N1, b)),
                )
            )
    return cells


def load_punctured(path: Path):
    cells = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cells.append(
                (
                    int(row["band0_weight"]),
                    int(row["band1_weight"]),
                    math.log(int(row["pair_count"]))
                    - math.log(int(row["denominator"])),
                )
            )
    return cells


def optimize(full_cells, punctured_cells, support_cap: int):
    # point = f0[1..42], f1[0..43], g0[1..41], g1[0..43], log(q).
    off_f1 = FULL_N0
    off_g0 = off_f1 + N1 + 1
    off_g1 = off_g0 + PUNCTURED_N0
    off_q = off_g1 + N1 + 1
    variables = off_q + 1

    log_binomial_f0 = np.array(
        [math.log(math.comb(FULL_N0, a)) for a in range(FULL_N0 + 1)]
    )
    log_binomial_g0 = np.array(
        [math.log(math.comb(PUNCTURED_N0, a)) for a in range(PUNCTURED_N0 + 1)]
    )
    log_binomial_1 = np.array(
        [math.log(math.comb(N1, b)) for b in range(N1 + 1)]
    )
    weights_f0 = np.arange(FULL_N0 + 1, dtype=np.float64)
    weights_g0 = np.arange(PUNCTURED_N0 + 1, dtype=np.float64)
    weights1 = np.arange(N1 + 1, dtype=np.float64)

    def unpack(point):
        f0 = np.concatenate(([0.0], point[:off_f1]))
        f1 = point[off_f1:off_g0]
        g0 = np.concatenate(([0.0], point[off_g0:off_g1]))
        g1 = point[off_g1:off_q]
        return f0, f1, g0, g1, float(point[off_q])

    def scores(point):
        f0, f1, g0, g1, log_q = unpack(point)
        return (
            log_binomial_f0 + weights_f0 * log_q + 8.0 * f0,
            log_binomial_1 + weights1 * log_q + 8.0 * f1,
            log_binomial_g0 + weights_g0 * log_q + 8.0 * g0,
            log_binomial_1 + weights1 * log_q + 8.0 * g1,
        )

    full_fraction = FULL_FACTOR_WEIGHT / SCALE
    punctured_fraction = PUNCTURED_FACTOR_WEIGHT / SCALE

    def objective(point):
        sf0, sf1, sg0, sg1 = scores(point)
        return (
            full_fraction * (logsumexp(sf0) + logsumexp(sf1))
            + punctured_fraction * (logsumexp(sg0) + logsumexp(sg1))
            - (support_cap / SCALE) * float(point[off_q])
        )

    def gradient(point):
        sf0, sf1, sg0, sg1 = scores(point)
        pf0 = np.exp(sf0 - logsumexp(sf0))
        pf1 = np.exp(sf1 - logsumexp(sf1))
        pg0 = np.exp(sg0 - logsumexp(sg0))
        pg1 = np.exp(sg1 - logsumexp(sg1))
        q_gradient = (
            full_fraction
            * (float(np.dot(weights_f0, pf0)) + float(np.dot(weights1, pf1)))
            + punctured_fraction
            * (float(np.dot(weights_g0, pg0)) + float(np.dot(weights1, pg1)))
            - support_cap / SCALE
        )
        return np.concatenate(
            (
                8.0 * full_fraction * pf0[1:],
                8.0 * full_fraction * pf1,
                8.0 * punctured_fraction * pg0[1:],
                8.0 * punctured_fraction * pg1,
                [q_gradient],
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
    projected_variables = 8 * 256 * (FULL_N0 + N1)
    initial_q = math.log(max(support_cap / projected_variables, 1e-12))
    initial = np.concatenate(
        (
            np.zeros(FULL_N0),
            initial_f1,
            np.zeros(PUNCTURED_N0),
            initial_g1,
            [initial_q],
        )
    )
    bounds = [(-40.0, 40.0)] * variables
    bounds[off_q] = (-40.0, 0.0)
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
    raw_full_slack = float(np.min(constraints(point)[: len(full_cells)]))
    raw_punctured_slack = float(np.min(constraints(point)[len(full_cells) :]))
    repair_full = max(0.0, -raw_full_slack) + 1e-12
    repair_punctured = max(0.0, -raw_punctured_slack) + 1e-12
    point[off_f1:off_g0] += repair_full
    point[off_g1:off_q] += repair_punctured
    minimum_slack = float(np.min(constraints(point)))

    sf0, sf1, sg0, sg1 = scores(point)
    log_q = float(point[off_q])
    full_log = (
        FULL_FACTOR_WEIGHT * (logsumexp(sf0) + logsumexp(sf1))
        + PUNCTURED_FACTOR_WEIGHT * (logsumexp(sg0) + logsumexp(sg1))
    )
    nonzero_log = (
        full_log + math.log1p(-math.exp(-full_log))
        if full_log > 50.0
        else math.log(math.expm1(full_log))
    )
    cumulative_log = nonzero_log - support_cap * log_q
    return (
        result,
        unpack(point),
        raw_full_slack,
        raw_punctured_slack,
        repair_full,
        repair_punctured,
        minimum_slack,
        cumulative_log,
    )


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
        "--punctured-spectrum",
        type=Path,
        default=root / "out" / "ebch84_punctured_band01_split_spectrum.csv",
    )
    args = parser.parse_args()
    if args.support_cap <= 0:
        raise SystemExit("punctured structured: support cap must be positive")
    result = optimize(
        load_full(args.full_spectrum),
        load_punctured(args.punctured_spectrum),
        args.support_cap,
    )
    (
        optimizer,
        factors,
        raw_full,
        raw_punctured,
        repair_full,
        repair_punctured,
        slack,
        cumulative_log,
    ) = result
    f0, f1, g0, g1, log_q = factors
    data_log2 = cumulative_log / math.log(2.0)
    print("constant-packet exact-length punctured structured upper probe")
    print(f"support_cap={args.support_cap}")
    print(f"optimizer_success={optimizer.success} message={optimizer.message}")
    print(f"raw_full_minimum_slack={raw_full:.12e}")
    print(f"raw_punctured_minimum_slack={raw_punctured:.12e}")
    print(f"full_uniform_repair={repair_full:.12e}")
    print(f"punctured_uniform_repair={repair_punctured:.12e}")
    print(f"minimum_constraint_slack={slack:.12e}")
    print(f"cumulative_nonzero_data_count_log2_upper={data_log2:.12f}")
    print(f"graph_match_log2_upper={-GRAPH_BITS}")
    print(f"cumulative_outer_with_graph_log2_upper={data_log2-GRAPH_BITS:.12f}")
    print(f"log_q={log_q:.12f}")
    print(f"f0_range={float(np.min(f0)):.9f},{float(np.max(f0)):.9f}")
    print(f"f1_range={float(np.min(f1)):.9f},{float(np.max(f1)):.9f}")
    print(f"g0_range={float(np.min(g0)):.9f},{float(np.max(g0)):.9f}")
    print(f"g1_range={float(np.min(g1)):.9f},{float(np.max(g1)):.9f}")
    print("model=ACTUAL_128_PUNCTURED_ROWS_DEGREE8_FINNER_PLUS_GRAPH_MATCH")
    print("status=DIAGNOSTIC_BINARY64_EXACT_LENGTH_CONSTANT_PACKET_UPPER")


if __name__ == "__main__":
    main()
