#!/usr/bin/env python3
"""g=4 packet-profile outer bound from exact EBCH band-pair spectra.

For band weights ``(a,b,c)``, a full EBCH word satisfies both projected-code
conditions on bands ``(0,1)`` and ``(1,2)``.  Exact split spectra give the
corresponding membership probabilities ``p01(a,b)`` and ``p12(b,c)``.
Positive rank-one factors satisfying

    log p01(a,b) <= u0[a] + u1[b],
    log p12(b,c) <= v1[b] + v2[c]

therefore dominate full membership by the geometric mean.  Degree-4 Finner
then separates the packet variables in each physical band.  The optimization
below is convex, subject to linear factor constraints.  Binary64 optimization
is diagnostic; a certificate needs frozen factors and outward evaluation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from packet_group_drive_stratified import profile_classes
from packet_group_outer_profile import K, N, atom_count, normalization_log2


GROUP_BITS = 4
CLASSES = GROUP_BITS + 1
BAND_SIZES = (42, 43, 43)
TILES = 256
ROWS_PER_BAND = 64 * TILES
FINNER_DEGREE = GROUP_BITS
ROW_FACTOR = ROWS_PER_BAND / FINNER_DEGREE


def load_cells(path: Path, left: int, right: int) -> list[tuple[int, int, float]]:
    cells = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            a = int(row["band0_weight"])
            b = int(row["band1_weight"])
            count = int(row["count"])
            if count <= 0:
                continue
            if not (0 <= a <= left and 0 <= b <= right):
                raise ValueError("g4 pair-spectrum outer: malformed spectrum cell")
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
        raise ValueError("g4 pair-spectrum outer: empty split spectrum")
    return cells


def optimize_profile(
    profile: list[int],
    cells01: list[tuple[int, int, float]],
    cells12: list[tuple[int, int, float]],
) -> tuple[object, np.ndarray, dict]:
    if len(profile) != CLASSES or sum(profile) != atom_count(GROUP_BITS):
        raise ValueError("g4 pair-spectrum outer: malformed packet profile")

    n0, n1, n2 = BAND_SIZES
    off_u0 = 0
    off_u1 = off_u0 + n0
    off_v1 = off_u1 + n1 + 1
    off_v2 = off_v1 + n1
    factor_variables = off_v2 + n2 + 1
    off_t = factor_variables
    variables = factor_variables + GROUP_BITS

    def unpack(point: np.ndarray):
        u0 = np.concatenate(([0.0], point[off_u0:off_u1]))
        u1 = point[off_u1:off_v1]
        v1 = np.concatenate(([0.0], point[off_v1:off_v2]))
        v2 = point[off_v2:factor_variables]
        log_t = np.concatenate(([0.0], point[off_t:]))
        return u0, u1, v1, v2, log_t

    rows = len(cells01) + len(cells12)
    constraint_matrix = np.zeros((rows, variables), dtype=np.float64)
    rhs = np.zeros(rows, dtype=np.float64)
    for row, (a, b, value) in enumerate(cells01):
        if a:
            constraint_matrix[row, off_u0 + a - 1] = 1.0
        constraint_matrix[row, off_u1 + b] = 1.0
        rhs[row] = value
    base = len(cells01)
    for local, (b, c, value) in enumerate(cells12):
        row = base + local
        if b:
            constraint_matrix[row, off_v1 + b - 1] = 1.0
        constraint_matrix[row, off_v2 + c] = 1.0
        rhs[row] = value

    log_one_coeff = np.full(CLASSES, -np.inf)
    log_zero_coeff = np.full(CLASSES, -np.inf)
    for weight in range(CLASSES):
        if weight:
            log_one_coeff[weight] = math.log(math.comb(GROUP_BITS - 1, weight - 1))
        if weight < GROUP_BITS:
            log_zero_coeff[weight] = math.log(math.comb(GROUP_BITS - 1, weight))
    log_binomials = [
        np.asarray([math.log(math.comb(size, weight)) for weight in range(size + 1)])
        for size in BAND_SIZES
    ]
    weight_vectors = [np.arange(size + 1, dtype=np.float64) for size in BAND_SIZES]
    profile_vector = np.asarray(profile, dtype=np.float64)

    def components(point: np.ndarray):
        u0, u1, v1, v2, log_t = unpack(point)
        log_a = logsumexp(log_one_coeff + log_t)
        log_c = logsumexp(log_zero_coeff + log_t)
        alpha = np.exp(log_one_coeff + log_t - log_a)
        beta = np.exp(log_zero_coeff + log_t - log_c)
        factors = (0.5 * u0, 0.5 * (u1 + v1), 0.5 * v2)
        scores = [
            log_binomial
            + weights * log_a
            + (size - weights) * log_c
            + FINNER_DEGREE * factor
            for size, weights, log_binomial, factor in zip(
                BAND_SIZES, weight_vectors, log_binomials, factors
            )
        ]
        return log_t, alpha, beta, scores

    def objective(point: np.ndarray) -> float:
        log_t, _alpha, _beta, scores = components(point)
        return ROW_FACTOR * sum(logsumexp(score) for score in scores) - float(
            profile_vector @ log_t
        )

    def gradient(point: np.ndarray) -> np.ndarray:
        u0, u1, v1, v2, _log_t = unpack(point)
        del u0, u1, v1, v2
        _log_t, alpha, beta, scores = components(point)
        result = np.zeros(variables, dtype=np.float64)
        probabilities = [np.exp(score - logsumexp(score)) for score in scores]
        # Each pair factor enters through the geometric mean and degree-4
        # Finner, hence multiplier two.
        result[off_u0:off_u1] = 2.0 * ROW_FACTOR * probabilities[0][1:]
        result[off_u1:off_v1] = 2.0 * ROW_FACTOR * probabilities[1]
        result[off_v1:off_v2] = 2.0 * ROW_FACTOR * probabilities[1][1:]
        result[off_v2:factor_variables] = 2.0 * ROW_FACTOR * probabilities[2]
        t_gradient = np.zeros(CLASSES, dtype=np.float64)
        for size, weights, probability in zip(BAND_SIZES, weight_vectors, probabilities):
            expected = float(weights @ probability)
            t_gradient += ROW_FACTOR * (
                expected * alpha + (size - expected) * beta
            )
        t_gradient -= profile_vector
        result[off_t:] = t_gradient[1:]
        return result

    def constraints(point: np.ndarray) -> np.ndarray:
        return constraint_matrix @ point - rhs

    initial_u1 = np.full(n1 + 1, -100.0)
    for _a, b, value in cells01:
        initial_u1[b] = max(initial_u1[b], value)
    initial_v2 = np.full(n2 + 1, -100.0)
    for _b, c, value in cells12:
        initial_v2[c] = max(initial_v2[c], value)
    empirical = np.full(CLASSES, -24.0)
    classes = profile_classes(GROUP_BITS)
    for weight, count in enumerate(profile):
        if count:
            empirical[weight] = math.log(count / classes[weight])
    empirical -= empirical[0]
    initial = np.concatenate(
        (np.zeros(n0), initial_u1, np.zeros(n1), initial_v2, empirical[1:])
    )
    result = minimize(
        objective,
        initial,
        jac=gradient,
        method="SLSQP",
        bounds=[(-60.0, 60.0)] * variables,
        constraints={
            "type": "ineq",
            "fun": constraints,
            "jac": lambda _point: constraint_matrix,
        },
        options={"maxiter": 8000, "ftol": 1e-10, "disp": False},
    )

    point = np.asarray(result.x, dtype=np.float64).copy()
    raw_slack = float(np.min(constraints(point)))
    repair = max(0.0, -raw_slack) + 1e-11
    point[off_u1:off_v1] += repair
    point[off_v2:factor_variables] += repair
    repaired_slack = float(np.min(constraints(point)))
    value = objective(point)
    u0, u1, v1, v2, log_t = unpack(point)
    return result, point, {
        "outer_log2": value / math.log(2.0),
        "log_variables": log_t.tolist(),
        "raw_constraint_slack": raw_slack,
        "uniform_factor_repair": repair,
        "minimum_constraint_slack": repaired_slack,
        "u0": u0.tolist(),
        "u1": u1.tolist(),
        "v1": v1.tolist(),
        "v2": v2.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument("--profile", required=True)
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
    parser.add_argument("--inner-log2", type=float)
    parser.add_argument("--target-log2", type=float, default=-111.41506501642425)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profile = [int(value) for value in args.profile.split(",")]
    cells01 = load_cells(args.spectrum01, BAND_SIZES[0], BAND_SIZES[1])
    cells12 = load_cells(args.spectrum12, BAND_SIZES[1], BAND_SIZES[2])
    result, _point, details = optimize_profile(profile, cells01, cells12)
    normalization = float(
        normalization_log2(GROUP_BITS, np.asarray(profile, dtype=np.float64))[0]
    )
    report = {
        "status": "DIAGNOSTIC_BINARY64_G4_PAIR_SPECTRUM_PROFILE_OUTER",
        "profile": profile,
        "outer_log2": details["outer_log2"],
        "normalization_log2": normalization,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_iterations": int(result.nit),
        "optimizer_evaluations": int(result.nfev),
        "spectrum01": str(args.spectrum01),
        "spectrum12": str(args.spectrum12),
        "details": details,
        "total_mass_sanity_gap_bits": details["outer_log2"] - K,
    }
    if args.inner_log2 is not None:
        combined = details["outer_log2"] + args.inner_log2
        report.update(
            {
                "inner_log2": args.inner_log2,
                "combined_log2": combined,
                "target_margin_bits": args.target_log2 - combined,
                "closes_target": bool(combined <= args.target_log2),
            }
        )
    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
