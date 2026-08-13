#!/usr/bin/env python3
"""One-conditioned-row packet-group outer bound with exact graph averaging.

Fix one BCH codeword row in each 64-row tile.  Conditional on that row,
apply the established three-band linear Brascamp--Lieb inequality to the
remaining 63 independent rows.  The resulting factor at a coordinate depends
only on the fixed row's bit.  Exact band-(0,1) and band-(1,2) split spectra,
combined by Cauchy--Schwarz, bound the sum over the fixed row.

In a tile containing a graph hole, choose the punctured data row as the fixed
row.  The hole factor then depends on the graph bit, not the punctured data
bit.  The exact averaged punctured band-(0,1) spectrum handles the other 84
coordinates, and the exact graph spectrum averages the 128 hole factors.

This script uses binary64 arithmetic for discovery.  A theorem-facing witness
requires frozen parameters and an independent outward implementation.
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

from certify_three_band_exact_length import load_graph_spectrum
from packet_group_drive_stratified import feasible_profile_count
from packet_group_outer_profile import N, atom_count, group_log_moments


GROUP_BITS = 4
CLASSES = GROUP_BITS + 1
ROWS = 64
TILES = 256
HOLE_TILES = 128
NORMAL_TILES = TILES - HOLE_TILES
GRAPH_DIMENSION = 24
BAND_SIZES = (42, 43, 43)


def load_split_spectrum(
    path: Path, *, count_field: str = "count", divisor: int = 1
) -> np.ndarray:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            count = int(row[count_field])
            if count:
                rows.append(
                    (
                        int(row["band0_weight"]),
                        int(row["band1_weight"]),
                        math.log(count) - math.log(divisor),
                    )
                )
    if not rows:
        raise ValueError(f"conditioned-row outer: empty spectrum {path}")
    return np.asarray(rows, dtype=np.float64)


def split_log_enumerator(table: np.ndarray, left: float, right: float) -> float:
    return float(logsumexp(table[:, 2] + table[:, 0] * left + table[:, 1] * right))


def pair_cauchy_log_enumerator(
    spectrum01: np.ndarray,
    spectrum12: np.ndarray,
    log_ratios: tuple[float, float, float],
    theta: float,
) -> float:
    """Bound one three-band enumerator by two exact pair enumerators."""

    r0, r1, r2 = log_ratios
    left = split_log_enumerator(spectrum01, 2.0 * r0, 2.0 * theta * r1)
    right = split_log_enumerator(
        spectrum12, 2.0 * (1.0 - theta) * r1, 2.0 * r2
    )
    return 0.5 * (left + right)


def conditioned_norms(
    log_variables: np.ndarray,
    coefficient: float,
    conditioned_rows: int,
    group_bits: int = GROUP_BITS,
) -> tuple[float, ...]:
    """Return log factors for every weight among the conditioned rows."""

    moments = group_log_moments(group_bits, log_variables)
    free_rows = ROWS - conditioned_rows
    binomial = np.asarray(
        [
            math.log(math.comb(free_rows, weight)) - free_rows * math.log(2.0)
            for weight in range(free_rows + 1)
        ]
    )
    return tuple(
        float(
            coefficient
            * logsumexp(
                binomial
                + moments[shift : shift + free_rows + 1] / coefficient
            )
        )
        for shift in range(conditioned_rows + 1)
    )


def evaluate_conditioned_outer(
    profile: list[int],
    log_variables: np.ndarray,
    band1_coefficient: float,
    theta: float,
    spectrum01: np.ndarray,
    punctured01: np.ndarray,
    spectrum12: np.ndarray,
    conditioned_rows: int = 1,
    band_coefficients: tuple[float, float, float] | None = None,
    group_bits: int = GROUP_BITS,
) -> tuple[float, dict]:
    classes = group_bits + 1
    if len(profile) != classes or sum(profile) != atom_count(group_bits):
        raise ValueError("conditioned-row outer: malformed profile")
    if log_variables.shape != (classes,) or log_variables[0] != 0.0:
        raise ValueError("conditioned-row outer: malformed packet variables")
    if not 0.0 <= theta <= 1.0:
        raise ValueError("conditioned-row outer: invalid Cauchy split")
    if not 1 <= conditioned_rows <= ROWS:
        raise ValueError("conditioned-row outer: invalid conditioned-row count")

    coefficients = band_coefficients or (
        1.0 - band1_coefficient,
        band1_coefficient,
        band1_coefficient,
    )
    if (
        len(coefficients) != 3
        or any(not 0.0 < value <= 1.0 for value in coefficients)
        or any(
            coefficients[left] + coefficients[right] < 1.0
            for left, right in ((0, 1), (0, 2), (1, 2))
        )
    ):
        raise ValueError("conditioned-row outer: invalid BL coefficient vector")
    factors = [
        conditioned_norms(log_variables, coefficient, conditioned_rows, group_bits)
        for coefficient in coefficients
    ]
    parity_factors = [
        (
            max(row[shift] for shift in range(0, conditioned_rows + 1, 2)),
            max(row[shift] for shift in range(1, conditioned_rows + 1, 2)),
        )
        for row in factors
    ]
    zeros = tuple(row[0] for row in parity_factors)
    ratios = tuple(row[1] - row[0] for row in parity_factors)
    normal_enumerator = pair_cauchy_log_enumerator(
        spectrum01, spectrum12, ratios, theta
    )
    punctured_enumerator = pair_cauchy_log_enumerator(
        punctured01, spectrum12, ratios, theta
    )

    # Change variables from the conditioned codewords to their XOR parity
    # codeword and r-1 free codewords.  Together with the 64-r unconditioned
    # messages, this always leaves 63 free 64-bit messages.
    free_message_log = (ROWS - 1) * ROWS * math.log(2.0)
    normal_tile_log = (
        free_message_log
        + sum(size * value for size, value in zip(BAND_SIZES, zeros))
        + normal_enumerator
    )
    punctured_common_log = (
        free_message_log
        + (BAND_SIZES[0] - 1) * zeros[0]
        + BAND_SIZES[1] * zeros[1]
        + BAND_SIZES[2] * zeros[2]
        + punctured_enumerator
    )
    # At the puncture, the graph bit replaces one conditioned-row bit.  After
    # maximizing over the other r-1 conditioned bits, parity no longer
    # constrains this coordinate.
    hole_tile_logs = tuple(
        punctured_common_log
        + max(factors[0][shift] for shift in range(graph_bit, graph_bit + conditioned_rows))
        for graph_bit in (0, 1)
    )
    graph_terms = [
        math.log(count)
        - GRAPH_DIMENSION * math.log(2.0)
        + (HOLE_TILES - weight) * hole_tile_logs[0]
        + weight * hole_tile_logs[1]
        for weight, count in enumerate(load_graph_spectrum())
        if count
    ]
    natural = (
        NORMAL_TILES * normal_tile_log
        + logsumexp(graph_terms)
        - float(np.asarray(profile, dtype=np.float64) @ log_variables)
    )
    return natural / math.log(2.0), {
        "group_bits": group_bits,
        "band_coefficients": list(coefficients),
        "conditioned_rows": conditioned_rows,
        "conditioned_log_factors": [list(row) for row in factors],
        "conditioned_parity_log_factors": [list(row) for row in parity_factors],
        "conditioned_log_ratios": list(ratios),
        "normal_pair_enumerator_log": normal_enumerator,
        "punctured_pair_enumerator_log": punctured_enumerator,
        "normal_tile_log": normal_tile_log,
        "hole0_tile_log": hole_tile_logs[0],
        "hole1_tile_log": hole_tile_logs[1],
        "graph_average_log": float(logsumexp(graph_terms)),
    }


def optimize_conditioned_outer(
    profile: list[int],
    initial_log_variables: np.ndarray,
    initial_band1: float,
    initial_theta: float,
    spectrum01: np.ndarray,
    punctured01: np.ndarray,
    spectrum12: np.ndarray,
    conditioned_rows: int = 1,
    group_bits: int = GROUP_BITS,
) -> tuple[float, np.ndarray, float, float, object, dict]:
    def objective(point: np.ndarray) -> float:
        variables = np.concatenate(([0.0], point[:group_bits]))
        value, _details = evaluate_conditioned_outer(
            profile,
            variables,
            float(point[-2]),
            float(point[-1]),
            spectrum01,
            punctured01,
            spectrum12,
            conditioned_rows,
            None,
            group_bits,
        )
        return value

    starts = [
        np.concatenate((initial_log_variables[1:], [initial_band1, initial_theta])),
        np.concatenate((initial_log_variables[1:], [initial_band1, 0.5])),
        np.concatenate((np.zeros(group_bits), [0.6, 0.5])),
    ]
    bounds = [(-40.0, 40.0)] * group_bits + [(0.500001, 0.999), (0.0, 1.0)]
    candidates = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-9},
        )
        candidates.append((objective(result.x), result))
    value, result = min(candidates, key=lambda row: row[0])
    variables = np.concatenate(([0.0], result.x[:group_bits]))
    band1 = float(result.x[-2])
    theta = float(result.x[-1])
    checked, details = evaluate_conditioned_outer(
        profile,
        variables,
        band1,
        theta,
        spectrum01,
        punctured01,
        spectrum12,
        conditioned_rows,
        None,
        group_bits,
    )
    return checked, variables, band1, theta, result, details


def optimize_conditioned_outer_asymmetric(
    profile: list[int],
    initial_log_variables: np.ndarray,
    initial_band1: float,
    initial_band2: float,
    initial_theta: float,
    spectrum01: np.ndarray,
    punctured01: np.ndarray,
    spectrum12: np.ndarray,
    group_bits: int = GROUP_BITS,
) -> tuple[float, np.ndarray, tuple[float, float, float], float, object, dict]:
    """Optimize the boundary p0+p1=1 with p2>=p1."""

    def unpack(point: np.ndarray):
        band1 = float(point[group_bits])
        fraction = float(point[group_bits + 1])
        band2 = band1 + (0.999 - band1) * fraction
        return (
            np.concatenate(([0.0], point[:group_bits])),
            (1.0 - band1, band1, band2),
            float(point[group_bits + 2]),
        )

    def objective(point: np.ndarray) -> float:
        variables, coefficients, theta = unpack(point)
        value, _details = evaluate_conditioned_outer(
            profile,
            variables,
            coefficients[1],
            theta,
            spectrum01,
            punctured01,
            spectrum12,
            1,
            coefficients,
            group_bits,
        )
        return value

    initial_fraction = (
        0.0
        if initial_band2 <= initial_band1
        else min(1.0, (initial_band2 - initial_band1) / (0.999 - initial_band1))
    )
    starts = [
        np.concatenate(
            (
                initial_log_variables[1:],
                [initial_band1, initial_fraction, initial_theta],
            )
        ),
        np.concatenate((initial_log_variables[1:], [initial_band1, 0.95, initial_theta])),
        np.asarray([0.0] * group_bits + [0.6, 0.95, 0.5]),
    ]
    bounds = (
        [(-40.0, 40.0)] * group_bits
        + [(0.500001, 0.998), (0.0, 1.0), (0.0, 1.0)]
    )
    candidates = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-9},
        )
        candidates.append((objective(result.x), result))
    value, result = min(candidates, key=lambda row: row[0])
    variables, coefficients, theta = unpack(result.x)
    checked, details = evaluate_conditioned_outer(
        profile,
        variables,
        coefficients[1],
        theta,
        spectrum01,
        punctured01,
        spectrum12,
        1,
        coefficients,
        group_bits,
    )
    return checked, variables, coefficients, theta, result, details


def parse_vector(text: str, expected: int, cast=float) -> list:
    values = [cast(value) for value in text.split(",")]
    if len(values) != expected:
        raise ValueError(f"expected {expected} comma-separated values")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument("--group-bits", type=int, default=GROUP_BITS)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--log-variables", required=True)
    parser.add_argument("--band1", type=float, required=True)
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--conditioned-rows", type=int, default=1)
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--inner-log2", type=float)
    parser.add_argument(
        "--target-log2",
        type=float,
        help="per-profile target; defaults to -40-log2(number of feasible profiles)",
    )
    parser.add_argument(
        "--spectrum01", type=Path, default=root / "out" / "ebch85_band01_split_spectrum.csv"
    )
    parser.add_argument(
        "--punctured01",
        type=Path,
        default=root / "out" / "ebch84_punctured_band01_split_spectrum.csv",
    )
    parser.add_argument(
        "--spectrum12", type=Path, default=root / "out" / "ebch86_band12_split_spectrum.csv"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    target_log2 = (
        args.target_log2
        if args.target_log2 is not None
        else -40.0 - math.log2(feasible_profile_count(args.group_bits, N, 21))
    )

    classes = args.group_bits + 1
    profile = parse_vector(args.profile, classes, int)
    variables = np.asarray(parse_vector(args.log_variables, classes), dtype=np.float64)
    spectrum01 = load_split_spectrum(args.spectrum01)
    punctured01 = load_split_spectrum(
        args.punctured01, count_field="pair_count", divisor=42
    )
    spectrum12 = load_split_spectrum(args.spectrum12)
    if args.optimize:
        value, variables, band1, theta, result, details = optimize_conditioned_outer(
            profile,
            variables,
            args.band1,
            args.theta,
            spectrum01,
            punctured01,
            spectrum12,
            args.conditioned_rows,
            args.group_bits,
        )
        optimizer = {
            "success": bool(result.success),
            "message": str(result.message),
            "iterations": int(result.nit),
            "evaluations": int(result.nfev),
        }
    else:
        band1 = args.band1
        theta = args.theta
        value, details = evaluate_conditioned_outer(
            profile,
            variables,
            band1,
            theta,
            spectrum01,
            punctured01,
            spectrum12,
            args.conditioned_rows,
            None,
            args.group_bits,
        )
        optimizer = None
    report = {
        "status": "DIAGNOSTIC_BINARY64_PACKET_GROUP_ONE_CONDITIONED_ROW_EXACT_GRAPH_OUTER",
        "group_bits": args.group_bits,
        "profile": profile,
        "log_variables": variables.tolist(),
        "band1_coefficient": band1,
        "pair_cauchy_theta": theta,
        "conditioned_rows": args.conditioned_rows,
        "outer_log2": value,
        "optimizer": optimizer,
        "details": details,
        "sources": {
            "spectrum01": str(args.spectrum01),
            "punctured01": str(args.punctured01),
            "spectrum12": str(args.spectrum12),
        },
        "scope": "binary64 diagnostic; requires independent outward hardening",
    }
    if args.inner_log2 is not None:
        combined = value + args.inner_log2
        report.update(
            {
                "inner_log2": args.inner_log2,
                "combined_log2": combined,
                "target_log2": target_log2,
                "target_margin_bits": target_log2 - combined,
                "closes_target": bool(combined <= target_log2),
            }
        )
    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
