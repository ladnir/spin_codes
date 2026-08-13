#!/usr/bin/env python3
"""Exact graph-spectrum/puncture-averaged multivariate linear-BL probe.

A normal g=4 inner group has packet polynomial P^16.  Each of the 128 graph
holes lies in a distinct band-zero group and replaces one data coordinate.
Conditional on graph bit b, the affected group has polynomial P^15 P_b,
where P_b enumerates the other three bits of that atom and fixes the hole bit
to b.  The remaining 63 data bits use their exact binomial slice.  Averaging
the 128 hole-group factors with the exact graph24 weight spectrum removes the
adversarial per-coordinate replacement ratio.

This is a binary64 diagnostic formula.  Its independence/layout premises and
all transcendental operations require outward certification before theorem
use.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from certify_three_band_exact_length import load_graph_spectrum
from packet_group_outer_profile import (
    BINOMIAL_LOG_PROBABILITIES,
    K,
    atom_count,
    evaluate_linear_outer,
    group_log_moments,
    GROUPS,
)
from punctured_ebch_outer import full_spectrum


GROUP_BITS = 4
ATOMS = 16
TILES = 256
BAND_ZERO_GROUPS = TILES * 42
BAND_ONE_GROUPS = TILES * 86
HOLES = 128
GRAPH_DIMENSION = 24


def log_convolve(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.full(len(left) + len(right) - 1, -np.inf)
    for index, value in enumerate(left):
        result[index : index + len(right)] = np.logaddexp(
            result[index : index + len(right)], value + right
        )
    return result


def hole_group_log_moments(log_variables: np.ndarray, graph_bit: int) -> np.ndarray:
    atom = np.asarray(
        [math.log(math.comb(GROUP_BITS, weight)) + log_variables[weight] for weight in range(5)]
    )
    power = np.asarray([0.0])
    for _ in range(ATOMS - 1):
        power = log_convolve(power, atom)
    hole_atom = np.asarray(
        [
            math.log(math.comb(GROUP_BITS - 1, other_weight))
            + log_variables[other_weight + graph_bit]
            for other_weight in range(GROUP_BITS)
        ]
    )
    polynomial = log_convolve(power, hole_atom)
    return np.asarray(
        [
            polynomial[weight]
            - math.log(math.comb(63, weight))
            for weight in range(64)
        ]
    )


def exact_graph_linear_outer(
    profile: list[int], log_variables: np.ndarray, band1: float
) -> tuple[float, dict]:
    if len(profile) != 5 or sum(profile) != atom_count(GROUP_BITS):
        raise ValueError("exact graph linear BL: malformed g=4 profile")
    band0 = 1.0 - band1
    standard_moments = group_log_moments(GROUP_BITS, log_variables)
    standard0 = band0 * logsumexp(
        BINOMIAL_LOG_PROBABILITIES + standard_moments / band0
    )
    standard1 = band1 * logsumexp(
        BINOMIAL_LOG_PROBABILITIES + standard_moments / band1
    )
    remaining_binomial = np.asarray(
        [math.log(math.comb(63, weight)) - 63 * math.log(2.0) for weight in range(64)]
    )
    hole_norms = []
    for graph_bit in (0, 1):
        moments = hole_group_log_moments(log_variables, graph_bit)
        hole_norms.append(
            band0 * logsumexp(remaining_binomial + moments / band0)
        )
    graph_terms = []
    for weight, count in enumerate(load_graph_spectrum()):
        if count:
            graph_terms.append(
                math.log(count)
                - GRAPH_DIMENSION * math.log(2.0)
                + (HOLES - weight) * hole_norms[0]
                + weight * hole_norms[1]
            )
    natural = (
        K * math.log(2.0)
        + (BAND_ZERO_GROUPS - HOLES) * standard0
        + BAND_ONE_GROUPS * standard1
        + logsumexp(graph_terms)
        - float(np.asarray(profile, dtype=np.float64) @ log_variables)
    )
    adversarial_unpunctured = evaluate_linear_outer(
        GROUP_BITS, profile, log_variables, band1
    )
    return natural / math.log(2.0), {
        "standard_band0_norm": standard0,
        "standard_band1_norm": standard1,
        "hole0_band0_norm": hole_norms[0],
        "hole1_band0_norm": hole_norms[1],
        "unpunctured_linear_bl_log2": adversarial_unpunctured / math.log(2.0),
        "graph_spectrum_average_log": logsumexp(graph_terms),
    }


def optimize_exact_graph_linear_outer(
    profile: list[int],
    *,
    initial_log_variables: np.ndarray | None = None,
    initial_band1: float | None = None,
    log_variable_bound: float = 40.0,
) -> tuple[float, np.ndarray, float, dict]:
    """Minimize the exact graph-averaged linear-BL branch.

    The packet variable for weight zero is fixed to one, so its logarithm is
    zero.  Each returned point still defines a globally reusable affine outer
    witness; optimization affects sharpness only.
    """

    if len(profile) != GROUP_BITS + 1 or sum(profile) != atom_count(GROUP_BITS):
        raise ValueError("exact graph linear optimizer: malformed g=4 profile")
    starts: list[np.ndarray] = []
    zero = np.zeros(GROUP_BITS, dtype=np.float64)
    if initial_log_variables is not None:
        initial = np.asarray(initial_log_variables, dtype=np.float64)
        if initial.shape != (GROUP_BITS + 1,) or initial[0] != 0.0:
            raise ValueError("exact graph linear optimizer: invalid initial variables")
        base = initial[1:]
        bands = (initial_band1,) if initial_band1 is not None else (0.5, 0.75, 0.9)
        starts.extend(np.concatenate((base, [float(band)])) for band in bands)
    for band in (0.5, 0.6, 0.75, 0.9):
        starts.append(np.concatenate((zero, [band])))

    def objective(point: np.ndarray) -> float:
        variables = np.concatenate(([0.0], point[:GROUP_BITS]))
        return exact_graph_linear_outer(profile, variables, float(point[-1]))[0]

    bounds = [(-log_variable_bound, log_variable_bound)] * GROUP_BITS + [
        (0.5, 0.999)
    ]
    candidates = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-9},
        )
        value = objective(result.x)
        candidates.append((value, result))
    value, result = min(candidates, key=lambda row: row[0])
    variables = np.concatenate(([0.0], result.x[:GROUP_BITS]))
    band1 = float(result.x[-1])
    checked, details = exact_graph_linear_outer(profile, variables, band1)
    details = {
        **details,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_iterations": int(result.nit),
        "optimizer_evaluations": int(result.nfev),
    }
    return checked, variables, band1, details


def exact_graph_total_spectrum_outer(
    profile: list[int], log_variables: np.ndarray, log_beta: float
) -> tuple[float, dict]:
    """Total-spectrum branch with exact hole-group graph averaging."""

    normal_moments = group_log_moments(GROUP_BITS, log_variables)
    weights = np.arange(65, dtype=np.float64)
    log_alpha = float(np.max(normal_moments - weights * log_beta))

    atom = np.asarray(
        [math.log(math.comb(GROUP_BITS, weight)) + log_variables[weight] for weight in range(5)]
    )
    power = np.asarray([0.0])
    for _ in range(ATOMS - 1):
        power = log_convolve(power, atom)
    hole_alphas = []
    for graph_bit in (0, 1):
        hole_atom = np.full(5, -np.inf)
        for old_bit in (0, 1):
            for other_weight in range(4):
                unpunctured_weight = old_bit + other_weight
                term = (
                    math.log(math.comb(3, other_weight))
                    + log_variables[other_weight + graph_bit]
                )
                hole_atom[unpunctured_weight] = np.logaddexp(
                    hole_atom[unpunctured_weight], term
                )
        polynomial = log_convolve(power, hole_atom)
        moments = np.asarray(
            [
                polynomial[weight] - math.log(math.comb(64, weight))
                for weight in range(65)
            ]
        )
        hole_alphas.append(float(np.max(moments - weights * log_beta)))

    spectrum = full_spectrum()
    enumerator = logsumexp(
        [math.log(count) + weight * log_beta for weight, count in enumerate(spectrum) if count]
    )
    graph = logsumexp(
        [
            math.log(count)
            - GRAPH_DIMENSION * math.log(2.0)
            + (HOLES - weight) * hole_alphas[0]
            + weight * hole_alphas[1]
            for weight, count in enumerate(load_graph_spectrum())
            if count
        ]
    )
    natural = (
        (GROUPS - HOLES) * log_alpha
        + (K // 64) * enumerator
        + graph
        - float(np.asarray(profile, dtype=np.float64) @ log_variables)
    )
    return natural / math.log(2.0), {
        "normal_log_alpha": log_alpha,
        "hole0_log_alpha": hole_alphas[0],
        "hole1_log_alpha": hole_alphas[1],
        "log_enumerator": enumerator,
        "graph_spectrum_average_log": graph,
        "log_beta": log_beta,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--log-variables", required=True)
    parser.add_argument("--band1", type=float, required=True)
    parser.add_argument("--log-beta", type=float)
    parser.add_argument("--old-outer-log2", type=float)
    parser.add_argument("--inner-log2", type=float)
    parser.add_argument("--target-log2", type=float, default=-111.41506501642425)
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profile = [int(value) for value in args.profile.split(",")]
    variables = np.asarray([float(value) for value in args.log_variables.split(",")])
    if args.optimize:
        value, variables, band1, details = optimize_exact_graph_linear_outer(
            profile,
            initial_log_variables=variables,
            initial_band1=args.band1,
        )
    else:
        band1 = args.band1
        value, details = exact_graph_linear_outer(profile, variables, band1)
    report = {
        "status": "DIAGNOSTIC_BINARY64_G4_EXACT_GRAPH_LINEAR_BL",
        "profile": profile,
        "log_variables": variables.tolist(),
        "band1_coefficient": band1,
        "outer_log2": value,
        "details": details,
        "lemma_candidate": (
            "128 distinct band-zero hole groups; P^15 P_b group moments; "
            "exact graph24 spectrum average"
        ),
    }
    if args.old_outer_log2 is not None:
        report["improvement_bits"] = args.old_outer_log2 - value
    if args.inner_log2 is not None:
        combined = value + args.inner_log2
        report["combined_log2"] = combined
        report["target_margin_bits"] = args.target_log2 - combined
        report["closes_target"] = bool(combined <= args.target_log2)
    if args.log_beta is not None:
        spectrum_value, spectrum_details = exact_graph_total_spectrum_outer(
            profile, variables, args.log_beta
        )
        spectrum_report = {
            "outer_log2": spectrum_value,
            "details": spectrum_details,
        }
        if args.old_outer_log2 is not None:
            spectrum_report["improvement_bits"] = args.old_outer_log2 - spectrum_value
        if args.inner_log2 is not None:
            combined = spectrum_value + args.inner_log2
            spectrum_report["combined_log2"] = combined
            spectrum_report["target_margin_bits"] = args.target_log2 - combined
            spectrum_report["closes_target"] = bool(combined <= args.target_log2)
        report["exact_graph_total_spectrum"] = spectrum_report
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
