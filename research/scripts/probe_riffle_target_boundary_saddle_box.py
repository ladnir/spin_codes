#!/usr/bin/env python3
"""Calibrate boundary saddles at the target Riffle g=4 packet length.

This is a numerical diagnostic, not a proof certificate.  It uses the Perron
limit of the five-state boundary transfer matrix, so its cost is independent
of the target packet count.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eig
from scipy.optimize import brentq, minimize


PACKET_WIDTH = 4
DEFAULT_PACKET_COUNT = 524_352
DEFAULT_DISTANCES = (76_000, 77_000, 188_743)


def component_matrix(packet_weight: int) -> np.ndarray:
    matrix = np.zeros((5, 5), dtype=np.float64)
    for old in range(5):
        for new in range(5):
            if old + new == packet_weight and packet_weight <= PACKET_WIDTH:
                matrix[old, new] = math.comb(PACKET_WIDTH - old, new)
    return matrix


COMPONENTS = tuple(component_matrix(weight) for weight in range(5))


def perron_value_and_gradient(parameters: np.ndarray) -> tuple[float, np.ndarray]:
    matrix = COMPONENTS[0].copy()
    weighted = []
    for packet_weight, parameter in enumerate(parameters, start=1):
        component = math.exp(float(parameter)) * COMPONENTS[packet_weight]
        matrix += component
        weighted.append(component)

    values, left, right = eig(matrix, left=True, right=True)
    index = int(np.argmax(values.real))
    rho = float(values[index].real)
    if not rho > 0.0 or abs(float(values[index].imag)) > 1e-9:
        raise RuntimeError("failed to isolate the positive Perron root")
    left_vector = left[:, index]
    right_vector = right[:, index]
    denominator = np.vdot(left_vector, right_vector)
    gradient = np.asarray(
        [
            (np.vdot(left_vector, component @ right_vector) / denominator).real
            / rho
            for component in weighted
        ],
        dtype=np.float64,
    )
    return math.log(rho), gradient


def optimize_saddle(packet_fractions: np.ndarray) -> dict[str, object]:
    target = np.asarray(packet_fractions[1:], dtype=np.float64)

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        log_rho, gradient = perron_value_and_gradient(parameters)
        return log_rho - float(np.dot(target, parameters)), gradient - target

    result = minimize(
        objective,
        np.full(4, -2.0, dtype=np.float64),
        method="L-BFGS-B",
        jac=True,
        bounds=[(-40.0, 20.0)] * 4,
        options={"maxiter": 2000, "ftol": 1e-15, "gtol": 1e-12},
    )
    _log_rho, gradient = perron_value_and_gradient(np.asarray(result.x))
    error = float(np.max(np.abs(gradient - target)))
    return {
        "parameters_natural_log_x1_through_x4": [float(x) for x in result.x],
        "maximum_fraction_error": error,
        "optimizer_reported_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_success": error < 2e-8,
    }


def uniform_weight_subset_packet_fractions(
    packet_count: int, input_weight: int
) -> np.ndarray:
    """Exact one-packet hypergeometric marginal of a uniform weight subset."""
    total_bits = PACKET_WIDTH * packet_count
    probabilities = np.empty(5, dtype=np.float64)
    probabilities[0] = math.prod(
        (total_bits - input_weight - index) / (total_bits - index)
        for index in range(PACKET_WIDTH)
    )
    for weight in range(PACKET_WIDTH):
        probabilities[weight + 1] = probabilities[weight] * (
            (PACKET_WIDTH - weight)
            / (weight + 1)
            * (input_weight - weight)
            / (total_bits - input_weight - PACKET_WIDTH + weight + 1)
        )
    probabilities /= np.sum(probabilities)
    return probabilities


def fixed_coordinate_packet_fractions(
    packet_count: int,
    input_weight: int,
    fixed_weight: int,
    fixed_count: float,
) -> np.ndarray:
    """Maximum-entropy packet law after fixing one packet-weight count."""
    remaining_packets = packet_count - fixed_count
    remaining_weight = input_weight - fixed_weight * fixed_count
    target_mean = remaining_weight / remaining_packets
    remaining_weights = [weight for weight in range(5) if weight != fixed_weight]

    def mean_at(log_x: float) -> float:
        terms = np.asarray(
            [
                math.comb(4, weight) * math.exp(log_x * weight)
                for weight in remaining_weights
            ],
            dtype=np.float64,
        )
        return float(np.dot(remaining_weights, terms) / np.sum(terms))

    log_x = brentq(lambda value: mean_at(value) - target_mean, -40.0, 40.0)
    terms = np.asarray(
        [
            math.comb(4, weight) * math.exp(log_x * weight)
            for weight in remaining_weights
        ],
        dtype=np.float64,
    )
    fractions = np.zeros(5, dtype=np.float64)
    for weight, term in zip(remaining_weights, terms):
        fractions[weight] = remaining_packets / packet_count * term / np.sum(terms)
    fractions[fixed_weight] = fixed_count / packet_count
    return fractions


def distance_row(packet_count: int, distance: int) -> dict[str, object]:
    input_weight = 2 * distance
    bit_density = input_weight / (4.0 * packet_count)
    center = uniform_weight_subset_packet_fractions(packet_count, input_weight)
    expected_counts = center * packet_count
    center_saddle = optimize_saddle(center)

    # The marginal scale is sufficient for this diagnostic rare-coordinate
    # sweep.  It is not asserted to be a simultaneous confidence interval.
    weight4_mean = float(expected_counts[4])
    weight4_sd = math.sqrt(packet_count * center[4] * (1.0 - center[4]))
    candidates = {1.0, 2.0, 4.0, 8.0, float(max(1, round(weight4_mean)))}
    candidates.add(float(max(1, round(weight4_mean + 6.0 * weight4_sd))))
    sweep = []
    for count in sorted(candidates):
        fractions = fixed_coordinate_packet_fractions(
            packet_count, input_weight, 4, count
        )
        saddle = optimize_saddle(fractions)
        sweep.append(
            {
                "weight4_count": count,
                "packet_fractions_h0_through_h4": [float(x) for x in fractions],
                **saddle,
            }
        )

    one_count_sweep = []
    for fixed_weight in range(1, 5):
        fractions = fixed_coordinate_packet_fractions(
            packet_count, input_weight, fixed_weight, 1.0
        )
        one_count_sweep.append(
            {
                "fixed_packet_weight": fixed_weight,
                "fixed_count": 1.0,
                "packet_fractions_h0_through_h4": [float(x) for x in fractions],
                **optimize_saddle(fractions),
            }
        )

    return {
        "distance": distance,
        "boundary_input_weight": input_weight,
        "bit_density": bit_density,
        "uniform_weight_subset_center_counts_h0_through_h4": [
            float(x) for x in expected_counts
        ],
        "center_saddle": center_saddle,
        "weight4_marginal_standard_deviation": weight4_sd,
        "fixed_weight4_count_sweep": sweep,
        "one_count_each_nonzero_packet_weight": one_count_sweep,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-count", type=int, default=DEFAULT_PACKET_COUNT)
    parser.add_argument(
        "--distances", type=int, nargs="+", default=list(DEFAULT_DISTANCES)
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.packet_count <= 0:
        raise ValueError("packet count must be positive")
    if any(not 0 < 2 * distance < 4 * args.packet_count for distance in args.distances):
        raise ValueError("every distance must give an interior boundary weight")

    rows = [distance_row(args.packet_count, distance) for distance in args.distances]
    parameters = [
        parameter
        for row in rows
        for sweep_row in (
            row["fixed_weight4_count_sweep"]
            + row["one_count_each_nonzero_packet_weight"]
        )
        for parameter in sweep_row["parameters_natural_log_x1_through_x4"]
    ]
    payload = {
        "schema": "riffle-target-boundary-saddle-calibration-v1",
        "evidence_label": "NUMERICAL_PERRON_LIMIT_DIAGNOSTIC",
        "packet_width": PACKET_WIDTH,
        "packet_count": args.packet_count,
        "rows": rows,
        "observed_parameter_range_over_sweeps": [min(parameters), max(parameters)],
        "scope": (
            "The center is the packet histogram of a uniformly random binary "
            "word conditional on boundary input weight. It is not the proved "
            "histogram law of the double-parity outer code. Perron-limit "
            "saddles omit finite-length bridge endpoint corrections."
        ),
    }
    output = args.output or Path(
        "constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/"
        "receipts/goal26_target_boundary_saddle_calibration.json"
    )
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(output),
                "observed_parameter_range": payload[
                    "observed_parameter_range_over_sweeps"
                ],
                "all_optimizers_succeeded": all(
                    row["center_saddle"]["optimizer_success"]
                    and all(
                        sweep_row["optimizer_success"]
                        for sweep_row in (
                            row["fixed_weight4_count_sweep"]
                            + row["one_count_each_nonzero_packet_weight"]
                        )
                    )
                    for row in rows
                ),
                "centers": [
                    {
                        "distance": row["distance"],
                        "counts": row[
                            "uniform_weight_subset_center_counts_h0_through_h4"
                        ],
                        "parameters": row["center_saddle"][
                            "parameters_natural_log_x1_through_x4"
                        ],
                    }
                    for row in rows
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
