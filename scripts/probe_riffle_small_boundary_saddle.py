#!/usr/bin/env python3
"""Compare exact boundary coefficients with face-aware saddle estimates."""

from __future__ import annotations

import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from analyze_riffle_small_exact_vs_type_bound import (
    histogram_sequence_count,
    log2_fraction,
    outer_histogram_distribution,
    target_accumulator_enumerators,
)
from probe_riffle_small_typewise_bound import DEFAULT_RECEIPTS


PACKET_WIDTH = 4


def boundary_transition(old: int, new: int, packet: int) -> int:
    if old + new != packet or packet > PACKET_WIDTH:
        return 0
    return math.comb(PACKET_WIDTH - old, new)


def log_mass_moments(
    active: tuple[int, ...],
    parameters: np.ndarray,
    packet_positions: int,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return log mass, gradient, and Hessian for one boundary face."""
    variables = active[1:]
    dimension = len(variables)
    matrices = []
    matrix = np.zeros((5, 5), dtype=np.float64)
    parameter_by_packet = {packet: parameters[i] for i, packet in enumerate(variables)}
    for packet in active:
        component = np.zeros((5, 5), dtype=np.float64)
        factor = math.exp(float(parameter_by_packet.get(packet, 0.0)))
        for old in range(5):
            for new in range(5):
                multiplicity = boundary_transition(old, new, packet)
                if multiplicity:
                    component[old, new] = factor * multiplicity
        matrix += component
        if packet:
            matrices.append(component)

    value = np.zeros(5, dtype=np.float64)
    value[0] = 1.0
    first = np.zeros((dimension, 5), dtype=np.float64)
    second = np.zeros((dimension, dimension, 5), dtype=np.float64)
    log_scale = 0.0
    for _ in range(packet_positions):
        next_value = value @ matrix
        next_first = np.empty_like(first)
        next_second = np.empty_like(second)
        for i in range(dimension):
            next_first[i] = first[i] @ matrix + value @ matrices[i]
        for i in range(dimension):
            for j in range(dimension):
                next_second[i, j] = (
                    second[i, j] @ matrix
                    + first[i] @ matrices[j]
                    + first[j] @ matrices[i]
                    + (value @ matrices[i] if i == j else 0.0)
                )
        scale = float(np.max(next_value))
        if not scale > 0.0:
            raise RuntimeError("the boundary transfer lost all mass")
        value = next_value / scale
        first = next_first / scale
        second = next_second / scale
        log_scale += math.log(scale)

    endpoint = float(value[0])
    gradient = first[:, 0] / endpoint
    hessian = second[:, :, 0] / endpoint - np.outer(gradient, gradient)
    return log_scale + math.log(endpoint), gradient, hessian


def optimize_face(
    histogram: tuple[int, ...], packet_positions: int
) -> tuple[float, np.ndarray, np.ndarray, bool]:
    active = tuple(packet for packet, count in enumerate(histogram) if count)
    if not active or active[0] != 0:
        raise ValueError("the current probe requires a positive zero-packet count")
    target = np.asarray([histogram[packet] for packet in active[1:]], dtype=np.float64)

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        log_mass, gradient, _hessian = log_mass_moments(
            active, parameters, packet_positions
        )
        return log_mass - float(np.dot(target, parameters)), gradient - target

    result = minimize(
        objective,
        np.zeros(len(target), dtype=np.float64),
        method="L-BFGS-B",
        jac=True,
        bounds=[(-80.0, 80.0)] * len(target),
        options={"maxiter": 3000, "ftol": 1e-14, "gtol": 1e-10},
    )
    log_mass, gradient, hessian = log_mass_moments(
        active, np.asarray(result.x), packet_positions
    )
    gradient_error = float(np.max(np.abs(gradient - target)))
    success = bool(result.success) and gradient_error < 1e-6
    bound_log = log_mass - float(np.dot(target, result.x))
    return bound_log, np.asarray(result.x), hessian, success


@lru_cache(maxsize=None)
def reachable_face(
    active: tuple[int, ...], packet_positions: int
) -> tuple[tuple[int, ...], ...]:
    """Return reachable nonzero-packet counts for zero-to-zero paths."""
    variables = active[1:]
    index = {packet: i for i, packet in enumerate(variables)}
    zero = (0,) * len(variables)
    current: list[set[tuple[int, ...]]] = [set() for _ in range(5)]
    current[0].add(zero)
    for _ in range(packet_positions):
        following: list[set[tuple[int, ...]]] = [set() for _ in range(5)]
        for old in range(5):
            for counts in current[old]:
                used_nonzero = sum(counts)
                for packet in active:
                    if packet == 0 and used_nonzero == packet_positions:
                        continue
                    for new in range(5):
                        if not boundary_transition(old, new, packet):
                            continue
                        next_counts = list(counts)
                        if packet:
                            next_counts[index[packet]] += 1
                        if sum(next_counts) <= packet_positions:
                            following[new].add(tuple(next_counts))
        current = following
    return tuple(sorted(current[0]))


def lattice_index(active: tuple[int, ...], packet_positions: int) -> int:
    points = reachable_face(active, packet_positions)
    dimension = len(active) - 1
    if dimension == 0:
        return 1
    origin = np.asarray(points[0], dtype=np.int64)
    differences = [
        tuple(int(value) for value in np.asarray(point, dtype=np.int64) - origin)
        for point in points[1:]
        if point != points[0]
    ]
    smith = smith_normal_form(Matrix(differences).T, domain=ZZ)
    diagonal = [abs(int(smith[i, i])) for i in range(dimension)]
    if any(value == 0 for value in diagonal):
        raise RuntimeError("the reachable face does not have full affine rank")
    return math.prod(diagonal)


def exact_source_rows(
    data_blocks: int,
    parity_symbols: int,
    input_weight: int,
    distance: int,
) -> tuple[list[dict[str, object]], int]:
    """Build exact boundary rows when a Goal 21 receipt is unavailable."""
    outer, _signatures = outer_histogram_distribution(data_blocks, parity_symbols)
    targets = {
        histogram
        for histogram in outer
        if sum(packet * count for packet, count in enumerate(histogram))
        == input_weight
    }
    packet_positions = 2 * (data_blocks + parity_symbols)
    accumulator = target_accumulator_enumerators(
        packet_positions, targets, distance
    )
    rows = []
    for histogram in sorted(targets):
        exact_count = accumulator[histogram].get(distance, 0)
        if not exact_count:
            continue
        exact_contribution = (
            outer[histogram]
            * exact_count
            / histogram_sequence_count(histogram)
        )
        rows.append(
            {
                "histogram_h0_through_h4": list(histogram),
                "exact_bad_sequences": str(exact_count),
                "exact_expected_contribution_log2": log2_fraction(
                    exact_contribution
                ),
            }
        )
    return rows, 4 * packet_positions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=4)
    parser.add_argument("--parity-symbols", type=int, default=2)
    parser.add_argument("--distance", type=int, default=8)
    parser.add_argument("--input-weight", type=int, default=16)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    source = DEFAULT_RECEIPTS / (
        f"goal21_one_type_loss_b{args.data_blocks}_p{args.parity_symbols}"
        f"_h{args.input_weight}_d{args.distance}.json"
    )
    packet_positions = 2 * (args.data_blocks + args.parity_symbols)
    if source.exists():
        goal21 = json.loads(source.read_text(encoding="utf-8"))
        source_rows = goal21["type_rows"]
        binary_output_length = goal21["binary_output_length"]
        exact_source = str(source)
    else:
        source_rows, binary_output_length = exact_source_rows(
            args.data_blocks,
            args.parity_symbols,
            args.input_weight,
            args.distance,
        )
        exact_source = "computed_by_boundary_probe"
    rows = []
    exact_logs = []
    raw_bound_logs = []
    approximate_logs = []
    failures = 0
    for source_row in source_rows:
        histogram = tuple(source_row["histogram_h0_through_h4"])
        active = tuple(packet for packet, count in enumerate(histogram) if count)
        bound_log, parameters, covariance, success = optimize_face(
            histogram, packet_positions
        )
        failures += int(not success)
        sign, log_determinant = np.linalg.slogdet(covariance)
        if not sign > 0.0:
            raise RuntimeError("the saddle covariance is not positive definite")
        index = lattice_index(active, packet_positions)
        dimension = len(active) - 1
        gaussian_log_probability = (
            math.log(index)
            - 0.5 * dimension * math.log(2.0 * math.pi)
            - 0.5 * float(log_determinant)
        )
        exact_log2 = float(source_row["exact_expected_contribution_log2"])
        exact_coefficient = int(source_row["exact_bad_sequences"])
        face_bound_loss_bits = bound_log / math.log(2.0) - math.log2(exact_coefficient)
        approximation_error_bits = (
            face_bound_loss_bits + gaussian_log_probability / math.log(2.0)
        )
        approximate_log2 = exact_log2 + approximation_error_bits
        exact_logs.append(exact_log2 * math.log(2.0))
        raw_bound_logs.append(
            (exact_log2 + face_bound_loss_bits) * math.log(2.0)
        )
        approximate_logs.append(approximate_log2 * math.log(2.0))
        rows.append(
            {
                "histogram_h0_through_h4": list(histogram),
                "active_packet_weights": list(active),
                "saddle_dimension": dimension,
                "lattice_index": index,
                "face_coefficient_bound_loss_bits": face_bound_loss_bits,
                "gaussian_log_probability_bits": (
                    gaussian_log_probability / math.log(2.0)
                ),
                "saddle_approximation_error_bits": approximation_error_bits,
                "exact_expected_contribution_log2": exact_log2,
                "saddle_expected_contribution_log2": approximate_log2,
                "parameters": [float(value) for value in parameters],
                "covariance_log_determinant": float(log_determinant),
                "optimizer_success": success,
            }
        )

    exact_aggregate = float(np.logaddexp.reduce(exact_logs)) / math.log(2.0)
    raw_bound_aggregate = (
        float(np.logaddexp.reduce(raw_bound_logs)) / math.log(2.0)
    )
    approximate_aggregate = (
        float(np.logaddexp.reduce(approximate_logs)) / math.log(2.0)
    )
    errors = [row["saddle_approximation_error_bits"] for row in rows]
    payload = {
        "schema": "riffle-small-boundary-saddle-v1",
        "evidence_label": "EXACT_COEFFICIENTS_AND_NUMERICAL_SADDLE_APPROXIMATION",
        "data_blocks": args.data_blocks,
        "parity_symbols": args.parity_symbols,
        "binary_output_length": binary_output_length,
        "packet_positions": packet_positions,
        "distance": args.distance,
        "input_weight": args.input_weight,
        "packet_types": len(rows),
        "exact_source": exact_source,
        "optimizer_failures": failures,
        "exact_shell_expected_count_log2": exact_aggregate,
        "face_coefficient_bound_shell_log2": raw_bound_aggregate,
        "face_coefficient_bound_loss_bits": (
            raw_bound_aggregate - exact_aggregate
        ),
        "saddle_shell_expected_count_log2": approximate_aggregate,
        "aggregate_saddle_error_bits": approximate_aggregate - exact_aggregate,
        "minimum_type_error_bits": min(errors),
        "maximum_type_error_bits": max(errors),
        "type_rows": rows,
        "scope": (
            "The source coefficients are exact integers. The saddle points, "
            "covariance determinants, and Gaussian local factors are numerical "
            "approximations. The Gaussian factors are not proved bounds."
        ),
    }
    output = args.output or (
        DEFAULT_RECEIPTS
        / (
            f"goal22_boundary_saddle_b{args.data_blocks}_p{args.parity_symbols}"
            f"_h{args.input_weight}_d{args.distance}.json"
        )
    )
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "packet_types": len(rows),
                "optimizer_failures": failures,
                "exact_log2": exact_aggregate,
                "face_bound_log2": raw_bound_aggregate,
                "face_bound_loss_bits": raw_bound_aggregate - exact_aggregate,
                "saddle_log2": approximate_aggregate,
                "aggregate_error_bits": payload["aggregate_saddle_error_bits"],
                "type_error_range_bits": [min(errors), max(errors)],
                "lattice_indices": sorted({row["lattice_index"] for row in rows}),
                "status": "PASS" if not failures else "CHECK",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
