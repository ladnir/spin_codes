#!/usr/bin/env python3
"""Decompose the fixed-type loss on the outer-weight-16 shell at D=8."""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_riffle_small_exact_vs_type_bound import (  # noqa: E402
    histogram_sequence_count,
    log2_fraction,
    log_fraction,
    optimize_coefficient_bound,
    outer_histogram_distribution,
    target_accumulator_enumerators,
)
from probe_riffle_small_clustered_tilts import optimize_rows  # noqa: E402
from probe_riffle_small_typewise_bound import DEFAULT_RECEIPTS  # noqa: E402


def logsumexp_or_negative_infinity(values: list[float]) -> float:
    return float(logsumexp(values)) if values else -math.inf


def finite_or_none(value: float) -> float | None:
    """Keep receipt JSON standards-compliant when an exact mass is zero."""
    return value if math.isfinite(value) else None


def boundary_log_mass(log_x1234: np.ndarray, packet_positions: int) -> float:
    """Return log(e0^T L(x)^N e0) in the gauge x0=1."""
    log_x = np.array([0.0, *log_x1234], dtype=np.float64)
    matrix = np.zeros((5, 5), dtype=np.float64)
    for old in range(5):
        for new in range(5 - old):
            packet = old + new
            matrix[old, new] = math.comb(4 - old, new) * math.exp(log_x[packet])

    vector = np.zeros(5, dtype=np.float64)
    vector[0] = 1.0
    log_scale = 0.0
    for _ in range(packet_positions):
        vector = vector @ matrix
        scale = float(np.max(vector))
        if not scale > 0.0:
            return math.inf
        vector /= scale
        log_scale += math.log(scale)
    return log_scale + math.log(float(vector[0]))


def optimize_boundary_coefficient(
    histogram: tuple[int, ...],
    packet_positions: int,
    cached_parameters: np.ndarray,
) -> tuple[float, np.ndarray, bool]:
    """Optimize the boundary coefficient bound after removing two gauges."""

    def objective(parameters: np.ndarray) -> float:
        log_x = np.array([0.0, *parameters], dtype=np.float64)
        return boundary_log_mass(parameters, packet_positions) - float(
            np.dot(np.asarray(histogram, dtype=np.float64), log_x)
        )

    # In the z->0 boundary limit, x'_k=x_k z^(k/2).  Transform the cached
    # full-bound parameters into these four nonredundant boundary variables.
    log_z = -math.exp(float(cached_parameters[0]))
    boundary_start = np.array(
        [cached_parameters[k] + 0.5 * k * log_z for k in range(1, 5)],
        dtype=np.float64,
    )
    starts = [np.zeros(4, dtype=np.float64), boundary_start]
    results = [
        minimize(
            objective,
            start,
            method="Nelder-Mead",
            bounds=[(-120.0, 120.0)] * 4,
            options={"maxiter": 6000, "xatol": 1e-11, "fatol": 1e-11},
        )
        for start in starts
    ]
    best = min(results, key=lambda result: float(result.fun))
    return float(best.fun), np.asarray(best.x), bool(best.success)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=4)
    parser.add_argument("--parity-symbols", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--distance", type=int, default=8)
    parser.add_argument("--input-weight", type=int, default=16)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    clustered_path = DEFAULT_RECEIPTS / (
        f"goal20_clustered_b{args.data_blocks}_p{args.parity_symbols}"
        f"_d{args.distance}.json"
    )
    clustered = json.loads(clustered_path.read_text(encoding="utf-8"))
    cached_by_histogram = {
        tuple(row["histogram_h0_through_h4"]): np.asarray(
            row["parameters"], dtype=float
        )
        for row in clustered["typewise_parameters"]
    }

    outer, _signatures = outer_histogram_distribution(
        args.data_blocks, args.parity_symbols
    )
    targets = {
        histogram
        for histogram in outer
        if sum(index * count for index, count in enumerate(histogram))
        == args.input_weight
    }
    if not targets:
        raise RuntimeError("the requested outer input shell is empty")
    packet_positions = 2 * (args.data_blocks + args.parity_symbols)
    accumulator = target_accumulator_enumerators(packet_positions, targets, None)

    rows = tuple(
        (
            histogram,
            log_fraction(outer[histogram]),
            math.log(histogram_sequence_count(histogram)),
        )
        for histogram in sorted(targets)
    )
    _common_log, common_parameters, common_success = optimize_coefficient_bound(
        args.distance, packet_positions, rows, None
    )

    exact_shell = Fraction(0)
    multivariate_logs = []
    positive_multivariate_logs = []
    impossible_multivariate_logs = []
    boundary_logs = []
    type_rows = []
    optimizer_failures = 0
    boundary_optimizer_failures = 0
    for histogram, log_outer, _log_total in rows:
        total_sequences = histogram_sequence_count(histogram)
        enumerator = accumulator[histogram]
        if sum(enumerator.values()) != total_sequences:
            raise RuntimeError("the complete conditional enumerator has the wrong mass")
        minimum_output_weight = min(enumerator)
        exact_bad_sequences = sum(
            count for weight, count in enumerator.items() if weight <= args.distance
        )
        exact_probability = Fraction(exact_bad_sequences, total_sequences)
        exact_contribution = outer[histogram] * exact_probability
        exact_shell += exact_contribution

        cached_parameters = cached_by_histogram[histogram]
        bound_log_without_outer, parameters, success = optimize_rows(
            args.distance,
            packet_positions,
            ((histogram, 0.0, math.log(total_sequences)),),
            [cached_parameters, common_parameters],
        )
        optimizer_failures += int(not success)
        multivariate_log = log_outer + bound_log_without_outer
        multivariate_logs.append(multivariate_log)
        if exact_bad_sequences:
            positive_multivariate_logs.append(multivariate_log)
        else:
            impossible_multivariate_logs.append(multivariate_log)

        boundary_log_without_outer, boundary_parameters, boundary_success = (
            optimize_boundary_coefficient(
                histogram, packet_positions, cached_parameters
            )
        )
        boundary_log = log_outer - math.log(total_sequences) + boundary_log_without_outer
        boundary_logs.append(boundary_log)
        boundary_optimizer_failures += int(not boundary_success)

        # For H=2D, the deterministic H<=2W inequality makes W=D the only
        # possible bad output.  The scalar Chernoff infimum as z->0 is exact.
        scalar_chernoff_is_exact = args.input_weight == 2 * args.distance
        if scalar_chernoff_is_exact:
            if exact_bad_sequences and minimum_output_weight != args.distance:
                raise RuntimeError("the boundary-shell minimum is not D")
            scalar_log2 = log2_fraction(exact_contribution)
        else:
            scalar_log2 = None

        type_rows.append(
            {
                "histogram_h0_through_h4": list(histogram),
                "minimum_conditional_output_weight": minimum_output_weight,
                "exact_bad_sequences": str(exact_bad_sequences),
                "total_conditional_sequences": str(total_sequences),
                "exact_expected_contribution_log2": finite_or_none(
                    log2_fraction(exact_contribution)
                ),
                "scalar_chernoff_expected_contribution_log2": (
                    finite_or_none(scalar_log2)
                    if scalar_log2 is not None
                    else None
                ),
                "multivariate_bound_expected_contribution_log2": (
                    multivariate_log / math.log(2.0)
                ),
                "multivariate_loss_over_exact_bits": (
                    multivariate_log / math.log(2.0)
                    - log2_fraction(exact_contribution)
                    if exact_contribution
                    else None
                ),
                "multivariate_parameters": [float(value) for value in parameters],
                "optimizer_success": success,
                "boundary_coefficient_bound_expected_contribution_log2": (
                    boundary_log / math.log(2.0)
                ),
                "boundary_coefficient_loss_over_exact_bits": (
                    boundary_log / math.log(2.0)
                    - log2_fraction(exact_contribution)
                ),
                "boundary_parameters_log_x1_x2_x3_x4": [
                    float(value) for value in boundary_parameters
                ],
                "boundary_optimizer_success": boundary_success,
            }
        )

    exact_log2 = log2_fraction(exact_shell)
    multivariate_log = logsumexp_or_negative_infinity(multivariate_logs)
    positive_multivariate_log = logsumexp_or_negative_infinity(
        positive_multivariate_logs
    )
    impossible_multivariate_log = logsumexp_or_negative_infinity(
        impossible_multivariate_logs
    )
    boundary_log = logsumexp_or_negative_infinity(boundary_logs)
    payload = {
        "schema": "riffle-small-one-type-loss-decomposition-v1",
        "evidence_label": "EXACT_CONDITIONAL_ENUMERATOR_AND_NUMERICAL_COEFFICIENT_BOUND",
        "data_blocks": args.data_blocks,
        "parity_symbols": args.parity_symbols,
        "binary_output_length": 4 * packet_positions,
        "distance": args.distance,
        "outer_input_weight": args.input_weight,
        "packet_types": len(rows),
        "types_with_nonzero_exact_bad_tail": sum(
            int(row["exact_bad_sequences"] != "0") for row in type_rows
        ),
        "types_with_zero_exact_bad_tail": sum(
            int(row["exact_bad_sequences"] == "0") for row in type_rows
        ),
        "exact_shell_expected_count_log2": exact_log2,
        "exact_coefficient_plus_scalar_chernoff_log2": (
            exact_log2 if args.input_weight == 2 * args.distance else None
        ),
        "scalar_chernoff_loss_bits": (
            0.0 if args.input_weight == 2 * args.distance else None
        ),
        "multivariate_coefficient_bound_log2": multivariate_log / math.log(2.0),
        "multivariate_coefficient_loss_bits": (
            multivariate_log / math.log(2.0) - exact_log2
        ),
        "boundary_coefficient_bound_log2": boundary_log / math.log(2.0),
        "boundary_coefficient_loss_bits": (
            boundary_log / math.log(2.0) - exact_log2
        ),
        "multivariate_bound_on_exact_positive_types_log2": (
            positive_multivariate_log / math.log(2.0)
        ),
        "multivariate_bound_on_exact_impossible_types_log2": (
            finite_or_none(impossible_multivariate_log / math.log(2.0))
        ),
        "common_optimizer_success": common_success,
        "type_optimizer_failures": optimizer_failures,
        "boundary_optimizer_failures": boundary_optimizer_failures,
        "type_rows": type_rows,
        "scope": (
            "The conditional enumerators and exact shell count use integer and "
            "rational arithmetic. Scalar Chernoff is exact here because H=2D. "
            "The full and boundary multivariate coefficient values are numerical "
            "upper bounds and are not certified against floating-point error."
        ),
    }
    output = args.output or (
        DEFAULT_RECEIPTS
        / (
            f"goal21_one_type_loss_b{args.data_blocks}_p{args.parity_symbols}"
            f"_h{args.input_weight}_d{args.distance}.json"
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "packet_types": len(rows),
                "positive_types": payload["types_with_nonzero_exact_bad_tail"],
                "impossible_types": payload["types_with_zero_exact_bad_tail"],
                "exact_log2": exact_log2,
                "scalar_chernoff_log2": payload[
                    "exact_coefficient_plus_scalar_chernoff_log2"
                ],
                "multivariate_log2": payload[
                    "multivariate_coefficient_bound_log2"
                ],
                "coefficient_loss_bits": payload[
                    "multivariate_coefficient_loss_bits"
                ],
                "boundary_coefficient_log2": payload[
                    "boundary_coefficient_bound_log2"
                ],
                "boundary_coefficient_loss_bits": payload[
                    "boundary_coefficient_loss_bits"
                ],
                "optimizer_failures": optimizer_failures,
                "boundary_optimizer_failures": boundary_optimizer_failures,
                "status": "PASS",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
