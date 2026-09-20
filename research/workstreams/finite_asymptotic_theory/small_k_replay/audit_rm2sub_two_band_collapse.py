#!/usr/bin/env python3
"""Audit the two-binomial collapse by direct finite enumeration."""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np

from probe_rm2sub_two_band_bridge import (
    two_binomial_region_log_matrices,
)


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "rm2sub_two_band_collapse_audit.json"


def direct_region_matrix(
    impulses: np.ndarray,
    total_positions: int,
    epoch_bits: int,
    first_count: int,
    second_count: int,
    first_probability: float,
    second_probability: float,
) -> np.ndarray:
    positions = range(total_positions)
    result = np.zeros((2, 2), dtype=np.float64)
    assignment_count = 0
    for first_positions in itertools.combinations(positions, first_count):
        first_set = set(first_positions)
        remaining = [position for position in positions if position not in first_set]
        for second_positions in itertools.combinations(remaining, second_count):
            assignment_count += 1
            typed = [
                *((position, first_probability) for position in first_positions),
                *((position, second_probability) for position in second_positions),
            ]
            for bits in itertools.product((0, 1), repeat=len(typed)):
                probability = 1.0
                epoch_weights = [0] * (total_positions // epoch_bits)
                for (position, live_probability), bit in zip(typed, bits):
                    probability *= (
                        live_probability if bit else 1.0 - live_probability
                    )
                    if bit:
                        epoch_weights[position // epoch_bits] += 1
                product = np.eye(2, dtype=np.float64)
                for weight in epoch_weights:
                    product = product @ impulses[weight]
                result += probability * product
    return result / assignment_count


def main() -> None:
    rng = np.random.default_rng(20260904)
    epoch_bits = 2
    epochs_per_region = 3
    total_positions = epoch_bits * epochs_per_region
    maximum_occupation = 4
    first_probability = 0.23
    second_probability = 0.71
    impulses = rng.uniform(0.02, 0.8, size=(epoch_bits + 1, 2, 2))
    with np.errstate(divide="ignore"):
        forced_epoch_logs = np.log(impulses)
    from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
        regular_region_log_matrices,
    )

    forced_region_logs = regular_region_log_matrices(
        forced_epoch_logs,
        epoch_bits,
        epochs_per_region,
        maximum_occupation,
    )
    collapsed = two_binomial_region_log_matrices(
        forced_region_logs,
        first_probability,
        second_probability,
        maximum_occupation,
    )
    common_probability = 0.37
    homogeneous = two_binomial_region_log_matrices(
        forced_region_logs,
        common_probability,
        common_probability,
        maximum_occupation,
    )
    from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
        binomial_transform,
        candidate_epoch_logs,
    )

    homogeneous_epoch_logs = candidate_epoch_logs(
        impulses, binomial_transform(epoch_bits, common_probability)
    )
    established_homogeneous = regular_region_log_matrices(
        homogeneous_epoch_logs,
        epoch_bits,
        epochs_per_region,
        maximum_occupation,
    )
    homogeneous_error = 0.0
    for first_count in range(maximum_occupation + 1):
        for second_count in range(maximum_occupation - first_count + 1):
            difference = np.abs(
                np.exp(homogeneous[first_count, second_count])
                - np.exp(established_homogeneous[first_count + second_count])
            )
            homogeneous_error = max(
                homogeneous_error, float(np.max(difference))
            )
    rows = []
    maximum_error = 0.0
    for first_count in range(maximum_occupation + 1):
        for second_count in range(maximum_occupation - first_count + 1):
            direct = direct_region_matrix(
                impulses,
                total_positions,
                epoch_bits,
                first_count,
                second_count,
                first_probability,
                second_probability,
            )
            recovered = np.exp(collapsed[first_count, second_count])
            error = float(np.max(np.abs(direct - recovered)))
            maximum_error = max(maximum_error, error)
            rows.append(
                {
                    "first_count": first_count,
                    "second_count": second_count,
                    "maximum_absolute_error": error,
                }
            )
    if maximum_error > 2e-14 or homogeneous_error > 2e-14:
        raise AssertionError(f"two-band collapse audit failed: {maximum_error}")
    payload = {
        "schema": "rm2sub-two-band-collapse-audit-v1",
        "status": "PASS",
        "parameters": {
            "epoch_bits": epoch_bits,
            "epochs_per_region": epochs_per_region,
            "total_positions": total_positions,
            "maximum_occupation": maximum_occupation,
            "first_probability": first_probability,
            "second_probability": second_probability,
        },
        "maximum_absolute_error": maximum_error,
        "homogeneous_reduction_maximum_absolute_error": homogeneous_error,
        "rows": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"maximum_absolute_error={maximum_error:.17g}")
    print(f"homogeneous_reduction_error={homogeneous_error:.17g}")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
