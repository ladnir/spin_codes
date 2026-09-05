#!/usr/bin/env python3
"""Audit the three-group collapse and consolidated middle-region receipt."""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parent.parent.parent
SCRIPTS = REPOSITORY / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    log_matrix_power_moments_batch,
    regular_region_log_matrices,
)
from probe_rm2sub_three_group_bridge import (
    GROUPS,
    compositions_three,
    log_matrix_power_moments_binary,
    three_binomial_distributions,
)


RECEIPT = HERE / "rm2sub_three_group_bridge_q30_q52_probe_d100.json"
OUTPUT = HERE / "rm2sub_three_group_bridge_audit.json"
LOG2 = math.log(2.0)


def direct_region_matrix(
    impulses: np.ndarray,
    total_positions: int,
    epoch_bits: int,
    counts: tuple[int, int, int],
    probabilities: tuple[float, float, float],
) -> np.ndarray:
    positions = tuple(range(total_positions))
    result = np.zeros((2, 2), dtype=np.float64)
    assignment_count = 0
    for first in itertools.combinations(positions, counts[0]):
        remaining_first = tuple(position for position in positions if position not in first)
        for second in itertools.combinations(remaining_first, counts[1]):
            remaining_second = tuple(
                position for position in remaining_first if position not in second
            )
            for third in itertools.combinations(remaining_second, counts[2]):
                assignment_count += 1
                typed = [
                    *((position, probabilities[0]) for position in first),
                    *((position, probabilities[1]) for position in second),
                    *((position, probabilities[2]) for position in third),
                ]
                for bits in itertools.product((0, 1), repeat=sum(counts)):
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


def close(first: float, second: float, tolerance: float = 2e-9) -> None:
    if not math.isclose(first, second, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"aggregation mismatch: {first} != {second}")


def main() -> None:
    rng = np.random.default_rng(20260904)
    epoch_bits = 2
    epochs_per_region = 3
    total_positions = epoch_bits * epochs_per_region
    maximum_occupation = 3
    probabilities = (0.21, 0.47, 0.79)
    impulses = rng.uniform(0.02, 0.8, size=(epoch_bits + 1, 2, 2))
    with np.errstate(divide="ignore"):
        forced_epoch_logs = np.log(impulses)
    forced_region_logs = regular_region_log_matrices(
        forced_epoch_logs,
        epoch_bits,
        epochs_per_region,
        maximum_occupation,
    )
    forced_region = np.exp(forced_region_logs).reshape(
        (maximum_occupation + 1, 4)
    )
    compositions = [
        composition
        for occupation in range(maximum_occupation + 1)
        for composition in compositions_three(occupation)
    ]
    distributions = three_binomial_distributions(
        compositions, maximum_occupation, probabilities
    )
    collapsed = (distributions @ forced_region).reshape((-1, 2, 2))
    direct_error = 0.0
    for composition, recovered in zip(compositions, collapsed):
        direct = direct_region_matrix(
            impulses,
            total_positions,
            epoch_bits,
            composition,
            probabilities,
        )
        direct_error = max(
            direct_error, float(np.max(np.abs(direct - recovered)))
        )
    if direct_error > 2e-14:
        raise AssertionError(f"three-group direct audit failed: {direct_error}")

    matrices = rng.uniform(0.01, 0.95, size=(41, 2, 2))
    with np.errstate(divide="ignore"):
        matrix_logs = np.log(matrices)
    binary = log_matrix_power_moments_binary(matrix_logs, 37)
    iterative = log_matrix_power_moments_batch(matrix_logs, 37)
    power_error = float(np.max(np.abs(binary - iterative)))
    if power_error > 2e-13:
        raise AssertionError(f"binary power audit failed: {power_error}")

    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    if receipt["schema"] != "rm2sub-fixed-rm-three-group-sparse-bridge-v1":
        raise AssertionError("unexpected receipt schema")
    parameters = receipt["parameters"]
    expected_parameters = {
        "message_bits": 1 << 16,
        "output_bits": 1 << 17,
        "outer_rows": 256,
        "outer_block_bits": 512,
        "step_bits": 64,
        "state_bits": 14,
        "bad_weight": 13107,
        "occupation_interval": [30, 52],
        "live_model": "support-averaged-preaddmul",
        "change_of_measure": "pointwise-density-envelope",
        "exact_zero_state_destination_split": True,
        "exact_nonactivation_verified_from_kernel_shell_counts": True,
    }
    for key, value in expected_parameters.items():
        if parameters[key] != value:
            raise AssertionError(f"parameter mismatch for {key}")
    groups = receipt["groups"]
    if [
        (
            row["name"],
            int(row["lower"]),
            int(row["upper"]),
            float(row["reference_probability"]),
        )
        for row in groups
    ] != list(GROUPS):
        raise AssertionError("group partition or reference schedule changed")

    receipt_compositions = receipt["composition_rows"]
    expected_compositions = [
        composition
        for occupation in range(30, 53)
        for composition in compositions_three(occupation)
    ]
    observed_compositions = [
        (
            int(row["low_rows"]),
            int(row["central_rows"]),
            int(row["high_rows"]),
        )
        for row in receipt_compositions
    ]
    if observed_compositions != expected_compositions:
        raise AssertionError("composition coverage or ordering changed")
    if len(observed_compositions) != 21275:
        raise AssertionError("unexpected composition count")

    by_composition = {
        composition: float(row["log2_upper_diagnostic"]) * LOG2
        for composition, row in zip(observed_compositions, receipt_compositions)
    }
    occupation_values = []
    for row in receipt["occupation_rows"]:
        occupation = int(row["occupation"])
        values = np.asarray(
            [
                by_composition[composition]
                for composition in expected_compositions
                if sum(composition) == occupation
            ]
        )
        aggregate = float(logsumexp(values))
        close(aggregate / LOG2, float(row["log2_upper_diagnostic"]))
        if -aggregate / LOG2 <= 40.0:
            raise AssertionError("an occupation fell below the comparison line")
        occupation_values.append(aggregate)
    interval = float(logsumexp(np.asarray(occupation_values)))
    close(interval / LOG2, float(receipt["interval_log2_upper_diagnostic"]))
    interval_margin = -interval / LOG2
    close(interval_margin, 1228.6379667587958)

    payload = {
        "schema": "rm2sub-three-group-bridge-audit-v1",
        "status": "PASS",
        "source_receipt": str(RECEIPT),
        "direct_collapse_maximum_absolute_error": direct_error,
        "binary_power_maximum_absolute_log_error": power_error,
        "checked_composition_count": len(observed_compositions),
        "checked_occupation_count": len(occupation_values),
        "interval_margin_bits_diagnostic": interval_margin,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"direct_collapse_error={direct_error:.17g}")
    print(f"binary_power_log_error={power_error:.17g}")
    print(f"interval_margin_bits_diagnostic={interval_margin:.12f}")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
