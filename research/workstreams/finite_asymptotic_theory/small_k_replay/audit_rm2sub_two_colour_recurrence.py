#!/usr/bin/env python3
"""Check the two-colour region recurrence against direct enumeration."""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np

from evaluate_rm2sub_q_ladder import (
    mixed_epoch_log_matrices,
    mixed_region_log_matrices,
)


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "rm2sub_two_colour_recurrence_audit.json"


def direct_average(
    impulses: np.ndarray,
    probability: float,
    positions: int,
    epoch_bits: int,
    regular: int,
    forced: int,
) -> np.ndarray:
    """Average directly over disjoint placements and regular bit values."""
    result = np.zeros((2, 2), dtype=np.float64)
    placement_count = 0
    universe = range(positions)
    for regular_positions in itertools.combinations(universe, regular):
        remaining = [item for item in universe if item not in regular_positions]
        for forced_positions in itertools.combinations(remaining, forced):
            placement_count += 1
            regular_set = set(regular_positions)
            forced_set = set(forced_positions)
            for bits in itertools.product((0, 1), repeat=regular):
                live_regular = {
                    position
                    for position, bit in zip(regular_positions, bits)
                    if bit
                }
                bit_mass = (
                    probability ** len(live_regular)
                    * (1.0 - probability) ** (regular - len(live_regular))
                )
                product = np.eye(2, dtype=np.float64)
                for start in range(0, positions, epoch_bits):
                    epoch = set(range(start, start + epoch_bits))
                    active = len(epoch & live_regular) + len(epoch & forced_set)
                    product = product @ impulses[active]
                result += bit_mass * product
    return result / placement_count


def main() -> None:
    epoch_bits = 3
    epochs_per_region = 2
    positions = epoch_bits * epochs_per_region
    maximum_occupation = 4
    probability = 0.37
    impulses = np.asarray(
        [
            [[0.91, 0.07], [0.03, 0.83]],
            [[0.61, 0.29], [0.11, 0.67]],
            [[0.43, 0.37], [0.19, 0.53]],
            [[0.31, 0.41], [0.23, 0.47]],
            [[0.27, 0.43], [0.29, 0.41]],
        ],
        dtype=np.float64,
    )
    epoch_logs = mixed_epoch_log_matrices(
        impulses, probability, maximum_occupation
    )
    recovered = mixed_region_log_matrices(
        epoch_logs, epoch_bits, epochs_per_region, maximum_occupation
    )
    rows = []
    maximum_error = 0.0
    for regular in range(maximum_occupation + 1):
        for forced in range(maximum_occupation - regular + 1):
            reference = direct_average(
                impulses,
                probability,
                positions,
                epoch_bits,
                regular,
                forced,
            )
            candidate = np.exp(recovered[regular, forced])
            error = float(np.max(np.abs(reference - candidate)))
            maximum_error = max(maximum_error, error)
            rows.append(
                {
                    "regular": regular,
                    "forced": forced,
                    "maximum_absolute_error": error,
                }
            )
    if maximum_error > 2e-14:
        raise ArithmeticError(f"recurrence audit failed: {maximum_error}")
    payload = {
        "schema": "rm2sub-two-colour-recurrence-audit-v1",
        "status": "PASSED_BINARY64_IDENTITY_CHECK",
        "parameters": {
            "positions": positions,
            "epoch_bits": epoch_bits,
            "epochs_per_region": epochs_per_region,
            "maximum_occupation": maximum_occupation,
            "regular_bit_probability": probability,
        },
        "checked_pairs": len(rows),
        "maximum_absolute_error": maximum_error,
        "rows": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"checked_pairs={len(rows)}")
    print(f"maximum_absolute_error={maximum_error:.3e}")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
