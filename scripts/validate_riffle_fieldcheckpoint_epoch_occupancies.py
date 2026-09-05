#!/usr/bin/env python3
"""Validate all target-size FieldCheckpointAccumulate epoch transfers."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_fieldcheckpoint_accumulate_oneblock import (
    occupancy_epoch_matrices,
    self_test,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal02_epoch_occupancy_validation.json"
)


def matrix_checks(matrices: list[np.ndarray]) -> dict[str, object]:
    stacked = np.stack(matrices)
    return {
        "matrix_count": len(matrices),
        "all_finite": bool(np.all(np.isfinite(stacked))),
        "minimum_entry": float(np.min(stacked)),
        "maximum_entry": float(np.max(stacked)),
        "maximum_row_sum": float(np.max(stacked.sum(axis=2))),
        "minimum_row_sum": float(np.min(stacked.sum(axis=2))),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=1024)
    parser.add_argument("--log-surprisal", type=float, default=-7.75)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    stochastic = occupancy_epoch_matrices(args.state_bits, args.epoch_bits, 1.0)
    z = math.exp(-math.exp(args.log_surprisal))
    tilted = occupancy_epoch_matrices(args.state_bits, args.epoch_bits, z)
    stochastic_rows = np.stack(stochastic).sum(axis=2)
    payload = {
        "schema": "riffle-fieldcheckpoint-epoch-occupancy-validation-v1",
        "parameters": {
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "visits_per_lane": args.epoch_bits // args.state_bits,
            "log_surprisal": args.log_surprisal,
            "z": z,
        },
        "small_exhaustive_tests": self_test(),
        "z_equals_one": {
            **matrix_checks(stochastic),
            "maximum_row_sum_error": float(
                np.max(np.abs(stochastic_rows - 1.0))
            ),
            "r2_zero_to_zero": float(stochastic[2][0, 0]),
            "r2_live_to_zero": float(stochastic[2][1, 0]),
        },
        "tilted": matrix_checks(tilted),
        "scope": (
            "Numerical validation of the exact coefficient formulas. "
            "This receipt is not an outward-rounded distance certificate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"maximum_row_sum_error,{payload['z_equals_one']['maximum_row_sum_error']:.17g}")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
