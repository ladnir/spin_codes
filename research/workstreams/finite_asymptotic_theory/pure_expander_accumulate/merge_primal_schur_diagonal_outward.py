#!/usr/bin/env python3
"""Merge outward diagonal bands and certify requested accumulator shells."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from certify_primal_schur_diagonal_outward import shell_ratio_upper


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "primal_schur_shells_outward.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, action="append", required=True)
    parser.add_argument("--shell-weight", type=int, action="append", required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    parameters: tuple[int, int, int] | None = None
    values: dict[tuple[int, int], np.float64] = {}
    sources = []
    for path in args.receipt:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") != "OUTWARD_BINARY64_UPPER_BOUND":
            raise ValueError(f"{path} is not an outward diagonal receipt")
        row_parameters = payload["parameters"]
        current = (
            int(row_parameters["message_bits"]),
            int(row_parameters["output_bits"]),
            int(row_parameters["right_degree"]),
        )
        if parameters is None:
            parameters = current
        elif current != parameters:
            raise ValueError("receipt parameters do not agree")
        for row in payload["entries"]:
            key = (int(row["sector"]), int(row["level"]))
            value = np.float64(row["scaled_diagonal_upper"])
            if key in values and value != values[key]:
                raise ValueError(f"conflicting duplicate diagonal {key}")
            values[key] = value
        sources.append(str(path.resolve()))

    assert parameters is not None
    message_bits, output_bits, right_degree = parameters
    shell_rows = []
    for shell_weight in sorted(set(args.shell_weight)):
        maximum_level = max(2 * shell_weight, 2 * (output_bits - shell_weight) + 1)
        maximum_level = min(maximum_level, output_bits // 2)
        # The low-side shells used by this workstream have support through 2w.
        if shell_weight <= output_bits // 2:
            maximum_level = 2 * shell_weight
        required = {
            (sector, level)
            for level in range(1, maximum_level + 1)
            for sector in range(min(2, level) + 1)
        }
        missing = sorted(required - values.keys())
        if missing:
            raise ValueError(
                f"shell {shell_weight} is missing {len(missing)} entries; first={missing[0]}"
            )
        ratio, mean_lower = shell_ratio_upper(
            message_bits,
            output_bits,
            right_degree,
            shell_weight,
            values,
        )
        shell_rows.append(
            {
                "shell_weight": shell_weight,
                "active_level_maximum": maximum_level,
                "mean_lower": float(mean_lower),
                "mean_log2_lower": math.log2(float(mean_lower)),
                "variance_to_mean_upper": float(ratio),
                "variance_to_mean_log2_upper": math.log2(float(ratio)),
                "passes_factor_512": bool(ratio <= 512),
            }
        )

    payload = {
        "schema": "pure-ea-primal-schur-shell-merge-outward-v1",
        "status": "OUTWARD_BINARY64_SHELL_CERTIFICATE",
        "parameters": {
            "message_bits": message_bits,
            "output_bits": output_bits,
            "right_degree": right_degree,
        },
        "source_receipts": sources,
        "shells": shell_rows,
        "claim": {
            "all_requested_shells_pass_factor_512": all(
                row["passes_factor_512"] for row in shell_rows
            ),
            "maximum_variance_to_mean_upper": max(
                row["variance_to_mean_upper"] for row in shell_rows
            ),
        },
        "scope": [
            "This receipt composes outward diagonal upper bounds with exact accumulator run counts.",
            "Complementary high shells require a separately stated symmetry argument before this receipt is used for them.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
