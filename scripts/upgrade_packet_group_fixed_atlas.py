#!/usr/bin/env python3
"""Upgrade full-support fixed witnesses with the total-weight outer branch."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from packet_group_outer_profile import fixed_total_weight_outer, normalization_log2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    upgraded = []
    group_bits = None
    for path in args.input:
        source = json.loads(path.read_text(encoding="utf-8"))
        rows = source if isinstance(source, list) else source["rows"]
        for row in rows:
            local_group = int(row["group_bits"])
            if group_bits is None:
                group_bits = local_group
            if local_group != group_bits:
                raise SystemExit("fixed-atlas upgrade: group mismatch")
            if "fugacities" not in row or "inner_constant_log2" not in row:
                raise SystemExit("fixed-atlas upgrade: full-support witness required")
            profile = np.asarray(row["profile"], dtype=np.float64)
            physical_weight = int(
                sum(index * int(count) for index, count in enumerate(row["profile"]))
            )
            outer_value, outer_details = fixed_total_weight_outer(physical_weight)
            outer_charge = (
                np.arange(local_group + 1, dtype=np.float64)
                * float(outer_details["charge_per_bit"])
            )
            inner_charge = np.log2(np.asarray(row["fugacities"], dtype=np.float64))
            charge = outer_charge + inner_charge
            constant = float(outer_details["constant_log2"]) + float(
                row["inner_constant_log2"]
            )
            normalization = float(normalization_log2(local_group, profile)[0])
            combined = constant - float(profile @ charge) - normalization
            candidate = dict(row)
            candidate.update(
                {
                    "name": f"{row['name']}_total_weight",
                    "source": str(path),
                    "outer_type": "total_weight",
                    "outer_charge": outer_charge.tolist(),
                    "outer_details": outer_details,
                    "outer_constant_log2": float(outer_details["constant_log2"]),
                    "outer_log2": outer_value,
                    "constant_log2": constant,
                    "charge": charge.tolist(),
                    "combined_log2": combined,
                    "margin_bits": float(row["target_log2"]) - combined,
                }
            )
            upgraded.append(candidate if combined < float(row["combined_log2"]) else row)

    report = {
        "status": "DIAGNOSTIC_BINARY64_TOTAL_WEIGHT_UPGRADED_FIXED_PROFILE_ATLAS",
        "group_bits": group_bits,
        "rows": upgraded,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"upgraded_witnesses={len(upgraded)} output={args.output}")


if __name__ == "__main__":
    main()
