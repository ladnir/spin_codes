#!/usr/bin/env python3
"""Convert a binary-edge scan into reusable fixed combined witnesses."""

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
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--base-class", type=int, required=True)
    parser.add_argument("--active-class", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for path in args.input:
        source = json.loads(path.read_text(encoding="utf-8"))
        for index, row in enumerate(source):
            if "linear_bl_outer_details" not in row:
                continue
            profile = np.asarray(row["profile"], dtype=np.float64)
            normalization = float(normalization_log2(args.group, profile)[0])
            inner_charge = np.zeros(args.group + 1)
            inner_charge[args.active_class] = math.log2(float(row["ratio"]))
            inner_constant = (
                float(row["inner_probability_log2"])
                + float(profile @ inner_charge)
                + normalization
            )
            outer_details = row["linear_bl_outer_details"]
            linear_charge = (
                np.asarray(outer_details["log_variables"], dtype=np.float64)
                / math.log(2.0)
            )
            linear_constant = (
                float(row["linear_bl_outer_log2"])
                + float(profile @ linear_charge)
            )
            physical_weight = int(
                sum(index * int(count) for index, count in enumerate(row["profile"]))
            )
            total_value, total_details = fixed_total_weight_outer(physical_weight)
            if total_value < float(row["linear_bl_outer_log2"]):
                outer_type = "total_weight"
                outer_charge = (
                    np.arange(args.group + 1, dtype=np.float64)
                    * float(total_details["charge_per_bit"])
                )
                outer_constant = float(total_details["constant_log2"])
                outer_anchor = total_value
            else:
                outer_type = "linear_bl"
                outer_charge = linear_charge
                outer_constant = linear_constant
                outer_anchor = float(row["linear_bl_outer_log2"])
            rows.append(
                {
                    "name": f"edge_{args.base_class}_{args.active_class}_{int(row['active_atoms'])}",
                    "source": str(path),
                    "source_index": index,
                    "profile": row["profile"],
                    "support": sorted((args.base_class, args.active_class)),
                    "constant_log2": outer_constant + inner_constant,
                    "charge": (outer_charge + inner_charge).tolist(),
                    "anchor_combined_log2": outer_anchor
                    + float(row["inner_probability_log2"]),
                    "outer_type": outer_type,
                    "pole": row["pole"],
                    "ratio": row["ratio"],
                    "outer_constant_log2": outer_constant,
                    "inner_constant_log2": inner_constant,
                }
            )
    report = {
        "status": "DIAGNOSTIC_BINARY64_GENERALIZED_FIXED_COMBINED_WITNESSES",
        "group_bits": args.group,
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"fixed_combined_witnesses={len(rows)} output={args.output}")


if __name__ == "__main__":
    main()
