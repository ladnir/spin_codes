#!/usr/bin/env python3
"""Clone g=4 fixed witnesses with exact graph-averaged outer branches."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from probe_packet_group_exact_graph_linear_bl import (
    exact_graph_linear_outer,
    exact_graph_total_spectrum_outer,
)


def rows(path: Path) -> list[dict]:
    value = json.loads(path.read_text(encoding="utf-8"))
    result = value if isinstance(value, list) else value.get("rows", [])
    if not isinstance(result, list):
        raise ValueError(f"exact graph upgrader: no rows in {path}")
    return result


def upgrade(row: dict) -> dict | None:
    if int(row.get("group_bits", -1)) != 4 or "fugacities" not in row:
        return None
    outer_type = row.get("outer_type")
    details = dict(row.get("outer_details", {}))
    variables = details.get("log_variables")
    if variables is None:
        return None
    log_variables = np.asarray(variables, dtype=np.float64)
    profile = [int(value) for value in row["profile"]]
    if outer_type == "linear_bl":
        band1 = float(details.get("band1_coefficient", 0.5))
        outer_value, exact_details = exact_graph_linear_outer(
            profile, log_variables, band1
        )
        new_type = "exact_graph_linear_bl"
    elif outer_type == "total_spectrum":
        outer_point = details.get("outer_point")
        if not isinstance(outer_point, list) or len(outer_point) < 5:
            return None
        log_beta = float(outer_point[4])
        outer_value, exact_details = exact_graph_total_spectrum_outer(
            profile, log_variables, log_beta
        )
        new_type = "exact_graph_total_spectrum"
    else:
        return None

    outer_charge = log_variables / math.log(2.0)
    outer_constant = outer_value + float(np.asarray(profile) @ outer_charge)
    fugacities = np.asarray(row["fugacities"], dtype=np.float64)
    inner_charge = np.log2(fugacities)
    charge = outer_charge + inner_charge
    constant = outer_constant + float(row["inner_constant_log2"])
    combined = (
        constant
        - float(np.asarray(profile) @ charge)
        - float(row["normalization_log2"])
    )
    target = float(row["target_log2"])
    return {
        **row,
        "name": str(row.get("name", "witness")) + "__exact_graph",
        "outer_type": new_type,
        "outer_details": {**details, **exact_details},
        "outer_log2": outer_value,
        "outer_constant_log2": outer_constant,
        "outer_charge": outer_charge.tolist(),
        "constant_log2": constant,
        "charge": charge.tolist(),
        "combined_log2": combined,
        "margin_bits": target - combined,
        "status": "DIAGNOSTIC_BINARY64_EXACT_GRAPH_UPGRADED_FIXED_PROFILE_WITNESS",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    upgraded = []
    seen = set()
    for artifact in args.artifact:
        for row in rows(artifact):
            result = upgrade(row)
            if result is None:
                continue
            key = (tuple(result["profile"]), result["outer_type"], tuple(result["charge"]))
            if key not in seen:
                seen.add(key)
                upgraded.append(result)
    report = {
        "status": "DIAGNOSTIC_BINARY64_G4_EXACT_GRAPH_UPGRADED_ATLAS",
        "group_bits": 4,
        "complete": True,
        "source_artifacts": [str(path) for path in args.artifact],
        "rows": upgraded,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "upgraded_witnesses": len(upgraded),
                "minimum_margin_bits": min(
                    (row["margin_bits"] for row in upgraded), default=None
                ),
                "closing_witnesses": int(
                    sum(bool(row["margin_bits"] >= 0) for row in upgraded)
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
