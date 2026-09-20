#!/usr/bin/env python3
"""Consolidate the 10%-distance RM2Sub first-tranche screening receipts."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
Q1_PATH = HERE / "rm2sub_primary_tranche_d100.json"
Q1_SUPPLEMENT_PATH = HERE / "rm2sub_primary_q1_rm49_t128_s14_d100.json"
OUTPUT_JSON = HERE / "rm2sub_primary_frontier_summary_d100.json"
OUTPUT_CSV = HERE / "rm2sub_primary_frontier_summary_d100.csv"

Q2_RECEIPTS = {
    "t128_s13": (
        HERE / "rm2sub_primary_frontier_rm49_t128_s13_complete_grid_d100.json",
    ),
    "t128_s14": (HERE / "rm2sub_primary_frontier_rm49_t128_s14_d100.json",),
    "t128_s15": (HERE / "rm2sub_primary_frontier_rm49_t128_s15_d100.json",),
    "t128_s17": (HERE / "rm2sub_primary_frontier_rm49_t128_s17_d100.json",),
    "t128_s19": (HERE / "rm2sub_primary_frontier_rm49_t128_s19_d100.json",),
    "t256_s12": (HERE / "rm2sub_primary_frontier_rm49_t256_s12_d100.json",),
    "t256_s14": (HERE / "rm2sub_primary_frontier_rm49_t256_s14_d100.json",),
    "t256_s16": (HERE / "rm2sub_primary_frontier_rm49_t256_s16_d100.json",),
    "t256_s18": (HERE / "rm2sub_primary_frontier_rm49_t256_s18_d100.json",),
}


def structured_rm_case(payload: dict[str, object], occupation: int) -> dict[str, object]:
    matches = [
        row
        for row in payload["cases"]
        if row["family"] == "structured"
        and row["constituent"] == "RM(4,9) [512,256,32]"
        and int(row["occupation"]) == occupation
    ]
    if len(matches) != 1:
        raise AssertionError(f"expected one RM Q={occupation} case, found {len(matches)}")
    return matches[0]


def main() -> None:
    q1_payload = json.loads(Q1_PATH.read_text(encoding="utf-8"))
    if q1_payload["status"] != "BINARY64_DIAGNOSTIC":
        raise AssertionError("unexpected Q1 receipt status")
    if q1_payload["parameters"]["q2_enabled"]:
        raise AssertionError("principal Q1 receipt unexpectedly contains Q2")

    q1_supplement = json.loads(Q1_SUPPLEMENT_PATH.read_text(encoding="utf-8"))
    precise_q1 = {
        row["configuration"]: row
        for row in [*q1_payload["cases"], *q1_supplement["cases"]]
        if row["family"] == "structured"
        and row["constituent"] == "RM(4,9) [512,256,32]"
        and int(row["message_exponent"]) == 16
    }
    rows = []
    for configuration, paths in Q2_RECEIPTS.items():
        candidates = []
        grids = []
        for path in paths:
            payload = json.loads(path.read_text(encoding="utf-8"))
            row = structured_rm_case(payload, 2)
            candidates.append(row)
            parameters = payload["parameters"]
            grids.append(
                {
                    "receipt": str(path),
                    "tilt_count": parameters["tilt_count"],
                }
            )
        q2 = max(candidates, key=lambda row: float(row["margin_bits"]))
        q1 = precise_q1[configuration]
        rows.append(
            {
                "configuration": configuration,
                "step_bits": int(q1["step_bits"]),
                "state_bits": int(q1["state_bits"]),
                "persistence_exponent": float(q1["persistence_exponent"]),
                "a_minimum_distance": int(q1["a_minimum_distance"]),
                "kernel_minimum_distance": int(q1["kernel_minimum_distance"]),
                "kernel_weight_four": int(q1["kernel_weight_four"]),
                "q1_margin_bits": float(q1["margin_bits"]),
                "q1_tilt_count": q1_payload["parameters"]["tilt_count"],
                "q2_screen_margin_bits": float(q2["margin_bits"]),
                "q2_dominant_first_weight": int(q2["dominant_first_weight"]),
                "q2_dominant_second_weight": int(q2["dominant_second_weight"]),
                "q2_screen_grids": grids,
                "passes_q1_40_bit_screen": float(q1["margin_bits"]) >= 40.0,
                "passes_q2_40_bit_screen": float(q2["margin_bits"]) >= 40.0,
            }
        )

    payload = {
        "schema": "rm2sub-primary-frontier-summary-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "construction": (
            "one fixed RM(4,9) [512,256,32] constituent repeated at k=2^16, "
            "uniform routing, and one fixed audited RM2Sub A/B pair"
        ),
        "target": {
            "relative_distance": 0.10,
            "bad_weight": 13107,
            "output_bits": 131072,
            "screen_margin_bits": 40,
        },
        "q1_sources": [str(Q1_PATH), str(Q1_SUPPLEMENT_PATH)],
        "cases": rows,
        "limitations": [
            "Q1 uses a 121-point binary64 tilt grid; it is not outward rounded.",
            "Q2 is a coarse elimination screen over the stated finite witness sets.",
            "Occupations three through 256 are not covered.",
            "The rows are not distance certificates and do not include implementation cost.",
        ],
    }
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    fields = [
        "configuration",
        "step_bits",
        "state_bits",
        "persistence_exponent",
        "a_minimum_distance",
        "kernel_minimum_distance",
        "kernel_weight_four",
        "q1_margin_bits",
        "q1_tilt_count",
        "q2_screen_margin_bits",
        "q2_dominant_first_weight",
        "q2_dominant_second_weight",
        "passes_q1_40_bit_screen",
        "passes_q2_40_bit_screen",
    ]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"json={OUTPUT_JSON}")
    print(f"csv={OUTPUT_CSV}")


if __name__ == "__main__":
    main()
