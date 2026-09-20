#!/usr/bin/env python3
"""Combine the strongest fixed-RM RM2Sub diagnostic bound at every Q."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
Q1_PATH = HERE / "rm2sub_epoch_geometry_rm49_t64_s14_q1_d100.json"
Q2_PATH = HERE / "rm2sub_epoch_geometry_rm49_t64_s14_q2_d100.json"
Q34_PATH = HERE / "rm2sub_q_ladder_rm49_t64_s14_d100.json"
Q5_29_PATH = HERE / "rm2sub_three_group_bridge_q03_q52_probe_d100.json"
Q30_256_PATH = HERE / "rm2sub_refined_band_bridge_q30_q256_diagnostic.json"
Q30_256_AUDIT_PATH = HERE / "rm2sub_refined_band_bridge_audit.json"
OUTPUT = HERE / "rm2sub_full_occupation_q1_q256_diagnostic.json"
AUDIT_OUTPUT = HERE / "rm2sub_full_occupation_q1_q256_audit.json"
LOG2 = math.log(2.0)


def close(first: float, second: float, tolerance: float = 3e-9) -> None:
    if not math.isclose(first, second, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"value mismatch: {first} != {second}")


def main() -> None:
    q1 = json.loads(Q1_PATH.read_text(encoding="utf-8"))
    q2 = json.loads(Q2_PATH.read_text(encoding="utf-8"))
    q34 = json.loads(Q34_PATH.read_text(encoding="utf-8"))
    q5_29 = json.loads(Q5_29_PATH.read_text(encoding="utf-8"))
    q30_256 = json.loads(Q30_256_PATH.read_text(encoding="utf-8"))
    q30_audit = json.loads(Q30_256_AUDIT_PATH.read_text(encoding="utf-8"))
    if q30_audit["status"] != "PASS":
        raise AssertionError("the refined-band bridge audit did not pass")

    expected_case_parameters = {
        "step_bits": 64,
        "state_bits": 14,
        "message_bits": 1 << 16,
        "outer_rows": 256,
        "output_bits": 1 << 17,
        "bad_weight": 13107,
        "outer_model": "fixed-exact-spectrum",
    }
    expected_receipt_parameters = {
        "step_bits": 64,
        "state_bits": 14,
        "message_bits": 1 << 16,
        "outer_rows": 256,
        "output_bits": 1 << 17,
        "bad_weight": 13107,
    }

    selected_rows = []
    q1_case = next(
        row
        for row in q1["cases"]
        if row["family"] == "structured" and int(row["occupation"]) == 1
    )
    for key, value in expected_case_parameters.items():
        if q1_case[key] != value:
            raise AssertionError(f"Q1 parameter mismatch for {key}")
    selected_rows.append(
        {
            "occupation": 1,
            "margin_bits_diagnostic": float(q1_case["margin_bits"]),
            "log2_upper_diagnostic": -float(q1_case["margin_bits"]),
            "source": Q1_PATH.name,
        }
    )
    q2_case = next(
        row
        for row in q2["cases"]
        if row["family"] == "structured" and int(row["occupation"]) == 2
    )
    for key, value in expected_case_parameters.items():
        if q2_case[key] != value:
            raise AssertionError(f"Q2 parameter mismatch for {key}")
    for receipt_name, receipt in (("Q3--4", q34), ("Q5--29", q5_29)):
        for key, value in expected_receipt_parameters.items():
            if receipt["parameters"][key] != value:
                raise AssertionError(
                    f"{receipt_name} parameter mismatch for {key}"
                )
    if q30_256["bad_event"] != (
        "some message of outer occupation 30 through 256 encodes to weight at most 13107"
    ):
        raise AssertionError("Q30--256 bad event changed")
    selected_rows.append(
        {
            "occupation": 2,
            "margin_bits_diagnostic": float(q2_case["margin_bits"]),
            "log2_upper_diagnostic": -float(q2_case["margin_bits"]),
            "source": Q2_PATH.name,
        }
    )

    q34_by_q = {
        int(row["occupation"]): row for row in q34["occupation_rows"]
    }
    for occupation in (3, 4):
        row = q34_by_q[occupation]
        selected_rows.append(
            {
                "occupation": occupation,
                "margin_bits_diagnostic": float(row["margin_bits_diagnostic"]),
                "log2_upper_diagnostic": float(row["log2_upper_diagnostic"]),
                "source": Q34_PATH.name,
            }
        )

    compositions_by_q: dict[int, list[float]] = {}
    for row in q5_29["composition_rows"]:
        occupation = int(row["occupation"])
        if 5 <= occupation <= 29:
            compositions_by_q.setdefault(occupation, []).append(
                float(row["log2_upper_diagnostic"]) * LOG2
            )
    bridge_rows = {
        int(row["occupation"]): row for row in q5_29["occupation_rows"]
    }
    checked_middle_compositions = 0
    for occupation in range(5, 30):
        values = compositions_by_q[occupation]
        expected_count = (occupation + 1) * (occupation + 2) // 2
        if len(values) != expected_count:
            raise AssertionError(f"incomplete Q={occupation} composition set")
        aggregate = float(logsumexp(np.asarray(values)))
        source_row = bridge_rows[occupation]
        close(aggregate / LOG2, float(source_row["log2_upper_diagnostic"]))
        selected_rows.append(
            {
                "occupation": occupation,
                "margin_bits_diagnostic": -aggregate / LOG2,
                "log2_upper_diagnostic": aggregate / LOG2,
                "source": Q5_29_PATH.name,
            }
        )
        checked_middle_compositions += expected_count

    for source_row in q30_256["occupation_rows"]:
        selected_rows.append(
            {
                "occupation": int(source_row["occupation"]),
                "margin_bits_diagnostic": float(
                    source_row["margin_bits_diagnostic"]
                ),
                "log2_upper_diagnostic": float(
                    source_row["log2_upper_diagnostic"]
                ),
                "source": Q30_256_PATH.name,
            }
        )
    if [row["occupation"] for row in selected_rows] != list(range(1, 257)):
        raise AssertionError("occupation coverage is not exactly Q=1 through Q=256")

    logs = np.asarray(
        [float(row["log2_upper_diagnostic"]) * LOG2 for row in selected_rows]
    )
    aggregate = float(logsumexp(logs))
    margin = -aggregate / LOG2
    weakest = min(
        selected_rows, key=lambda row: float(row["margin_bits_diagnostic"])
    )
    if margin <= 40.0:
        raise AssertionError("the full diagnostic union does not clear 40 bits")

    payload = {
        "schema": "rm2sub-fixed-rm-full-occupation-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "construction": (
            "one fixed RM(4,9) constituent repeated in 256 rows, uniform "
            "routing, and the fixed audited RM2Sub t=64,s=14 A/B pair"
        ),
        "parameters": {
            "message_bits": 1 << 16,
            "output_bits": 1 << 17,
            "bad_weight": 13107,
            "distance_fraction": 0.10,
            "occupation_interval": [1, 256],
        },
        "source_schedule": [
            {"occupation_interval": [1, 1], "receipt": Q1_PATH.name},
            {"occupation_interval": [2, 2], "receipt": Q2_PATH.name},
            {"occupation_interval": [3, 4], "receipt": Q34_PATH.name},
            {"occupation_interval": [5, 29], "receipt": Q5_29_PATH.name},
            {"occupation_interval": [30, 256], "receipt": Q30_256_PATH.name},
        ],
        "occupation_rows": selected_rows,
        "union_log2_upper_diagnostic": aggregate / LOG2,
        "union_margin_bits_diagnostic": margin,
        "union_gap_to_40_bits": margin - 40.0,
        "weakest_occupation": {
            "occupation": int(weakest["occupation"]),
            "margin_bits_diagnostic": float(weakest["margin_bits_diagnostic"]),
        },
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "The source receipts use finite Chernoff and Holder witness grids.",
            "This receipt is not a formal distance certificate.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    audit = {
        "schema": "rm2sub-fixed-rm-full-occupation-audit-v1",
        "status": "PASS",
        "consolidated_receipt": str(OUTPUT),
        "covered_occupation_interval": [1, 256],
        "checked_occupation_count": len(selected_rows),
        "checked_q5_q29_composition_count": checked_middle_compositions,
        "checked_q30_q256_composition_count": int(
            q30_audit["checked_composition_count"]
        ),
        "union_margin_bits_diagnostic": margin,
        "weakest_occupation": int(weakest["occupation"]),
    }
    AUDIT_OUTPUT.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"union_margin_bits_diagnostic={margin:.12f}")
    print(f"gap_to_40_bits={margin - 40.0:.12f}")
    print(f"weakest_occupation={weakest['occupation']}")
    print(f"output={OUTPUT}")
    print(f"audit={AUDIT_OUTPUT}")


if __name__ == "__main__":
    main()
