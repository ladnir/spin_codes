#!/usr/bin/env python3
"""Audit the fixed-RM RM2Sub occupation-ladder receipts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "rm2sub_q_ladder_audit.json"
CASES = (
    ("rm2sub_q_ladder_rm49_t64_s14_d100.json", 64, 14, list(range(3, 9))),
    ("rm2sub_q_ladder_rm49_t128_s13_d100.json", 128, 13, list(range(3, 9))),
    ("rm2sub_q_ladder_rm49_t256_s12_d100.json", 256, 12, list(range(3, 9))),
    ("rm2sub_q_ladder_rm49_t64_s14_q09_q16_d100.json", 64, 14, list(range(9, 17))),
    ("rm2sub_q_ladder_rm49_t64_s14_q17_q32_d100.json", 64, 14, list(range(17, 33))),
    ("rm2sub_q_ladder_rm49_t64_s14_q30_q32_fine_d100.json", 64, 14, list(range(30, 33))),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    recurrence_path = HERE / "rm2sub_two_colour_recurrence_audit.json"
    recurrence = json.loads(recurrence_path.read_text(encoding="utf-8"))
    if recurrence["status"] != "PASSED_BINARY64_IDENTITY_CHECK":
        raise AssertionError("two-colour recurrence audit did not pass")
    if int(recurrence["checked_pairs"]) != 15:
        raise AssertionError("unexpected recurrence audit coverage")
    if float(recurrence["maximum_absolute_error"]) > 2e-14:
        raise AssertionError("recurrence audit error exceeds tolerance")

    audited = []
    for filename, step_bits, state_bits, occupations in CASES:
        path = HERE / filename
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload["schema"] != "rm2sub-fixed-occupation-rm-holder-v1":
            raise AssertionError(f"schema mismatch in {filename}")
        parameters = payload["parameters"]
        expected = {
            "message_bits": 1 << 16,
            "output_bits": 1 << 17,
            "outer_rows": 256,
            "outer_block_bits": 512,
            "outer_dimension": 256,
            "bad_weight": 13107,
            "step_bits": step_bits,
            "state_bits": state_bits,
            "persistence_exponent": 20,
            "occupations": occupations,
        }
        for key, value in expected.items():
            if parameters[key] != value:
                raise AssertionError(f"{filename}: {key} mismatch")
        rows = payload["occupation_rows"]
        if [int(row["occupation"]) for row in rows] != occupations:
            raise AssertionError(f"{filename}: occupation order mismatch")
        for row in rows:
            occupation = int(row["occupation"])
            if len(row["mixture_rows"]) != occupation + 1:
                raise AssertionError(f"{filename}: incomplete all-one mixture")
            if int(row["dominant_all_one_rows"]) != 0:
                raise AssertionError(f"{filename}: unexpected dominant face")
        audited.append(
            {
                "path": str(path),
                "sha256": sha256(path),
                "step_bits": step_bits,
                "state_bits": state_bits,
                "occupations": occupations,
                "minimum_margin_bits": min(
                    float(row["margin_bits_diagnostic"]) for row in rows
                ),
            }
        )

    matched = [
        json.loads((HERE / filename).read_text(encoding="utf-8"))
        for filename, _step, _state, _occupations in CASES[:3]
    ]
    for index in range(6):
        margins = [
            float(payload["occupation_rows"][index]["margin_bits_diagnostic"])
            for payload in matched
        ]
        if not margins[0] > margins[1] > margins[2]:
            raise AssertionError("matched-persistence ordering mismatch")

    fine = json.loads((HERE / CASES[-1][0]).read_text(encoding="utf-8"))
    fine_margins = {
        int(row["occupation"]): float(row["margin_bits_diagnostic"])
        for row in fine["occupation_rows"]
    }
    if fine_margins[30] >= 0.0:
        raise AssertionError("fine Q=30 boundary unexpectedly closed")

    payload = {
        "schema": "rm2sub-q-ladder-audit-v1",
        "status": "PASS",
        "scope": (
            "receipt schemas and parameters, complete all-one mixtures, "
            "matched-persistence ordering, boundary sign, and independent "
            "two-colour recurrence identity check"
        ),
        "recurrence_audit": {
            "path": str(recurrence_path),
            "sha256": sha256(recurrence_path),
            "maximum_absolute_error": recurrence["maximum_absolute_error"],
        },
        "receipts": audited,
        "limitations": [
            "This audit does not outward-round the numerical margins.",
            "It checks the recorded finite witness grids, not witness optimality.",
            "It does not close occupations 30 through 256.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"receipt={OUTPUT}")


if __name__ == "__main__":
    main()
