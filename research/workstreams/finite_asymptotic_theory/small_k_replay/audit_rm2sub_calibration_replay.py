#!/usr/bin/env python3
"""Audit RM2Sub calibration constituents and replay receipts."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
GENERATED = HERE / "rm2sub_calibration_constituents"
OUTPUT = HERE / "rm2sub_calibration_replay_audit.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def krawtchouk(length: int, degree: int, point: int) -> int:
    return sum(
        (-1 if intersection & 1 else 1)
        * math.comb(point, intersection)
        * math.comb(length - point, degree - intersection)
        for intersection in range(
            max(0, degree - (length - point)), min(degree, point) + 1
        )
    )


def audit_constituent(stem: str) -> dict[str, object]:
    selection_path = GENERATED / f"{stem}_selection.json"
    a_path = GENERATED / f"{stem}_a_spectrum.json"
    b_path = GENERATED / f"{stem}_b_kernel_spectrum.json"
    selection_payload = json.loads(selection_path.read_text(encoding="utf-8"))
    a_payload = json.loads(a_path.read_text(encoding="utf-8"))
    b_payload = json.loads(b_path.read_text(encoding="utf-8"))
    parameters = selection_payload["parameters"]
    selected = selection_payload["selected"]
    length = int(parameters["step_bits"])
    state_bits = int(parameters["state_bits"])
    generators = [int(word, 16) for word in selected["A_generator_words_hex"]]
    columns = [int(column, 16) for column in selected["B_columns_hex"]]
    if len(generators) != state_bits or len(columns) != length:
        raise AssertionError(f"{stem}: generator/column shape mismatch")
    for generator in generators:
        image = 0
        word = generator
        while word:
            low = word & -word
            image ^= columns[low.bit_length() - 1]
            word ^= low
        if image:
            raise AssertionError(f"{stem}: BA is nonzero")
    if len(set(columns)) != length or 0 in columns:
        raise AssertionError(f"{stem}: B columns are not distinct nonzero")

    a = [0] * (length + 1)
    for row in a_payload["spectrum"]:
        a[int(row["weight"])] = int(row["count"])
    if sum(a) != 1 << state_bits or a[0] != 1:
        raise AssertionError(f"{stem}: A spectrum mass mismatch")
    a_distance = next(weight for weight, count in enumerate(a[1:], 1) if count)
    if a_distance != int(selected["minimum_A_distance"]):
        raise AssertionError(f"{stem}: A distance mismatch")

    b = [0] * (length + 1)
    for row in b_payload["by_total_weight"]:
        b[int(row["total_weight"])] = int(row["kernel_words"])
    if sum(b) != 1 << (length - state_bits) or b[0] != 1:
        raise AssertionError(f"{stem}: kernel spectrum mass mismatch")
    b_distance = next(weight for weight, count in enumerate(b[1:], 1) if count)
    if b_distance != int(selected["minimum_kernel_distance"]):
        raise AssertionError(f"{stem}: kernel distance mismatch")
    if b[4] != int(selected["weight_four_kernel_words"]):
        raise AssertionError(f"{stem}: weight-four kernel count mismatch")

    support = [(weight, count) for weight, count in enumerate(a) if count]
    denominator = 1 << state_bits
    for weight, claimed in enumerate(b):
        numerator = sum(
            count * krawtchouk(length, weight, dual_weight)
            for dual_weight, count in support
        )
        if numerator % denominator or numerator // denominator != claimed:
            raise AssertionError(f"{stem}: MacWilliams mismatch at weight {weight}")
    return {
        "stem": stem,
        "step_bits": length,
        "state_bits": state_bits,
        "a_minimum_distance": a_distance,
        "kernel_minimum_distance": b_distance,
        "kernel_weight_four": b[4],
        "BA_zero": True,
        "macwilliams_exact": True,
        "files": {
            str(path.relative_to(HERE)): sha256(path)
            for path in (selection_path, a_path, b_path)
        },
    }


def audit_json_csv(json_path: Path, csv_path: Path) -> dict[str, object]:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    with csv_path.open(newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    if len(rows) != len(payload["cases"]):
        raise AssertionError(f"{json_path.name}: JSON/CSV row-count mismatch")
    for left, right in zip(payload["cases"], rows):
        for key in ("inner", "constituent", "message_exponent", "margin_bits"):
            if key in left and key in right and str(left[key]) != right[key]:
                raise AssertionError(
                    f"{json_path.name}: JSON/CSV mismatch for {key}"
                )
    return {
        "json": json_path.name,
        "json_sha256": sha256(json_path),
        "csv": csv_path.name,
        "csv_sha256": sha256(csv_path),
        "case_count": len(rows),
    }


def find_case(cases: list[dict[str, object]], **wanted: object) -> dict[str, object]:
    matches = [
        row for row in cases
        if all(str(row.get(key)) == str(value) for key, value in wanted.items())
    ]
    if len(matches) != 1:
        raise AssertionError(f"case lookup {wanted} returned {len(matches)} rows")
    return matches[0]


def main() -> None:
    stems = [f"t64_s{state_bits}" for state_bits in range(7, 17)]
    stems += ["t64_s18", "t64_s19", "t64_s20", "t256_s18"]
    constituents = [audit_constituent(stem) for stem in stems]

    pairs = [
        (
            HERE / "rm2sub_epoch_calibration_bch_q1_d100.json",
            HERE / "rm2sub_epoch_calibration_bch_q1_d100.csv",
        ),
        (
            HERE / "rm2sub_t64_state_calibration_bch_q1_d100.json",
            HERE / "rm2sub_t64_state_calibration_bch_q1_d100.csv",
        ),
        (
            HERE / "rm2sub_q2_calibration_bch_k13_d100.json",
            HERE / "rm2sub_q2_calibration_bch_k13_d100.csv",
        ),
        (
            HERE / "rate_half_family_k_margin_d100_rm2sub_t64_matched.json",
            HERE / "rate_half_family_k_margin_d100_rm2sub_t64_matched.csv",
        ),
    ]
    receipts = [audit_json_csv(*pair) for pair in pairs]

    family = json.loads(pairs[-1][0].read_text(encoding="utf-8"))
    cases = family["cases"]
    for row in cases:
        exponent = int(row["message_exponent"])
        if int(row["state_bits"]) != max(7, exponent - 4):
            raise AssertionError("family replay state schedule mismatch")
        if int(row["step_bits"]) != 64:
            raise AssertionError("family replay epoch mismatch")
        if int(row["outer_rows"]) % 64:
            raise AssertionError("family replay includes an inadmissible region")
        if int(row["output_bits"]) != 2 * int(row["message_bits"]):
            raise AssertionError("family replay rate-half identity mismatch")
        if int(row["bad_weight"]) != math.floor(0.1 * int(row["output_bits"])):
            raise AssertionError("family replay bad-weight convention mismatch")

    anchors = [
        (
            {
                "series": "RM(4,9) [512,256,32] exact",
                "message_exponent": 16,
            },
            41.46156261012271,
        ),
        (
            {
                "series": "extended BCH [128,64,22] exact",
                "message_exponent": 16,
            },
            33.96633639389293,
        ),
        (
            {
                "series": "random full-rank [256,128] reused",
                "message_exponent": 16,
            },
            62.31811059396261,
        ),
    ]
    checked_anchors = []
    for key, expected in anchors:
        row = find_case(cases, **key)
        actual = float(row["margin_bits"])
        if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-11):
            raise AssertionError(f"anchor mismatch for {key}: {actual}")
        checked_anchors.append({**key, "margin_bits": actual})

    payload = {
        "schema": "rm2sub-calibration-replay-audit-v1",
        "status": "PASS",
        "exact_checks": {
            "selected_constituents": constituents,
            "family_schedule_and_admissibility": True,
            "json_csv_consistency": receipts,
        },
        "binary64_anchor_checks": checked_anchors,
        "scope": (
            "Exact integer checks authenticate the selected A/B pairs and "
            "receipt structure. Matching stored binary64 values does not make "
            "the transfer calculations outward rounded."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,PASS")
    print(f"output,{OUTPUT}")


if __name__ == "__main__":
    main()
