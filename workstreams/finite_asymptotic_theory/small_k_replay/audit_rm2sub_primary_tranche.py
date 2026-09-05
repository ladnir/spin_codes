#!/usr/bin/env python3
"""Audit exact RM2Sub maps and receipt consistency for the primary tranche."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from generate_rm2sub_calibration_constituent import (
    enumerate_spectrum,
    generator_words,
    macwilliams_kernel,
    rank,
)


HERE = Path(__file__).resolve().parent
MAP_DIR = HERE / "rm2sub_calibration_constituents"
SUMMARY_JSON = HERE / "rm2sub_primary_frontier_summary_d100.json"
SUMMARY_CSV = HERE / "rm2sub_primary_frontier_summary_d100.csv"
OUTPUT = HERE / "rm2sub_primary_tranche_audit.json"
NESTED_DIR = HERE / "nested_rm2sub" / "t128_chain0"
CONFIGURATIONS = (
    "t128_s13",
    "t128_s14",
    "t128_s15",
    "t128_s17",
    "t128_s19",
    "t256_s12",
    "t256_s14",
    "t256_s16",
    "t256_s18",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def spectrum_vector(payload: dict[str, object], key: str, length: int) -> list[int]:
    result = [0] * (length + 1)
    for row in payload[key]:
        weight_key = "weight" if "weight" in row else "total_weight"
        count_key = "count" if "count" in row else "kernel_words"
        result[int(row[weight_key])] = int(row[count_key])
    return result


def audit_map(stem: str, map_dir: Path = MAP_DIR) -> dict[str, object]:
    selection_path = map_dir / f"{stem}_selection.json"
    a_path = map_dir / f"{stem}_a_spectrum.json"
    b_path = map_dir / f"{stem}_b_kernel_spectrum.json"
    selection_payload = json.loads(selection_path.read_text(encoding="utf-8"))
    selected = selection_payload["selected"]
    parameters = selection_payload["parameters"]
    step_bits = int(parameters["step_bits"])
    state_bits = int(parameters["state_bits"])
    if stem != f"t{step_bits}_s{state_bits}":
        raise AssertionError(f"name/parameter mismatch for {stem}")

    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    generators = [int(value, 16) for value in selected["A_generator_words_hex"]]
    if len(columns) != step_bits or len(generators) != state_bits:
        raise AssertionError(f"shape mismatch for {stem}")
    if generator_words(columns, state_bits) != generators:
        raise AssertionError(f"A/B transpose mismatch for {stem}")
    if rank(generators) != state_bits:
        raise AssertionError(f"rank failure for {stem}")
    if any((left & right).bit_count() & 1 for left in generators for right in generators):
        raise AssertionError(f"BA != 0 for {stem}")

    exact_a = enumerate_spectrum(generators, state_bits, step_bits)
    a_payload = json.loads(a_path.read_text(encoding="utf-8"))
    if exact_a != spectrum_vector(a_payload, "spectrum", step_bits):
        raise AssertionError(f"A spectrum mismatch for {stem}")
    exact_kernel = macwilliams_kernel(exact_a, state_bits)
    b_payload = json.loads(b_path.read_text(encoding="utf-8"))
    if exact_kernel != spectrum_vector(b_payload, "by_total_weight", step_bits):
        raise AssertionError(f"kernel spectrum mismatch for {stem}")

    a_distance = next(weight for weight, count in enumerate(exact_a[1:], 1) if count)
    kernel_distance = next(
        weight for weight, count in enumerate(exact_kernel[1:], 1) if count
    )
    if a_distance != int(selected["minimum_A_distance"]):
        raise AssertionError(f"A distance mismatch for {stem}")
    if kernel_distance != int(selected["minimum_kernel_distance"]):
        raise AssertionError(f"kernel distance mismatch for {stem}")
    if exact_kernel[4] != int(selected["weight_four_kernel_words"]):
        raise AssertionError(f"weight-four mismatch for {stem}")

    return {
        "configuration": stem,
        "step_bits": step_bits,
        "state_bits": state_bits,
        "a_minimum_distance": a_distance,
        "kernel_minimum_distance": kernel_distance,
        "kernel_weight_four": exact_kernel[4],
        "selection_sha256": sha256(selection_path),
        "a_spectrum_sha256": sha256(a_path),
        "kernel_spectrum_sha256": sha256(b_path),
    }


def main() -> None:
    maps = [audit_map(stem) for stem in CONFIGURATIONS]
    nested_maps = [audit_map(f"t128_s{state_bits}", NESTED_DIR) for state_bits in range(12, 17)]
    previous = None
    for state_bits in range(12, 17):
        selection = json.loads(
            (NESTED_DIR / f"t128_s{state_bits}_selection.json").read_text(
                encoding="utf-8"
            )
        )["selected"]
        generators = [int(value, 16) for value in selection["A_generator_words_hex"]]
        if previous is not None and generators[:-1] != previous:
            raise AssertionError(f"nested prefix mismatch at s={state_bits}")
        previous = generators
    summary = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))
    summary_rows = summary["cases"]
    if [row["configuration"] for row in summary_rows] != list(CONFIGURATIONS):
        raise AssertionError("summary configuration order mismatch")
    by_name = {row["configuration"]: row for row in maps}
    for row in summary_rows:
        exact = by_name[row["configuration"]]
        for field in (
            "step_bits",
            "state_bits",
            "a_minimum_distance",
            "kernel_minimum_distance",
            "kernel_weight_four",
        ):
            if row[field] != exact[field]:
                raise AssertionError(f"summary {field} mismatch for {row['configuration']}")

    with SUMMARY_CSV.open(newline="", encoding="utf-8") as source:
        csv_rows = list(csv.DictReader(source))
    if [row["configuration"] for row in csv_rows] != list(CONFIGURATIONS):
        raise AssertionError("CSV configuration order mismatch")
    for json_row, csv_row in zip(summary_rows, csv_rows, strict=True):
        if abs(float(json_row["q1_margin_bits"]) - float(csv_row["q1_margin_bits"])) > 1e-12:
            raise AssertionError("Q1 JSON/CSV mismatch")
        if abs(
            float(json_row["q2_screen_margin_bits"])
            - float(csv_row["q2_screen_margin_bits"])
        ) > 1e-12:
            raise AssertionError("Q2 JSON/CSV mismatch")

    payload = {
        "schema": "rm2sub-primary-tranche-audit-v1",
        "status": "PASS",
        "scope": (
            "exact A/B transpose, rank, BA=0, A spectrum, MacWilliams kernel "
            "spectrum, selected distances, and frontier JSON/CSV consistency"
        ),
        "maps": maps,
        "nested_maps": nested_maps,
        "nested_prefixes_verified": True,
        "frontier_summary": {
            "json": str(SUMMARY_JSON),
            "json_sha256": sha256(SUMMARY_JSON),
            "csv": str(SUMMARY_CSV),
            "csv_sha256": sha256(SUMMARY_CSV),
        },
        "limitations": [
            "This audit does not outward-round the numerical margins.",
            "This audit does not cover occupations three and above.",
            "This audit does not prove optimality of the selected A/B maps.",
            "The nested-family audit covers one unselected deterministic chain.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status=PASS")
    print(f"receipt={OUTPUT}")


if __name__ == "__main__":
    main()
