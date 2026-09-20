#!/usr/bin/env python3
"""Audit and compact the complete refined-band RM2Sub bridge receipts."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
SOURCE_NAMES = (
    "rm2sub_three_group_bridge_q30_q52_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q053_q096_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q097_q128_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q129_q160_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q161_q176_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q177_q192_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q193_q208_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q209_q224_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q225_q232_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q233_q240_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q241_q248_refined_bands_probe_d100.json",
    "rm2sub_three_group_bridge_q249_q256_refined_bands_probe_d100.json",
)
OUTPUT = HERE / "rm2sub_refined_band_bridge_q30_q256_diagnostic.json"
AUDIT_OUTPUT = HERE / "rm2sub_refined_band_bridge_audit.json"
LOG2 = math.log(2.0)
EXPECTED_GROUPS = (
    ("low", 1, 95, 0.25),
    ("central", 96, 416, 0.5),
    ("high", 417, 512, 0.8677722630069483),
)


def close(first: float, second: float, tolerance: float = 3e-9) -> None:
    if not math.isclose(first, second, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"aggregation mismatch: {first} != {second}")


def main() -> None:
    occupation_rows = []
    source_rows = []
    seen_occupations = []
    checked_compositions = 0

    for name in SOURCE_NAMES:
        path = HERE / name
        raw = path.read_bytes()
        receipt = json.loads(raw)
        if receipt["schema"] != "rm2sub-fixed-rm-three-group-sparse-bridge-v1":
            raise AssertionError(f"unexpected schema in {name}")
        parameters = receipt["parameters"]
        expected_parameters = {
            "message_bits": 1 << 16,
            "output_bits": 1 << 17,
            "outer_rows": 256,
            "outer_block_bits": 512,
            "step_bits": 64,
            "state_bits": 14,
            "bad_weight": 13107,
            "live_model": "support-averaged-preaddmul",
            "change_of_measure": "pointwise-density-envelope",
            "exact_zero_state_destination_split": True,
            "exact_nonactivation_verified_from_kernel_shell_counts": True,
        }
        for key, value in expected_parameters.items():
            if parameters[key] != value:
                raise AssertionError(f"parameter mismatch for {key} in {name}")
        observed_groups = tuple(
            (
                row["name"],
                int(row["lower"]),
                int(row["upper"]),
                float(row["reference_probability"]),
            )
            for row in receipt["groups"]
        )
        for observed, expected in zip(observed_groups, EXPECTED_GROUPS):
            if observed[:3] != expected[:3]:
                raise AssertionError(f"group partition mismatch in {name}")
            close(observed[3], expected[3], 2e-12)

        composition_rows = receipt["composition_rows"]
        cursor = 0
        chunk_logs = []
        for row in receipt["occupation_rows"]:
            occupation = int(row["occupation"])
            expected_count = (occupation + 1) * (occupation + 2) // 2
            selected = composition_rows[cursor : cursor + expected_count]
            if len(selected) != expected_count:
                raise AssertionError(f"incomplete composition block at Q={occupation}")
            local_index = 0
            for low in range(occupation + 1):
                for central in range(occupation - low + 1):
                    high = occupation - low - central
                    observed = selected[local_index]
                    if (
                        int(observed["low_rows"]),
                        int(observed["central_rows"]),
                        int(observed["high_rows"]),
                    ) != (low, central, high):
                        raise AssertionError(
                            f"composition ordering mismatch at Q={occupation}"
                        )
                    local_index += 1
            block = composition_rows[cursor : cursor + expected_count]
            values = np.asarray(
                [float(item["log2_upper_diagnostic"]) * LOG2 for item in block]
            )
            aggregate = float(logsumexp(values))
            close(aggregate / LOG2, float(row["log2_upper_diagnostic"]))
            occupation_rows.append(row)
            seen_occupations.append(occupation)
            chunk_logs.append(aggregate)
            cursor += expected_count
            checked_compositions += expected_count
        if cursor != len(composition_rows):
            raise AssertionError(f"unexpected trailing compositions in {name}")
        chunk_aggregate = float(logsumexp(np.asarray(chunk_logs)))
        close(
            chunk_aggregate / LOG2,
            float(receipt["interval_log2_upper_diagnostic"]),
        )
        source_rows.append(
            {
                "name": name,
                "sha256": hashlib.sha256(raw).hexdigest(),
                "occupation_interval": parameters["occupation_interval"],
                "composition_count": cursor,
                "interval_margin_bits_diagnostic": -chunk_aggregate / LOG2,
            }
        )

    if seen_occupations != list(range(30, 257)):
        raise AssertionError("source receipts do not cover Q=30 through Q=256 exactly")
    all_logs = np.asarray(
        [float(row["log2_upper_diagnostic"]) * LOG2 for row in occupation_rows]
    )
    aggregate = float(logsumexp(all_logs))
    margin = -aggregate / LOG2
    weakest = min(
        occupation_rows, key=lambda row: float(row["margin_bits_diagnostic"])
    )
    if margin <= 40.0:
        raise AssertionError("the refined-band union does not clear 40 bits")

    payload = {
        "schema": "rm2sub-refined-band-bridge-q30-q256-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "construction": (
            "one fixed RM(4,9) constituent repeated in 256 rows, uniform "
            "routing, and the fixed audited RM2Sub t=64,s=14 A/B pair"
        ),
        "probability_space": {
            "outer": "one fixed authenticated RM(4,9) constituent",
            "routing": "independent uniform row-coordinate and region permutations",
            "inner": "independent nonzero field scalar in every RM2Sub epoch",
        },
        "bad_event": "some message of outer occupation 30 through 256 encodes to weight at most 13107",
        "groups": [
            {
                "name": name,
                "lower": lower,
                "upper": upper,
                "reference_probability": probability,
            }
            for name, lower, upper, probability in EXPECTED_GROUPS
        ],
        "source_receipts": source_rows,
        "checked_composition_count": checked_compositions,
        "checked_occupation_count": len(occupation_rows),
        "occupation_rows": occupation_rows,
        "union_log2_upper_diagnostic": aggregate / LOG2,
        "union_margin_bits_diagnostic": margin,
        "union_gap_to_40_bits": margin - 40.0,
        "weakest_occupation": {
            "occupation": int(weakest["occupation"]),
            "margin_bits_diagnostic": float(weakest["margin_bits_diagnostic"]),
        },
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "The finite Chernoff grid is not asserted optimal.",
            "Occupations one through 29 are covered by separate receipts.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    audit = {
        "schema": "rm2sub-refined-band-bridge-audit-v1",
        "status": "PASS",
        "consolidated_receipt": str(OUTPUT),
        "checked_source_receipts": len(SOURCE_NAMES),
        "checked_composition_count": checked_compositions,
        "checked_occupation_count": len(occupation_rows),
        "covered_occupation_interval": [30, 256],
        "union_margin_bits_diagnostic": margin,
    }
    AUDIT_OUTPUT.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"checked_compositions={checked_compositions}")
    print(f"union_margin_bits_diagnostic={margin:.12f}")
    print(f"weakest_occupation={weakest['occupation']}")
    print(f"output={OUTPUT}")
    print(f"audit={AUDIT_OUTPUT}")


if __name__ == "__main__":
    main()
