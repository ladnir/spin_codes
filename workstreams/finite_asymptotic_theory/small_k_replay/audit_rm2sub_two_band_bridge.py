#!/usr/bin/env python3
"""Audit structure and aggregation of the consolidated two-band receipt."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
RECEIPT = (
    HERE
    / "rm2sub_two_band_bridge_all_partners_holder_split_reference_probe_d100.json"
)
COLLAPSE_AUDIT = HERE / "rm2sub_two_band_collapse_audit.json"
OUTPUT = HERE / "rm2sub_two_band_bridge_audit.json"
LOG2 = math.log(2.0)


def close(first: float, second: float, tolerance: float = 2e-9) -> None:
    if not math.isclose(first, second, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"aggregation mismatch: {first} != {second}")


def main() -> None:
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    collapse = json.loads(COLLAPSE_AUDIT.read_text(encoding="utf-8"))
    if receipt["schema"] != "rm2sub-fixed-rm-two-band-sparse-bridge-v1":
        raise AssertionError("unexpected bridge schema")
    if collapse["status"] != "PASS":
        raise AssertionError("collapse audit did not pass")
    parameters = receipt["parameters"]
    expected = {
        "message_bits": 1 << 16,
        "output_bits": 1 << 17,
        "outer_rows": 256,
        "outer_block_bits": 512,
        "step_bits": 64,
        "state_bits": 14,
        "bad_weight": 13107,
        "occupation_interval": [30, 52],
        "live_model": "support-averaged-preaddmul",
        "change_of_measure": "holder",
        "exact_zero_state_destination_split": True,
        "exact_nonactivation_verified_from_kernel_shell_counts": True,
    }
    for key, value in expected.items():
        if parameters[key] != value:
            raise AssertionError(f"parameter mismatch for {key}")
    selected = receipt["selected_bands"]
    if [int(row["index"]) for row in selected] != list(range(9)):
        raise AssertionError("selected band indexes are incomplete")
    if [float(row["probability"]) for row in selected] != [
        0.25,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.75,
    ]:
        raise AssertionError("reference probability schedule changed")
    if [(int(row["lower"]), int(row["upper"])) for row in selected] != [
        (1, 47),
        (48, 63),
        (64, 95),
        (96, 159),
        (160, 191),
        (192, 320),
        (321, 352),
        (353, 416),
        (417, 512),
    ]:
        raise AssertionError("band partition changed")

    global_values = []
    pair_margins = []
    for pair in receipt["pair_rows"]:
        if [int(row["occupation"]) for row in pair["occupation_rows"]] != list(
            range(30, 53)
        ):
            raise AssertionError("occupation interval is incomplete")
        occupation_values = []
        for occupation in pair["occupation_rows"]:
            mixed = [
                float(row["log2_upper_diagnostic"]) * LOG2
                for row in occupation["split_rows"]
                if row["is_mixed"]
            ]
            aggregate = float(logsumexp(np.asarray(mixed)))
            close(
                aggregate / LOG2,
                float(occupation["mixed_log2_upper_diagnostic"]),
            )
            occupation_values.append(aggregate)
            global_values.append(aggregate)
        pair_aggregate = float(logsumexp(np.asarray(occupation_values)))
        close(
            pair_aggregate / LOG2,
            float(pair["mixed_interval_log2_upper_diagnostic"]),
        )
        margin = -pair_aggregate / LOG2
        pair_margins.append(margin)
        if margin <= 40.0:
            raise AssertionError("a displayed two-band family fell below 40 bits")
    global_aggregate = float(logsumexp(np.asarray(global_values)))
    close(
        global_aggregate / LOG2,
        float(
            receipt["union_of_displayed_mixed_pairs"]["log2_upper_diagnostic"]
        ),
    )
    global_margin = -global_aggregate / LOG2
    close(global_margin, 1172.8820317434236)
    if min(pair_margins) != pair_margins[0]:
        raise AssertionError("the recorded weakest partner changed")

    payload = {
        "schema": "rm2sub-two-band-bridge-audit-v1",
        "status": "PASS",
        "source_receipt": str(RECEIPT),
        "collapse_audit": str(COLLAPSE_AUDIT),
        "checked_pair_count": len(pair_margins),
        "checked_occupations_per_pair": 23,
        "union_margin_bits_diagnostic": global_margin,
        "minimum_pair_margin_bits_diagnostic": min(pair_margins),
        "collapse_maximum_absolute_error": collapse["maximum_absolute_error"],
        "homogeneous_reduction_maximum_absolute_error": collapse[
            "homogeneous_reduction_maximum_absolute_error"
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"union_margin_bits_diagnostic={global_margin:.12f}")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
