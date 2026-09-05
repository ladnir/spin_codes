#!/usr/bin/env python3
"""Audit the two floating-point receipts for the packet s=256 certificate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


LOG2 = math.log(2.0)
ROOT = Path(
    "constructions/riffle_bchperm_transpose_packetshuffle_fieldcheckpoint_g4_s256"
)


def close(left: float, right: float, tolerance: float = 2e-9) -> None:
    if abs(left - right) > tolerance:
        raise AssertionError(f"receipt mismatch: {left} != {right}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--regular",
        type=Path,
        default=ROOT / "receipts/all_profiles_regular_support_relaxation_delta09.json",
    )
    parser.add_argument(
        "--combined",
        type=Path,
        default=ROOT / "receipts/all_one_completion_delta09.json",
    )
    parser.add_argument("--target-margin-bits", type=float, default=40.0)
    args = parser.parse_args()

    regular = json.loads(args.regular.read_text(encoding="utf-8"))
    combined = json.loads(args.combined.read_text(encoding="utf-8"))
    if regular["schema"] != "riffle-packet4-fieldcheckpoint-s256-regular-support-v1":
        raise AssertionError("unexpected regular receipt schema")
    if combined["schema"] != "riffle-packet4-fieldcheckpoint-s256-all-one-completion-v1":
        raise AssertionError("unexpected combined receipt schema")
    if regular["candidate"] != combined["candidate"]:
        raise AssertionError("candidate names differ")

    regular_rows = regular["occupation_rows"]
    if [row["active_regular_outer_blocks"] for row in regular_rows] != list(
        range(1, 8193)
    ):
        raise AssertionError("regular rows do not cover active counts 1 through 8192")
    regular_recomputed = float(
        logsumexp(
            np.asarray(
                [float(row["pointwise_log2_upper"]) * LOG2 for row in regular_rows]
            )
        )
        / LOG2
    )
    close(regular_recomputed, float(regular["regular_log2_upper"]))

    mixed_rows = combined["mixed_rows"]
    if [row["regular_active_blocks"] for row in mixed_rows] != list(range(8192)):
        raise AssertionError("mixed rows do not cover regular counts 0 through 8191")
    mixed_recomputed = float(
        logsumexp(
            np.asarray(
                [float(row["pointwise_log2_upper"]) * LOG2 for row in mixed_rows]
            )
        )
        / LOG2
    )
    close(mixed_recomputed, float(combined["all_one_class_log2_upper"]))
    close(float(regular["regular_log2_upper"]), float(combined["regular_class_log2_upper"]))
    total_recomputed = float(
        np.logaddexp(regular_recomputed * LOG2, mixed_recomputed * LOG2) / LOG2
    )
    close(total_recomputed, float(combined["combined_log2_upper"]))

    low_end = int(combined["parameters"]["low_end"])
    high_start = int(combined["parameters"]["high_start"])
    for row in mixed_rows:
        active = int(row["regular_active_blocks"])
        method = str(row["method"])
        if active <= low_end:
            expected = "one forced packet; all regular inputs deleted"
        elif active < high_start:
            expected = "one forced packet plus pessimistically packed regular packets"
        else:
            expected = "all all-one inputs deleted; regular bound reused"
        if method != expected:
            raise AssertionError(f"wrong mixed method at active count {active}")

    margin = -total_recomputed
    if margin < args.target_margin_bits:
        raise AssertionError(
            f"combined margin {margin:.12f} is below {args.target_margin_bits}"
        )
    print("certificate_audit,PASS")
    print(f"regular_lambda_bits,{-regular_recomputed:.12f}")
    print(f"all_one_lambda_bits,{-mixed_recomputed:.12f}")
    print(f"combined_lambda_bits,{margin:.12f}")


if __name__ == "__main__":
    main()
