#!/usr/bin/env python3
"""Audit a fixed B32 activation profile against the recurrence target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


DEFAULT_FIXED = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_b32_kernel_spectrum.json"
)
DEFAULT_TARGET = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/"
    "receipts/zero_state_activation_table.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_b32_activation_domination.json"
)


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    fixed = json.loads(args.fixed.read_text(encoding="utf-8"))
    target = json.loads(args.target.read_text(encoding="utf-8"))
    fixed_by_weight = {
        int(row["total_weight"]): float(row["fixed_nonactivation_probability"])
        for row in fixed["by_total_weight"]
    }
    target_by_weight = {
        int(row["total_weight"]): float(
            row["uniform_support_average_distinct_upper_bound"]
        )
        for row in target["by_total_weight"]
    }
    rows = []
    for weight in range(1, 129):
        fixed_value = fixed_by_weight[weight]
        target_value = target_by_weight[weight]
        ratio = fixed_value / target_value if target_value else 0.0
        rows.append(
            {
                "total_weight": weight,
                "fixed_nonactivation_probability": fixed_value,
                "target_upper_bound": target_value,
                "fixed_to_target_ratio": ratio,
                "dominated": fixed_value <= target_value,
            }
        )
    worst = max(rows, key=lambda row: row["fixed_to_target_ratio"])
    return {
        "schema": "riffle-splitstate-fixed-b32-activation-domination-v1",
        "candidate": "Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul t=128 s=32",
        "all_nonzero_shells_dominated": all(row["dominated"] for row in rows),
        "worst_ratio_row": worst,
        "by_total_weight": rows,
        "consequence": (
            "Every zero-to-zero entry used by the completed recurrence is "
            "at least its exact value for this fixed B. Since all transfer "
            "entries are nonnegative, the existing end-to-end upper bounds "
            "remain valid after substituting the fixed B."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed", type=Path, default=DEFAULT_FIXED)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    if not payload["all_nonzero_shells_dominated"]:
        raise SystemExit("fixed B does not dominate the target profile")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    row = payload["worst_ratio_row"]
    print(f"all_nonzero_shells_dominated,1")
    print(
        f"worst_weight,{row['total_weight']},"
        f"ratio,{row['fixed_to_target_ratio']:.12g}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()

