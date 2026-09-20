#!/usr/bin/env python3
"""Audit a complete fixed-A spectrum and its coordinate moment identities."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


DEFAULT_SELECTION = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_a32_selection.json"
)
DEFAULT_HISTOGRAM = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_a32_histogram.csv"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_a32_spectrum_audit.json"
)


def falling(value: int, degree: int) -> int:
    result = 1
    for offset in range(degree):
        result *= value - offset
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    selection = json.loads(args.selection.read_text(encoding="utf-8"))
    spectrum = [0] * (args.step_bits + 1)
    with args.histogram.open(newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            spectrum[int(row["weight"])] = int(row["count"])
    mass = sum(spectrum)
    if mass != 1 << args.state_bits:
        raise ArithmeticError("fixed A spectrum mass mismatch")
    minimum_distance = next(
        weight for weight in range(1, args.step_bits + 1) if spectrum[weight]
    )
    moment_rows = []
    for degree in (1, 2, 3):
        observed = sum(
            count * falling(weight, degree)
            for weight, count in enumerate(spectrum)
        )
        expected = falling(args.step_bits, degree) * (
            1 << (args.state_bits - degree)
        )
        moment_rows.append(
            {
                "degree": degree,
                "observed_including_zero": observed,
                "expected_from_coordinate_independence": expected,
                "matches": observed == expected,
            }
        )
    coordinate = selection["selected"]["coordinate_audit"]
    return {
        "schema": "riffle-splitstate-fixed-a-spectrum-audit-v1",
        "candidate": args.candidate,
        "checks": {
            "spectrum_mass": mass,
            "minimum_distance": minimum_distance,
            "ba_zero": bool(selection["selected"]["ba_zero"]),
            "zero_coordinate_forms": coordinate["zero_coordinate_forms"],
            "distinct_coordinate_forms": coordinate["distinct_coordinate_forms"],
            "triple_dependencies": coordinate["triple_dependencies"],
            "first_three_factorial_moments_match": all(
                row["matches"] for row in moment_rows
            ),
        },
        "factorial_moments": moment_rows,
        "spectrum": [
            {"weight": weight, "count": count}
            for weight, count in enumerate(spectrum)
            if count
        ],
        "scope": (
            f"Exact audit of the complete 2^{args.state_bits}-word histogram "
            "produced by the Gray-code enumerator."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--histogram", type=Path, default=DEFAULT_HISTOGRAM)
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--state-bits", type=int, default=32)
    parser.add_argument(
        "--candidate",
        default="Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul t=128 s=32",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    checks = payload["checks"]
    if not (
        checks["ba_zero"]
        and checks["minimum_distance"] >= 20
        and checks["zero_coordinate_forms"] == 0
        and checks["distinct_coordinate_forms"] == args.step_bits
        and checks["triple_dependencies"] == 0
        and checks["first_three_factorial_moments_match"]
    ):
        raise SystemExit("fixed A failed the required proof profile")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"minimum_distance,{checks['minimum_distance']}")
    print(f"ba_zero,{int(checks['ba_zero'])}")
    print("coordinate_profile,1")
    print("factorial_moments,1")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
