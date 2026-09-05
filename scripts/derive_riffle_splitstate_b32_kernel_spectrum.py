#!/usr/bin/env python3
"""Apply the exact MacWilliams transform to a fixed-B dual spectrum."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


DEFAULT_DUAL = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_b32_dual_histogram.csv"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_b32_kernel_spectrum.json"
)


def krawtchouk(length: int, degree: int, point: int) -> int:
    result = 0
    for intersection in range(
        max(0, degree - (length - point)), min(degree, point) + 1
    ):
        term = math.comb(point, intersection) * math.comb(
            length - point, degree - intersection
        )
        result += -term if intersection & 1 else term
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dual = [0] * (args.length + 1)
    with args.dual.open(newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            dual[int(row["weight"])] = int(row["count"])
    dual_mass = sum(dual)
    expected_dual_mass = 1 << args.redundancy
    if dual_mass != expected_dual_mass:
        raise ArithmeticError("dual spectrum mass mismatch")

    kernel = []
    for weight in range(args.length + 1):
        numerator = sum(
            count * krawtchouk(args.length, weight, dual_weight)
            for dual_weight, count in enumerate(dual)
        )
        if numerator % expected_dual_mass:
            raise ArithmeticError(
                f"nonintegral MacWilliams coefficient at weight {weight}"
            )
        kernel.append(numerator // expected_dual_mass)
    expected_kernel_mass = 1 << (args.length - args.redundancy)
    if sum(kernel) != expected_kernel_mass:
        raise ArithmeticError("kernel spectrum mass mismatch")
    if any(value < 0 for value in kernel):
        raise ArithmeticError("negative kernel spectrum coefficient")

    rows = []
    for weight, count in enumerate(kernel):
        shell = math.comb(args.length, weight)
        rows.append(
            {
                "total_weight": weight,
                "kernel_words": count,
                "shell_size": shell,
                "fixed_nonactivation_probability": count / shell,
                "uniform_support_average_distinct_upper_bound": count / shell,
            }
        )
    return {
        "schema": "riffle-splitstate-fixed-b-kernel-spectrum-v1",
        "candidate": args.candidate,
        "parameters": {
            "length": args.length,
            "redundancy": args.redundancy,
            "kernel_dimension": args.length - args.redundancy,
        },
        "checks": {
            "dual_mass": dual_mass,
            "kernel_mass": sum(kernel),
            "minimum_kernel_distance": next(
                weight for weight in range(1, len(kernel)) if kernel[weight]
            ),
            "complement_symmetric": all(
                kernel[weight] == kernel[args.length - weight]
                for weight in range(args.length + 1)
            ),
        },
        "by_total_weight": rows,
        "scope": (
            "Exact integer MacWilliams transform of the recorded complete "
            "dual histogram. Each probability is the exact fraction of one "
            "weight shell contained in the kernel of the fixed B."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dual", type=Path, default=DEFAULT_DUAL)
    parser.add_argument("--length", type=int, default=128)
    parser.add_argument("--redundancy", type=int, default=32)
    parser.add_argument(
        "--candidate",
        default="Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul t=128 s=32",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"minimum_kernel_distance,{payload['checks']['minimum_kernel_distance']}"
    )
    for weight in (4, 6, 8, 16, 32, 64):
        row = payload["by_total_weight"][weight]
        print(
            f"weight,{weight},kernel_words,{row['kernel_words']},"
            f"probability,{row['fixed_nonactivation_probability']:.12g}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
