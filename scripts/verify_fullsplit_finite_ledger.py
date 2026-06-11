#!/usr/bin/env python3
"""Verify the current finite full-split dense-inner first-moment ledger.

This is a manifest-style checker for the finite ``N=2^21, delta=.09`` RM/full
split certificate recorded in innerDense.tex.  The tiny prefix is still backed
by explicit finite episode CSVs, while the post-prefix and high tail use the
theorem-facing all-episode wrapper.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent


HIGH_INTERVALS = [
    ("2001--7858", 0.02, 0.01, -588.081877),
    ("7859--20550", 0.05, 0.03, -2657.429602),
    ("20551--75000", 0.2, 0.1, -6027.058821),
    ("75001--250000", 0.5, 0.3, -12846.353611),
    ("250001--350000", 1.2, 0.5, -28478.125315),
    ("350001--400000", 1.2, 1.0, -81557.079312),
    ("400001--450000", 1.2, 1.0, -91423.257760),
    ("450001--550000", 1.2, 1.0, -82083.446102),
    ("550001--650000", 1.2, 1.0, -1519.978235),
    ("650001--725000", 1.2, 2.0, -116329.618720),
    ("725001--750000", 1.2, 3.0, -169889.387695),
    ("750001--850000", 1.2, 3.0, -134525.068420),
    ("850001--950000", 1.2, 3.0, -117536.328304),
    ("950001--1050000", 1.2, 3.0, -285054.035442),
    ("1050001--1148736", 1.2, 3.0, -562862.201412),
]


@dataclass(frozen=True)
class CsvLedger:
    name: str
    path: Path
    column: str
    expected_log2: float


def log2add(a: float, b: float) -> float:
    if a == float("-inf"):
        return b
    if b == float("-inf"):
        return a
    if a < b:
        a, b = b, a
    return a + math.log2(1.0 + 2.0 ** (b - a))


def csv_logsum(path: Path, column: str) -> tuple[float, int, dict[str, str] | None]:
    total = float("-inf")
    rows = 0
    peak: dict[str, str] | None = None
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows += 1
            value = float(row[column])
            total = log2add(total, value)
            if peak is None or value > float(peak[column]):
                peak = row
    return total, rows, peak


def check_close(name: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise SystemExit(
            f"{name}: expected {expected:.6f}, got {actual:.6f}; "
            f"delta={actual - expected:.6g}"
        )


def print_peak(prefix: str, peak: dict[str, str] | None, column: str) -> None:
    if peak is None:
        return
    keys = [
        "outer_weight",
        "first_r",
        "remaining_ones",
        "gap_min",
        "gap_max",
        "inner_mode",
        column,
    ]
    fields = [f"{key}={peak[key]}" for key in keys if key in peak]
    print(f"{prefix}_peak," + ",".join(fields))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prefix-e-le8-csv",
        type=Path,
        default=ROOT / "fullsplit_piecewise_h32_500_csv.csv",
    )
    parser.add_argument(
        "--prefix-e-ge9-tail-csv",
        type=Path,
        default=ROOT / "fullsplit_turnoff_tail_h32_500_r1_64_emin9.csv",
    )
    parser.add_argument(
        "--postprefix-eall-hsummary-csv",
        type=Path,
        default=ROOT / "fullsplit_piecewise_h501_2000_eall_hsummary.csv",
    )
    parser.add_argument("--tolerance", type=float, default=5e-6)
    args = parser.parse_args()

    csv_ledgers = [
        CsvLedger(
            "prefix_32_500_e_le8",
            args.prefix_e_le8_csv,
            "term_log2",
            -34.767174,
        ),
        CsvLedger(
            "prefix_32_500_e_ge9_tail",
            args.prefix_e_ge9_tail_csv,
            "term_log2",
            -269.335258,
        ),
        CsvLedger(
            "postprefix_501_2000_eall",
            args.postprefix_eall_hsummary_csv,
            "aggregate_log2",
            -182.259739,
        ),
    ]

    print("Full-split finite ledger")
    parts: dict[str, float] = {}
    for item in csv_ledgers:
        total, rows, peak = csv_logsum(item.path, item.column)
        check_close(item.name, total, item.expected_log2, args.tolerance)
        parts[item.name] = total
        print(f"{item.name}_rows,{rows}")
        print(f"{item.name}_log2,{total:.6f}")
        print_peak(item.name, peak, item.column)

    prefix_all = log2add(parts["prefix_32_500_e_le8"], parts["prefix_32_500_e_ge9_tail"])
    print(f"prefix_32_500_all_e_log2,{prefix_all:.6f}")

    high_total = float("-inf")
    for label, lam, rho, value in HIGH_INTERVALS:
        high_total = log2add(high_total, value)
        print(f"interval_{label}_lambda,{lam:g}")
        print(f"interval_{label}_rho,{rho:g}")
        print(f"interval_{label}_log2,{value:.6f}")
    print(f"interval_2001_1148736_total_log2,{high_total:.6f}")
    print("feasible_far_bucket_cutoff_h,1148736")

    total = float("-inf")
    for value in [prefix_all, parts["postprefix_501_2000_eall"], high_total]:
        total = log2add(total, value)
    print(f"finite_ledger_total_log2,{total:.6f}")
    print(f"finite_ledger_margin_bits,{-total:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
