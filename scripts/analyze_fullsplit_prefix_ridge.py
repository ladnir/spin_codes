#!/usr/bin/env python3
"""Summarize the tiny-prefix ridge in the full-split finite ledger.

The finite ``h=32..500`` prefix ledger is currently the tightest part of the
RM/full-split certificate.  This helper keeps the bottleneck decomposition
visible: which gap bucket, first-block occupancy, and outer weight range carry
the mass.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def log2add(a: float, b: float) -> float:
    if a == float("-inf"):
        return b
    if b == float("-inf"):
        return a
    if a < b:
        a, b = b, a
    return a + math.log2(1.0 + 2.0 ** (b - a))


def log2sum(values: list[float]) -> float:
    total = float("-inf")
    for value in values:
        total = log2add(total, value)
    return total


def share(log2_value: float, total_log2: float) -> float:
    if log2_value == float("-inf"):
        return 0.0
    return 2.0 ** (log2_value - total_log2)


def parse_h_cutoffs(text: str) -> list[int]:
    if not text:
        return []
    return [int(part) for part in text.split(",") if part]


def row_matches_ridge(row: dict[str, str | int | float], *, gap: str, r: int) -> bool:
    return int(row["_r"]) == r and str(row["_gap"]) == gap


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv",
        type=Path,
        default=ROOT / "fullsplit_piecewise_h32_500_csv.csv",
        help="Prefix row ledger produced by sum_fullsplit_piecewise_certificate.py.",
    )
    parser.add_argument(
        "--expected-total-log2",
        type=float,
        default=None,
        help="Optional expected log2 total; fail if the audited CSV disagrees.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=5e-7,
        help="Tolerance for --expected-total-log2.",
    )
    parser.add_argument("--top-rows", type=int, default=8)
    parser.add_argument("--top-h", type=int, default=12)
    parser.add_argument("--top-hr", type=int, default=8)
    parser.add_argument(
        "--h-cutoffs",
        default="64,80,96,128,160",
        help="Comma-separated cumulative h cutoffs to report.",
    )
    parser.add_argument("--ridge-gap", default="1-4000")
    parser.add_argument("--ridge-r", type=int, default=1)
    args = parser.parse_args()

    rows: list[dict[str, str | int | float]] = []
    with args.csv.open(newline="") as f:
        for row in csv.DictReader(f):
            h = int(row["outer_weight"])
            r = int(row["first_r"])
            gap_min = int(row["gap_min"])
            gap_max = int(row["gap_max"])
            term = float(row["term_log2"])
            row["_h"] = h
            row["_r"] = r
            row["_gap"] = f"{gap_min}-{gap_max}"
            row["_term"] = term
            rows.append(row)

    total = log2sum([float(row["_term"]) for row in rows])
    if args.expected_total_log2 is not None:
        diff = abs(total - args.expected_total_log2)
        if diff > args.tolerance:
            raise SystemExit(
                f"{args.csv}: total {total:.12f} differs from expected "
                f"{args.expected_total_log2:.12f} by {diff:.3g}"
            )

    by_gap: defaultdict[str, float] = defaultdict(lambda: float("-inf"))
    by_r: defaultdict[int, float] = defaultdict(lambda: float("-inf"))
    by_h: defaultdict[int, float] = defaultdict(lambda: float("-inf"))
    by_hr: defaultdict[tuple[int, int], float] = defaultdict(lambda: float("-inf"))

    for row in rows:
        h = int(row["_h"])
        r = int(row["_r"])
        gap = str(row["_gap"])
        term = float(row["_term"])
        by_gap[gap] = log2add(by_gap[gap], term)
        by_r[r] = log2add(by_r[r], term)
        by_h[h] = log2add(by_h[h], term)
        by_hr[(h, r)] = log2add(by_hr[(h, r)], term)

    sorted_rows = sorted(rows, key=lambda row: float(row["_term"]), reverse=True)
    peak = sorted_rows[0]
    print(f"prefix_ridge_rows,{len(rows)}")
    print(f"prefix_ridge_total_log2,{total:.9f}")
    print(
        "prefix_ridge_peak,"
        f"h={peak['_h']},r={peak['_r']},gap={peak['_gap']},"
        f"term_log2={float(peak['_term']):.9f}"
    )
    print(f"prefix_ridge_peak_to_total_bits,{total - float(peak['_term']):.9f}")

    print("prefix_ridge_gap_totals")
    for gap, value in sorted(by_gap.items()):
        print(f"gap,{gap},log2,{value:.9f},share,{share(value, total):.12g}")

    print("prefix_ridge_r_totals")
    for r, value in sorted(by_r.items(), key=lambda item: item[1], reverse=True)[: args.top_hr]:
        print(f"r,{r},log2,{value:.9f},share,{share(value, total):.12g}")

    print("prefix_ridge_h_totals")
    for h, value in sorted(by_h.items(), key=lambda item: item[1], reverse=True)[: args.top_h]:
        print(f"h,{h},log2,{value:.9f},share,{share(value, total):.12g}")

    print("prefix_ridge_hr_totals")
    for (h, r), value in sorted(by_hr.items(), key=lambda item: item[1], reverse=True)[: args.top_hr]:
        print(f"h,{h},r,{r},log2,{value:.9f},share,{share(value, total):.12g}")

    print("prefix_ridge_h_cumulative")
    cumulative = float("-inf")
    cutoffs = set(parse_h_cutoffs(args.h_cutoffs))
    for h in sorted(by_h):
        cumulative = log2add(cumulative, by_h[h])
        if h in cutoffs:
            print(f"h_le,{h},log2,{cumulative:.9f},share,{share(cumulative, total):.12g}")

    print("prefix_ridge_adjacent_ratios_r1")
    ridge_rows = sorted(
        [row for row in rows if row_matches_ridge(row, gap=args.ridge_gap, r=args.ridge_r)],
        key=lambda row: int(row["_h"]),
    )
    previous_h: int | None = None
    previous_value: float | None = None
    max_ratio_all = (-1.0, None, None)
    max_ratio_after_first = (-1.0, None, None)
    for row in ridge_rows:
        h = int(row["_h"])
        value = float(row["_term"])
        if previous_h is not None and previous_value is not None:
            ratio = 2.0 ** (value - previous_value)
            if ratio > max_ratio_all[0]:
                max_ratio_all = (ratio, previous_h, h)
            if previous_h >= int(ridge_rows[1]["_h"]) and ratio > max_ratio_after_first[0]:
                max_ratio_after_first = (ratio, previous_h, h)
        previous_h = h
        previous_value = value
    ratio, h0, h1 = max_ratio_all
    print(f"max_ratio_all,{ratio:.12g},from,{h0},to,{h1}")
    ratio, h0, h1 = max_ratio_after_first
    print(f"max_ratio_after_33,{ratio:.12g},from,{h0},to,{h1}")

    if len(ridge_rows) >= 2:
        ridge_exact = log2sum([float(row["_term"]) for row in ridge_rows])
        remainder_exact = log2sum(
            [
                float(row["_term"])
                for row in rows
                if not row_matches_ridge(row, gap=args.ridge_gap, r=args.ridge_r)
            ]
        )
        ridge_peak = float(ridge_rows[0]["_term"])
        first_ratio = 2.0 ** (float(ridge_rows[1]["_term"]) - ridge_peak)
        tail_ratio = max_ratio_after_first[0]
        finite_factor = 1.0 + first_ratio * (1.0 - tail_ratio ** (len(ridge_rows) - 1)) / (1.0 - tail_ratio)
        infinite_factor = 1.0 + first_ratio / (1.0 - tail_ratio)
        finite_geometric = ridge_peak + math.log2(finite_factor)
        infinite_geometric = ridge_peak + math.log2(infinite_factor)
        finite_total_bound = log2add(finite_geometric, remainder_exact)
        infinite_total_bound = log2add(infinite_geometric, remainder_exact)
        inner_values = [float(row["inner_log2"]) for row in ridge_rows]
        print("prefix_ridge_ratio_certificate")
        print(
            f"ridge_subset,r,{args.ridge_r},gap,{args.ridge_gap},"
            f"h_min,{ridge_rows[0]['_h']},h_max,{ridge_rows[-1]['_h']},rows,{len(ridge_rows)}"
        )
        print(f"ridge_exact_log2,{ridge_exact:.9f}")
        print(f"ridge_remainder_exact_log2,{remainder_exact:.9f}")
        print(f"ridge_first_ratio,{first_ratio:.12g}")
        print(f"ridge_tail_ratio,{tail_ratio:.12g}")
        print(f"ridge_geometric_finite_log2,{finite_geometric:.9f}")
        print(f"ridge_geometric_infinite_log2,{infinite_geometric:.9f}")
        print(f"ridge_geometric_slack_bits,{finite_geometric - ridge_exact:.12g}")
        print(f"ridge_total_ratio_bound_log2,{finite_total_bound:.9f}")
        print(f"ridge_total_ratio_bound_slack_bits,{finite_total_bound - total:.12g}")
        print(f"ridge_total_infinite_bound_log2,{infinite_total_bound:.9f}")
        print(f"ridge_inner_log2_span,{max(inner_values) - min(inner_values):.12g}")

    print("prefix_ridge_top_rows")
    for index, row in enumerate(sorted_rows[: args.top_rows], start=1):
        print(
            f"row,{index},h,{row['_h']},r,{row['_r']},gap,{row['_gap']},"
            f"term_log2,{float(row['_term']):.9f},"
            f"share,{share(float(row['_term']), total):.12g},"
            f"inner_log2,{float(row['inner_log2']):.9f},"
            f"placement_log2,{float(row['placement_log2']):.9f}"
        )


if __name__ == "__main__":
    main()
