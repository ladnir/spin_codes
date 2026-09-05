#!/usr/bin/env python3
"""Check endpoint safety for the full-split exact tiny-prefix grids.

The finite ledger uses exact ``e <= 8`` CSV grids for the tiny prefix and then
reuses one grid across each placement-gap bucket.  This helper checks the
current endpoint audit: for every recorded remaining weight H, the right
endpoint total must be no larger than the left endpoint total.

This is an audit guard for the CSV artifact, not a replacement for the
eventual analytic monotonicity lemma over every interior T.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class EndpointPair:
    name: str
    left: Path
    right: Path


DEFAULT_PAIRS = [
    EndpointPair(
        "gap_1_4000",
        ROOT / "fast_fullsplit_e08_T5949_H0_499_all.csv",
        ROOT / "fast_fullsplit_e08_T9948_H0_499_all.csv",
    ),
    EndpointPair(
        "gap_4001_8000",
        ROOT / "fast_fullsplit_e08_T9949_H0_499_all.csv",
        ROOT / "fast_fullsplit_e08_T13948_H0_499_all.csv",
    ),
    EndpointPair(
        "gap_8001_12000",
        ROOT / "fast_fullsplit_e08_T13949_H0_499_all.csv",
        ROOT / "fast_fullsplit_e08_T17948_H0_499_all.csv",
    ),
]


def parse_log2(text: str) -> float:
    if text == "-inf":
        return float("-inf")
    return float(text)


def load_grid(path: Path) -> dict[int, dict[str, str]]:
    with path.open(newline="") as f:
        return {int(row["H"]): row for row in csv.DictReader(f)}


def parse_pair(text: str) -> EndpointPair:
    fields = [part.strip() for part in text.split(",")]
    if len(fields) != 3:
        raise ValueError(f"bad pair {text!r}; expected name,left.csv,right.csv")
    return EndpointPair(fields[0], Path(fields[1]), Path(fields[2]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pair",
        action="append",
        help="Endpoint pair as name,left.csv,right.csv. May be repeated.",
    )
    parser.add_argument("--h-min", type=int, default=0)
    parser.add_argument("--h-max", type=int, default=499)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()

    pairs = [parse_pair(item) for item in args.pair] if args.pair else DEFAULT_PAIRS
    columns = [f"e{idx}_log2" for idx in range(9)] + ["total_log2"]

    print("Full-split exact-grid endpoint audit")
    print(
        "bucket,H_min,H_max,max_total_right_minus_left,H_total,"
        "left_total,right_total,max_column_right_minus_left,column,H_column"
    )

    failed = False
    for pair in pairs:
        left = load_grid(pair.left)
        right = load_grid(pair.right)
        h_values = list(range(args.h_min, args.h_max + 1))
        missing = [h for h in h_values if h not in left or h not in right]
        if missing:
            raise SystemExit(f"{pair.name}: missing H values, first missing H={missing[0]}")

        worst_total = (float("-inf"), -1, float("nan"), float("nan"))
        worst_column = (float("-inf"), "", -1)
        for h in h_values:
            left_total = parse_log2(left[h]["total_log2"])
            right_total = parse_log2(right[h]["total_log2"])
            total_diff = right_total - left_total
            if total_diff > worst_total[0]:
                worst_total = (total_diff, h, left_total, right_total)

            for column in columns:
                diff = parse_log2(right[h][column]) - parse_log2(left[h][column])
                if diff > worst_column[0]:
                    worst_column = (diff, column, h)

        if worst_total[0] > args.tolerance:
            failed = True

        print(
            f"{pair.name},{args.h_min},{args.h_max},"
            f"{worst_total[0]:.12g},{worst_total[1]},"
            f"{worst_total[2]:.12g},{worst_total[3]:.12g},"
            f"{worst_column[0]:.12g},{worst_column[1]},{worst_column[2]}"
        )

    if failed:
        raise SystemExit("endpoint audit failed: a right endpoint total exceeds its left endpoint")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
