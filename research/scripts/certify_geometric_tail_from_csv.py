#!/usr/bin/env python3
"""Post-process a first-moment CSV with a conditional geometric tail certificate.

The global episode-cover evaluator produces pointwise log2 contributions

    t_h = A_h^out p_h^in(delta).

This helper checks an observed finite window for adjacent decay and computes the
residual that would follow from a proved geometric-ratio lemma beyond the last
computed h.  It deliberately keeps the proof obligation explicit: the script can
verify the ratio on rows present in the CSV, but a paper proof must still supply
the ratio for the uncomputed tail.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def log2add(x: float, y: float) -> float:
    if x == float("-inf"):
        return y
    if y == float("-inf"):
        return x
    hi = max(x, y)
    lo = min(x, y)
    return hi + math.log2(1.0 + 2.0 ** (lo - hi))


def format_log2(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.6f}"


def parse_log2(text: str) -> float:
    text = text.strip()
    if text == "-inf":
        return float("-inf")
    return float(text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, help="Pointwise CSV from certify_global_episode_cover.py")
    parser.add_argument("--sigma", type=int, default=None, help="Optional sigma filter.")
    parser.add_argument(
        "--check-from",
        type=int,
        required=True,
        help="First h where adjacent ratios are checked on the available CSV rows.",
    )
    parser.add_argument(
        "--tail-after",
        type=int,
        default=None,
        help="Last computed h; residual starts at h=tail_after+1. Defaults to the largest checked h.",
    )
    parser.add_argument(
        "--assume-log2-ratio",
        type=float,
        default=None,
        help="Use this log2 adjacent ratio for the uncomputed residual. Defaults to the observed worst ratio.",
    )
    parser.add_argument("--out", default=None, help="Optional one-row summary CSV.")
    args = parser.parse_args()

    rows: list[tuple[int, float]] = []
    with open(args.csv, newline="") as f:
        for row in csv.DictReader(f):
            if args.sigma is not None and int(row["sigma"]) != args.sigma:
                continue
            rows.append((int(row["h"]), parse_log2(row["term_log2"])))
    if not rows:
        raise SystemExit("no rows matched")
    rows.sort()
    terms = dict(rows)
    h_min = rows[0][0]
    h_max = rows[-1][0]
    tail_after = args.tail_after if args.tail_after is not None else h_max
    if tail_after > h_max:
        raise SystemExit("--tail-after exceeds largest h in CSV")
    if args.check_from >= h_max:
        raise SystemExit("--check-from must be smaller than largest h in CSV")

    observed_worst = float("-inf")
    observed_at = -1
    checked_pairs = 0
    prev_h = None
    prev_term = None
    for h, term in rows:
        if h < args.check_from:
            continue
        if prev_h is not None and h == prev_h + 1:
            diff = term - prev_term
            checked_pairs += 1
            if diff > observed_worst:
                observed_worst = diff
                observed_at = prev_h
        prev_h = h
        prev_term = term
    if checked_pairs == 0:
        raise SystemExit("no adjacent pairs checked")

    assumed = observed_worst if args.assume_log2_ratio is None else args.assume_log2_ratio
    if assumed >= 0.0:
        residual = float("inf")
        certified_total = float("inf")
    else:
        q = 2.0 ** assumed
        last = terms[tail_after]
        residual = last + assumed - math.log2(1.0 - q)
        prefix = float("-inf")
        for h, term in rows:
            if h <= tail_after:
                prefix = log2add(prefix, term)
        certified_total = log2add(prefix, residual)

    row = {
        "source_csv": str(Path(args.csv)),
        "sigma": "" if args.sigma is None else args.sigma,
        "h_min": h_min,
        "h_max": h_max,
        "check_from": args.check_from,
        "tail_after": tail_after,
        "observed_worst_log2_ratio": format_log2(observed_worst),
        "observed_worst_at_h": observed_at,
        "assumed_log2_ratio": format_log2(assumed),
        "tail_first_h": tail_after + 1,
        "tail_residual_log2": format_log2(residual),
        "certified_total_log2_if_ratio_holds": format_log2(certified_total),
    }

    for key, value in row.items():
        print(f"{key}={value}")

    if args.out is not None:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
