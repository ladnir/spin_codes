#!/usr/bin/env python3
"""Sweep h-dependent full-split episode bounds and outer terms.

This wraps bound_fullsplit_episode_gaps.py in the first-moment shape:

    A_h^out * Pr[first active block has occupancy r in a gap bucket]
              * p_inner_episode(h-r, live_after_first).

It is intended for coarse prefix-certificate development, not final theorem
generation.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from bound_fullsplit_episode_gaps import (
    episode_gap_bound_log2,
    load_spectrum,
    precompute_survival_log2,
    split_entries,
)
from certify_fullsplit_inner_bucket_sum import first_active_bucket_r_log2, load_outer
from dense_largek_eval import log2add
from probe_block_recursive_inner import parse_int_list


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-prefix-csv", type=Path, required=True)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--h-values", required=True)
    parser.add_argument("--first-r-values", default="1")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=1)
    parser.add_argument("--gap-stop", type=int, default=12000)
    parser.add_argument("--gap-step", type=int, default=4000)
    parser.add_argument("--e-max", type=int, default=4)
    parser.add_argument("--extra-turnoff-log2", type=float)
    parser.add_argument("--lambda-min", type=float, default=0.001)
    parser.add_argument("--lambda-max", type=float, default=2.0)
    parser.add_argument("--lambda-step", type=float, default=0.002)
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    d = math.floor(args.distance_delta * args.N)
    b = args.block_bits
    outer = load_outer(args.outer_prefix_csv)
    entries = split_entries(load_spectrum(args.spectrum), b)
    p0 = sum(p for _j, q, p in entries if q == 0)
    turnoff_log2 = math.log2(p0)
    if args.extra_turnoff_log2 is not None:
        turnoff_log2 = log2add(turnoff_log2, args.extra_turnoff_log2)

    gap_ranges: list[tuple[int, int, int]] = []
    for gap_min in range(args.gap_start, args.gap_stop + 1, args.gap_step):
        gap_max = min(args.gap_stop, gap_min + args.gap_step - 1)
        live_after = args.late_blocks + gap_min - 1
        gap_ranges.append((gap_min, gap_max, live_after))

    survival_by_live = {
        live_after: precompute_survival_log2(
            entries=entries,
            max_live=live_after,
            distance=d,
            lambda_min=args.lambda_min,
            lambda_max=args.lambda_max,
            lambda_step=args.lambda_step,
        )
        for _gap_min, _gap_max, live_after in gap_ranges
    }

    rows = []
    total = float("-inf")
    peak = None
    for h in parse_int_list(args.h_values):
        if h not in outer:
            continue
        for r in parse_int_list(args.first_r_values):
            if r < 1 or r > h:
                continue
            H = h - r
            for gap_min, gap_max, live_after in gap_ranges:
                placement = first_active_bucket_r_log2(
                    n=args.N,
                    b=b,
                    h=h,
                    late_blocks=args.late_blocks,
                    gap_min=gap_min,
                    gap_max=gap_max,
                    first_r=r,
                )
                if placement == float("-inf"):
                    continue
                inner, inner_peak = episode_gap_bound_log2(
                    b=b,
                    remaining_blocks=live_after,
                    remaining_ones=H,
                    p0_log2=turnoff_log2,
                    survival_log2=survival_by_live[live_after],
                    e_max=args.e_max,
                )
                term = outer[h] + placement + min(0.0, inner)
                total = log2add(total, term)
                row = {
                    "h": h,
                    "first_r": r,
                    "gap_min": gap_min,
                    "gap_max": gap_max,
                    "live_after_first": live_after,
                    "outer_log2": outer[h],
                    "placement_log2": placement,
                    "inner_log2": inner,
                    "term_log2": term,
                    "turnoff_log2": turnoff_log2,
                    "inner_peak_occupied": "" if inner_peak is None else inner_peak[0],
                    "inner_peak_selected": "" if inner_peak is None else inner_peak[1],
                    "inner_peak_skipped": "" if inner_peak is None else inner_peak[2],
                    "inner_peak_term_log2": "" if inner_peak is None else inner_peak[3],
                }
                rows.append(row)
                if peak is None or term > peak["term_log2"]:
                    peak = row

    fields = [
        "h",
        "first_r",
        "gap_min",
        "gap_max",
        "live_after_first",
        "outer_log2",
        "placement_log2",
        "inner_log2",
        "term_log2",
        "turnoff_log2",
        "inner_peak_occupied",
        "inner_peak_selected",
        "inner_peak_skipped",
        "inner_peak_term_log2",
    ]
    if args.output_csv:
        with args.output_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    print("Full-split episode first-moment sweep")
    print(f"N={args.N}, b={b}, d={d}, e_max={args.e_max}, turnoff_log2={turnoff_log2:.6f}")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    if peak is not None:
        print(f"peak_h,{peak['h']}")
        print(f"peak_first_r,{peak['first_r']}")
        print(f"peak_gap_min,{peak['gap_min']}")
        print(f"peak_gap_max,{peak['gap_max']}")
        print(f"peak_term_log2,{peak['term_log2']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
