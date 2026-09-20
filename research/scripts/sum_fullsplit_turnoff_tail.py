#!/usr/bin/env python3
"""Sum a crude multi-turnoff tail for the full-split episode ledger.

After the first active block, if H input ones remain, then at most H later
blocks are occupied and hence at most H+1 zero gaps can be selected for
turnoff/skipping.  The event that at least e_min gaps terminate is therefore
bounded by

    sum_{e>=e_min} C(H+1,e) p_term^e,

without using any survival or output-weight charge.  This is intentionally
pessimistic, but it is a theorem-shaped way to control the all-e tail after
the explicit e=0,1 ledger.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from certify_fullsplit_inner_bucket_sum import load_outer
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list
from sum_fullsplit_episode_envelope import precompute_gap_sums


def turnoff_tail_log2(*, gaps: int, e_min: int, turnoff_log2: float) -> float:
    total = float("-inf")
    for e in range(e_min, gaps + 1):
        total = log2add(total, log2_binom(gaps, e) + e * turnoff_log2)
    return total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-prefix-csv", type=Path, required=True)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--h-values", default="32:500")
    parser.add_argument("--first-r-values", default="1:64")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=1)
    parser.add_argument("--gap-stop", type=int, default=12000)
    parser.add_argument("--gap-step", type=int, default=4000)
    parser.add_argument("--turnoff-log2", type=float, required=True)
    parser.add_argument("--e-min", type=int, default=2)
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    outer = load_outer(args.outer_prefix_csv)
    gap_ranges: list[tuple[int, int]] = []
    for gap_min in range(args.gap_start, args.gap_stop + 1, args.gap_step):
        gap_max = min(args.gap_stop, gap_min + args.gap_step - 1)
        gap_ranges.append((gap_min, gap_max))

    h_values = parse_int_list(args.h_values)
    r_values = parse_int_list(args.first_r_values)
    max_h_after_first = max((h - 1 for h in h_values if h in outer), default=0)
    gap_sums = precompute_gap_sums(
        b=args.block_bits,
        late_blocks=args.late_blocks,
        gap_ranges=gap_ranges,
        max_h_after_first=max_h_after_first,
    )

    rows = []
    total = float("-inf")
    by_gap = [float("-inf") for _ in gap_ranges]
    peak = None
    tail_cache: dict[int, float] = {}
    for h in h_values:
        if h not in outer:
            continue
        for r in r_values:
            if r < 1 or r > min(h, args.block_bits):
                continue
            H = h - r
            tail = tail_cache.get(H)
            if tail is None:
                tail = turnoff_tail_log2(
                    gaps=H + 1,
                    e_min=args.e_min,
                    turnoff_log2=args.turnoff_log2,
                )
                tail_cache[H] = tail
            if tail == float("-inf"):
                continue
            for idx, (gap_min, gap_max) in enumerate(gap_ranges):
                gap_sum = gap_sums.get((H, idx), float("-inf"))
                if gap_sum == float("-inf"):
                    continue
                placement = (
                    log2_binom(args.block_bits, r)
                    + gap_sum
                    - log2_binom(args.N, h)
                )
                term = outer[h] + placement + tail
                total = log2add(total, term)
                by_gap[idx] = log2add(by_gap[idx], term)
                row = {
                    "outer_weight": h,
                    "first_r": r,
                    "remaining_ones": H,
                    "gap_min": gap_min,
                    "gap_max": gap_max,
                    "outer_log2": outer[h],
                    "placement_log2": placement,
                    "tail_log2": tail,
                    "term_log2": term,
                }
                rows.append(row)
                if peak is None or term > peak["term_log2"]:
                    peak = row

    if args.output_csv:
        fields = [
            "outer_weight",
            "first_r",
            "remaining_ones",
            "gap_min",
            "gap_max",
            "outer_log2",
            "placement_log2",
            "tail_log2",
            "term_log2",
        ]
        with args.output_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    print("Full-split multi-turnoff tail sum")
    print(f"e_min,{args.e_min}")
    print(f"turnoff_log2,{args.turnoff_log2:.6f}")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    for (gap_min, gap_max), val in zip(gap_ranges, by_gap):
        print(f"gap_{gap_min}_{gap_max}_log2,{val:.6f}")
    if peak is not None:
        print(f"peak_h,{peak['outer_weight']}")
        print(f"peak_first_r,{peak['first_r']}")
        print(f"peak_gap_min,{peak['gap_min']}")
        print(f"peak_gap_max,{peak['gap_max']}")
        print(f"peak_tail_log2,{peak['tail_log2']:.6f}")
        print(f"peak_term_log2,{peak['term_log2']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
