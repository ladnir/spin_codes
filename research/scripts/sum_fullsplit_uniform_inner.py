#!/usr/bin/env python3
"""Sum outer/placement terms using supplied full-split inner bucket bounds.

This is a bookkeeping helper for the full-codeword random-split inner.  It does
not run the BCH transfer.  Instead it takes one inner log2 bound per gap bucket
and sums

    A_h^out * Pr[first active block has occupancy r in bucket] * p_inner(bucket)

over the requested outer prefix.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from certify_fullsplit_inner_bucket_sum import first_active_bucket_r_log2, load_outer
from dense_largek_eval import log2_binom, log2add


def parse_float_list(text: str) -> list[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def log2sub(log_a: float, log_b: float) -> float:
    if log_b == float("-inf"):
        return log_a
    if log_a < log_b:
        raise ValueError("log2sub requires a >= b")
    if log_a == log_b:
        return float("-inf")
    return log_a + math.log2(-math.expm1((log_b - log_a) * math.log(2.0)))


def first_active_bucket_any_log2(
    *,
    n: int,
    b: int,
    h: int,
    late_blocks: int,
    gap_min: int,
    gap_max: int,
) -> float:
    """Probability that the first active block is anywhere in the bucket."""

    hi_coords = b * (late_blocks + gap_max)
    lo_coords = b * (late_blocks + gap_min - 1)
    if h > hi_coords:
        return float("-inf")
    hi = log2_binom(hi_coords, h)
    lo = log2_binom(lo_coords, h) if h <= lo_coords else float("-inf")
    return log2sub(hi, lo) - log2_binom(n, h)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-prefix-csv", type=Path, required=True)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--h-max", type=int, default=500)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=1)
    parser.add_argument("--gap-stop", type=int, default=12000)
    parser.add_argument("--gap-step", type=int, default=4000)
    parser.add_argument(
        "--inner-log2-by-gap",
        required=True,
        help="Comma-separated inner log2 bounds, one per gap bucket.",
    )
    parser.add_argument("--r-max", type=int, default=64)
    args = parser.parse_args()

    outer = load_outer(args.outer_prefix_csv)
    inner_by_gap = parse_float_list(args.inner_log2_by_gap)
    gap_ranges: list[tuple[int, int]] = []
    for gap_min in range(args.gap_start, args.gap_stop + 1, args.gap_step):
        gap_max = min(args.gap_stop, gap_min + args.gap_step - 1)
        gap_ranges.append((gap_min, gap_max))
    if len(inner_by_gap) != len(gap_ranges):
        raise ValueError("inner-log2-by-gap count must match the number of gap buckets")

    total = float("-inf")
    by_gap = [float("-inf") for _ in gap_ranges]
    peak: tuple[int, int, int, int, float] | None = None

    sum_all_r = args.r_max >= args.block_bits
    for h, outer_log2 in outer.items():
        if h > args.h_max:
            continue
        if sum_all_r:
            for idx, (gap_min, gap_max) in enumerate(gap_ranges):
                placement = first_active_bucket_any_log2(
                    n=args.N,
                    b=args.block_bits,
                    h=h,
                    late_blocks=args.late_blocks,
                    gap_min=gap_min,
                    gap_max=gap_max,
                )
                if placement == float("-inf"):
                    continue
                inner = min(0.0, inner_by_gap[idx])
                term = outer_log2 + placement + inner
                total = log2add(total, term)
                by_gap[idx] = log2add(by_gap[idx], term)
                if peak is None or term > peak[4]:
                    peak = (h, 0, gap_min, gap_max, term)
            continue
        for r in range(1, min(args.r_max, args.block_bits, h) + 1):
            for idx, (gap_min, gap_max) in enumerate(gap_ranges):
                placement = first_active_bucket_r_log2(
                    n=args.N,
                    b=args.block_bits,
                    h=h,
                    late_blocks=args.late_blocks,
                    gap_min=gap_min,
                    gap_max=gap_max,
                    first_r=r,
                )
                if placement == float("-inf"):
                    continue
                inner = min(0.0, inner_by_gap[idx])
                term = outer_log2 + placement + inner
                total = log2add(total, term)
                by_gap[idx] = log2add(by_gap[idx], term)
                if peak is None or term > peak[4]:
                    peak = (h, r, gap_min, gap_max, term)

    print("Full-split uniform-inner prefix sum")
    print(
        f"N={args.N}, b={args.block_bits}, h_max={args.h_max}, "
        f"late_blocks={args.late_blocks}, r_max={args.r_max}"
    )
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    for (gap_min, gap_max), val in zip(gap_ranges, by_gap):
        print(f"gap_{gap_min}_{gap_max}_log2,{val:.6f}")
    if peak is not None:
        h, r, gap_min, gap_max, term = peak
        print(f"peak_h,{h}")
        print(f"peak_r,{r}")
        print(f"peak_gap_min,{gap_min}")
        print(f"peak_gap_max,{gap_max}")
        print(f"peak_term_log2,{term:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
