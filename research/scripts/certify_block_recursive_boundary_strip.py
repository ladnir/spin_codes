#!/usr/bin/env python3
"""First-moment mass of a block-recursive boundary strip.

The late/buffer branch covers inputs whose first active block is in the final
L_late blocks.  This script measures the adjacent strip just before that
window:

    first active block is g blocks before the late window, g_min <= g <= g_max.

No inner-distance gain is claimed here.  The output is a bookkeeping
certificate for the proof split: it tells us how much first-moment mass the
future cancellation/column-distance lemma must control.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from dense_largek_eval import log2_binom, log2add


def load_outer(path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if "h" in row and "outer_log2" in row:
                val = row["outer_log2"]
                if val and val != "-inf":
                    out[int(row["h"])] = float(val)
    return out


def first_block_gap_log2(*, n: int, b: int, h: int, late_blocks: int, gap: int, first_r: int) -> float:
    """Probability of a fixed first-active-block gap and first occupancy."""

    after_coords = b * (late_blocks + gap - 1)
    if first_r < 1 or first_r > min(b, h):
        return float("-inf")
    if h - first_r > after_coords:
        return float("-inf")
    return log2_binom(b, first_r) + log2_binom(after_coords, h - first_r) - log2_binom(n, h)


def log2diff(a: float, b: float) -> float:
    if b == float("-inf"):
        return a
    if b > a:
        raise ValueError("log2diff requires a >= b")
    if a == b:
        return float("-inf")
    return a + math.log2(1.0 - 2.0 ** (b - a))


def first_active_strip_log2(*, n: int, b: int, h: int, late_blocks: int, gap_count: int) -> float:
    enlarged = log2_binom((late_blocks + gap_count) * b, h)
    late = log2_binom(late_blocks * b, h)
    return log2diff(enlarged, late) - log2_binom(n, h)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-min", type=int, default=1)
    parser.add_argument("--gap-max", type=int, default=51)
    parser.add_argument("--outer-prefix-csv", type=Path, required=True)
    parser.add_argument("--h-max", type=int, default=500)
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()

    if args.N % args.block_bits:
        raise ValueError("N must be divisible by block size")
    outer = load_outer(args.outer_prefix_csv)
    total = float("-inf")
    peak: tuple[int, int, int, float, float, float] | None = None

    print("Block-recursive boundary-strip mass")
    print(
        f"N={args.N}, block_bits={args.block_bits}, late_blocks={args.late_blocks}, "
        f"gap_min={args.gap_min}, gap_max={args.gap_max}, h_max={args.h_max}"
    )
    if not args.summary_only:
        print("h,gap,first_r,outer_log2,strip_log2,term_log2")
    if args.summary_only:
        gap_count = args.gap_max - args.gap_min + 1
        if args.gap_min != 1:
            raise ValueError("summary-only fast path currently requires --gap-min 1")
        for h in range(1, args.h_max + 1):
            if h not in outer:
                continue
            strip = first_active_strip_log2(
                n=args.N,
                b=args.block_bits,
                h=h,
                late_blocks=args.late_blocks,
                gap_count=gap_count,
            )
            term = outer[h] + strip
            total = log2add(total, term)
            if peak is None or term > peak[5]:
                # Peak gap/first_r are not resolved on the fast aggregate path.
                peak = (h, args.gap_max, -1, outer[h], strip, term)
    else:
        for h in range(1, args.h_max + 1):
            if h not in outer:
                continue
            max_r = min(args.block_bits, h)
            for gap in range(args.gap_min, args.gap_max + 1):
                for first_r in range(1, max_r + 1):
                    strip = first_block_gap_log2(
                        n=args.N,
                        b=args.block_bits,
                        h=h,
                        late_blocks=args.late_blocks,
                        gap=gap,
                        first_r=first_r,
                    )
                    if strip == float("-inf"):
                        continue
                    term = outer[h] + strip
                    total = log2add(total, term)
                    if peak is None or term > peak[5]:
                        peak = (h, gap, first_r, outer[h], strip, term)
                    print(f"{h},{gap},{first_r},{outer[h]:.6f},{strip:.6f},{term:.6f}")

    print("summary")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    if peak is not None:
        h, gap, first_r, outer_log, strip_log, term_log = peak
        print(f"peak_h,{h}")
        print(f"peak_gap,{gap}")
        print(f"peak_first_r,{first_r}")
        print(f"peak_outer_log2,{outer_log:.6f}")
        print(f"peak_strip_log2,{strip_log:.6f}")
        print(f"peak_term_log2,{term_log:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
