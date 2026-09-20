#!/usr/bin/env python3
"""First-moment contribution of the block-recursive late/buffer branch.

For a block-recursive inner with M=N/b blocks, the late/buffer branch is the
event that the first active input block lies in the final L blocks.  This event
has an exact combinatorial probability for a uniform weight-h input.  This
script multiplies that probability by a block-outer spectrum prefix.

This certifies only the late/buffer branch.  The early branch still needs the
column-distance/cancellation ledger.
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


def first_active_late_log2(*, n: int, b: int, h: int, late_blocks: int) -> float:
    return log2_binom(late_blocks * b, h) - log2_binom(n, h)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--late-blocks", type=int, default=6050)
    parser.add_argument("--outer-prefix-csv", type=Path, required=True)
    parser.add_argument("--h-max", type=int, default=500)
    args = parser.parse_args()

    if args.N % args.block_bits:
        raise ValueError("N must be divisible by block size")
    outer = load_outer(args.outer_prefix_csv)
    total = float("-inf")
    peak: tuple[int, float, float, float] | None = None

    print("Block-recursive late/buffer branch certificate")
    print(f"N={args.N}, block_bits={args.block_bits}, late_blocks={args.late_blocks}, h_max={args.h_max}")
    print("h,outer_log2,late_log2,term_log2")
    for h in range(1, args.h_max + 1):
        if h not in outer:
            continue
        late = first_active_late_log2(n=args.N, b=args.block_bits, h=h, late_blocks=args.late_blocks)
        term = outer[h] + late
        total = log2add(total, term)
        if peak is None or term > peak[3]:
            peak = (h, outer[h], late, term)
        print(f"{h},{outer[h]:.6f},{late:.6f},{term:.6f}")

    print("summary")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    if peak is not None:
        h, outer_log, late_log, term_log = peak
        print(f"peak_h,{h}")
        print(f"peak_outer_log2,{outer_log:.6f}")
        print(f"peak_late_log2,{late_log:.6f}")
        print(f"peak_term_log2,{term_log:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
