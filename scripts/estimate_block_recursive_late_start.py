#!/usr/bin/env python3
"""Late-start balance estimates for the block-recursive BCH inner.

This is a deliberately simple first-moment sanity check.  If a live
block-recursive episode emits about avg_output ones per zero-input block, then
an h-support can have low output only if its first active block lies in roughly
the final floor(d/avg_output) blocks.  The corresponding placement factor is

    C(L*b, h) / C(N, h).

The estimate ignores cancellation and local output fluctuations; it is not a
certificate.  It is useful for asking whether a proposed outer floor has enough
room to reach a target margin.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from dense_largek_eval import log2_binom


def load_outer_logs(summary_csv: Path | None, prefix_csv: Path | None) -> dict[int, float]:
    logs: dict[int, float] = {}
    if prefix_csv is not None:
        with prefix_csv.open(newline="") as f:
            for row in csv.DictReader(f):
                if "h" not in row or "outer_log2" not in row:
                    continue
                logs[int(row["h"])] = float(row["outer_log2"])
    return logs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--avg-output", type=float, default=32.0)
    parser.add_argument("--h-values", default="4,8,16,22,24,25,32,40,48,64")
    parser.add_argument("--outer-prefix-csv", type=Path, default=None)
    parser.add_argument("--target-bits", type=float, default=40.0)
    args = parser.parse_args()

    if args.N % args.block_bits != 0:
        raise ValueError("N must be divisible by block size")
    d = math.floor(args.delta * args.N)
    live_blocks = min(args.N // args.block_bits, math.floor(d / args.avg_output))
    live_coords = live_blocks * args.block_bits
    outer_logs = load_outer_logs(None, args.outer_prefix_csv)

    print("Block-recursive late-start estimate")
    print(
        f"N={args.N}, b={args.block_bits}, delta={args.delta}, d={d}, "
        f"avg_output={args.avg_output}, live_blocks={live_blocks}, live_coords={live_coords}"
    )
    print("h,late_log2,outer_log2,term_log2,margin_bits,has_outer")
    for part in args.h_values.split(","):
        part = part.strip()
        if not part:
            continue
        h = int(part)
        late = log2_binom(live_coords, h) - log2_binom(args.N, h)
        if h in outer_logs:
            outer = outer_logs[h]
            term = outer + late
            margin = -term
            has_outer = 1
        else:
            outer = float("nan")
            term = float("nan")
            margin = float("nan")
            has_outer = 0
        def fmt(x: float) -> str:
            if math.isnan(x):
                return "nan"
            if x == float("inf"):
                return "inf"
            if x == float("-inf"):
                return "-inf"
            return f"{x:.6f}"

        print(f"{h},{fmt(late)},{fmt(outer)},{fmt(term)},{fmt(margin)},{has_outer}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
