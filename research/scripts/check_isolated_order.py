#!/usr/bin/env python3
"""Check exact slice probabilities against isolated order-statistic events."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from check_q1_law import parse_enum
from dense_largek_eval import isolated_full_run_prob, isolated_order_late_prob


def exact_slice_prob(rows: list[list[float]], n: int, w: int, cut: int) -> float:
    row = rows[w]
    num = sum((0.0 if v == float("-inf") else 2.0**v) for v in row[: cut + 1])
    return num / math.comb(n, w)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enum", type=Path, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--delta", type=float, required=True)
    ap.add_argument("--r", type=int, default=2)
    ap.add_argument("--h-max", type=int, default=10)
    args = ap.parse_args()

    rows = parse_enum(args.enum)
    cut = math.floor(args.delta * args.n)
    l = math.ceil(2.0 * args.delta * args.n)

    print(f"enum   : {args.enum}")
    print(f"n      : {args.n}")
    print(f"delta  : {args.delta}")
    print(f"cut    : {cut}")
    print(f"L      : {l}")
    print(f"r      : {args.r}")
    print()
    print("h, exact, order_bound, gap_bits")
    for h in range(max(1, args.r), min(args.h_max, len(rows) - 1) + 1):
        exact = exact_slice_prob(rows, args.n, h, cut)
        bound = isolated_full_run_prob(args.n, h) * isolated_order_late_prob(args.n, h, args.r, l)
        bound += 1.0 - isolated_full_run_prob(args.n, h)
        gap_bits = math.log2(bound / exact)
        print(f"{h}, {exact:.6e}, {bound:.6e}, {gap_bits:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
