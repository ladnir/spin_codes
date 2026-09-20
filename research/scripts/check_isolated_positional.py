#!/usr/bin/env python3
"""Check the position-sensitive isolated one-run model against exact inner-only tables."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from check_q1_law import parse_enum
from dense_largek_eval import isolated_position_single_model_prob


def exact_slice_prob(rows: list[list[float]], n: int, w: int, cut: int) -> float:
    row = rows[w]
    num = sum((0.0 if v == float("-inf") else 2.0**v) for v in row[: cut + 1])
    return num / math.comb(n, w)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enum", type=Path, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--delta", type=float, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--h-max", type=int, default=10)
    args = ap.parse_args()

    rows = parse_enum(args.enum)
    cut = math.floor(args.delta * args.n)

    print(f"enum   : {args.enum}")
    print(f"n      : {args.n}")
    print(f"delta  : {args.delta}")
    print(f"cut    : {cut}")
    print(f"m      : {args.m}")
    print()
    print("h, exact, positional_model, gap_bits")
    for h in range(1, min(args.h_max, len(rows) - 1) + 1):
        exact = exact_slice_prob(rows, args.n, h, cut)
        model = isolated_position_single_model_prob(args.n, h, cut, args.m)
        gap_bits = math.log2(model / exact)
        print(f"{h}, {exact:.6e}, {model:.6e}, {gap_bits:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
