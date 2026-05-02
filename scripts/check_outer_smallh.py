#!/usr/bin/env python3
"""Compare exact systematic-outer low-weight counts against the current formula."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from check_inner_smallw import parse_enum
from dense_largek_eval import outer_small_h_exact_log2


def exact_outer_weight(rows: list[list[float]], h: int) -> float:
    total = 0.0
    for row in rows:
        if h < len(row):
            v = row[h]
            if v != float("-inf"):
                total += 2.0**v
    return total


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enum", type=Path, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--sigma", type=int, required=True)
    ap.add_argument("--h-max", type=int, default=12)
    args = ap.parse_args()

    rows = parse_enum(args.enum)
    print(f"enum   : {args.enum}")
    print(f"k      : {args.k}")
    print(f"sigma  : {args.sigma}")
    print()
    print("h, exact, model, model/exact, log2(model/exact)")

    for h in range(1, args.h_max + 1):
        exact = exact_outer_weight(rows, h)
        log2_model = outer_small_h_exact_log2(args.k, args.sigma, h)
        model = 0.0 if log2_model == float("-inf") else 2.0**log2_model
        ratio = model / exact if exact > 0.0 else float("inf")
        bits = math.log2(ratio) if ratio > 0.0 else float("-inf")
        print(f"{h}, {exact:.6e}, {model:.6e}, {ratio:.6e}, {bits:.6f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
