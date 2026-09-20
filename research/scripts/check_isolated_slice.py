#!/usr/bin/env python3
"""Compare exact small-weight inner slices against the isolated late-placement base."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from check_q1_law import parse_enum
from dense_largek_eval import late_start_base_prob


def exact_slice_prob(rows: list[list[float]], n: int, w: int, cut: int) -> float:
    row = rows[w]
    num = sum((0.0 if v == float("-inf") else 2.0**v) for v in row[: cut + 1])
    return num / math.comb(n, w)


def isolated_full_run_prob(n: int, h: int) -> float:
    return math.comb(n - h + 1, h) / math.comb(n, h)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enum", type=Path, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--delta", type=float, required=True)
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
    print()
    print("h, exact, isolated_base, gap_bits, gap_bits_per_pair, Pr[R=h]")
    for h in range(1, min(args.h_max, len(rows) - 1) + 1):
        exact = exact_slice_prob(rows, args.n, h, cut)
        p_full = isolated_full_run_prob(args.n, h)
        base = p_full * late_start_base_prob(args.n, h, h, l)
        if h == 1:
            gap_bits = 0.0 if exact == 0.0 else math.log2(exact / base)
            gap_bits_per_pair = 0.0
        else:
            gap_bits = math.log2(exact / base)
            gap_bits_per_pair = gap_bits / (h * (h - 1) / 2.0)
        print(
            f"{h}, {exact:.6e}, {base:.6e}, {gap_bits:.6f}, {gap_bits_per_pair:.6f}, {p_full:.6e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
