#!/usr/bin/env python3
"""Inspect the deterministic cluster count forced by gaps shorter than sigma."""

from __future__ import annotations

import argparse
import math


def forced_cluster_mean(n: int, w: int, s: int, m: int) -> float:
    if s <= 1:
        return 1.0 if w > 0 else 0.0
    M = n - w - (s - 1)
    if M < 0:
        return float("nan")
    num = math.comb(M - m + s + 1, s) if M - m + 1 >= 0 else 0
    den = math.comb(M + s, s)
    p_long = num / den
    return 1.0 + (s - 1) * p_long


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--w-max", type=int, default=10)
    args = ap.parse_args()

    print(f"n = {args.n}")
    print(f"m = {args.m}")
    print("w, s, E[C_m | w,s]")
    for w in range(1, args.w_max + 1):
        for s in range(1, w + 1):
            print(f"{w}, {s}, {forced_cluster_mean(args.n, w, s, args.m):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
