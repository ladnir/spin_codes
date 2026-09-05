#!/usr/bin/env python3
"""Small utilities for the gap-survival law behind episode clustering."""

from __future__ import annotations

import argparse
import math


def no_m_zero_run_prob(g: int, m: int) -> float:
    """Probability a fair-coin string of length g has no run of m zeros."""
    if g < 0:
        return 0.0
    if m <= 0:
        return 0.0
    if g < m:
        return 1.0

    # dp[r] = probability mass of prefixes whose trailing zero-run has length r,
    # for r=0..m-1, while never having hit an m-zero run.
    dp = [0.0] * m
    dp[0] = 1.0
    for _ in range(g):
        nxt = [0.0] * m
        # append a 1: reset trailing zero count
        nxt[0] += 0.5 * sum(dp)
        # append a 0: increase trailing zero count if still < m
        for r in range(m - 1):
            nxt[r + 1] += 0.5 * dp[r]
        dp = nxt
    return sum(dp)


def mean_internal_gap_survival(n: int, w: int, m: int) -> float:
    """Average gap-survival probability for the s=2 run model at weight w.

    This averages over the exact internal gap distribution when a weight-w word
    has exactly two runs of ones. The gap length g ranges from 1 up to n-w,
    and the number of zero-gap placements with a fixed internal gap g is
    (n-w-g+1), corresponding to the split of the remaining boundary zeros.
    """
    z = n - w
    if z < 1:
        return 0.0
    total = 0.0
    denom = 0.0
    for g in range(1, z + 1):
        mult = z - g + 1
        total += mult * no_m_zero_run_prob(g, m)
        denom += mult
    return total / denom if denom > 0 else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--g-max", type=int, default=24)
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--w", type=int, default=None)
    args = ap.parse_args()

    print(f"gap survival for m={args.m}")
    for g in range(args.g_max + 1):
        exact = no_m_zero_run_prob(g, args.m)
        ub = max(0.0, 1.0 - g * (2.0 ** (-args.m)))
        print(f"  g={g:2d}: exact={exact:.6f}  union-lb={ub:.6f}")

    if args.n is not None and args.w is not None:
        avg = mean_internal_gap_survival(args.n, args.w, args.m)
        print()
        print(f"mean internal-gap survival for n={args.n}, w={args.w}, s=2: {avg:.6f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
