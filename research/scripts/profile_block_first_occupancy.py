#!/usr/bin/env python3
"""Exact first-active-block occupancy probabilities for uniform weight-h inputs."""

from __future__ import annotations

import argparse
import math

from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--h-values", default="4,8,16,22,24,25,32,40,64")
    parser.add_argument("--late-blocks", type=int, default=None)
    args = parser.parse_args()

    if args.N % args.block_bits != 0:
        raise ValueError("N must be divisible by block size")
    blocks = args.N // args.block_bits
    b = args.block_bits
    late_blocks = args.late_blocks if args.late_blocks is not None else blocks

    print("First active block occupancy")
    print(f"N={args.N}, block_bits={b}, blocks={blocks}, late_blocks={late_blocks}")
    print("h,r,log2_prob,prob,log2_prob_and_late,prob_and_late")
    for h in parse_int_list(args.h_values):
        denom = log2_binom(args.N, h)
        for r in range(1, min(b, h) + 1):
            total = float("-inf")
            late_total = float("-inf")
            for t in range(blocks):
                later = (blocks - t - 1) * b
                if h - r > later:
                    continue
                term = log2_binom(b, r) + log2_binom(later, h - r) - denom
                total = log2add(total, term)
                if t >= blocks - late_blocks:
                    late_total = log2add(late_total, term)
            if total == float("-inf"):
                continue
            print(
                f"{h},{r},{total:.6f},{2**total:.12g},"
                f"{late_total:.6f},{0.0 if late_total == float('-inf') else 2**late_total:.12g}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
