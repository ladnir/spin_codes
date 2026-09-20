#!/usr/bin/env python3
"""Scan positive-coefficient blocks for a binary biregular EC tail."""

from __future__ import annotations

import argparse
import math

from scipy.special import logsumexp

from binary_biregular_diagnostic import (
    geometric_blocks,
    saddle_biregular_block_logterm,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=1_048_575)
    parser.add_argument("--left-degree", type=int, default=6)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--cutoff", type=int, default=230_729)
    parser.add_argument("--memory", type=int, default=80)
    parser.add_argument("--support-start", type=int, default=384)
    parser.add_argument("--support-limit", type=int, default=1_048_574)
    parser.add_argument("--relative-width", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    blocks = geometric_blocks(
        args.support_start, args.support_limit, args.relative_width
    )
    terms = []
    for index, (lo, hi) in enumerate(blocks, start=1):
        term = saddle_biregular_block_logterm(
            code="ec",
            k=args.k,
            left_degree=args.left_degree,
            right_degree=args.right_degree,
            cutoff=args.cutoff,
            memory=args.memory,
            support_start=lo,
            support_limit=hi,
        )
        terms.append(term)
        print(
            f"{index}/{len(blocks)}\t{lo}..{hi}\t"
            f"x={term.input_marker:.12g}\tz={term.output_marker:.12g}\t"
            f"log2={term.log2_bound:.9f}",
            flush=True,
        )
    logs = [term.log2_bound * math.log(2.0) for term in terms]
    worst = max(terms, key=lambda term: term.log2_bound)
    print(f"block_count={len(terms)}")
    print(f"worst_block={worst.support_start}..{worst.support_limit}")
    print(f"worst_log2_bound={worst.log2_bound:.12f}")
    print(f"summed_log2_bound={float(logsumexp(logs) / math.log(2.0)):.12f}")


if __name__ == "__main__":
    main()
