#!/usr/bin/env python3
"""Conditioned early-branch probe for the block-recursive BCH inner.

This script forces the first active block to occur before the late-start
window and measures how much later input support cancels the zero-input
trajectory.  It is a proof-development diagnostic: the output identifies the
size of the cancellation lemma needed for the early-first-active branch.
"""

from __future__ import annotations

import argparse
import math
import random

from probe_block_recursive_inner import (
    build_output_systematic_basis,
    build_parity_tables,
    parity_state,
    parse_int_list,
)


def random_mask(rng: random.Random, b: int, weight: int) -> int:
    out = 0
    for bit in rng.sample(range(b), weight):
        out |= 1 << bit
    return out


def sample_later_blocks(
    rng: random.Random,
    *,
    b: int,
    blocks_n: int,
    first_block: int,
    remaining_weight: int,
) -> dict[int, int]:
    later_coords = (blocks_n - first_block - 1) * b
    if remaining_weight > later_coords:
        raise ValueError("not enough later coordinates")
    blocks: dict[int, int] = {}
    for pos in rng.sample(range(later_coords), remaining_weight):
        block = first_block + 1 + pos // b
        off = pos % b
        blocks[block] = blocks.get(block, 0) ^ (1 << off)
    return blocks


def simulate_from_blocks(
    *,
    parity_tables: list[list[int]],
    b: int,
    blocks_n: int,
    input_blocks: dict[int, int],
    first_block: int,
    terminal_state: bool,
    stop_after: int | None = None,
) -> tuple[int, int]:
    state = 0
    total = 0
    end = blocks_n if stop_after is None else min(blocks_n, first_block + stop_after)
    for i in range(first_block, end):
        v = input_blocks.get(i, 0) ^ state
        total += v.bit_count()
        state = parity_state(parity_tables, 16, v) if v else 0
    if terminal_state:
        total += state.bit_count()
    return total, state.bit_count()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--h", type=int, default=32)
    parser.add_argument("--first-r-values", default="1,2")
    parser.add_argument(
        "--gaps-before-late",
        default="1,10,100,500,1000,2000",
        help="Force first block to be this many blocks before the late threshold.",
    )
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    parity_tables = build_parity_tables(sys_basis, b)
    if args.N % b:
        raise ValueError("N must be divisible by block size")
    blocks_n = args.N // b
    d = math.floor(args.distance_delta * args.N)
    rng = random.Random(args.seed)

    print("Block-recursive early-branch cancellation probe")
    print(
        f"m={args.m}, designed_delta={args.delta_bch}, b={b}, N={args.N}, "
        f"blocks={blocks_n}, h={args.h}, d={d}, late_blocks={args.late_blocks}"
    )
    print(
        "first_r,gap,first_block,trials,failures,min_total,avg_total,"
        "min_first_only,avg_first_only,min_later_only,avg_later_only,"
        "max_cancellation,avg_cancellation,min_excess_first"
    )
    late_start = blocks_n - args.late_blocks
    for first_r in parse_int_list(args.first_r_values):
        if first_r > args.h:
            continue
        for gap in parse_int_list(args.gaps_before_late):
            first_block = late_start - gap
            if first_block < 0:
                continue
            failures = 0
            min_total = 10**18
            sum_total = 0
            min_first = 10**18
            sum_first = 0
            min_later = 10**18
            sum_later = 0
            max_cancel = 0
            sum_cancel = 0
            min_excess_first = 10**18
            for _ in range(args.trials):
                first_mask = random_mask(rng, b, first_r)
                later = sample_later_blocks(
                    rng,
                    b=b,
                    blocks_n=blocks_n,
                    first_block=first_block,
                    remaining_weight=args.h - first_r,
                )
                combined = dict(later)
                combined[first_block] = first_mask
                first_only = {first_block: first_mask}

                total, _ = simulate_from_blocks(
                    parity_tables=parity_tables,
                    b=b,
                    blocks_n=blocks_n,
                    input_blocks=combined,
                    first_block=first_block,
                    terminal_state=True,
                )
                first_total, _ = simulate_from_blocks(
                    parity_tables=parity_tables,
                    b=b,
                    blocks_n=blocks_n,
                    input_blocks=first_only,
                    first_block=first_block,
                    terminal_state=True,
                )
                later_first = min(later) if later else blocks_n
                later_total, _ = simulate_from_blocks(
                    parity_tables=parity_tables,
                    b=b,
                    blocks_n=blocks_n,
                    input_blocks=later,
                    first_block=later_first,
                    terminal_state=True,
                ) if later else (0, 0)
                cancellation = first_total + later_total - total

                failures += int(total <= d)
                min_total = min(min_total, total)
                sum_total += total
                min_first = min(min_first, first_total)
                sum_first += first_total
                min_later = min(min_later, later_total)
                sum_later += later_total
                max_cancel = max(max_cancel, cancellation)
                sum_cancel += cancellation
                min_excess_first = min(min_excess_first, first_total - d)
            print(
                f"{first_r},{gap},{first_block},{args.trials},{failures},"
                f"{min_total},{sum_total / args.trials:.6f},"
                f"{min_first},{sum_first / args.trials:.6f},"
                f"{min_later},{sum_later / args.trials:.6f},"
                f"{max_cancel},{sum_cancel / args.trials:.6f},{min_excess_first}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
