#!/usr/bin/env python3
"""Exact-turnoff ledger for the block-recursive BCH inner.

After the first active block, the cheapest way to suppress the zero-input
trajectory is for a later input block to equal the current state exactly,
turning the recursion off.  This script bounds that event by a union bound over
candidate turnoff blocks.

The calculation is conditional on:

* total inner input weight h;
* first active block is g blocks before the late threshold;
* first active block occupancy is r.

For each first-block mask of weight r, we follow the zero-input trajectory.  A
turnoff at later offset t is only counted as dangerous if the accumulated output
before that turnoff is still at most d.  Given remaining weight H=h-r uniformly
placed in later coordinates, the probability that block t equals a prescribed
state of weight q and has no extra ones is

    C(after_t, H-q) / C(total_later, H).

This is not the final early-branch proof: near-turnoffs and later restarts still
need handling.  It isolates the cleanest structural obstruction.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from check_bch_generator_inner import iter_weight_masks
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import (
    build_output_systematic_basis,
    build_parity_tables,
    parity_state,
    parse_int_list,
)


def log2_subprob_exact_turnoff(
    *,
    total_later_coords: int,
    after_coords: int,
    remaining_weight: int,
    state_weight: int,
) -> float:
    if state_weight > remaining_weight:
        return float("-inf")
    if remaining_weight - state_weight > after_coords:
        return float("-inf")
    return log2_binom(after_coords, remaining_weight - state_weight) - log2_binom(
        total_later_coords, remaining_weight
    )


def first_block_gap_log2(*, n: int, b: int, h: int, late_blocks: int, gap: int, first_r: int) -> float:
    after_coords = b * (late_blocks + gap - 1)
    if first_r < 1 or first_r > min(b, h):
        return float("-inf")
    if h - first_r > after_coords:
        return float("-inf")
    return log2_binom(b, first_r) + log2_binom(after_coords, h - first_r) - log2_binom(n, h)


def dangerous_turnoff_log2(
    *,
    tables: list[list[int]],
    b: int,
    first_mask: int,
    live_blocks: int,
    remaining_weight: int,
    d: int,
) -> tuple[float, int, int, int, int, int, float]:
    """Return log2 union bound and a few diagnostics for one first mask."""

    total_later = (live_blocks - 1) * b
    state = first_mask
    prefix = 0
    total = float("-inf")
    counted = 0
    last_t = 0
    min_state_q = b + 1
    peak_t = 0
    peak_prefix = 0
    peak_q = 0
    peak_log = float("-inf")

    # t=0 is the first active block output.  A later exact turnoff at offset t
    # emits zero at block t, so its weight is sum outputs at offsets < t.
    for t in range(0, live_blocks):
        if t > 0:
            q = state.bit_count()
            if prefix <= d and q <= remaining_weight:
                after = (live_blocks - t - 1) * b
                term = log2_subprob_exact_turnoff(
                    total_later_coords=total_later,
                    after_coords=after,
                    remaining_weight=remaining_weight,
                    state_weight=q,
                )
                total = log2add(total, term)
                counted += int(term != float("-inf"))
                last_t = t
                min_state_q = min(min_state_q, q)
                if term > peak_log:
                    peak_log = term
                    peak_t = t
                    peak_prefix = prefix
                    peak_q = q
        prefix += state.bit_count()
        if prefix > d and t > 0:
            break
        state = parity_state(tables, 16, state) if state else 0
    if min_state_q == b + 1:
        min_state_q = -1
    return total, counted, last_t, min_state_q, peak_t, peak_prefix, peak_q, peak_log


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--h-values", default="32")
    parser.add_argument("--first-r-values", default="1,2")
    parser.add_argument("--gap-values", default="102,200,500,1000")
    parser.add_argument("--exact-mask-limit", type=int, default=5000)
    parser.add_argument("--outer-log2", type=float, default=None)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    tables = build_parity_tables(sys_basis, b)
    if args.N % b:
        raise ValueError("N must be divisible by block size")
    blocks_n = args.N // b
    d = math.floor(args.distance_delta * args.N)

    print("Block-recursive exact-turnoff ledger")
    print(
        f"m={args.m}, designed_delta={args.delta_bch}, b={b}, local_length={ext_n}, "
        f"N={args.N}, blocks={blocks_n}, d={d}, late_blocks={args.late_blocks}, row_d0={row_d0}"
    )
    print(
        "h,first_r,gap,live_blocks,masks,first_log2,worst_cond_turnoff_log2,"
        "term_log2,with_outer_log2,worst_counted,worst_last_t,worst_min_state_q,"
        "peak_t,peak_prefix,peak_q,peak_cond_log2"
    )

    late_start = blocks_n - args.late_blocks
    for h in parse_int_list(args.h_values):
        for first_r in parse_int_list(args.first_r_values):
            if first_r > h:
                continue
            mask_count = math.comb(b, first_r)
            if mask_count > args.exact_mask_limit:
                print(f"# skipping h={h}, first_r={first_r}: {mask_count} masks exceed exact limit")
                continue
            masks = list(iter_weight_masks(b, first_r))
            for gap in parse_int_list(args.gap_values):
                first_block = late_start - gap
                if first_block < 0:
                    continue
                live_blocks = blocks_n - first_block
                first_log = first_block_gap_log2(
                    n=args.N,
                    b=b,
                    h=h,
                    late_blocks=args.late_blocks,
                    gap=gap,
                    first_r=first_r,
                )
                worst = float("-inf")
                worst_counted = 0
                worst_last_t = 0
                worst_min_q = -1
                worst_peak_t = 0
                worst_peak_prefix = 0
                worst_peak_q = 0
                worst_peak_log = float("-inf")
                for mask in masks:
                    cond, counted, last_t, min_q, peak_t, peak_prefix, peak_q, peak_log = dangerous_turnoff_log2(
                        tables=tables,
                        b=b,
                        first_mask=mask,
                        live_blocks=live_blocks,
                        remaining_weight=h - first_r,
                        d=d,
                    )
                    if cond > worst:
                        worst = cond
                        worst_counted = counted
                        worst_last_t = last_t
                        worst_min_q = min_q
                        worst_peak_t = peak_t
                        worst_peak_prefix = peak_prefix
                        worst_peak_q = peak_q
                        worst_peak_log = peak_log
                term = first_log + worst
                with_outer = term + args.outer_log2 if args.outer_log2 is not None else float("nan")
                print(
                    f"{h},{first_r},{gap},{live_blocks},{len(masks)},"
                    f"{first_log:.6f},{worst:.6f},{term:.6f},{with_outer:.6f},"
                    f"{worst_counted},{worst_last_t},{worst_min_q},"
                    f"{worst_peak_t},{worst_peak_prefix},{worst_peak_q},{worst_peak_log:.6f}"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
