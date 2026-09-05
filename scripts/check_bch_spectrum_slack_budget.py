#!/usr/bin/env python3
"""Estimate the finite first-moment slack budget for the native BCH block.

The random-like [511,256,61]-shaped diagnostic has substantial headroom, but
the proof needs a finite upper spectrum envelope.  This script asks how many
uniform bits of multiplicative loss in the nonzero local spectrum can be
tolerated before the block-outer first moment loses a requested margin.

This remains a diagnostic, because the input spectrum is random-like.  Its
output is useful as a target for a future certified BCH parent/subcode
spectrum envelope.
"""

from __future__ import annotations

import argparse
import math

from block_outer_upgrade_probe import format_log2, log2_sub_one, z_grid
from certify_global_episode_cover import (
    binom_cdf_half_be_log2,
    global_episode_inner_log2,
    survivor_only_exact_prefix_log2,
)
from dense_largek_eval import log2_binom, log2add


def parse_float_list(text: str) -> list[float]:
    out: list[float] = []
    for part in text.split(","):
        part = part.strip()
        if part:
            out.append(float(part))
    return out


def local_random_nonzero_log2(block_bits: int, local_length: int, d0: int, z: float) -> float:
    """Random-linear expected nonzero local enumerator above a hard floor."""
    log_rate = block_bits - local_length
    total = float("-inf")
    log_z = math.log2(z)
    for j in range(d0, local_length + 1):
        total = log2add(total, log_rate + log2_binom(local_length, j) + j * log_z)
    return total


def local_slack_total_log2(
    *,
    block_bits: int,
    local_length: int,
    d0: int,
    z: float,
    slack_bits: float,
) -> float:
    nonzero = local_random_nonzero_log2(block_bits, local_length, d0, z) + slack_bits
    if nonzero < -60.0:
        return math.log1p(2.0**nonzero) / math.log(2.0)
    return math.log2(1.0 + 2.0**nonzero)


def outer_with_slack(
    *,
    blocks: int,
    block_bits: int,
    local_length: int,
    d0: int,
    h_max: int,
    zs: list[float],
    slack_bits: float,
) -> tuple[list[float], list[float]]:
    vals = [float("-inf")] * (h_max + 1)
    best_z = [float("nan")] * (h_max + 1)
    local_logs = []
    for z in zs:
        local = local_slack_total_log2(
            block_bits=block_bits,
            local_length=local_length,
            d0=d0,
            z=z,
            slack_bits=slack_bits,
        )
        local_logs.append((z, log2_sub_one(blocks * local)))

    for h in range(d0, h_max + 1):
        best = float("inf")
        best_here = float("nan")
        for z, global_log in local_logs:
            if global_log == float("-inf"):
                continue
            candidate = global_log - h * math.log2(z)
            if candidate < best:
                best = candidate
                best_here = z
        vals[h] = best if best != float("inf") else float("-inf")
        best_z[h] = best_here
    return vals, best_z


def total_for_slack(
    *,
    blocks: int,
    block_bits: int,
    local_length: int,
    d0: int,
    h_max: int,
    zs: list[float],
    inner_logs: list[float],
    slack_bits: float,
) -> tuple[float, int, float, float, float]:
    outer_logs, _ = outer_with_slack(
        blocks=blocks,
        block_bits=block_bits,
        local_length=local_length,
        d0=d0,
        h_max=h_max,
        zs=zs,
        slack_bits=slack_bits,
    )
    total = float("-inf")
    peak_h = -1
    peak_term = float("-inf")
    peak_outer = float("-inf")
    peak_inner = float("-inf")
    for h in range(1, h_max + 1):
        outer = outer_logs[h]
        inner = inner_logs[h]
        term = outer + inner if outer != float("-inf") and inner != float("-inf") else float("-inf")
        total = log2add(total, term)
        if term > peak_term:
            peak_h = h
            peak_term = term
            peak_outer = outer
            peak_inner = inner
    return total, peak_h, peak_term, peak_outer, peak_inner


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=2**20)
    parser.add_argument("--block-bits", type=int, default=256)
    parser.add_argument("--local-length", type=int, default=511)
    parser.add_argument("--d0", type=int, default=61)
    parser.add_argument("--sigma", type=int, default=32)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--h-max", type=int, default=220)
    parser.add_argument("--r-max", type=int, default=12)
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    parser.add_argument("--block-ratio", type=float, default=1.02)
    parser.add_argument("--suffix-block-ratio", type=float, default=1.1)
    parser.add_argument("--slacks", default="0,3,10,20,25,28,29,30")
    parser.add_argument("--target-margin", type=float, default=40.0)
    args = parser.parse_args()

    blocks = (args.k + args.block_bits - 1) // args.block_bits
    n = blocks * args.local_length
    d = math.floor(args.delta * n)
    zs = z_grid(args.z_min, args.z_max, args.z_count)
    slacks = parse_float_list(args.slacks)

    print("BCH random-like spectrum slack budget")
    print(f"k={args.k}, blocks={blocks}, local=[{args.local_length},{args.block_bits},>={args.d0}], N={n}")
    print(f"sigma={args.sigma}, delta={args.delta}, d={d}, h_max={args.h_max}")
    print("precomputing inner terms", flush=True)
    tail_cache = [binom_cdf_half_be_log2(u, d) for u in range(n + 1)]
    exact_survivor = survivor_only_exact_prefix_log2(n=n, h_max=args.h_max, tail_cache=tail_cache)
    inner_logs = [float("-inf")] * (args.h_max + 1)
    for h in range(1, args.h_max + 1):
        inner, _, _, _ = global_episode_inner_log2(
            n=n,
            h=h,
            sigma=args.sigma,
            delta=args.delta,
            r_max=args.r_max,
            block_ratio=args.block_ratio,
            suffix_block_ratio=args.suffix_block_ratio,
            stop_gap_bits=60.0,
            tail_confirm=8,
            suffix_survivor=True,
            tail_mode="be",
            exact_survivor=True,
            tail_cache=tail_cache,
            exact_survivor_value=exact_survivor[h],
        )
        inner_logs[h] = inner

    print("slack_bits,total_log2,margin_bits,peak_h,peak_term_log2,peak_outer_log2,peak_inner_log2")
    for slack in slacks:
        total, peak_h, peak_term, peak_outer, peak_inner = total_for_slack(
            blocks=blocks,
            block_bits=args.block_bits,
            local_length=args.local_length,
            d0=args.d0,
            h_max=args.h_max,
            zs=zs,
            inner_logs=inner_logs,
            slack_bits=slack,
        )
        print(
            f"{slack:.3f},{format_log2(total)},{format_log2(-total)},{peak_h},"
            f"{format_log2(peak_term)},{format_log2(peak_outer)},{format_log2(peak_inner)}"
        )

    lo = 0.0
    hi = max(max(slacks) if slacks else 0.0, 1.0)
    while True:
        total, *_ = total_for_slack(
            blocks=blocks,
            block_bits=args.block_bits,
            local_length=args.local_length,
            d0=args.d0,
            h_max=args.h_max,
            zs=zs,
            inner_logs=inner_logs,
            slack_bits=hi,
        )
        if -total < args.target_margin:
            break
        lo = hi
        hi *= 2.0
    for _ in range(36):
        mid = 0.5 * (lo + hi)
        total, *_ = total_for_slack(
            blocks=blocks,
            block_bits=args.block_bits,
            local_length=args.local_length,
            d0=args.d0,
            h_max=args.h_max,
            zs=zs,
            inner_logs=inner_logs,
            slack_bits=mid,
        )
        if -total >= args.target_margin:
            lo = mid
        else:
            hi = mid
    total, peak_h, peak_term, peak_outer, peak_inner = total_for_slack(
        blocks=blocks,
        block_bits=args.block_bits,
        local_length=args.local_length,
        d0=args.d0,
        h_max=args.h_max,
        zs=zs,
        inner_logs=inner_logs,
        slack_bits=lo,
    )
    print()
    print(
        f"slack budget for {args.target_margin:.3f} bits: {lo:.6f} bits "
        f"(total={format_log2(total)}, peak_h={peak_h}, peak={format_log2(peak_term)}, "
        f"outer={format_log2(peak_outer)}, inner={format_log2(peak_inner)})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
