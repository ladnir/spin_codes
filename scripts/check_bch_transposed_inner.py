#!/usr/bin/env python3
"""Low-weight diagnostics for a transposed-BCH block-recursive inner.

For a primitive narrow-sense BCH code of length n=2^m-1, extend by one parity
coordinate.  The extended [n+1,k] code has a parity-check matrix H with
redundancy b=n+1-k.  When k=b this is exactly the proposed block recursion

    s_t = y_t = H (s_{t-1} | u_t),

with b-bit state/input/output blocks.

This script enumerates low-weight x=(s,u) and measures the syndrome/output
weight wt(Hx).  The zero-syndrome count is the local termination barrier; the
small nonzero syndrome weights are the possible "live but low-output" leakage.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random

from bch_candidate_params import bch_dimension
from dense_largek_eval import log2_binom
from exact_bch_spectrum_small import parity_rows
from check_bch_boundary_smallfield import find_primitive_poly


def row_rank(rows: list[int]) -> int:
    rows = [r for r in rows if r]
    rank = 0
    for col in range(max(rows).bit_length() if rows else 0):
        pivot = None
        mask = 1 << col
        for i in range(rank, len(rows)):
            if rows[i] & mask:
                pivot = i
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] & mask):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def extended_bch_check_rows(m: int, delta: int) -> tuple[list[int], int, int, int]:
    n, dim, redundancy = bch_dimension(m, delta)
    poly = find_primitive_poly(m)
    rows = parity_rows(m, delta, poly)
    ext_n = n + 1
    rows = rows + [(1 << ext_n) - 1]
    rank = row_rank(rows)
    return rows, ext_n, dim, rank


def syndrome(rows: list[int], x: int) -> int:
    out = 0
    for i, row in enumerate(rows):
        if (row & x).bit_count() & 1:
            out |= 1 << i
    return out


def log2_cdf_bin_half(width: int, cutoff: int) -> float:
    if cutoff < 0:
        return float("-inf")
    if cutoff >= width:
        return 0.0
    total = sum(math.comb(width, j) for j in range(cutoff + 1))
    return math.log2(total) - width


def iter_masks(n: int, weight: int):
    for combo in itertools.combinations(range(n), weight):
        x = 0
        for bit in combo:
            x |= 1 << bit
        yield x


def total_weight_table(rows: list[int], ext_n: int, rank: int, h_max: int, low_cuts: list[int]) -> None:
    print()
    print("total-weight syndrome table")
    header = ["r", "patterns", "zero", "log2_zero", "rand_zero_log2"]
    for cut in low_cuts:
        header.extend([f"le{cut}", f"log2_le{cut}", f"rand_le{cut}_log2"])
    print(",".join(header))
    for r in range(1, h_max + 1):
        counts = {cut: 0 for cut in low_cuts}
        zero = 0
        total = 0
        for x in iter_masks(ext_n, r):
            total += 1
            sw = syndrome(rows, x).bit_count()
            if sw == 0:
                zero += 1
            for cut in low_cuts:
                if sw <= cut:
                    counts[cut] += 1
        rand_zero = log2_binom(ext_n, r) - rank
        fields = [
            str(r),
            str(total),
            str(zero),
            f"{math.log2(zero):.6f}" if zero else "-inf",
            f"{rand_zero:.6f}",
        ]
        for cut in low_cuts:
            val = counts[cut]
            rand = log2_binom(ext_n, r) + log2_cdf_bin_half(rank, cut)
            fields.extend([
                str(val),
                f"{math.log2(val):.6f}" if val else "-inf",
                f"{rand:.6f}",
            ])
        print(",".join(fields))


def sampled_total_weight_table(
    rows: list[int],
    ext_n: int,
    rank: int,
    weights: list[int],
    low_cuts: list[int],
    samples: int,
    seed: int,
) -> None:
    rng = random.Random(seed)
    print()
    print("sampled total-weight syndrome table")
    header = ["r", "samples", "zero", "zero_rate_log2", "mean_syndrome_w"]
    for cut in low_cuts:
        header.extend([f"le{cut}", f"rate_le{cut}_log2", f"rand_le{cut}_log2"])
    print(",".join(header))
    for r in weights:
        counts = {cut: 0 for cut in low_cuts}
        zero = 0
        syndrome_sum = 0
        for _ in range(samples):
            x = 0
            for bit in rng.sample(range(ext_n), r):
                x |= 1 << bit
            sw = syndrome(rows, x).bit_count()
            syndrome_sum += sw
            if sw == 0:
                zero += 1
            for cut in low_cuts:
                if sw <= cut:
                    counts[cut] += 1
        fields = [
            str(r),
            str(samples),
            str(zero),
            f"{math.log2(zero / samples):.6f}" if zero else "-inf",
            f"{syndrome_sum / samples:.6f}",
        ]
        for cut in low_cuts:
            val = counts[cut]
            rand = log2_cdf_bin_half(rank, cut)
            fields.extend([
                str(val),
                f"{math.log2(val / samples):.6f}" if val else "-inf",
                f"{rand:.6f}",
            ])
        print(",".join(fields))


def split_weight_table(rows: list[int], block_bits: int, rank: int, h_max: int, low_cut: int) -> None:
    print()
    print("split-weight syndrome table")
    print("state_w,input_w,total,zero,le_cut,cut,mean_syndrome_w")
    left_masks = [[0] for _ in range(h_max + 1)]
    right_masks = [[0] for _ in range(h_max + 1)]
    for w in range(1, h_max + 1):
        left_masks[w] = list(iter_masks(block_bits, w))
        right_masks[w] = [x << block_bits for x in iter_masks(block_bits, w)]

    for state_w in range(0, h_max + 1):
        for input_w in range(0, h_max + 1 - state_w):
            if state_w + input_w == 0:
                continue
            total = 0
            zero = 0
            low = 0
            syndrome_sum = 0
            for lm in left_masks[state_w]:
                for rm in right_masks[input_w]:
                    total += 1
                    sw = syndrome(rows, lm | rm).bit_count()
                    syndrome_sum += sw
                    if sw == 0:
                        zero += 1
                    if sw <= low_cut:
                        low += 1
            print(
                f"{state_w},{input_w},{total},{zero},{low},{low_cut},"
                f"{syndrome_sum / total:.6f}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--delta", type=int, default=7)
    parser.add_argument("--h-max", type=int, default=8)
    parser.add_argument("--low-cuts", default="0,1,2,4")
    parser.add_argument("--split", action="store_true")
    parser.add_argument("--split-cut", type=int, default=4)
    parser.add_argument(
        "--sample-weights",
        default=None,
        help="Comma-separated total weights to sample instead of exhaustive enumeration.",
    )
    parser.add_argument("--samples-per-weight", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    rows, ext_n, dim, rank = extended_bch_check_rows(args.m, args.delta)
    block_bits = ext_n - dim
    low_cuts = [int(x) for x in args.low_cuts.split(",") if x.strip()]
    if ext_n != 2 * block_bits:
        print(
            "warning: extended BCH is not rate half, so ext_n != 2*redundancy "
            f"({ext_n} != {2 * block_bits})"
        )

    print("Transposed extended-BCH inner diagnostic")
    print(f"m={args.m}, primitive_n={ext_n - 1}, extended_n={ext_n}, designed_delta={args.delta}")
    print(f"extended_dimension={dim}, parity_rank={rank}, block_bits={block_bits}")
    print(f"rate_half_shape={ext_n == 2 * block_bits}")
    if ext_n == 2 * block_bits:
        left_rank = row_rank([row & ((1 << block_bits) - 1) for row in rows])
        right_rank = row_rank([(row >> block_bits) & ((1 << block_bits) - 1) for row in rows])
        print(f"left_half_rank={left_rank}, right_half_rank={right_rank}")
    if args.sample_weights is None:
        total_weight_table(rows, ext_n, rank, args.h_max, low_cuts)
    else:
        weights = [int(x) for x in args.sample_weights.split(",") if x.strip()]
        sampled_total_weight_table(
            rows,
            ext_n,
            rank,
            weights,
            low_cuts,
            args.samples_per_weight,
            args.seed,
        )
    if args.split and ext_n == 2 * block_bits:
        split_weight_table(rows, block_bits, rank, args.h_max, args.split_cut)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
