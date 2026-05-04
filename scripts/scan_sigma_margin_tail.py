#!/usr/bin/env python3
"""Scan sigma versus first-moment margin for the terminated-tail model."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

from certify_global_episode_cover import (
    binom_cdf_half_be_log2,
    global_episode_inner_log2,
    survivor_only_exact_prefix_log2,
)
from dense_largek_eval import log2add, outer_small_h_prefix


def format_log2(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.6f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=2**20)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--sigma-start", type=int, default=24)
    parser.add_argument("--sigma-step", type=int, default=8)
    parser.add_argument("--sigma-max", type=int, default=200)
    parser.add_argument("--target-bits", type=float, default=40.0)
    parser.add_argument("--h-max", type=int, default=30)
    parser.add_argument("--r-max", type=int, default=12)
    parser.add_argument("--block-ratio", type=float, default=1.02)
    parser.add_argument("--suffix-block-ratio", type=float, default=1.1)
    parser.add_argument("--out-prefix", default=None)
    args = parser.parse_args()

    if args.out_prefix is None:
        out_base = Path(__file__).with_name(
            f"sigma_margin_tail_k{args.k}_delta{args.delta:g}_step{args.sigma_step}_h{args.h_max}"
        )
    else:
        out_base = Path(args.out_prefix)
        if not out_base.is_absolute():
            out_base = Path.cwd() / out_base

    rows: list[dict[str, str | int | float]] = []
    for sigma in range(args.sigma_start, args.sigma_max + 1, args.sigma_step):
        parity_n = args.k + sigma
        n = args.k + parity_n
        d = math.floor(args.delta * n)
        print(f"sigma={sigma}: N={n}, d={d}", flush=True)

        tail_cache = [binom_cdf_half_be_log2(u, d) for u in range(n + 1)]
        exact_survivor = survivor_only_exact_prefix_log2(n=n, h_max=args.h_max, tail_cache=tail_cache)
        outer = outer_small_h_prefix(
            args.k,
            sigma,
            args.h_max,
            outer_mode="banded-fixedtap",
            parity_n=parity_n,
        )

        total = float("-inf")
        peak_h = -1
        peak_term = float("-inf")
        peak_outer = float("-inf")
        peak_inner = float("-inf")
        for h in range(1, args.h_max + 1):
            inner, _, _, _ = global_episode_inner_log2(
                n=n,
                h=h,
                sigma=sigma,
                delta=args.delta,
                r_max=args.r_max,
                block_ratio=args.block_ratio,
                suffix_block_ratio=args.suffix_block_ratio,
                stop_gap_bits=60,
                tail_confirm=8,
                suffix_survivor=True,
                tail_mode="be",
                exact_survivor=True,
                tail_cache=tail_cache,
                exact_survivor_value=exact_survivor[h],
            )
            term = outer[h] + inner if outer[h] != float("-inf") and inner != float("-inf") else float("-inf")
            total = log2add(total, term)
            if term > peak_term:
                peak_h = h
                peak_term = term
                peak_outer = outer[h]
                peak_inner = inner

        margin = -total
        rows.append(
            {
                "sigma": sigma,
                "offset": sigma - math.ceil(math.log2(args.k)),
                "k": args.k,
                "N": n,
                "delta": args.delta,
                "d": d,
                "h_max": args.h_max,
                "r_max": args.r_max,
                "total_log2_hle": format_log2(total),
                "margin_bits_hle": format_log2(margin),
                "peak_h": peak_h,
                "peak_term_log2": format_log2(peak_term),
                "peak_outer_log2": format_log2(peak_outer),
                "peak_inner_log2": format_log2(peak_inner),
            }
        )
        print(f"  margin={margin:.3f}, peak_h={peak_h}, peak={peak_term:.3f}", flush=True)
        if margin >= args.target_bits:
            break

    csv_path = out_base.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {csv_path}")

    xs = [int(row["sigma"]) for row in rows]
    ys = [float(row["margin_bits_hle"]) for row in rows]
    peaks = [int(row["peak_h"]) for row in rows]
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.plot(xs, ys, marker="o", linewidth=2)
    ax.axhline(args.target_bits, color="black", linestyle="--", linewidth=1)
    for x, y, h in zip(xs, ys, peaks):
        ax.annotate(f"h={h}", (x, y), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
    ax.set_xlabel("sigma = m")
    ax.set_ylabel(f"-log2 first moment, h <= {args.h_max} (bits)")
    ax.set_title(f"Terminated-tail dense+dense margin at delta={args.delta:g}, k={args.k}")
    ax.grid(True, alpha=0.3)
    png_path = out_base.with_suffix(".png")
    fig.savefig(png_path, dpi=180)
    print(f"wrote {png_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
