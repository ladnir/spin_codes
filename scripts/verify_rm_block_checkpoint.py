#!/usr/bin/env python3
"""Verify the RM(4,9) block-outer finite checkpoint artifacts.

This checker consumes the summary artifacts from

    block_outer_upgrade_probe.py --model spectrum-csv

and the geometric tail post-processor.  It verifies the published local
RM(4,9) spectrum was copied consistently, and that the finite prefix plus the
conditional post-prefix tail meet the requested first-moment margin.

The final tail is intentionally conditional: the script checks the observed
finite ratio and resulting residual, while the paper still needs the matching
post-prefix domination lemma.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


DEFAULT_SUMMARY = Path(__file__).with_name("block_outer_probe_k1048576_rm512_256_sig32_d009_h500_exact_summary.csv")
DEFAULT_TAILCERT = Path(__file__).with_name("block_outer_probe_k1048576_rm512_256_sig32_d009_h500_exact_tailcert.csv")
DEFAULT_SPECTRUM = Path(__file__).with_name("rm512_256_spectrum.csv")
DEFAULT_PREFIX = Path(__file__).with_name("block_outer_probe_k1048576_rm512_256_sig32_d009_h500_exact.csv")


def read_one(path: Path) -> dict[str, str]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 1:
        raise SystemExit(f"{path} must contain exactly one row, found {len(rows)}")
    return rows[0]


def read_spectrum_sum(path: Path) -> tuple[int, int, int]:
    total = 0
    min_nonzero = None
    max_weight = 0
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            weight = int(row["weight"])
            count = int(row["count"])
            total += count
            max_weight = max(max_weight, weight)
            if weight > 0 and count > 0:
                min_nonzero = weight if min_nonzero is None else min(min_nonzero, weight)
    if min_nonzero is None:
        raise SystemExit(f"{path} has no nonzero local codewords")
    return total, min_nonzero, max_weight


def parse_log2(text: str) -> float:
    text = text.strip()
    if text == "-inf":
        return float("-inf")
    return float(text)


def read_term(path: Path, h: int) -> float:
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if int(row["h"]) == h:
                return parse_log2(row["term_log2"])
    raise SystemExit(f"{path} has no row h={h}")


def log2add(x: float, y: float) -> float:
    if x == float("-inf"):
        return y
    if y == float("-inf"):
        return x
    hi = max(x, y)
    lo = min(x, y)
    return hi + math.log2(1.0 + 2.0 ** (lo - hi))


def check(name: str, passed: bool, detail: str) -> bool:
    print(f"{name}: {'PASS' if passed else 'FAIL'} ({detail})")
    return passed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--tailcert", type=Path, default=DEFAULT_TAILCERT)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--prefix-csv", type=Path, default=DEFAULT_PREFIX)
    parser.add_argument("--k", type=int, default=2**20)
    parser.add_argument("--n", type=int, default=2**21)
    parser.add_argument("--delta", type=float, default=0.09)
    parser.add_argument("--d", type=int, default=188743)
    parser.add_argument("--sigma", type=int, default=32)
    parser.add_argument("--block-bits", type=int, default=256)
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--local-length", type=int, default=512)
    parser.add_argument("--local-dim", type=int, default=256)
    parser.add_argument("--local-distance", type=int, default=32)
    parser.add_argument("--h-max", type=int, default=500)
    parser.add_argument("--peak-h", type=int, default=32)
    parser.add_argument("--max-total-log2", type=float, default=-40.5)
    parser.add_argument("--max-tail-log2", type=float, default=-500.0)
    parser.add_argument("--max-ratio-log2", type=float, default=-0.9)
    parser.add_argument("--analytic-ratio-log2", type=float, default=-0.119385964068)
    parser.add_argument("--max-analytic-tail-log2", type=float, default=-500.0)
    args = parser.parse_args()

    summary = read_one(args.summary)
    tail = read_one(args.tailcert)
    spectrum_total, spectrum_min, spectrum_max = read_spectrum_sum(args.spectrum)

    total = float(summary["total_log2"])
    tail_residual = float(tail["tail_residual_log2"])
    certified_total = float(tail["certified_total_log2_if_ratio_holds"])
    observed_ratio = float(tail["observed_worst_log2_ratio"])
    tail_after = int(tail["tail_after"])
    last_term = read_term(args.prefix_csv, tail_after)
    q = 2.0 ** args.analytic_ratio_log2
    analytic_tail = last_term + args.analytic_ratio_log2 - math.log2(1.0 - q)
    total_with_tail = log2add(total, tail_residual)
    total_with_analytic_tail = log2add(total, analytic_tail)

    print("RM(4,9) block-outer checkpoint verification")
    print(f"summary = {args.summary}")
    print(f"tailcert = {args.tailcert}")
    print(f"spectrum = {args.spectrum}")
    print(f"prefix csv = {args.prefix_csv}")
    print(f"prefix total log2 = {total:.6f}")
    print(f"tail residual log2 = {tail_residual:.6f}")
    print(f"analytic-ratio tail residual log2 = {analytic_tail:.6f}")
    print(f"certified total log2 = {certified_total:.6f}")
    print(f"prefix plus residual log2 = {total_with_tail:.6f}")
    print(f"prefix plus analytic-ratio residual log2 = {total_with_analytic_tail:.6f}")
    print(f"observed finite tail ratio log2 = {observed_ratio:.6f}")
    print()

    ok = True
    ok &= check("spectrum total", spectrum_total == (1 << args.local_dim), f"sum={spectrum_total}")
    ok &= check("spectrum min distance", spectrum_min == args.local_distance, f"d0={spectrum_min}")
    ok &= check("spectrum length", spectrum_max == args.local_length, f"max_weight={spectrum_max}")
    ok &= check("k", summary["k"] == str(args.k), summary["k"])
    ok &= check("k_eff", summary["k_eff"] == str(args.k), summary["k_eff"])
    ok &= check("N", summary["N"] == str(args.n), summary["N"])
    ok &= check("delta", abs(float(summary["delta"]) - args.delta) <= 1e-15, summary["delta"])
    ok &= check("d", summary["d"] == str(args.d), summary["d"])
    ok &= check("sigma", summary["sigma"] == str(args.sigma), summary["sigma"])
    ok &= check("block bits", summary["block_bits"] == str(args.block_bits), summary["block_bits"])
    ok &= check("blocks", summary["blocks"] == str(args.blocks), summary["blocks"])
    ok &= check("local d0", summary["d0"] == str(args.local_distance), summary["d0"])
    ok &= check("h_max", summary["h_max"] == str(args.h_max), summary["h_max"])
    ok &= check("model", summary["model"] == "spectrum-csv", summary["model"])
    ok &= check("peak", summary["peak_h"] == str(args.peak_h), summary["peak_h"])
    ok &= check("prefix margin", total <= args.max_total_log2, f"{total:.6f} <= {args.max_total_log2:.6f}")
    ok &= check(
        "tail ratio",
        observed_ratio <= args.max_ratio_log2,
        f"{observed_ratio:.6f} <= {args.max_ratio_log2:.6f}",
    )
    ok &= check(
        "tail residual",
        tail_residual <= args.max_tail_log2,
        f"{tail_residual:.6f} <= {args.max_tail_log2:.6f}",
    )
    ok &= check(
        "analytic-ratio tail residual",
        analytic_tail <= args.max_analytic_tail_log2,
        f"{analytic_tail:.6f} <= {args.max_analytic_tail_log2:.6f}",
    )
    ok &= check(
        "certified total",
        certified_total <= args.max_total_log2
        and total_with_tail <= args.max_total_log2
        and total_with_analytic_tail <= args.max_total_log2,
        f"{certified_total:.6f}, {total_with_tail:.6f}, {total_with_analytic_tail:.6f} <= {args.max_total_log2:.6f}",
    )
    print(f"STATUS: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
