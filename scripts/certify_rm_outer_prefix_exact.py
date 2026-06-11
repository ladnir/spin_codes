#!/usr/bin/env python3
"""Exact RM direct-sum outer prefix coefficients and ledger reweighting.

The RM block outer is a direct sum of identical local blocks.  For the finite
``h <= 500`` prefix we can compute exact global coefficients of

    W_local(z)^blocks

by truncated polynomial exponentiation.  This removes the Cauchy/GF smoothing
artifact that assigns mass to off-support weights such as ``h=33``.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def log2add(a: float, b: float) -> float:
    if a == float("-inf"):
        return b
    if b == float("-inf"):
        return a
    if a < b:
        a, b = b, a
    return a + math.log2(1.0 + 2.0 ** (b - a))


def log2_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2.0)


def load_local_spectrum(path: Path, h_max: int) -> list[int]:
    coeffs = [0] * (h_max + 1)
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            weight = int(row["weight"])
            if weight <= h_max:
                coeffs[weight] = int(row["count"])
    return coeffs


def truncated_convolution(a: list[int], b: list[int], h_max: int) -> list[int]:
    out = [0] * (h_max + 1)
    a_items = [(i, value) for i, value in enumerate(a) if value]
    b_items = [(i, value) for i, value in enumerate(b) if value]
    for i, x in a_items:
        limit = h_max - i
        for j, y in b_items:
            if j > limit:
                break
            out[i + j] += x * y
    return out


def direct_sum_coefficients(local: list[int], *, blocks: int, h_max: int) -> list[int]:
    poly = [1] + [0] * h_max
    base = local
    exponent = blocks
    while exponent:
        if exponent & 1:
            poly = truncated_convolution(poly, base, h_max)
        exponent >>= 1
        if exponent:
            base = truncated_convolution(base, base, h_max)
    return poly


@dataclass(frozen=True)
class ReweightedLedgerSummary:
    rows: int
    live_rows: int
    cauchy_log2: float
    exact_log2: float
    split_h: int
    split_log2: float
    above_split_log2: float
    peak_term_log2: float
    peak_h: int
    peak_first_r: int
    peak_gap_min: int
    peak_gap_max: int
    peak_outer_exact_log2: float


def exact_outer_log2(coeffs: list[int], h: int) -> float:
    if h < 0 or h >= len(coeffs) or coeffs[h] == 0:
        return float("-inf")
    return math.log2(coeffs[h])


def positive_support(coeffs: list[int], h_max: int | None = None) -> list[int]:
    if h_max is None:
        h_max = len(coeffs) - 1
    return [h for h, coeff in enumerate(coeffs[: h_max + 1]) if h > 0 and coeff]


def old_outer_column(row: dict[str, str]) -> str:
    if "outer_log2_bound" in row:
        return "outer_log2_bound"
    if "outer_log2" in row:
        return "outer_log2"
    raise KeyError("row has neither outer_log2_bound nor outer_log2")


def reweighted_ledger_sum(path: Path, coeffs: list[int], *, split_h: int | None = None) -> ReweightedLedgerSummary:
    if split_h is None:
        support = positive_support(coeffs)
        if not support:
            raise SystemExit("exact outer coefficients have no positive support")
        split_h = support[0]
    rows = 0
    live_rows = 0
    cauchy_total = float("-inf")
    exact_total = float("-inf")
    split_total = float("-inf")
    above_split_total = float("-inf")
    peak: tuple[float, dict[str, str], float] | None = None
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows += 1
            cauchy_term = float(row["term_log2"])
            cauchy_total = log2add(cauchy_total, cauchy_term)
            h = int(row["outer_weight"])
            outer = exact_outer_log2(coeffs, h)
            if outer == float("-inf"):
                continue
            column = old_outer_column(row)
            old_outer = float(row[column])
            exact_term = cauchy_term - old_outer + outer
            exact_total = log2add(exact_total, exact_term)
            if h == split_h:
                split_total = log2add(split_total, exact_term)
            elif h > split_h:
                above_split_total = log2add(above_split_total, exact_term)
            live_rows += 1
            if peak is None or exact_term > peak[0]:
                peak = (exact_term, row, outer)
    if peak is None:
        raise SystemExit(f"{path}: no live exact-outer rows")
    peak_term, peak_row, peak_outer = peak
    return ReweightedLedgerSummary(
        rows=rows,
        live_rows=live_rows,
        cauchy_log2=cauchy_total,
        exact_log2=exact_total,
        split_h=split_h,
        split_log2=split_total,
        above_split_log2=above_split_total,
        peak_term_log2=peak_term,
        peak_h=int(peak_row["outer_weight"]),
        peak_first_r=int(peak_row["first_r"]),
        peak_gap_min=int(peak_row["gap_min"]),
        peak_gap_max=int(peak_row["gap_max"]),
        peak_outer_exact_log2=peak_outer,
    )


@dataclass(frozen=True)
class LatePrefixSummary:
    total_log2: float
    rows: int
    split_h: int
    split_log2: float
    above_split_log2: float
    peak_h: int
    peak_outer_log2: float
    peak_late_log2: float
    peak_term_log2: float


def exact_late_prefix_sum(
    coeffs: list[int],
    *,
    n: int,
    b: int,
    late_blocks: int,
    h_max: int,
    split_h: int | None = None,
) -> LatePrefixSummary:
    if split_h is None:
        support = positive_support(coeffs, h_max)
        if not support:
            raise SystemExit("exact late prefix: no positive support")
        split_h = support[0]
    total = float("-inf")
    split_total = float("-inf")
    above_split_total = float("-inf")
    rows = 0
    peak: tuple[float, int, float, float] | None = None
    for h in range(1, min(h_max, len(coeffs) - 1) + 1):
        outer = exact_outer_log2(coeffs, h)
        if outer == float("-inf"):
            continue
        late = log2_binom(late_blocks * b, h) - log2_binom(n, h)
        term = outer + late
        total = log2add(total, term)
        if h == split_h:
            split_total = log2add(split_total, term)
        elif h > split_h:
            above_split_total = log2add(above_split_total, term)
        rows += 1
        if peak is None or term > peak[0]:
            peak = (term, h, outer, late)
    if peak is None:
        raise SystemExit("exact late prefix: no live rows")
    peak_term, peak_h, peak_outer, peak_late = peak
    return LatePrefixSummary(
        total_log2=total,
        rows=rows,
        split_h=split_h,
        split_log2=split_total,
        above_split_log2=above_split_total,
        peak_h=peak_h,
        peak_outer_log2=peak_outer,
        peak_late_log2=peak_late,
        peak_term_log2=peak_term,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-spectrum-csv", type=Path, default=ROOT / "rm512_256_spectrum.csv")
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--h-max", type=int, default=500)
    parser.add_argument(
        "--ledger-csv",
        action="append",
        type=Path,
        default=[],
        help="Ledger CSV to reweight from Cauchy/GF outer rows to exact outer coefficients.",
    )
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--late-blocks", type=int, default=5949)
    args = parser.parse_args()

    local = load_local_spectrum(args.local_spectrum_csv, args.h_max)
    coeffs = direct_sum_coefficients(local, blocks=args.blocks, h_max=args.h_max)
    support = [h for h, coeff in enumerate(coeffs) if coeff]
    positive = positive_support(coeffs, args.h_max)
    print("Exact RM outer prefix")
    print(f"local_spectrum_csv,{args.local_spectrum_csv}")
    print(f"blocks,{args.blocks}")
    print(f"h_max,{args.h_max}")
    print(f"support_count,{len(support)}")
    print(f"support_positive_count,{len(positive)}")
    print(f"support_min_positive,{positive[0]}")
    if len(positive) > 1:
        print(f"support_next_positive_after_min,{positive[1]}")
    print(f"support_first_positive,{';'.join(str(h) for h in positive if h <= 80)}")

    late = exact_late_prefix_sum(
        coeffs,
        n=args.N,
        b=args.block_bits,
        late_blocks=args.late_blocks,
        h_max=args.h_max,
    )
    print(f"late_prefix_exact_rows,{late.rows}")
    print(f"late_prefix_exact_log2,{late.total_log2:.6f}")
    print(f"late_prefix_exact_split_h,{late.split_h}")
    print(f"late_prefix_exact_split_log2,{late.split_log2:.6f}")
    print(f"late_prefix_exact_above_split_log2,{late.above_split_log2:.6f}")
    print(
        "late_prefix_exact_peak,"
        f"h={late.peak_h},outer_log2={late.peak_outer_log2:.6f},"
        f"late_log2={late.peak_late_log2:.6f},term_log2={late.peak_term_log2:.6f}"
    )

    for path in args.ledger_csv:
        summary = reweighted_ledger_sum(path, coeffs)
        print(f"ledger,{path}")
        print(f"ledger_rows,{summary.rows}")
        print(f"ledger_live_rows,{summary.live_rows}")
        print(f"ledger_cauchy_log2,{summary.cauchy_log2:.6f}")
        print(f"ledger_exact_outer_log2,{summary.exact_log2:.6f}")
        print(f"ledger_exact_split_h,{summary.split_h}")
        print(f"ledger_exact_split_log2,{summary.split_log2:.6f}")
        print(f"ledger_exact_above_split_log2,{summary.above_split_log2:.6f}")
        print(
            "ledger_exact_peak,"
            f"h={summary.peak_h},r={summary.peak_first_r},"
            f"gap={summary.peak_gap_min}--{summary.peak_gap_max},"
            f"outer_log2={summary.peak_outer_exact_log2:.6f},"
            f"term_log2={summary.peak_term_log2:.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
