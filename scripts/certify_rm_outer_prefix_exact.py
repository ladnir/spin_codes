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


@dataclass(frozen=True)
class CoefficientBuildTrace:
    exponent: int
    multiply_steps: int
    square_steps: int
    max_coefficient_bits: int
    max_coefficient_weight: int


@dataclass(frozen=True)
class PrefixCoefficientCertificate:
    coeffs: list[int]
    blocks: int
    h_max: int
    local_nonzero_terms: int
    local_min_positive: int
    local_first_positive_le_80: list[int]
    trace: CoefficientBuildTrace
    support_count: int
    positive_count: int
    min_positive: int
    next_positive_after_min: int
    first_positive_le_80: list[int]
    zero_prefix_stop: int
    gap_after_min_start: int
    gap_after_min_stop: int
    split_h: int
    split_coefficient_bits: int
    next_coefficient_bits: int


def direct_sum_coefficients_with_trace(
    local: list[int],
    *,
    blocks: int,
    h_max: int,
) -> tuple[list[int], CoefficientBuildTrace]:
    poly = [1] + [0] * h_max
    base = local
    exponent = blocks
    multiply_steps = 0
    square_steps = 0
    while exponent:
        if exponent & 1:
            poly = truncated_convolution(poly, base, h_max)
            multiply_steps += 1
        exponent >>= 1
        if exponent:
            base = truncated_convolution(base, base, h_max)
            square_steps += 1
    max_weight, max_coeff = max(enumerate(poly), key=lambda item: item[1].bit_length())
    trace = CoefficientBuildTrace(
        exponent=blocks,
        multiply_steps=multiply_steps,
        square_steps=square_steps,
        max_coefficient_bits=max_coeff.bit_length(),
        max_coefficient_weight=max_weight,
    )
    return poly, trace


def direct_sum_coefficients(local: list[int], *, blocks: int, h_max: int) -> list[int]:
    coeffs, _trace = direct_sum_coefficients_with_trace(local, blocks=blocks, h_max=h_max)
    return coeffs


def prefix_coefficient_certificate(
    local: list[int],
    *,
    blocks: int,
    h_max: int,
) -> PrefixCoefficientCertificate:
    coeffs, trace = direct_sum_coefficients_with_trace(local, blocks=blocks, h_max=h_max)
    positive = positive_support(coeffs, h_max)
    if len(positive) < 2:
        raise SystemExit("exact RM prefix: need at least two positive weights")
    local_positive = positive_support(local, h_max)
    if not local_positive:
        raise SystemExit("exact RM prefix: local spectrum has no positive support")
    return PrefixCoefficientCertificate(
        coeffs=coeffs,
        blocks=blocks,
        h_max=h_max,
        local_nonzero_terms=sum(1 for coeff in local if coeff),
        local_min_positive=local_positive[0],
        local_first_positive_le_80=[h for h in local_positive if h <= 80],
        trace=trace,
        support_count=sum(1 for coeff in coeffs if coeff),
        positive_count=len(positive),
        min_positive=positive[0],
        next_positive_after_min=positive[1],
        first_positive_le_80=[h for h in positive if h <= 80],
        zero_prefix_stop=positive[0] - 1,
        gap_after_min_start=positive[0] + 1,
        gap_after_min_stop=positive[1] - 1,
        split_h=positive[0],
        split_coefficient_bits=coeffs[positive[0]].bit_length(),
        next_coefficient_bits=coeffs[positive[1]].bit_length(),
    )


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
    cert = prefix_coefficient_certificate(local, blocks=args.blocks, h_max=args.h_max)
    coeffs = cert.coeffs
    print("Exact RM outer prefix")
    print(f"local_spectrum_csv,{args.local_spectrum_csv}")
    print(f"blocks,{args.blocks}")
    print(f"h_max,{args.h_max}")
    print(f"local_nonzero_terms,{cert.local_nonzero_terms}")
    print(f"local_min_positive,{cert.local_min_positive}")
    print(f"local_first_positive,{';'.join(str(h) for h in cert.local_first_positive_le_80)}")
    print(f"coefficient_build_exponent,{cert.trace.exponent}")
    print(f"coefficient_build_multiply_steps,{cert.trace.multiply_steps}")
    print(f"coefficient_build_square_steps,{cert.trace.square_steps}")
    print(f"coefficient_build_max_coefficient_bits,{cert.trace.max_coefficient_bits}")
    print(f"coefficient_build_max_coefficient_weight,{cert.trace.max_coefficient_weight}")
    print(f"support_count,{cert.support_count}")
    print(f"support_positive_count,{cert.positive_count}")
    print(f"support_min_positive,{cert.min_positive}")
    print(f"support_next_positive_after_min,{cert.next_positive_after_min}")
    print(f"support_zero_prefix,1--{cert.zero_prefix_stop}")
    print(f"support_gap_after_min,{cert.gap_after_min_start}--{cert.gap_after_min_stop}")
    print(f"support_first_positive,{';'.join(str(h) for h in cert.first_positive_le_80)}")
    print(f"support_split_h,{cert.split_h}")
    print(f"support_split_coefficient_bits,{cert.split_coefficient_bits}")
    print(f"support_next_coefficient_bits,{cert.next_coefficient_bits}")

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
