#!/usr/bin/env python3
"""Verify the current finite full-split dense-inner first-moment ledger.

This is a manifest-style checker for the finite ``N=2^21, delta=.09`` RM/full
split certificate recorded in innerDense.tex.  The tiny prefix is still backed
by explicit finite episode CSVs, while the post-prefix and high tail use the
theorem-facing all-episode wrapper.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from certify_prefix_placement_ratio import compute_sums, decimal_fraction, ratio_float
from certify_rm_outer_prefix_exact import (
    direct_sum_coefficients,
    exact_late_prefix_sum,
    load_local_spectrum,
    reweighted_ledger_sum,
)


ROOT = Path(__file__).resolve().parent


GLOBAL_TURNOFF_LOG2 = -62.4078758
HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS = 1.484774


@dataclass(frozen=True)
class HighInterval:
    start: int
    stop: int
    lam: float
    rho: float
    raw_log2: float

    @property
    def label(self) -> str:
        return f"{self.start}--{self.stop}"

    @property
    def adjusted_log2(self) -> float:
        return self.raw_log2 + HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS

    def command(self) -> list[str]:
        return [
            sys.executable,
            str(ROOT / "sum_fullsplit_piecewise_certificate.py"),
            "--h-values",
            f"{self.start}:{self.stop}",
            "--inner-mode-by-gap",
            "cap,eallratio,eallratio",
            "--gap-sum-mode",
            "endpoint",
            "--inner-T-by-gap",
            "feasible-min",
            "--turnoff-log2",
            f"{GLOBAL_TURNOFF_LOG2:.7f}",
            "--lambdas",
            f"{self.lam:g}",
            "--rhos",
            f"{self.rho:g}",
        ]

    def display_command(self) -> str:
        script = r"scripts\sum_fullsplit_piecewise_certificate.py"
        return (
            f"python {script} --h-values {self.start}:{self.stop} "
            "--inner-mode-by-gap cap,eallratio,eallratio --gap-sum-mode endpoint "
            f"--inner-T-by-gap feasible-min --turnoff-log2 {GLOBAL_TURNOFF_LOG2:.7f} "
            f"--lambdas {self.lam:g} --rhos {self.rho:g}"
        )


@dataclass(frozen=True)
class ComplementHighInterval:
    start: int
    stop: int
    rho: float
    log2_value: float

    @property
    def label(self) -> str:
        return f"{self.start}--{self.stop}"

    def command(self) -> list[str]:
        return [
            sys.executable,
            str(ROOT / "certify_fullsplit_high_density_interval.py"),
            "--intervals",
            f"{self.start}:{self.stop}:{self.rho:g}",
            "--lambda-value",
            "2.3",
        ]

    def display_command(self) -> str:
        script = r"scripts\certify_fullsplit_high_density_interval.py"
        return (
            f"python {script} --intervals {self.start}:{self.stop}:{self.rho:g} "
            "--lambda-value 2.3"
        )


HIGH_INTERVALS_RAW = [
    # Raw interval rows generated with the sharper finite-prefix p_term.  The
    # theorem-facing wrapper uses the global Bernstein p_term; the checked
    # adjustment above covers this wider atom.
    HighInterval(2001, 7858, 0.02, 0.01, -588.081877),
    HighInterval(7859, 20550, 0.05, 0.03, -2657.429602),
    HighInterval(20551, 75000, 0.2, 0.1, -6027.058821),
    HighInterval(75001, 250000, 0.5, 0.3, -12846.353611),
    HighInterval(250001, 350000, 1.2, 0.5, -28478.125315),
    HighInterval(350001, 400000, 1.2, 1.0, -81557.079312),
    HighInterval(400001, 450000, 1.2, 1.0, -91423.257760),
    HighInterval(450001, 550000, 1.2, 1.0, -82083.446102),
    HighInterval(550001, 650000, 1.2, 1.0, -1519.978235),
    HighInterval(650001, 725000, 1.2, 2.0, -116329.618720),
    HighInterval(725001, 750000, 1.2, 3.0, -169889.387695),
    HighInterval(750001, 850000, 1.2, 3.0, -134525.068420),
    HighInterval(850001, 950000, 1.2, 3.0, -117536.328304),
    HighInterval(950001, 1050000, 1.2, 3.0, -285054.035442),
    HighInterval(1050001, 1148736, 1.2, 3.0, -562862.201412),
]


HIGH_INTERVALS = [
    (interval.label, interval.lam, interval.rho, interval.adjusted_log2)
    for interval in HIGH_INTERVALS_RAW
]


COMPLEMENT_HIGH_INTERVALS = [
    ComplementHighInterval(1048577, 1148736, 1.0, -128915.315746),
    ComplementHighInterval(1148737, 1300000, 1.0, -124889.651755),
    ComplementHighInterval(1300001, 1500000, 1.0, -113710.357653),
    ComplementHighInterval(1500001, 1700000, 1.0, -35150.097366),
    ComplementHighInterval(1700001, 1900000, 10.0, -524628.366199),
    ComplementHighInterval(1900001, 2097089, 10.0, -894140.350713),
    ComplementHighInterval(2097090, 2097152, 10.0, -893792.285010),
]


@dataclass(frozen=True)
class CsvLedger:
    name: str
    path: Path
    column: str
    expected_log2: float


@dataclass(frozen=True)
class InteriorAuditSummary:
    rows: int
    worst: dict[str, str]
    worst_repro: dict[str, str]


@dataclass(frozen=True)
class PlacementRatioSummary:
    peak_h: int
    peak_ratio: float
    endpoint_ratio_at_h_min: float
    endpoint_slack_factor: float


def log2_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2.0)


def log2add(a: float, b: float) -> float:
    if a == float("-inf"):
        return b
    if b == float("-inf"):
        return a
    if a < b:
        a, b = b, a
    return a + math.log2(1.0 + 2.0 ** (b - a))


def late_prefix_cap_log2(
    path: Path,
    *,
    n: int,
    b: int,
    late_blocks: int,
    h_max: int,
) -> tuple[float, int, dict[str, float] | None]:
    total = float("-inf")
    rows = 0
    peak: dict[str, float] | None = None
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if "h" not in row or "outer_log2" not in row:
                continue
            outer_text = row["outer_log2"]
            if not outer_text or outer_text == "-inf":
                continue
            h = int(row["h"])
            if h < 1 or h > h_max:
                continue
            outer = float(outer_text)
            late = log2_binom(late_blocks * b, h) - log2_binom(n, h)
            term = outer + late
            rows += 1
            total = log2add(total, term)
            if peak is None or term > peak["term_log2"]:
                peak = {
                    "h": float(h),
                    "outer_log2": outer,
                    "late_log2": late,
                    "term_log2": term,
                }
    return total, rows, peak


def csv_logsum(path: Path, column: str) -> tuple[float, int, dict[str, str] | None]:
    total = float("-inf")
    rows = 0
    peak: dict[str, str] | None = None
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows += 1
            value = float(row[column])
            total = log2add(total, value)
            if peak is None or value > float(peak[column]):
                peak = row
    return total, rows, peak


@dataclass(frozen=True)
class PrefixRidgeRatioSummary:
    rows: int
    ridge_rows: int
    ridge_exact_log2: float
    remainder_exact_log2: float
    peak_log2: float
    first_ratio: float
    tail_ratio: float
    geometric_log2: float
    total_bound_log2: float
    total_exact_log2: float
    inner_span_bits: float
    first_outer_ratio: float
    first_placement_ratio: float
    first_inner_ratio: float
    tail_outer_ratio: float
    tail_placement_ratio: float
    tail_inner_ratio: float


def prefix_ridge_ratio_summary(
    path: Path,
    *,
    gap_min: int = 1,
    gap_max: int = 4000,
    first_r: int = 1,
) -> PrefixRidgeRatioSummary:
    rows = 0
    total = float("-inf")
    ridge_exact = float("-inf")
    remainder_exact = float("-inf")
    ridge: list[tuple[int, float, float, float, float]] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows += 1
            h = int(row["outer_weight"])
            r = int(row["first_r"])
            g0 = int(row["gap_min"])
            g1 = int(row["gap_max"])
            term = float(row["term_log2"])
            total = log2add(total, term)
            if r == first_r and g0 == gap_min and g1 == gap_max:
                ridge_exact = log2add(ridge_exact, term)
                ridge.append(
                    (
                        h,
                        term,
                        float(row["outer_log2_bound"]),
                        float(row["placement_log2"]),
                        float(row["inner_log2"]),
                    )
                )
            else:
                remainder_exact = log2add(remainder_exact, term)

    if len(ridge) < 2:
        raise SystemExit("prefix ridge ratio: need at least two ridge rows")
    ridge.sort(key=lambda item: item[0])
    for (prev_h, *_), (h, *_) in zip(ridge, ridge[1:]):
        if h != prev_h + 1:
            raise SystemExit(f"prefix ridge ratio: non-consecutive h values {prev_h}, {h}")

    peak_log2 = ridge[0][1]
    first_ratio = 2.0 ** (ridge[1][1] - peak_log2)
    tail_ratio = max(2.0 ** (ridge[i][1] - ridge[i - 1][1]) for i in range(2, len(ridge)))
    first_outer_ratio = 2.0 ** (ridge[1][2] - ridge[0][2])
    first_placement_ratio = 2.0 ** (ridge[1][3] - ridge[0][3])
    first_inner_ratio = 2.0 ** (ridge[1][4] - ridge[0][4])
    tail_outer_ratio = max(2.0 ** (ridge[i][2] - ridge[i - 1][2]) for i in range(2, len(ridge)))
    tail_placement_ratio = max(2.0 ** (ridge[i][3] - ridge[i - 1][3]) for i in range(2, len(ridge)))
    tail_inner_ratio = max(2.0 ** (ridge[i][4] - ridge[i - 1][4]) for i in range(2, len(ridge)))
    factor = 1.0 + first_ratio / (1.0 - tail_ratio)
    geometric = peak_log2 + math.log2(factor)
    total_bound = log2add(geometric, remainder_exact)
    inner_values = [inner for *_, inner in ridge]
    return PrefixRidgeRatioSummary(
        rows=rows,
        ridge_rows=len(ridge),
        ridge_exact_log2=ridge_exact,
        remainder_exact_log2=remainder_exact,
        peak_log2=peak_log2,
        first_ratio=first_ratio,
        tail_ratio=tail_ratio,
        geometric_log2=geometric,
        total_bound_log2=total_bound,
        total_exact_log2=total,
        inner_span_bits=max(inner_values) - min(inner_values),
        first_outer_ratio=first_outer_ratio,
        first_placement_ratio=first_placement_ratio,
        first_inner_ratio=first_inner_ratio,
        tail_outer_ratio=tail_outer_ratio,
        tail_placement_ratio=tail_placement_ratio,
        tail_inner_ratio=tail_inner_ratio,
    )


def exact_placement_ratio_summary(
    *,
    n: int,
    b: int,
    late_blocks: int,
    gap_min: int,
    gap_max: int,
    h_min: int,
    h_max: int,
    threshold_text: str,
) -> PlacementRatioSummary:
    threshold = decimal_fraction(threshold_text)
    sums = compute_sums(
        b=b,
        late_blocks=late_blocks,
        gap_min=gap_min,
        gap_max=gap_max,
        h_max=h_max,
    )
    max_num = 0
    max_den = 1
    max_h = -1
    for h in range(h_min, h_max):
        numerator = sums[h] * (h + 1)
        denominator = sums[h - 1] * (n - h)
        if numerator * threshold.denominator > denominator * threshold.numerator:
            raise SystemExit(f"prefix placement ratio: h={h}->{h + 1} exceeds {threshold_text}")
        if numerator * max_den > max_num * denominator:
            max_num = numerator
            max_den = denominator
            max_h = h

    n_max = b * (late_blocks + gap_max - 1)
    endpoint_num = (n_max - h_min + 1) * (h_min + 1)
    endpoint_den = h_min * (n - h_min)
    return PlacementRatioSummary(
        peak_h=max_h,
        peak_ratio=ratio_float(max_num, max_den),
        endpoint_ratio_at_h_min=ratio_float(endpoint_num, endpoint_den),
        endpoint_slack_factor=ratio_float(endpoint_num * max_den, endpoint_den * max_num),
    )


def read_interior_audit_summary(path: Path) -> InteriorAuditSummary:
    rows: list[dict[str, str]]
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f"{path}: no interior-audit rows")
    worst = max(rows, key=lambda row: float(row["max_log2_minus_left"]))
    worst_repro = max(rows, key=lambda row: abs(float(row["left_reproduction_delta"])))
    return InteriorAuditSummary(rows=len(rows), worst=worst, worst_repro=worst_repro)


def check_interior_audit(
    path: Path,
    *,
    expected_rows: int,
    tolerance: float,
    reproduction_tolerance: float,
) -> InteriorAuditSummary:
    summary = read_interior_audit_summary(path)
    if summary.rows != expected_rows:
        raise SystemExit(f"prefix interior audit: expected {expected_rows} rows, got {summary.rows}")
    worst_diff = float(summary.worst["max_log2_minus_left"])
    if worst_diff > tolerance:
        raise SystemExit(
            f"prefix interior audit: worst interior increase {worst_diff:.12g} "
            f"exceeds tolerance {tolerance:.12g}"
        )
    worst_repro = abs(float(summary.worst_repro["left_reproduction_delta"]))
    if worst_repro > reproduction_tolerance:
        raise SystemExit(
            f"prefix interior audit: worst left-endpoint reproduction error "
            f"{worst_repro:.12g} exceeds tolerance {reproduction_tolerance:.12g}"
        )
    return summary


def check_close(name: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise SystemExit(
            f"{name}: expected {expected:.6f}, got {actual:.6f}; "
            f"delta={actual - expected:.6g}"
        )


def print_peak(prefix: str, peak: dict[str, str] | None, column: str) -> None:
    if peak is None:
        return
    keys = [
        "outer_weight",
        "first_r",
        "remaining_ones",
        "gap_min",
        "gap_max",
        "inner_mode",
        column,
    ]
    fields = [f"{key}={peak[key]}" for key in keys if key in peak]
    print(f"{prefix}_peak," + ",".join(fields))


def parse_piecewise_total(stdout: str) -> float:
    for line in stdout.splitlines():
        if line.startswith("total_log2,"):
            return float(line.split(",", 1)[1])
    raise ValueError("piecewise certificate output did not contain total_log2")


def selected_high_intervals(selection: str) -> list[HighInterval]:
    if selection.lower() == "all":
        return list(HIGH_INTERVALS_RAW)
    labels = {part.strip() for part in selection.split(",") if part.strip()}
    by_label = {interval.label: interval for interval in HIGH_INTERVALS_RAW}
    missing = sorted(labels - set(by_label))
    if missing:
        raise SystemExit(f"unknown high interval label(s): {', '.join(missing)}")
    return [interval for interval in HIGH_INTERVALS_RAW if interval.label in labels]


def selected_complement_high_intervals(selection: str) -> list[ComplementHighInterval]:
    if selection.lower() == "all":
        return list(COMPLEMENT_HIGH_INTERVALS)
    labels = {part.strip() for part in selection.split(",") if part.strip()}
    by_label = {interval.label: interval for interval in COMPLEMENT_HIGH_INTERVALS}
    missing = sorted(labels - set(by_label))
    if missing:
        raise SystemExit(f"unknown complement high interval label(s): {', '.join(missing)}")
    return [interval for interval in COMPLEMENT_HIGH_INTERVALS if interval.label in labels]


def check_high_interval(interval: HighInterval, tolerance: float) -> float:
    result = subprocess.run(
        interval.command(),
        cwd=ROOT.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    actual = parse_piecewise_total(result.stdout)
    check_close(f"interval_{interval.label}_recompute", actual, interval.adjusted_log2, tolerance)
    return actual


def check_complement_high_interval(interval: ComplementHighInterval, tolerance: float) -> float:
    result = subprocess.run(
        interval.command(),
        cwd=ROOT.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    actual = parse_piecewise_total(result.stdout)
    check_close(f"complement_interval_{interval.label}_recompute", actual, interval.log2_value, tolerance)
    return actual


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prefix-e-le8-csv",
        type=Path,
        default=ROOT / "fullsplit_piecewise_h32_500_csv.csv",
    )
    parser.add_argument(
        "--prefix-e-ge9-tail-csv",
        type=Path,
        default=ROOT / "fullsplit_turnoff_tail_h32_500_r1_64_emin9.csv",
    )
    parser.add_argument(
        "--early-prefix-e-le16-csv",
        type=Path,
        default=ROOT / "fullsplit_piecewise_early_h32_500_e16_uniformsurv.csv",
    )
    parser.add_argument(
        "--early-prefix-e-ge17-tail-csv",
        type=Path,
        default=ROOT / "fullsplit_turnoff_tail_early_h32_500_r1_64_emin17.csv",
    )
    parser.add_argument(
        "--postprefix-eall-hsummary-csv",
        type=Path,
        default=ROOT / "fullsplit_piecewise_h501_2000_eall_hsummary.csv",
    )
    parser.add_argument(
        "--prefix-interior-audit-csv",
        type=Path,
        default=ROOT / "fullsplit_exact_interior_h0_499_allT.csv",
    )
    parser.add_argument(
        "--outer-prefix-csv",
        type=Path,
        default=ROOT / "block_outer_probe_k1048576_rm512_256_sig32_d009_h500_exact.csv",
    )
    parser.add_argument("--exact-outer-local-spectrum-csv", type=Path, default=ROOT / "rm512_256_spectrum.csv")
    parser.add_argument("--exact-outer-blocks", type=int, default=4096)
    parser.add_argument("--exact-outer-h-max", type=int, default=500)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--late-prefix-h-max", type=int, default=500)
    parser.add_argument("--skip-prefix-interior-audit", action="store_true")
    parser.add_argument("--prefix-interior-rows", type=int, default=1500)
    parser.add_argument("--prefix-interior-tolerance", type=float, default=1e-7)
    parser.add_argument("--prefix-interior-reproduction-tolerance", type=float, default=5e-6)
    parser.add_argument("--prefix-ridge-first-ratio-max", type=float, default=0.852700)
    parser.add_argument("--prefix-ridge-tail-ratio-max", type=float, default=0.833795)
    parser.add_argument("--prefix-ridge-total-bound-max-log2", type=float, default=-34.76713)
    parser.add_argument("--prefix-ridge-inner-span-tolerance", type=float, default=1e-9)
    parser.add_argument("--prefix-ridge-first-outer-ratio-max", type=float, default=2.809)
    parser.add_argument("--prefix-ridge-tail-outer-ratio-max", type=float, default=2.747)
    parser.add_argument("--prefix-ridge-placement-ratio-max", type=float, default=0.303594)
    parser.add_argument("--prefix-ridge-inner-ratio-max", type=float, default=1.000000001)
    parser.add_argument("--skip-prefix-placement-ratio-cert", action="store_true")
    parser.add_argument("--prefix-placement-ratio-threshold", default="0.303594")
    parser.add_argument("--prefix-placement-gap-min", type=int, default=1)
    parser.add_argument("--prefix-placement-gap-max", type=int, default=4000)
    parser.add_argument("--prefix-placement-h-min", type=int, default=32)
    parser.add_argument("--prefix-placement-h-max", type=int, default=500)
    parser.add_argument(
        "--print-high-interval-commands",
        action="store_true",
        help="Print the theorem-facing recomputation command for every high interval.",
    )
    parser.add_argument(
        "--check-high-intervals",
        help=(
            "Comma-separated high-interval labels to recompute, or 'all'. "
            "Example: 2001--7858,550001--650000. This is opt-in because it is slower."
        ),
    )
    parser.add_argument(
        "--check-complement-high-intervals",
        help=(
            "Comma-separated complement-high interval labels to recompute, or 'all'. "
            "Example: 1048577--1148736,1500001--1700000."
        ),
    )
    parser.add_argument("--tolerance", type=float, default=5e-6)
    args = parser.parse_args()

    if args.print_high_interval_commands:
        print("High-interval recomputation commands")
        for interval in HIGH_INTERVALS_RAW:
            print(f"interval_{interval.label}_command,{interval.display_command()}")
        print("Complement-high interval recomputation commands")
        for interval in COMPLEMENT_HIGH_INTERVALS:
            print(f"complement_interval_{interval.label}_command,{interval.display_command()}")
        if not args.check_high_intervals and not args.check_complement_high_intervals:
            return 0

    csv_ledgers = [
        CsvLedger(
            "prefix_32_500_e_le8",
            args.prefix_e_le8_csv,
            "term_log2",
            -34.767174,
        ),
        CsvLedger(
            "prefix_32_500_e_ge9_tail",
            args.prefix_e_ge9_tail_csv,
            "term_log2",
            -269.335258,
        ),
        CsvLedger(
            "early_32_500_e_le16_uniformsurv",
            args.early_prefix_e_le16_csv,
            "term_log2",
            -85.403338,
        ),
        CsvLedger(
            "early_32_500_e_ge17_tail",
            args.early_prefix_e_ge17_tail_csv,
            "term_log2",
            -305.738816,
        ),
        CsvLedger(
            "postprefix_501_2000_eall",
            args.postprefix_eall_hsummary_csv,
            "aggregate_log2",
            -182.259739,
        ),
    ]

    print("Full-split finite ledger")
    if not args.skip_prefix_interior_audit:
        interior = check_interior_audit(
            args.prefix_interior_audit_csv,
            expected_rows=args.prefix_interior_rows,
            tolerance=args.prefix_interior_tolerance,
            reproduction_tolerance=args.prefix_interior_reproduction_tolerance,
        )
        print(f"prefix_interior_audit_rows,{interior.rows}")
        print(
            "prefix_interior_audit_worst,"
            f"bucket={interior.worst['bucket']},"
            f"H={interior.worst['H']},"
            f"diff={interior.worst['max_log2_minus_left']},"
            f"T={interior.worst['T_at_max']}"
        )
        print(
            "prefix_interior_audit_worst_repro,"
            f"bucket={interior.worst_repro['bucket']},"
            f"H={interior.worst_repro['H']},"
            f"delta={interior.worst_repro['left_reproduction_delta']}"
        )

    late_prefix, late_rows, late_peak = late_prefix_cap_log2(
        args.outer_prefix_csv,
        n=args.N,
        b=args.block_bits,
        late_blocks=args.late_blocks,
        h_max=args.late_prefix_h_max,
    )
    check_close("late_prefix_T_lt_5949_cap", late_prefix, -40.115431, args.tolerance)
    print(f"late_prefix_T_lt_{args.late_blocks}_cap_rows,{late_rows}")
    print(f"late_prefix_T_lt_{args.late_blocks}_cap_log2,{late_prefix:.6f}")
    if late_peak is not None:
        print(
            f"late_prefix_T_lt_{args.late_blocks}_cap_peak,"
            f"h={int(late_peak['h'])},"
            f"outer_log2={late_peak['outer_log2']:.6f},"
            f"late_log2={late_peak['late_log2']:.6f},"
            f"term_log2={late_peak['term_log2']:.6f}"
        )

    exact_outer_coeffs = direct_sum_coefficients(
        load_local_spectrum(args.exact_outer_local_spectrum_csv, args.exact_outer_h_max),
        blocks=args.exact_outer_blocks,
        h_max=args.exact_outer_h_max,
    )
    exact_late_prefix = exact_late_prefix_sum(
        exact_outer_coeffs,
        n=args.N,
        b=args.block_bits,
        late_blocks=args.late_blocks,
        h_max=args.late_prefix_h_max,
    )
    check_close("late_prefix_T_lt_5949_exact_outer", exact_late_prefix.total_log2, -41.113442, args.tolerance)
    print(f"late_prefix_T_lt_{args.late_blocks}_exact_outer_rows,{exact_late_prefix.rows}")
    print(f"late_prefix_T_lt_{args.late_blocks}_exact_outer_log2,{exact_late_prefix.total_log2:.6f}")
    print(
        f"late_prefix_T_lt_{args.late_blocks}_exact_outer_peak,"
        f"h={exact_late_prefix.peak_h},"
        f"outer_log2={exact_late_prefix.peak_outer_log2:.6f},"
        f"late_log2={exact_late_prefix.peak_late_log2:.6f},"
        f"term_log2={exact_late_prefix.peak_term_log2:.6f}"
    )

    parts: dict[str, float] = {}
    for item in csv_ledgers:
        total, rows, peak = csv_logsum(item.path, item.column)
        check_close(item.name, total, item.expected_log2, args.tolerance)
        parts[item.name] = total
        print(f"{item.name}_rows,{rows}")
        print(f"{item.name}_log2,{total:.6f}")
        print_peak(item.name, peak, item.column)

    exact_outer_expected = {
        "prefix_32_500_e_le8": -37.383345,
        "prefix_32_500_e_ge9_tail": -284.004805,
        "early_32_500_e_le16_uniformsurv": -86.910456,
        "early_32_500_e_ge17_tail": -319.977808,
    }
    exact_outer_parts: dict[str, float] = {}
    for item in csv_ledgers:
        if item.name not in exact_outer_expected:
            continue
        summary = reweighted_ledger_sum(item.path, exact_outer_coeffs)
        check_close(f"{item.name}_exact_outer", summary.exact_log2, exact_outer_expected[item.name], args.tolerance)
        exact_outer_parts[item.name] = summary.exact_log2
        print(f"{item.name}_exact_outer_live_rows,{summary.live_rows}")
        print(f"{item.name}_exact_outer_log2,{summary.exact_log2:.6f}")
        print(
            f"{item.name}_exact_outer_peak,"
            f"outer_weight={summary.peak_h},first_r={summary.peak_first_r},"
            f"gap_min={summary.peak_gap_min},gap_max={summary.peak_gap_max},"
            f"outer_log2={summary.peak_outer_exact_log2:.6f},"
            f"term_log2={summary.peak_term_log2:.6f}"
        )

    ridge_ratio = prefix_ridge_ratio_summary(args.prefix_e_le8_csv)
    if ridge_ratio.first_ratio > args.prefix_ridge_first_ratio_max:
        raise SystemExit(
            f"prefix ridge: first ratio {ridge_ratio.first_ratio:.12g} exceeds "
            f"{args.prefix_ridge_first_ratio_max:.12g}"
        )
    if ridge_ratio.tail_ratio > args.prefix_ridge_tail_ratio_max:
        raise SystemExit(
            f"prefix ridge: tail ratio {ridge_ratio.tail_ratio:.12g} exceeds "
            f"{args.prefix_ridge_tail_ratio_max:.12g}"
        )
    if ridge_ratio.total_bound_log2 > args.prefix_ridge_total_bound_max_log2:
        raise SystemExit(
            f"prefix ridge: ratio bound {ridge_ratio.total_bound_log2:.9f} exceeds "
            f"{args.prefix_ridge_total_bound_max_log2:.9f}"
        )
    if ridge_ratio.inner_span_bits > args.prefix_ridge_inner_span_tolerance:
        raise SystemExit(
            f"prefix ridge: inner span {ridge_ratio.inner_span_bits:.12g} exceeds "
            f"{args.prefix_ridge_inner_span_tolerance:.12g}"
        )
    if ridge_ratio.first_outer_ratio > args.prefix_ridge_first_outer_ratio_max:
        raise SystemExit(
            f"prefix ridge: first outer ratio {ridge_ratio.first_outer_ratio:.12g} exceeds "
            f"{args.prefix_ridge_first_outer_ratio_max:.12g}"
        )
    if ridge_ratio.tail_outer_ratio > args.prefix_ridge_tail_outer_ratio_max:
        raise SystemExit(
            f"prefix ridge: tail outer ratio {ridge_ratio.tail_outer_ratio:.12g} exceeds "
            f"{args.prefix_ridge_tail_outer_ratio_max:.12g}"
        )
    if max(ridge_ratio.first_placement_ratio, ridge_ratio.tail_placement_ratio) > args.prefix_ridge_placement_ratio_max:
        raise SystemExit(
            "prefix ridge: placement ratio exceeds "
            f"{args.prefix_ridge_placement_ratio_max:.12g}"
        )
    if max(ridge_ratio.first_inner_ratio, ridge_ratio.tail_inner_ratio) > args.prefix_ridge_inner_ratio_max:
        raise SystemExit(
            "prefix ridge: inner ratio exceeds "
            f"{args.prefix_ridge_inner_ratio_max:.12g}"
        )
    print(f"prefix_ridge_ratio_rows,{ridge_ratio.ridge_rows}")
    print(f"prefix_ridge_exact_log2,{ridge_ratio.ridge_exact_log2:.6f}")
    print(f"prefix_ridge_remainder_exact_log2,{ridge_ratio.remainder_exact_log2:.6f}")
    print(f"prefix_ridge_first_ratio,{ridge_ratio.first_ratio:.12g}")
    print(f"prefix_ridge_tail_ratio,{ridge_ratio.tail_ratio:.12g}")
    print(f"prefix_ridge_first_outer_ratio,{ridge_ratio.first_outer_ratio:.12g}")
    print(f"prefix_ridge_first_placement_ratio,{ridge_ratio.first_placement_ratio:.12g}")
    print(f"prefix_ridge_first_inner_ratio,{ridge_ratio.first_inner_ratio:.12g}")
    print(f"prefix_ridge_tail_outer_ratio,{ridge_ratio.tail_outer_ratio:.12g}")
    print(f"prefix_ridge_tail_placement_ratio,{ridge_ratio.tail_placement_ratio:.12g}")
    print(f"prefix_ridge_tail_inner_ratio,{ridge_ratio.tail_inner_ratio:.12g}")
    print(f"prefix_ridge_geometric_bound_log2,{ridge_ratio.geometric_log2:.6f}")
    print(f"prefix_ridge_total_ratio_bound_log2,{ridge_ratio.total_bound_log2:.6f}")
    print(f"prefix_ridge_total_ratio_bound_slack_bits,{ridge_ratio.total_bound_log2 - ridge_ratio.total_exact_log2:.12g}")
    print(f"prefix_ridge_inner_span_bits,{ridge_ratio.inner_span_bits:.12g}")

    if not args.skip_prefix_placement_ratio_cert:
        placement = exact_placement_ratio_summary(
            n=args.N,
            b=args.block_bits,
            late_blocks=args.late_blocks,
            gap_min=args.prefix_placement_gap_min,
            gap_max=args.prefix_placement_gap_max,
            h_min=args.prefix_placement_h_min,
            h_max=args.prefix_placement_h_max,
            threshold_text=args.prefix_placement_ratio_threshold,
        )
        print(f"prefix_placement_ratio_threshold,{args.prefix_placement_ratio_threshold}")
        print(f"prefix_placement_ratio_peak_h,{placement.peak_h}")
        print(f"prefix_placement_ratio_peak,{placement.peak_ratio:.12g}")
        print(f"prefix_placement_endpoint_bound_at_h_min,{placement.endpoint_ratio_at_h_min:.12g}")
        print(f"prefix_placement_endpoint_slack_factor,{placement.endpoint_slack_factor:.12g}")

    prefix_all = log2add(
        exact_outer_parts["prefix_32_500_e_le8"],
        exact_outer_parts["prefix_32_500_e_ge9_tail"],
    )
    print(f"prefix_32_500_all_e_exact_outer_log2,{prefix_all:.6f}")
    early_all = log2add(
        exact_outer_parts["early_32_500_e_le16_uniformsurv"],
        exact_outer_parts["early_32_500_e_ge17_tail"],
    )
    print(f"early_32_500_all_e_exact_outer_log2,{early_all:.6f}")

    high_total = float("-inf")
    for label, lam, rho, value in HIGH_INTERVALS:
        high_total = log2add(high_total, value)
        print(f"interval_{label}_lambda,{lam:g}")
        print(f"interval_{label}_rho,{rho:g}")
        print(f"interval_{label}_log2,{value:.6f}")
    print(f"interval_2001_1148736_total_log2,{high_total:.6f}")
    print("feasible_far_bucket_cutoff_h,1148736")

    complement_high_total = float("-inf")
    for interval in COMPLEMENT_HIGH_INTERVALS:
        complement_high_total = log2add(complement_high_total, interval.log2_value)
        print(f"complement_interval_{interval.label}_rho,{interval.rho:g}")
        print(f"complement_interval_{interval.label}_log2,{interval.log2_value:.6f}")
    print(f"complement_interval_1048577_2097152_total_log2,{complement_high_total:.6f}")

    if args.check_high_intervals:
        for interval in selected_high_intervals(args.check_high_intervals):
            actual = check_high_interval(interval, args.tolerance)
            print(f"interval_{interval.label}_recomputed_log2,{actual:.6f}")
            print(f"interval_{interval.label}_recompute_status,PASS")

    if args.check_complement_high_intervals:
        for interval in selected_complement_high_intervals(args.check_complement_high_intervals):
            actual = check_complement_high_interval(interval, args.tolerance)
            print(f"complement_interval_{interval.label}_recomputed_log2,{actual:.6f}")
            print(f"complement_interval_{interval.label}_recompute_status,PASS")

    total = float("-inf")
    for value in [prefix_all, early_all, parts["postprefix_501_2000_eall"], high_total, complement_high_total]:
        total = log2add(total, value)
    print(f"finite_ledger_total_log2,{total:.6f}")
    print(f"finite_ledger_margin_bits,{-total:.6f}")
    late_plus_window = log2add(exact_late_prefix.total_log2, total)
    print(f"late_plus_window_total_log2,{late_plus_window:.6f}")
    print(f"late_plus_window_margin_bits,{-late_plus_window:.6f}")
    h32_500_all_positions = log2add(exact_late_prefix.total_log2, log2add(prefix_all, early_all))
    print(f"h32_500_all_first_active_positions_log2,{h32_500_all_positions:.6f}")
    print(f"h32_500_all_first_active_positions_margin_bits,{-h32_500_all_positions:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
