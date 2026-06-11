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
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--late-prefix-h-max", type=int, default=500)
    parser.add_argument("--skip-prefix-interior-audit", action="store_true")
    parser.add_argument("--prefix-interior-rows", type=int, default=1500)
    parser.add_argument("--prefix-interior-tolerance", type=float, default=1e-7)
    parser.add_argument("--prefix-interior-reproduction-tolerance", type=float, default=5e-6)
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
    parser.add_argument("--tolerance", type=float, default=5e-6)
    args = parser.parse_args()

    if args.print_high_interval_commands:
        print("High-interval recomputation commands")
        for interval in HIGH_INTERVALS_RAW:
            print(f"interval_{interval.label}_command,{interval.display_command()}")
        if not args.check_high_intervals:
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

    parts: dict[str, float] = {}
    for item in csv_ledgers:
        total, rows, peak = csv_logsum(item.path, item.column)
        check_close(item.name, total, item.expected_log2, args.tolerance)
        parts[item.name] = total
        print(f"{item.name}_rows,{rows}")
        print(f"{item.name}_log2,{total:.6f}")
        print_peak(item.name, peak, item.column)

    prefix_all = log2add(parts["prefix_32_500_e_le8"], parts["prefix_32_500_e_ge9_tail"])
    print(f"prefix_32_500_all_e_log2,{prefix_all:.6f}")
    early_all = log2add(
        parts["early_32_500_e_le16_uniformsurv"],
        parts["early_32_500_e_ge17_tail"],
    )
    print(f"early_32_500_all_e_log2,{early_all:.6f}")

    high_total = float("-inf")
    for label, lam, rho, value in HIGH_INTERVALS:
        high_total = log2add(high_total, value)
        print(f"interval_{label}_lambda,{lam:g}")
        print(f"interval_{label}_rho,{rho:g}")
        print(f"interval_{label}_log2,{value:.6f}")
    print(f"interval_2001_1148736_total_log2,{high_total:.6f}")
    print("feasible_far_bucket_cutoff_h,1148736")

    if args.check_high_intervals:
        for interval in selected_high_intervals(args.check_high_intervals):
            actual = check_high_interval(interval, args.tolerance)
            print(f"interval_{interval.label}_recomputed_log2,{actual:.6f}")
            print(f"interval_{interval.label}_recompute_status,PASS")

    total = float("-inf")
    for value in [prefix_all, early_all, parts["postprefix_501_2000_eall"], high_total]:
        total = log2add(total, value)
    print(f"finite_ledger_total_log2,{total:.6f}")
    print(f"finite_ledger_margin_bits,{-total:.6f}")
    late_plus_window = log2add(late_prefix, total)
    print(f"late_plus_window_total_log2,{late_plus_window:.6f}")
    print(f"late_plus_window_margin_bits,{-late_plus_window:.6f}")
    h32_500_all_positions = log2add(late_prefix, log2add(prefix_all, early_all))
    print(f"h32_500_all_first_active_positions_log2,{h32_500_all_positions:.6f}")
    print(f"h32_500_all_first_active_positions_margin_bits,{-h32_500_all_positions:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
