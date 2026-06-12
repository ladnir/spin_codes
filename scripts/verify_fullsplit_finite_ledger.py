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
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from bound_fullsplit_episode_gaps import load_spectrum as load_inner_spectrum
from bound_fullsplit_episode_gaps import split_entries as split_inner_entries
from check_fullsplit_T_monotonicity import (
    DEFAULT_INTERVALS as T_MONOTONICITY_INTERVALS,
    sufficient_reduction_rows,
)
from certify_prefix_placement_ratio import certify_placement_ratio
from certify_rm_outer_prefix_exact import (
    exact_late_prefix_sum,
    load_local_spectrum,
    prefix_coefficient_certificate,
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


@dataclass(frozen=True)
class EarlyPostprefixInterval:
    start: int
    stop: int
    log2_value: float
    lambdas: str | None = None
    rhos: str | None = None

    @property
    def label(self) -> str:
        return f"{self.start}--{self.stop}"

    def command(self) -> list[str]:
        cmd = [
            sys.executable,
            str(ROOT / "sum_fullsplit_piecewise_certificate.py"),
            "--h-values",
            f"{self.start}:{self.stop}",
            "--gap-start",
            "12001",
            "--gap-stop",
            "26819",
            "--gap-step",
            "5000",
            "--inner-mode-by-gap",
            "eallratio,eallratio,eallratio",
            "--gap-sum-mode",
            "endpoint",
            "--inner-T-by-gap",
            "feasible-min",
            "--turnoff-log2",
            f"{GLOBAL_TURNOFF_LOG2:.7f}",
        ]
        if self.lambdas is not None:
            cmd.extend(["--lambdas", self.lambdas])
        if self.rhos is not None:
            cmd.extend(["--rhos", self.rhos])
        return cmd

    def display_command(self) -> str:
        script = r"scripts\sum_fullsplit_piecewise_certificate.py"
        cmd = (
            f"python {script} --h-values {self.start}:{self.stop} "
            "--gap-start 12001 --gap-stop 26819 --gap-step 5000 "
            "--inner-mode-by-gap eallratio,eallratio,eallratio "
            "--gap-sum-mode endpoint --inner-T-by-gap feasible-min "
            f"--turnoff-log2 {GLOBAL_TURNOFF_LOG2:.7f}"
        )
        if self.lambdas is not None:
            cmd += f" --lambdas {self.lambdas}"
        if self.rhos is not None:
            cmd += f" --rhos {self.rhos}"
        return cmd


@dataclass(frozen=True)
class EarlyAcceleratedInterval:
    start: int
    stop: int
    lam: float
    rho: float
    log2_value: float
    script_name: str

    @property
    def label(self) -> str:
        return f"{self.start}--{self.stop}"

    @property
    def kind(self) -> str:
        if "paired" in self.script_name:
            return "paired"
        return "endpoint"

    def command(self) -> list[str]:
        return [
            sys.executable,
            str(ROOT / self.script_name),
            "--intervals",
            f"{self.start}:{self.stop}:{self.lam:g}:{self.rho:g}",
        ]

    def display_command(self) -> str:
        script = rf"scripts\{self.script_name}"
        return f"python {script} --intervals {self.start}:{self.stop}:{self.lam:g}:{self.rho:g}"


@dataclass(frozen=True)
class EarlyAcceleratedCoverSummary:
    count: int
    cover_start: int
    cover_stop: int
    endpoint_start: int
    endpoint_stop: int
    endpoint_count: int
    paired_start: int
    paired_stop: int
    paired_count: int


@dataclass(frozen=True)
class IntervalCoverSummary:
    count: int
    cover_start: int
    cover_stop: int


@dataclass(frozen=True)
class HighFeasibleMinSummary:
    bucket_cutoffs: list[tuple[str, int, int, int]]
    far_cutoff_h: int


@dataclass(frozen=True)
class LatePostprefixShapeSummary:
    interval_count: int
    peak_h_min: int
    peak_h_max: int
    left_endpoint_peaks: int
    right_endpoint_peaks: int
    critical_peaks: int
    rows: list[dict[str, object]]


@dataclass(frozen=True)
class ComplementHighShapeSummary:
    interval_count: int
    min_convexity_growth_bits: float
    rows: list[dict[str, object]]


@dataclass(frozen=True)
class TMonotonicitySummary:
    interval_count: int
    bucket_count: int
    row_count: int
    min_slack_bits: float
    worst_row: dict[str, object]
    rows: list[dict[str, object]]


@dataclass(frozen=True)
class LatePostprefixInterval:
    start: int
    stop: int
    z: float
    log2_value: float

    @property
    def label(self) -> str:
        return f"{self.start}--{self.stop}"

    def command(self) -> list[str]:
        return [
            sys.executable,
            str(ROOT / "certify_fullsplit_late_prefix_interval.py"),
            "--intervals",
            f"{self.start}:{self.stop}:{self.z:.17g}",
        ]

    def display_command(self) -> str:
        script = r"scripts\certify_fullsplit_late_prefix_interval.py"
        return f"python {script} --intervals {self.start}:{self.stop}:{self.z:.17g}"


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


EARLY_POSTPREFIX_INTERVALS = [
    EarlyPostprefixInterval(501, 2000, -380.829185),
    EarlyPostprefixInterval(2001, 7858, -852.998772, "0.02", "0.01"),
    EarlyPostprefixInterval(7859, 20550, -3712.365087, "0.05", "0.03"),
    EarlyPostprefixInterval(20551, 30000, -8614.382763, "0.2", "0.1"),
    EarlyPostprefixInterval(30001, 50000, -21332.131033, "0.2", "0.1"),
    EarlyPostprefixInterval(50001, 75000, -34847.929043, "0.2", "0.1"),
]


EARLY_ACCELERATED_INTERVALS = [
    EarlyAcceleratedInterval(75001, 90000, 0.2, 0.1, -28574.869802, "certify_fullsplit_early_endpoint_interval.py"),
    EarlyAcceleratedInterval(90001, 100000, 0.2, 0.1, -22149.918612, "certify_fullsplit_early_endpoint_interval.py"),
    EarlyAcceleratedInterval(100001, 110000, 0.5, 0.3, -1056.187188, "certify_fullsplit_early_endpoint_interval.py"),
    EarlyAcceleratedInterval(110001, 125000, 0.5, 0.3, -10362.025949, "certify_fullsplit_early_endpoint_interval.py"),
    EarlyAcceleratedInterval(125001, 160000, 0.8, 0.3, -45589.980049, "certify_fullsplit_early_endpoint_interval.py"),
    EarlyAcceleratedInterval(160001, 250000, 1.2, 0.5, -60293.751685, "certify_fullsplit_early_endpoint_interval.py"),
    EarlyAcceleratedInterval(250001, 350000, 1.2, 0.5, -94682.265141, "certify_fullsplit_early_endpoint_interval.py"),
    EarlyAcceleratedInterval(350001, 524288, 2.0, 0.8, -181215.853042, "certify_fullsplit_early_paired_interval.py"),
    EarlyAcceleratedInterval(524289, 750000, 2.0, 0.8, -194565.849984, "certify_fullsplit_early_paired_interval.py"),
    EarlyAcceleratedInterval(750001, 1048576, 2.0, 0.8, -98386.449256, "certify_fullsplit_early_paired_interval.py"),
]


LATE_POSTPREFIX_Z = 0.39605985943459426
LATE_POSTPREFIX_INTERVALS = [
    LatePostprefixInterval(501, 2000, LATE_POSTPREFIX_Z, -545.690526),
    LatePostprefixInterval(2001, 7858, LATE_POSTPREFIX_Z, -2237.593414),
    LatePostprefixInterval(7859, 20550, LATE_POSTPREFIX_Z, -8919.160723),
    LatePostprefixInterval(20551, 75000, LATE_POSTPREFIX_Z, -23772.761808),
    LatePostprefixInterval(75001, 150000, LATE_POSTPREFIX_Z, -93856.032305),
    LatePostprefixInterval(150001, 250000, LATE_POSTPREFIX_Z, -210536.528296),
    LatePostprefixInterval(250001, 380736, LATE_POSTPREFIX_Z, -417970.896784),
]


@dataclass(frozen=True)
class CsvLedger:
    name: str
    path: Path
    column: str
    expected_log2: float


@dataclass(frozen=True)
class InteriorAuditBucket:
    name: str
    t_min: int
    t_max: int

    @property
    def t_count(self) -> int:
        return self.t_max - self.t_min + 1

    @property
    def label(self) -> str:
        return f"{self.name}:{self.t_min}--{self.t_max}:{self.t_count}"


INTERIOR_AUDIT_BUCKETS = (
    InteriorAuditBucket("gap_1_4000", 5949, 9948),
    InteriorAuditBucket("gap_4001_8000", 9949, 13948),
    InteriorAuditBucket("gap_8001_12000", 13949, 17948),
)
INTERIOR_AUDIT_BUCKET_BY_NAME = {bucket.name: bucket for bucket in INTERIOR_AUDIT_BUCKETS}


@dataclass(frozen=True)
class InteriorAuditSummary:
    rows: int
    h_min: int
    h_max: int
    h_count: int
    bucket_counts: dict[str, int]
    bucket_t_counts: dict[str, int]
    left_endpoint_maxima: bool
    worst_diff: float
    worst_diff_slack: float
    worst_reproduction_abs: float
    worst_reproduction_slack: float
    worst: dict[str, str]
    worst_repro: dict[str, str]


@dataclass(frozen=True)
class PlacementRatioSummary:
    peak_h: int
    peak_to_h: int
    peak_ratio: float
    peak_log2: float
    threshold_num: int
    threshold_den: int
    exact_comparisons: int
    peak_num_bits: int
    peak_den_bits: int
    threshold_slack_bits: int
    threshold_gap_log2: float
    endpoint_ratio_at_h_min: float
    endpoint_ratio_at_h_min_log2: float
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


def finite_float(value: float) -> float | None:
    if math.isinf(value) or math.isnan(value):
        return None
    return value


def interval_range(start: int, stop: int) -> list[int]:
    return [start, stop]


def manifest_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.parent))
    except ValueError:
        return str(path)


def add_row(
    manifest: dict[str, object],
    *,
    name: str,
    log2_value: float,
    role: str,
    h_range: tuple[int, int] | None = None,
    rows: int | None = None,
    command: str | None = None,
    notes: str | None = None,
    **extra: object,
) -> None:
    row: dict[str, object] = {
        "name": name,
        "role": role,
        "log2": finite_float(log2_value),
    }
    if h_range is not None:
        row["h_range"] = interval_range(*h_range)
    if rows is not None:
        row["rows"] = rows
    if command is not None:
        row["command"] = command
    if notes is not None:
        row["notes"] = notes
    row.update(extra)
    manifest["row_families"].append(row)  # type: ignore[index,union-attr]


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
    cert = certify_placement_ratio(
        N=n,
        b=b,
        late_blocks=late_blocks,
        gap_min=gap_min,
        gap_max=gap_max,
        h_min=h_min,
        h_max=h_max,
        threshold_text=threshold_text,
    )
    return PlacementRatioSummary(
        peak_h=cert.peak_h,
        peak_to_h=cert.peak_to_h,
        peak_ratio=cert.peak_ratio,
        peak_log2=cert.peak_log2,
        threshold_num=cert.threshold.numerator,
        threshold_den=cert.threshold.denominator,
        exact_comparisons=cert.comparisons,
        peak_num_bits=cert.peak_num.bit_length(),
        peak_den_bits=cert.peak_den.bit_length(),
        threshold_slack_bits=cert.threshold_slack.bit_length(),
        threshold_gap_log2=cert.threshold_gap_log2,
        endpoint_ratio_at_h_min=cert.endpoint_ratio,
        endpoint_ratio_at_h_min_log2=cert.endpoint_log2,
        endpoint_slack_factor=cert.endpoint_slack_factor,
    )


def read_interior_audit_summary(
    path: Path,
    *,
    expected_h_min: int,
    expected_h_max: int,
    tolerance: float,
    reproduction_tolerance: float,
) -> InteriorAuditSummary:
    rows: list[dict[str, str]]
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f"{path}: no interior-audit rows")

    required_columns = {
        "bucket",
        "H",
        "T_count",
        "max_log2_minus_left",
        "T_at_max",
        "left_log2",
        "max_log2",
        "left_reproduction_delta",
    }
    missing_columns = required_columns - set(rows[0])
    if missing_columns:
        raise SystemExit(f"{path}: missing interior-audit columns {sorted(missing_columns)}")

    expected_h_values = set(range(expected_h_min, expected_h_max + 1))
    bucket_counts = {bucket.name: 0 for bucket in INTERIOR_AUDIT_BUCKETS}
    bucket_t_counts: dict[str, int] = {}
    seen: set[tuple[str, int]] = set()
    left_endpoint_maxima = True

    for row in rows:
        bucket_name = row["bucket"]
        bucket = INTERIOR_AUDIT_BUCKET_BY_NAME.get(bucket_name)
        if bucket is None:
            raise SystemExit(f"{path}: unexpected interior-audit bucket {bucket_name!r}")
        try:
            h = int(row["H"])
            t_count = int(row["T_count"])
            t_at_max = int(row["T_at_max"])
            float(row["left_log2"])
            float(row["max_log2"])
            float(row["max_log2_minus_left"])
            float(row["left_reproduction_delta"])
        except ValueError as exc:
            raise SystemExit(f"{path}: malformed interior-audit row {row}") from exc
        if h not in expected_h_values:
            raise SystemExit(
                f"{path}: interior-audit H={h} is outside "
                f"{expected_h_min}--{expected_h_max}"
            )
        key = (bucket_name, h)
        if key in seen:
            raise SystemExit(f"{path}: duplicate interior-audit row for bucket={bucket_name}, H={h}")
        seen.add(key)
        bucket_counts[bucket_name] += 1
        if t_count != bucket.t_count:
            raise SystemExit(
                f"{path}: bucket={bucket_name}, H={h} has T_count={t_count}, "
                f"expected {bucket.t_count}"
            )
        bucket_t_counts[bucket_name] = t_count
        if not (bucket.t_min <= t_at_max <= bucket.t_max):
            raise SystemExit(
                f"{path}: bucket={bucket_name}, H={h} has T_at_max={t_at_max}, "
                f"outside {bucket.t_min}--{bucket.t_max}"
            )
        if t_at_max != bucket.t_min:
            left_endpoint_maxima = False

    missing_pairs = [
        f"{bucket.name}:H={h}"
        for bucket in INTERIOR_AUDIT_BUCKETS
        for h in range(expected_h_min, expected_h_max + 1)
        if (bucket.name, h) not in seen
    ]
    if missing_pairs:
        preview = ", ".join(missing_pairs[:5])
        suffix = "" if len(missing_pairs) <= 5 else f", ... ({len(missing_pairs)} total)"
        raise SystemExit(f"{path}: missing interior-audit rows {preview}{suffix}")

    worst = max(rows, key=lambda row: float(row["max_log2_minus_left"]))
    worst_repro = max(rows, key=lambda row: abs(float(row["left_reproduction_delta"])))
    worst_diff = float(worst["max_log2_minus_left"])
    worst_reproduction_abs = abs(float(worst_repro["left_reproduction_delta"]))
    return InteriorAuditSummary(
        rows=len(rows),
        h_min=expected_h_min,
        h_max=expected_h_max,
        h_count=expected_h_max - expected_h_min + 1,
        bucket_counts=bucket_counts,
        bucket_t_counts=bucket_t_counts,
        left_endpoint_maxima=left_endpoint_maxima,
        worst_diff=worst_diff,
        worst_diff_slack=tolerance - worst_diff,
        worst_reproduction_abs=worst_reproduction_abs,
        worst_reproduction_slack=reproduction_tolerance - worst_reproduction_abs,
        worst=worst,
        worst_repro=worst_repro,
    )


def check_interior_audit(
    path: Path,
    *,
    expected_rows: int,
    expected_h_min: int,
    expected_h_max: int,
    tolerance: float,
    reproduction_tolerance: float,
) -> InteriorAuditSummary:
    summary = read_interior_audit_summary(
        path,
        expected_h_min=expected_h_min,
        expected_h_max=expected_h_max,
        tolerance=tolerance,
        reproduction_tolerance=reproduction_tolerance,
    )
    if summary.rows != expected_rows:
        raise SystemExit(f"prefix interior audit: expected {expected_rows} rows, got {summary.rows}")
    expected_rows_from_grid = len(INTERIOR_AUDIT_BUCKETS) * summary.h_count
    if summary.rows != expected_rows_from_grid:
        raise SystemExit(
            f"prefix interior audit: expected {expected_rows_from_grid} rows from "
            f"{summary.h_count} H-values and {len(INTERIOR_AUDIT_BUCKETS)} buckets, got {summary.rows}"
        )
    if not summary.left_endpoint_maxima:
        raise SystemExit("prefix interior audit: an interior row is maximized away from the left endpoint")
    if summary.worst_diff > tolerance:
        raise SystemExit(
            f"prefix interior audit: worst interior increase {summary.worst_diff:.12g} "
            f"exceeds tolerance {tolerance:.12g}"
        )
    if summary.worst_reproduction_abs > reproduction_tolerance:
        raise SystemExit(
            f"prefix interior audit: worst left-endpoint reproduction error "
            f"{summary.worst_reproduction_abs:.12g} exceeds tolerance {reproduction_tolerance:.12g}"
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


def selected_early_postprefix_intervals(selection: str) -> list[EarlyPostprefixInterval]:
    if selection.lower() == "all":
        return list(EARLY_POSTPREFIX_INTERVALS)
    labels = {part.strip() for part in selection.split(",") if part.strip()}
    by_label = {interval.label: interval for interval in EARLY_POSTPREFIX_INTERVALS}
    missing = sorted(labels - set(by_label))
    if missing:
        raise SystemExit(f"unknown early-postprefix interval label(s): {', '.join(missing)}")
    return [interval for interval in EARLY_POSTPREFIX_INTERVALS if interval.label in labels]


def selected_early_accelerated_intervals(selection: str) -> list[EarlyAcceleratedInterval]:
    if selection.lower() == "all":
        return list(EARLY_ACCELERATED_INTERVALS)
    labels = {part.strip() for part in selection.split(",") if part.strip()}
    by_label = {interval.label: interval for interval in EARLY_ACCELERATED_INTERVALS}
    missing = sorted(labels - set(by_label))
    if missing:
        raise SystemExit(f"unknown early-accelerated interval label(s): {', '.join(missing)}")
    return [interval for interval in EARLY_ACCELERATED_INTERVALS if interval.label in labels]


def check_early_accelerated_cover(
    intervals: list[EarlyAcceleratedInterval],
    *,
    expected_start: int,
    expected_stop: int,
    expected_endpoint_stop: int,
    expected_paired_start: int,
) -> EarlyAcceleratedCoverSummary:
    if not intervals:
        raise SystemExit("early accelerated cover: no intervals configured")
    ordered = sorted(intervals, key=lambda interval: interval.start)
    if ordered[0].start != expected_start or ordered[-1].stop != expected_stop:
        raise SystemExit(
            "early accelerated cover: expected "
            f"{expected_start}--{expected_stop}, got {ordered[0].start}--{ordered[-1].stop}"
        )
    previous_stop = expected_start - 1
    for interval in ordered:
        if interval.start != previous_stop + 1:
            raise SystemExit(
                "early accelerated cover: gap or overlap before "
                f"{interval.start}--{interval.stop}; previous stop was {previous_stop}"
            )
        previous_stop = interval.stop
        expected_helper = (
            "certify_fullsplit_early_endpoint_interval.py"
            if interval.kind == "endpoint"
            else "certify_fullsplit_early_paired_interval.py"
        )
        if interval.script_name != expected_helper:
            raise SystemExit(
                f"early accelerated cover: interval {interval.label} has helper "
                f"{interval.script_name}, expected {expected_helper}"
            )
        if interval.lam <= 0.0 or interval.rho <= 0.0:
            raise SystemExit(f"early accelerated cover: interval {interval.label} has nonpositive pole")

    endpoint = [interval for interval in ordered if interval.kind == "endpoint"]
    paired = [interval for interval in ordered if interval.kind == "paired"]
    if not endpoint or not paired:
        raise SystemExit("early accelerated cover: both endpoint and paired regions are required")
    if endpoint[0].start != expected_start or endpoint[-1].stop != expected_endpoint_stop:
        raise SystemExit(
            "early accelerated endpoint cover: expected "
            f"{expected_start}--{expected_endpoint_stop}, got {endpoint[0].start}--{endpoint[-1].stop}"
        )
    if paired[0].start != expected_paired_start or paired[-1].stop != expected_stop:
        raise SystemExit(
            "early accelerated paired cover: expected "
            f"{expected_paired_start}--{expected_stop}, got {paired[0].start}--{paired[-1].stop}"
        )
    if expected_paired_start != expected_endpoint_stop + 1:
        raise SystemExit("early accelerated cover: endpoint/paired split is not adjacent")
    if any(interval.kind != "endpoint" for interval in ordered[: len(endpoint)]):
        raise SystemExit("early accelerated cover: endpoint intervals are not a prefix")
    if any(interval.kind != "paired" for interval in ordered[len(endpoint) :]):
        raise SystemExit("early accelerated cover: paired intervals are not a suffix")

    return EarlyAcceleratedCoverSummary(
        count=len(ordered),
        cover_start=ordered[0].start,
        cover_stop=ordered[-1].stop,
        endpoint_start=endpoint[0].start,
        endpoint_stop=endpoint[-1].stop,
        endpoint_count=len(endpoint),
        paired_start=paired[0].start,
        paired_stop=paired[-1].stop,
        paired_count=len(paired),
    )


def selected_late_postprefix_intervals(selection: str) -> list[LatePostprefixInterval]:
    if selection.lower() == "all":
        return list(LATE_POSTPREFIX_INTERVALS)
    labels = {part.strip() for part in selection.split(",") if part.strip()}
    by_label = {interval.label: interval for interval in LATE_POSTPREFIX_INTERVALS}
    missing = sorted(labels - set(by_label))
    if missing:
        raise SystemExit(f"unknown late-postprefix interval label(s): {', '.join(missing)}")
    return [interval for interval in LATE_POSTPREFIX_INTERVALS if interval.label in labels]


def check_contiguous_cover(
    name: str,
    intervals,
    *,
    expected_start: int,
    expected_stop: int,
) -> IntervalCoverSummary:
    if not intervals:
        raise SystemExit(f"{name}: no intervals configured")
    ordered = sorted(intervals, key=lambda interval: interval.start)
    if ordered[0].start != expected_start or ordered[-1].stop != expected_stop:
        raise SystemExit(
            f"{name}: expected cover {expected_start}--{expected_stop}, "
            f"got {ordered[0].start}--{ordered[-1].stop}"
        )
    previous_stop = expected_start - 1
    for interval in ordered:
        if interval.start != previous_stop + 1:
            raise SystemExit(
                f"{name}: gap or overlap before {interval.start}--{interval.stop}; "
                f"previous stop was {previous_stop}"
            )
        if interval.stop < interval.start:
            raise SystemExit(f"{name}: malformed interval {interval.start}--{interval.stop}")
        previous_stop = interval.stop
    return IntervalCoverSummary(
        count=len(ordered),
        cover_start=ordered[0].start,
        cover_stop=ordered[-1].stop,
    )


def interval_overlap(
    left: IntervalCoverSummary,
    right: IntervalCoverSummary,
) -> tuple[int, int] | None:
    start = max(left.cover_start, right.cover_start)
    stop = min(left.cover_stop, right.cover_stop)
    if start > stop:
        return None
    return (start, stop)


def check_high_feasible_min_t(
    *,
    b: int,
    high_cover: IntervalCoverSummary,
) -> HighFeasibleMinSummary:
    bucket_cutoffs: list[tuple[str, int, int, int]] = []
    for bucket in INTERIOR_AUDIT_BUCKETS:
        cutoff_h = b * bucket.t_max + b
        bucket_cutoffs.append((bucket.name, bucket.t_min, bucket.t_max, cutoff_h))
    far_cutoff_h = max(cutoff for _, _, _, cutoff in bucket_cutoffs)
    if high_cover.cover_stop != far_cutoff_h:
        raise SystemExit(
            "high feasible-min T audit: expected high cover to stop at "
            f"far cutoff {far_cutoff_h}, got {high_cover.cover_stop}"
        )
    return HighFeasibleMinSummary(bucket_cutoffs=bucket_cutoffs, far_cutoff_h=far_cutoff_h)


def late_postprefix_reduced_log2(*, n: int, m: int, h: int, z: float) -> float:
    return -h * math.log2(z) + log2_binom(m, h) - log2_binom(n, h)


def late_postprefix_adjacent_ratio_log2(*, n: int, m: int, h: int, z: float) -> float:
    if h < 0 or h >= min(m, n):
        raise ValueError(f"late post-prefix adjacent ratio is undefined at h={h}")
    return math.log2(m - h) - math.log2(n - h) - math.log2(z)


def check_late_postprefix_shapes(
    intervals: list[LatePostprefixInterval],
    *,
    n: int,
    m: int,
    tolerance: float = 1e-10,
) -> LatePostprefixShapeSummary:
    rows: list[dict[str, object]] = []
    left_endpoint_peaks = 0
    right_endpoint_peaks = 0
    critical_peaks = 0
    peak_h_values: list[int] = []
    for interval in intervals:
        if not (0.0 < interval.z < 1.0):
            raise SystemExit(f"late post-prefix shape: {interval.label} has z outside (0,1)")
        if interval.start < 0 or interval.stop > min(n, m):
            raise SystemExit(
                f"late post-prefix shape: {interval.label} is outside feasible placement range 0--{min(n, m)}"
            )
        critical = (m - interval.z * n) / (1.0 - interval.z)
        candidates = {interval.start, interval.stop}
        for h in (math.floor(critical), math.ceil(critical)):
            if interval.start <= h <= interval.stop:
                candidates.add(h)
        best_h = max(
            candidates,
            key=lambda h: late_postprefix_reduced_log2(n=n, m=m, h=h, z=interval.z),
        )
        if best_h > interval.start:
            left_ratio = late_postprefix_adjacent_ratio_log2(n=n, m=m, h=best_h - 1, z=interval.z)
            if left_ratio < -tolerance:
                raise SystemExit(
                    f"late post-prefix shape: {interval.label} peak {best_h} is below its left neighbor"
                )
        if best_h < interval.stop:
            right_ratio = late_postprefix_adjacent_ratio_log2(n=n, m=m, h=best_h, z=interval.z)
            if right_ratio > tolerance:
                raise SystemExit(
                    f"late post-prefix shape: {interval.label} peak {best_h} is below its right neighbor"
                )
        if best_h == interval.start:
            peak_source = "left_endpoint"
            left_endpoint_peaks += 1
        elif best_h == interval.stop:
            peak_source = "right_endpoint"
            right_endpoint_peaks += 1
        else:
            peak_source = "critical"
            critical_peaks += 1
        peak_h_values.append(best_h)
        rows.append(
            {
                "range": [interval.start, interval.stop],
                "z": interval.z,
                "critical_h": critical,
                "peak_h": best_h,
                "peak_source": peak_source,
            }
        )
    return LatePostprefixShapeSummary(
        interval_count=len(intervals),
        peak_h_min=min(peak_h_values),
        peak_h_max=max(peak_h_values),
        left_endpoint_peaks=left_endpoint_peaks,
        right_endpoint_peaks=right_endpoint_peaks,
        critical_peaks=critical_peaks,
        rows=rows,
    )


def check_complement_high_shapes(
    intervals: list[ComplementHighInterval],
    *,
    n: int,
) -> ComplementHighShapeSummary:
    rows: list[dict[str, object]] = []
    min_growth = float("inf")
    for interval in intervals:
        if interval.start <= n // 2 or interval.stop > n:
            raise SystemExit(f"complement high shape: {interval.label} is not inside N/2<h<=N")
        if interval.rho <= 0.0:
            raise SystemExit(f"complement high shape: {interval.label} has nonpositive rho")
        if interval.stop > interval.start:
            first_ratio = (n - interval.start) / (interval.start + 1)
            last_ratio = (n - (interval.stop - 1)) / interval.stop
            if last_ratio > first_ratio:
                raise SystemExit(
                    f"complement high shape: binomial adjacent ratio increases on {interval.label}"
                )
            growth_bits = math.log2(first_ratio / last_ratio)
            min_growth = min(min_growth, growth_bits)
        else:
            growth_bits = float("inf")
        rows.append(
            {
                "range": [interval.start, interval.stop],
                "rho": interval.rho,
                "e0_endpoint": interval.start,
                "ege1_outer_volume_maxima": "endpoints",
                "convexity_growth_bits": finite_float(growth_bits),
            }
        )
    return ComplementHighShapeSummary(
        interval_count=len(intervals),
        min_convexity_growth_bits=min_growth,
        rows=rows,
    )


def check_t_monotonicity_sufficient_reduction(
    *,
    inner_spectrum: Path,
    b: int,
    turnoff_log2: float,
) -> TMonotonicitySummary:
    ordered = sorted(T_MONOTONICITY_INTERVALS, key=lambda interval: interval.h_min)
    previous_stop = 2000
    for interval in ordered:
        if interval.h_min != previous_stop + 1:
            raise SystemExit(
                "T monotonicity audit: interval cover gap or overlap before "
                f"{interval.h_min}--{interval.h_max}"
            )
        if interval.h_max < interval.h_min:
            raise SystemExit(f"T monotonicity audit: malformed interval {interval.h_min}--{interval.h_max}")
        if interval.lam <= 0.0 or interval.rho <= 0.0:
            raise SystemExit(f"T monotonicity audit: nonpositive pole on {interval.h_min}--{interval.h_max}")
        previous_stop = interval.h_max
    if previous_stop != 1148736:
        raise SystemExit(f"T monotonicity audit: expected cover through 1148736, got {previous_stop}")

    buckets = [
        (bucket.t_min, bucket.t_max)
        for bucket in INTERIOR_AUDIT_BUCKETS
        if bucket.name in {"gap_4001_8000", "gap_8001_12000"}
    ]
    bucket_names = {(bucket.t_min, bucket.t_max): bucket.name for bucket in INTERIOR_AUDIT_BUCKETS}
    entries = split_inner_entries(load_inner_spectrum(inner_spectrum), b)
    raw_rows = sufficient_reduction_rows(
        intervals=list(T_MONOTONICITY_INTERVALS),
        buckets=buckets,
        entries=entries,
        b=b,
        turnoff_log2=turnoff_log2,
    )
    if len(raw_rows) != len(T_MONOTONICITY_INTERVALS) * len(buckets):
        raise SystemExit(
            "T monotonicity audit: expected "
            f"{len(T_MONOTONICITY_INTERVALS) * len(buckets)} rows, got {len(raw_rows)}"
        )

    rows: list[dict[str, object]] = []
    for row in raw_rows:
        if row.min_slack <= 0.0:
            raise SystemExit(
                "T monotonicity audit: nonpositive slack on "
                f"{row.h_min}--{row.h_max}, T={row.bucket_t_min}--{row.bucket_t_max}: {row.min_slack}"
            )
        rows.append(
            {
                "h_range": [row.h_min, row.h_max],
                "lambda": row.lam,
                "rho": row.rho,
                "bucket": bucket_names[(row.bucket_t_min, row.bucket_t_max)],
                "T_range": [row.bucket_t_min, row.bucket_t_max],
                "no_new_slack_bits": finite_float(row.no_new_slack),
                "new_slack_bits": finite_float(row.new_slack),
                "min_slack_bits": finite_float(row.min_slack),
                "H_new": row.h_new,
                "T_new": row.t_new,
            }
        )
    finite_rows = [row for row in rows if row["min_slack_bits"] is not None]
    if not finite_rows:
        raise SystemExit("T monotonicity audit: no finite slack row found")
    worst = min(finite_rows, key=lambda row: float(row["min_slack_bits"]))
    return TMonotonicitySummary(
        interval_count=len(T_MONOTONICITY_INTERVALS),
        bucket_count=len(buckets),
        row_count=len(rows),
        min_slack_bits=float(worst["min_slack_bits"]),
        worst_row=worst,
        rows=rows,
    )


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


def check_early_postprefix_interval(interval: EarlyPostprefixInterval, tolerance: float) -> float:
    result = subprocess.run(
        interval.command(),
        cwd=ROOT.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    actual = parse_piecewise_total(result.stdout)
    check_close(f"early_postprefix_interval_{interval.label}_recompute", actual, interval.log2_value, tolerance)
    return actual


def check_early_accelerated_interval(interval: EarlyAcceleratedInterval, tolerance: float) -> float:
    result = subprocess.run(
        interval.command(),
        cwd=ROOT.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    actual = parse_piecewise_total(result.stdout)
    check_close(f"early_accelerated_interval_{interval.label}_recompute", actual, interval.log2_value, tolerance)
    return actual


def check_late_postprefix_interval(interval: LatePostprefixInterval, tolerance: float) -> float:
    result = subprocess.run(
        interval.command(),
        cwd=ROOT.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    actual = parse_piecewise_total(result.stdout)
    check_close(f"late_postprefix_interval_{interval.label}_recompute", actual, interval.log2_value, tolerance)
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
    parser.add_argument("--inner-spectrum", type=Path, default=ROOT / "EBCH128_64.wd")
    parser.add_argument("--exact-outer-blocks", type=int, default=4096)
    parser.add_argument("--exact-outer-h-max", type=int, default=500)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--late-prefix-h-max", type=int, default=500)
    parser.add_argument("--skip-prefix-interior-audit", action="store_true")
    parser.add_argument("--prefix-interior-rows", type=int, default=1500)
    parser.add_argument("--prefix-interior-h-min", type=int, default=0)
    parser.add_argument("--prefix-interior-h-max", type=int, default=499)
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
        "--print-early-postprefix-commands",
        action="store_true",
        help="Print the recomputation command for every checked early post-prefix interval.",
    )
    parser.add_argument(
        "--print-early-accelerated-commands",
        action="store_true",
        help="Print the recomputation command for every accelerated early interval.",
    )
    parser.add_argument(
        "--print-late-postprefix-commands",
        action="store_true",
        help="Print the recomputation command for every ultra-late post-prefix interval.",
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
    parser.add_argument(
        "--check-early-postprefix-intervals",
        help=(
            "Comma-separated early post-prefix interval labels to recompute, or 'all'. "
            "Example: 501--2000,2001--7858."
        ),
    )
    parser.add_argument(
        "--check-early-accelerated-intervals",
        help=(
            "Comma-separated accelerated early interval labels to recompute, or 'all'. "
            "Example: 100001--110000,750001--1048576."
        ),
    )
    parser.add_argument(
        "--check-late-postprefix-intervals",
        help=(
            "Comma-separated ultra-late post-prefix interval labels to recompute, or 'all'. "
            "Example: 501--2000,250001--380736."
        ),
    )
    parser.add_argument(
        "--write-manifest-json",
        type=Path,
        help="Write a machine-readable manifest of the checked row families, constants, and totals.",
    )
    parser.add_argument("--tolerance", type=float, default=5e-6)
    args = parser.parse_args()

    manifest: dict[str, object] = {
        "schema": "fullsplit_finite_ledger.v1",
        "construction": "RM512_256 block outer with full-split dense inner",
        "target": {
            "N": args.N,
            "delta": 0.09,
            "quantity": "first-moment upper bound for checked finite dense+dense rows",
        },
        "parameters": {
            "block_bits": args.block_bits,
            "late_blocks": args.late_blocks,
            "late_coordinates": args.block_bits * args.late_blocks,
            "exact_outer_blocks": args.exact_outer_blocks,
            "exact_outer_h_max": args.exact_outer_h_max,
            "global_turnoff_log2": GLOBAL_TURNOFF_LOG2,
            "high_interval_turnoff_adjustment_bits": HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS,
            "exact_outer_local_spectrum_csv": manifest_path(args.exact_outer_local_spectrum_csv),
            "inner_spectrum": manifest_path(args.inner_spectrum),
        },
        "checks": {},
        "row_families": [],
        "totals": {},
    }

    if (
        args.print_high_interval_commands
        or args.print_early_postprefix_commands
        or args.print_early_accelerated_commands
        or args.print_late_postprefix_commands
    ):
        if args.print_early_postprefix_commands:
            print("Early post-prefix interval recomputation commands")
            for interval in EARLY_POSTPREFIX_INTERVALS:
                print(f"early_postprefix_interval_{interval.label}_command,{interval.display_command()}")
        if args.print_early_accelerated_commands:
            print("Accelerated early interval recomputation commands")
            for interval in EARLY_ACCELERATED_INTERVALS:
                print(f"early_accelerated_interval_{interval.label}_command,{interval.display_command()}")
        if args.print_late_postprefix_commands:
            print("Ultra-late post-prefix interval recomputation commands")
            for interval in LATE_POSTPREFIX_INTERVALS:
                print(f"late_postprefix_interval_{interval.label}_command,{interval.display_command()}")
        if not args.print_high_interval_commands:
            if (
                not args.check_early_postprefix_intervals
                and not args.check_early_accelerated_intervals
                and not args.check_late_postprefix_intervals
            ):
                return 0
        else:
            print("High-interval recomputation commands")
            for interval in HIGH_INTERVALS_RAW:
                print(f"interval_{interval.label}_command,{interval.display_command()}")
            print("Complement-high interval recomputation commands")
            for interval in COMPLEMENT_HIGH_INTERVALS:
                print(f"complement_interval_{interval.label}_command,{interval.display_command()}")
        if (
            not args.check_high_intervals
            and not args.check_complement_high_intervals
            and not args.check_early_postprefix_intervals
            and not args.check_early_accelerated_intervals
            and not args.check_late_postprefix_intervals
        ):
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
            expected_h_min=args.prefix_interior_h_min,
            expected_h_max=args.prefix_interior_h_max,
            tolerance=args.prefix_interior_tolerance,
            reproduction_tolerance=args.prefix_interior_reproduction_tolerance,
        )
        print(f"prefix_interior_audit_rows,{interior.rows}")
        print(f"prefix_interior_audit_h_range,{interior.h_min}--{interior.h_max}")
        print(
            "prefix_interior_audit_buckets,"
            + ";".join(bucket.label for bucket in INTERIOR_AUDIT_BUCKETS)
        )
        print("prefix_interior_audit_left_endpoint_maxima,PASS")
        print(f"prefix_interior_audit_tolerance_bits,{args.prefix_interior_tolerance:.12g}")
        print(f"prefix_interior_audit_worst_slack_bits,{interior.worst_diff_slack:.12g}")
        print(
            "prefix_interior_audit_reproduction_tolerance_bits,"
            f"{args.prefix_interior_reproduction_tolerance:.12g}"
        )
        print(
            "prefix_interior_audit_worst_repro_slack_bits,"
            f"{interior.worst_reproduction_slack:.12g}"
        )
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
        manifest["checks"]["prefix_interior_audit"] = {  # type: ignore[index]
            "rows": interior.rows,
            "h_range": [interior.h_min, interior.h_max],
            "h_count": interior.h_count,
            "buckets": [
                {
                    "name": bucket.name,
                    "T_range": [bucket.t_min, bucket.t_max],
                    "T_count": bucket.t_count,
                    "rows": interior.bucket_counts[bucket.name],
                }
                for bucket in INTERIOR_AUDIT_BUCKETS
            ],
            "left_endpoint_maxima": interior.left_endpoint_maxima,
            "tolerance_bits": args.prefix_interior_tolerance,
            "worst_slack_bits": interior.worst_diff_slack,
            "reproduction_tolerance_bits": args.prefix_interior_reproduction_tolerance,
            "worst_reproduction_slack_bits": interior.worst_reproduction_slack,
            "worst": interior.worst,
            "worst_reproduction": interior.worst_repro,
        }

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
    add_row(
        manifest,
        name=f"late_prefix_T_lt_{args.late_blocks}_cap_h_le_{args.late_prefix_h_max}",
        role="diagnostic_cap",
        h_range=(1, args.late_prefix_h_max),
        rows=late_rows,
        log2_value=late_prefix,
        notes="Placement-only ultra-late cap using the outer-prefix CSV coefficients.",
        peak=late_peak,
    )

    exact_outer_certificate = prefix_coefficient_certificate(
        load_local_spectrum(args.exact_outer_local_spectrum_csv, args.exact_outer_h_max),
        blocks=args.exact_outer_blocks,
        h_max=args.exact_outer_h_max,
    )
    exact_outer_coeffs = exact_outer_certificate.coeffs
    exact_support = exact_outer_certificate.first_positive_le_80
    if exact_outer_certificate.min_positive != 32 or exact_outer_certificate.next_positive_after_min != 48:
        raise SystemExit(
            "exact RM support: expected first positive weights 32,48; "
            f"got {exact_outer_certificate.min_positive},"
            f"{exact_outer_certificate.next_positive_after_min}"
        )
    if exact_outer_certificate.zero_prefix_stop != 31:
        raise SystemExit(
            "exact RM support: expected zero prefix through h=31; "
            f"got h={exact_outer_certificate.zero_prefix_stop}"
        )
    if (exact_outer_certificate.gap_after_min_start, exact_outer_certificate.gap_after_min_stop) != (33, 47):
        raise SystemExit(
            "exact RM support: expected gap h=33--47; "
            f"got {exact_outer_certificate.gap_after_min_start}--"
            f"{exact_outer_certificate.gap_after_min_stop}"
        )
    manifest["checks"]["exact_outer_support"] = {  # type: ignore[index]
        "local_nonzero_terms": exact_outer_certificate.local_nonzero_terms,
        "local_min_positive": exact_outer_certificate.local_min_positive,
        "local_first_positive_le_80": exact_outer_certificate.local_first_positive_le_80,
        "coefficient_build_exponent": exact_outer_certificate.trace.exponent,
        "coefficient_build_multiply_steps": exact_outer_certificate.trace.multiply_steps,
        "coefficient_build_square_steps": exact_outer_certificate.trace.square_steps,
        "coefficient_build_max_coefficient_bits": exact_outer_certificate.trace.max_coefficient_bits,
        "coefficient_build_max_coefficient_weight": exact_outer_certificate.trace.max_coefficient_weight,
        "support_count": exact_outer_certificate.support_count,
        "positive_count": exact_outer_certificate.positive_count,
        "min_positive": exact_outer_certificate.min_positive,
        "next_positive_after_min": exact_outer_certificate.next_positive_after_min,
        "zero_prefix": [1, exact_outer_certificate.zero_prefix_stop],
        "gap_after_min": [
            exact_outer_certificate.gap_after_min_start,
            exact_outer_certificate.gap_after_min_stop,
        ],
        "first_positive_le_80": exact_outer_certificate.first_positive_le_80,
        "split_h": exact_outer_certificate.split_h,
        "split_coefficient_bits": exact_outer_certificate.split_coefficient_bits,
        "next_coefficient_bits": exact_outer_certificate.next_coefficient_bits,
    }
    print(f"exact_outer_local_nonzero_terms,{exact_outer_certificate.local_nonzero_terms}")
    print(f"exact_outer_local_min_positive,{exact_outer_certificate.local_min_positive}")
    print(
        "exact_outer_local_first_positive,"
        f"{';'.join(str(h) for h in exact_outer_certificate.local_first_positive_le_80)}"
    )
    print(f"exact_outer_coefficient_build_exponent,{exact_outer_certificate.trace.exponent}")
    print(f"exact_outer_coefficient_build_multiply_steps,{exact_outer_certificate.trace.multiply_steps}")
    print(f"exact_outer_coefficient_build_square_steps,{exact_outer_certificate.trace.square_steps}")
    print(
        "exact_outer_coefficient_build_max_coefficient_bits,"
        f"{exact_outer_certificate.trace.max_coefficient_bits}"
    )
    print(
        "exact_outer_coefficient_build_max_coefficient_weight,"
        f"{exact_outer_certificate.trace.max_coefficient_weight}"
    )
    print(f"exact_outer_support_count,{exact_outer_certificate.support_count}")
    print(f"exact_outer_support_positive_count,{exact_outer_certificate.positive_count}")
    print(f"exact_outer_support_min_positive,{exact_outer_certificate.min_positive}")
    print(f"exact_outer_support_next_positive_after_min,{exact_outer_certificate.next_positive_after_min}")
    print(f"exact_outer_support_zero_prefix,1--{exact_outer_certificate.zero_prefix_stop}")
    print(
        "exact_outer_support_gap_after_min,"
        f"{exact_outer_certificate.gap_after_min_start}--{exact_outer_certificate.gap_after_min_stop}"
    )
    print(f"exact_outer_support_first_positive,{';'.join(str(h) for h in exact_support)}")
    print(f"exact_outer_support_split_h,{exact_outer_certificate.split_h}")
    print(f"exact_outer_support_split_coefficient_bits,{exact_outer_certificate.split_coefficient_bits}")
    print(f"exact_outer_support_next_coefficient_bits,{exact_outer_certificate.next_coefficient_bits}")
    exact_late_prefix = exact_late_prefix_sum(
        exact_outer_coeffs,
        n=args.N,
        b=args.block_bits,
        late_blocks=args.late_blocks,
        h_max=args.late_prefix_h_max,
    )
    check_close("late_prefix_T_lt_5949_exact_outer", exact_late_prefix.total_log2, -41.113442, args.tolerance)
    check_close("late_prefix_T_lt_5949_exact_outer_split", exact_late_prefix.split_log2, -41.113442, args.tolerance)
    check_close(
        "late_prefix_T_lt_5949_exact_outer_above_split",
        exact_late_prefix.above_split_log2,
        -66.416502,
        args.tolerance,
    )
    print(f"late_prefix_T_lt_{args.late_blocks}_exact_outer_rows,{exact_late_prefix.rows}")
    print(f"late_prefix_T_lt_{args.late_blocks}_exact_outer_log2,{exact_late_prefix.total_log2:.6f}")
    print(f"late_prefix_T_lt_{args.late_blocks}_exact_outer_split_h,{exact_late_prefix.split_h}")
    print(f"late_prefix_T_lt_{args.late_blocks}_exact_outer_split_log2,{exact_late_prefix.split_log2:.6f}")
    print(
        f"late_prefix_T_lt_{args.late_blocks}_exact_outer_above_split_log2,"
        f"{exact_late_prefix.above_split_log2:.6f}"
    )
    print(
        f"late_prefix_T_lt_{args.late_blocks}_exact_outer_peak,"
        f"h={exact_late_prefix.peak_h},"
        f"outer_log2={exact_late_prefix.peak_outer_log2:.6f},"
        f"late_log2={exact_late_prefix.peak_late_log2:.6f},"
        f"term_log2={exact_late_prefix.peak_term_log2:.6f}"
    )
    add_row(
        manifest,
        name=f"late_prefix_T_lt_{args.late_blocks}_exact_outer_h_le_{args.late_prefix_h_max}",
        role="active_h_le_500",
        h_range=(1, args.late_prefix_h_max),
        rows=exact_late_prefix.rows,
        log2_value=exact_late_prefix.total_log2,
        notes="Exact RM direct-sum outer coefficients for the ultra-late prefix.",
        split_h=exact_late_prefix.split_h,
        split_log2=exact_late_prefix.split_log2,
        above_split_log2=exact_late_prefix.above_split_log2,
        peak={
            "h": exact_late_prefix.peak_h,
            "outer_log2": exact_late_prefix.peak_outer_log2,
            "late_log2": exact_late_prefix.peak_late_log2,
            "term_log2": exact_late_prefix.peak_term_log2,
        },
    )

    parts: dict[str, float] = {}
    for item in csv_ledgers:
        total, rows, peak = csv_logsum(item.path, item.column)
        check_close(item.name, total, item.expected_log2, args.tolerance)
        parts[item.name] = total
        print(f"{item.name}_rows,{rows}")
        print(f"{item.name}_log2,{total:.6f}")
        print_peak(item.name, peak, item.column)
        add_row(
            manifest,
            name=item.name,
            role="raw_csv_check",
            rows=rows,
            log2_value=total,
            source_csv=manifest_path(item.path),
            column=item.column,
            expected_log2=item.expected_log2,
            peak=peak,
        )

    exact_outer_expected = {
        "prefix_32_500_e_le8": -37.383345,
        "prefix_32_500_e_ge9_tail": -284.004805,
        "early_32_500_e_le16_uniformsurv": -86.910456,
        "early_32_500_e_ge17_tail": -319.977808,
    }
    exact_outer_split_expected = {
        "prefix_32_500_e_le8": (-37.383482, -50.738253),
        "prefix_32_500_e_ge9_tail": (-540.423048, -284.004805),
        "early_32_500_e_le16_uniformsurv": (-86.910458, -106.352289),
        "early_32_500_e_ge17_tail": (-1019.441595, -319.977808),
    }
    dominant_prefix_expected = {
        "prefix_32_500_e_le8": {
            "selector": (32, 1, 31, 1, 4000, 5949),
            "term_log2": -37.3857668415005,
            "total_minus_peak_bits": 0.002422317512923655,
            "total_remainder_log2": -46.602718187648264,
            "split_minus_peak_bits": 0.0022846069119495382,
            "split_remainder_log2": -46.68722908748783,
            "above_split_gap_bits": 13.352485666805116,
        },
    }
    exact_outer_parts: dict[str, float] = {}
    for item in csv_ledgers:
        if item.name not in exact_outer_expected:
            continue
        summary = reweighted_ledger_sum(item.path, exact_outer_coeffs)
        check_close(f"{item.name}_exact_outer", summary.exact_log2, exact_outer_expected[item.name], args.tolerance)
        split_expected, above_split_expected = exact_outer_split_expected[item.name]
        check_close(f"{item.name}_exact_outer_split", summary.split_log2, split_expected, args.tolerance)
        check_close(
            f"{item.name}_exact_outer_above_split",
            summary.above_split_log2,
            above_split_expected,
            args.tolerance,
        )
        exact_outer_parts[item.name] = summary.exact_log2
        print(f"{item.name}_exact_outer_live_rows,{summary.live_rows}")
        print(f"{item.name}_exact_outer_log2,{summary.exact_log2:.6f}")
        print(f"{item.name}_exact_outer_split_h,{summary.split_h}")
        print(f"{item.name}_exact_outer_split_log2,{summary.split_log2:.6f}")
        print(f"{item.name}_exact_outer_above_split_log2,{summary.above_split_log2:.6f}")
        print(
            f"{item.name}_exact_outer_peak,"
            f"outer_weight={summary.peak_h},first_r={summary.peak_first_r},"
            f"gap_min={summary.peak_gap_min},gap_max={summary.peak_gap_max},"
            f"outer_log2={summary.peak_outer_exact_log2:.6f},"
            f"term_log2={summary.peak_term_log2:.6f}"
        )
        dom = summary.dominant
        dominant_payload: dict[str, object] | None = None
        if item.name in dominant_prefix_expected:
            expected = dominant_prefix_expected[item.name]
            selector = (dom.h, dom.first_r, dom.remaining_ones, dom.gap_min, dom.gap_max, dom.bucket_T)
            if selector != expected["selector"]:
                raise SystemExit(
                    f"{item.name}_dominant: expected selector {expected['selector']}, got {selector}"
                )
            for field in [
                "term_log2",
                "total_minus_peak_bits",
                "total_remainder_log2",
                "split_minus_peak_bits",
                "split_remainder_log2",
                "above_split_gap_bits",
            ]:
                check_close(f"{item.name}_dominant_{field}", getattr(dom, field), expected[field], args.tolerance)
            print(
                f"{item.name}_dominant_row,"
                f"outer_weight={dom.h},first_r={dom.first_r},"
                f"remaining_ones={dom.remaining_ones},"
                f"gap_min={dom.gap_min},gap_max={dom.gap_max},"
                f"bucket_T={dom.bucket_T},inner_mode={dom.inner_mode}"
            )
            print(f"{item.name}_dominant_old_outer_log2,{dom.old_outer_log2:.6f}")
            print(f"{item.name}_dominant_exact_outer_log2,{dom.exact_outer_log2:.6f}")
            print(f"{item.name}_dominant_placement_log2,{dom.placement_log2:.6f}")
            print(f"{item.name}_dominant_inner_log2,{dom.inner_log2:.6f}")
            print(f"{item.name}_dominant_term_log2,{dom.term_log2:.6f}")
            print(f"{item.name}_dominant_total_minus_peak_bits,{dom.total_minus_peak_bits:.6f}")
            print(f"{item.name}_dominant_total_remainder_log2,{dom.total_remainder_log2:.6f}")
            print(f"{item.name}_dominant_split_minus_peak_bits,{dom.split_minus_peak_bits:.6f}")
            print(f"{item.name}_dominant_split_remainder_log2,{dom.split_remainder_log2:.6f}")
            print(f"{item.name}_dominant_above_split_gap_bits,{dom.above_split_gap_bits:.6f}")
            dominant_payload = {
                "outer_weight": dom.h,
                "first_r": dom.first_r,
                "remaining_ones": dom.remaining_ones,
                "gap_min": dom.gap_min,
                "gap_max": dom.gap_max,
                "bucket_T": dom.bucket_T,
                "inner_mode": dom.inner_mode,
                "old_outer_log2": dom.old_outer_log2,
                "exact_outer_log2": dom.exact_outer_log2,
                "placement_log2": dom.placement_log2,
                "inner_log2": dom.inner_log2,
                "term_log2": dom.term_log2,
                "total_minus_peak_bits": dom.total_minus_peak_bits,
                "total_remainder_log2": dom.total_remainder_log2,
                "split_minus_peak_bits": dom.split_minus_peak_bits,
                "split_remainder_log2": dom.split_remainder_log2,
                "above_split_gap_bits": dom.above_split_gap_bits,
            }
        add_row_extra = {}
        if dominant_payload is not None:
            add_row_extra["dominant"] = dominant_payload
        add_row(
            manifest,
            name=f"{item.name}_exact_outer",
            role="active_h_le_500",
            rows=summary.live_rows,
            log2_value=summary.exact_log2,
            source_csv=manifest_path(item.path),
            split_h=summary.split_h,
            split_log2=summary.split_log2,
            above_split_log2=summary.above_split_log2,
            peak={
                "outer_weight": summary.peak_h,
                "first_r": summary.peak_first_r,
                "gap_min": summary.peak_gap_min,
                "gap_max": summary.peak_gap_max,
                "outer_log2": summary.peak_outer_exact_log2,
                "term_log2": summary.peak_term_log2,
            },
            **add_row_extra,
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
    manifest["checks"]["prefix_ridge_ratio"] = {  # type: ignore[index]
        "rows": ridge_ratio.rows,
        "ridge_rows": ridge_ratio.ridge_rows,
        "ridge_exact_log2": ridge_ratio.ridge_exact_log2,
        "remainder_exact_log2": ridge_ratio.remainder_exact_log2,
        "first_ratio": ridge_ratio.first_ratio,
        "tail_ratio": ridge_ratio.tail_ratio,
        "first_outer_ratio": ridge_ratio.first_outer_ratio,
        "first_placement_ratio": ridge_ratio.first_placement_ratio,
        "first_inner_ratio": ridge_ratio.first_inner_ratio,
        "tail_outer_ratio": ridge_ratio.tail_outer_ratio,
        "tail_placement_ratio": ridge_ratio.tail_placement_ratio,
        "tail_inner_ratio": ridge_ratio.tail_inner_ratio,
        "geometric_bound_log2": ridge_ratio.geometric_log2,
        "total_ratio_bound_log2": ridge_ratio.total_bound_log2,
        "total_ratio_bound_slack_bits": ridge_ratio.total_bound_log2 - ridge_ratio.total_exact_log2,
        "inner_span_bits": ridge_ratio.inner_span_bits,
    }

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
        print(f"prefix_placement_ratio_threshold_num,{placement.threshold_num}")
        print(f"prefix_placement_ratio_threshold_den,{placement.threshold_den}")
        print("prefix_placement_ratio_exact_cross_multiply_status,PASS")
        print(f"prefix_placement_ratio_exact_comparisons,{placement.exact_comparisons}")
        print(f"prefix_placement_ratio_peak_h,{placement.peak_h}")
        print(f"prefix_placement_ratio_peak_to_h,{placement.peak_to_h}")
        print(f"prefix_placement_ratio_peak,{placement.peak_ratio:.12g}")
        print(f"prefix_placement_ratio_peak_log2,{placement.peak_log2:.12g}")
        print(f"prefix_placement_ratio_peak_num_bits,{placement.peak_num_bits}")
        print(f"prefix_placement_ratio_peak_den_bits,{placement.peak_den_bits}")
        print(f"prefix_placement_ratio_threshold_slack_bits,{placement.threshold_slack_bits}")
        print(f"prefix_placement_ratio_threshold_gap_log2,{placement.threshold_gap_log2:.12g}")
        print(f"prefix_placement_endpoint_bound_at_h_min,{placement.endpoint_ratio_at_h_min:.12g}")
        print(f"prefix_placement_endpoint_bound_at_h_min_log2,{placement.endpoint_ratio_at_h_min_log2:.12g}")
        print(f"prefix_placement_endpoint_slack_factor,{placement.endpoint_slack_factor:.12g}")
        manifest["checks"]["prefix_placement_ratio"] = {  # type: ignore[index]
            "threshold": args.prefix_placement_ratio_threshold,
            "threshold_num": placement.threshold_num,
            "threshold_den": placement.threshold_den,
            "exact_cross_multiply_status": "PASS",
            "exact_comparisons": placement.exact_comparisons,
            "gap_range": [args.prefix_placement_gap_min, args.prefix_placement_gap_max],
            "h_range": [args.prefix_placement_h_min, args.prefix_placement_h_max],
            "peak_h": placement.peak_h,
            "peak_to_h": placement.peak_to_h,
            "peak_ratio": placement.peak_ratio,
            "peak_log2": placement.peak_log2,
            "peak_num_bits": placement.peak_num_bits,
            "peak_den_bits": placement.peak_den_bits,
            "threshold_slack_bits": placement.threshold_slack_bits,
            "threshold_gap_log2": placement.threshold_gap_log2,
            "endpoint_bound_at_h_min": placement.endpoint_ratio_at_h_min,
            "endpoint_bound_at_h_min_log2": placement.endpoint_ratio_at_h_min_log2,
            "endpoint_slack_factor": placement.endpoint_slack_factor,
        }

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

    high_cover = check_contiguous_cover(
        "high interval cover",
        HIGH_INTERVALS_RAW,
        expected_start=2001,
        expected_stop=1148736,
    )
    for interval in HIGH_INTERVALS_RAW:
        if interval.lam <= 0.0 or interval.rho <= 0.0:
            raise SystemExit(f"high interval cover: interval {interval.label} has nonpositive pole")
        check_close(
            f"interval_{interval.label}_turnoff_adjustment",
            interval.adjusted_log2 - interval.raw_log2,
            HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS,
            1e-9,
        )
    print("high_interval_cover_status,PASS")
    print(f"high_interval_cover_range,{high_cover.cover_start}--{high_cover.cover_stop}")
    print(f"high_interval_count,{high_cover.count}")
    print(f"high_interval_turnoff_adjustment_bits,{HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS:.6f}")
    manifest["checks"]["high_interval_cover"] = {  # type: ignore[index]
        "status": "PASS",
        "cover_range": [high_cover.cover_start, high_cover.cover_stop],
        "interval_count": high_cover.count,
        "turnoff_adjustment_bits": HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS,
        "intervals": [
            {
                "range": [interval.start, interval.stop],
                "lambda": interval.lam,
                "rho": interval.rho,
                "raw_log2": interval.raw_log2,
                "adjusted_log2": interval.adjusted_log2,
                "command": interval.display_command(),
            }
            for interval in HIGH_INTERVALS_RAW
        ],
    }
    high_feasible_t = check_high_feasible_min_t(b=args.block_bits, high_cover=high_cover)
    print("high_feasible_min_T_status,PASS")
    print(f"high_feasible_min_T_far_cutoff_h,{high_feasible_t.far_cutoff_h}")
    for name, t_min, t_max, cutoff_h in high_feasible_t.bucket_cutoffs:
        print(f"high_feasible_min_T_{name}_range,{t_min}--{t_max}")
        print(f"high_feasible_min_T_{name}_cutoff_h,{cutoff_h}")
    manifest["checks"]["high_feasible_min_T"] = {  # type: ignore[index]
        "status": "PASS",
        "far_cutoff_h": high_feasible_t.far_cutoff_h,
        "rule": "T_eff=max(T_min,ceil(H/block_bits)); h cutoff is block_bits*T_max+block_bits",
        "buckets": [
            {
                "name": name,
                "T_range": [t_min, t_max],
                "cutoff_h": cutoff_h,
            }
            for name, t_min, t_max, cutoff_h in high_feasible_t.bucket_cutoffs
        ],
    }
    t_mono = check_t_monotonicity_sufficient_reduction(
        inner_spectrum=args.inner_spectrum,
        b=args.block_bits,
        turnoff_log2=GLOBAL_TURNOFF_LOG2,
    )
    worst_t = t_mono.worst_row["T_range"]
    worst_h = t_mono.worst_row["h_range"]
    print("t_monotonicity_sufficient_status,PASS")
    print(f"t_monotonicity_sufficient_rows,{t_mono.row_count}")
    print(f"t_monotonicity_sufficient_min_slack_bits,{t_mono.min_slack_bits:.6f}")
    print(
        "t_monotonicity_sufficient_worst,"
        f"h={worst_h[0]}--{worst_h[1]},"
        f"bucket={t_mono.worst_row['bucket']},"
        f"T={worst_t[0]}--{worst_t[1]},"
        f"lambda={t_mono.worst_row['lambda']:g},"
        f"rho={t_mono.worst_row['rho']:g},"
        f"slack={t_mono.min_slack_bits:.6f}"
    )
    manifest["checks"]["t_monotonicity_sufficient_reduction"] = {  # type: ignore[index]
        "status": "PASS",
        "interval_count": t_mono.interval_count,
        "bucket_count": t_mono.bucket_count,
        "row_count": t_mono.row_count,
        "min_slack_bits": t_mono.min_slack_bits,
        "worst_row": t_mono.worst_row,
        "rule": "two-case sufficient endpoint reduction for fixed-pole e>=1 envelope monotonicity in T",
        "rows": t_mono.rows,
    }

    high_total = float("-inf")
    for label, lam, rho, value in HIGH_INTERVALS:
        high_total = log2add(high_total, value)
        print(f"interval_{label}_lambda,{lam:g}")
        print(f"interval_{label}_rho,{rho:g}")
        print(f"interval_{label}_log2,{value:.6f}")
        add_row(
            manifest,
            name=f"interval_{label}",
            role="active_postprefix_late_window_fixed_pole",
            h_range=tuple(int(x) for x in label.split("--")),  # type: ignore[arg-type]
            log2_value=value,
            command=next(interval.display_command() for interval in HIGH_INTERVALS_RAW if interval.label == label),
            lambda_value=lam,
            rho=rho,
            raw_log2=next(interval.raw_log2 for interval in HIGH_INTERVALS_RAW if interval.label == label),
            turnoff_adjustment_bits=HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS,
        )
    print(f"interval_2001_1148736_total_log2,{high_total:.6f}")
    print("feasible_far_bucket_cutoff_h,1148736")

    complement_cover = check_contiguous_cover(
        "complement high interval cover",
        COMPLEMENT_HIGH_INTERVALS,
        expected_start=args.N // 2 + 1,
        expected_stop=args.N,
    )
    for interval in COMPLEMENT_HIGH_INTERVALS:
        if interval.rho <= 0.0:
            raise SystemExit(f"complement high interval cover: interval {interval.label} has nonpositive rho")
    high_complement_overlap = interval_overlap(high_cover, complement_cover)
    if high_complement_overlap != (args.N // 2 + 1, high_cover.cover_stop):
        raise SystemExit(
            "high/complement cover: expected deliberate overlap "
            f"{args.N // 2 + 1}--{high_cover.cover_stop}, got {high_complement_overlap}"
        )
    overlap_count = high_complement_overlap[1] - high_complement_overlap[0] + 1
    print("complement_high_cover_status,PASS")
    print(f"complement_high_cover_range,{complement_cover.cover_start}--{complement_cover.cover_stop}")
    print(f"complement_high_interval_count,{complement_cover.count}")
    print(f"high_complement_overlap_range,{high_complement_overlap[0]}--{high_complement_overlap[1]}")
    print(f"high_complement_overlap_count,{overlap_count}")
    manifest["checks"]["complement_high_cover"] = {  # type: ignore[index]
        "status": "PASS",
        "cover_range": [complement_cover.cover_start, complement_cover.cover_stop],
        "interval_count": complement_cover.count,
        "high_overlap_range": [high_complement_overlap[0], high_complement_overlap[1]],
        "high_overlap_count": overlap_count,
        "overlap_treatment": "deliberate union-bound overcount",
        "intervals": [
            {
                "range": [interval.start, interval.stop],
                "rho": interval.rho,
                "log2": interval.log2_value,
                "command": interval.display_command(),
            }
            for interval in COMPLEMENT_HIGH_INTERVALS
        ],
    }
    complement_shape = check_complement_high_shapes(COMPLEMENT_HIGH_INTERVALS, n=args.N)
    print("complement_high_shape_status,PASS")
    print(f"complement_high_shape_interval_count,{complement_shape.interval_count}")
    print(f"complement_high_shape_min_convexity_growth_bits,{complement_shape.min_convexity_growth_bits:.12g}")
    manifest["checks"]["complement_high_shape"] = {  # type: ignore[index]
        "status": "PASS",
        "interval_count": complement_shape.interval_count,
        "e0_outer_endpoint_rule": "z<1 makes complement-side e=0 outer bound largest at h_min",
        "ege1_outer_volume_rule": "discrete convexity makes the outer/volume/rho term largest at an endpoint",
        "min_convexity_growth_bits": complement_shape.min_convexity_growth_bits,
        "intervals": complement_shape.rows,
    }

    complement_high_total = float("-inf")
    for interval in COMPLEMENT_HIGH_INTERVALS:
        complement_high_total = log2add(complement_high_total, interval.log2_value)
        print(f"complement_interval_{interval.label}_rho,{interval.rho:g}")
        print(f"complement_interval_{interval.label}_log2,{interval.log2_value:.6f}")
        add_row(
            manifest,
            name=f"complement_interval_{interval.label}",
            role="active_complement_high_density",
            h_range=(interval.start, interval.stop),
            log2_value=interval.log2_value,
            command=interval.display_command(),
            rho=interval.rho,
        )
    print(f"complement_interval_1048577_2097152_total_log2,{complement_high_total:.6f}")

    early_postprefix_total = float("-inf")
    for interval in EARLY_POSTPREFIX_INTERVALS:
        early_postprefix_total = log2add(early_postprefix_total, interval.log2_value)
        if interval.lambdas is not None:
            print(f"early_postprefix_interval_{interval.label}_lambda,{interval.lambdas}")
        if interval.rhos is not None:
            print(f"early_postprefix_interval_{interval.label}_rho,{interval.rhos}")
        print(f"early_postprefix_interval_{interval.label}_log2,{interval.log2_value:.6f}")
        add_row(
            manifest,
            name=f"early_postprefix_interval_{interval.label}",
            role="active_postprefix_early_endpoint",
            h_range=(interval.start, interval.stop),
            log2_value=interval.log2_value,
            command=interval.display_command(),
            lambdas=interval.lambdas,
            rhos=interval.rhos,
        )
    check_close("early_postprefix_501_75000", early_postprefix_total, -380.829185, args.tolerance)
    print(f"early_postprefix_501_75000_total_log2,{early_postprefix_total:.6f}")

    early_accelerated_cover = check_early_accelerated_cover(
        EARLY_ACCELERATED_INTERVALS,
        expected_start=75001,
        expected_stop=1048576,
        expected_endpoint_stop=350000,
        expected_paired_start=350001,
    )
    print("early_accelerated_cover_status,PASS")
    print(
        "early_accelerated_cover_range,"
        f"{early_accelerated_cover.cover_start}--{early_accelerated_cover.cover_stop}"
    )
    print(f"early_accelerated_interval_count,{early_accelerated_cover.count}")
    print(
        "early_accelerated_endpoint_cover,"
        f"{early_accelerated_cover.endpoint_start}--{early_accelerated_cover.endpoint_stop}"
    )
    print(f"early_accelerated_endpoint_interval_count,{early_accelerated_cover.endpoint_count}")
    print(
        "early_accelerated_paired_cover,"
        f"{early_accelerated_cover.paired_start}--{early_accelerated_cover.paired_stop}"
    )
    print(f"early_accelerated_paired_interval_count,{early_accelerated_cover.paired_count}")
    manifest["checks"]["early_accelerated_cover"] = {  # type: ignore[index]
        "status": "PASS",
        "cover_range": [early_accelerated_cover.cover_start, early_accelerated_cover.cover_stop],
        "interval_count": early_accelerated_cover.count,
        "endpoint_range": [early_accelerated_cover.endpoint_start, early_accelerated_cover.endpoint_stop],
        "endpoint_interval_count": early_accelerated_cover.endpoint_count,
        "paired_range": [early_accelerated_cover.paired_start, early_accelerated_cover.paired_stop],
        "paired_interval_count": early_accelerated_cover.paired_count,
        "intervals": [
            {
                "range": [interval.start, interval.stop],
                "kind": interval.kind,
                "helper": interval.script_name,
                "lambda": interval.lam,
                "rho": interval.rho,
                "log2": interval.log2_value,
                "command": interval.display_command(),
            }
            for interval in EARLY_ACCELERATED_INTERVALS
        ],
    }

    early_accelerated_total = float("-inf")
    endpoint_total = float("-inf")
    paired_total = float("-inf")
    for interval in EARLY_ACCELERATED_INTERVALS:
        early_accelerated_total = log2add(early_accelerated_total, interval.log2_value)
        if interval.kind == "endpoint":
            endpoint_total = log2add(endpoint_total, interval.log2_value)
        else:
            paired_total = log2add(paired_total, interval.log2_value)
        print(f"early_accelerated_interval_{interval.label}_kind,{interval.kind}")
        print(f"early_accelerated_interval_{interval.label}_lambda,{interval.lam:g}")
        print(f"early_accelerated_interval_{interval.label}_rho,{interval.rho:g}")
        print(f"early_accelerated_interval_{interval.label}_log2,{interval.log2_value:.6f}")
        add_row(
            manifest,
            name=f"early_accelerated_interval_{interval.label}",
            role=f"active_postprefix_early_{interval.kind}",
            h_range=(interval.start, interval.stop),
            log2_value=interval.log2_value,
            command=interval.display_command(),
            lambda_value=interval.lam,
            rho=interval.rho,
            helper=interval.script_name,
        )
    check_close("early_accelerated_75001_1048576", early_accelerated_total, -1056.187188, args.tolerance)
    print(f"early_accelerated_endpoint_75001_350000_total_log2,{endpoint_total:.6f}")
    print(f"early_accelerated_paired_350001_1048576_total_log2,{paired_total:.6f}")
    print(f"early_accelerated_75001_1048576_total_log2,{early_accelerated_total:.6f}")

    late_postprefix_cover = check_contiguous_cover(
        "late post-prefix cover",
        LATE_POSTPREFIX_INTERVALS,
        expected_start=501,
        expected_stop=args.block_bits * args.late_blocks,
    )
    for interval in LATE_POSTPREFIX_INTERVALS:
        if interval.z != LATE_POSTPREFIX_Z:
            raise SystemExit(f"late post-prefix cover: interval {interval.label} has z={interval.z}")
    print("late_postprefix_cover_status,PASS")
    print(f"late_postprefix_cover_range,{late_postprefix_cover.cover_start}--{late_postprefix_cover.cover_stop}")
    print(f"late_postprefix_interval_count,{late_postprefix_cover.count}")
    print(f"late_postprefix_feasible_cutoff_h,{args.block_bits * args.late_blocks}")
    print(f"late_postprefix_z,{LATE_POSTPREFIX_Z:.17g}")
    manifest["checks"]["late_postprefix_cover"] = {  # type: ignore[index]
        "status": "PASS",
        "cover_range": [late_postprefix_cover.cover_start, late_postprefix_cover.cover_stop],
        "interval_count": late_postprefix_cover.count,
        "feasible_cutoff_h": args.block_bits * args.late_blocks,
        "cutoff_reason": "block_bits * late_blocks",
        "z": LATE_POSTPREFIX_Z,
        "intervals": [
            {
                "range": [interval.start, interval.stop],
                "z": interval.z,
                "log2": interval.log2_value,
                "command": interval.display_command(),
            }
            for interval in LATE_POSTPREFIX_INTERVALS
        ],
    }
    late_shape = check_late_postprefix_shapes(
        LATE_POSTPREFIX_INTERVALS,
        n=args.N,
        m=args.block_bits * args.late_blocks,
    )
    print("late_postprefix_shape_status,PASS")
    print(f"late_postprefix_shape_interval_count,{late_shape.interval_count}")
    print(f"late_postprefix_shape_peak_h_range,{late_shape.peak_h_min}--{late_shape.peak_h_max}")
    print(f"late_postprefix_shape_left_endpoint_peaks,{late_shape.left_endpoint_peaks}")
    print(f"late_postprefix_shape_right_endpoint_peaks,{late_shape.right_endpoint_peaks}")
    print(f"late_postprefix_shape_critical_peaks,{late_shape.critical_peaks}")
    manifest["checks"]["late_postprefix_shape"] = {  # type: ignore[index]
        "status": "PASS",
        "interval_count": late_shape.interval_count,
        "peak_h_range": [late_shape.peak_h_min, late_shape.peak_h_max],
        "left_endpoint_peaks": late_shape.left_endpoint_peaks,
        "right_endpoint_peaks": late_shape.right_endpoint_peaks,
        "critical_peaks": late_shape.critical_peaks,
        "shape_rule": "adjacent ratio ((m-h)/(N-h))/z is decreasing, so each row is checked at endpoints and the critical point",
        "intervals": late_shape.rows,
    }

    late_postprefix_total = float("-inf")
    for interval in LATE_POSTPREFIX_INTERVALS:
        late_postprefix_total = log2add(late_postprefix_total, interval.log2_value)
        print(f"late_postprefix_interval_{interval.label}_z,{interval.z:.17g}")
        print(f"late_postprefix_interval_{interval.label}_log2,{interval.log2_value:.6f}")
        add_row(
            manifest,
            name=f"late_postprefix_interval_{interval.label}",
            role="active_postprefix_ultralate_placement",
            h_range=(interval.start, interval.stop),
            log2_value=interval.log2_value,
            command=interval.display_command(),
            z=interval.z,
        )
    check_close("late_postprefix_501_380736", late_postprefix_total, -545.690526, args.tolerance)
    print(f"late_postprefix_501_380736_total_log2,{late_postprefix_total:.6f}")

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

    if args.check_early_postprefix_intervals:
        for interval in selected_early_postprefix_intervals(args.check_early_postprefix_intervals):
            actual = check_early_postprefix_interval(interval, args.tolerance)
            print(f"early_postprefix_interval_{interval.label}_recomputed_log2,{actual:.6f}")
            print(f"early_postprefix_interval_{interval.label}_recompute_status,PASS")

    if args.check_early_accelerated_intervals:
        for interval in selected_early_accelerated_intervals(args.check_early_accelerated_intervals):
            actual = check_early_accelerated_interval(interval, args.tolerance)
            print(f"early_accelerated_interval_{interval.label}_recomputed_log2,{actual:.6f}")
            print(f"early_accelerated_interval_{interval.label}_recompute_status,PASS")

    if args.check_late_postprefix_intervals:
        for interval in selected_late_postprefix_intervals(args.check_late_postprefix_intervals):
            actual = check_late_postprefix_interval(interval, args.tolerance)
            print(f"late_postprefix_interval_{interval.label}_recomputed_log2,{actual:.6f}")
            print(f"late_postprefix_interval_{interval.label}_recompute_status,PASS")

    h501_plus_checked = float("-inf")
    for value in [
        parts["postprefix_501_2000_eall"],
        late_postprefix_total,
        early_postprefix_total,
        early_accelerated_total,
        high_total,
        complement_high_total,
    ]:
        h501_plus_checked = log2add(h501_plus_checked, value)
    check_close("h501_plus_checked_rows", h501_plus_checked, -182.259739, args.tolerance)
    print(f"h501_plus_checked_rows_log2,{h501_plus_checked:.6f}")

    total = float("-inf")
    for value in [prefix_all, early_all, h501_plus_checked]:
        total = log2add(total, value)
    print(f"finite_ledger_total_log2,{total:.6f}")
    print(f"finite_ledger_margin_bits,{-total:.6f}")
    late_plus_window = log2add(exact_late_prefix.total_log2, total)
    print(f"late_plus_window_total_log2,{late_plus_window:.6f}")
    print(f"late_plus_window_margin_bits,{-late_plus_window:.6f}")
    h32_500_all_positions = log2add(exact_late_prefix.total_log2, log2add(prefix_all, early_all))
    print(f"h32_500_all_first_active_positions_log2,{h32_500_all_positions:.6f}")
    print(f"h32_500_all_first_active_positions_margin_bits,{-h32_500_all_positions:.6f}")
    current_checked = log2add(h32_500_all_positions, h501_plus_checked)
    check_close("current_checked_ledger", current_checked, -37.278528, args.tolerance)
    print(f"current_checked_ledger_total_log2,{current_checked:.6f}")
    print(f"current_checked_ledger_margin_bits,{-current_checked:.6f}")
    manifest["totals"] = {  # type: ignore[index]
        "prefix_32_500_all_e_exact_outer_log2": prefix_all,
        "early_32_500_all_e_exact_outer_log2": early_all,
        "interval_2001_1148736_total_log2": high_total,
        "complement_interval_1048577_2097152_total_log2": complement_high_total,
        "early_postprefix_501_75000_total_log2": early_postprefix_total,
        "early_accelerated_endpoint_75001_350000_total_log2": endpoint_total,
        "early_accelerated_paired_350001_1048576_total_log2": paired_total,
        "early_accelerated_75001_1048576_total_log2": early_accelerated_total,
        "late_postprefix_501_380736_total_log2": late_postprefix_total,
        "h501_plus_checked_rows_log2": h501_plus_checked,
        "finite_ledger_total_log2": total,
        "finite_ledger_margin_bits": -total,
        "late_plus_window_total_log2": late_plus_window,
        "late_plus_window_margin_bits": -late_plus_window,
        "h32_500_all_first_active_positions_log2": h32_500_all_positions,
        "h32_500_all_first_active_positions_margin_bits": -h32_500_all_positions,
        "current_checked_ledger_total_log2": current_checked,
        "current_checked_ledger_margin_bits": -current_checked,
    }
    if args.write_manifest_json:
        output_manifest_path = args.write_manifest_json
        if not output_manifest_path.is_absolute():
            output_manifest_path = (Path.cwd() / output_manifest_path).resolve()
        output_manifest_path.parent.mkdir(parents=True, exist_ok=True)
        output_manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(f"manifest_json,{output_manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
