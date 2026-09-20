#!/usr/bin/env python3
"""Compare outer block models against the current full-split EBCH inner ledger.

This is a projection/audit utility, not a theorem upgrader.  It keeps the
current inner full-split ledger fixed and swaps only the local outer spectrum
model.  The RM row uses the exact local RM spectrum CSV.  The BCH-like rows use
the random-like envelope described in the report and are labelled heuristic.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent

N_DEFAULT = 2**21
K_DEFAULT = 2**20
INNER_B_DEFAULT = 64
LATE_BLOCKS_DEFAULT = 5949
H_PREFIX_MAX_DEFAULT = 500
H_MAX_DEFAULT = 2000
CHECKED_BASELINE_LOG2 = -37.278528
RM_GT2000_EXISTING_THEOREM_LOG2 = -586.597103


def log2add(a: float, b: float) -> float:
    if a == float("-inf"):
        return b
    if b == float("-inf"):
        return a
    if a < b:
        a, b = b, a
    return a + math.log2(1.0 + 2.0 ** (b - a))


def log2sum(values: Iterable[float]) -> float:
    total = float("-inf")
    for value in values:
        total = log2add(total, value)
    return total


def log2_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2.0)


def fmt_log2(value: float) -> str:
    if value == float("-inf"):
        return "-inf"
    if math.isnan(value):
        return "NA"
    return f"{value:.6f}"


def fmt_margin(value: float) -> str:
    if value == float("-inf"):
        return "inf"
    if math.isnan(value):
        return "NA"
    return f"{-value:.6f}"


def parse_float(text: str) -> float:
    if text.lower() == "nan":
        return math.nan
    return float(text)


@dataclass(frozen=True)
class OuterModeSpec:
    name: str
    status: str
    local_name: str
    local_length: int
    local_dimension: int
    local_distance: int
    blocks: int
    spectrum_kind: str
    source: str
    inflation_bits: int = 0

    @property
    def family(self) -> str:
        return self.name

    @property
    def display_name(self) -> str:
        if self.spectrum_kind == "exact_csv":
            return self.name
        return f"{self.name}_plus{self.inflation_bits}"

    @property
    def cost_proxy(self) -> float:
        return self.blocks * self.local_length * math.log2(self.local_length)


@dataclass(frozen=True)
class Term:
    component: str
    h: int
    log2_value: float
    first_r: str = ""
    gap_min: str = ""
    gap_max: str = ""
    bucket_T: str = ""
    inner_mode: str = ""
    source: str = ""

    @property
    def regime(self) -> str:
        if self.component == "ultra_late_T_lt_5949":
            return "ultra_late_prefix"
        first = f"r={self.first_r}" if self.first_r else "r=NA"
        gap = f"gap={self.gap_min}..{self.gap_max}" if self.gap_min else "gap=NA"
        return f"{self.component}; {first}; {gap}"


@dataclass(frozen=True)
class ModeResult:
    spec: OuterModeSpec
    total_log2: float
    h_lt32_log2: float
    h_32_500_log2: float
    h_501_2000_log2: float
    h_gt2000_log2: float
    h_gt2000_status: str
    dominant: Term
    local_min_positive: int
    local_support_le_hmax: int


def direct_sum_by_active_blocks(
    local: list[float],
    *,
    blocks: int,
    h_max: int,
    min_positive: int,
) -> list[float]:
    """Return log2 coefficients of (1 + L_+(z))**blocks through h_max.

    The direct sum only needs at most ``h_max // min_positive`` active local
    blocks.  That is far smaller than binary powering for the present
    low-weight projection, especially when we evaluate several sensitivity
    lanes.
    """
    local_items = [(i, value) for i, value in enumerate(local) if i > 0 and value != float("-inf")]
    if not local_items:
        raise ValueError("local spectrum has no positive terms")
    result = [float("-inf")] * (h_max + 1)
    result[0] = 0.0
    power = [float("-inf")] * (h_max + 1)
    power[0] = 0.0
    max_active = min(blocks, h_max // min_positive)
    local_max = local_items[-1][0]
    for active in range(1, max_active + 1):
        new_power = [float("-inf")] * (h_max + 1)
        prev_min = (active - 1) * min_positive
        prev_max = min(h_max, (active - 1) * local_max)
        for prev_h in range(prev_min, prev_max + 1):
            prev_value = power[prev_h]
            if prev_value == float("-inf"):
                continue
            limit = h_max - prev_h
            for local_h, local_value in local_items:
                if local_h > limit:
                    break
                target = prev_h + local_h
                new_power[target] = log2add(new_power[target], prev_value + local_value)
        power = new_power
        choose = log2_binom(blocks, active)
        for h, value in enumerate(power):
            if value != float("-inf"):
                result[h] = log2add(result[h], choose + value)
    return result


def exact_local_log_spectrum(path: Path, h_max: int) -> list[float]:
    coeffs = [float("-inf")] * (h_max + 1)
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            weight = int(row["weight"])
            if weight <= h_max:
                count = int(row["count"])
                if count > 0:
                    coeffs[weight] = math.log2(count)
    if coeffs[0] == float("-inf"):
        coeffs[0] = 0.0
    return coeffs


def random_like_local_log_spectrum(spec: OuterModeSpec, h_max: int) -> list[float]:
    coeffs = [float("-inf")] * (h_max + 1)
    coeffs[0] = 0.0
    for w in range(spec.local_distance, min(spec.local_length, h_max) + 1):
        value = spec.local_dimension - spec.local_length + log2_binom(spec.local_length, w)
        if w <= spec.local_distance + 32:
            value += spec.inflation_bits
        coeffs[w] = min(value, log2_binom(spec.local_length, w))
    return coeffs


def local_log_spectrum(spec: OuterModeSpec, h_max: int) -> list[float]:
    if spec.spectrum_kind == "exact_csv":
        return exact_local_log_spectrum(ROOT / spec.source, h_max)
    if spec.spectrum_kind == "random_like_floor":
        return random_like_local_log_spectrum(spec, h_max)
    raise ValueError(f"unknown spectrum kind: {spec.spectrum_kind}")


def positive_support(coeffs: list[float]) -> list[int]:
    return [i for i, value in enumerate(coeffs) if i > 0 and value != float("-inf")]


def old_outer_column(row: dict[str, str]) -> str:
    if "outer_log2_bound" in row:
        return "outer_log2_bound"
    if "outer_log2" in row:
        return "outer_log2"
    if "peak_outer_log2_bound" in row:
        return "peak_outer_log2_bound"
    raise KeyError("ledger row has no outer log column")


def add_term(
    terms: list[Term],
    *,
    component: str,
    h: int,
    log2_value: float,
    row: dict[str, str] | None = None,
    source: str,
) -> None:
    if log2_value == float("-inf"):
        return
    row = row or {}
    terms.append(
        Term(
            component=component,
            h=h,
            log2_value=log2_value,
            first_r=row.get("first_r", ""),
            gap_min=row.get("gap_min", ""),
            gap_max=row.get("gap_max", ""),
            bucket_T=row.get("bucket_T", ""),
            inner_mode=row.get("inner_mode", ""),
            source=source,
        )
    )


def add_ultra_late_terms(
    terms: list[Term],
    outer_log: list[float],
    *,
    n: int,
    b: int,
    late_blocks: int,
    h_prefix_max: int,
) -> None:
    for h in range(1, min(h_prefix_max, len(outer_log) - 1) + 1):
        outer = outer_log[h]
        if outer == float("-inf"):
            continue
        late = log2_binom(late_blocks * b, h) - log2_binom(n, h)
        add_term(
            terms,
            component=f"ultra_late_T_lt_{late_blocks}",
            h=h,
            log2_value=outer + late,
            source="exact late-placement formula",
        )


def add_row_csv_terms(
    terms: list[Term],
    outer_log: list[float],
    *,
    path: Path,
    component: str,
) -> None:
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            h = int(row["outer_weight"])
            if h >= len(outer_log):
                continue
            outer = outer_log[h]
            if outer == float("-inf"):
                continue
            old_outer = parse_float(row[old_outer_column(row)])
            raw_term = parse_float(row["term_log2"])
            term = raw_term - old_outer + outer
            add_term(terms, component=component, h=h, log2_value=term, row=row, source=path.name)


def add_hsummary_terms(
    terms: list[Term],
    outer_log: list[float],
    *,
    path: Path,
    component: str,
) -> None:
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            h = int(row["outer_weight"])
            if h >= len(outer_log):
                continue
            outer = outer_log[h]
            if outer == float("-inf"):
                continue
            old_outer = parse_float(row["peak_outer_log2_bound"])
            aggregate = parse_float(row["aggregate_log2"])
            term = aggregate - old_outer + outer
            add_term(terms, component=component, h=h, log2_value=term, row=row, source=path.name)


def range_logsum(terms: list[Term], start: int, stop: int) -> float:
    return log2sum(term.log2_value for term in terms if start <= term.h <= stop)


def evaluate_mode(spec: OuterModeSpec, args: argparse.Namespace) -> ModeResult:
    if spec.blocks * spec.local_dimension != args.K:
        raise SystemExit(
            f"{spec.display_name}: dimension mismatch, "
            f"{spec.blocks} * {spec.local_dimension} != {args.K}"
        )
    if spec.blocks * spec.local_length != args.N:
        raise SystemExit(
            f"{spec.display_name}: length mismatch, "
            f"{spec.blocks} * {spec.local_length} != {args.N}"
        )
    local = local_log_spectrum(spec, args.h_max)
    support = positive_support(local)
    if not support:
        raise SystemExit(f"{spec.display_name}: local model has no positive support")
    local_min_positive = support[0]
    if local_min_positive < spec.local_distance:
        raise SystemExit(
            f"{spec.display_name}: local model has mass below distance floor "
            f"{spec.local_distance}: first positive {local_min_positive}"
        )
    if spec.status == "heuristic" and "heuristic" not in spec.name:
        raise SystemExit(f"{spec.display_name}: heuristic row is not labelled heuristic")

    outer_log = direct_sum_by_active_blocks(
        local,
        blocks=spec.blocks,
        h_max=args.h_max,
        min_positive=local_min_positive,
    )
    terms: list[Term] = []
    add_ultra_late_terms(
        terms,
        outer_log,
        n=args.N,
        b=args.inner_b,
        late_blocks=args.late_blocks,
        h_prefix_max=args.h_prefix_max,
    )
    add_row_csv_terms(
        terms,
        outer_log,
        path=ROOT / "fullsplit_piecewise_h32_500_csv.csv",
        component="window_e_le8",
    )
    add_row_csv_terms(
        terms,
        outer_log,
        path=ROOT / "fullsplit_turnoff_tail_h32_500_r1_64_emin9.csv",
        component="window_e_ge9_tail",
    )
    add_row_csv_terms(
        terms,
        outer_log,
        path=ROOT / "fullsplit_piecewise_early_h32_500_e16_uniformsurv.csv",
        component="early_e_le16_uniformsurv",
    )
    add_row_csv_terms(
        terms,
        outer_log,
        path=ROOT / "fullsplit_turnoff_tail_early_h32_500_r1_64_emin17.csv",
        component="early_e_ge17_tail",
    )
    add_hsummary_terms(
        terms,
        outer_log,
        path=ROOT / "fullsplit_piecewise_h501_2000_eall_hsummary.csv",
        component="postprefix_501_2000_eall_projection",
    )
    if not terms:
        raise SystemExit(f"{spec.display_name}: no live ledger terms")
    dominant = max(terms, key=lambda term: term.log2_value)
    gt2000_log2 = (
        RM_GT2000_EXISTING_THEOREM_LOG2
        if spec.name == "rm512_exact"
        else math.nan
    )
    total_log2 = log2sum(term.log2_value for term in terms)
    if not math.isnan(gt2000_log2):
        total_log2 = log2add(total_log2, gt2000_log2)
    return ModeResult(
        spec=spec,
        total_log2=total_log2,
        h_lt32_log2=range_logsum(terms, 1, 31),
        h_32_500_log2=range_logsum(terms, 32, 500),
        h_501_2000_log2=range_logsum(terms, 501, 2000),
        h_gt2000_log2=gt2000_log2,
        h_gt2000_status=(
            "covered by existing RM theorem; not recomputed by this projection driver"
            if spec.name == "rm512_exact"
            else "not recomputed for heuristic projection"
        ),
        dominant=dominant,
        local_min_positive=local_min_positive,
        local_support_le_hmax=len(support),
    )


def outer_specs() -> list[OuterModeSpec]:
    specs = [
        OuterModeSpec(
            name="rm512_exact",
            status="proved",
            local_name="RM(4,9) [512,256,32]",
            local_length=512,
            local_dimension=256,
            local_distance=32,
            blocks=4096,
            spectrum_kind="exact_csv",
            source="rm512_256_spectrum.csv",
        )
    ]
    for inflation in [0, 10, 20, 40]:
        specs.append(
            OuterModeSpec(
                name="bch256_heuristic",
                status="heuristic",
                local_name="BCH-like [256,128,38]",
                local_length=256,
                local_dimension=128,
                local_distance=38,
                blocks=8192,
                spectrum_kind="random_like_floor",
                source="random-like hard-floor envelope; CodeTables construction record",
                inflation_bits=inflation,
            )
        )
    for inflation in [0, 10, 20, 40]:
        specs.append(
            OuterModeSpec(
                name="bch512_heuristic",
                status="heuristic",
                local_name="extended BCH-like [512,256,>=62]",
                local_length=512,
                local_dimension=256,
                local_distance=62,
                blocks=4096,
                spectrum_kind="random_like_floor",
                source="random-like hard-floor envelope; primitive BCH designed-distance projection",
                inflation_bits=inflation,
            )
        )
    return specs


def check_monotone_inflation(results: list[ModeResult], family: str) -> None:
    rows = [result for result in results if result.spec.name == family]
    rows.sort(key=lambda result: result.spec.inflation_bits)
    for prev, cur in zip(rows, rows[1:]):
        if cur.total_log2 + 1e-9 < prev.total_log2:
            raise SystemExit(
                f"{family}: inflation monotonicity failed: "
                f"+{cur.spec.inflation_bits} total {cur.total_log2:.9f} < "
                f"+{prev.spec.inflation_bits} total {prev.total_log2:.9f}"
            )


def check_baseline(results: list[ModeResult], tolerance: float) -> None:
    rm = next(result for result in results if result.spec.name == "rm512_exact")
    if abs(rm.total_log2 - CHECKED_BASELINE_LOG2) > tolerance:
        raise SystemExit(
            "rm512_exact baseline reproduction failed: "
            f"{rm.total_log2:.9f} differs from {CHECKED_BASELINE_LOG2:.9f}"
        )


def write_csv(results: list[ModeResult], path: Path, args: argparse.Namespace) -> None:
    fieldnames = [
        "mode",
        "status",
        "local_outer",
        "local_length",
        "local_dimension",
        "local_distance_floor",
        "blocks",
        "spectrum_model",
        "sensitivity_bits",
        "source",
        "N",
        "K",
        "delta",
        "d",
        "projected_total_log2",
        "margin_bits",
        "dominant_h",
        "dominant_component",
        "dominant_first_active_regime",
        "dominant_first_r",
        "dominant_gap_min",
        "dominant_gap_max",
        "dominant_bucket_T",
        "dominant_term_log2",
        "h_lt32_log2",
        "h_32_500_log2",
        "h_501_2000_log2",
        "h_gt2000_log2",
        "h_gt2000_status",
        "local_min_positive",
        "local_support_le_hmax",
        "cost_blocks",
        "cost_local_length",
        "cost_n_log2n",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            spec = result.spec
            dom = result.dominant
            writer.writerow(
                {
                    "mode": spec.display_name,
                    "status": spec.status,
                    "local_outer": spec.local_name,
                    "local_length": spec.local_length,
                    "local_dimension": spec.local_dimension,
                    "local_distance_floor": spec.local_distance,
                    "blocks": spec.blocks,
                    "spectrum_model": spec.spectrum_kind,
                    "sensitivity_bits": spec.inflation_bits,
                    "source": spec.source,
                    "N": args.N,
                    "K": args.K,
                    "delta": f"{args.delta:.12g}",
                    "d": args.d,
                    "projected_total_log2": fmt_log2(result.total_log2),
                    "margin_bits": fmt_margin(result.total_log2),
                    "dominant_h": dom.h,
                    "dominant_component": dom.component,
                    "dominant_first_active_regime": dom.regime,
                    "dominant_first_r": dom.first_r,
                    "dominant_gap_min": dom.gap_min,
                    "dominant_gap_max": dom.gap_max,
                    "dominant_bucket_T": dom.bucket_T,
                    "dominant_term_log2": fmt_log2(dom.log2_value),
                    "h_lt32_log2": fmt_log2(result.h_lt32_log2),
                    "h_32_500_log2": fmt_log2(result.h_32_500_log2),
                    "h_501_2000_log2": fmt_log2(result.h_501_2000_log2),
                    "h_gt2000_log2": fmt_log2(result.h_gt2000_log2),
                    "h_gt2000_status": result.h_gt2000_status,
                    "local_min_positive": result.local_min_positive,
                    "local_support_le_hmax": result.local_support_le_hmax,
                    "cost_blocks": spec.blocks,
                    "cost_local_length": spec.local_length,
                    "cost_n_log2n": f"{spec.cost_proxy:.0f}",
                }
            )


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |"]
    out.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def write_markdown(results: list[ModeResult], path: Path, args: argparse.Namespace) -> None:
    baseline_rows = []
    sensitivity_rows = []
    cumulative_rows = []
    for result in results:
        spec = result.spec
        dom = result.dominant
        if spec.inflation_bits == 0 or spec.spectrum_kind == "exact_csv":
            baseline_rows.append(
                [
                    spec.display_name,
                    spec.status,
                    spec.local_name,
                    str(spec.blocks),
                    str(spec.local_distance),
                    fmt_log2(result.total_log2),
                    fmt_margin(result.total_log2),
                    str(dom.h),
                    dom.regime,
                    f"{spec.cost_proxy:.0f}",
                ]
            )
        sensitivity_rows.append(
            [
                spec.display_name,
                spec.status,
                f"+{spec.inflation_bits}",
                fmt_log2(result.total_log2),
                fmt_margin(result.total_log2),
                str(dom.h),
                fmt_log2(dom.log2_value),
            ]
        )
        cumulative_rows.append(
            [
                spec.display_name,
                fmt_log2(result.h_lt32_log2),
                fmt_log2(result.h_32_500_log2),
                fmt_log2(result.h_501_2000_log2),
                fmt_log2(result.h_gt2000_log2),
            ]
        )

    text = [
        "# Outer Mode Comparison at delta = 0.09",
        "",
        f"Generated by `scripts/compare_outer_modes_fullsplit.py` on {datetime.now().isoformat(timespec='seconds')}.",
        "",
        f"Checkpoint: `N={args.N}`, `K={args.K}`, `delta={args.delta}`, `d={args.d}`.  The inner model is the current full-split EBCH `[128,64,22]` ledger with block size `b={args.inner_b}`.",
        "",
        "Important: only `rm512_exact` is a proved/certificate row.  The BCH rows are spectrum-model projections using a hard distance floor and random-like local weight counts.  They do not change the current theorem.",
        "",
        "## Baseline Rows",
        "",
        markdown_table(
            [
                "mode",
                "status",
                "local outer",
                "blocks",
                "d0",
                "log2 mu",
                "margin bits",
                "dom h",
                "dominant first-active regime",
                "cost B n log2 n",
            ],
            baseline_rows,
        ),
        "",
        "## Low-Weight Sensitivity",
        "",
        "For BCH projections, `+s` means adding `s` bits to every local low-weight count in the window `d0 <= w <= d0+32` before direct-sum convolution.",
        "",
        markdown_table(
            ["mode", "status", "low-weight inflation", "log2 mu", "margin bits", "dom h", "dom term"],
            sensitivity_rows,
        ),
        "",
        "## Cumulative Contribution by Outer Weight",
        "",
        markdown_table(
            ["mode", "h<32", "32..500", "501..2000", ">2000"],
            cumulative_rows,
        ),
        "",
        "The `>2000` range is not isolated in this projection driver.  For the proved RM row it is already covered by the existing finite theorem and is far below the displayed mass; for BCH projections it must be recomputed with a real spectrum/envelope before becoming proof text.",
        "",
        "## Spectrum Notes",
        "",
        "- `rm512_exact`: exact RM `[512,256,32]` local spectrum from `scripts/rm512_256_spectrum.csv`.",
        "- `bch256_heuristic`: BCH-like `[256,128,38]`, motivated by the CodeTables `[256,128]` construction record; modeled as `2^(k-n) binom(n,w)` above the hard floor.",
        "- `bch512_heuristic`: extended BCH-like `[512,256,>=62]`, motivated by primitive BCH length `511`, designed distance `61`, extension, and subcode to dimension `256`; modeled the same way.",
        "",
        "Sources: [CodeTables `[256,128]`](https://www.codetables.de/BKLC/BKLC.php?k=128&n=256&q=2), [Kusaka weight-distribution index](https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/), and [BCH designed-distance bound](https://errorcorrectionzoo.org/c/bch).",
        "",
    ]
    path.write_text("\n".join(text), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=N_DEFAULT)
    parser.add_argument("--K", type=int, default=K_DEFAULT)
    parser.add_argument("--delta", type=float, default=0.09)
    parser.add_argument("--d", type=int, default=None)
    parser.add_argument("--inner-b", type=int, default=INNER_B_DEFAULT)
    parser.add_argument("--late-blocks", type=int, default=LATE_BLOCKS_DEFAULT)
    parser.add_argument("--h-prefix-max", type=int, default=H_PREFIX_MAX_DEFAULT)
    parser.add_argument("--h-max", type=int, default=H_MAX_DEFAULT)
    parser.add_argument("--baseline-tolerance", type=float, default=5e-6)
    parser.add_argument("--out-csv", type=Path, default=ROOT / "outer_mode_comparison_delta009.csv")
    parser.add_argument("--out-md", type=Path, default=ROOT / "outer_mode_comparison_delta009.md")
    parser.add_argument("--skip-checks", action="store_true")
    args = parser.parse_args()
    if args.d is None:
        args.d = math.floor(args.delta * args.N)

    results = [evaluate_mode(spec, args) for spec in outer_specs()]
    if not args.skip_checks:
        check_baseline(results, args.baseline_tolerance)
        check_monotone_inflation(results, "bch256_heuristic")
        check_monotone_inflation(results, "bch512_heuristic")

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    write_csv(results, args.out_csv, args)
    write_markdown(results, args.out_md, args)

    print(f"wrote_csv,{args.out_csv}")
    print(f"wrote_markdown,{args.out_md}")
    for result in results:
        print(
            f"{result.spec.display_name},status={result.spec.status},"
            f"log2_mu={result.total_log2:.6f},margin_bits={-result.total_log2:.6f},"
            f"dominant_h={result.dominant.h},dominant={result.dominant.regime}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
