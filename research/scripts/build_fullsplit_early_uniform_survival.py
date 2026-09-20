#!/usr/bin/env python3
"""Build bucket-uniform early-prefix inner CSVs for the full-split certificate.

For the early first-active buckets the exact left-endpoint grids are strong, but
auditing every interior T is expensive.  This helper builds a direct bucket
upper bound for all T in [T_min,T_max].

For fixed remaining input weight H, occupied remaining blocks x, and selected
termination gaps e, it uses

  C(b T, H) >= C(b T_min, H)

and, after writing live = T - skipped,

  C(T-live+e-1,e-1) C(live-e,x-e)
    <= C(T_max-live+e-1,e-1) C(live-e,x-e).

The Chernoff survival factor is kept as survival(live).  Thus each bucket gets
one H-indexed CSV that upper-bounds the exact e<=E episode ledger throughout the
bucket, without assuming endpoint monotonicity.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

from bound_fullsplit_episode_gaps import (
    load_spectrum,
    precompute_survival_log2,
    split_entries,
)
from dense_largek_eval import log2_binom, log2add
from fast_fullsplit_episode_e01 import precompute_nonempty_block_logcoeff


ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Bucket:
    name: str
    t_min: int
    t_max: int
    output_name: str


EARLY_BUCKETS = [
    Bucket(
        "gap_12001_17000",
        17949,
        22948,
        "fullsplit_early_{e_tag}_uniformsurv_T17949_22948_H0_499.csv",
    ),
    Bucket(
        "gap_17001_22000",
        22949,
        27948,
        "fullsplit_early_{e_tag}_uniformsurv_T22949_27948_H0_499.csv",
    ),
    Bucket(
        "gap_22001_26819",
        27949,
        32767,
        "fullsplit_early_{e_tag}_uniformsurv_T27949_32767_H0_499.csv",
    ),
]


def precompute_logc_table(*, n_max: int, k_max: int) -> list[list[float]]:
    table: list[list[float]] = []
    for k in range(k_max + 1):
        row = [float("-inf")] * (n_max + 1)
        if k == 0:
            for n in range(n_max + 1):
                row[n] = 0.0
        elif k <= n_max:
            row[k] = 0.0
            for n in range(k + 1, n_max + 1):
                row[n] = row[n - 1] + math.log2(n / (n - k))
        table.append(row)
    return table


def check_survival_monotone(survival: list[float], tolerance: float = 1e-10) -> None:
    for live in range(1, len(survival)):
        if survival[live] > survival[live - 1] + tolerance:
            raise ValueError(
                f"survival is not nonincreasing at live={live}: "
                f"{survival[live - 1]} -> {survival[live]}"
            )


def uniform_suffix_bounds(
    *,
    bucket: Bucket,
    h_max: int,
    e_max: int,
    survival: list[float],
    logc: list[list[float]],
) -> list[list[float]]:
    suffix = [[float("-inf")] * (h_max + 1) for _ in range(e_max + 1)]
    t_max = bucket.t_max
    for e in range(1, e_max + 1):
        choose_gap_count = logc[e]
        left_row = logc[e - 1]
        for x in range(h_max + 1):
            if e > x + 1:
                continue
            if e == x + 1:
                suffix[e][x] = logc[x][t_max] + survival[x]
                continue

            right_row = logc[x - e]
            total = float("-inf")
            for live in range(max(x, e), t_max + 1):
                left = left_row[t_max - live + e - 1]
                right = right_row[live - e]
                if left == float("-inf") or right == float("-inf"):
                    continue
                total = log2add(total, left + right + survival[live])
            suffix[e][x] = choose_gap_count[x + 1] + total
    return suffix


def build_rows(
    *,
    bucket: Bucket,
    b: int,
    h_max: int,
    e_max: int,
    turnoff_log2: float,
    survival: list[float],
    logc: list[list[float]],
    log_coeff: list[list[float]],
) -> list[dict[str, float | int]]:
    suffix = uniform_suffix_bounds(
        bucket=bucket,
        h_max=h_max,
        e_max=e_max,
        survival=survival,
        logc=logc,
    )
    rows: list[dict[str, float | int]] = []
    for H in range(h_max + 1):
        denom = log2_binom(b * bucket.t_min, H)
        total = survival[bucket.t_min]
        x_min = 0 if H == 0 else max(1, (H + b - 1) // b)
        x_max = min(H, bucket.t_max)
        for x in range(x_min, x_max + 1):
            coeff = log_coeff[x][H] if x <= h_max else float("-inf")
            if coeff == float("-inf"):
                continue
            base = coeff - denom
            for e in range(1, min(e_max, x + 1) + 1):
                total = log2add(total, base + e * turnoff_log2 + suffix[e][x])
        rows.append({"H": H, "total_log2": min(0.0, total)})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=ROOT / "EBCH128_64.wd")
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--h-max", type=int, default=499)
    parser.add_argument("--e-max", type=int, default=8)
    parser.add_argument("--turnoff-log2", type=float, default=-63.8926492)
    parser.add_argument("--lambda-min", type=float, default=0.001)
    parser.add_argument("--lambda-max", type=float, default=2.0)
    parser.add_argument("--lambda-step", type=float, default=0.004)
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    args = parser.parse_args()

    d = math.floor(args.distance_delta * args.N)
    buckets = EARLY_BUCKETS
    max_live = max(bucket.t_max for bucket in buckets)
    entries = split_entries(load_spectrum(args.spectrum), args.block_bits)
    survival = precompute_survival_log2(
        entries=entries,
        max_live=max_live,
        distance=d,
        lambda_min=args.lambda_min,
        lambda_max=args.lambda_max,
        lambda_step=args.lambda_step,
    )
    check_survival_monotone(survival)
    logc = precompute_logc_table(n_max=max_live + args.e_max, k_max=args.h_max + args.e_max + 1)
    log_coeff = precompute_nonempty_block_logcoeff(
        b=args.block_bits,
        h_max=args.h_max,
        x_max=args.h_max,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print("Full-split early uniform-survival builder")
    print(f"N,{args.N}")
    print(f"delta,{args.distance_delta}")
    print(f"d,{d}")
    print(f"block_bits,{args.block_bits}")
    print(f"h_max,{args.h_max}")
    print(f"e_max,{args.e_max}")
    e_tag = f"e{args.e_max:02d}"
    for bucket in buckets:
        rows = build_rows(
            bucket=bucket,
            b=args.block_bits,
            h_max=args.h_max,
            e_max=args.e_max,
            turnoff_log2=args.turnoff_log2,
            survival=survival,
            logc=logc,
            log_coeff=log_coeff,
        )
        out_path = args.output_dir / bucket.output_name.format(e_tag=e_tag)
        with out_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["H", "total_log2"])
            writer.writeheader()
            writer.writerows(rows)
        sample_low = min(31, args.h_max)
        sample_high = min(499, args.h_max)
        print(f"{bucket.name}_csv,{out_path}")
        print(f"{bucket.name}_H{sample_low}_log2,{rows[sample_low]['total_log2']:.6f}")
        print(f"{bucket.name}_H{sample_high}_log2,{rows[sample_high]['total_log2']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
