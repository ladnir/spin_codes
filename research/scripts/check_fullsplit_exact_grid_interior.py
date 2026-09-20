#!/usr/bin/env python3
"""Interior-T audit for the full-split exact tiny-prefix grids.

For a fixed remaining input weight H, the exact e<=E episode formula can be
swept over T by maintaining short cumulative sums for each selected-gap count.
This helper compares every T in the current placement buckets against the saved
left-endpoint CSV total.

The implementation is intended as a finite audit tool for the tiny prefix.  It
does not replace the cleaner analytic monotonicity lemma we ultimately want.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from bound_fullsplit_episode_gaps import (
    load_spectrum,
    precompute_survival_log2,
    split_entries,
)
from dense_largek_eval import log2_binom, log2add
from fast_fullsplit_episode_e01 import precompute_nonempty_block_logcoeff
from probe_block_recursive_inner import parse_int_list


ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Bucket:
    name: str
    t_min: int
    t_max: int
    left_csv: Path


WINDOW_BUCKETS = [
    Bucket("gap_1_4000", 5949, 9948, ROOT / "fast_fullsplit_e08_T5949_H0_499_all.csv"),
    Bucket("gap_4001_8000", 9949, 13948, ROOT / "fast_fullsplit_e08_T9949_H0_499_all.csv"),
    Bucket("gap_8001_12000", 13949, 17948, ROOT / "fast_fullsplit_e08_T13949_H0_499_all.csv"),
]

EARLY_BUCKETS = [
    Bucket("gap_12001_17000", 17949, 22948, ROOT / "fast_fullsplit_e08_T17949_H0_499_all.csv"),
    Bucket("gap_17001_22000", 22949, 27948, ROOT / "fast_fullsplit_e08_T22949_H0_499_all.csv"),
    Bucket("gap_22001_26819", 27949, 32767, ROOT / "fast_fullsplit_e08_T27949_H0_499_all.csv"),
]

DEFAULT_BUCKETS = WINDOW_BUCKETS


def select_buckets(bucket_set: str) -> list[Bucket]:
    if bucket_set == "window":
        return list(WINDOW_BUCKETS)
    if bucket_set == "early":
        return list(EARLY_BUCKETS)
    if bucket_set == "window-plus-early":
        return list(WINDOW_BUCKETS) + list(EARLY_BUCKETS)
    raise ValueError(f"unknown bucket set {bucket_set!r}")


def logadd2_array(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    out = np.logaddexp(a * math.log(2.0), b * math.log(2.0)) / math.log(2.0)
    return out


def logsumexp2_axis0(mat: np.ndarray) -> np.ndarray:
    max_vals = np.max(mat, axis=0)
    out = np.full(max_vals.shape, float("-inf"), dtype=np.float64)
    finite = np.isfinite(max_vals)
    if np.any(finite):
        shifted = np.exp2(mat[:, finite] - max_vals[finite])
        out[finite] = max_vals[finite] + np.log2(np.sum(shifted, axis=0))
    return out


def parse_log2(text: str) -> float:
    if text == "-inf":
        return float("-inf")
    return float(text)


def load_left_totals(path: Path) -> dict[int, float]:
    with path.open(newline="") as f:
        return {int(row["H"]): parse_log2(row["total_log2"]) for row in csv.DictReader(f)}


def precompute_logc_table(*, n_max: int, k_max: int) -> np.ndarray:
    out = np.full((k_max + 1, n_max + 1), float("-inf"), dtype=np.float64)
    out[0, :] = 0.0
    for k in range(1, k_max + 1):
        out[k, k] = 0.0
        for n in range(k + 1, n_max + 1):
            out[k, n] = out[k, n - 1] + math.log2(n / (n - k))
    return out


def sweep_h(
    *,
    H: int,
    b: int,
    e_max: int,
    t_max: int,
    turnoff_log2: float,
    survival: list[float],
    logc: np.ndarray,
    log_coeff_by_x: list[float],
    watch_t: set[int],
) -> dict[int, float]:
    """Return exact e<=e_max totals for one H at requested T values."""

    states = {
        e: np.full((e + 1, H + 1), float("-inf"), dtype=np.float64)
        for e in range(1, e_max + 1)
    }
    x_min = 0 if H == 0 else max(1, (H + b - 1) // b)
    coeff = np.array(log_coeff_by_x, dtype=np.float64)
    out: dict[int, float] = {}

    for T in range(1, t_max + 1):
        max_x = min(H, T)
        for e in range(1, min(e_max, max_x) + 1):
            xs = np.arange(e, max_x + 1, dtype=np.int64)
            bvals = logc[xs - e, T - e] + survival[T]
            st = states[e]
            st[1, xs] = logadd2_array(st[1, xs], bvals)
            for order in range(2, e + 1):
                st[order, xs] = logadd2_array(st[order, xs], st[order - 1, xs])

        if T not in watch_t:
            continue

        denom = log2_binom(b * T, H)
        total = survival[T]
        for e in range(1, e_max + 1):
            # Endpoint case e=x+1, where all gaps are selected.
            x_endpoint = e - 1
            if x_min <= x_endpoint <= min(H, T):
                c = coeff[x_endpoint]
                if c != float("-inf"):
                    term = c - denom + e * turnoff_log2 + log2_binom(T, x_endpoint) + survival[x_endpoint]
                    total = log2add(total, term)

            if e > max_x:
                continue
            xs = np.arange(max(e, x_min), max_x + 1, dtype=np.int64)
            if xs.size == 0:
                continue
            suffix = np.array([log2_binom(int(x) + 1, e) for x in xs]) + states[e][e, xs]
            vals = coeff[xs] - denom + e * turnoff_log2 + suffix
            vals = vals[np.isfinite(vals)]
            if vals.size:
                term = float(np.logaddexp.reduce(vals * math.log(2.0)) / math.log(2.0))
                total = log2add(total, term)
        out[T] = total
    return out


def update_denom_vector(*, denom: np.ndarray, h_values: np.ndarray, n_start: int, b: int) -> None:
    """Update log2 binom(n,H) from n_start to n_start+b in place."""

    for offset in range(1, b + 1):
        n = n_start + offset
        birth = h_values == n
        denom[birth] = 0.0
        valid = (h_values > 0) & (h_values < n)
        denom[valid] += np.log2(n / (n - h_values[valid]))


def sweep_combined(
    *,
    H_values: list[int],
    buckets: list[Bucket],
    b: int,
    e_max: int,
    t_max: int,
    turnoff_log2: float,
    survival: list[float],
    logc: np.ndarray,
    log_coeff: list[list[float]],
    sample_step: int,
) -> dict[tuple[str, int], tuple[float, int, float]]:
    """Return best interior totals for all requested H values.

    The result maps (bucket_name,H) to (best_log2, T_at_best, reproduced_left).
    """

    h_max = max(H_values)
    h_indices = np.array(H_values, dtype=np.int64)
    coeff = np.array(log_coeff, dtype=np.float64)[:, h_indices]
    survival_arr = np.array(survival, dtype=np.float64)
    states = {
        e: np.full((e + 1, h_max + 1), float("-inf"), dtype=np.float64)
        for e in range(1, e_max + 1)
    }
    choose_xe = {
        e: np.array(
            [log2_binom(x + 1, e) if x >= e else float("-inf") for x in range(h_max + 1)],
            dtype=np.float64,
        )
        for e in range(1, e_max + 1)
    }
    x_mins = np.array([0 if H == 0 else max(1, (H + b - 1) // b) for H in H_values], dtype=np.int64)
    denom = np.full(len(H_values), float("-inf"), dtype=np.float64)
    denom[h_indices == 0] = 0.0

    best: dict[str, np.ndarray] = {
        bucket.name: np.full(len(H_values), float("-inf"), dtype=np.float64)
        for bucket in buckets
    }
    best_t: dict[str, np.ndarray] = {
        bucket.name: np.full(len(H_values), -1, dtype=np.int64)
        for bucket in buckets
    }
    left_repro: dict[str, np.ndarray] = {
        bucket.name: np.full(len(H_values), float("nan"), dtype=np.float64)
        for bucket in buckets
    }

    bucket_by_t: dict[int, Bucket] = {}
    for bucket in buckets:
        for T in range(bucket.t_min, bucket.t_max + 1, sample_step):
            bucket_by_t[T] = bucket
        bucket_by_t[bucket.t_min] = bucket
        bucket_by_t[bucket.t_max] = bucket

    for T in range(1, t_max + 1):
        update_denom_vector(denom=denom, h_values=h_indices, n_start=b * (T - 1), b=b)

        max_x = min(h_max, T)
        for e in range(1, min(e_max, max_x) + 1):
            xs = np.arange(e, max_x + 1, dtype=np.int64)
            bvals = logc[xs - e, T - e] + survival_arr[T]
            st = states[e]
            st[1, xs] = logadd2_array(st[1, xs], bvals)
            for order in range(2, e + 1):
                st[order, xs] = logadd2_array(st[order, xs], st[order - 1, xs])

        bucket = bucket_by_t.get(T)
        if bucket is None:
            continue

        total = np.full(len(H_values), survival_arr[T], dtype=np.float64)
        x_grid = np.arange(h_max + 1, dtype=np.int64)
        for e in range(1, e_max + 1):
            x_endpoint = e - 1
            if x_endpoint <= min(h_max, T):
                endpoint_valid = x_mins <= x_endpoint
                endpoint_vals = (
                    coeff[x_endpoint, :]
                    - denom
                    + e * turnoff_log2
                    + log2_binom(T, x_endpoint)
                    + survival_arr[x_endpoint]
                )
                endpoint_vals = np.where(endpoint_valid, endpoint_vals, float("-inf"))
                total = logadd2_array(total, endpoint_vals)

            if e > max_x:
                continue
            suffix = choose_xe[e] + states[e][e, :]
            suffix = np.where(x_grid >= e, suffix, float("-inf"))
            vals = coeff + suffix[:, None]
            reduced = logsumexp2_axis0(vals) - denom + e * turnoff_log2
            total = logadd2_array(total, reduced)

        if T == bucket.t_min:
            left_repro[bucket.name] = total.copy()
        improve = total > best[bucket.name]
        best[bucket.name][improve] = total[improve]
        best_t[bucket.name][improve] = T

    result: dict[tuple[str, int], tuple[float, int, float]] = {}
    for bucket in buckets:
        for idx, H in enumerate(H_values):
            result[(bucket.name, H)] = (
                float(best[bucket.name][idx]),
                int(best_t[bucket.name][idx]),
                float(left_repro[bucket.name][idx]),
            )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=ROOT / "EBCH128_64.wd")
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--turnoff-log2", type=float, default=-63.8926492)
    parser.add_argument("--h-values", default="0:80")
    parser.add_argument("--e-max", type=int, default=8)
    parser.add_argument("--lambda-min", type=float, default=0.001)
    parser.add_argument("--lambda-max", type=float, default=2.0)
    parser.add_argument("--lambda-step", type=float, default=0.004)
    parser.add_argument("--tolerance", type=float, default=1e-7)
    parser.add_argument("--sample-step", type=int, default=1)
    parser.add_argument("--method", choices=("combined", "single"), default="combined")
    parser.add_argument(
        "--bucket-set",
        choices=("window", "early", "window-plus-early"),
        default="window",
        help="Which saved left-endpoint bucket family to audit.",
    )
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    H_values = parse_int_list(args.h_values)
    h_max = max(H_values)
    buckets = select_buckets(args.bucket_set)
    t_max = max(bucket.t_max for bucket in buckets)
    d = math.floor(args.distance_delta * args.N)
    entries = split_entries(load_spectrum(args.spectrum), args.block_bits)
    survival = precompute_survival_log2(
        entries=entries,
        max_live=t_max,
        distance=d,
        lambda_min=args.lambda_min,
        lambda_max=args.lambda_max,
        lambda_step=args.lambda_step,
    )
    logc = precompute_logc_table(n_max=t_max, k_max=h_max)
    log_coeff = precompute_nonempty_block_logcoeff(b=args.block_bits, h_max=h_max, x_max=h_max)
    left_by_bucket = {bucket.name: load_left_totals(bucket.left_csv) for bucket in buckets}

    fields = [
        "bucket",
        "H",
        "T_count",
        "max_log2_minus_left",
        "T_at_max",
        "left_log2",
        "max_log2",
        "left_reproduction_delta",
    ]
    print("Full-split exact-grid interior audit")
    if not args.summary_only:
        print(",".join(fields))
    output_rows: list[dict[str, str]] = []
    failed = False
    if args.method == "combined":
        combined = sweep_combined(
            H_values=H_values,
            buckets=buckets,
            b=args.block_bits,
            e_max=args.e_max,
            t_max=t_max,
            turnoff_log2=args.turnoff_log2,
            survival=survival,
            logc=logc,
            log_coeff=log_coeff,
            sample_step=args.sample_step,
        )
        for H in H_values:
            for bucket in buckets:
                left = left_by_bucket[bucket.name][H]
                best, best_t, left_repro = combined[(bucket.name, H)]
                diff = best - left
                repro_delta = left_repro - left
                count = len(range(bucket.t_min, bucket.t_max + 1, args.sample_step))
                if (bucket.t_max - bucket.t_min) % args.sample_step:
                    count += 1
                if diff > args.tolerance or abs(repro_delta) > 5e-6:
                    failed = True
                row = {
                    "bucket": bucket.name,
                    "H": str(H),
                    "T_count": str(count),
                    "max_log2_minus_left": f"{diff:.12g}",
                    "T_at_max": str(best_t),
                    "left_log2": f"{left:.12g}",
                    "max_log2": f"{best:.12g}",
                    "left_reproduction_delta": f"{repro_delta:.12g}",
                }
                output_rows.append(row)
                if not args.summary_only:
                    print(",".join(row[field] for field in fields))
        if args.output_csv:
            with args.output_csv.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerows(output_rows)
        worst = max(output_rows, key=lambda row: float(row["max_log2_minus_left"]))
        worst_repro = max(output_rows, key=lambda row: abs(float(row["left_reproduction_delta"])))
        print(
            "summary,"
            f"rows={len(output_rows)},"
            f"worst_bucket={worst['bucket']},"
            f"worst_H={worst['H']},"
            f"worst_diff={worst['max_log2_minus_left']},"
            f"worst_T={worst['T_at_max']},"
            f"worst_repro_bucket={worst_repro['bucket']},"
            f"worst_repro_H={worst_repro['H']},"
            f"worst_repro_delta={worst_repro['left_reproduction_delta']}"
        )
        if failed:
            raise SystemExit("interior audit failed")
        return 0

    for H in H_values:
        watch_t: set[int] = set()
        for bucket in buckets:
            watch_t.update(range(bucket.t_min, bucket.t_max + 1, args.sample_step))
            watch_t.add(bucket.t_min)
            watch_t.add(bucket.t_max)
        totals = sweep_h(
            H=H,
            b=args.block_bits,
            e_max=args.e_max,
            t_max=t_max,
            turnoff_log2=args.turnoff_log2,
            survival=survival,
            logc=logc,
            log_coeff_by_x=[log_coeff[x][H] for x in range(h_max + 1)],
            watch_t=watch_t,
        )
        for bucket in buckets:
            left = left_by_bucket[bucket.name][H]
            left_repro = totals[bucket.t_min]
            best_t = bucket.t_min
            best = float("-inf")
            count = 0
            for T in range(bucket.t_min, bucket.t_max + 1, args.sample_step):
                value = totals[T]
                count += 1
                if value > best:
                    best = value
                    best_t = T
            # Ensure the right endpoint is included when sample_step does not land on it.
            if (bucket.t_max - bucket.t_min) % args.sample_step:
                count += 1
                value = totals[bucket.t_max]
                if value > best:
                    best = value
                    best_t = bucket.t_max
            diff = best - left
            repro_delta = left_repro - left
            if diff > args.tolerance or abs(repro_delta) > 5e-6:
                failed = True
            row = {
                "bucket": bucket.name,
                "H": str(H),
                "T_count": str(count),
                "max_log2_minus_left": f"{diff:.12g}",
                "T_at_max": str(best_t),
                "left_log2": f"{left:.12g}",
                "max_log2": f"{best:.12g}",
                "left_reproduction_delta": f"{repro_delta:.12g}",
            }
            output_rows.append(row)
            if not args.summary_only:
                print(",".join(row[field] for field in fields))

    if args.output_csv:
        with args.output_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(output_rows)
    if args.summary_only:
        worst = max(output_rows, key=lambda row: float(row["max_log2_minus_left"]))
        worst_repro = max(output_rows, key=lambda row: abs(float(row["left_reproduction_delta"])))
        print(
            "summary,"
            f"rows={len(output_rows)},"
            f"worst_bucket={worst['bucket']},"
            f"worst_H={worst['H']},"
            f"worst_diff={worst['max_log2_minus_left']},"
            f"worst_T={worst['T_at_max']},"
            f"worst_repro_bucket={worst_repro['bucket']},"
            f"worst_repro_H={worst_repro['H']},"
            f"worst_repro_delta={worst_repro['left_reproduction_delta']}"
        )
    if failed:
        raise SystemExit("interior audit failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
