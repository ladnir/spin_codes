#!/usr/bin/env python3
"""Piecewise first-moment certificate for the full-split inner.

This is the theorem-interface companion to the finite-grid diagnostics.  The
outer code is supplied by a local spectrum generating-function bound, while the
inner contribution is selected per first-active-position bucket:

  cap     : use the trivial inner bound 1;
  cauchy  : use the coefficient-free all-episode envelope from innerDense.tex;
  e01cauchy: use a sharper coefficient-free envelope for e=0 and e=1 only;
  eallratio: use the e=0,1 envelope plus the fixed-pole e>=2 ratio lemma;
  csv     : use an audited inner CSV indexed by H.

The point is not to replace exact finite audits, but to expose the summable
proof split: short remaining spans can be charged to outer spectrum and
placement, while middle/far spans use the all-episode inner suppression.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from analyze_fullsplit_all_episode_gf import (
    all_episode_cauchy_occupancy_log2,
    parse_float_list,
)
from block_outer_upgrade_probe import load_local_spectrum, outer_block_gf_bounds, z_grid
from bound_fullsplit_episode_gaps import load_spectrum, split_entries
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list
from sum_fullsplit_episode_envelope import envelope_value, load_knot_csv, precompute_gap_sums


def lambda_grid(lo: float, hi: float, step: float) -> list[float]:
    n = int(math.floor((hi - lo) / step + 1e-12))
    return [lo + i * step for i in range(n + 1)]


def bucket_ranges(start: int, stop: int, step: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for gap_min in range(start, stop + 1, step):
        out.append((gap_min, min(stop, gap_min + step - 1)))
    return out


def wrapper_default_lambdas() -> list[float]:
    return [0.001, 0.003, 0.01, 0.02, 0.05, 0.08, 0.12, 0.2, 0.3, 0.5, 0.8, 1.2]


def wrapper_default_alphas() -> list[float]:
    return [0.5, 0.75, 0.9, 0.95, 0.98, 0.99, 0.99683772234, 0.999]


def wrapper_default_rhos() -> list[float]:
    return [
        1e-6,
        3e-6,
        1e-5,
        3e-5,
        1e-4,
        3e-4,
        5.6234132519e-4,
        1e-3,
        3e-3,
        1e-2,
        3e-2,
        1e-1,
        3e-1,
        5e-1,
        1.0,
        3.0,
        10.0,
        30.0,
        100.0,
        300.0,
        1000.0,
        3000.0,
        10000.0,
    ]


def log2_expm1_pos(x: float) -> float:
    if x <= 0.0:
        return float("-inf")
    if x < 50.0:
        return math.log2(math.expm1(x))
    return (x + math.log1p(-math.exp(-x))) / math.log(2.0)


def log2_arith_geom_sum(log2_r: float, lo: int, hi: int) -> float:
    """Return log2(sum_{x=lo}^hi (x+1) r^x)."""

    if lo > hi:
        return float("-inf")
    n = hi - lo
    if abs(log2_r) < 1e-8:
        count = n + 1
        return math.log2(count * (lo + hi + 2) / 2.0)

    ln_r = log2_r * math.log(2.0)

    def finite_sums(q: float) -> tuple[float, float]:
        # A=sum q^i, B=sum i q^i, i=0..n, for 0<q<1.
        if q == 0.0:
            return 1.0, 0.0
        qn1 = 0.0 if (n + 1) * math.log(q) < -745.0 else q ** (n + 1)
        qn2 = qn1 * q
        one_minus = 1.0 - q
        A = (1.0 - qn1) / one_minus
        B = (q - (n + 1) * qn1 + n * qn2) / (one_minus * one_minus)
        return A, B

    if log2_r < 0.0:
        r = math.exp(ln_r)
        A, B = finite_sums(r)
        inner = (lo + 1) * A + B
        return lo * log2_r + math.log2(inner)

    q = math.exp(-ln_r)
    A, B = finite_sums(q)
    inner = (hi + 1) * A - B
    return hi * log2_r + math.log2(inner)


def log2_one_plus_pow2(x: float) -> float:
    if x == float("-inf"):
        return 0.0
    if x > 40.0:
        return x + math.log2(1.0 + 2.0 ** (-x))
    return math.log2(1.0 + 2.0 ** x)


def log2_ht_tail_envelope(*, H: int, T: int, log_q: float, cutoff_bits: float = 120.0) -> float:
    """Return log2(sum_{k>=1} binom(H,k) binom(T,k)/(k+1) q^k)."""

    total = float("-inf")
    for k in range(1, min(H, T) + 1):
        term = (
            log2_binom(H, k)
            + log2_binom(T, k)
            - math.log2(k + 1)
            + k * log_q
        )
        total = log2add(total, term)
        if term < total - cutoff_bits and k > 4:
            break
    return total


def e01_cauchy_log2(
    *,
    b: int,
    T: int,
    H: int,
    distance: int,
    term_log2: float,
    entries: list[tuple[int, int, float]],
    lambdas: list[float],
    rhos: list[float],
) -> tuple[float, float, float, float]:
    x_min = 0 if H == 0 else max(1, (H + b - 1) // b)
    x_max = min(H, T)
    if x_min > x_max:
        return float("-inf"), float("nan"), float("nan"), float("nan")

    denom = log2_binom(b * T, H)
    best_e0 = float("inf")
    best_e1 = float("inf")
    best_lam = float("nan")
    best_rho = float("nan")
    for lam in lambdas:
        M = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
        if not (0.0 < M < 1.0):
            continue
        pref = lam * distance / math.log(2.0)
        e0 = pref + T * math.log2(M)
        if e0 < best_e0:
            best_e0 = e0
        log_m_over = math.log2(M / (1.0 - M))
        for rho in rhos:
            log_block = log2_expm1_pos(b * math.log1p(rho))
            log_g = log_block + log_m_over
            e1 = (
                term_log2
                + pref
                - denom
                - H * math.log2(rho)
                + log2_arith_geom_sum(log_g, max(1, x_min), x_max)
            )
            if e1 < best_e1:
                best_e1 = e1
                best_lam = lam
                best_rho = rho
    return min(0.0, log2add(best_e0, best_e1)), best_lam, float("nan"), best_rho


def eall_ratio_log2(
    *,
    b: int,
    T: int,
    H: int,
    distance: int,
    term_log2: float,
    entries: list[tuple[int, int, float]],
    lambdas: list[float],
    rhos: list[float],
) -> tuple[float, float, float, float]:
    x_min = 0 if H == 0 else max(1, (H + b - 1) // b)
    x_max = min(H, T)
    if x_min > x_max:
        return float("-inf"), float("nan"), float("nan"), float("nan")

    denom = log2_binom(b * T, H)
    best_e0 = float("inf")
    best_tail = float("inf")
    best_lam = float("nan")
    best_rho = float("nan")
    for lam in lambdas:
        M = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
        if not (0.0 < M < 1.0):
            continue
        pref = lam * distance / math.log(2.0)
        e0 = pref + T * math.log2(M)
        if e0 < best_e0:
            best_e0 = e0

        log_m_over = math.log2(M / (1.0 - M))
        tail_log2 = log2_ht_tail_envelope(
            H=H,
            T=T,
            log_q=term_log2 + math.log2(1.0 - M),
        )
        tail_factor_log2 = log2_one_plus_pow2(tail_log2)
        for rho in rhos:
            log_block = log2_expm1_pos(b * math.log1p(rho))
            log_g = log_block + log_m_over
            e1 = (
                term_log2
                + pref
                - denom
                - H * math.log2(rho)
                + log2_arith_geom_sum(log_g, max(1, x_min), x_max)
            )
            e_ge1 = e1 + tail_factor_log2
            if e_ge1 < best_tail:
                best_tail = e_ge1
                best_lam = lam
                best_rho = rho
    return min(0.0, log2add(best_e0, best_tail)), best_lam, float("nan"), best_rho


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-spectrum-csv", type=Path, default=Path(__file__).with_name("rm512_256_spectrum.csv"))
    parser.add_argument("--outer-blocks", type=int, default=4096)
    parser.add_argument("--outer-block-bits", type=int, default=256)
    parser.add_argument("--local-length", type=int, default=512)
    parser.add_argument("--local-distance", type=int, default=32)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--inner-spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--inner-block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--turnoff-log2", type=float, default=-62.4078758)
    parser.add_argument("--h-values", default="32:500")
    parser.add_argument("--first-r-values", default="1:64")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=1)
    parser.add_argument("--gap-stop", type=int, default=12000)
    parser.add_argument("--gap-step", type=int, default=4000)
    parser.add_argument(
        "--gap-sum-mode",
        choices=("exact", "endpoint"),
        default="exact",
        help=(
            "How to sum placements inside a gap bucket. endpoint uses the rigorous "
            "upper bound bucket_size*C(n_max,H), which is much faster for high H."
        ),
    )
    parser.add_argument(
        "--inner-T-by-gap",
        choices=("min", "max", "feasible-min"),
        default="min",
        help=(
            "Use the first, last, or earliest H-feasible T in each gap bucket when "
            "applying bucket-level inner bounds."
        ),
    )
    parser.add_argument(
        "--inner-mode-by-gap",
        default="cap,cauchy,cauchy",
        help="Comma-separated modes, one per gap bucket: cap, cauchy, e01cauchy, eallratio, or csv.",
    )
    parser.add_argument(
        "--inner-knot-csvs",
        help=(
            "Semicolon-separated CSV paths, one per gap bucket. Use none for "
            "non-csv buckets. Required for any csv bucket."
        ),
    )
    parser.add_argument("--knot-h-column", default="H")
    parser.add_argument("--knot-value-column", default="total_log2")
    parser.add_argument("--require-knot-coverage", action="store_true")
    parser.add_argument(
        "--use-lambda-range",
        action="store_true",
        help="Use --lambda-min/--lambda-max/--lambda-step instead of the default certificate grid.",
    )
    parser.add_argument("--lambda-min", type=float, default=0.001)
    parser.add_argument("--lambda-max", type=float, default=0.2)
    parser.add_argument("--lambda-step", type=float, default=0.05)
    parser.add_argument(
        "--lambdas",
        help=(
            "Comma list or lo:hi:step lambda values. If omitted, use the wrapper "
            "certificate grid unless --use-lambda-range is set."
        ),
    )
    parser.add_argument("--alphas", help="Comma list or lo:hi:step alpha fractions for a=alpha*M.")
    parser.add_argument("--rhos", help="Comma list or lo:hi:step rho values for the occupancy Cauchy bound.")
    parser.add_argument(
        "--cauchy-prune-cap-below-log2",
        type=float,
        help=(
            "If a row's cap-only outer+placement term is already below this log2 value, "
            "skip Cauchy optimization and keep the cap bound for that row."
        ),
    )
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    parser.add_argument("--output-csv", type=Path)
    parser.add_argument("--output-h-summary-csv", type=Path)
    args = parser.parse_args()

    h_values = parse_int_list(args.h_values)
    r_values = parse_int_list(args.first_r_values)
    h_max = max(h_values)
    d = math.floor(args.distance_delta * args.N)

    spectrum = load_local_spectrum(str(args.local_spectrum_csv))
    zs = z_grid(args.z_min, args.z_max, args.z_count)
    outer_logs, outer_z = outer_block_gf_bounds(
        blocks=args.outer_blocks,
        block_bits=args.outer_block_bits,
        local_length=args.local_length,
        d0=args.local_distance,
        h_max=h_max,
        zs=zs,
        model="spectrum-csv",
        spectrum=spectrum,
    )

    gaps = bucket_ranges(args.gap_start, args.gap_stop, args.gap_step)
    modes = [part.strip().lower() for part in args.inner_mode_by_gap.split(",") if part.strip()]
    if len(modes) != len(gaps):
        raise ValueError("--inner-mode-by-gap must have one mode per gap bucket")
    for mode in modes:
        if mode not in {"cap", "cauchy", "e01cauchy", "eallratio", "csv"}:
            raise ValueError(f"unknown inner mode {mode!r}")
    knot_lists: list[list[tuple[int, float]] | None] = [None for _ in gaps]
    if any(mode == "csv" for mode in modes):
        if not args.inner_knot_csvs:
            raise ValueError("--inner-knot-csvs is required when any bucket uses csv mode")
        parts = [part.strip() for part in args.inner_knot_csvs.split(";")]
        if len(parts) != len(gaps):
            raise ValueError("--inner-knot-csvs must have one semicolon-separated entry per gap bucket")
        for idx, (mode, part) in enumerate(zip(modes, parts)):
            if mode != "csv":
                continue
            if not part or part.lower() == "none":
                raise ValueError(f"bucket {idx} uses csv mode but has no CSV path")
            knot_lists[idx] = load_knot_csv(
                Path(part),
                h_column=args.knot_h_column,
                value_column=args.knot_value_column,
            )

    max_h_after_first = max((h - 1 for h in h_values), default=0)
    gap_sums = precompute_gap_sums(
        b=args.inner_block_bits,
        late_blocks=args.late_blocks,
        gap_ranges=gaps,
        max_h_after_first=max_h_after_first,
        mode=args.gap_sum_mode,
    )

    entries = split_entries(load_spectrum(args.inner_spectrum), args.inner_block_bits)
    if args.lambdas:
        lambdas = parse_float_list(args.lambdas)
    elif args.use_lambda_range:
        lambdas = lambda_grid(args.lambda_min, args.lambda_max, args.lambda_step)
    else:
        lambdas = wrapper_default_lambdas()
    alphas = parse_float_list(args.alphas) if args.alphas else wrapper_default_alphas()
    rhos = parse_float_list(args.rhos) if args.rhos else wrapper_default_rhos()
    inner_cache: dict[tuple[str, int, int], tuple[float, float, float, float]] = {}

    def inner_bound(bucket_idx: int, H: int) -> tuple[float, float, float, float]:
        mode = modes[bucket_idx]
        if mode == "cap":
            return 0.0, float("nan"), float("nan"), float("nan")
        if mode == "csv":
            knots = knot_lists[bucket_idx]
            assert knots is not None
            return (
                min(0.0, envelope_value(knots, H, require_coverage=args.require_knot_coverage)),
                float("nan"),
                float("nan"),
                float("nan"),
            )
        gap_min, gap_max = gaps[bucket_idx]
        if args.inner_T_by_gap == "max":
            T = args.late_blocks + gap_max - 1
        elif args.inner_T_by_gap == "feasible-min":
            T = max(args.late_blocks + gap_min - 1, (H + args.inner_block_bits - 1) // args.inner_block_bits)
            if T > args.late_blocks + gap_max - 1:
                return float("-inf"), float("nan"), float("nan"), float("nan")
        else:
            T = args.late_blocks + gap_min - 1
        key = (mode, T, H)
        cached = inner_cache.get(key)
        if cached is not None:
            return cached
        if mode == "e01cauchy":
            cached = e01_cauchy_log2(
                b=args.inner_block_bits,
                T=T,
                H=H,
                distance=d,
                term_log2=args.turnoff_log2,
                entries=entries,
                lambdas=lambdas,
                rhos=rhos,
            )
            inner_cache[key] = cached
            return cached
        if mode == "eallratio":
            cached = eall_ratio_log2(
                b=args.inner_block_bits,
                T=T,
                H=H,
                distance=d,
                term_log2=args.turnoff_log2,
                entries=entries,
                lambdas=lambdas,
                rhos=rhos,
            )
            inner_cache[key] = cached
            return cached
        best = float("inf")
        best_lam = float("nan")
        best_alpha = float("nan")
        best_rho = float("nan")
        for lam in lambdas:
            for alpha in alphas:
                for rho in rhos:
                    val = all_episode_cauchy_occupancy_log2(
                        b=args.inner_block_bits,
                        T=T,
                        H=H,
                        distance=d,
                        term_log2=args.turnoff_log2,
                        entries=entries,
                        lam=lam,
                        alpha=alpha,
                        rho=rho,
                    )
                    if val < best:
                        best = val
                        best_lam = lam
                        best_alpha = alpha
                        best_rho = rho
        cached = (best, best_lam, best_alpha, best_rho)
        inner_cache[key] = cached
        return cached

    rows = []
    total = float("-inf")
    by_gap = [float("-inf") for _ in gaps]
    by_mode: dict[str, float] = {}
    by_h: dict[int, float] = {}
    peak_by_h: dict[int, dict] = {}
    peak = None

    for h in h_values:
        outer = outer_logs[h] if h < len(outer_logs) else float("-inf")
        if outer == float("-inf"):
            continue
        for r in r_values:
            if r < 1 or r > min(h, args.inner_block_bits):
                continue
            H = h - r
            for idx, (gap_min, gap_max) in enumerate(gaps):
                gap_sum = gap_sums.get((H, idx), float("-inf"))
                if gap_sum == float("-inf"):
                    continue
                placement = (
                    log2_binom(args.inner_block_bits, r)
                    + gap_sum
                    - log2_binom(args.N, h)
                )
                if (
                    modes[idx] == "cauchy"
                    and args.cauchy_prune_cap_below_log2 is not None
                    and outer + placement <= args.cauchy_prune_cap_below_log2
                ):
                    inner, best_lam, best_alpha, best_rho = 0.0, float("nan"), float("nan"), float("nan")
                else:
                    inner, best_lam, best_alpha, best_rho = inner_bound(idx, H)
                term = outer + placement + inner
                total = log2add(total, term)
                by_gap[idx] = log2add(by_gap[idx], term)
                by_mode[modes[idx]] = log2add(by_mode.get(modes[idx], float("-inf")), term)
                by_h[h] = log2add(by_h.get(h, float("-inf")), term)
                row = {
                    "outer_weight": h,
                    "first_r": r,
                    "remaining_ones": H,
                    "gap_min": gap_min,
                    "gap_max": gap_max,
                    "bucket_T": (
                        args.late_blocks + gap_max - 1
                        if args.inner_T_by_gap == "max"
                        else max(
                            args.late_blocks + gap_min - 1,
                            (H + args.inner_block_bits - 1) // args.inner_block_bits,
                        )
                        if args.inner_T_by_gap == "feasible-min"
                        else args.late_blocks + gap_min - 1
                    ),
                    "inner_mode": modes[idx],
                    "outer_log2_bound": outer,
                    "outer_z": outer_z[h],
                    "placement_log2": placement,
                    "inner_log2": inner,
                    "best_lambda": best_lam,
                    "best_alpha": best_alpha,
                    "best_rho": best_rho,
                    "term_log2": term,
                }
                if args.output_csv:
                    rows.append(row)
                if peak is None or term > peak["term_log2"]:
                    peak = row
                h_peak = peak_by_h.get(h)
                if h_peak is None or term > h_peak["term_log2"]:
                    peak_by_h[h] = row

    if args.output_csv:
        fields = [
            "outer_weight",
            "first_r",
            "remaining_ones",
            "gap_min",
            "gap_max",
            "bucket_T",
            "inner_mode",
            "outer_log2_bound",
            "outer_z",
            "placement_log2",
            "inner_log2",
            "best_lambda",
            "best_alpha",
            "best_rho",
            "term_log2",
        ]
        with args.output_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    if args.output_h_summary_csv:
        fields = [
            "outer_weight",
            "aggregate_log2",
            "peak_first_r",
            "peak_remaining_ones",
            "peak_gap_min",
            "peak_gap_max",
            "peak_inner_mode",
            "peak_outer_log2_bound",
            "peak_placement_log2",
            "peak_inner_log2",
            "peak_best_lambda",
            "peak_best_alpha",
            "peak_best_rho",
            "peak_term_log2",
        ]
        with args.output_h_summary_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for h in sorted(by_h):
                h_peak = peak_by_h[h]
                writer.writerow(
                    {
                        "outer_weight": h,
                        "aggregate_log2": by_h[h],
                        "peak_first_r": h_peak["first_r"],
                        "peak_remaining_ones": h_peak["remaining_ones"],
                        "peak_gap_min": h_peak["gap_min"],
                        "peak_gap_max": h_peak["gap_max"],
                        "peak_inner_mode": h_peak["inner_mode"],
                        "peak_outer_log2_bound": h_peak["outer_log2_bound"],
                        "peak_placement_log2": h_peak["placement_log2"],
                        "peak_inner_log2": h_peak["inner_log2"],
                        "peak_best_lambda": h_peak["best_lambda"],
                        "peak_best_alpha": h_peak["best_alpha"],
                        "peak_best_rho": h_peak["best_rho"],
                        "peak_term_log2": h_peak["term_log2"],
                    }
                )

    print("Full-split piecewise certificate sum")
    print(f"N,{args.N}")
    print(f"delta,{args.distance_delta}")
    print(f"d,{d}")
    print(f"h_values,{args.h_values}")
    print(f"r_values,{args.first_r_values}")
    print(f"inner_modes,{','.join(modes)}")
    print(f"lambda_count,{len(lambdas)}")
    print(f"alpha_count,{len(alphas)}")
    print(f"rho_count,{len(rhos)}")
    print(f"inner_cache_entries,{len(inner_cache)}")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    for (gap_min, gap_max), mode, val in zip(gaps, modes, by_gap):
        print(f"gap_{gap_min}_{gap_max}_{mode}_log2,{val:.6f}")
    for mode, val in sorted(by_mode.items()):
        print(f"mode_{mode}_log2,{val:.6f}")
    if by_h:
        h_peak, h_val = max(by_h.items(), key=lambda kv: kv[1])
        print(f"peak_aggregate_h,{h_peak}")
        print(f"peak_aggregate_h_log2,{h_val:.6f}")
    if peak is not None:
        print(f"peak_h,{peak['outer_weight']}")
        print(f"peak_first_r,{peak['first_r']}")
        print(f"peak_gap_min,{peak['gap_min']}")
        print(f"peak_gap_max,{peak['gap_max']}")
        print(f"peak_inner_mode,{peak['inner_mode']}")
        print(f"peak_outer_log2_bound,{peak['outer_log2_bound']:.6f}")
        print(f"peak_placement_log2,{peak['placement_log2']:.6f}")
        print(f"peak_inner_log2,{peak['inner_log2']:.6f}")
        print(f"peak_term_log2,{peak['term_log2']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
