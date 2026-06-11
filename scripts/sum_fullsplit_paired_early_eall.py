#!/usr/bin/env python3
"""Paired-T early all-episode first-moment sum for the full-split inner.

This is the theorem-facing version of ``probe_fullsplit_early_paired_survival``.
For each first-active bucket it keeps the placement span T paired with the
inner bound, but it sums the all-episode wrapper by branches:

  e=0    : binom(bT,H) * exp(lambda d) * M(lambda)^T,
  e>=1  : the Cauchy e=1 envelope, with the same-lambda e>=2 multiplier.

The e>=1 branch cancels the occupancy denominator binom(bT,H), so its T-sum is
an explicit arithmetic-geometric endpoint-dominated sum.  This is much faster
than calling the per-T inner bound for every r and also exposes the proof split
we need for the early high-density certificate.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from block_outer_upgrade_probe import (
    load_local_spectrum,
    local_spectrum_is_complement_symmetric,
    outer_block_gf_bounds_for_weights,
    z_grid,
)
from bound_fullsplit_episode_gaps import load_spectrum, split_entries
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list
from sum_fullsplit_piecewise_certificate import (
    log2_arith_geom_sum,
    log2_expm1_pos,
    log2_ht_tail_envelope,
    log2_one_plus_pow2,
    parse_float_list,
)


def bucket_ranges(start: int, stop: int, step: int) -> list[tuple[int, int]]:
    return [
        (gap_min, min(stop, gap_min + step - 1))
        for gap_min in range(start, stop + 1, step)
    ]


def log2_sum_e0_T(
    *,
    b: int,
    H: int,
    T_min: int,
    T_max: int,
    pref: float,
    log_m: float,
) -> tuple[float, int, float]:
    total = float("-inf")
    peak_T = -1
    peak = float("-inf")
    for T in range(T_min, T_max + 1):
        if H > b * T:
            continue
        val = log2_binom(b * T, H) + pref + T * log_m
        total = log2add(total, val)
        if val > peak:
            peak = val
            peak_T = T
    return total, peak_T, peak


def log2_sum_ege1_T(
    *,
    b: int,
    H: int,
    T_min: int,
    T_max: int,
    term_log2: float,
    pref: float,
    log_rho: float,
    log_g: float,
    log_q: float,
) -> tuple[float, int, float]:
    x_min = 0 if H == 0 else max(1, (H + b - 1) // b)
    total = float("-inf")
    peak_T = -1
    peak = float("-inf")
    for T in range(T_min, T_max + 1):
        x_max = min(H, T)
        if x_min > x_max:
            continue
        s_t = log2_arith_geom_sum(log_g, max(1, x_min), x_max)
        tail = log2_one_plus_pow2(log2_ht_tail_envelope(H=H, T=T, log_q=log_q))
        val = term_log2 + pref - H * log_rho + s_t + tail
        total = log2add(total, val)
        if val > peak:
            peak = val
            peak_T = T
    return total, peak_T, peak


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-spectrum-csv", type=Path, default=Path(__file__).with_name("rm512_256_spectrum.csv"))
    parser.add_argument("--outer-blocks", type=int, default=4096)
    parser.add_argument("--outer-block-bits", type=int, default=256)
    parser.add_argument("--local-length", type=int, default=512)
    parser.add_argument("--local-distance", type=int, default=32)
    parser.add_argument("--inner-spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--turnoff-log2", type=float, default=-62.4078758)
    parser.add_argument("--h-values", required=True)
    parser.add_argument("--first-r-values", default="1:64")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=12001)
    parser.add_argument("--gap-stop", type=int, default=26819)
    parser.add_argument("--gap-step", type=int, default=5000)
    parser.add_argument("--lambdas", default="2.3")
    parser.add_argument("--rhos", default="1")
    parser.add_argument("--outer-complement-symmetry", action="store_true")
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    args = parser.parse_args()

    b = args.block_bits
    h_values = parse_int_list(args.h_values)
    r_values = parse_int_list(args.first_r_values)
    local_spectrum = load_local_spectrum(str(args.local_spectrum_csv))
    if args.outer_complement_symmetry:
        if args.outer_blocks * args.local_length != args.N:
            raise ValueError("--outer-complement-symmetry requires outer_blocks*local_length == N")
        if not local_spectrum_is_complement_symmetric(local_spectrum, args.local_length):
            raise ValueError("--outer-complement-symmetry requested, but local spectrum is not symmetric")
        outer_weights = sorted({min(h, args.N - h) for h in h_values})
    else:
        outer_weights = sorted(set(h_values))

    outer_logs, _outer_z = outer_block_gf_bounds_for_weights(
        blocks=args.outer_blocks,
        block_bits=args.outer_block_bits,
        local_length=args.local_length,
        d0=args.local_distance,
        weights=outer_weights,
        zs=z_grid(args.z_min, args.z_max, args.z_count),
        model="spectrum-csv",
        spectrum=local_spectrum,
    )

    d = math.floor(args.distance_delta * args.N)
    entries = split_entries(load_spectrum(args.inner_spectrum), b)
    lambdas = parse_float_list(args.lambdas)
    rhos = parse_float_list(args.rhos)
    buckets = bucket_ranges(args.gap_start, args.gap_stop, args.gap_step)

    mgf_by_lam: dict[float, float] = {}
    for lam in lambdas:
        M = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
        if 0.0 < M < 1.0:
            mgf_by_lam[lam] = M
    if not mgf_by_lam:
        raise SystemExit("no valid lambda values")

    total = float("-inf")
    by_h: dict[int, float] = {}
    by_bucket = [float("-inf") for _ in buckets]
    peak = None
    branch_totals = {"e0": float("-inf"), "ege1": float("-inf")}

    for h in h_values:
        outer_h = min(h, args.N - h) if args.outer_complement_symmetry else h
        outer = outer_logs.get(outer_h, float("-inf"))
        if outer == float("-inf"):
            continue
        common_h = outer - log2_binom(args.N, h)
        h_total = float("-inf")
        for r in r_values:
            if r < 1 or r > min(h, b):
                continue
            H = h - r
            common = common_h + log2_binom(b, r)
            for idx, (gap_min, gap_max) in enumerate(buckets):
                T_min = args.late_blocks + gap_min - 1
                T_max = args.late_blocks + gap_max - 1
                if H > b * T_max:
                    continue

                best_e0 = float("inf")
                best_e0_meta = None
                best_ege1 = float("inf")
                best_ege1_meta = None
                for lam, M in mgf_by_lam.items():
                    pref = lam * d / math.log(2.0)
                    log_m = math.log2(M)
                    e0_sum, e0_T, e0_peak = log2_sum_e0_T(
                        b=b,
                        H=H,
                        T_min=T_min,
                        T_max=T_max,
                        pref=pref,
                        log_m=log_m,
                    )
                    if e0_sum < best_e0:
                        best_e0 = e0_sum
                        best_e0_meta = (lam, float("nan"), e0_T, e0_peak)

                    log_q = args.turnoff_log2 + math.log2(1.0 - M)
                    log_m_over = math.log2(M / (1.0 - M))
                    for rho in rhos:
                        log_block = log2_expm1_pos(b * math.log1p(rho))
                        log_g = log_block + log_m_over
                        ege1_sum, ege1_T, ege1_peak = log2_sum_ege1_T(
                            b=b,
                            H=H,
                            T_min=T_min,
                            T_max=T_max,
                            term_log2=args.turnoff_log2,
                            pref=pref,
                            log_rho=math.log2(rho),
                            log_g=log_g,
                            log_q=log_q,
                        )
                        if ege1_sum < best_ege1:
                            best_ege1 = ege1_sum
                            best_ege1_meta = (lam, rho, ege1_T, ege1_peak)

                e0_term = common + best_e0
                ege1_term = common + best_ege1
                term = log2add(e0_term, ege1_term)
                total = log2add(total, term)
                h_total = log2add(h_total, term)
                by_bucket[idx] = log2add(by_bucket[idx], term)
                branch_totals["e0"] = log2add(branch_totals["e0"], e0_term)
                branch_totals["ege1"] = log2add(branch_totals["ege1"], ege1_term)
                if peak is None or term > peak["term_log2"]:
                    peak_branch = "e0" if e0_term >= ege1_term else "ege1"
                    peak_meta = best_e0_meta if peak_branch == "e0" else best_ege1_meta
                    peak = {
                        "term_log2": term,
                        "h": h,
                        "r": r,
                        "H": H,
                        "gap_min": gap_min,
                        "gap_max": gap_max,
                        "outer": outer,
                        "common": common,
                        "e0_term": e0_term,
                        "ege1_term": ege1_term,
                        "branch": peak_branch,
                        "lambda": peak_meta[0] if peak_meta else float("nan"),
                        "rho": peak_meta[1] if peak_meta else float("nan"),
                        "peak_T": peak_meta[2] if peak_meta else -1,
                        "branch_peak_log2": peak_meta[3] if peak_meta else float("nan"),
                    }
        by_h[h] = h_total

    print("Full-split paired early all-episode split sum")
    print(f"N,{args.N}")
    print(f"delta,{args.distance_delta}")
    print(f"d,{d}")
    print(f"h_values,{args.h_values}")
    print(f"r_values,{args.first_r_values}")
    print(f"outer_complement_symmetry,{int(args.outer_complement_symmetry)}")
    print(f"lambda_count,{len(mgf_by_lam)}")
    print(f"rho_count,{len(rhos)}")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    print(f"branch_e0_log2,{branch_totals['e0']:.6f}")
    print(f"branch_ege1_log2,{branch_totals['ege1']:.6f}")
    for idx, (gap_min, gap_max) in enumerate(buckets):
        print(f"gap_{gap_min}_{gap_max}_log2,{by_bucket[idx]:.6f}")
    for h in h_values:
        print(f"h_{h}_log2,{by_h.get(h, float('-inf')):.6f}")
    if peak:
        print(f"peak_h,{peak['h']}")
        print(f"peak_first_r,{peak['r']}")
        print(f"peak_remaining_ones,{peak['H']}")
        print(f"peak_gap_min,{peak['gap_min']}")
        print(f"peak_gap_max,{peak['gap_max']}")
        print(f"peak_branch,{peak['branch']}")
        print(f"peak_lambda,{peak['lambda']:.6g}")
        print(f"peak_rho,{peak['rho']:.6g}")
        print(f"peak_T,{peak['peak_T']}")
        print(f"peak_outer_log2,{peak['outer']:.6f}")
        print(f"peak_common_log2,{peak['common']:.6f}")
        print(f"peak_e0_term_log2,{peak['e0_term']:.6f}")
        print(f"peak_ege1_term_log2,{peak['ege1_term']:.6f}")
        print(f"peak_term_log2,{peak['term_log2']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
