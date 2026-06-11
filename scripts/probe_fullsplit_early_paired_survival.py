#!/usr/bin/env python3
"""Probe paired-T early high-weight full-split bounds.

The standard piecewise helper can combine endpoint placement at a bucket's
right endpoint with an inner bound at the left endpoint.  That is rigorous when
the inner bound is monotone and the placement loss is small, but it is far too
pessimistic for high-density early buckets: placement is dominated by large T,
where the survival exponent is also much stronger.

This diagnostic keeps T paired:

    sum_T binom(b,r) binom(bT,H) / binom(N,h) * survival(T).

It is deliberately a probe, not the final certificate: an optional
``episode-slack-bits`` parameter can be added to the survival exponent to test
how much room remains for a later all-episode multiplier.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from block_outer_upgrade_probe import (
    load_local_spectrum,
    local_spectrum_is_complement_symmetric,
    outer_block_gf_bounds,
    z_grid,
)
from bound_fullsplit_episode_gaps import load_spectrum, split_entries
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list
from sum_fullsplit_piecewise_certificate import parse_float_list


def bucket_ranges(start: int, stop: int, step: int) -> list[tuple[int, int]]:
    return [
        (gap_min, min(stop, gap_min + step - 1))
        for gap_min in range(start, stop + 1, step)
    ]


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
    parser.add_argument("--h-values", required=True)
    parser.add_argument("--first-r-values", default="1:64")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=12001)
    parser.add_argument("--gap-stop", type=int, default=26819)
    parser.add_argument("--gap-step", type=int, default=5000)
    parser.add_argument("--lambdas", default="0.05:20:0.05")
    parser.add_argument("--episode-slack-bits", type=float, default=0.0)
    parser.add_argument(
        "--outer-complement-symmetry",
        action="store_true",
        help="Use A_h=A_{N-h} for complement-symmetric local spectra, avoiding high-weight outer tables.",
    )
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    args = parser.parse_args()

    h_values = parse_int_list(args.h_values)
    r_values = parse_int_list(args.first_r_values)
    if args.outer_complement_symmetry:
        if args.outer_blocks * args.local_length != args.N:
            raise ValueError("--outer-complement-symmetry requires outer_blocks*local_length == N")
        local_spectrum = load_local_spectrum(str(args.local_spectrum_csv))
        if not local_spectrum_is_complement_symmetric(local_spectrum, args.local_length):
            raise ValueError("--outer-complement-symmetry requested, but local spectrum is not symmetric")
        h_max = max(min(h, args.N - h) for h in h_values)
    else:
        local_spectrum = load_local_spectrum(str(args.local_spectrum_csv))
        h_max = max(h_values)
    d = math.floor(args.distance_delta * args.N)

    outer_logs, _outer_z = outer_block_gf_bounds(
        blocks=args.outer_blocks,
        block_bits=args.outer_block_bits,
        local_length=args.local_length,
        d0=args.local_distance,
        h_max=h_max,
        zs=z_grid(args.z_min, args.z_max, args.z_count),
        model="spectrum-csv",
        spectrum=local_spectrum,
    )

    entries = split_entries(load_spectrum(args.inner_spectrum), args.block_bits)
    lambdas = parse_float_list(args.lambdas)
    log_mgfs = []
    for lam in lambdas:
        mgf = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
        if 0.0 < mgf < 1.0:
            log_mgfs.append((lam, math.log2(mgf)))
    if not log_mgfs:
        raise SystemExit("no valid lambda values")

    buckets = bucket_ranges(args.gap_start, args.gap_stop, args.gap_step)
    t_values: list[tuple[int, int, int]] = []
    for idx, (gap_min, gap_max) in enumerate(buckets):
        for gap in range(gap_min, gap_max + 1):
            t_values.append((idx, gap, args.late_blocks + gap - 1))
    survival: dict[int, tuple[float, float]] = {}
    inv_log2 = 1.0 / math.log(2.0)
    for _idx, _gap, T in t_values:
        best = float("inf")
        best_lam = float("nan")
        for lam, log_mgf in log_mgfs:
            val = lam * d * inv_log2 + T * log_mgf
            if val < best:
                best = val
                best_lam = lam
        survival[T] = (min(0.0, best) + args.episode_slack_bits, best_lam)

    total = float("-inf")
    by_bucket = [float("-inf") for _ in buckets]
    peak = None
    rows_by_h: dict[int, float] = {}
    for h in h_values:
        outer_h = min(h, args.N - h) if args.outer_complement_symmetry else h
        outer = outer_logs[outer_h] if 0 <= outer_h < len(outer_logs) else float("-inf")
        if outer == float("-inf"):
            continue
        denom = log2_binom(args.N, h)
        h_total = float("-inf")
        for r in r_values:
            if r < 1 or r > min(h, args.block_bits):
                continue
            H = h - r
            choose_r = log2_binom(args.block_bits, r)
            for idx, gap, T in t_values:
                if H > args.block_bits * T:
                    continue
                inner, lam = survival[T]
                placement = choose_r + log2_binom(args.block_bits * T, H) - denom
                term = outer + placement + inner
                total = log2add(total, term)
                h_total = log2add(h_total, term)
                by_bucket[idx] = log2add(by_bucket[idx], term)
                if peak is None or term > peak["term_log2"]:
                    peak = {
                        "term_log2": term,
                        "h": h,
                        "r": r,
                        "H": H,
                        "gap": gap,
                        "T": T,
                        "outer": outer,
                        "placement": placement,
                        "inner": inner,
                        "lambda": lam,
                    }
        rows_by_h[h] = h_total

    print("Full-split early paired survival probe")
    print(f"N,{args.N}")
    print(f"delta,{args.distance_delta}")
    print(f"d,{d}")
    print(f"h_values,{args.h_values}")
    print(f"outer_complement_symmetry,{int(args.outer_complement_symmetry)}")
    print(f"episode_slack_bits,{args.episode_slack_bits:.6f}")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    for idx, (gap_min, gap_max) in enumerate(buckets):
        print(f"gap_{gap_min}_{gap_max}_log2,{by_bucket[idx]:.6f}")
    for h in h_values:
        print(f"h_{h}_log2,{rows_by_h.get(h, float('-inf')):.6f}")
    if peak:
        print(f"peak_h,{peak['h']}")
        print(f"peak_first_r,{peak['r']}")
        print(f"peak_remaining_ones,{peak['H']}")
        print(f"peak_gap,{peak['gap']}")
        print(f"peak_T,{peak['T']}")
        print(f"peak_lambda,{peak['lambda']:.6g}")
        print(f"peak_outer_log2,{peak['outer']:.6f}")
        print(f"peak_placement_log2,{peak['placement']:.6f}")
        print(f"peak_inner_log2,{peak['inner']:.6f}")
        print(f"peak_term_log2,{peak['term_log2']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
