#!/usr/bin/env python3
"""Interval certificate for the full-split high-density early branch.

This helper packages the proof reductions used after the paired-T diagnostics:

* use complement symmetry of the RM512 direct-sum outer;
* freeze one Chernoff pole z for an h-interval, so the outer/placement term
  is bounded by interval endpoints;
* use right-endpoint dominance in T for the paired e>=1 branch;
* sum all first-block weights by the crude exact factor sum_r binom(b,r).

For rho != 1 the factor rho^{-H}=rho^{-h+r} is split into a convex-in-h
interval term rho^{-h} and the first-block sum sum_r binom(b,r)rho^r.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

from block_outer_upgrade_probe import (
    load_local_spectrum,
    local_spectrum_is_complement_symmetric,
    local_spectrum_log2,
    log2_sub_one,
    z_grid,
)
from bound_fullsplit_episode_gaps import load_spectrum, split_entries
from dense_largek_eval import log2_binom, log2add
from sum_fullsplit_piecewise_certificate import (
    log2_arith_geom_sum,
    log2_expm1_pos,
    log2_ht_tail_envelope,
    log2_one_plus_pow2,
)


@dataclass(frozen=True)
class Interval:
    h_min: int
    h_max: int
    rho: float | None = None


def parse_intervals(text: str) -> list[Interval]:
    out: list[Interval] = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        fields = part.split(":")
        if len(fields) not in (2, 3):
            raise ValueError(f"bad interval {part!r}; expected h0:h1[:rho]")
        lo, hi = [int(x) for x in fields[:2]]
        if lo > hi:
            raise ValueError(f"bad interval {part!r}")
        rho = None if len(fields) == 2 else float(fields[2])
        out.append(Interval(lo, hi, rho))
    return out


def fixed_z_outer_logs(
    *,
    blocks: int,
    spectrum: list[tuple[int, int]],
    zs: list[float],
) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    for z in zs:
        local = local_spectrum_log2(spectrum, z)
        out.append((z, log2_sub_one(blocks * local)))
    return out


def max_outer_minus_volume_on_interval(
    *,
    N: int,
    interval: Interval,
    z: float,
    global_log: float,
    log_rho: float,
) -> float:
    log_z = math.log2(z)

    def value(h: int) -> float:
        w = N - h
        return global_log - w * log_z - log2_binom(N, h) - h * log_rho

    # For fixed z,rho this is linear in w=N-h minus log binom(N,w), hence
    # convex in w.
    return max(value(interval.h_min), value(interval.h_max))


def max_outer_on_interval(*, N: int, interval: Interval, z: float, global_log: float) -> float:
    log_z = math.log2(z)
    # z<1, so the complement-side outer bound is largest at largest complement
    # weight, i.e. the left endpoint of the h-interval.
    return global_log - (N - interval.h_min) * log_z


def best_fixed_z_outer_minus_volume(
    *,
    N: int,
    interval: Interval,
    local_logs: list[tuple[float, float]],
    log_rho: float,
) -> tuple[float, float]:
    best = float("inf")
    best_z = float("nan")
    for z, global_log in local_logs:
        if global_log == float("-inf"):
            continue
        val = max_outer_minus_volume_on_interval(
            N=N,
            interval=interval,
            z=z,
            global_log=global_log,
            log_rho=log_rho,
        )
        if val < best:
            best = val
            best_z = z
    return best, best_z


def best_fixed_z_outer(
    *,
    N: int,
    interval: Interval,
    local_logs: list[tuple[float, float]],
) -> tuple[float, float]:
    best = float("inf")
    best_z = float("nan")
    for z, global_log in local_logs:
        if global_log == float("-inf"):
            continue
        val = max_outer_on_interval(N=N, interval=interval, z=z, global_log=global_log)
        if val < best:
            best = val
            best_z = z
    return best, best_z


def log2_one_minus_pow2(log_x: float) -> float:
    if log_x >= 0.0:
        return float("nan")
    if log_x < -40.0:
        return -2.0**log_x / math.log(2.0)
    return math.log2(1.0 - 2.0**log_x)


def endpoint_geom_overhead(log_ratio: float) -> float:
    if log_ratio >= 0.0:
        return float("inf")
    return -log2_one_minus_pow2(log_ratio)


def log2_weighted_first_block_sum(*, b: int, rho: float) -> float:
    total = float("-inf")
    log_rho = math.log2(rho)
    for r in range(1, b + 1):
        total = log2add(total, log2_binom(b, r) + r * log_rho)
    return total


@dataclass(frozen=True)
class RhoParams:
    rho: float
    log_rho: float
    log_g: float
    ege1_t_overhead: float
    r_sum_log2: float
    s_full: float


def rho_params(
    *,
    b: int,
    T_max: int,
    M: float,
    rho: float,
) -> RhoParams:
    log_g = log2_expm1_pos(b * math.log1p(rho)) + math.log2(M / (1.0 - M))
    if log_g <= 1.0:
        raise ValueError(f"paired endpoint dominance needs G>2; rho={rho:g} gives log2 G={log_g:.6f}")
    return RhoParams(
        rho=rho,
        log_rho=math.log2(rho),
        log_g=log_g,
        ege1_t_overhead=endpoint_geom_overhead(-math.log2(2.0**log_g - 1.0)),
        r_sum_log2=log2_weighted_first_block_sum(b=b, rho=rho),
        s_full=log2_arith_geom_sum(log_g, 1, T_max),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-spectrum-csv", type=Path, default=Path(__file__).with_name("rm512_256_spectrum.csv"))
    parser.add_argument("--inner-spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--outer-blocks", type=int, default=4096)
    parser.add_argument("--local-length", type=int, default=512)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--turnoff-log2", type=float, default=-62.4078758)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=22001)
    parser.add_argument("--gap-stop", type=int, default=26819)
    parser.add_argument(
        "--intervals",
        default="1115101:1115201",
        help="Semicolon-separated h0:h1 or h0:h1:rho rows. Per-row rho overrides --rho.",
    )
    parser.add_argument("--lambda-value", type=float, default=2.3)
    parser.add_argument("--rho", type=float, default=1.0)
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    args = parser.parse_args()

    if args.outer_blocks * args.local_length != args.N:
        raise ValueError("expected outer_blocks*local_length == N")

    local_spectrum = load_local_spectrum(str(args.local_spectrum_csv))
    if not local_spectrum_is_complement_symmetric(local_spectrum, args.local_length):
        raise ValueError("local outer spectrum is not complement-symmetric")

    b = args.block_bits
    d = math.floor(args.distance_delta * args.N)
    T_min = args.late_blocks + args.gap_start - 1
    T_max = args.late_blocks + args.gap_stop - 1
    n_endpoint = b * T_max
    intervals = parse_intervals(args.intervals)

    entries = split_entries(load_spectrum(args.inner_spectrum), b)
    lam = args.lambda_value
    M = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
    if not (0.0 < M < 1.0):
        raise ValueError("invalid lambda: M_+(lambda) is not in (0,1)")

    pref = lam * d / math.log(2.0)
    log_m = math.log2(M)
    log_q = args.turnoff_log2 + math.log2(1.0 - M)

    local_logs = fixed_z_outer_logs(
        blocks=args.outer_blocks,
        spectrum=local_spectrum,
        zs=z_grid(args.z_min, args.z_max, args.z_count),
    )

    total = float("-inf")
    print("Full-split high-density interval certificate")
    print(f"N,{args.N}")
    print(f"delta,{args.distance_delta}")
    print(f"d,{d}")
    print(f"T_min,{T_min}")
    print(f"T_max,{T_max}")
    print(f"lambda,{lam:.12g}")
    print(f"default_rho,{args.rho:.12g}")
    print(f"M,{M:.12g}")
    print(
        "h_min,h_max,rho,count,e0_log2,ege1_log2,total_log2,"
        "outer_z_e0,outer_z_ege1,log2_g,e0_T_overhead,ege1_T_overhead,ege1_tail_log2"
    )

    rho_cache: dict[float, RhoParams] = {}
    for interval in intervals:
        rho_value = args.rho if interval.rho is None else interval.rho
        params = rho_cache.get(rho_value)
        if params is None:
            params = rho_params(b=b, T_max=T_max, M=M, rho=rho_value)
            rho_cache[rho_value] = params

        if interval.h_min <= args.N // 2:
            raise ValueError("interval must be on high-density complement side")
        if interval.h_max > args.N:
            raise ValueError("interval must not exceed N")

        h_count_log2 = math.log2(interval.h_max - interval.h_min + 1)
        H_min = interval.h_min - b
        H_max = interval.h_max - 1

        outer_e0, z_e0 = best_fixed_z_outer(
            N=args.N,
            interval=interval,
            local_logs=local_logs,
        )
        outer_vol, z_ege1 = best_fixed_z_outer_minus_volume(
            N=args.N,
            interval=interval,
            local_logs=local_logs,
            log_rho=params.log_rho,
        )

        e0_ratio_log2 = -log_m + b * math.log2(1.0 - H_min / n_endpoint)
        e0_t_overhead = endpoint_geom_overhead(e0_ratio_log2)
        if e0_t_overhead == float("inf"):
            raise ValueError(
                f"e=0 endpoint dominance failed on {interval.h_min}:{interval.h_max}: "
                f"log2 ratio {e0_ratio_log2:.6f}"
            )

        e0 = h_count_log2 + outer_e0 + pref + T_max * log_m + e0_t_overhead

        tail = log2_one_plus_pow2(log2_ht_tail_envelope(H=H_max, T=T_max, log_q=log_q))
        ege1_inner = args.turnoff_log2 + pref + params.s_full + tail
        ege1 = h_count_log2 + outer_vol + params.r_sum_log2 + ege1_inner + params.ege1_t_overhead

        interval_total = log2add(e0, ege1)
        total = log2add(total, interval_total)
        print(
            f"{interval.h_min},{interval.h_max},{rho_value:.12g},{interval.h_max - interval.h_min + 1},"
            f"{e0:.6f},{ege1:.6f},{interval_total:.6f},"
            f"{z_e0:.12g},{z_ege1:.12g},{params.log_g:.6f},"
            f"{e0_t_overhead:.6f},{params.ege1_t_overhead:.6f},{tail:.12g}"
        )

    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
