#!/usr/bin/env python3
"""Paired-T interval certificate for the early full-split branch.

This helper packages the proof-shaped bound for the remaining early
middle/high-density range.  Unlike endpoint-placement rows, it keeps the
placement length T paired with the survival bound before summing over T.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

from block_outer_upgrade_probe import load_local_spectrum, local_spectrum_log2, log2_sub_one, z_grid
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
    start: int
    stop: int
    lam: float
    rho: float

    @property
    def label(self) -> str:
        return f"{self.start}--{self.stop}"


def parse_intervals(text: str) -> list[Interval]:
    out: list[Interval] = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        fields = part.split(":")
        if len(fields) != 4:
            raise ValueError(f"bad interval {part!r}; expected h0:h1:lambda:rho")
        start, stop = int(fields[0]), int(fields[1])
        if start > stop:
            raise ValueError(f"bad interval {part!r}: start > stop")
        out.append(Interval(start, stop, float(fields[2]), float(fields[3])))
    return out


def log2add_many(values: list[float]) -> float:
    total = float("-inf")
    for value in values:
        total = log2add(total, value)
    return total


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


def fixed_z_outer_logs(
    *,
    blocks: int,
    spectrum: list[tuple[int, int]],
    zs: list[float],
) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    for z in zs:
        local = local_spectrum_log2(spectrum, z)
        global_log = log2_sub_one(blocks * local)
        if global_log != float("-inf"):
            out.append((z, global_log))
    return out


def best_e0_volume_for_t(
    *,
    n_total: int,
    interval: Interval,
    n_live_plus_first: int,
    local_logs: list[tuple[float, float]],
) -> tuple[float, float]:
    best = float("inf")
    best_z = float("nan")
    for z, global_log in local_logs:
        log_z = math.log2(z)
        candidates = [interval.start, interval.stop]
        critical = (n_live_plus_first - z * n_total) / (1.0 - z)
        for h in [math.floor(critical), math.ceil(critical)]:
            if interval.start <= h <= interval.stop:
                candidates.append(h)

        worst = float("-inf")
        for h in candidates:
            if h < 0 or h > n_live_plus_first:
                continue
            value = (
                global_log
                - h * log_z
                + log2_binom(n_live_plus_first, h)
                - log2_binom(n_total, h)
            )
            if value > worst:
                worst = value
        if worst < best:
            best = worst
            best_z = z
    return best, best_z


def e0_paired_interval_bound(
    *,
    b: int,
    n_total: int,
    interval: Interval,
    t_min: int,
    t_max: int,
    pref: float,
    log_m: float,
    local_logs: list[tuple[float, float]],
) -> tuple[float, float]:
    total = float("-inf")
    best_z_at_peak = float("nan")
    peak = float("-inf")
    for t in range(t_min, t_max + 1):
        value, z = best_e0_volume_for_t(
            n_total=n_total,
            interval=interval,
            n_live_plus_first=b * t + b,
            local_logs=local_logs,
        )
        term = value + pref + t * log_m
        total = log2add(total, term)
        if term > peak:
            peak = term
            best_z_at_peak = z
    return math.log2(interval.stop - interval.start + 1) + total, best_z_at_peak


def best_low_outer_minus_volume(
    *,
    n_total: int,
    interval: Interval,
    local_logs: list[tuple[float, float]],
    log_rho: float,
) -> tuple[float, float]:
    best = float("inf")
    best_z = float("nan")
    for z, global_log in local_logs:
        log_z = math.log2(z)

        def value(h: int) -> float:
            return global_log - h * log_z - log2_binom(n_total, h) - h * log_rho

        # For fixed z and rho this is linear in h minus the concave log-volume,
        # hence convex; the interval maximum is at an endpoint.
        val = max(value(interval.start), value(interval.stop))
        if val < best:
            best = val
            best_z = z
    return best, best_z


def log2_weighted_first_block_sum(*, b: int, rho: float) -> float:
    total = float("-inf")
    log_rho = math.log2(rho)
    for r in range(1, b + 1):
        total = log2add(total, log2_binom(b, r) + r * log_rho)
    return total


def interval_total(
    *,
    interval: Interval,
    b: int,
    n_total: int,
    distance: int,
    t_min: int,
    t_max: int,
    local_logs: list[tuple[float, float]],
    entries: list[tuple[int, int, float]],
    turnoff_log2: float,
) -> tuple[float, dict[str, float]]:
    lam = interval.lam
    rho = interval.rho
    m_plus = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
    if not (0.0 < m_plus < 1.0):
        raise ValueError(f"{interval.label}: invalid lambda {lam:g}")
    pref = lam * distance / math.log(2.0)
    log_m = math.log2(m_plus)
    log_q = turnoff_log2 + math.log2(1.0 - m_plus)
    log_rho = math.log2(rho)
    log_g = log2_expm1_pos(b * math.log1p(rho)) + math.log2(m_plus / (1.0 - m_plus))
    if log_g <= 1.0:
        raise ValueError(
            f"{interval.label}: paired endpoint dominance needs G>2; "
            f"lambda={lam:g}, rho={rho:g} gives log2 G={log_g:.6f}"
        )

    h_count = math.log2(interval.stop - interval.start + 1)
    e0, z_e0 = e0_paired_interval_bound(
        b=b,
        n_total=n_total,
        interval=interval,
        t_min=t_min,
        t_max=t_max,
        pref=pref,
        log_m=log_m,
        local_logs=local_logs,
    )

    outer_vol, z_ege1 = best_low_outer_minus_volume(
        n_total=n_total,
        interval=interval,
        local_logs=local_logs,
        log_rho=log_rho,
    )
    t_overhead = endpoint_geom_overhead(-math.log2(2.0**log_g - 1.0))
    r_sum = log2_weighted_first_block_sum(b=b, rho=rho)
    s_full = log2_arith_geom_sum(log_g, 1, t_max)
    tail = log2_one_plus_pow2(
        log2_ht_tail_envelope(H=interval.stop - 1, T=t_max, log_q=log_q)
    )
    ege1 = h_count + outer_vol + r_sum + turnoff_log2 + pref + s_full + tail + t_overhead
    return log2add(e0, ege1), {
        "e0": e0,
        "ege1": ege1,
        "outer_z_e0": z_e0,
        "outer_z_ege1": z_ege1,
        "log2_g": log_g,
        "t_overhead": t_overhead,
        "tail": tail,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-spectrum-csv", type=Path, default=Path(__file__).with_name("rm512_256_spectrum.csv"))
    parser.add_argument("--inner-spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--outer-blocks", type=int, default=4096)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--turnoff-log2", type=float, default=-62.4078758)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=12001)
    parser.add_argument("--gap-stop", type=int, default=26819)
    parser.add_argument(
        "--intervals",
        required=True,
        help="Semicolon-separated h0:h1:lambda:rho rows.",
    )
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    args = parser.parse_args()

    intervals = parse_intervals(args.intervals)
    t_min = args.late_blocks + args.gap_start - 1
    t_max = args.late_blocks + args.gap_stop - 1
    d = math.floor(args.distance_delta * args.N)
    local_spectrum = load_local_spectrum(str(args.local_spectrum_csv))
    local_logs = fixed_z_outer_logs(
        blocks=args.outer_blocks,
        spectrum=local_spectrum,
        zs=z_grid(args.z_min, args.z_max, args.z_count),
    )
    entries = split_entries(load_spectrum(args.inner_spectrum), args.block_bits)

    total = float("-inf")
    print("Full-split early paired-T interval certificate")
    print(f"N,{args.N}")
    print(f"delta,{args.distance_delta}")
    print(f"d,{d}")
    print(f"T_min,{t_min}")
    print(f"T_max,{t_max}")
    print("h_min,h_max,lambda,rho,e0_log2,ege1_log2,total_log2,outer_z_e0,outer_z_ege1,log2_g,T_overhead,tail_log2")
    for interval in intervals:
        value, meta = interval_total(
            interval=interval,
            b=args.block_bits,
            n_total=args.N,
            distance=d,
            t_min=t_min,
            t_max=t_max,
            local_logs=local_logs,
            entries=entries,
            turnoff_log2=args.turnoff_log2,
        )
        total = log2add(total, value)
        print(
            f"{interval.start},{interval.stop},{interval.lam:.12g},{interval.rho:.12g},"
            f"{meta['e0']:.6f},{meta['ege1']:.6f},{value:.6f},"
            f"{meta['outer_z_e0']:.12g},{meta['outer_z_ege1']:.12g},"
            f"{meta['log2_g']:.6f},{meta['t_overhead']:.6f},{meta['tail']:.12g}"
        )
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
