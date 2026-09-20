#!/usr/bin/env python3
"""Vectorized interval certificate for early endpoint full-split rows.

This is a faster companion to ``sum_fullsplit_piecewise_certificate.py`` for the
early first-active region ``T > 17948``.  It evaluates the same conservative
endpoint-placement, same-pole ``eallratio`` bound, but over whole h-intervals
using vectorized binomial recurrences.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from block_outer_upgrade_probe import load_local_spectrum, local_spectrum_log2, log2_sub_one, z_grid
from bound_fullsplit_episode_gaps import load_spectrum, split_entries
from dense_largek_eval import log2_binom, log2add
from sum_fullsplit_piecewise_certificate import (
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


def log2_sum_exp(values: np.ndarray) -> float:
    if values.size == 0:
        return float("-inf")
    peak = float(np.max(values))
    if not math.isfinite(peak):
        return float("-inf")
    return peak + math.log2(float(np.sum(np.exp2(values - peak))))


def log2add_many(values: list[float]) -> float:
    total = float("-inf")
    for value in values:
        total = log2add(total, value)
    return total


def log2_binom_range(n: int, k_min: int, k_max: int) -> np.ndarray:
    if k_min > k_max:
        return np.empty(0, dtype=np.float64)
    vals = np.empty(k_max - k_min + 1, dtype=np.float64)
    vals[0] = log2_binom(n, k_min)
    for idx, k in enumerate(range(k_min, k_max), start=1):
        vals[idx] = vals[idx - 1] + math.log2((n - k) / (k + 1))
    return vals


def log2_arith_geom_sum_vec(log2_r: float, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Vectorized log2(sum_{x=lo}^hi (x+1) r^x)."""

    n = hi - lo
    out = np.full(lo.shape, float("-inf"), dtype=np.float64)
    valid = n >= 0
    if not np.any(valid):
        return out
    lo_v = lo[valid].astype(np.float64)
    hi_v = hi[valid].astype(np.float64)
    n_v = n[valid].astype(np.float64)

    if abs(log2_r) < 1e-8:
        count = n_v + 1.0
        out[valid] = np.log2(count * (lo_v + hi_v + 2.0) / 2.0)
        return out

    if log2_r < 0.0:
        r = 2.0**log2_r
        qn1 = np.power(r, n_v + 1.0)
        qn2 = qn1 * r
        one_minus = 1.0 - r
        a = (1.0 - qn1) / one_minus
        b = (r - (n_v + 1.0) * qn1 + n_v * qn2) / (one_minus * one_minus)
        inner = (lo_v + 1.0) * a + b
        out[valid] = lo_v * log2_r + np.log2(inner)
        return out

    q = 2.0 ** (-log2_r)
    qn1 = np.power(q, n_v + 1.0)
    qn2 = qn1 * q
    one_minus = 1.0 - q
    a = (1.0 - qn1) / one_minus
    b = (q - (n_v + 1.0) * qn1 + n_v * qn2) / (one_minus * one_minus)
    inner = (hi_v + 1.0) * a - b
    out[valid] = hi_v * log2_r + np.log2(inner)
    return out


def fixed_z_outer_logs(
    *,
    blocks: int,
    spectrum: list[tuple[int, int]],
    zs: list[float],
) -> list[tuple[float, float, float]]:
    out: list[tuple[float, float, float]] = []
    for z in zs:
        local = local_spectrum_log2(spectrum, z)
        global_log = log2_sub_one(blocks * local)
        if global_log != float("-inf"):
            out.append((z, math.log2(z), global_log))
    return out


def outer_bounds_for_h(h_values: np.ndarray, local_logs: list[tuple[float, float, float]]) -> tuple[np.ndarray, np.ndarray]:
    best = np.full(h_values.shape, float("inf"), dtype=np.float64)
    best_z = np.full(h_values.shape, float("nan"), dtype=np.float64)
    h_float = h_values.astype(np.float64)
    for z, log_z, global_log in local_logs:
        candidate = global_log - h_float * log_z
        mask = candidate < best
        best[mask] = candidate[mask]
        best_z[mask] = z
    best[~np.isfinite(best)] = float("-inf")
    return best, best_z


def interval_total(
    *,
    interval: Interval,
    b: int,
    n_total: int,
    distance: int,
    late_blocks: int,
    gaps: list[tuple[int, int]],
    local_logs: list[tuple[float, float, float]],
    entries: list[tuple[int, int, float]],
    turnoff_log2: float,
) -> tuple[float, dict[str, float | int]]:
    h_values = np.arange(interval.start, interval.stop + 1, dtype=np.int64)
    outer, outer_z = outer_bounds_for_h(h_values, local_logs)
    log_c_n_h = log2_binom_range(n_total, interval.start, interval.stop)

    lam = interval.lam
    rho = interval.rho
    m_plus = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
    if not (0.0 < m_plus < 1.0):
        raise ValueError(f"{interval.label}: invalid lambda {lam:g}")
    pref = lam * distance / math.log(2.0)
    log_m = math.log2(m_plus)
    log_rho = math.log2(rho)
    log_q = turnoff_log2 + math.log2(1.0 - m_plus)
    log_g = log2_expm1_pos(b * math.log1p(rho)) + math.log2(m_plus / (1.0 - m_plus))

    terms: list[float] = []
    peak_term = float("-inf")
    peak_h = -1
    peak_r = -1
    peak_gap_min = -1
    peak_branch = ""
    peak_outer_z = float("nan")

    for gap_idx, (gap_min, gap_max) in enumerate(gaps):
        t_min = late_blocks + gap_min - 1
        t_max = late_blocks + gap_max - 1
        n_min = b * t_min
        n_end = b * t_max
        gap_log = math.log2(gap_max - gap_min + 1)
        e0_inner = pref + t_min * log_m

        h_lo = interval.start - b
        h_hi = interval.stop - 1
        log_c_nmin = log2_binom_range(n_min, h_lo, h_hi)
        log_c_nend = log2_binom_range(n_end, h_lo, h_hi)
        tail_factor = log2_one_plus_pow2(
            log2_ht_tail_envelope(H=interval.stop - 1, T=t_min, log_q=log_q)
        )

        for r in range(1, b + 1):
            H = h_values - r
            valid = (H >= 0) & (H <= n_end)
            if not np.any(valid):
                continue
            h_valid = h_values[valid]
            H_valid = H[valid]
            h_index = h_valid - interval.start
            H_index = H_valid - h_lo
            common = (
                outer[valid]
                - log_c_n_h[h_index]
                + log2_binom(b, r)
                + gap_log
                + log_c_nend[H_index]
            )
            e0_term = common + e0_inner

            inner_ege1 = np.full(common.shape, float("-inf"), dtype=np.float64)
            feasible = H_valid <= n_min
            if np.any(feasible):
                H_feas = H_valid[feasible]
                x_min = np.maximum(1, (H_feas + b - 1) // b)
                x_max = np.minimum(H_feas, t_min)
                s_values = log2_arith_geom_sum_vec(log_g, x_min, x_max)
                inner_ege1[feasible] = (
                    turnoff_log2
                    + pref
                    - log_c_nmin[H_index[feasible]]
                    - H_feas.astype(np.float64) * log_rho
                    + s_values
                    + tail_factor
                )
            ege1_term = common + inner_ege1
            # min(0, inner) in the row enumerator means common + min(0, raw).
            raw_inner = np.logaddexp2(e0_inner, inner_ege1)
            row_terms = common + np.minimum(0.0, raw_inner)
            total_here = log2_sum_exp(row_terms)
            terms.append(total_here)

            local_idx = int(np.argmax(row_terms))
            local_peak = float(row_terms[local_idx])
            if local_peak > peak_term:
                peak_term = local_peak
                peak_h = int(h_valid[local_idx])
                peak_r = r
                peak_gap_min = gap_min
                branch = "e0" if float(e0_term[local_idx]) >= float(ege1_term[local_idx]) else "ege1"
                peak_branch = branch if raw_inner[local_idx] <= 0.0 else "cap"
                peak_outer_z = float(outer_z[valid][local_idx])

    total = log2add_many(terms)
    meta: dict[str, float | int] = {
        "peak_h": peak_h,
        "peak_r": peak_r,
        "peak_gap_min": peak_gap_min,
        "peak_term_log2": peak_term,
        "peak_branch": peak_branch,
        "peak_outer_z": peak_outer_z,
        "lambda": lam,
        "rho": rho,
        "log2_g": log_g,
    }
    return total, meta


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
    parser.add_argument("--gap-step", type=int, default=5000)
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
    gaps = [
        (gap_min, min(args.gap_stop, gap_min + args.gap_step - 1))
        for gap_min in range(args.gap_start, args.gap_stop + 1, args.gap_step)
    ]
    d = math.floor(args.distance_delta * args.N)
    local_spectrum = load_local_spectrum(str(args.local_spectrum_csv))
    local_logs = fixed_z_outer_logs(
        blocks=args.outer_blocks,
        spectrum=local_spectrum,
        zs=z_grid(args.z_min, args.z_max, args.z_count),
    )
    entries = split_entries(load_spectrum(args.inner_spectrum), args.block_bits)

    total = float("-inf")
    print("Full-split early endpoint interval certificate")
    print(f"N,{args.N}")
    print(f"delta,{args.distance_delta}")
    print(f"d,{d}")
    print("h_min,h_max,lambda,rho,total_log2,peak_h,peak_r,peak_gap_min,peak_branch,peak_term_log2,log2_g,peak_outer_z")
    for interval in intervals:
        value, meta = interval_total(
            interval=interval,
            b=args.block_bits,
            n_total=args.N,
            distance=d,
            late_blocks=args.late_blocks,
            gaps=gaps,
            local_logs=local_logs,
            entries=entries,
            turnoff_log2=args.turnoff_log2,
        )
        total = log2add(total, value)
        print(
            f"{interval.start},{interval.stop},{interval.lam:.12g},{interval.rho:.12g},"
            f"{value:.6f},{meta['peak_h']},{meta['peak_r']},{meta['peak_gap_min']},"
            f"{meta['peak_branch']},{meta['peak_term_log2']:.6f},{meta['log2_g']:.6f},"
            f"{meta['peak_outer_z']:.12g}"
        )
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
