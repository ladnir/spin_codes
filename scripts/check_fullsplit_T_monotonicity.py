#!/usr/bin/env python3
"""Audit T-monotonicity for the full-split fixed-pole envelope.

For fixed H and pole (lambda,rho), the e>=1 part of the theorem-facing
envelope has the form

    const * S_T(H;lambda,rho) / binom(bT,H) * (1 + R_{H,T}(lambda)).

This helper checks the one-step slack

    log2 binom(b(T+1),H)/binom(bT,H)
      - log2 S_{T+1}/S_T
      - log2 (1+R_{H,T+1})/(1+R_{H,T}).

Positive slack means the fixed-pole envelope is nonincreasing from T to T+1.
The default intervals are the current RM/full-split certificate ranges.

With --sufficient-reduction, the helper instead checks the two-case endpoint
reduction used in the manuscript: T>=H has no new occupancy term, while T<H
uses H=max(H_min,T+1) and the bound

    S_{T+1}/S_T <= 1 + ((T+2)/(T+1))*G.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

from bound_fullsplit_episode_gaps import load_spectrum, split_entries
from dense_largek_eval import log2_binom
from probe_block_recursive_inner import parse_int_list
from sum_fullsplit_piecewise_certificate import (
    log2_arith_geom_sum,
    log2_ht_tail_envelope,
    log2_one_plus_pow2,
)


@dataclass(frozen=True)
class Interval:
    h_min: int
    h_max: int
    lam: float
    rho: float


DEFAULT_INTERVALS = [
    Interval(2001, 7858, 0.02, 0.01),
    Interval(7859, 20550, 0.05, 0.03),
    Interval(20551, 75000, 0.2, 0.1),
    Interval(75001, 250000, 0.5, 0.3),
    Interval(250001, 350000, 1.2, 0.5),
    Interval(350001, 650000, 1.2, 1.0),
    Interval(650001, 725000, 1.2, 2.0),
    Interval(725001, 950000, 1.2, 3.0),
    Interval(950001, 1148736, 1.2, 3.0),
]


def parse_intervals(text: str) -> list[Interval]:
    out: list[Interval] = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        fields = [x.strip() for x in part.split(",")]
        if len(fields) != 3:
            raise ValueError(f"bad interval {part!r}; expected h0:h1,lambda,rho")
        h0, h1 = [int(x) for x in fields[0].split(":")]
        out.append(Interval(h0, h1, float(fields[1]), float(fields[2])))
    return out


def arith_geom_log_ratio(*, log_g: float, x_min: int, t: int, h: int) -> float:
    x_max = min(h, t)
    x_next = min(h, t + 1)
    if x_next == x_max:
        return 0.0
    log_s = log2_arith_geom_sum(log_g, x_min, x_max)
    log_s_next = log2_arith_geom_sum(log_g, x_min, x_next)
    return log_s_next - log_s


def t_values_for_bucket(*, t_min: int, t_max: int, h: int, b: int, mode: str, samples: int) -> list[int]:
    lo = max(t_min, (h + b - 1) // b)
    hi = t_max - 1
    if lo > hi:
        return []
    if mode == "endpoints":
        return sorted(set([lo, hi]))
    if mode == "sample":
        if samples <= 2:
            return sorted(set([lo, hi]))
        vals = {lo, hi}
        for i in range(1, samples - 1):
            vals.add(lo + round(i * (hi - lo) / (samples - 1)))
        return sorted(vals)
    if mode == "all":
        return list(range(lo, hi + 1))
    raise ValueError(f"unknown t mode {mode!r}")


def denominator_log_ratio(*, b: int, H: int, T: int) -> float:
    return log2_binom(b * (T + 1), H) - log2_binom(b * T, H)


def print_sufficient_reduction(
    *,
    intervals: list[Interval],
    buckets: list[tuple[int, int]],
    entries: list[tuple[int, int, float]],
    b: int,
    turnoff_log2: float,
) -> None:
    print("Full-split T-monotonicity sufficient endpoint reduction", flush=True)
    print(
        "h_min,h_max,lambda,rho,bucket_T_min,bucket_T_max,"
        "no_new_slack,new_slack,min_slack,H_new,T_new",
        flush=True,
    )
    for interval in intervals:
        M = sum(p * math.exp(-interval.lam * j) for j, q, p in entries if q > 0)
        G = math.expm1(b * math.log1p(interval.rho)) * M / (1.0 - M)
        log_q = turnoff_log2 + math.log2(1.0 - M)
        H_min = max(1, interval.h_min - b)
        H_max = interval.h_max - 1
        for t_min, t_max in buckets:
            log_tail_bound = log2_one_plus_pow2(log2_ht_tail_envelope(H=H_max, T=t_max, log_q=log_q))

            no_new_slack = float("inf")
            t_den = t_max - 1
            if H_min <= min(H_max, t_den):
                no_new_slack = denominator_log_ratio(b=b, H=H_min, T=t_den) - log_tail_bound

            new_slack = float("inf")
            new_arg: tuple[int, int] | None = None
            for T in range(t_min, t_max):
                H = max(H_min, T + 1)
                if H > H_max or H > b * T:
                    continue
                log_s_bound = math.log2(1.0 + ((T + 2) / (T + 1)) * G)
                slack = denominator_log_ratio(b=b, H=H, T=T) - log_s_bound - log_tail_bound
                if slack < new_slack:
                    new_slack = slack
                    new_arg = (H, T)

            min_slack = min(no_new_slack, new_slack)
            H_new = "" if new_arg is None else str(new_arg[0])
            T_new = "" if new_arg is None else str(new_arg[1])
            print(
                f"{interval.h_min},{interval.h_max},{interval.lam:g},{interval.rho:g},"
                f"{t_min},{t_max},{no_new_slack:.6f},{new_slack:.6f},"
                f"{min_slack:.6f},{H_new},{T_new}",
                flush=True,
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--turnoff-log2", type=float, default=-63.8926492)
    parser.add_argument("--intervals", help="Semicolon-separated h0:h1,lambda,rho entries.")
    parser.add_argument("--r-values", default="1:64")
    parser.add_argument("--buckets", default="4001:8000,8001:12000")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--h-step", type=int, default=100)
    parser.add_argument("--t-mode", choices=("endpoints", "sample", "all"), default="endpoints")
    parser.add_argument("--t-samples", type=int, default=9)
    parser.add_argument(
        "--sufficient-reduction",
        action="store_true",
        help="Check the two-case endpoint sufficient condition used in the manuscript.",
    )
    args = parser.parse_args()

    intervals = parse_intervals(args.intervals) if args.intervals else DEFAULT_INTERVALS
    r_values = parse_int_list(args.r_values)
    buckets: list[tuple[int, int]] = []
    for part in args.buckets.split(","):
        lo, hi = [int(x) for x in part.split(":")]
        buckets.append((args.late_blocks + lo - 1, args.late_blocks + hi - 1))

    entries = split_entries(load_spectrum(args.spectrum), args.block_bits)

    if args.sufficient_reduction:
        print_sufficient_reduction(
            intervals=intervals,
            buckets=buckets,
            entries=entries,
            b=args.block_bits,
            turnoff_log2=args.turnoff_log2,
        )
        return 0

    print("Full-split T-monotonicity audit", flush=True)
    print("h_min,h_max,lambda,rho,bucket_T_min,bucket_T_max,min_slack,H,T,logD,logS,logTail", flush=True)
    for interval in intervals:
        M = sum(p * math.exp(-interval.lam * j) for j, q, p in entries if q > 0)
        A = math.expm1(args.block_bits * math.log1p(interval.rho))
        log_g = math.log2(A * M / (1.0 - M))
        log_q = args.turnoff_log2 + math.log2(1.0 - M)

        h_values = list(range(interval.h_min, interval.h_max + 1, args.h_step))
        if h_values[-1] != interval.h_max:
            h_values.append(interval.h_max)

        for t_min, t_max in buckets:
            worst = (float("inf"), None)
            for h in h_values:
                for r in r_values:
                    if r < 1 or r > min(h, args.block_bits):
                        continue
                    H = h - r
                    x_min = max(1, (H + args.block_bits - 1) // args.block_bits)
                    for T in t_values_for_bucket(
                        t_min=t_min,
                        t_max=t_max,
                        h=H,
                        b=args.block_bits,
                        mode=args.t_mode,
                        samples=args.t_samples,
                    ):
                        log_d = denominator_log_ratio(b=args.block_bits, H=H, T=T)
                        log_s = arith_geom_log_ratio(log_g=log_g, x_min=x_min, t=T, h=H)
                        log_tail = log2_one_plus_pow2(log2_ht_tail_envelope(H=H, T=T + 1, log_q=log_q))
                        log_tail -= log2_one_plus_pow2(log2_ht_tail_envelope(H=H, T=T, log_q=log_q))
                        slack = log_d - log_s - log_tail
                        if slack < worst[0]:
                            worst = (slack, (H, T, log_d, log_s, log_tail))
            if worst[1] is None:
                print(
                    f"{interval.h_min},{interval.h_max},{interval.lam:g},{interval.rho:g},"
                    f"{t_min},{t_max},inf,,,,,"
                )
                continue
            H, T, log_d, log_s, log_tail = worst[1]
            print(
                f"{interval.h_min},{interval.h_max},{interval.lam:g},{interval.rho:g},"
                f"{t_min},{t_max},{worst[0]:.6f},{H},{T},"
                f"{log_d:.6f},{log_s:.6f},{log_tail:.12g}"
                ,
                flush=True,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
