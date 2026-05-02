#!/usr/bin/env python3
"""Compare exact inner slice probabilities against current analytic envelopes."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from dense_largek_eval import binom_cdf_half_upper, exact_smallw_inner_bound


def parse_enum(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        vals: list[float] = []
        for x in line.split(","):
            x = x.strip()
            if not x:
                continue
            vals.append(float("-inf") if x == "-inf" else float(x))
        if vals:
            rows.append(vals)
    return rows


def exact_slice_prob(rows: list[list[float]], n: int, w: int, cut: int) -> float:
    row = rows[w]
    num = sum((0.0 if v == float("-inf") else 2.0**v) for v in row[: cut + 1])
    return num / math.comb(n, w)


def global_inner_bound(n: int, w: int, sigma: int, delta: float, steps: int = 400) -> float:
    best = 1.0
    off = n * (2.0 ** (-sigma))
    cut = math.floor(delta * n)
    for i in range(1, steps):
        eps = i / steps
        rho = 1.0 - eps
        l = math.floor(rho * n)
        q = binom_cdf_half_upper(l, cut)
        val = (rho**w) + off + q
        if val < best:
            best = val
    return min(1.0, best)


def q_power_model(q_single: float, w: int) -> float:
    return q_single**w


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enum", type=Path, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--sigma", type=int, required=True)
    ap.add_argument("--delta", type=float, required=True)
    ap.add_argument("--w-max", type=int, default=12)
    ap.add_argument("--xi", type=float, default=8.0)
    args = ap.parse_args()

    rows = parse_enum(args.enum)
    cut = math.floor(args.delta * args.n)

    print(f"enum   : {args.enum}")
    print(f"n      : {args.n}")
    print(f"sigma  : {args.sigma}")
    print(f"delta  : {args.delta}")
    print(f"cut    : {cut}")
    print()
    q_single = min(1.0, 2.0 * args.delta)
    print(f"q_single = min(1, 2*delta) = {q_single:.6f}")
    print()
    print("w, exact, run_bound, global_bound, q^w, run/exact, global/exact, q^w/exact")

    for w in range(1, min(args.w_max, len(rows) - 1) + 1):
        exact = exact_slice_prob(rows, args.n, w, cut)
        runb = exact_smallw_inner_bound(args.n, w, args.sigma, args.delta, args.xi)
        glob = global_inner_bound(args.n, w, args.sigma, args.delta)
        qpow = q_power_model(q_single, w)
        rr = (runb / exact) if exact > 0.0 else float("inf")
        gr = (glob / exact) if exact > 0.0 else float("inf")
        qr = (qpow / exact) if exact > 0.0 else float("inf")
        print(
            f"{w}, {exact:.6e}, {runb:.6e}, {glob:.6e}, {qpow:.6e}, {rr:.3e}, {gr:.3e}, {qr:.3e}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
