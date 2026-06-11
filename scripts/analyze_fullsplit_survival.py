#!/usr/bin/env python3
"""Survival-conditioned Chernoff diagnostics for the full-split BCH inner.

The unconditional transfer can be dominated by the q'=0 self-turnoff atom.
This helper separates that atom and reports the Chernoff exponent available on
the class where every nonzero branch keeps q'>0 for T consecutive live blocks.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from import_wd_spectrum import parse_wd
from probe_block_recursive_inner import parse_int_list


def load_spectrum(path: Path) -> list[tuple[int, int]]:
    if path.suffix.lower() == ".wd":
        return parse_wd(path.read_text())
    raise ValueError("only .wd spectra are supported")


def split_entries(spectrum: list[tuple[int, int]], b: int) -> list[tuple[int, int, float]]:
    nonzero_total = sum(c for w, c in spectrum if w > 0)
    out: list[tuple[int, int, float]] = []
    for w, count in spectrum:
        if w <= 0 or count <= 0:
            continue
        p_w = count / nonzero_total
        den = math.comb(2 * b, w)
        lo = max(0, w - b)
        hi = min(b, w)
        for q in range(lo, hi + 1):
            j = w - q
            p = p_w * math.comb(b, q) * math.comb(b, j) / den
            out.append((j, q, p))
    return out


def best_survival_chernoff(
    *,
    entries: list[tuple[int, int, float]],
    live_blocks: int,
    distance: int,
    lambda_min: float,
    lambda_max: float,
    lambda_step: float,
) -> tuple[float, float]:
    best = (float("inf"), lambda_min)
    steps = int(round((lambda_max - lambda_min) / lambda_step))
    for i in range(steps + 1):
        lam = lambda_min + i * lambda_step
        mgf = 0.0
        for j, q, p in entries:
            if q > 0:
                mgf += p * math.exp(-lam * j)
        if mgf <= 0.0:
            continue
        val = (lam * distance + live_blocks * math.log(mgf)) / math.log(2)
        if val < best[0]:
            best = (val, lam)
    return best


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--live-blocks", default="5949,9949,13949")
    parser.add_argument("--lambda-min", type=float, default=0.001)
    parser.add_argument("--lambda-max", type=float, default=2.0)
    parser.add_argument("--lambda-step", type=float, default=0.001)
    args = parser.parse_args()

    d = math.floor(args.distance_delta * args.N)
    entries = split_entries(load_spectrum(args.spectrum), args.block_bits)
    p_self = sum(p for _j, q, p in entries if q == 0)
    p_survive = sum(p for _j, q, p in entries if q > 0)
    mean_j = sum(j * p for j, _q, p in entries)
    mean_j_survive_mass = sum(j * p for j, q, p in entries if q > 0)

    print("Full-split survival-conditioned branch diagnostic")
    print(f"N={args.N}, b={args.block_bits}, d={d}, spectrum={args.spectrum}")
    print(f"self_turnoff_log2,{math.log2(p_self):.6f}")
    print(f"survival_mass_log2,{math.log2(p_survive):.6f}")
    print(f"mean_output,{mean_j:.6f}")
    print(f"mean_output_survival_cond,{mean_j_survive_mass / p_survive:.6f}")
    print("live_blocks,best_log2,best_lambda")
    for live_blocks in parse_int_list(args.live_blocks):
        val, lam = best_survival_chernoff(
            entries=entries,
            live_blocks=live_blocks,
            distance=d,
            lambda_min=args.lambda_min,
            lambda_max=args.lambda_max,
            lambda_step=args.lambda_step,
        )
        print(f"{live_blocks},{val:.6f},{lam:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
