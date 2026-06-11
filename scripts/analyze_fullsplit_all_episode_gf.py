#!/usr/bin/env python3
"""Evaluate the all-episode generating-function bound for the full-split inner.

This implements Lemma fullsplit-all-episode-gf-bound from innerDense.tex.  It
does not truncate the selected-gap count e; instead it uses the Cauchy bound

  [t^(T-x)] (1/(1-t) + p_term/(1-t/M))^(x+1)

with a = alpha*M, 0 < alpha < 1.  The default path also applies a second
Cauchy bound to the occupancy coefficient [z^H]((1+z)^b-1)^x, which is the
version intended for the theorem interface.  The exact coefficient path is kept
for small audit cases.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from bound_fullsplit_episode_gaps import load_spectrum, split_entries
from dense_largek_eval import log2_binom, log2add
from fast_fullsplit_episode_e01 import precompute_nonempty_block_logcoeff
from probe_block_recursive_inner import parse_int_list


def parse_float_list(text: str) -> list[float]:
    out: list[float] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            vals = [float(x) for x in part.split(":")]
            if len(vals) == 3:
                lo, hi, step = vals
                n = int(math.floor((hi - lo) / step + 1e-12))
                out.extend(lo + i * step for i in range(n + 1))
            else:
                raise ValueError(f"bad float range {part!r}; use lo:hi:step")
        else:
            out.append(float(part))
    return out


def default_alphas() -> list[float]:
    vals = [0.25, 0.5, 0.75, 0.9, 0.95, 0.98, 0.99]
    vals.extend(1.0 - 10.0 ** (-k / 2) for k in range(3, 17))
    vals.extend(1.0 - 2.0 ** (-k) for k in range(8, 25, 2))
    return sorted({x for x in vals if 0.0 < x < 1.0})


def default_rhos() -> list[float]:
    vals = [10.0 ** (k / 4) for k in range(-24, 9)]
    vals.extend([0.002, 0.003, 0.005, 0.007, 0.02, 0.03, 0.05, 0.07])
    return sorted({x for x in vals if x > 0.0})


def log2_expm1(x: float) -> float:
    """Return log2(exp(x)-1), stable for the small rho regime."""

    if x <= 0.0:
        return float("-inf")
    if x < 50.0:
        return math.log(math.expm1(x), 2.0)
    return (x + math.log1p(-math.exp(-x))) / math.log(2.0)


def log2_geom_range_from_logg(log_g: float, lo: int, hi: int) -> float:
    """Return log2(sum_{x=lo}^hi g^x) from log2(g)."""

    if lo > hi:
        return float("-inf")
    n = hi - lo + 1
    if abs(log_g) < 1e-12:
        return math.log2(n)

    ln_g = log_g * math.log(2.0)
    if log_g < 0.0:
        tail = 0.0 if n * ln_g < -745.0 else math.exp(n * ln_g)
        return (lo * ln_g + math.log1p(-tail) - math.log1p(-math.exp(ln_g))) / math.log(2.0)

    tail = 0.0 if -n * ln_g < -745.0 else math.exp(-n * ln_g)
    return (hi * ln_g + math.log1p(-tail) - math.log1p(-math.exp(-ln_g))) / math.log(2.0)


def all_episode_exact_occupancy_log2(
    *,
    b: int,
    T: int,
    H: int,
    distance: int,
    term_log2: float,
    entries: list[tuple[int, int, float]],
    coeff_log2: list[list[float]],
    lam: float,
    alpha: float,
) -> float:
    m = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
    if not (0.0 < m < 1.0):
        return float("inf")
    a = alpha * m
    if not (0.0 < a < m):
        return float("inf")
    p_term = 2.0 ** term_log2
    base = lam * distance / math.log(2.0) + T * math.log2(m) - log2_binom(b * T, H)
    gap_factor = (1.0 / (1.0 - a)) + (p_term / (1.0 - a / m))
    if gap_factor <= 0.0:
        return float("inf")
    log_gap = math.log2(gap_factor)
    log_a = math.log2(a)

    total = float("-inf")
    x_min = 0 if H == 0 else max(1, (H + b - 1) // b)
    x_max = min(H, T)
    for x in range(x_min, x_max + 1):
        cx = coeff_log2[x][H]
        if cx == float("-inf"):
            continue
        term = cx + base - (T - x) * log_a + (x + 1) * log_gap
        total = log2add(total, term)
    return min(0.0, total)


def all_episode_cauchy_occupancy_log2(
    *,
    b: int,
    T: int,
    H: int,
    distance: int,
    term_log2: float,
    entries: list[tuple[int, int, float]],
    lam: float,
    alpha: float,
    rho: float,
) -> float:
    m = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
    if not (0.0 < m < 1.0):
        return float("inf")
    a = alpha * m
    if not (0.0 < a < m):
        return float("inf")
    p_term = 2.0 ** term_log2
    gap_factor = (1.0 / (1.0 - a)) + (p_term / (1.0 - a / m))
    if gap_factor <= 0.0:
        return float("inf")

    x_min = 0 if H == 0 else max(1, (H + b - 1) // b)
    x_max = min(H, T)
    if x_min > x_max:
        return float("-inf")

    log_a = math.log2(a)
    log_gap = math.log2(gap_factor)
    log_a_gap = log_a + log_gap
    log_block_poly = log2_expm1(b * math.log1p(rho))
    log_geom = log2_geom_range_from_logg(log_a_gap + log_block_poly, x_min, x_max)
    total = (
        lam * distance / math.log(2.0)
        + T * math.log2(m)
        - log2_binom(b * T, H)
        + log_gap
        - T * log_a
        - H * math.log2(rho)
        + log_geom
    )
    return min(0.0, total)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--remaining-blocks", type=int, required=True)
    parser.add_argument("--remaining-ones", required=True)
    parser.add_argument("--turnoff-log2", type=float, default=-63.8926492)
    parser.add_argument("--lambda-min", type=float, default=0.001)
    parser.add_argument("--lambda-max", type=float, default=2.0)
    parser.add_argument("--lambda-step", type=float, default=0.01)
    parser.add_argument("--alphas", help="Comma list or lo:hi:step alpha fractions for a=alpha*M.")
    parser.add_argument("--rhos", help="Comma list or lo:hi:step positive rho values for the occupancy Cauchy bound.")
    parser.add_argument(
        "--occupancy-bound",
        choices=("cauchy", "exact"),
        default="cauchy",
        help="Use a second Cauchy bound over [z^H] coefficients, or exact coefficient tables.",
    )
    args = parser.parse_args()

    b = args.block_bits
    T = args.remaining_blocks
    H_values = parse_int_list(args.remaining_ones)
    d = math.floor(args.distance_delta * args.N)
    entries = split_entries(load_spectrum(args.spectrum), b)
    h_max = max(H_values)
    coeff_log2 = None
    if args.occupancy_bound == "exact":
        coeff_log2 = precompute_nonempty_block_logcoeff(b=b, h_max=h_max, x_max=min(h_max, T))
    lambdas = [
        args.lambda_min + i * args.lambda_step
        for i in range(int(round((args.lambda_max - args.lambda_min) / args.lambda_step)) + 1)
    ]
    alphas = parse_float_list(args.alphas) if args.alphas else default_alphas()
    rhos = parse_float_list(args.rhos) if args.rhos else default_rhos()

    print("Full-split all-episode GF bound")
    print(f"N={args.N}, b={b}, T={T}, d={d}, turnoff_log2={args.turnoff_log2:.6f}")
    print(f"occupancy_bound,{args.occupancy_bound}")
    print("H,bound_log2,best_lambda,best_alpha,best_rho")
    for H in H_values:
        best = float("inf")
        best_lam = float("nan")
        best_alpha = float("nan")
        best_rho = float("nan")
        for lam in lambdas:
            for alpha in alphas:
                if args.occupancy_bound == "exact":
                    assert coeff_log2 is not None
                    val = all_episode_exact_occupancy_log2(
                        b=b,
                        T=T,
                        H=H,
                        distance=d,
                        term_log2=args.turnoff_log2,
                        entries=entries,
                        coeff_log2=coeff_log2,
                        lam=lam,
                        alpha=alpha,
                    )
                    if val < best:
                        best = val
                        best_lam = lam
                        best_alpha = alpha
                        best_rho = float("nan")
                    continue
                for rho in rhos:
                    val = all_episode_cauchy_occupancy_log2(
                        b=b,
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
        print(f"{H},{best:.6f},{best_lam:.6f},{best_alpha:.12f},{best_rho:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
