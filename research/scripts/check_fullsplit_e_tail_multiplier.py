#!/usr/bin/env python3
"""Audit the fixed-pole e>=2 tail multiplier for the full-split inner.

For fixed lambda, M=M_+(lambda), p=p_term, T remaining blocks, and x occupied
later blocks, the proposed ratio lemma bounds the e>=2 contribution relative
to the fixed-lambda e=1 envelope by

    sum_{e=2}^{x+1} C(x+1,e)/(x+1) * C(T-x+e-1,e-1)
        * (p*(1-M))^(e-1).

This helper maximizes that multiplier over the possible x-range for each H,
and also reports the simpler H,T envelope

    sum_{k>=1} C(H,k) C(T,k)/(k+1) * (p*(1-M))^k.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from bound_fullsplit_episode_gaps import load_spectrum, split_entries
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list
from sum_fullsplit_piecewise_certificate import parse_float_list


def log2_ratio_tail_for_x(*, T: int, x: int, log_q: float, cutoff_bits: float) -> float:
    total = float("-inf")
    # k=e-1.  The endpoint e=x+1 is k=x and is covered by this formula.
    for k in range(1, x + 1):
        # C(x+1,k+1)/(x+1) = C(x,k+1)/(k+1) for k+1<=x,
        # but the direct form below is clear and handles k=x.
        term = (
            log2_binom(x + 1, k + 1)
            - math.log2(x + 1)
            + log2_binom(T - x + k, k)
            + k * log_q
        )
        total = log2add(total, term)
        if term < total - cutoff_bits and k > 4:
            break
    return total


def log2_ht_envelope(*, H: int, T: int, log_q: float, cutoff_bits: float) -> float:
    total = float("-inf")
    k_max = min(H, T)
    for k in range(1, k_max + 1):
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--remaining-blocks", type=int, required=True)
    parser.add_argument("--remaining-ones", required=True)
    parser.add_argument("--turnoff-log2", type=float, default=-62.4078758)
    parser.add_argument("--lambdas", default="0.001,0.003,0.01,0.02,0.05")
    parser.add_argument("--cutoff-bits", type=float, default=120.0)
    args = parser.parse_args()

    entries = split_entries(load_spectrum(args.spectrum), args.block_bits)
    H_values = parse_int_list(args.remaining_ones)
    lambdas = parse_float_list(args.lambdas)
    T = args.remaining_blocks

    print("Full-split e>=2 tail multiplier audit")
    print(f"T,{T}")
    print(f"turnoff_log2,{args.turnoff_log2:.12f}")
    print("H,lambda,M,log2_1_minus_M,sharp_x,sharp_log2,ht_log2,first_order_ht_log2")
    for H in H_values:
        x_min = 0 if H == 0 else max(1, (H + args.block_bits - 1) // args.block_bits)
        x_max = min(H, T)
        for lam in lambdas:
            M = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
            if not (0.0 < M < 1.0):
                continue
            log_one_minus = math.log2(1.0 - M)
            log_q = args.turnoff_log2 + log_one_minus
            sharp = float("-inf")
            sharp_x = -1
            for x in range(x_min, x_max + 1):
                val = log2_ratio_tail_for_x(T=T, x=x, log_q=log_q, cutoff_bits=args.cutoff_bits)
                if val > sharp:
                    sharp = val
                    sharp_x = x
            ht = log2_ht_envelope(H=H, T=T, log_q=log_q, cutoff_bits=args.cutoff_bits)
            first = (
                log2_binom(H, 1)
                + log2_binom(T, 1)
                - 1.0
                + log_q
                if H >= 1 and T >= 1
                else float("-inf")
            )
            print(
                f"{H},{lam:.9g},{M:.12g},{log_one_minus:.6f},"
                f"{sharp_x},{sharp:.6f},{ht:.6f},{first:.6f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
