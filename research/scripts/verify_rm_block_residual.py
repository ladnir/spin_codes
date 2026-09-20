#!/usr/bin/env python3
"""Verify residual constants for the RM(4,9) block-outer checkpoint.

This is the RM-block analogue of the dense-band residual checker.  The outer
enters only through the product generating function

    W_out(z) = W_RM(4,9)(z)^4096 - 1.

The checked constants are designed to support a finite-prefix theorem with
H=500: bounded-r low-slot geometric decay, bounded-r high-slot exponent, a
large-r low-slot ledger, and the linear-window handoff.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from block_outer_upgrade_probe import load_local_spectrum, local_spectrum_log2
from dense_largek_eval import log2add
from verify_checkpoint_residual import (
    bounded_r_log2_ratio,
    high_slot_exponent,
    h_candidates_for_fixed_rc,
    ledger_value,
)
from verify_dense_claims import ETA_CRIT, scan_gap


def rm_product_log2_w(*, spectrum_path: Path, blocks: int, z: float) -> float:
    local = local_spectrum_log2(load_local_spectrum(str(spectrum_path)), z)
    all_log = blocks * local
    if all_log == float("-inf"):
        return float("-inf")
    if all_log <= 1e-10:
        val = math.expm1(all_log * math.log(2.0))
        return math.log2(val) if val > 0.0 else float("-inf")
    return all_log + math.log2(1.0 - 2.0 ** (-all_log))


def large_r_low_slot_max_rm(
    *,
    log2_w: float,
    n: int,
    offset_s: int,
    h_min: int,
    r_cap: int,
    eta: float,
    theta: float,
    z: float,
) -> tuple[float, int, int, int, int]:
    h_high_global = math.floor(eta * n)
    m_slot = math.floor(theta * n)
    best_value = float("-inf")
    best_h = -1
    best_r = -1
    best_c = -1

    for r in range(r_cap + 1, h_high_global // 2 + 1):
        for c in (2 * r, 2 * r + 1):
            h_low = max(h_min, c)
            if h_low > h_high_global:
                continue
            for h in h_candidates_for_fixed_rc(n=n, h_low=h_low, h_high=h_high_global, c=c, m_slot=m_slot, z=z):
                value = ledger_value(
                    log2_w=log2_w,
                    n=n,
                    offset_s=offset_s,
                    z=z,
                    h=h,
                    r=r,
                    c=c,
                    m_slot=m_slot,
                )
                if value > best_value:
                    best_value = value
                    best_h = h
                    best_r = r
                    best_c = c
    return best_value, best_h, best_r, best_c, m_slot


def check(name: str, passed: bool, detail: str) -> bool:
    print(f"{name}: {'PASS' if passed else 'FAIL'} ({detail})")
    return passed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("rm512_256_spectrum.csv"))
    parser.add_argument("--k", type=int, default=2**20)
    parser.add_argument("--n", type=int, default=2**21)
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--sigma", type=int, default=32)
    parser.add_argument("--delta", type=float, default=0.09)
    parser.add_argument("--h-min", type=int, default=500)
    parser.add_argument("--r-cap", type=int, default=12)
    parser.add_argument("--eta", type=float, default=0.25)
    parser.add_argument("--linear-eta-hi", type=float, default=ETA_CRIT)
    parser.add_argument("--theta-slot", type=float, default=0.35)
    parser.add_argument("--z", type=float, default=0.4)
    parser.add_argument("--run-theta", type=float, default=0.001)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument("--gap-step", type=float, default=1e-4)
    parser.add_argument("--min-ratio-slack", type=float, default=0.10)
    parser.add_argument("--min-high-slot-exponent", type=float, default=0.05)
    parser.add_argument("--max-low-slot-ledger", type=float, default=-250.0)
    parser.add_argument("--max-union-bound", type=float, default=-160.0)
    parser.add_argument("--min-linear-gap", type=float, default=0.002)
    args = parser.parse_args()

    offset_s = args.sigma - math.ceil(math.log2(args.k))
    eta_lo = args.h_min / args.n
    log2_w = rm_product_log2_w(spectrum_path=args.spectrum, blocks=args.blocks, z=args.z)

    log2_q = bounded_r_log2_ratio(
        n=args.n,
        h=args.h_min,
        r_cap=args.r_cap,
        theta=args.theta_slot,
        z=args.z,
    )
    high_slot, high_slot_x = high_slot_exponent(
        delta=args.delta,
        theta=args.theta_slot,
        eta=args.eta,
        z=args.z,
    )
    low_value, low_h, low_r, low_c, m_slot = large_r_low_slot_max_rm(
        log2_w=log2_w,
        n=args.n,
        offset_s=offset_s,
        h_min=args.h_min,
        r_cap=args.r_cap,
        eta=args.eta,
        theta=args.theta_slot,
        z=args.z,
    )
    low_union = math.log2(3.0) + 4.0 * math.log2(args.n) + low_value
    gap = scan_gap(
        eta_lo=eta_lo,
        eta_hi=args.linear_eta_hi,
        theta=args.run_theta,
        delta=args.delta,
        xi=args.xi,
        step=args.gap_step,
    )

    print("RM(4,9) block-outer residual verification")
    print(f"k = {args.k}")
    print(f"N = {args.n}")
    print(f"sigma = {args.sigma}")
    print(f"offset s = {offset_s}")
    print(f"delta = {args.delta}")
    print(f"H = {args.h_min}")
    print(f"R = {args.r_cap}")
    print(f"z = {args.z}")
    print(f"log2 W_out(z) = {log2_w:.12f}")
    print(f"bounded-r low-slot log2 q = {log2_q:.12f}")
    print(f"bounded-r high-slot exponent = {high_slot:.12f} at x={high_slot_x:.12f}")
    print(
        "large-r low-slot ledger max = "
        f"{low_value:.12f} at h={low_h}, r={low_r}, c={low_c}, M={m_slot}"
    )
    print(f"large-r low-slot after 3*N^4 union log2 = {low_union:.12f}")
    print(f"linear handoff eta window = [{eta_lo:.12f}, {args.linear_eta_hi:.12f}]")
    print(f"linear handoff worst gap = {gap.worst_gap:.12f} at eta={gap.worst_eta:.12f}")
    print()

    ok = True
    ok &= check("outer generating value", log2_w < 40.0, f"{log2_w:.6f} < 40")
    ok &= check(
        "bounded-r ratio",
        log2_q <= -args.min_ratio_slack,
        f"{log2_q:.6f} <= {-args.min_ratio_slack:.6f}",
    )
    ok &= check(
        "bounded-r high-slot exponent",
        high_slot >= args.min_high_slot_exponent,
        f"{high_slot:.6f} >= {args.min_high_slot_exponent:.6f}",
    )
    ok &= check(
        "large-r low-slot ledger",
        low_value <= args.max_low_slot_ledger,
        f"{low_value:.6f} <= {args.max_low_slot_ledger:.6f}",
    )
    ok &= check(
        "large-r union bound",
        low_union <= args.max_union_bound,
        f"{low_union:.6f} <= {args.max_union_bound:.6f}",
    )
    ok &= check(
        "linear handoff gap",
        gap.worst_gap <= -args.min_linear_gap,
        f"{gap.worst_gap:.6f} <= {-args.min_linear_gap:.6f}",
    )
    print(f"STATUS: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
