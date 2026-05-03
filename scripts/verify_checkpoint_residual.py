#!/usr/bin/env python3
"""Verify the numeric residual constants for the dense+dense checkpoint.

This script checks the four quantitative constants used in the
``Checkpoint residual split at delta=0.106'' lemma:

* bounded-r low-slot adjacent ratio;
* bounded-r high-slot exponent;
* large-r low-slot ledger maximum;
* the resulting coarse N^4 union slack;
* the low-linear handoff gap on the sampled grid.

It intentionally verifies only the residual constants.  The finite prefix
through H is produced by the global episode-cover evaluator.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from dense_largek_eval import fixedtap_banded_outer_gf_log2
from verify_dense_claims import ETA_CRIT, scan_gap


@dataclass
class LedgerMax:
    value: float
    h: int
    r: int
    c: int
    m_slot: int


def bounded_r_log2_ratio(*, n: int, h: int, r_cap: int, theta: float, z: float) -> float:
    q = (1.0 / z) * ((theta * n - h + 2 * r_cap + 1) / (h - 2 * r_cap)) * ((h + 1) / (n - h))
    return math.log2(q)


def binary_entropy(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


def high_slot_exponent(*, delta: float, theta: float, eta: float, z: float) -> tuple[float, float]:
    """Minimize the high-slot exponent from the bounded-r suppression lemma."""

    def exponent(x: float) -> float:
        return x * (1.0 - binary_entropy(delta / x)) - eta * max(0.0, math.log2(x / z))

    def derivative_high_piece(x: float) -> float:
        return 1.0 + math.log2(1.0 - delta / x) - eta / (x * math.log(2.0))

    candidates = {theta, 1.0}
    if theta <= z <= 1.0:
        candidates.add(z)

    lo = max(theta, z)
    hi = 1.0
    if lo < hi and derivative_high_piece(lo) < 0.0 and derivative_high_piece(hi) > 0.0:
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if derivative_high_piece(mid) < 0.0:
                lo = mid
            else:
                hi = mid
        candidates.add(0.5 * (lo + hi))

    best_x = min(candidates, key=exponent)
    return exponent(best_x), best_x


def large_r_candidate_tail_log2(*, eta: float, offset_s: int, r_cap: int) -> float:
    base = math.e * eta * (2.0 ** (2 - offset_s))
    if base >= r_cap + 2:
        return float("inf")
    first = base / (r_cap + 1)
    denom = 1.0 - base / (r_cap + 2)
    return (r_cap + 1) * math.log2(first) - math.log2(denom)


def ledger_value(
    *,
    log2_w: float,
    n: int,
    offset_s: int,
    z: float,
    h: int,
    r: int,
    c: int,
    m_slot: int,
) -> float:
    return (
        log2_w
        + r * math.log2((math.e**2) * (2.0 ** (2 - offset_s)) * (m_slot + 2 * r + 2) / (r * r))
        + c * math.log2(h / (n - h))
        + (h - c) * math.log2(m_slot / n)
        - h * math.log2(z)
    )


def h_candidates_for_fixed_rc(*, n: int, h_low: int, h_high: int, c: int, m_slot: int, z: float) -> set[int]:
    out = {h_low, h_high}
    slope = math.log2(m_slot / n) - math.log2(z)
    if slope < 0.0 and c > 0:
        # d/dh = c*N/(ln(2)*h*(N-h)) + slope.  On h <= N/4 this is decreasing,
        # so any interior maximum occurs at the lower derivative root.
        target = c * n / (math.log(2.0) * (-slope))
        disc = n * n - 4.0 * target
        if disc >= 0.0:
            root = 0.5 * (n - math.sqrt(disc))
            for h in range(math.floor(root) - 2, math.floor(root) + 4):
                if h_low <= h <= h_high:
                    out.add(h)
    return out


def large_r_low_slot_max(
    *,
    k: int,
    sigma: int,
    offset_s: int,
    h_min: int,
    r_cap: int,
    eta: float,
    theta: float,
    z: float,
) -> LedgerMax:
    n = 2 * k
    h_high_global = math.floor(eta * n)
    m_slot = math.floor(theta * n)
    log2_w = fixedtap_banded_outer_gf_log2(k, k, sigma, z)
    best = LedgerMax(value=float("-inf"), h=-1, r=-1, c=-1, m_slot=m_slot)

    r_high = h_high_global // 2
    for r in range(r_cap + 1, r_high + 1):
        for c in (2 * r, 2 * r + 1):
            h_low = max(h_min, c)
            if h_low > h_high_global:
                continue
            for h in h_candidates_for_fixed_rc(n=n, h_low=h_low, h_high=h_high_global, c=c, m_slot=m_slot, z=z):
                val = ledger_value(
                    log2_w=log2_w,
                    n=n,
                    offset_s=offset_s,
                    z=z,
                    h=h,
                    r=r,
                    c=c,
                    m_slot=m_slot,
                )
                if val > best.value:
                    best = LedgerMax(value=val, h=h, r=r, c=c, m_slot=m_slot)
    return best


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=2**20)
    parser.add_argument("--sigma", type=int, default=25)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--h-min", type=int, default=2000)
    parser.add_argument("--r-cap", type=int, default=20)
    parser.add_argument("--eta", type=float, default=0.25)
    parser.add_argument("--linear-eta-hi", type=float, default=ETA_CRIT)
    parser.add_argument("--theta-slot", type=float, default=0.35)
    parser.add_argument("--z", type=float, default=0.370884)
    parser.add_argument("--run-theta", type=float, default=0.001)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument("--eta-lo", type=float, default=0.0009)
    parser.add_argument("--gap-step", type=float, default=1e-4)
    parser.add_argument("--min-ratio-slack", type=float, default=0.05)
    parser.add_argument("--min-high-slot-exponent", type=float, default=0.018)
    parser.add_argument("--max-low-slot-ledger", type=float, default=-170.0)
    parser.add_argument("--max-union-bound", type=float, default=-85.0)
    parser.add_argument("--min-linear-gap", type=float, default=0.009)
    args = parser.parse_args()

    n = 2 * args.k
    offset_s = args.sigma - math.ceil(math.log2(args.k))

    log2_q = bounded_r_log2_ratio(n=n, h=args.h_min, r_cap=args.r_cap, theta=args.theta_slot, z=args.z)
    high_slot, high_slot_x = high_slot_exponent(
        delta=args.delta,
        theta=args.theta_slot,
        eta=args.eta,
        z=args.z,
    )
    candidate_tail = large_r_candidate_tail_log2(eta=args.eta, offset_s=offset_s, r_cap=args.r_cap)
    ledger = large_r_low_slot_max(
        k=args.k,
        sigma=args.sigma,
        offset_s=offset_s,
        h_min=args.h_min,
        r_cap=args.r_cap,
        eta=args.eta,
        theta=args.theta_slot,
        z=args.z,
    )
    union_log2 = math.log2(3.0) + 4.0 * math.log2(n) + ledger.value
    gap = scan_gap(
        eta_lo=args.eta_lo,
        eta_hi=args.linear_eta_hi,
        theta=args.run_theta,
        delta=args.delta,
        xi=args.xi,
        step=args.gap_step,
    )

    print("Dense+dense checkpoint residual verification")
    print(f"k = {args.k}")
    print(f"N = {n}")
    print(f"sigma = {args.sigma}")
    print(f"offset s = {offset_s}")
    print(f"delta = {args.delta}")
    print(f"H = {args.h_min}")
    print(f"R = {args.r_cap}")
    print(f"theta_slot = {args.theta_slot}")
    print(f"z = {args.z}")
    print()
    print(f"bounded-r low-slot log2 q = {log2_q:.12f}")
    print(f"bounded-r high-slot exponent = {high_slot:.12f} at x={high_slot_x:.12f}")
    print(f"large-r candidate tail log2 = {candidate_tail:.12f}")
    print(
        "large-r low-slot ledger max = "
        f"{ledger.value:.12f} at h={ledger.h}, r={ledger.r}, c={ledger.c}, M={ledger.m_slot}"
    )
    print(f"large-r low-slot after 3*N^4 union log2 = {union_log2:.12f}")
    print()
    print(f"linear handoff eta window = [{args.eta_lo:.12f}, {args.linear_eta_hi:.12f}]")
    print(f"linear handoff worst gap = {gap.worst_gap:.12f} at eta={gap.worst_eta:.12f}")
    print(f"linear handoff outer = {gap.outer_at_worst:.12f}")
    print(f"linear handoff inner = {gap.inner_at_worst:.12f}")

    ok = True
    checks = [
        ("bounded-r ratio", log2_q <= -args.min_ratio_slack),
        ("bounded-r high-slot exponent", high_slot >= args.min_high_slot_exponent),
        ("large-r low-slot ledger", ledger.value <= args.max_low_slot_ledger),
        ("large-r union bound", union_log2 <= args.max_union_bound),
        ("linear handoff gap", gap.worst_gap <= -args.min_linear_gap),
    ]
    for name, passed in checks:
        print(f"{name}: {'PASS' if passed else 'FAIL'}")
        ok = ok and passed
    print(f"STATUS: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
