#!/usr/bin/env python3
"""Verify the effective one-step termination atom for the full-split inner.

The finite prefix certificate uses a termination atom that includes both

  * the local BCH split self-turnoff event q=0, and
  * exact cancellation of a nonzero next-state pattern by the next input block.

For the current tiny prefix, H<=499 input ones remain after the first active
block and the relevant placement buckets have T>=5949 remaining blocks.  For
fixed H and q, the exact-cancellation probability

    C(b(T-1), H-q) / C(bT, H)

is decreasing in T: enlarging the outside universe while keeping H fixed only
decreases the chance that the distinguished block contains exactly the fixed
q-subset and no other selected coordinates.  Hence the finite-prefix supremum
is attained at T=5949 and can be checked by scanning H=0..499.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from analyze_fullsplit_exact_cancellation import exact_hit_log2, load_spectrum, split_q_law
from dense_largek_eval import log2add


ROOT = Path(__file__).resolve().parent


def check_close(name: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise SystemExit(
            f"{name}: expected {expected:.12f}, got {actual:.12f}; "
            f"delta={actual - expected:.6g}"
        )


def check_leq_log2(name: str, actual: float, upper: float, tol: float) -> None:
    if actual > upper + tol:
        raise SystemExit(
            f"{name}: expected <= {upper:.12f}, got {actual:.12f}; "
            f"excess={actual - upper:.6g}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=ROOT / "EBCH128_64.wd")
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--h-min", type=int, default=0)
    parser.add_argument("--h-max", type=int, default=499)
    parser.add_argument("--t-min", type=int, default=5949)
    parser.add_argument("--t-max", type=int, default=17948)
    parser.add_argument("--expected-p0-log2", type=float, default=-64.00457432724919)
    parser.add_argument("--expected-eta-log2", type=float, default=-67.63641076470455)
    parser.add_argument("--expected-computed-pterm-log2", type=float, default=-63.89264922047382)
    parser.add_argument(
        "--pterm-upper-log2",
        type=float,
        default=-63.8926492,
        help="Conservative rounded log2 upper bound used by finite certificate scripts.",
    )
    parser.add_argument("--tolerance", type=float, default=5e-10)
    args = parser.parse_args()

    if args.h_min < 0 or args.h_max < args.h_min:
        raise SystemExit("invalid H range")
    if args.t_min < 1 or args.t_max < args.t_min:
        raise SystemExit("invalid T range")

    b = args.block_bits
    pq = split_q_law(load_spectrum(args.spectrum), b)
    p0_log2 = math.log2(pq[0])

    eta_log2 = float("-inf")
    eta_h = -1
    for H in range(args.h_min, args.h_max + 1):
        value = exact_hit_log2(
            pq=pq,
            b=b,
            remaining_ones=H,
            remaining_blocks=args.t_min,
        )
        if value > eta_log2:
            eta_log2 = value
            eta_h = H

    pterm_log2 = log2add(p0_log2, eta_log2)
    extra_upper_log2 = math.log2(2.0**args.pterm_upper_log2 - 2.0**p0_log2)

    check_close("p0_log2", p0_log2, args.expected_p0_log2, args.tolerance)
    check_close("eta_log2", eta_log2, args.expected_eta_log2, args.tolerance)
    check_close(
        "computed_pterm_log2",
        pterm_log2,
        args.expected_computed_pterm_log2,
        args.tolerance,
    )
    check_leq_log2("computed_pterm_log2", pterm_log2, args.pterm_upper_log2, args.tolerance)

    print("Full-split termination atom verifier")
    print(f"spectrum,{args.spectrum}")
    print(f"block_bits,{b}")
    print(f"finite_prefix_H_range,{args.h_min}..{args.h_max}")
    print(f"finite_prefix_T_range,{args.t_min}..{args.t_max}")
    print("T_reduction,monotone_decreasing_so_T_min_suffices")
    print(f"self_turnoff_p0_log2,{p0_log2:.14f}")
    print(f"eta_sup_H,{eta_h}")
    print(f"eta_sup_T,{args.t_min}")
    print(f"eta_sup_log2,{eta_log2:.14f}")
    print(f"computed_pterm_log2,{pterm_log2:.14f}")
    print(f"certified_pterm_upper_log2,{args.pterm_upper_log2:.14f}")
    print(f"extra_turnoff_for_certified_upper_log2,{extra_upper_log2:.14f}")
    print("status,PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
