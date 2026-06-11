#!/usr/bin/env python3
"""Certify the tiny-prefix first-gap placement ratio.

For the dominant prefix ridge we need the placement ratio

    P_{h+1}/P_h = (S_h / S_{h-1}) * (h+1)/(N-h),
    S_H = sum_{g=gap_min}^{gap_max} C(b*(late_blocks+g-1), H).

This helper computes the S_H as exact Python integers and verifies the ratio
against a rational decimal threshold by cross multiplication.  It is deliberately
small and independent of the larger CSV ledger, so it can serve as an audit
object for the placement-slope sublemma.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction


def decimal_fraction(text: str) -> Fraction:
    return Fraction(text)


def compute_sums(*, b: int, late_blocks: int, gap_min: int, gap_max: int, h_max: int) -> list[int]:
    sums = [0] * h_max
    for gap in range(gap_min, gap_max + 1):
        n = b * (late_blocks + gap - 1)
        coeff = 1
        sums[0] += coeff
        for h in range(1, h_max):
            coeff = (coeff * (n - h + 1)) // h
            sums[h] += coeff
    return sums


def ratio_float(num: int, den: int) -> float:
    return math.exp(math.log(num) - math.log(den))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-min", type=int, default=1)
    parser.add_argument("--gap-max", type=int, default=4000)
    parser.add_argument("--h-min", type=int, default=32)
    parser.add_argument("--h-max", type=int, default=500)
    parser.add_argument("--max-placement-ratio", default="0.303594")
    args = parser.parse_args()

    if args.h_min < 1 or args.h_max <= args.h_min:
        raise SystemExit("--h-min must be positive and --h-max must be larger")
    threshold = decimal_fraction(args.max_placement_ratio)
    sums = compute_sums(
        b=args.block_bits,
        late_blocks=args.late_blocks,
        gap_min=args.gap_min,
        gap_max=args.gap_max,
        h_max=args.h_max,
    )

    max_num = 0
    max_den = 1
    max_h = -1
    for h in range(args.h_min, args.h_max):
        numerator = sums[h] * (h + 1)
        denominator = sums[h - 1] * (args.N - h)
        if numerator * threshold.denominator > denominator * threshold.numerator:
            raise SystemExit(
                f"placement ratio h={h}->{h + 1} exceeds {args.max_placement_ratio}"
            )
        if numerator * max_den > max_num * denominator:
            max_num = numerator
            max_den = denominator
            max_h = h

    n_max = args.block_bits * (args.late_blocks + args.gap_max - 1)
    endpoint_num = (n_max - args.h_min + 1) * (args.h_min + 1)
    endpoint_den = args.h_min * (args.N - args.h_min)

    print("Prefix first-gap placement ratio certificate")
    print(f"N,{args.N}")
    print(f"block_bits,{args.block_bits}")
    print(f"late_blocks,{args.late_blocks}")
    print(f"gap_min,{args.gap_min}")
    print(f"gap_max,{args.gap_max}")
    print(f"h_min,{args.h_min}")
    print(f"h_max,{args.h_max}")
    print(f"placement_ratio_threshold,{args.max_placement_ratio}")
    print(f"placement_ratio_peak_h,{max_h}")
    print(f"placement_ratio_peak_to_h,{max_h + 1}")
    print(f"placement_ratio_peak,{ratio_float(max_num, max_den):.12g}")
    print(f"placement_ratio_peak_log2,{math.log2(max_num) - math.log2(max_den):.12g}")
    print(f"endpoint_bound_at_h_min,{ratio_float(endpoint_num, endpoint_den):.12g}")
    print(f"endpoint_bound_at_h_min_log2,{math.log2(endpoint_num) - math.log2(endpoint_den):.12g}")
    print(f"endpoint_bound_slack_factor,{ratio_float(endpoint_num * max_den, endpoint_den * max_num):.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
