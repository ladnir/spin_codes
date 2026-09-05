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
from dataclasses import dataclass
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


@dataclass(frozen=True)
class PlacementRatioCertificate:
    threshold: Fraction
    comparisons: int
    peak_h: int
    peak_num: int
    peak_den: int
    threshold_slack: int
    endpoint_num: int
    endpoint_den: int

    @property
    def peak_to_h(self) -> int:
        return self.peak_h + 1

    @property
    def peak_ratio(self) -> float:
        return ratio_float(self.peak_num, self.peak_den)

    @property
    def peak_log2(self) -> float:
        return math.log2(self.peak_num) - math.log2(self.peak_den)

    @property
    def threshold_gap_log2(self) -> float:
        gap_num = self.threshold_slack
        gap_den = self.threshold.denominator * self.peak_den
        return math.log2(gap_num) - math.log2(gap_den)

    @property
    def endpoint_ratio(self) -> float:
        return ratio_float(self.endpoint_num, self.endpoint_den)

    @property
    def endpoint_log2(self) -> float:
        return math.log2(self.endpoint_num) - math.log2(self.endpoint_den)

    @property
    def endpoint_slack_factor(self) -> float:
        return ratio_float(self.endpoint_num * self.peak_den, self.endpoint_den * self.peak_num)


def certify_placement_ratio(
    *,
    N: int,
    b: int,
    late_blocks: int,
    gap_min: int,
    gap_max: int,
    h_min: int,
    h_max: int,
    threshold_text: str,
) -> PlacementRatioCertificate:
    if h_min < 1 or h_max <= h_min:
        raise ValueError("h_min must be positive and h_max must be larger")
    threshold = decimal_fraction(threshold_text)
    sums = compute_sums(
        b=b,
        late_blocks=late_blocks,
        gap_min=gap_min,
        gap_max=gap_max,
        h_max=h_max,
    )

    max_num = 0
    max_den = 1
    max_h = -1
    comparisons = 0
    for h in range(h_min, h_max):
        numerator = sums[h] * (h + 1)
        denominator = sums[h - 1] * (N - h)
        comparisons += 1
        if numerator * threshold.denominator > denominator * threshold.numerator:
            raise SystemExit(f"placement ratio h={h}->{h + 1} exceeds {threshold_text}")
        if numerator * max_den > max_num * denominator:
            max_num = numerator
            max_den = denominator
            max_h = h

    n_max = b * (late_blocks + gap_max - 1)
    endpoint_num = (n_max - h_min + 1) * (h_min + 1)
    endpoint_den = h_min * (N - h_min)
    threshold_slack = max_den * threshold.numerator - max_num * threshold.denominator
    if threshold_slack <= 0:
        raise SystemExit("placement ratio peak has no positive threshold slack")
    return PlacementRatioCertificate(
        threshold=threshold,
        comparisons=comparisons,
        peak_h=max_h,
        peak_num=max_num,
        peak_den=max_den,
        threshold_slack=threshold_slack,
        endpoint_num=endpoint_num,
        endpoint_den=endpoint_den,
    )


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

    cert = certify_placement_ratio(
        N=args.N,
        b=args.block_bits,
        late_blocks=args.late_blocks,
        gap_min=args.gap_min,
        gap_max=args.gap_max,
        h_min=args.h_min,
        h_max=args.h_max,
        threshold_text=args.max_placement_ratio,
    )

    print("Prefix first-gap placement ratio certificate")
    print(f"N,{args.N}")
    print(f"block_bits,{args.block_bits}")
    print(f"late_blocks,{args.late_blocks}")
    print(f"gap_min,{args.gap_min}")
    print(f"gap_max,{args.gap_max}")
    print(f"h_min,{args.h_min}")
    print(f"h_max,{args.h_max}")
    print(f"placement_ratio_threshold,{args.max_placement_ratio}")
    print(f"placement_ratio_threshold_num,{cert.threshold.numerator}")
    print(f"placement_ratio_threshold_den,{cert.threshold.denominator}")
    print("placement_ratio_exact_cross_multiply_status,PASS")
    print(f"placement_ratio_exact_comparisons,{cert.comparisons}")
    print(f"placement_ratio_peak_h,{cert.peak_h}")
    print(f"placement_ratio_peak_to_h,{cert.peak_to_h}")
    print(f"placement_ratio_peak,{cert.peak_ratio:.12g}")
    print(f"placement_ratio_peak_log2,{cert.peak_log2:.12g}")
    print(f"placement_ratio_peak_num_bits,{cert.peak_num.bit_length()}")
    print(f"placement_ratio_peak_den_bits,{cert.peak_den.bit_length()}")
    print(f"placement_ratio_threshold_slack_bits,{cert.threshold_slack.bit_length()}")
    print(f"placement_ratio_threshold_gap_log2,{cert.threshold_gap_log2:.12g}")
    print(f"endpoint_bound_at_h_min,{cert.endpoint_ratio:.12g}")
    print(f"endpoint_bound_at_h_min_log2,{cert.endpoint_log2:.12g}")
    print(f"endpoint_bound_slack_factor,{cert.endpoint_slack_factor:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
