#!/usr/bin/env python3
"""Budget higher-interpolation savings for the BCH boundary problem.

For normalized h=61 boundary words, degree-30 monic g is determined by 30
agreement points with y^255.  If we count using s>30 points instead, a valid
boundary word is counted binom(h,s) times, while only a subset of the
binom(511,s) possible s-sets are consistent with some monic degree-30 g.

This script prints how many bits of consistency saving among s-sets are needed
to reach the certificate allowance.  One independent F_512 equation is worth
9 bits.
"""

from __future__ import annotations

import argparse
import math

from dense_largek_eval import log2_binom


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=511)
    parser.add_argument("--h", type=int, default=61)
    parser.add_argument("--degree", type=int, default=30)
    parser.add_argument("--allowed-codewords-log2", type=float, default=67.130781)
    parser.add_argument("--s-min", type=int, default=30)
    parser.add_argument("--s-max", type=int, default=45)
    parser.add_argument("--field-bits", type=float, default=9.0)
    args = parser.parse_args()

    orbit = math.log2(args.n)
    print("BCH boundary higher-interpolation budget")
    print(f"n={args.n}, h={args.h}, degree={args.degree}, orbit_log2={orbit:.6f}")
    print(f"allowed_codewords_log2={args.allowed_codewords_log2:.6f}")
    print("s,trivial_sset_bound_log2,excess_bits,needed_independent_Fq_conditions")
    for s in range(args.s_min, args.s_max + 1):
        trivial = log2_binom(args.n, s) - log2_binom(args.h, s) + orbit
        excess = trivial - args.allowed_codewords_log2
        needed = max(0.0, excess / args.field_bits)
        print(f"{s},{trivial:.6f},{excess:.6f},{needed:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
