#!/usr/bin/env python3
"""Finite interpolation bounds for the BCH boundary locator family.

For the normalized h=61 boundary family, roots satisfy

    x g(x^2) = 1.

Since squaring is a permutation of F_512^*, this is equivalent to saying that
the monic degree-30 polynomial g agrees with the function y -> y^255 on h
points.  Any 30 agreement points determine g.  Therefore the number of valid
g's is at most binom(511,30)/binom(h,30).

This is rigorous but still too weak for the certificate; it is best understood
as the first finite splitting/list-decoding bound.
"""

from __future__ import annotations

import argparse
import math

from dense_largek_eval import log2_binom


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=511)
    parser.add_argument("--h", type=int, default=61)
    parser.add_argument("--degree", type=int, default=None)
    parser.add_argument("--allowed-codewords-log2", type=float, default=67.130781)
    args = parser.parse_args()

    degree = args.degree if args.degree is not None else (args.h - 1) // 2
    # A monic degree-degree polynomial has degree free coefficients, hence
    # degree agreement points determine it.
    determining_points = degree
    normalized = log2_binom(args.n, determining_points) - log2_binom(args.h, determining_points)
    with_orbit = normalized + math.log2(args.n)
    print("BCH boundary interpolation bound")
    print(f"n={args.n}, h={args.h}, degree={degree}")
    print(f"determining_points={determining_points}")
    print(f"normalized_count_bound_log2={normalized:.6f}")
    print(f"with_scaling_orbit_log2={with_orbit:.6f}")
    print(f"allowed_codewords_log2={args.allowed_codewords_log2:.6f}")
    print(f"excess_over_allowed_bits={with_orbit - args.allowed_codewords_log2:.6f}")
    print(f"rs_agreement_view=degree-{degree} monic polynomials agreeing with y^255 on {args.h} of {args.n} points")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
