#!/usr/bin/env python3
"""Heuristic splitting probability for BCH boundary locator polynomials.

The normalized h=61 boundary family has q^((h-1)/2) possible monic degree-30
polynomials g, giving f(X)=X g(X^2)+1.  If these behaved like random monic
degree-h polynomials, the probability of splitting into h distinct nonzero
roots would be binom(q-1,h)/q^h.  This script compares that heuristic with the
certificate allowance.
"""

from __future__ import annotations

import argparse
import math

from dense_largek_eval import log2_binom


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, default=512)
    parser.add_argument("--h", type=int, default=61)
    parser.add_argument("--allowed-codewords-log2", type=float, default=67.130781)
    args = parser.parse_args()

    n = args.q - 1
    free_coeffs = (args.h - 1) // 2
    field_bits = math.log2(args.q)
    locator_family = free_coeffs * field_bits
    split_prob = log2_binom(n, args.h) - args.h * field_bits
    normalized_split = locator_family + split_prob
    codeword_split = math.log2(n) + normalized_split
    print("BCH boundary random-splitting heuristic")
    print(f"q={args.q}, n={n}, h={args.h}")
    print(f"locator_family_log2={locator_family:.6f}")
    print(f"random_monic_split_probability_log2={split_prob:.6f}")
    print(f"expected_normalized_split_count_log2={normalized_split:.6f}")
    print(f"expected_codeword_count_with_scaling_log2={codeword_split:.6f}")
    print(f"allowed_codewords_log2={args.allowed_codewords_log2:.6f}")
    print(f"heuristic_margin_bits={args.allowed_codewords_log2 - codeword_split:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
