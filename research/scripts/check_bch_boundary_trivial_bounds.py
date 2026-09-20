#!/usr/bin/env python3
"""Compare trivial locator boundary counts with the BCH spectrum target."""

from __future__ import annotations

import argparse
import math

from dense_largek_eval import log2_binom


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=511)
    parser.add_argument("--h", type=int, default=61)
    parser.add_argument("--parent-dim", type=int, default=259)
    parser.add_argument("--subcode-dim", type=int, default=256)
    parser.add_argument("--slack-bits", type=float, default=56.75)
    parser.add_argument("--field-bits", type=int, default=9)
    args = parser.parse_args()

    redundancy = args.n - args.parent_dim
    free_even_coeffs = (args.h - 1) // 2
    normalized_polys = args.field_bits * free_even_coeffs
    orbit_factor = math.log2(args.n)
    trivial_codewords = normalized_polys + orbit_factor
    randomlike_parent = log2_binom(args.n, args.h) - redundancy
    allowed_parent = log2_binom(args.n, args.h) - (args.n - args.subcode_dim) + args.slack_bits
    print("BCH boundary trivial-bound audit")
    print(f"n={args.n}, h={args.h}, parent_dim={args.parent_dim}, subcode_dim={args.subcode_dim}")
    print(f"randomlike_parent_log2_Ah={randomlike_parent:.6f}")
    print(f"allowed_parent_log2_Ah_for_subcode_certificate={allowed_parent:.6f}")
    print(f"normalized_locator_polynomial_log2_count={normalized_polys:.6f}")
    print(f"with_scaling_orbit_log2_count={trivial_codewords:.6f}")
    print(f"excess_over_allowed_bits={trivial_codewords - allowed_parent:.6f}")
    print(f"needed_splitting_savings_bits={max(0.0, trivial_codewords - allowed_parent):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
