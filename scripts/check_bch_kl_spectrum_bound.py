#!/usr/bin/env python3
"""Check the finite usefulness of the Krasikov--Litsyn BCH spectrum bound.

Krasikov--Litsyn give an asymptotic binomial approximation for primitive BCH
spectra.  For our native BCH candidate, the relevant parent is length 511 with
designed distance 61 and dimension 259.  This script evaluates the explicit
error factor from their Theorem 3 at the weights controlling our first moment.

The point of this checker is diagnostic: at n=511 the asymptotic error term is
too large to use directly as a finite certificate.
"""

from __future__ import annotations

import argparse
import math

from bch_candidate_params import bch_dimension


def parse_weights(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            vals = [int(x) for x in part.split(":")]
            if len(vals) == 2:
                lo, hi = vals
                step = 1
            elif len(vals) == 3:
                lo, hi, step = vals
            else:
                raise ValueError(f"bad weight range {part!r}")
            out.extend(range(lo, hi + 1, step))
        else:
            out.append(int(part))
    return out


def log2_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) / math.log(2.0) - math.lgamma(k + 1) / math.log(2.0) - math.lgamma(n - k + 1) / math.log(2.0)


def kl_log2_error_bound(n: int, t: int, weight: int) -> float:
    ell = (weight + 1) // 2
    # Theorem 3 error factor:
    # sqrt(2) (2 ell)^t exp(2(t-1)^2 - ell(n-2ell)/n) n^(ell-t) (1+o(1)).
    # We report the displayed finite part, without claiming control of o(1).
    return (
        0.5
        + t * math.log2(2 * ell)
        + (2 * (t - 1) ** 2 - ell * (n - 2 * ell) / n) / math.log(2.0)
        + (ell - t) * math.log2(n)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=9)
    parser.add_argument("--designed-distance", type=int, default=61)
    parser.add_argument(
        "--kl-t",
        type=int,
        default=None,
        help="The t parameter in the KL denominator (n+1)^t. Default: generator_degree/m.",
    )
    parser.add_argument("--subcode-dim", type=int, default=256)
    parser.add_argument("--weights", default="61,70,100,150,220")
    args = parser.parse_args()

    n, dim, degree = bch_dimension(args.m, args.designed_distance)
    t = args.kl_t if args.kl_t is not None else degree // args.m
    weights = parse_weights(args.weights)
    print("Native BCH spectrum-bound diagnostic")
    print(f"m = {args.m}")
    print(f"n = {n}")
    print(f"designed distance = {args.designed_distance}")
    print(f"parent dimension = {dim}")
    print(f"generator degree = {degree}")
    print(f"KL t = {t}")
    print(f"subcode dimension = {args.subcode_dim}")
    print()
    print("weight,ell,randomlike_subcode_log2_Ah,parent_main_log2_Ah,kl_log2_error_bound")
    for weight in weights:
        ell = (weight + 1) // 2
        randomlike = log2_binom(n, weight) + args.subcode_dim - n
        parent_main = log2_binom(n, weight) - t * math.log2(n + 1)
        kl_error = kl_log2_error_bound(n, t, weight)
        print(f"{weight},{ell},{randomlike:.6f},{parent_main:.6f},{kl_error:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
