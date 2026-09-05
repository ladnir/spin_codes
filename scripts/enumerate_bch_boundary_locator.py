#!/usr/bin/env python3
"""Enumerate BCH boundary words via the locator-polynomial family.

For odd boundary weight h=Delta, the BCH locator reduction gives

    f(X) = X g(X^2) + c,

where g is monic of degree (h-1)/2 and c is the first free power sum / product
term.  When gcd(h, 2^m-1)=1, scaling normalizes c=1.  Otherwise c ranges over
representatives of the quotient F_q^* / (F_q^*)^h.

This enumerator is intended for small m where all g can still be enumerated.
It reports the number of boundary codewords after restoring scaling orbits and
compares naturally against exact spectra from exact_bch_spectrum_small.py.
"""

from __future__ import annotations

import argparse
import math

from check_bch_boundary_smallfield import find_primitive_poly, gf_mul
from check_bch_newton_boundary import defining_set


def gf_pow(a: int, e: int, poly: int, m: int) -> int:
    out = 1
    base = a
    while e:
        if e & 1:
            out = gf_mul(out, base, poly, m)
        e >>= 1
        if e:
            base = gf_mul(base, base, poly, m)
    return out


def eval_g(coeffs: list[int], y: int, poly: int, m: int) -> int:
    acc = 1
    for c in reversed(coeffs):
        acc = gf_mul(acc, y, poly, m) ^ c
    return acc


def constant_reps(q: int, h: int, poly: int, m: int) -> list[int]:
    n = q - 1
    alpha = 2
    step = math.gcd(h, n)
    return [gf_pow(alpha, i, poly, m) for i in range(step)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--delta", type=int, default=7)
    parser.add_argument("--h", type=int, default=None)
    parser.add_argument("--max-family-log2", type=float, default=32.0)
    args = parser.parse_args()

    q = 1 << args.m
    n = q - 1
    h = args.h if args.h is not None else args.delta
    if h % 2 != 1:
        raise ValueError("boundary h must be odd")
    d = (h - 1) // 2
    family_log2 = args.m * d
    if family_log2 > args.max_family_log2:
        raise SystemExit(f"family size 2^{family_log2} exceeds --max-family-log2 {args.max_family_log2}")
    poly = find_primitive_poly(args.m)
    field = list(range(1, q))
    constants = constant_reps(q, h, poly, args.m)
    zeros = defining_set(n, args.delta)
    if any(i not in zeros for i in range(1, args.delta)):
        raise AssertionError("defining set closure sanity failed")

    normalized_hits = 0
    root_hist: dict[int, int] = {}
    total_family = q**d
    for c in constants:
        for idx in range(total_family):
            x = idx
            coeffs = []
            for _ in range(d):
                coeffs.append(x & (q - 1))
                x >>= args.m
            roots = []
            for xval in field:
                y = gf_mul(xval, xval, poly, args.m)
                if gf_mul(xval, eval_g(coeffs, y, poly, args.m), poly, args.m) ^ c == 0:
                    roots.append(xval)
            root_hist[len(roots)] = root_hist.get(len(roots), 0) + 1
            if len(roots) == h:
                normalized_hits += 1

    scaling_classes = math.gcd(h, n)
    orbit_size = n // scaling_classes
    boundary_codewords = normalized_hits * orbit_size
    print("BCH boundary locator enumeration")
    print(f"m={args.m}, q={q}, n={n}, delta={args.delta}, h={h}, degree={d}, primitive_poly=0x{poly:x}")
    print(f"family_log2={family_log2:.6f}, constants={len(constants)}, scaling_classes={scaling_classes}")
    print(f"normalized_hits={normalized_hits}")
    print(f"orbit_size={orbit_size}")
    print(f"boundary_codewords={boundary_codewords}")
    print("root_hist=" + ",".join(f"{k}:{root_hist[k]}" for k in sorted(root_hist)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
