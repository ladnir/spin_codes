#!/usr/bin/env python3
"""Exact small-field counts for the boundary splitting family.

This is a sanity-check model for the q=512 boundary target.  For small
q=2^m and odd h=2d+1, enumerate all monic degree-d polynomials g and count how
often

    F(Y) = Y g(Y)^2 + 1

has exactly h nonzero roots.  The q=512 target is too large for enumeration,
but small fields help detect obvious structural bias in the lacunary family.
"""

from __future__ import annotations

import argparse
import math

from dense_largek_eval import log2_binom


def gf_mul(a: int, b: int, poly: int, m: int) -> int:
    out = 0
    aa = a
    bb = b
    while bb:
        if bb & 1:
            out ^= aa
        bb >>= 1
        aa <<= 1
        if aa & (1 << m):
            aa ^= poly
    return out & ((1 << m) - 1)


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


def is_irreducible(poly: int, m: int) -> bool:
    # Brute force divisibility by monic polynomials of degree 1..m//2.
    for deg in range(1, m // 2 + 1):
        for div_low in range(1 << deg):
            div = (1 << deg) | div_low
            rem = poly
            while rem.bit_length() - 1 >= deg:
                rem ^= div << (rem.bit_length() - 1 - deg)
            if rem == 0:
                return False
    return True


def find_primitive_poly(m: int) -> int:
    n = (1 << m) - 1
    for low in range(1 << m):
        poly = (1 << m) | low
        if poly & 1 == 0 or not is_irreducible(poly, m):
            continue
        x = 2
        if gf_pow(x, n, poly, m) != 1:
            continue
        primitive = True
        for p in prime_factors(n):
            if gf_pow(x, n // p, poly, m) == 1:
                primitive = False
                break
        if primitive:
            return poly
    raise ValueError(f"no primitive polynomial found for m={m}")


def prime_factors(n: int) -> list[int]:
    out = []
    d = 2
    x = n
    while d * d <= x:
        if x % d == 0:
            out.append(d)
            while x % d == 0:
                x //= d
        d += 1
    if x > 1:
        out.append(x)
    return out


def eval_g(coeffs: list[int], y: int, poly: int, m: int) -> int:
    acc = 1  # monic leading coefficient
    for c in reversed(coeffs):
        acc = gf_mul(acc, y, poly, m) ^ c
    return acc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--h", type=int, default=5)
    args = parser.parse_args()

    q = 1 << args.m
    n = q - 1
    if args.h % 2 != 1:
        raise ValueError("--h must be odd")
    d = (args.h - 1) // 2
    poly = find_primitive_poly(args.m)
    field = list(range(1, q))
    root_hist: dict[int, int] = {}
    split_count = 0
    total = q**d
    for idx in range(total):
        x = idx
        coeffs = []
        for _ in range(d):
            coeffs.append(x & (q - 1))
            x >>= args.m
        roots = 0
        for y in field:
            gy = eval_g(coeffs, y, poly, args.m)
            val = gf_mul(y, gf_mul(gy, gy, poly, args.m), poly, args.m) ^ 1
            if val == 0:
                roots += 1
        root_hist[roots] = root_hist.get(roots, 0) + 1
        if roots == args.h:
            split_count += 1

    random_split = log2_binom(n, args.h) - args.h * math.log2(q)
    expected = math.log2(total) + random_split
    actual = math.log2(split_count) if split_count else float("-inf")
    print("Small-field BCH boundary splitting")
    print(f"m={args.m}, q={q}, n={n}, h={args.h}, d={d}, primitive_poly=0x{poly:x}")
    print(f"total_g={total}")
    print(f"split_count={split_count}")
    print(f"actual_split_count_log2={actual:.6f}")
    print(f"random_expected_split_count_log2={expected:.6f}")
    print("root_hist=" + ",".join(f"{k}:{root_hist[k]}" for k in sorted(root_hist)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
