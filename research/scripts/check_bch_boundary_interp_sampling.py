#!/usr/bin/env python3
"""Sample overdetermined interpolation consistency for the BCH boundary target.

For q=512, h=61, and degree d=30, the normalized boundary problem asks for
many points y where a monic degree-d polynomial g agrees with y^255.

For a random s-subset S of F_q^*, consistency means that there is a monic
degree-d g agreeing with y^255 on every y in S.  If the overdetermined
conditions behave randomly, the probability should be about q^{-(s-d)}.

This sampler is only a diagnostic; it tests whether small excesses s-d show
obvious bias before we try to prove a finite rank/equidistribution lemma.
"""

from __future__ import annotations

import argparse
import random


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


def gf_inv(a: int, poly: int, m: int) -> int:
    if a == 0:
        raise ZeroDivisionError
    return gf_pow(a, (1 << m) - 2, poly, m)


def solve_linear_square(a: list[list[int]], b: list[int], inv_table: list[int], mul_table: list[list[int]]) -> list[int] | None:
    n = len(b)
    mat = [row[:] + [rhs] for row, rhs in zip(a, b)]
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if mat[r][col]:
                pivot = r
                break
        if pivot is None:
            return None
        if pivot != col:
            mat[col], mat[pivot] = mat[pivot], mat[col]
        inv = inv_table[mat[col][col]]
        for c in range(col, n + 1):
            mat[col][c] = mul_table[mat[col][c]][inv]
        for r in range(n):
            if r == col or mat[r][col] == 0:
                continue
            factor = mat[r][col]
            for c in range(col, n + 1):
                mat[r][c] ^= mul_table[factor][mat[col][c]]
    return [mat[i][n] for i in range(n)]


def eval_linear(coeffs: list[int], powers: list[int], mul_table: list[list[int]], leading_value: int) -> int:
    acc = leading_value
    for c, yp in zip(coeffs, powers):
        if c:
            acc ^= mul_table[c][yp]
    return acc


def consistent(
    points: list[int],
    degree: int,
    target_power: int,
    inv_table: list[int],
    mul_table: list[list[int]],
    powers: list[list[int]],
) -> bool:
    base = points[:degree]
    a: list[list[int]] = []
    b: list[int] = []
    for y in base:
        row = powers[y][:degree]
        a.append(row)
        b.append(powers[y][target_power] ^ powers[y][degree])
    coeffs = solve_linear_square(a, b, inv_table, mul_table)
    if coeffs is None:
        return False
    for y in points[degree:]:
        if eval_linear(coeffs, powers[y][:degree], mul_table, powers[y][degree]) != powers[y][target_power]:
            return False
    return True


def parse_int_list(text: str) -> list[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=9)
    parser.add_argument("--poly", type=lambda x: int(x, 0), default=0x211)
    parser.add_argument("--degree", type=int, default=30)
    parser.add_argument("--target-power", type=int, default=255)
    parser.add_argument("--s-values", default="31,32")
    parser.add_argument("--trials", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    q = 1 << args.m
    field = list(range(1, q))
    print("precomputing field tables", flush=True)
    mul_table = [[0] * q for _ in range(q)]
    for a in range(q):
        for b in range(q):
            mul_table[a][b] = gf_mul(a, b, args.poly, args.m)
    inv_table = [0] * q
    for a in range(1, q):
        inv_table[a] = gf_inv(a, args.poly, args.m)
    max_power = max(args.target_power, args.degree)
    powers = [[0] * (max_power + 1) for _ in range(q)]
    powers[0][0] = 1
    for y in field:
        powers[y][0] = 1
        for e in range(1, max_power + 1):
            powers[y][e] = mul_table[powers[y][e - 1]][y]
    rng = random.Random(args.seed)
    print("BCH boundary interpolation consistency sampler")
    print(f"m={args.m}, q={q}, poly=0x{args.poly:x}, degree={args.degree}, target_power={args.target_power}")
    print(f"trials={args.trials}, seed={args.seed}")
    print("s,successes,empirical_probability,random_probability")
    for s in parse_int_list(args.s_values):
        successes = 0
        for _ in range(args.trials):
            pts = rng.sample(field, s)
            if consistent(pts, args.degree, args.target_power, inv_table, mul_table, powers):
                successes += 1
        empirical = successes / args.trials
        random_prob = q ** -(s - args.degree)
        print(f"{s},{successes},{empirical:.12g},{random_prob:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
