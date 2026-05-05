#!/usr/bin/env python3
"""Check the divided-difference form of the BCH boundary interpolation condition.

For the normalized h=61 boundary problem, a monic degree-d polynomial g agrees
with y^r on a (d+1)-set S iff the top interpolation coefficient of y^r on S
is one.  For a monomial, that top coefficient is the complete homogeneous
symmetric polynomial h_{r-d}(S).

For q=512, d=30, r=255, the first consistency condition is therefore

    h_225(S) = 1

for a 31-set S.  This script verifies the identity against direct interpolation
and samples the condition.
"""

from __future__ import annotations

import argparse
import random

from check_bch_boundary_interp_sampling import (
    gf_inv,
    gf_mul,
    solve_linear_square,
)


def precompute_tables(q: int, poly: int, m: int, max_power: int) -> tuple[list[list[int]], list[int], list[list[int]]]:
    mul_table = [[0] * q for _ in range(q)]
    for a in range(q):
        for b in range(q):
            mul_table[a][b] = gf_mul(a, b, poly, m)
    inv_table = [0] * q
    for a in range(1, q):
        inv_table[a] = gf_inv(a, poly, m)
    powers = [[0] * (max_power + 1) for _ in range(q)]
    for y in range(q):
        powers[y][0] = 1
        for e in range(1, max_power + 1):
            powers[y][e] = mul_table[powers[y][e - 1]][y]
    return mul_table, inv_table, powers


def complete_homogeneous_fast(points: list[int], degree: int, powers: list[list[int]], mul_table: list[list[int]]) -> int:
    h = [0] * (degree + 1)
    h[0] = 1
    for y in points:
        old = h[:]
        for t in range(1, degree + 1):
            acc = old[t]
            for j in range(1, t + 1):
                acc ^= mul_table[old[t - j]][powers[y][j]]
            h[t] = acc
    return h[degree]


def interpolation_top_coeff(
    points: list[int],
    d: int,
    r: int,
    inv_table: list[int],
    mul_table: list[list[int]],
    powers: list[list[int]],
) -> int | None:
    # Solve for degree <= d interpolation of y^r on d+1 points and return top coeff.
    a = []
    b = []
    for y in points:
        a.append(powers[y][: d + 1])
        b.append(powers[y][r])
    coeffs = solve_linear_square(a, b, inv_table, mul_table)
    return None if coeffs is None else coeffs[d]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=9)
    parser.add_argument("--poly", type=lambda x: int(x, 0), default=0x211)
    parser.add_argument("--degree", type=int, default=30)
    parser.add_argument("--target-power", type=int, default=255)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--check-trials", type=int, default=20)
    args = parser.parse_args()

    q = 1 << args.m
    field = list(range(1, q))
    sym_degree = args.target_power - args.degree
    print("precomputing field tables", flush=True)
    mul_table, inv_table, powers = precompute_tables(q, args.poly, args.m, args.target_power)
    rng = random.Random(args.seed)

    mismatches = 0
    for _ in range(args.check_trials):
        pts = rng.sample(field, args.degree + 1)
        top = interpolation_top_coeff(pts, args.degree, args.target_power, inv_table, mul_table, powers)
        sym = complete_homogeneous_fast(pts, sym_degree, powers, mul_table)
        if top != sym:
            mismatches += 1

    successes = 0
    for _ in range(args.trials):
        pts = rng.sample(field, args.degree + 1)
        sym = complete_homogeneous_fast(pts, sym_degree, powers, mul_table)
        if sym == 1:
            successes += 1

    print("BCH boundary divided-difference sampler")
    print(f"m={args.m}, q={q}, degree={args.degree}, target_power={args.target_power}")
    print(f"condition=h_{sym_degree}(S)=1 for |S|={args.degree + 1}")
    print(f"identity_check_trials={args.check_trials}, mismatches={mismatches}")
    print(f"sample_trials={args.trials}, successes={successes}")
    print(f"empirical_probability={successes / args.trials:.12g}")
    print(f"random_probability={1 / q:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
