#!/usr/bin/env python3
"""Audit BCH boundary interpolation-rank constraints.

This is a Route-B tool for the BCH boundary program.  For q=2^m and odd
boundary h=2d+1, normalized boundary supports are exactly the h-sets
S in F_q^* on which a monic degree-d polynomial g agrees with

    y^(q/2 - 1).

Equivalently, the overdetermined interpolation system with d unknown lower
coefficients is consistent on S.  Exact h-set enumeration is only intended for
small fields, but it validates the rank formulation against the m=5 exact BCH
boundary counts.  The sampling mode probes the same rank conditions for larger
fields and can search for split boundary polynomials without enumerating all
q^d locator families.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from itertools import combinations
from pathlib import Path

from check_bch_boundary_smallfield import find_primitive_poly, gf_mul


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


def precompute_tables(q: int, poly: int, m: int, max_power: int) -> tuple[list[list[int]], list[int], list[list[int]]]:
    mul_table = [[0] * q for _ in range(q)]
    for a in range(q):
        row = mul_table[a]
        for b in range(q):
            row[b] = gf_mul(a, b, poly, m)

    inv_table = [0] * q
    for a in range(1, q):
        inv_table[a] = gf_pow(a, q - 2, poly, m)

    powers = [[0] * (max_power + 1) for _ in range(q)]
    for y in range(q):
        powers[y][0] = 1
        prow = powers[y]
        for e in range(1, max_power + 1):
            prow[e] = mul_table[prow[e - 1]][y]
    return mul_table, inv_table, powers


def solve_square(
    a: list[list[int]],
    b: list[int],
    inv_table: list[int],
    mul_table: list[list[int]],
) -> list[int] | None:
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
        pivot_row = mat[col]
        for c in range(col, n + 1):
            pivot_row[c] = mul_table[pivot_row[c]][inv]
        for r in range(n):
            if r == col or mat[r][col] == 0:
                continue
            factor = mat[r][col]
            row = mat[r]
            for c in range(col, n + 1):
                row[c] ^= mul_table[factor][pivot_row[c]]
    return [mat[i][n] for i in range(n)]


def interpolate_monic(
    points: tuple[int, ...],
    degree: int,
    target_power: int,
    inv_table: list[int],
    mul_table: list[list[int]],
    powers: list[list[int]],
) -> list[int] | None:
    base = points[:degree]
    a: list[list[int]] = []
    b: list[int] = []
    for y in base:
        a.append(powers[y][:degree])
        b.append(powers[y][target_power] ^ powers[y][degree])
    return solve_square(a, b, inv_table, mul_table)


def eval_monic(coeffs: list[int], y: int, degree: int, mul_table: list[list[int]], powers: list[list[int]]) -> int:
    out = powers[y][degree]
    for i, c in enumerate(coeffs):
        if c:
            out ^= mul_table[c][powers[y][i]]
    return out


def consistent_set(
    points: tuple[int, ...],
    degree: int,
    target_power: int,
    inv_table: list[int],
    mul_table: list[list[int]],
    powers: list[list[int]],
) -> tuple[bool, tuple[int, ...] | None]:
    coeffs = interpolate_monic(points, degree, target_power, inv_table, mul_table, powers)
    if coeffs is None:
        return False, None
    for y in points[degree:]:
        if eval_monic(coeffs, y, degree, mul_table, powers) != powers[y][target_power]:
            return False, None
    return True, tuple(coeffs)


def count_roots(
    coeffs: tuple[int, ...],
    field: list[int],
    degree: int,
    target_power: int,
    mul_table: list[list[int]],
    powers: list[list[int]],
) -> int:
    coeff_list = list(coeffs)
    roots = 0
    for y in field:
        if eval_monic(coeff_list, y, degree, mul_table, powers) == powers[y][target_power]:
            roots += 1
    return roots


def read_spectrum_count(path: Path, h: int) -> int | None:
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if int(row["weight"]) == h:
                return int(row["count"])
    return None


def exact_hsets(args: argparse.Namespace) -> int:
    q = 1 << args.m
    n = q - 1
    h = args.h
    degree = args.degree
    target_power = args.target_power
    total = math.comb(n, h)
    if total > args.max_combinations:
        raise SystemExit(
            f"exact h-set enumeration would visit {total} sets; raise --max-combinations deliberately"
        )

    poly = find_primitive_poly(args.m) if args.poly is None else args.poly
    field = list(range(1, q))
    mul_table, inv_table, powers = precompute_tables(q, poly, args.m, max(target_power, degree))

    consistent = 0
    root_hist: dict[int, int] = {}
    split_polys: set[tuple[int, ...]] = set()
    for pts in combinations(field, h):
        ok, coeffs = consistent_set(pts, degree, target_power, inv_table, mul_table, powers)
        if not ok or coeffs is None:
            continue
        consistent += 1
        roots = count_roots(coeffs, field, degree, target_power, mul_table, powers)
        root_hist[roots] = root_hist.get(roots, 0) + 1
        if roots == h:
            split_polys.add(coeffs)

    orbit_size = n // math.gcd(h, n)
    print("BCH boundary interpolation-rank exact h-set audit")
    print(f"m={args.m}, q={q}, n={n}, h={h}, degree={degree}, target_power={target_power}, poly=0x{poly:x}")
    print(f"hsets_tested={total}")
    print(f"consistent_hsets={consistent}")
    print(f"split_polynomials={len(split_polys)}")
    print(f"implied_boundary_codewords={len(split_polys) * orbit_size}")
    if args.compare_spectrum is not None:
        exact = read_spectrum_count(args.compare_spectrum, h)
        print(f"exact_spectrum_boundary_count={exact}")
        print(f"matches_exact_spectrum={exact == len(split_polys) * orbit_size}")
    print("consistent_root_hist=" + ",".join(f"{k}:{root_hist[k]}" for k in sorted(root_hist)))
    return 0


def sample_sets(args: argparse.Namespace) -> int:
    q = 1 << args.m
    n = q - 1
    s = args.s
    degree = args.degree
    target_power = args.target_power
    poly = find_primitive_poly(args.m) if args.poly is None else args.poly
    field = list(range(1, q))
    mul_table, inv_table, powers = precompute_tables(q, poly, args.m, max(target_power, degree))
    rng = random.Random(args.seed)

    consistent = 0
    split_polys: set[tuple[int, ...]] = set()
    root_hist: dict[int, int] = {}
    for _ in range(args.trials):
        pts = tuple(rng.sample(field, s))
        ok, coeffs = consistent_set(pts, degree, target_power, inv_table, mul_table, powers)
        if not ok or coeffs is None:
            continue
        consistent += 1
        if args.scan_roots:
            roots = count_roots(coeffs, field, degree, target_power, mul_table, powers)
            root_hist[roots] = root_hist.get(roots, 0) + 1
            if roots == args.h:
                split_polys.add(coeffs)

    random_prob = q ** -max(0, s - degree)
    print("BCH boundary interpolation-rank sampler")
    print(f"m={args.m}, q={q}, n={n}, h={args.h}, s={s}, degree={degree}, target_power={target_power}, poly=0x{poly:x}")
    print(f"trials={args.trials}, seed={args.seed}")
    print(f"consistent_sets={consistent}")
    print(f"empirical_probability={consistent / args.trials:.12g}")
    print(f"random_rank_probability={random_prob:.12g}")
    if args.scan_roots:
        orbit_size = n // math.gcd(args.h, n)
        print(f"unique_split_polynomials_found={len(split_polys)}")
        print(f"found_boundary_codewords={len(split_polys) * orbit_size}")
        if args.compare_spectrum is not None:
            exact = read_spectrum_count(args.compare_spectrum, args.h)
            print(f"exact_spectrum_boundary_count={exact}")
        print("consistent_root_hist=" + ",".join(f"{k}:{root_hist[k]}" for k in sorted(root_hist)))
    return 0


def recover_polys(args: argparse.Namespace) -> int:
    q = 1 << args.m
    n = q - 1
    degree = args.degree
    target_power = args.target_power
    poly = find_primitive_poly(args.m) if args.poly is None else args.poly
    field = list(range(1, q))
    mul_table, inv_table, powers = precompute_tables(q, poly, args.m, max(target_power, degree))
    rng = random.Random(args.seed)

    orbit_size = n // math.gcd(args.h, n)
    target_polys = args.target_polys
    exact = None
    if args.compare_spectrum is not None:
        exact = read_spectrum_count(args.compare_spectrum, args.h)
        if exact is not None:
            if exact % orbit_size != 0:
                raise SystemExit(f"exact boundary count {exact} is not divisible by orbit size {orbit_size}")
            target_polys = exact // orbit_size

    split_polys: set[tuple[int, ...]] = set()
    split_hits = 0
    root_hist: dict[int, int] = {}
    for trial in range(1, args.trials + 1):
        pts = tuple(rng.sample(field, degree))
        coeffs = interpolate_monic(pts, degree, target_power, inv_table, mul_table, powers)
        if coeffs is None:
            continue
        coeff_key = tuple(coeffs)
        roots = count_roots(coeff_key, field, degree, target_power, mul_table, powers)
        root_hist[roots] = root_hist.get(roots, 0) + 1
        if roots == args.h:
            split_hits += 1
            split_polys.add(coeff_key)
            if target_polys is not None and len(split_polys) >= target_polys:
                break

    print("BCH boundary split-polynomial recovery sampler")
    print(f"m={args.m}, q={q}, n={n}, h={args.h}, degree={degree}, target_power={target_power}, poly=0x{poly:x}")
    print(f"trials_limit={args.trials}, seed={args.seed}")
    print(f"unique_split_polynomials_found={len(split_polys)}")
    print(f"split_hits_in_samples={split_hits}")
    print(f"found_boundary_codewords={len(split_polys) * orbit_size}")
    if exact is not None:
        print(f"exact_spectrum_boundary_count={exact}")
        print(f"matches_exact_spectrum={exact == len(split_polys) * orbit_size}")
    if target_polys is not None:
        print(f"target_split_polynomials={target_polys}")
        print(f"target_reached={len(split_polys) >= target_polys}")
    print("root_hist=" + ",".join(f"{k}:{root_hist[k]}" for k in sorted(root_hist)))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--h", type=int, default=7)
    parser.add_argument("--degree", type=int, default=None)
    parser.add_argument("--target-power", type=int, default=None)
    parser.add_argument("--poly", type=lambda x: int(x, 0), default=None)
    parser.add_argument("--compare-spectrum", type=Path, default=None)
    sub = parser.add_subparsers(dest="mode", required=True)

    exact = sub.add_parser("exact-hsets", help="exactly enumerate h-sets; small fields only")
    exact.add_argument("--max-combinations", type=int, default=5_000_000)

    sample = sub.add_parser("sample-sets", help="sample s-sets for rank consistency")
    sample.add_argument("--s", type=int, default=None)
    sample.add_argument("--trials", type=int, default=5000)
    sample.add_argument("--seed", type=int, default=1)
    sample.add_argument("--scan-roots", action="store_true")

    recover = sub.add_parser("recover-polys", help="sample d-sets until split polynomials are recovered")
    recover.add_argument("--trials", type=int, default=200000)
    recover.add_argument("--seed", type=int, default=1)
    recover.add_argument("--target-polys", type=int, default=None)

    args = parser.parse_args()
    if args.h % 2 != 1:
        raise ValueError("--h must be odd")
    args.degree = (args.h - 1) // 2 if args.degree is None else args.degree
    args.target_power = ((1 << args.m) // 2 - 1) if args.target_power is None else args.target_power
    if args.mode == "sample-sets" and args.s is None:
        args.s = args.degree + 1

    if args.mode == "exact-hsets":
        return exact_hsets(args)
    if args.mode == "recover-polys":
        return recover_polys(args)
    return sample_sets(args)


if __name__ == "__main__":
    raise SystemExit(main())
