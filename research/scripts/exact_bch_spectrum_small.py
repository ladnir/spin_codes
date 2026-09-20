#!/usr/bin/env python3
"""Exact spectra for small primitive narrow-sense binary BCH codes.

This is the small-instance ground-truth tool for the BCH spectra program.  It
builds the binary parity-check matrix from the BCH power-sum constraints,
computes a binary nullspace basis, and enumerates all codewords by Gray code.

The default instance is the primitive [31,16,7] BCH code.  Larger dimensions
grow exponentially; use --max-dim deliberately.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from bch_candidate_params import bch_dimension, cyclotomic_coset
from check_bch_boundary_smallfield import find_primitive_poly, gf_mul
from dense_largek_eval import log2_binom


def coset_reps(n: int, delta: int) -> list[int]:
    seen: set[int] = set()
    reps: list[int] = []
    for e in range(1, delta):
        c = cyclotomic_coset(e, n)
        if not any(x in seen for x in c):
            reps.append(e)
            seen.update(c)
    return reps


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


def parity_rows(m: int, delta: int, poly: int) -> list[int]:
    n = (1 << m) - 1
    alpha = 2
    rows: list[int] = []
    reps = coset_reps(n, delta)
    for e in reps:
        bit_rows = [0] * m
        for pos in range(n):
            val = gf_pow(alpha, (pos * e) % n, poly, m)
            for bit in range(m):
                if (val >> bit) & 1:
                    bit_rows[bit] |= 1 << pos
        rows.extend(bit_rows)
    return rows


def nullspace_basis(rows: list[int], n: int) -> list[int]:
    rows = [r for r in rows if r]
    rank = 0
    pivots: list[int] = []
    for col in range(n):
        pivot = None
        mask = 1 << col
        for i in range(rank, len(rows)):
            if rows[i] & mask:
                pivot = i
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] & mask):
                rows[i] ^= rows[rank]
        pivots.append(col)
        rank += 1
        if rank == len(rows):
            break

    pivot_set = set(pivots)
    free_cols = [c for c in range(n) if c not in pivot_set]
    basis: list[int] = []
    pivot_rows = rows[:rank]
    for free in free_cols:
        vec = 1 << free
        free_mask = 1 << free
        for pivot_col, row in zip(pivots, pivot_rows):
            if row & free_mask:
                vec |= 1 << pivot_col
        basis.append(vec)
    return basis


def enumerate_spectrum(basis: list[int], n: int) -> list[int]:
    k = len(basis)
    spectrum = [0] * (n + 1)
    current = 0
    prev_gray = 0
    spectrum[0] = 1
    for t in range(1, 1 << k):
        gray = t ^ (t >> 1)
        diff = gray ^ prev_gray
        idx = diff.bit_length() - 1
        current ^= basis[idx]
        spectrum[current.bit_count()] += 1
        prev_gray = gray
    return spectrum


def write_spectrum(path: Path, spectrum: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["weight", "count"])
        for weight, count in enumerate(spectrum):
            if count:
                writer.writerow([weight, count])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--delta", type=int, default=7)
    parser.add_argument("--max-dim", type=int, default=24)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    n, dim, redundancy = bch_dimension(args.m, args.delta)
    if dim > args.max_dim:
        raise SystemExit(f"dimension {dim} exceeds --max-dim {args.max_dim}")
    poly = find_primitive_poly(args.m)
    rows = parity_rows(args.m, args.delta, poly)
    basis = nullspace_basis(rows, n)
    spectrum = enumerate_spectrum(basis, n)
    total = sum(spectrum)
    dmin = next(i for i, count in enumerate(spectrum) if i > 0 and count)

    print("Exact small BCH spectrum")
    print(f"m={args.m}, n={n}, delta={args.delta}, primitive_poly=0x{poly:x}")
    print(f"dimension_from_cosets={dim}, redundancy={redundancy}")
    print(f"dimension_from_nullspace={len(basis)}")
    print(f"total_codewords={total}")
    print(f"minimum_weight={dmin}")
    print("weight,count,randomlike_log2_count,loss_vs_randomlike_bits")
    for weight, count in enumerate(spectrum):
        if not count:
            continue
        randomlike = log2_binom(n, weight) + dim - n
        loss = math.log2(count) - randomlike if count > 0 else float("-inf")
        print(f"{weight},{count},{randomlike:.6f},{loss:.6f}")
    if args.out is not None:
        write_spectrum(args.out, spectrum)
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
