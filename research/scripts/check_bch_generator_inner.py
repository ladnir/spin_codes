#!/usr/bin/env python3
"""Diagnostics for the generator-side BCH block-recursive inner.

The candidate recursion is

    v_i = u_i + s_{i-1}
    (y_i, s_i) = E(v_i),

where E is a generator map for an extended primitive BCH [2b,b,d0] code.
Unlike the transposed/parity-check version, the behavior depends on the chosen
generator basis, not just on the code.  This script compares:

* the raw nullspace basis produced from the BCH parity checks;
* a systematic encoder on the left half:  (v, P v);
* a systematic encoder on the right half: (Q v, v).

It reports projection ranks and exact/sampled split-weight behavior for
low-weight driving words v.
"""

from __future__ import annotations

import argparse
import math
import random

from bch_candidate_params import bch_dimension
from check_bch_boundary_smallfield import find_primitive_poly
from check_bch_transposed_inner import extended_bch_check_rows, row_rank
from exact_bch_spectrum_small import nullspace_basis


def rref_rows(rows: list[int], n: int, preferred_pivots: list[int]) -> tuple[list[int], list[int]]:
    """Return RREF rows with pivots chosen from preferred columns when possible."""
    rows = [r for r in rows if r]
    rank = 0
    pivots: list[int] = []
    order = preferred_pivots + [c for c in range(n) if c not in set(preferred_pivots)]
    for col in order:
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
    return rows[:rank], pivots


def systematic_basis(code_basis: list[int], n: int, systematic_cols: list[int]) -> list[int]:
    """Change basis so the projection to systematic_cols is the identity.

    The input code_basis spans a b-dimensional code.  If projection onto the
    requested b columns is full rank, the returned basis vectors g_j satisfy
    g_j restricted to systematic_cols = e_j.
    """
    b = len(systematic_cols)
    rows: list[int] = []
    for i, word in enumerate(code_basis):
        proj = 0
        for j, col in enumerate(systematic_cols):
            if (word >> col) & 1:
                proj |= 1 << j
        rows.append(proj | (1 << (b + i)))

    rref, pivots = rref_rows(rows, 2 * b, list(range(b)))
    if pivots[:b] != list(range(b)):
        raise ValueError("requested systematic projection is not full rank")

    out = [0] * b
    for row in rref[:b]:
        pivot = (row & ((1 << b) - 1)).bit_length() - 1
        combo = row >> b
        word = 0
        x = combo
        while x:
            bit = x & -x
            idx = bit.bit_length() - 1
            word ^= code_basis[idx]
            x ^= bit
        out[pivot] = word
    return out


def encode(basis: list[int], msg: int) -> int:
    word = 0
    x = msg
    while x:
        bit = x & -x
        word ^= basis[bit.bit_length() - 1]
        x ^= bit
    return word


def iter_weight_masks(n: int, weight: int):
    if weight < 0 or weight > n:
        return
    if weight == 0:
        yield 0
        return
    combo = (1 << weight) - 1
    limit = 1 << n
    while combo < limit:
        yield combo
        c = combo & -combo
        r = combo + c
        combo = (((r ^ combo) >> 2) // c) | r


def split_weight(word: int, b: int) -> tuple[int, int]:
    y = word & ((1 << b) - 1)
    s = word >> b
    return y.bit_count(), s.bit_count()


def summarize_basis(name: str, basis: list[int], b: int, h_max: int, samples: int, seed: int) -> None:
    left_rank = row_rank([g & ((1 << b) - 1) for g in basis])
    right_rank = row_rank([g >> b for g in basis])
    row_y = [split_weight(g, b)[0] for g in basis]
    row_s = [split_weight(g, b)[1] for g in basis]
    print()
    print(f"encoder={name}")
    print(f"left_output_rank={left_rank}, right_state_rank={right_rank}")
    print(
        "row_weight_stats="
        f"y_min:{min(row_y)},y_avg:{sum(row_y)/len(row_y):.3f},y_max:{max(row_y)},"
        f"s_min:{min(row_s)},s_avg:{sum(row_s)/len(row_s):.3f},s_max:{max(row_s)}"
    )
    print("msg_w,count,min_y,avg_y,min_s,avg_s,min_total,avg_total,zero_state,zero_output")
    rng = random.Random(seed)
    for w in range(1, h_max + 1):
        total_patterns = math.comb(b, w)
        exact = total_patterns <= samples
        iterator = iter_weight_masks(b, w) if exact else (
            sum(1 << bit for bit in rng.sample(range(b), w)) for _ in range(samples)
        )
        count = 0
        min_y = b + 1
        min_s = b + 1
        min_total = 2 * b + 1
        sum_y = 0
        sum_s = 0
        sum_total = 0
        zero_s = 0
        zero_y = 0
        for msg in iterator:
            word = encode(basis, msg)
            y, s = split_weight(word, b)
            count += 1
            sum_y += y
            sum_s += s
            sum_total += y + s
            min_y = min(min_y, y)
            min_s = min(min_s, s)
            min_total = min(min_total, y + s)
            if s == 0:
                zero_s += 1
            if y == 0:
                zero_y += 1
        print(
            f"{w},{count},{min_y},{sum_y/count:.6f},{min_s},{sum_s/count:.6f},"
            f"{min_total},{sum_total/count:.6f},{zero_s},{zero_y}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--delta", type=int, default=7)
    parser.add_argument("--h-max", type=int, default=8)
    parser.add_argument("--samples", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    rows, ext_n, dim, rank = extended_bch_check_rows(args.m, args.delta)
    n, primitive_dim, redundancy = bch_dimension(args.m, args.delta)
    if ext_n != 2 * dim:
        raise SystemExit(f"extended code is not rate half: length={ext_n}, dim={dim}")
    b = dim
    basis = nullspace_basis(rows, ext_n)
    if len(basis) != b:
        raise SystemExit(f"unexpected nullspace dimension {len(basis)} != {b}")

    print("Generator-side extended-BCH inner diagnostic")
    print(
        f"m={args.m}, primitive_n={n}, extended_n={ext_n}, designed_delta={args.delta}, "
        f"dimension={dim}, parity_rank={rank}, block_bits={b}"
    )
    print(f"primitive_dimension={primitive_dim}, primitive_redundancy={redundancy}")

    left_cols = list(range(b))
    right_cols = list(range(b, 2 * b))
    variants = [("raw_nullspace", basis)]
    try:
        variants.append(("systematic_output_left", systematic_basis(basis, ext_n, left_cols)))
    except ValueError as exc:
        print(f"systematic_output_left unavailable: {exc}")
    try:
        variants.append(("systematic_state_right", systematic_basis(basis, ext_n, right_cols)))
    except ValueError as exc:
        print(f"systematic_state_right unavailable: {exc}")

    for name, variant in variants:
        summarize_basis(name, variant, b, args.h_max, args.samples, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
