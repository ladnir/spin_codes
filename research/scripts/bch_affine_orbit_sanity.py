#!/usr/bin/env python3
"""Sanity checks for affine orbits of extended BCH supports.

The length-128/256 BCH spectrum computations exploit the affine automorphism
group of the extended primitive BCH code.  This small-instance tool builds the
same coordinate model: primitive coordinates are F_q^*, the extension parity
coordinate is 0, and affine maps act by

    x -> a x + b,  a in F_q^*, b in F_q.

For exact-enumerable BCH instances, it enumerates codewords of selected
primitive weights, extends them by the overall parity bit, and partitions their
supports into affine orbits.  This is the first low-risk validation layer for
the affine-coset spectrum ladder; it is not intended to enumerate large codes.
"""

from __future__ import annotations

import argparse
from collections import Counter

from check_bch_boundary_smallfield import find_primitive_poly, gf_mul
from exact_bch_spectrum_small import enumerate_spectrum, nullspace_basis, parity_rows


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


def primitive_position_values(m: int, poly: int) -> list[int]:
    n = (1 << m) - 1
    alpha = 2
    values = [1] * n
    x = 1
    for pos in range(n):
        values[pos] = x
        x = gf_mul(x, alpha, poly, m)
    return values


def extended_support_mask(word: int, values: list[int], primitive_weight: int) -> int:
    mask = 0
    x = word
    while x:
        bit = x & -x
        pos = bit.bit_length() - 1
        mask |= 1 << values[pos]
        x ^= bit
    if primitive_weight & 1:
        mask |= 1
    return mask


def mask_points(mask: int) -> list[int]:
    points: list[int] = []
    x = mask
    while x:
        bit = x & -x
        points.append(bit.bit_length() - 1)
        x ^= bit
    return points


def affine_maps(q: int, poly: int, m: int) -> list[tuple[int, int]]:
    return [(a, b) for a in range(1, q) for b in range(q)]


def affine_image_mask(points: list[int], a: int, b: int, poly: int, m: int) -> int:
    out = 0
    for x in points:
        out |= 1 << (gf_mul(a, x, poly, m) ^ b)
    return out


def affine_canonical(mask: int, maps: list[tuple[int, int]], poly: int, m: int) -> tuple[int, int]:
    points = mask_points(mask)
    best: int | None = None
    stabilizer = 0
    for a, b in maps:
        image = affine_image_mask(points, a, b, poly, m)
        if image == mask:
            stabilizer += 1
        if best is None or image < best:
            best = image
    if best is None:
        raise AssertionError("empty affine map list")
    return best, stabilizer


def enumerate_codeword_weights(basis: list[int], wanted: set[int]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    current = 0
    prev_gray = 0
    if 0 in wanted:
        out.append((0, 0))
    for t in range(1, 1 << len(basis)):
        gray = t ^ (t >> 1)
        diff = gray ^ prev_gray
        idx = diff.bit_length() - 1
        current ^= basis[idx]
        weight = current.bit_count()
        if weight in wanted:
            out.append((current, weight))
        prev_gray = gray
    return out


def parse_weights(text: str | None, dmin: int) -> set[int]:
    if text is None:
        return {dmin}
    return {int(x.strip()) for x in text.split(",") if x.strip()}


def parse_extended_weights(text: str | None) -> set[int]:
    if text is None:
        return set()
    return {int(x.strip()) for x in text.split(",") if x.strip()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--delta", type=int, default=7)
    parser.add_argument("--weights", default=None, help="Primitive weights to inspect, comma-separated; default dmin")
    parser.add_argument(
        "--extended-weights",
        default=None,
        help="Extended weights to inspect; includes primitive weights W and W-1",
    )
    parser.add_argument("--max-dim", type=int, default=21)
    args = parser.parse_args()

    q = 1 << args.m
    n = q - 1
    poly = find_primitive_poly(args.m)
    rows = parity_rows(args.m, args.delta, poly)
    basis = nullspace_basis(rows, n)
    if len(basis) > args.max_dim:
        raise SystemExit(f"dimension {len(basis)} exceeds --max-dim {args.max_dim}")

    spectrum = enumerate_spectrum(basis, n)
    dmin = next(w for w, c in enumerate(spectrum) if w > 0 and c)
    wanted = parse_weights(args.weights, dmin)
    for extended_weight in parse_extended_weights(args.extended_weights):
        wanted.update([extended_weight - 1, extended_weight])
    wanted = {w for w in wanted if 0 <= w <= n}
    values = primitive_position_values(args.m, poly)
    maps = affine_maps(q, poly, args.m)

    orbit_counts: Counter[int] = Counter()
    stabilizers: Counter[int] = Counter()
    by_extended_weight: Counter[int] = Counter()
    for word, weight in enumerate_codeword_weights(basis, wanted):
        mask = extended_support_mask(word, values, weight)
        canonical, stabilizer = affine_canonical(mask, maps, poly, args.m)
        orbit_counts[canonical] += 1
        stabilizers[stabilizer] += 1
        by_extended_weight[mask.bit_count()] += 1

    print("BCH affine support-orbit sanity")
    print(f"m={args.m}, q={q}, n={n}, delta={args.delta}, primitive_poly=0x{poly:x}")
    print(f"dimension={len(basis)}, primitive_dmin={dmin}")
    print(f"weights={','.join(str(w) for w in sorted(wanted))}")
    print(f"codewords_examined={sum(orbit_counts.values())}")
    print(f"affine_group_size={len(maps)}")
    print(f"affine_orbits={len(orbit_counts)}")
    print("extended_weight_hist=" + ",".join(f"{k}:{by_extended_weight[k]}" for k in sorted(by_extended_weight)))
    print("orbit_size_hist=" + ",".join(f"{v}:{c}" for v, c in sorted(Counter(orbit_counts.values()).items())))
    print("stabilizer_hist=" + ",".join(f"{k}:{stabilizers[k]}" for k in sorted(stabilizers)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
