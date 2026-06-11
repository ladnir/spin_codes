#!/usr/bin/env python3
"""Exact input/state cancellation hazard for the full-split BCH inner.

After a nonzero full-split branch creates state weight q, exact cancellation in
the next block requires the remaining input support to equal that fixed q-set.
With H remaining ones over T remaining blocks, the exact-hit probability is

    C(b(T-1), H-q) / C(bT, H).

This helper averages that quantity over the BCH split state-weight law and
reports whether it is comparable to the q'=0 self-turnoff atom.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from dense_largek_eval import log2_binom, log2add
from import_wd_spectrum import parse_wd
from probe_block_recursive_inner import parse_int_list


def load_spectrum(path: Path) -> list[tuple[int, int]]:
    if path.suffix.lower() == ".wd":
        return parse_wd(path.read_text())
    raise ValueError("only .wd spectra are supported")


def split_q_law(spectrum: list[tuple[int, int]], b: int) -> list[float]:
    nonzero_total = sum(c for w, c in spectrum if w > 0)
    pq = [0.0 for _ in range(b + 1)]
    for w, count in spectrum:
        if w <= 0 or count <= 0:
            continue
        p_w = count / nonzero_total
        den = math.comb(2 * b, w)
        for q in range(max(0, w - b), min(b, w) + 1):
            pq[q] += p_w * math.comb(b, q) * math.comb(b, w - q) / den
    return pq


def exact_hit_log2(*, pq: list[float], b: int, remaining_ones: int, remaining_blocks: int) -> float:
    total = float("-inf")
    n = b * remaining_blocks
    rest = b * (remaining_blocks - 1)
    for q, p in enumerate(pq):
        if q == 0 or p <= 0.0 or q > remaining_ones:
            continue
        if remaining_ones - q > rest:
            continue
        hit = log2_binom(rest, remaining_ones - q) - log2_binom(n, remaining_ones)
        total = log2add(total, math.log2(p) + hit)
    return total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--remaining-ones", default="31,79,159,319,499")
    parser.add_argument("--remaining-blocks", default="5949,9949,13949")
    args = parser.parse_args()

    b = args.block_bits
    pq = split_q_law(load_spectrum(args.spectrum), b)
    p0_log2 = math.log2(pq[0])
    best: tuple[int, int, float] | None = None

    print("Full-split exact-cancellation diagnostic")
    print(f"b={b}, spectrum={args.spectrum}")
    print(f"self_turnoff_log2,{p0_log2:.6f}")
    print("remaining_ones,remaining_blocks,exact_hit_log2,union_over_blocks_log2,effective_atom_log2")
    for H in parse_int_list(args.remaining_ones):
        for T in parse_int_list(args.remaining_blocks):
            hit = exact_hit_log2(pq=pq, b=b, remaining_ones=H, remaining_blocks=T)
            union = hit + math.log2(T) if hit != float("-inf") else float("-inf")
            effective = log2add(p0_log2, hit)
            print(f"{H},{T},{hit:.6f},{union:.6f},{effective:.6f}")
            if best is None or hit > best[2]:
                best = (H, T, hit)
    if best is not None:
        H, T, hit = best
        print("summary")
        print(f"max_hit_remaining_ones,{H}")
        print(f"max_hit_remaining_blocks,{T}")
        print(f"max_exact_hit_log2,{hit:.6f}")
        print(f"max_effective_atom_log2,{log2add(p0_log2, hit):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
