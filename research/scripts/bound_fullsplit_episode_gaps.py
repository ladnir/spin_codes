#!/usr/bin/env python3
"""Episode-gap upper bound for the full-split BCH inner.

This diagnostic models the missing short-episode branch after the survival
split.  Fix the first active block and let T blocks remain after it.  If a
self-turnoff q'=0 occurs after an occupied/live block, the following zero gap is
off and can be skipped.  For a random set of remaining occupied blocks, this
script union-bounds over selected skipped gaps:

    p0^e * E[# selected e-gap subsets with total skipped length s]
          * Chernoff_survival(T - s).

The output is intentionally pessimistic: q'=0 branches pay only p0 and their
emitted weight is ignored.  Exact input/state cancellations are not included
yet; this isolates the self-turnoff/restart mechanism.
"""

from __future__ import annotations

import argparse
import math
from functools import lru_cache
from pathlib import Path

from dense_largek_eval import log2_binom, log2add
from import_wd_spectrum import parse_wd


def load_spectrum(path: Path) -> list[tuple[int, int]]:
    if path.suffix.lower() == ".wd":
        return parse_wd(path.read_text())
    raise ValueError("only .wd spectra are supported")


def split_entries(spectrum: list[tuple[int, int]], b: int) -> list[tuple[int, int, float]]:
    nonzero_total = sum(c for w, c in spectrum if w > 0)
    out: list[tuple[int, int, float]] = []
    for w, count in spectrum:
        if w <= 0 or count <= 0:
            continue
        p_w = count / nonzero_total
        den = math.comb(2 * b, w)
        for q in range(max(0, w - b), min(b, w) + 1):
            j = w - q
            p = p_w * math.comb(b, q) * math.comb(b, j) / den
            out.append((j, q, p))
    return out


def precompute_survival_log2(
    *,
    entries: list[tuple[int, int, float]],
    max_live: int,
    distance: int,
    lambda_min: float,
    lambda_max: float,
    lambda_step: float,
) -> list[float]:
    lambdas = [
        lambda_min + i * lambda_step
        for i in range(int(round((lambda_max - lambda_min) / lambda_step)) + 1)
    ]
    log2_mgf = []
    for lam in lambdas:
        mgf = sum(p * math.exp(-lam * j) for j, q, p in entries if q > 0)
        log2_mgf.append(math.log2(mgf))

    out = [0.0] * (max_live + 1)
    inv_log2 = 1.0 / math.log(2.0)
    for live in range(1, max_live + 1):
        best = float("inf")
        for lam, lm in zip(lambdas, log2_mgf):
            val = lam * distance * inv_log2 + live * lm
            if val < best:
                best = val
        out[live] = min(0.0, best)
    return out


@lru_cache(maxsize=None)
def nonempty_block_coeff(x: int, h: int, b: int) -> int:
    """Coefficient [z^h] ((1+z)^b - 1)^x as an exact integer."""

    if h < x or h > b * x:
        return 0
    total = 0
    for i in range(x + 1):
        n = b * (x - i)
        if h <= n:
            term = math.comb(x, i) * math.comb(n, h)
            total = total - term if i & 1 else total + term
    return total


def occupied_logprob(*, blocks: int, h: int, occupied: int, b: int) -> float:
    if occupied < 0 or occupied > min(blocks, h):
        return float("-inf")
    if h == 0:
        return 0.0 if occupied == 0 else float("-inf")
    coeff = nonempty_block_coeff(occupied, h, b)
    if coeff <= 0:
        return float("-inf")
    return log2_binom(blocks, occupied) + math.log2(coeff) - log2_binom(b * blocks, h)


def selected_gap_sum_logcount(*, gaps: int, empty: int, selected: int, skipped: int) -> float:
    """Expected number of selected gap subsets with a given skipped length."""

    if selected < 0 or selected > gaps or skipped < 0 or skipped > empty:
        return float("-inf")
    if selected == 0:
        return 0.0 if skipped == 0 else float("-inf")
    if selected == gaps:
        return 0.0 if skipped == empty else float("-inf")

    unselected = gaps - selected
    return (
        log2_binom(gaps, selected)
        + log2_binom(skipped + selected - 1, selected - 1)
        + log2_binom(empty - skipped + unselected - 1, unselected - 1)
        - log2_binom(empty + gaps - 1, gaps - 1)
    )


def episode_gap_bound_log2(
    *,
    b: int,
    remaining_blocks: int,
    remaining_ones: int,
    p0_log2: float,
    survival_log2: list[float],
    e_max: int,
) -> tuple[float, tuple[int, int, int, float] | None]:
    total = float("-inf")
    peak: tuple[int, int, int, float] | None = None

    for occupied in range(0, min(remaining_blocks, remaining_ones) + 1):
        occ_lp = occupied_logprob(
            blocks=remaining_blocks,
            h=remaining_ones,
            occupied=occupied,
            b=b,
        )
        if occ_lp == float("-inf"):
            continue
        gaps = occupied + 1  # first active block plus remaining occupied blocks.
        empty = remaining_blocks - occupied
        max_e = min(e_max, gaps)
        for selected in range(max_e + 1):
            selected_cost = selected * p0_log2
            if selected == 0:
                term = occ_lp + survival_log2[remaining_blocks]
                total = log2add(total, term)
                if peak is None or term > peak[3]:
                    peak = (occupied, selected, 0, term)
                continue
            if selected == gaps:
                skipped_iter = (empty,)
            else:
                skipped_iter = range(empty + 1)
            for skipped in skipped_iter:
                lg = selected_gap_sum_logcount(
                    gaps=gaps,
                    empty=empty,
                    selected=selected,
                    skipped=skipped,
                )
                if lg == float("-inf"):
                    continue
                live = remaining_blocks - skipped
                term = occ_lp + lg + selected_cost + survival_log2[live]
                total = log2add(total, term)
                if peak is None or term > peak[3]:
                    peak = (occupied, selected, skipped, term)
    return total, peak


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--remaining-blocks", type=int, required=True)
    parser.add_argument("--remaining-ones", type=int, required=True)
    parser.add_argument("--e-max", type=int, default=64)
    parser.add_argument(
        "--extra-turnoff-log2",
        type=float,
        help="Optional extra per-live-branch termination atom, log2 scale; added to q'=0.",
    )
    parser.add_argument("--lambda-min", type=float, default=0.001)
    parser.add_argument("--lambda-max", type=float, default=2.0)
    parser.add_argument("--lambda-step", type=float, default=0.001)
    args = parser.parse_args()

    d = math.floor(args.distance_delta * args.N)
    entries = split_entries(load_spectrum(args.spectrum), args.block_bits)
    p0 = sum(p for _j, q, p in entries if q == 0)
    turnoff_log2 = math.log2(p0)
    if args.extra_turnoff_log2 is not None:
        turnoff_log2 = log2add(turnoff_log2, args.extra_turnoff_log2)
    survival = precompute_survival_log2(
        entries=entries,
        max_live=args.remaining_blocks,
        distance=d,
        lambda_min=args.lambda_min,
        lambda_max=args.lambda_max,
        lambda_step=args.lambda_step,
    )
    bound, peak = episode_gap_bound_log2(
        b=args.block_bits,
        remaining_blocks=args.remaining_blocks,
        remaining_ones=args.remaining_ones,
        p0_log2=turnoff_log2,
        survival_log2=survival,
        e_max=args.e_max,
    )

    print("Full-split episode-gap diagnostic")
    print(
        f"N={args.N}, b={args.block_bits}, d={d}, "
        f"remaining_blocks={args.remaining_blocks}, remaining_ones={args.remaining_ones}, e_max={args.e_max}"
    )
    print(f"self_turnoff_log2,{math.log2(p0):.6f}")
    if args.extra_turnoff_log2 is not None:
        print(f"extra_turnoff_log2,{args.extra_turnoff_log2:.6f}")
    print(f"effective_turnoff_log2,{turnoff_log2:.6f}")
    print(f"bound_log2,{bound:.6f}")
    print(f"margin_bits,{-bound:.6f}")
    if peak is not None:
        occupied, selected, skipped, term = peak
        print(f"peak_occupied_remaining,{occupied}")
        print(f"peak_selected_gaps,{selected}")
        print(f"peak_skipped_blocks,{skipped}")
        print(f"peak_term_log2,{term:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
