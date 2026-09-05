#!/usr/bin/env python3
"""Bucketed certificate diagnostic for the full-codeword random-split inner.

Variant:

    V_i = U_i + S_{i-1}.

If V_i=0, emit zero and set state zero.  If V_i != 0, apply a random
invertible scrambler, encode with a rate-half BCH block code, randomly split
the 2b codeword coordinates into b output and b state coordinates, and continue.

For nonzero V_i, the scrambler makes the BCH codeword uniform over nonzero
codewords.  Therefore the local transition uses only the full BCH spectrum
A_w and the hypergeometric split:

    Pr[q'=s | w] = C(b,s) C(b,w-s) / C(2b,w),
    j = w - q'.

This script mirrors certify_permuted_inner_bucket_sum.py, but avoids the
systematic split profile B_{j,q'} entirely.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from dense_largek_eval import log2_binom, log2add
from import_wd_spectrum import parse_wd
from probe_block_recursive_inner import parse_int_list
from certify_permuted_inner_chernoff import hypergeom_block_probs


def load_outer(path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if "h" in row and "outer_log2" in row:
                val = row["outer_log2"]
                if val and val != "-inf":
                    out[int(row["h"])] = float(val)
    return out


def load_spectrum(path: Path) -> list[tuple[int, int]]:
    if path.suffix.lower() == ".wd":
        return parse_wd(path.read_text())
    rows: list[tuple[int, int]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            w_key = "weight" if "weight" in row else "w"
            c_key = "count" if "count" in row else "A_w"
            rows.append((int(row[w_key]), int(row[c_key])))
    return rows


def first_active_bucket_r_log2(
    *,
    n: int,
    b: int,
    h: int,
    late_blocks: int,
    gap_min: int,
    gap_max: int,
    first_r: int,
) -> float:
    if first_r < 1 or first_r > min(b, h):
        return float("-inf")
    total = float("-inf")
    for gap in range(gap_min, gap_max + 1):
        after_coords = b * (late_blocks + gap - 1)
        if h - first_r > after_coords:
            continue
        term = (
            log2_binom(b, first_r)
            + log2_binom(after_coords, h - first_r)
            - log2_binom(n, h)
        )
        total = log2add(total, term)
    return total


def fullsplit_distribution(
    *,
    spectrum: list[tuple[int, int]],
    b: int,
    q_bin_size: int,
) -> dict[tuple[int, int], float]:
    """Distribution of (j_output, q_state) for a nonzero local branch."""

    nonzero_total = sum(c for w, c in spectrum if w > 0)
    dist: dict[tuple[int, int], float] = defaultdict(float)
    for w, count in spectrum:
        if w <= 0 or count <= 0:
            continue
        p_w = count / nonzero_total
        den = math.comb(2 * b, w)
        lo = max(0, w - b)
        hi = min(b, w)
        for q in range(lo, hi + 1):
            ways = math.comb(b, q) * math.comb(b, w - q)
            q_bin = (q // q_bin_size) * q_bin_size if q_bin_size > 1 else q
            j = w - q
            dist[(j, q_bin)] += p_w * ways / den
    # Normalize away tiny floating drift.
    s = sum(dist.values())
    if s > 0.0:
        inv = 1.0 / s
        for key in list(dist):
            dist[key] *= inv
    return dict(dist)


def split_tail_summary(dist: dict[tuple[int, int], float], cuts: list[int]) -> list[tuple[int, float]]:
    out = []
    for cut in cuts:
        p = sum(prob for (_j, q), prob in dist.items() if q <= cut)
        out.append((cut, math.log2(p) if p > 0 else float("-inf")))
    return out


def weighted_branch_vector(
    *,
    b: int,
    branch_dist: dict[tuple[int, int], float],
    lam: float,
) -> list[tuple[int, float]]:
    """Weighted q'-law for a nonzero full-split branch.

    The full-split scrambler makes the BCH codeword law independent of
    wt(V_i), as long as V_i is nonzero.  The only q/r-specific part of the
    local transition is therefore the exact-cancellation event V_i=0.
    """

    weighted_branch: dict[int, float] = defaultdict(float)
    for (j, q_next), p in branch_dist.items():
        weighted_branch[q_next] += p * math.exp(-lam * j)
    return sorted(weighted_branch.items())


def weighted_branch_array(
    *,
    b: int,
    branch_dist: dict[tuple[int, int], float],
    lam: float,
) -> np.ndarray:
    out = np.zeros(b + 1, dtype=np.float64)
    for (j, q_next), p in branch_dist.items():
        out[q_next] += p * math.exp(-lam * j)
    return out


def continuation_moment_from_states(
    *,
    b: int,
    branch_weight: np.ndarray,
    binoms: list[int],
    live_blocks_after_first: int,
    states: dict[tuple[int, int], float],
    q_support: list[int],
) -> float:
    if not states:
        return 0.0

    max_h = max(H for H, _q in states)
    cur = np.zeros((max_h + 1, b + 1), dtype=np.float64)
    for (H, q), mass in states.items():
        cur[H][q] += mass
    q_idx = np.array(q_support, dtype=np.intp)

    for step in range(live_blocks_after_first):
        remaining_blocks = live_blocks_after_first - step
        nxt = np.zeros_like(cur)
        branch_scales = np.zeros(max_h + 1, dtype=np.float64)
        p_cache: dict[int, list[float]] = {}
        row_totals = cur[:, q_idx].sum(axis=1)
        for H in range(max_h + 1):
            row_total = row_totals[H]
            if row_total == 0.0:
                continue
            row = cur[H]
            p_r = p_cache.get(H)
            if p_r is None:
                p_r = hypergeom_block_probs(b=b, remaining_blocks=remaining_blocks, H=H)
                p_cache[H] = p_r
            for r, p_occ in enumerate(p_r):
                if p_occ == 0.0:
                    continue
                target_h = H - r
                cancel_mass = row[r] / binoms[r] if r < len(row) and row[r] else 0.0
                noncancel_mass = row_total - cancel_mass
                if cancel_mass:
                    nxt[target_h][0] += p_occ * cancel_mass
                if noncancel_mass:
                    branch_scales[target_h] += p_occ * noncancel_mass
        nxt += branch_scales[:, None] * branch_weight[None, :]
        cur = nxt
    return float(cur.sum())


def open_row_writer(path: Path | None, row_fields: list[str]) -> tuple[object | None, csv.DictWriter | None]:
    if path is None:
        return None, None
    f = path.open("w", newline="")
    writer = csv.DictWriter(f, fieldnames=row_fields)
    writer.writeheader()
    f.flush()
    return f, writer


def initial_bound(
    *,
    branch_dist: dict[tuple[int, int], float],
    branch_by_lambda: dict[float, np.ndarray],
    b: int,
    binoms: list[int],
    live_blocks_after_first: int,
    remaining_h: int,
    d: int,
    lambdas: list[float],
    q_support: list[int],
) -> tuple[float, int, int, float]:
    """Chernoff bound including first active split and continuation."""

    best = float("inf")
    best_lam = lambdas[0]
    peak_j = 0
    peak_q = 0
    for lam in lambdas:
        init_states: dict[tuple[int, int], float] = defaultdict(float)
        for (j, q), p_init in branch_dist.items():
            if p_init <= 0.0:
                continue
            init_states[(remaining_h, q)] += p_init * math.exp(-lam * j)
            if p_init > branch_dist.get((peak_j, peak_q), 0.0):
                peak_j = j
                peak_q = q
        moment = continuation_moment_from_states(
            b=b,
            branch_weight=branch_by_lambda[lam],
            binoms=binoms,
            live_blocks_after_first=live_blocks_after_first,
            states=dict(init_states),
            q_support=q_support,
        )
        bound = (lam * d + math.log(moment)) / math.log(2) if moment > 0.0 else float("-inf")
        if bound < best:
            best = bound
            best_lam = lam
    return best, peak_j, peak_q, best_lam


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-prefix-csv", type=Path, required=True)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--h-values", default="32")
    parser.add_argument("--first-r-values", default="1")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=4001)
    parser.add_argument("--gap-stop", type=int, default=8000)
    parser.add_argument("--gap-step", type=int, default=4000)
    parser.add_argument("--lambdas", default="0.0001,0.0002,0.0003")
    parser.add_argument("--q-bin-size", type=int, default=16)
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    b = args.block_bits
    if args.N % b:
        raise ValueError("N must be divisible by block size")
    d = math.floor(args.distance_delta * args.N)
    outer = load_outer(args.outer_prefix_csv)
    spectrum = load_spectrum(args.spectrum)
    total_words = sum(c for _, c in spectrum)
    d0 = min(w for w, c in spectrum if w > 0 and c > 0)
    raw_branch_dist = fullsplit_distribution(spectrum=spectrum, b=b, q_bin_size=1)
    branch_dist = fullsplit_distribution(spectrum=spectrum, b=b, q_bin_size=args.q_bin_size)
    q_support = sorted({0, *(q for _j, q in branch_dist)})
    lambdas = [float(x) for x in args.lambdas.split(",") if x.strip()]
    branch_by_lambda = {
        lam: weighted_branch_array(b=b, branch_dist=branch_dist, lam=lam)
        for lam in lambdas
    }
    binoms = [math.comb(b, i) for i in range(b + 1)]

    print("Full-codeword random-split BCH inner bucket diagnostic")
    print(
        f"N={args.N}, b={b}, d={d}, late_blocks={args.late_blocks}, "
        f"spectrum={args.spectrum}, spectrum_words={total_words}, d0={d0}, "
        f"q_bin_size={args.q_bin_size}"
    )
    print(
        "split_tail_log2="
        + ",".join(f"q<={cut}:{val:.6f}" for cut, val in split_tail_summary(raw_branch_dist, [0, 4, 8, 16]))
    )
    row_fields = [
        "h",
        "first_r",
        "gap_min",
        "gap_max",
        "outer_log2",
        "placement_log2",
        "inner_log2",
        "term_log2",
        "peak_j",
        "peak_q",
        "peak_lambda",
        "live_after_first",
    ]
    csv_file, csv_writer = open_row_writer(args.output_csv, row_fields)
    if not args.summary_only:
        print(",".join(row_fields))

    grand = float("-inf")
    peak: tuple[int, int, int, int, float] | None = None
    for h in parse_int_list(args.h_values):
        if h not in outer:
            continue
        for first_r in parse_int_list(args.first_r_values):
            if first_r > h:
                continue
            remaining_h = h - first_r
            for gap_min in range(args.gap_start, args.gap_stop + 1, args.gap_step):
                gap_max = min(args.gap_stop, gap_min + args.gap_step - 1)
                placement = first_active_bucket_r_log2(
                    n=args.N,
                    b=b,
                    h=h,
                    late_blocks=args.late_blocks,
                    gap_min=gap_min,
                    gap_max=gap_max,
                    first_r=first_r,
                )
                if placement == float("-inf"):
                    continue
                live_after = args.late_blocks + gap_min - 1
                inner, peak_j, peak_q, peak_lam = initial_bound(
                    branch_dist=branch_dist,
                    branch_by_lambda=branch_by_lambda,
                    b=b,
                    binoms=binoms,
                    live_blocks_after_first=live_after,
                    remaining_h=remaining_h,
                    d=d,
                    lambdas=lambdas,
                    q_support=q_support,
                )
                term = outer[h] + placement + min(0.0, inner)
                grand = log2add(grand, term)
                if peak is None or term > peak[4]:
                    peak = (h, first_r, gap_min, gap_max, term)
                row = {
                    "h": h,
                    "first_r": first_r,
                    "gap_min": gap_min,
                    "gap_max": gap_max,
                    "outer_log2": outer[h],
                    "placement_log2": placement,
                    "inner_log2": inner,
                    "term_log2": term,
                    "peak_j": peak_j,
                    "peak_q": peak_q,
                    "peak_lambda": peak_lam,
                    "live_after_first": live_after,
                }
                if csv_writer is not None and csv_file is not None:
                    csv_writer.writerow(row)
                    csv_file.flush()
                if not args.summary_only:
                    print(
                        f"{h},{first_r},{gap_min},{gap_max},{outer[h]:.6f},"
                        f"{placement:.6f},{inner:.6f},{term:.6f},"
                        f"{peak_j},{peak_q},{peak_lam:.8g},{live_after}"
                        ,
                        flush=True,
                    )

    print("summary")
    print(f"total_log2,{grand:.6f}")
    print(f"margin_bits,{-grand:.6f}")
    if peak is not None:
        h, first_r, gap_min, gap_max, term = peak
        print(f"peak_h,{h}")
        print(f"peak_first_r,{first_r}")
        print(f"peak_gap_min,{gap_min}")
        print(f"peak_gap_max,{gap_max}")
        print(f"peak_term_log2,{term:.6f}")
    if csv_file is not None:
        csv_file.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
