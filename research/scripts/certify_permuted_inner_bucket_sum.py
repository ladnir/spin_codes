#!/usr/bin/env python3
"""Bucketed conditional first-moment sum for the permuted BCH inner.

This diagnostic combines:

* exact first-active-block placement probabilities for a weight-h input;
* a conditional Chernoff continuation bound from the guessed/local profile;
* an outer spectrum coefficient.

For a bucket of gaps g_min..g_max before the late threshold L, the script sums
the exact placement mass in that bucket, but uses the continuation bound at the
near edge g_min.  That is the weakest point in the intended proof split and is
a conservative diagnostic if the continuation bound improves with live length.
The script prints this convention explicitly; the final certificate should
either prove monotonicity or use sufficiently fine buckets/endpoints.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from certify_permuted_inner_chernoff import (
    build_envelope_profile,
    build_output_systematic_basis,
    build_parity_tables,
    build_profile,
    chernoff_continuation_with_local_fast,
    chernoff_continuation_sparse_fast,
    precompute_local_transitions_prob,
    quantize_profile,
)
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list


def load_outer(path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if "h" in row and "outer_log2" in row:
                val = row["outer_log2"]
                if val and val != "-inf":
                    out[int(row["h"])] = float(val)
    return out


def log2diff(a: float, b: float) -> float:
    if b == float("-inf"):
        return a
    if b > a:
        raise ValueError("log2diff requires a >= b")
    if a == b:
        return float("-inf")
    return a + math.log2(1.0 - 2.0 ** (b - a))


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
    """Exact probability for first occupancy r and gap in [gap_min,gap_max]."""

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


def best_continuation(
    *,
    profile: list[dict[int, float]],
    b: int,
    live_blocks_after_first: int,
    remaining_h: int,
    remaining_d: int,
    q_values: list[int],
    lambdas: list[float],
    local_by_lambda: dict[float, list[list[list[tuple[int, float]]]]],
    sparse: bool,
) -> tuple[float, int, float]:
    worst_q = q_values[0]
    best_for_worst = float("-inf")
    best_lam_for_worst = lambdas[0]
    for q in q_values:
        best = float("inf")
        best_lam = lambdas[0]
        for lam in lambdas:
            fn = chernoff_continuation_sparse_fast if sparse else chernoff_continuation_with_local_fast
            val = fn(
                    local=local_by_lambda[lam],
                    b=b,
                    live_blocks_after_first=live_blocks_after_first,
                    remaining_h=remaining_h,
                    initial_q=q,
                    remaining_d=remaining_d,
                    lam=lam,
                )
            if val < best:
                best = val
                best_lam = lam
        # Need an upper bound uniform over q in this support, so take max.
        if best > best_for_worst:
            best_for_worst = best
            worst_q = q
            best_lam_for_worst = best_lam
    return best_for_worst, worst_q, best_lam_for_worst


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-prefix-csv", type=Path, required=True)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--h-values", default="32")
    parser.add_argument("--first-r-values", default="1")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=1)
    parser.add_argument("--gap-stop", type=int, default=12000)
    parser.add_argument("--gap-step", type=int, default=1000)
    parser.add_argument("--lambdas", default="0.0001,0.0002,0.0003")
    parser.add_argument("--profile-mode", choices=("exact-prefix", "random-tail", "entropy-tail", "floor"), default="random-tail")
    parser.add_argument("--profile-slack-bits", type=float, default=0.0)
    parser.add_argument("--exact-j-max", type=int, default=4)
    parser.add_argument("--q-bin-size", type=int, default=16)
    parser.add_argument("--floor-only-q", type=int, default=32)
    parser.add_argument("--q-start-mode", choices=("all", "min"), default="min")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--dense-dp", action="store_true")
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    if args.N % b:
        raise ValueError("N must be divisible by block size")
    blocks = args.N // b
    d = math.floor(args.distance_delta * args.N)
    outer = load_outer(args.outer_prefix_csv)
    lambdas = [float(x) for x in args.lambdas.split(",") if x.strip()]

    if args.profile_mode == "exact-prefix":
        tables = build_parity_tables(sys_basis, b)
        raw = build_profile(tables, b, args.exact_j_max)
        binoms = [math.comb(b, i) for i in range(b + 1)]
        profile = [{q: cnt / binoms[j] for q, cnt in raw[j].items()} for j in range(b + 1)]
    else:
        profile = build_envelope_profile(
            b=b,
            d0=row_d0,
            mode=args.profile_mode,
            slack_bits=args.profile_slack_bits,
            floor_only_q=args.floor_only_q,
        )
    profile = quantize_profile(profile, args.q_bin_size)
    local_by_lambda = {
        lam: precompute_local_transitions_prob(profile, b, lam)
        for lam in lambdas
    }

    print("Permuted block-recursive bucketed conditional first-moment diagnostic")
    print(
        f"N={args.N}, blocks={blocks}, b={b}, d={d}, late_blocks={args.late_blocks}, "
        f"profile_mode={args.profile_mode}, slack_bits={args.profile_slack_bits}, "
        f"q_bin_size={args.q_bin_size}, convention=continuation_at_bucket_near_edge"
    )
    if not args.summary_only:
        print(
            "h,first_r,gap_min,gap_max,outer_log2,placement_log2,cont_log2,"
            "term_log2,worst_q,best_lambda,live_after_first,remaining_d"
        )
    grand = float("-inf")
    peak: tuple[int, int, int, int, float] | None = None
    for h in parse_int_list(args.h_values):
        if h not in outer:
            continue
        for first_r in parse_int_list(args.first_r_values):
            if first_r > h:
                continue
            # Pessimistic support for q after the first block under the local profile.
            q_values = sorted(profile[first_r])
            if args.q_start_mode == "min" and q_values:
                q_values = [q_values[0]]
            remaining_h = h - first_r
            remaining_d = d - first_r
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
                cont, worst_q, best_lam = best_continuation(
                    profile=profile,
                    b=b,
                    live_blocks_after_first=live_after,
                    remaining_h=remaining_h,
                    remaining_d=remaining_d,
                    q_values=q_values,
                    lambdas=lambdas,
                    local_by_lambda=local_by_lambda,
                    sparse=not args.dense_dp,
                )
                term = outer[h] + placement + min(0.0, cont)
                grand = log2add(grand, term)
                if peak is None or term > peak[4]:
                    peak = (h, first_r, gap_min, gap_max, term)
                if not args.summary_only:
                    print(
                        f"{h},{first_r},{gap_min},{gap_max},{outer[h]:.6f},"
                        f"{placement:.6f},{cont:.6f},{term:.6f},{worst_q},{best_lam:.8g},"
                        f"{live_after},{remaining_d}"
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
