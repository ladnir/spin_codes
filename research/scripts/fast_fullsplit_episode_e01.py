#!/usr/bin/env python3
"""Fast episode-gap scan for the full-split BCH inner.

The full episode diagnostic sums over selected skipped-gap subsets.  This
helper uses the closed form for a fixed selected-gap count e.  For e=1,

  P[X=x] E[# one-gap choices with skipped length s | X=x]
    = coeff(H,x) / C(bT,H) * (x+1) * C(T-s-1,x-1),

where coeff(H,x) = [z^H]((1+z)^b-1)^x.

For general 1 <= e <= x, the corresponding factor is

  coeff(H,x) / C(bT,H) * C(x+1,e) C(s+e-1,e-1) C(T-s-e,x-e).

The endpoint e=x+1 is handled separately.  The script reports the e=0 term,
the e=1..E terms, and their log-sum for a list of remaining input weights H.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from bound_fullsplit_episode_gaps import (
    load_spectrum,
    precompute_survival_log2,
    split_entries,
)
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list


def precompute_logc_nk(*, n_max: int, k_max: int) -> list[list[float]]:
    table: list[list[float]] = []
    for k in range(k_max + 1):
        row = [float("-inf")] * (n_max + 1)
        if k == 0:
            for n in range(n_max + 1):
                row[n] = 0.0
        elif k <= n_max:
            row[k] = 0.0
            for n in range(k + 1, n_max + 1):
                row[n] = row[n - 1] + math.log2(n / (n - k))
        table.append(row)
    return table


def precompute_nonempty_block_logcoeff(
    *,
    b: int,
    h_max: int,
    x_max: int,
) -> list[list[float]]:
    """Log coefficients of ((1+z)^b-1)^x through z^h_max.

    The inclusion-exclusion integer formula is exact, but slow when auditing
    many H values.  This recurrence is the same coefficient identity, evaluated
    in log space and truncated to the H range needed by the certificate.
    """

    choose = [float("-inf")] + [log2_binom(b, r) for r in range(1, min(b, h_max) + 1)]
    coeff = [[float("-inf")] * (h_max + 1) for _ in range(x_max + 1)]
    coeff[0][0] = 0.0
    for x in range(1, x_max + 1):
        prev = coeff[x - 1]
        cur = coeff[x]
        h_hi = min(h_max, b * x)
        for h in range(x, h_hi + 1):
            total = float("-inf")
            r_lo = max(1, h - b * (x - 1))
            r_hi = min(b, h)
            for r in range(r_lo, r_hi + 1):
                prev_val = prev[h - r]
                if prev_val != float("-inf"):
                    total = log2add(total, prev_val + choose[r])
            cur[h] = total
    return coeff


def precompute_gap_suffix_logsum(
    *,
    T: int,
    x_max: int,
    e_max: int,
    survival_log2: list[float],
    logc: list[list[float]],
) -> list[list[float]]:
    """Precompute skipped-gap suffix sums for each selected-gap count.

    For fixed x and 1 <= e <= x, the selected-gap term contains

        C(x+1,e) sum_s C(s+e-1,e-1) C(T-s-e,x-e) survival(T-s).

    This depends on x, e, and T, but not on H.  For e=x+1 all gaps are
    selected, all empty blocks are skipped, and the C(T,x) occupancy factor no
    longer cancels, so the suffix is C(T,x) survival(x).
    """

    suffix = [[float("-inf")] * (x_max + 1) for _ in range(e_max + 1)]
    for x in range(0, x_max + 1):
        for e in range(1, min(e_max, x + 1) + 1):
            if e == x + 1:
                suffix[e][x] = log2_binom(T, x) + survival_log2[x]
                continue
            total = float("-inf")
            choose_gaps = log2_binom(x + 1, e)
            row_left = logc[e - 1]
            row_right = logc[x - e]
            for skipped in range(0, T - x + 1):
                left = row_left[skipped + e - 1]
                right = row_right[T - skipped - e]
                if left == float("-inf") or right == float("-inf"):
                    continue
                total = log2add(total, choose_gaps + left + right + survival_log2[T - skipped])
            suffix[e][x] = total
    return suffix


def episode_terms_log2(
    *,
    b: int,
    T: int,
    H: int,
    turnoff_log2: float,
    survival_log2: list[float],
    suffix_logsum: list[list[float]],
    log_coeff: list[list[float]],
) -> list[float]:
    if H < 0:
        return []
    denom = log2_binom(b * T, H)
    e_max = len(suffix_logsum) - 1
    terms = [float("-inf")] * (e_max + 1)
    terms[0] = survival_log2[T]
    x_min = 0 if H == 0 else max(1, (H + b - 1) // b)
    x_max = min(H, T)
    for x in range(x_min, x_max + 1):
        coeff_log2 = log_coeff[x][H]
        if coeff_log2 == float("-inf"):
            continue
        base = coeff_log2 - denom
        for e in range(1, min(e_max, x + 1) + 1):
            terms[e] = log2add(terms[e], base + e * turnoff_log2 + suffix_logsum[e][x])
    return terms


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("EBCH128_64.wd"))
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--remaining-blocks", type=int, required=True)
    parser.add_argument("--remaining-ones", required=True)
    parser.add_argument("--extra-turnoff-log2", type=float)
    parser.add_argument("--lambda-min", type=float, default=0.001)
    parser.add_argument("--lambda-max", type=float, default=2.0)
    parser.add_argument("--lambda-step", type=float, default=0.004)
    parser.add_argument("--e-max", type=int, default=1)
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    b = args.block_bits
    T = args.remaining_blocks
    H_values = parse_int_list(args.remaining_ones)
    d = math.floor(args.distance_delta * args.N)
    entries = split_entries(load_spectrum(args.spectrum), b)
    p0 = sum(p for _j, q, p in entries if q == 0)
    turnoff_log2 = math.log2(p0)
    if args.extra_turnoff_log2 is not None:
        turnoff_log2 = log2add(turnoff_log2, args.extra_turnoff_log2)
    survival = precompute_survival_log2(
        entries=entries,
        max_live=T,
        distance=d,
        lambda_min=args.lambda_min,
        lambda_max=args.lambda_max,
        lambda_step=args.lambda_step,
    )
    max_h = max(H_values)
    x_max = min(max_h, T)
    logc = precompute_logc_nk(n_max=T - 1, k_max=x_max)
    suffix_logsum = precompute_gap_suffix_logsum(
        T=T,
        x_max=x_max,
        e_max=args.e_max,
        survival_log2=survival,
        logc=logc,
    )
    log_coeff = precompute_nonempty_block_logcoeff(b=b, h_max=max_h, x_max=x_max)

    rows = []
    for H in H_values:
        e_terms = episode_terms_log2(
            b=b,
            T=T,
            H=H,
            turnoff_log2=turnoff_log2,
            survival_log2=survival,
            suffix_logsum=suffix_logsum,
            log_coeff=log_coeff,
        )
        total = float("-inf")
        for val in e_terms:
            total = log2add(total, val)
        row = {"H": H, "total_log2": total}
        for e, val in enumerate(e_terms):
            row[f"e{e}_log2"] = val
        rows.append(row)

    if args.output_csv:
        with args.output_csv.open("w", newline="") as f:
            fields = ["H"] + [f"e{e}_log2" for e in range(args.e_max + 1)] + ["total_log2"]
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    print("Fast full-split episode scan")
    print(f"N={args.N}, b={b}, T={T}, d={d}, e_max={args.e_max}, turnoff_log2={turnoff_log2:.6f}")
    print(",".join(["H"] + [f"e{e}_log2" for e in range(args.e_max + 1)] + ["total_log2"]))
    for row in rows:
        vals = [str(row["H"])]
        vals.extend(f"{row[f'e{e}_log2']:.6f}" for e in range(args.e_max + 1))
        vals.append(f"{row['total_log2']:.6f}")
        print(",".join(vals))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
