"""Corrected r=1 burst+late-tail diagnostic.

The older paired-cover envelope allows late-tail support after an early
terminated burst without charging the late episode. This diagnostic fixes that
for the representative r=1, fixed-T lane: after the early burst terminates, the
first late support position is charged by either a late fixed-tap termination or
by the binomial lower tail over its suffix length.
"""

from __future__ import annotations

import argparse
import csv
import math
from functools import lru_cache

from dense_largek_eval import (
    binom_cdf_half_entropy_log2,
    log2_binom,
    log2add,
    outer_small_h_prefix,
)


def format_log2(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.6f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--sigma", type=int, required=True)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--h-max", type=int, default=80)
    parser.add_argument("--small-xi", type=float, default=0.1)
    parser.add_argument("--t-factor", type=float, default=2.0)
    parser.add_argument("--outer-csv", default=None, help="Optional plot_fcib_peaks CSV with outer_log2 column.")
    parser.add_argument("--outer-smallh-mode", choices=("banded-fixedtap",), default="banded-fixedtap")
    args = parser.parse_args()

    n = 2 * args.k
    d = math.floor(args.delta * n)
    t = math.ceil(args.t_factor * d)
    l = math.ceil((2.0 + args.small_xi) * d)
    prefix = n - l
    if t > prefix:
        raise SystemExit("T exceeds early prefix")

    if args.outer_csv:
        outer = [float("-inf")] * (args.h_max + 1)
        with open(args.outer_csv, newline="") as f:
            for row in csv.DictReader(f):
                if int(row["sigma"]) == args.sigma:
                    h = int(row["h"])
                    if 0 <= h <= args.h_max:
                        outer[h] = float(row["outer_log2"])
    else:
        outer = outer_small_h_prefix(
            args.k,
            args.sigma,
            args.h_max,
            outer_mode=args.outer_smallh_mode,
            parity_n=args.k,
        )
    early_cost = -(args.sigma - 1) + binom_cdf_half_entropy_log2(max(0, t - args.sigma), d)

    @lru_cache(maxsize=None)
    def late_sum_log(j: int) -> float:
        if j == 0:
            return 0.0
        total = float("-inf")
        logc = 0.0  # C(s-1,j-1) at s=j.
        cand = (j - 1) * (2.0 ** (-(args.sigma - 1)))
        for s in range(j, l + 1):
            tail_log = binom_cdf_half_entropy_log2(s, d)
            tail = 1.0 if tail_log >= 0.0 else 2.0**tail_log
            late_prob = min(1.0, cand + tail)
            total = log2add(total, logc + math.log2(late_prob))
            if s < l:
                logc += math.log2(s) - math.log2(s - j + 1)
        return total

    cumulative = float("-inf")
    peak = (float("-inf"), 0)
    print(f"k={args.k} sigma={args.sigma} n={n} d={d} T={t} L={l}")
    print("h, outer_log2, corrected_inner_place_log2, term_log2, cumulative_log2")
    for h in range(2, args.h_max + 1):
        denom = log2_binom(n, h)
        inner_place = float("-inf")
        for j in range(0, h - 1):
            i = h - 2 - j
            if i < 0 or i > t - 2:
                continue
            base = (
                log2_binom(prefix - t + 1, 1)
                + log2_binom(t - 2, i)
                - denom
                + early_cost
            )
            inner_place = log2add(inner_place, base + late_sum_log(j))
        term = outer[h] + inner_place
        cumulative = log2add(cumulative, term)
        if term > peak[0]:
            peak = (term, h)
        print(f"{h}, {format_log2(outer[h])}, {format_log2(inner_place)}, {format_log2(term)}, {format_log2(cumulative)}")
    print(f"TOTAL_LOG2={format_log2(cumulative)}")
    print(f"PEAK_H={peak[1]} PEAK_LOG2={format_log2(peak[0])}")


if __name__ == "__main__":
    main()
