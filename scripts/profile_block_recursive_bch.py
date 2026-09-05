#!/usr/bin/env python3
"""Profile B_{j,q} for the output-systematic block-recursive BCH inner.

For the candidate recursion

    y_i = v_i,  s_i = P v_i,

the local object is

    B[j,q] = #{v : wt(v)=j, wt(P v)=q}.

This script computes B[j,q] exactly for small j when feasible and by sampling
for larger j.  The output is a compact CSV-style summary suitable for deciding
where a rigorous DP certificate should spend effort.
"""

from __future__ import annotations

import argparse
import math
import random

from check_bch_generator_inner import iter_weight_masks
from probe_block_recursive_inner import (
    build_output_systematic_basis,
    build_parity_tables,
    parity_state,
    parse_int_list,
)


def random_weight_mask(rng: random.Random, n: int, weight: int) -> int:
    out = 0
    for bit in rng.sample(range(n), weight):
        out |= 1 << bit
    return out


def quantile_from_hist(hist: dict[int, int], total: int, alpha: float) -> int:
    threshold = math.ceil(alpha * total)
    seen = 0
    for q in sorted(hist):
        seen += hist[q]
        if seen >= threshold:
            return q
    return max(hist)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--j-values", default="1:16,20,24,32,40,48,64")
    parser.add_argument("--exact-max-patterns", type=int, default=200000)
    parser.add_argument("--samples", type=int, default=200000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    tables = build_parity_tables(sys_basis, b)
    rng = random.Random(args.seed)

    print("Block-recursive BCH parity-expansion profile")
    print(f"m={args.m}, designed_delta={args.delta_bch}, local_length={ext_n}, block_bits={b}")
    print(
        "j,mode,count,min_q,q001,q01,q05,q10,median,q90,avg_q,max_q,"
        "le4,le8,le12,le16,le20,le24"
    )
    for j in parse_int_list(args.j_values):
        total_patterns = math.comb(b, j)
        exact = total_patterns <= args.exact_max_patterns
        count = total_patterns if exact else args.samples
        hist: dict[int, int] = {}
        q_sum = 0
        if exact:
            iterator = iter_weight_masks(b, j)
        else:
            iterator = (random_weight_mask(rng, b, j) for _ in range(args.samples))
        for v in iterator:
            q = parity_state(tables, 16, v).bit_count()
            hist[q] = hist.get(q, 0) + 1
            q_sum += q
        fields = [
            str(j),
            "exact" if exact else "sample",
            str(count),
            str(min(hist)),
            str(quantile_from_hist(hist, count, 0.001)),
            str(quantile_from_hist(hist, count, 0.01)),
            str(quantile_from_hist(hist, count, 0.05)),
            str(quantile_from_hist(hist, count, 0.10)),
            str(quantile_from_hist(hist, count, 0.50)),
            str(quantile_from_hist(hist, count, 0.90)),
            f"{q_sum / count:.6f}",
            str(max(hist)),
        ]
        for cut in [4, 8, 12, 16, 20, 24]:
            fields.append(str(sum(c for q, c in hist.items() if q <= cut)))
        print(",".join(fields))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
