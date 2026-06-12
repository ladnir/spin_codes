#!/usr/bin/env python3
"""Low-state recurrence profile for the block-recursive BCH inner.

Short recurrences such as wt(P^2 e_i)=1 create cheap exact-turnoff pairs.  This
script profiles how often a low-weight initial state returns to a low-weight
state after a small number of zero-input steps.  The output is intended to size
the pair/cluster counting lemma for the early branch.
"""

from __future__ import annotations

import argparse
import math

from check_bch_generator_inner import iter_weight_masks
from probe_block_recursive_inner import (
    build_output_systematic_basis,
    build_parity_tables,
    parity_state,
    parse_int_list,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--start-weights", default="1,2")
    parser.add_argument("--T-max", type=int, default=20)
    parser.add_argument("--q-cuts", default="1,2,4,8,16,20")
    parser.add_argument("--exact-mask-limit", type=int, default=200000)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    tables = build_parity_tables(sys_basis, b)
    cuts = parse_int_list(args.q_cuts)

    print("Block-recursive low-state recurrence profile")
    print(
        f"m={args.m}, designed_delta={args.delta_bch}, b={b}, local_length={ext_n}, "
        f"T_max={args.T_max}, row_d0={row_d0}"
    )
    header = ["start_w", "t", "count", "min_q", "avg_q", "max_q"]
    header.extend(f"le{cut}" for cut in cuts)
    print(",".join(header))

    for start_w in parse_int_list(args.start_weights):
        total = math.comb(b, start_w)
        if total > args.exact_mask_limit:
            print(f"# skipping start_w={start_w}: {total} masks exceed exact limit")
            continue
        states = list(iter_weight_masks(b, start_w))
        for t in range(args.T_max + 1):
            qs = [s.bit_count() for s in states]
            row = [
                str(start_w),
                str(t),
                str(len(states)),
                str(min(qs)),
                f"{sum(qs) / len(qs):.6f}",
                str(max(qs)),
            ]
            for cut in cuts:
                row.append(str(sum(q <= cut for q in qs)))
            print(",".join(row))
            states = [parity_state(tables, 16, s) if s else 0 for s in states]
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
