#!/usr/bin/env python3
"""Exact late-start thresholds for low-weight BCH recursive states.

For initial mismatch weights j small enough to enumerate, compute the largest
number of zero-input live blocks whose emitted weight can remain below a
distance budget d.  This calibrates the block-level late-start term using the
actual P trajectory instead of the crude b/2 average.
"""

from __future__ import annotations

import argparse
import math

from check_bch_generator_inner import iter_weight_masks
from probe_block_recursive_inner import build_output_systematic_basis, build_parity_tables, parity_state, parse_int_list


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--distance-deltas", default="0.106,0.09")
    parser.add_argument("--j-values", default="1,2")
    parser.add_argument("--max-live-blocks", type=int, default=20000)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    tables = build_parity_tables(sys_basis, b)
    if args.N % b != 0:
        raise ValueError("N must be divisible by block size")
    blocks = args.N // b
    deltas = [float(x) for x in args.distance_deltas.split(",") if x.strip()]

    print("Block-recursive exact low-weight late thresholds")
    print(f"m={args.m}, designed_delta={args.delta_bch}, local_length={ext_n}, block_bits={b}, N={args.N}, blocks={blocks}")
    print("j,count,delta,d,max_live_blocks_ok,live_fraction,log2_fraction")
    for j in parse_int_list(args.j_values):
        masks = list(iter_weight_masks(b, j))
        for delta in deltas:
            d = math.floor(delta * args.N)
            max_ok = 0
            for s0 in masks:
                s = s0
                acc = 0
                live = 0
                while live < args.max_live_blocks:
                    acc += s.bit_count()
                    if acc > d:
                        break
                    live += 1
                    s = parity_state(tables, 16, s)
                max_ok = max(max_ok, live)
            frac = max_ok / blocks
            print(f"{j},{len(masks)},{delta},{d},{max_ok},{frac:.12g},{math.log2(frac):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
