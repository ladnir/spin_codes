#!/usr/bin/env python3
"""Low-state recurrence profile with inter-block state permutations.

The unpermuted block-recursive BCH inner uses the deterministic trajectory

    s, P s, P^2 s, ...

and the current [128,64,22] systematic parity map has many unit-state returns
at P^2.  This diagnostic inserts fixed random state-coordinate permutations:

    s_{t+1} = pi_t P s_t.

It checks whether the low-power algebraic recurrence disappears while
preserving the local BCH expansion.
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


def random_perm(rng: random.Random, b: int) -> list[int]:
    perm = list(range(b))
    rng.shuffle(perm)
    return perm


def apply_perm(x: int, perm: list[int]) -> int:
    out = 0
    while x:
        bit = x & -x
        src = bit.bit_length() - 1
        out |= 1 << perm[src]
        x ^= bit
    return out


def summarize(states: list[int], cuts: list[int]) -> list[str]:
    qs = [s.bit_count() for s in states]
    row = [
        str(len(states)),
        str(min(qs)),
        f"{sum(qs) / len(qs):.6f}",
        str(max(qs)),
    ]
    for cut in cuts:
        row.append(str(sum(q <= cut for q in qs)))
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--start-weights", default="1,2")
    parser.add_argument("--T-max", type=int, default=20)
    parser.add_argument("--q-cuts", default="1,2,4,8,16,20")
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--exact-mask-limit", type=int, default=200000)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    tables = build_parity_tables(sys_basis, b)
    cuts = parse_int_list(args.q_cuts)
    rng = random.Random(args.seed)

    print("Block-recursive permuted-state low-recurrence profile")
    print(
        f"m={args.m}, designed_delta={args.delta_bch}, b={b}, local_length={ext_n}, "
        f"T_max={args.T_max}, trials={args.trials}, row_d0={row_d0}"
    )
    header = ["trial", "start_w", "t", "count", "min_q", "avg_q", "max_q"]
    header.extend(f"le{cut}" for cut in cuts)
    print(",".join(header))

    for trial in range(args.trials):
        perms = [random_perm(rng, b) for _ in range(args.T_max)]
        for start_w in parse_int_list(args.start_weights):
            total = math.comb(b, start_w)
            if total > args.exact_mask_limit:
                print(f"# skipping start_w={start_w}: {total} masks exceed exact limit")
                continue
            states = list(iter_weight_masks(b, start_w))
            for t in range(args.T_max + 1):
                print(",".join([str(trial), str(start_w), str(t)] + summarize(states, cuts)))
                if t < args.T_max:
                    states = [
                        apply_perm(parity_state(tables, 16, s), perms[t]) if s else 0
                        for s in states
                    ]
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
