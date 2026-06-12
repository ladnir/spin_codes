#!/usr/bin/env python3
"""Column-distance probes for the block-recursive BCH inner.

For the zero-input live recursion

    y_t = s_{t-1},  s_t = P s_{t-1},

the relevant growth statistic is

    D_L(P) = min_{s != 0} sum_{t=0}^{L-1} wt(P^t s).

This script computes D_L exactly for small b when 2^b is feasible and by
sampling for b=64.  It is a diagnostic for the proof obligation: a good
block-recursive inner needs a usable lower bound on D_L, not just the local
two-block distance wt(s)+wt(Ps).
"""

from __future__ import annotations

import argparse
import math
import random

from check_bch_generator_inner import iter_weight_masks
from probe_block_recursive_inner import build_output_systematic_basis, build_parity_tables, parity_state


def random_mask(rng: random.Random, b: int, mode: str, weight: int | None) -> int:
    if mode == "uniform":
        x = 0
        while x == 0:
            x = rng.getrandbits(b)
        return x
    if mode == "weight":
        if weight is None:
            raise ValueError("--sample-weight is required for weight mode")
        out = 0
        for bit in rng.sample(range(b), weight):
            out |= 1 << bit
        return out
    raise ValueError(f"unknown sample mode {mode!r}")


def trajectory_sums(tables: list[list[int]], b: int, s0: int, l_max: int) -> list[int]:
    out = [0] * (l_max + 1)
    s = s0
    acc = 0
    for ell in range(1, l_max + 1):
        acc += s.bit_count()
        out[ell] = acc
        s = parity_state(tables, 16, s)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--delta-bch", type=int, default=7)
    parser.add_argument("--L-max", type=int, default=32)
    parser.add_argument("--exact-max-states", type=int, default=1 << 20)
    parser.add_argument("--samples", type=int, default=200000)
    parser.add_argument("--sample-mode", choices=("uniform", "weight"), default="uniform")
    parser.add_argument("--sample-weight", type=int, default=None)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    tables = build_parity_tables(sys_basis, b)
    total_states = (1 << b) - 1
    if args.sample_mode == "weight" and args.sample_weight is not None:
        candidate_count = math.comb(b, args.sample_weight)
        exact = candidate_count <= args.exact_max_states
        count = candidate_count if exact else args.samples
    else:
        exact = total_states <= args.exact_max_states
        count = total_states if exact else args.samples
    rng = random.Random(args.seed)

    mins = [10**18] * (args.L_max + 1)
    sums = [0] * (args.L_max + 1)
    argmins = [0] * (args.L_max + 1)

    if exact and args.sample_mode == "weight" and args.sample_weight is not None:
        iterator = iter_weight_masks(b, args.sample_weight)
        mode = f"exact-weight-{args.sample_weight}"
    elif exact:
        iterator = range(1, 1 << b)
        mode = "exact"
    else:
        iterator = (random_mask(rng, b, args.sample_mode, args.sample_weight) for _ in range(args.samples))
        mode = f"sample-{args.sample_mode}"

    for idx, s0 in enumerate(iterator, 1):
        vals = trajectory_sums(tables, b, s0, args.L_max)
        for ell in range(1, args.L_max + 1):
            val = vals[ell]
            sums[ell] += val
            if val < mins[ell]:
                mins[ell] = val
                argmins[ell] = s0

    print("Block-recursive zero-input column profile")
    print(
        f"m={args.m}, designed_delta={args.delta_bch}, local_length={ext_n}, block_bits={b}, "
        f"mode={mode}, count={count}, sample_weight={args.sample_weight}"
    )
    print("L,min_sum,avg_sum,min_per_block,avg_per_block,argmin_weight,pair_bound")
    for ell in range(1, args.L_max + 1):
        pair_bound = (ell // 2) * row_d0
        if ell % 2:
            pair_bound += 1
        print(
            f"{ell},{mins[ell]},{sums[ell] / count:.6f},"
            f"{mins[ell] / ell:.6f},{(sums[ell] / count) / ell:.6f},"
            f"{argmins[ell].bit_count()},{pair_bound}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
