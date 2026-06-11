#!/usr/bin/env python3
"""Monte Carlo probe for the output-systematic block-recursive BCH inner.

The candidate inner is

    v_i = u_i + s_{i-1}
    y_i = v_i
    s_i = P v_i,

where (v, P v) is a systematic encoder for an extended primitive BCH
[2b,b,d0] code.  This script samples uniform weight-h inputs over N=M*b
coordinates and simulates the recursion.

With --state-permutation random, the state update is instead

    s_i = pi_i P v_i,

where pi_i is a fixed random coordinate permutation for block i.  This is meant
to test whether inter-block state permutations remove low-power recurrences of
the parity map.

This is evidence, not a certificate.  It is intended to answer whether the
small-h scalar fixed-tap bottleneck still looks dangerous after the BCH parity
state expansion.
"""

from __future__ import annotations

import argparse
import math
import random
from collections import defaultdict

from check_bch_generator_inner import encode, systematic_basis
from check_bch_transposed_inner import extended_bch_check_rows
from exact_bch_spectrum_small import nullspace_basis


def parse_int_list(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            vals = [int(x) for x in part.split(":")]
            if len(vals) == 2:
                lo, hi = vals
                step = 1
            elif len(vals) == 3:
                lo, hi, step = vals
            else:
                raise ValueError(f"bad range {part!r}")
            out.extend(range(lo, hi + 1, step))
        else:
            out.append(int(part))
    return out


def build_output_systematic_basis(m: int, delta: int) -> tuple[list[int], int, int, int]:
    rows, ext_n, dim, rank = extended_bch_check_rows(m, delta)
    if ext_n != 2 * dim:
        raise ValueError(f"extended code is not rate half: length={ext_n}, dim={dim}")
    basis = nullspace_basis(rows, ext_n)
    sys_basis = systematic_basis(basis, ext_n, list(range(dim)))
    d0 = min((encode(sys_basis, 1 << i).bit_count() for i in range(dim)), default=0)
    return sys_basis, dim, ext_n, d0


def build_parity_tables(sys_basis: list[int], b: int, chunk_bits: int = 16) -> list[list[int]]:
    cols = [word >> b for word in sys_basis]
    chunks = (b + chunk_bits - 1) // chunk_bits
    tables: list[list[int]] = []
    for chunk in range(chunks):
        width = min(chunk_bits, b - chunk * chunk_bits)
        table = [0] * (1 << width)
        for mask in range(1, 1 << width):
            bit = mask & -mask
            idx = bit.bit_length() - 1
            table[mask] = table[mask ^ bit] ^ cols[chunk * chunk_bits + idx]
        tables.append(table)
    return tables


def parity_state(tables: list[list[int]], chunk_bits: int, v: int) -> int:
    out = 0
    mask = (1 << chunk_bits) - 1
    for chunk, table in enumerate(tables):
        out ^= table[(v >> (chunk * chunk_bits)) & mask]
    return out


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


def sample_input_blocks(rng: random.Random, n: int, b: int, h: int) -> dict[int, int]:
    blocks: dict[int, int] = defaultdict(int)
    for pos in rng.sample(range(n), h):
        block = pos // b
        off = pos - block * b
        blocks[block] ^= 1 << off
    return blocks


def simulate_once(
    *,
    rng: random.Random,
    parity_tables: list[list[int]],
    b: int,
    blocks_n: int,
    h: int,
    d: int,
    terminal: str,
    perms: list[list[int]] | None,
) -> tuple[int, bool, int, int]:
    n = blocks_n * b
    input_blocks = sample_input_blocks(rng, n, b, h)
    if not input_blocks:
        return 0, True, blocks_n, 0

    state = 0
    total = 0
    first_active = min(input_blocks)
    active_blocks = 0
    for i in range(first_active, blocks_n):
        u = input_blocks.get(i, 0)
        v = u ^ state
        if v:
            active_blocks += 1
        total += v.bit_count()
        if total > d:
            return total, False, first_active, active_blocks
        state = parity_state(parity_tables, 16, v) if v else 0
        if state and perms is not None:
            state = apply_perm(state, perms[i])

    if terminal == "state":
        total += state.bit_count()
    elif terminal == "require-zero" and state:
        return total, False, first_active, active_blocks
    elif terminal != "none":
        raise ValueError(f"unknown terminal mode {terminal!r}")
    return total, total <= d, first_active, active_blocks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--distance-delta", type=float, default=0.106)
    parser.add_argument("--h-values", default="1,2,4,8,16,22,24,25,32,40,64")
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--terminal", choices=("none", "state", "require-zero"), default="state")
    parser.add_argument("--state-permutation", choices=("none", "random"), default="none")
    parser.add_argument("--perm-seed", type=int, default=12345)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    parity_tables = build_parity_tables(sys_basis, b)
    if args.N % b != 0:
        raise ValueError(f"N={args.N} must be divisible by b={b}")
    blocks_n = args.N // b
    d = math.floor(args.distance_delta * args.N)
    rng = random.Random(args.seed)
    perms = None
    if args.state_permutation == "random":
        perm_rng = random.Random(args.perm_seed)
        perms = [random_perm(perm_rng, b) for _ in range(blocks_n)]

    print("Block-recursive BCH inner Monte Carlo")
    print(
        f"m={args.m}, designed_delta={args.delta_bch}, local_length={ext_n}, "
        f"block_bits={b}, N={args.N}, blocks={blocks_n}, d={d}, terminal={args.terminal}, "
        f"state_permutation={args.state_permutation}, perm_seed={args.perm_seed}"
    )
    print("h,trials,failures,fail_rate,fail_log2,min_weight,avg_capped_weight,avg_first_active,avg_active_blocks")
    for h in parse_int_list(args.h_values):
        failures = 0
        min_weight = None
        weight_sum = 0
        first_sum = 0
        active_sum = 0
        for _ in range(args.trials):
            wt, ok, first, active = simulate_once(
                rng=rng,
                parity_tables=parity_tables,
                b=b,
                blocks_n=blocks_n,
                h=h,
                d=d,
                terminal=args.terminal,
                perms=perms,
            )
            failures += int(ok)
            min_weight = wt if min_weight is None else min(min_weight, wt)
            weight_sum += min(wt, d + 1)
            first_sum += first
            active_sum += active
        rate = failures / args.trials
        log_rate = math.log2(rate) if rate > 0 else float("-inf")
        print(
            f"{h},{args.trials},{failures},{rate:.8g},"
            f"{log_rate:.6f}" if rate > 0 else f"{h},{args.trials},{failures},0,-inf",
            end="",
        )
        print(
            f",{min_weight},{weight_sum / args.trials:.6f},"
            f"{first_sum / args.trials:.6f},{active_sum / args.trials:.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
