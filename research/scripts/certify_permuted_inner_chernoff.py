#!/usr/bin/env python3
"""Chernoff transfer certificate for the permuted block-recursive BCH inner.

For the permuted-state recursion

    V_i = U_i + S_{i-1},   Y_i = V_i,   S_i = pi_i P V_i,

the state support is uniform conditional on its weight.  This script builds the
exact weight transition kernel for one block:

    current state weight q,
    r input ones in the current block,
    ell overlaps with the state,
    j = q + r - 2 ell output ones,
    next state weight q' distributed by B[j,q'] / C(b,j).

It then computes an exponential-moment bound for a fixed total input weight h:

    Pr[wt(Y) <= d] <= exp(lambda d) E[exp(-lambda wt(Y))].

The input support is sampled without replacement across the N=Mb coordinates,
so the DP also tracks remaining input weight H and remaining blocks R.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict

from check_bch_generator_inner import iter_weight_masks
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import (
    build_output_systematic_basis,
    build_parity_tables,
    parity_state,
    parse_int_list,
)


def binary_entropy(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


def log2_random_tail(b: int, q: int) -> float:
    if q < 0:
        return float("-inf")
    if q >= b:
        return 0.0
    # Exact enough for b=64, and still cheap for moderate b.
    total = sum(math.comb(b, i) for i in range(q + 1))
    return math.log2(total) - b


def build_envelope_profile(
    *,
    b: int,
    d0: int,
    mode: str,
    slack_bits: float,
    floor_only_q: int,
) -> list[dict[int, float]]:
    """Return cumulative-tail envelope as a synthetic distribution.

    The returned profile[j] is a probability distribution over q_next that
    stochastically dominates the true low-tail behavior.  For an upper bound on
    E exp(-lambda q_next), we push as much mass as allowed by the cumulative
    envelope onto the smallest q values.
    """

    profiles: list[dict[int, float]] = []
    for j in range(b + 1):
        min_q = max(0, d0 - j) if j else 0
        cdf: list[float] = []
        for q in range(b + 1):
            if q < min_q:
                cdf.append(0.0)
            elif mode == "floor":
                cdf.append(1.0 if q >= floor_only_q else 0.0)
            elif mode == "random-tail":
                val = min(0.0, slack_bits + log2_random_tail(b, q))
                cdf.append(min(1.0, 2.0**val))
            elif mode == "entropy-tail":
                if q >= b:
                    cdf.append(1.0)
                else:
                    val = min(0.0, slack_bits - b * (1.0 - binary_entropy(q / b)))
                    cdf.append(min(1.0, 2.0**val))
            else:
                raise ValueError(f"unknown envelope mode {mode!r}")
        cdf[-1] = 1.0
        prev = 0.0
        dist: dict[int, float] = {}
        for q, cur in enumerate(cdf):
            cur = max(prev, min(1.0, cur))
            mass = cur - prev
            if mass > 0.0:
                dist[q] = mass
            prev = cur
        profiles.append(dist)
    return profiles


def quantize_profile(profile: list[dict[int, float]], bin_size: int) -> list[dict[int, float]]:
    if bin_size <= 1:
        return profile
    out: list[dict[int, float]] = []
    for row in profile:
        acc: dict[int, float] = defaultdict(float)
        for q, mass in row.items():
            # Round down: lower state weight is pessimistic for future output.
            acc[(q // bin_size) * bin_size] += mass
        out.append(dict(acc))
    return out


def build_profile(tables: list[list[int]], b: int, exact_j_max: int) -> list[dict[int, int]]:
    prof: list[dict[int, int]] = [defaultdict(int) for _ in range(b + 1)]
    for j in range(b + 1):
        if j <= exact_j_max or j == b:
            for mask in iter_weight_masks(b, j):
                q = parity_state(tables, 16, mask).bit_count() if mask else 0
                prof[j][q] += 1
        else:
            # Conservative fallback: forget expansion for unprofiled output
            # weights.  This preserves a valid upper bound but can be loose.
            prof[j][0] = math.comb(b, j)
    return [dict(row) for row in prof]


def block_input_logprob(*, b: int, remaining_blocks: int, H: int, r: int) -> float:
    if r < 0 or r > b or r > H:
        return float("-inf")
    rest = b * (remaining_blocks - 1)
    if H - r > rest:
        return float("-inf")
    return log2_binom(b, r) + log2_binom(rest, H - r) - log2_binom(b * remaining_blocks, H)


def overlap_logprob(*, b: int, q: int, r: int, ell: int) -> float:
    if ell < 0 or ell > q or ell > r:
        return float("-inf")
    if r - ell > b - q:
        return float("-inf")
    return log2_binom(q, ell) + log2_binom(b - q, r - ell) - log2_binom(b, r)


def log2_profile_prob(profile: list[dict[int, int]], b: int, j: int, q_next: int) -> float:
    cnt = profile[j].get(q_next, 0)
    if cnt == 0:
        return float("-inf")
    return math.log2(cnt) - log2_binom(b, j)


def precompute_local_transitions(
    profile: list[dict[int, int]],
    b: int,
    lam: float,
) -> list[list[list[tuple[int, int, float]]]]:
    """local[q][r] -> [(j, q_next, log2 prob including exp charge)]."""

    log2_e = math.log2(math.e)
    local: list[list[list[tuple[int, int, float]]]] = [[[] for _ in range(b + 1)] for _ in range(b + 1)]
    for q in range(b + 1):
        for r in range(b + 1):
            rows: list[tuple[int, int, float]] = []
            lo = max(0, r - (b - q))
            hi = min(q, r)
            for ell in range(lo, hi + 1):
                j = q + r - 2 * ell
                lp_ell = overlap_logprob(b=b, q=q, r=r, ell=ell)
                charge = -lam * j * log2_e
                for q_next in profile[j]:
                    lp_q = log2_profile_prob(profile, b, j, q_next)
                    rows.append((j, q_next, lp_ell + lp_q + charge))
            local[q][r] = rows
    return local


def chernoff_for_h(
    *,
    profile: list[dict[int, int]],
    b: int,
    blocks: int,
    h: int,
    d: int,
    lam: float,
    stop_after_blocks: int | None,
) -> float:
    # map (H_remaining, q_state) -> log2 weighted mass.
    states: dict[tuple[int, int], float] = {(h, 0): 0.0}
    max_blocks = blocks if stop_after_blocks is None else min(blocks, stop_after_blocks)
    local = precompute_local_transitions(profile, b, lam)
    for step in range(max_blocks):
        remaining_blocks = blocks - step
        nxt: dict[tuple[int, int], float] = {}
        for (H, q), base in states.items():
            if H == 0 and q == 0:
                key = (0, 0)
                nxt[key] = log2add(nxt.get(key, float("-inf")), base)
                continue
            max_r = min(b, H)
            for r in range(max_r + 1):
                lp_r = block_input_logprob(b=b, remaining_blocks=remaining_blocks, H=H, r=r)
                if lp_r == float("-inf"):
                    continue
                for _j, q_next, lp_local in local[q][r]:
                    key = (H - r, q_next)
                    val = base + lp_r + lp_local
                    nxt[key] = log2add(nxt.get(key, float("-inf")), val)
        states = nxt
    total = float("-inf")
    for (H, q), val in states.items():
        # If stopped early, pessimistically emit nothing later; this is only
        # useful for debugging.  For full runs H=q=0 terminal mass dominates.
        total = log2add(total, val)
    return lam * d * log2_e + total


def precompute_local_transitions_prob(
    profile: list[dict[int, float]],
    b: int,
    lam: float,
) -> list[list[list[tuple[int, float]]]]:
    local: list[list[list[tuple[int, float]]]] = [[[] for _ in range(b + 1)] for _ in range(b + 1)]
    binoms = [math.comb(b, i) for i in range(b + 1)]
    for q in range(b + 1):
        for r in range(b + 1):
            acc: dict[int, float] = defaultdict(float)
            lo = max(0, r - (b - q))
            hi = min(q, r)
            for ell in range(lo, hi + 1):
                j = q + r - 2 * ell
                ways_overlap = math.comb(q, ell) * math.comb(b - q, r - ell)
                p_overlap = ways_overlap / binoms[r]
                weight = math.exp(-lam * j)
                for q_next, prob in profile[j].items():
                    acc[q_next] += p_overlap * prob * weight
            local[q][r] = list(acc.items())
    return local


def hypergeom_block_probs(*, b: int, remaining_blocks: int, H: int) -> list[float]:
    """Distribution of input occupancy in the next b-bit block.

    This is Hypergeom(total=b*remaining_blocks, good=b, draws=H).  H is small in
    our certificate regime, so a recurrence avoids expensive huge binomials.
    """

    max_r = min(b, H)
    out = [0.0] * (max_r + 1)
    rest = b * (remaining_blocks - 1)
    min_r = max(0, H - rest)
    if min_r > max_r:
        return out
    total = b * remaining_blocks
    logp = (
        math.lgamma(b + 1)
        - math.lgamma(min_r + 1)
        - math.lgamma(b - min_r + 1)
        + math.lgamma(rest + 1)
        - math.lgamma(H - min_r + 1)
        - math.lgamma(rest - H + min_r + 1)
        - math.lgamma(total + 1)
        + math.lgamma(H + 1)
        + math.lgamma(total - H + 1)
    )
    p = math.exp(logp)
    out[min_r] = p
    for r in range(min_r, max_r):
        denom = (r + 1) * (rest - H + r + 1)
        if denom == 0:
            break
        p *= ((b - r) * (H - r)) / denom
        out[r + 1] = p
    # Tiny floating drift is harmless, but normalize when possible.
    s = sum(out)
    if s > 0.0:
        inv = 1.0 / s
        out = [x * inv for x in out]
    return out


def chernoff_for_h_fast(
    *,
    profile: list[dict[int, int]],
    b: int,
    blocks: int,
    h: int,
    d: int,
    lam: float,
    stop_after_blocks: int | None,
) -> float:
    local = precompute_local_transitions_prob(profile, b, lam)
    max_blocks = blocks if stop_after_blocks is None else min(blocks, stop_after_blocks)
    cur = [[0.0 for _ in range(b + 1)] for _ in range(h + 1)]
    cur[h][0] = 1.0
    for step in range(max_blocks):
        remaining_blocks = blocks - step
        nxt = [[0.0 for _ in range(b + 1)] for _ in range(h + 1)]
        for H in range(h + 1):
            row = cur[H]
            if not any(row):
                continue
            p_r = hypergeom_block_probs(b=b, remaining_blocks=remaining_blocks, H=H)
            for q, base in enumerate(row):
                if base == 0.0:
                    continue
                if H == 0 and q == 0:
                    nxt[0][0] += base
                    continue
                for r, pr in enumerate(p_r):
                    if pr == 0.0:
                        continue
                    target_H = H - r
                    scale = base * pr
                    for q_next, p_local in local[q][r]:
                        nxt[target_H][q_next] += scale * p_local
        cur = nxt
    moment = sum(sum(row) for row in cur)
    if moment <= 0.0:
        return float("-inf")
    return (lam * d + math.log(moment)) / math.log(2)


def chernoff_continuation_fast(
    *,
    profile: list[dict[int, int]],
    b: int,
    live_blocks_after_first: int,
    remaining_h: int,
    initial_q: int,
    remaining_d: int,
    lam: float,
) -> float:
    local = precompute_local_transitions_prob(profile, b, lam)
    h = remaining_h
    cur = [[0.0 for _ in range(b + 1)] for _ in range(h + 1)]
    cur[h][initial_q] = 1.0
    for step in range(live_blocks_after_first):
        remaining_blocks = live_blocks_after_first - step
        nxt = [[0.0 for _ in range(b + 1)] for _ in range(h + 1)]
        for H in range(h + 1):
            row = cur[H]
            if not any(row):
                continue
            denom = math.comb(b * remaining_blocks, H)
            max_r = min(b, H)
            p_r = [0.0] * (max_r + 1)
            rest = b * (remaining_blocks - 1)
            for r in range(max_r + 1):
                if H - r <= rest:
                    p_r[r] = math.comb(b, r) * math.comb(rest, H - r) / denom
            for q, base in enumerate(row):
                if base == 0.0:
                    continue
                if H == 0 and q == 0:
                    nxt[0][0] += base
                    continue
                for r, pr in enumerate(p_r):
                    if pr == 0.0:
                        continue
                    target_H = H - r
                    scale = base * pr
                    for q_next, p_local in local[q][r]:
                        nxt[target_H][q_next] += scale * p_local
        cur = nxt
    moment = sum(sum(row) for row in cur)
    if moment <= 0.0:
        return float("-inf")
    return (lam * remaining_d + math.log(moment)) / math.log(2)


def chernoff_continuation_with_local_fast(
    *,
    local: list[list[list[tuple[int, float]]]],
    b: int,
    live_blocks_after_first: int,
    remaining_h: int,
    initial_q: int,
    remaining_d: int,
    lam: float,
) -> float:
    h = remaining_h
    cur = [[0.0 for _ in range(b + 1)] for _ in range(h + 1)]
    cur[h][initial_q] = 1.0
    for step in range(live_blocks_after_first):
        remaining_blocks = live_blocks_after_first - step
        nxt = [[0.0 for _ in range(b + 1)] for _ in range(h + 1)]
        for H in range(h + 1):
            row = cur[H]
            if not any(row):
                continue
            denom = math.comb(b * remaining_blocks, H)
            max_r = min(b, H)
            p_r = [0.0] * (max_r + 1)
            rest = b * (remaining_blocks - 1)
            for r in range(max_r + 1):
                if H - r <= rest:
                    p_r[r] = math.comb(b, r) * math.comb(rest, H - r) / denom
            for q, base in enumerate(row):
                if base == 0.0:
                    continue
                if H == 0 and q == 0:
                    nxt[0][0] += base
                    continue
                for r, pr in enumerate(p_r):
                    if pr == 0.0:
                        continue
                    target_H = H - r
                    scale = base * pr
                    for q_next, p_local in local[q][r]:
                        nxt[target_H][q_next] += scale * p_local
        cur = nxt
    moment = sum(sum(row) for row in cur)
    if moment <= 0.0:
        return float("-inf")
    return (lam * remaining_d + math.log(moment)) / math.log(2)


def chernoff_continuation_sparse_fast(
    *,
    local: list[list[list[tuple[int, float]]]],
    b: int,
    live_blocks_after_first: int,
    remaining_h: int,
    initial_q: int,
    remaining_d: int,
    lam: float,
    prune_below: float = 0.0,
) -> float:
    states: dict[tuple[int, int], float] = {(remaining_h, initial_q): 1.0}
    for step in range(live_blocks_after_first):
        remaining_blocks = live_blocks_after_first - step
        nxt: dict[tuple[int, int], float] = defaultdict(float)
        p_cache: dict[int, list[float]] = {}
        for (H, q), base in states.items():
            if base <= prune_below:
                continue
            if H == 0 and q == 0:
                nxt[(0, 0)] += base
                continue
            p_r = p_cache.get(H)
            if p_r is None:
                p_r = hypergeom_block_probs(b=b, remaining_blocks=remaining_blocks, H=H)
                p_cache[H] = p_r
            for r, pr in enumerate(p_r):
                if pr == 0.0:
                    continue
                scale = base * pr
                target_H = H - r
                for q_next, p_local in local[q][r]:
                    nxt[(target_H, q_next)] += scale * p_local
        states = dict(nxt)
    moment = sum(states.values())
    if moment <= 0.0:
        return float("-inf")
    return (lam * remaining_d + math.log(moment)) / math.log(2)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=7)
    parser.add_argument("--delta-bch", type=int, default=21)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--distance-delta", type=float, default=0.09)
    parser.add_argument("--h-values", default="32")
    parser.add_argument("--lambdas", default="0.001,0.002,0.004,0.008,0.016")
    parser.add_argument("--stop-after-blocks", type=int, default=None)
    parser.add_argument("--exact-j-max", type=int, default=4)
    parser.add_argument("--log-domain", action="store_true")
    parser.add_argument("--profile-mode", choices=("exact-prefix", "random-tail", "entropy-tail", "floor"), default="exact-prefix")
    parser.add_argument("--profile-slack-bits", type=float, default=0.0)
    parser.add_argument("--floor-only-q", type=int, default=32)
    parser.add_argument("--q-bin-size", type=int, default=1)
    parser.add_argument("--continuation", action="store_true")
    parser.add_argument("--live-blocks-after-first", type=int, default=None)
    parser.add_argument("--remaining-h", type=int, default=None)
    parser.add_argument("--initial-q-values", default=None)
    parser.add_argument("--remaining-d", type=int, default=None)
    args = parser.parse_args()

    sys_basis, b, ext_n, row_d0 = build_output_systematic_basis(args.m, args.delta_bch)
    if args.N % b:
        raise ValueError("N must be divisible by block size")
    blocks = args.N // b
    d = math.floor(args.distance_delta * args.N)
    tables = build_parity_tables(sys_basis, b)
    if args.profile_mode == "exact-prefix":
        raw_profile = build_profile(tables, b, args.exact_j_max)
        binoms = [math.comb(b, i) for i in range(b + 1)]
        profile = [
            {q: cnt / binoms[j] for q, cnt in raw_profile[j].items()}
            for j in range(b + 1)
        ]
    else:
        profile = build_envelope_profile(
            b=b,
            d0=row_d0,
            mode=args.profile_mode,
            slack_bits=args.profile_slack_bits,
            floor_only_q=args.floor_only_q,
        )
    profile = quantize_profile(profile, args.q_bin_size)

    print("Permuted block-recursive BCH Chernoff certificate")
    print(
        f"m={args.m}, designed_delta={args.delta_bch}, b={b}, local_length={ext_n}, "
        f"N={args.N}, blocks={blocks}, d={d}, row_d0={row_d0}, "
        f"stop_after_blocks={args.stop_after_blocks}, exact_j_max={args.exact_j_max}, "
        f"profile_mode={args.profile_mode}, profile_slack_bits={args.profile_slack_bits}, "
        f"q_bin_size={args.q_bin_size}"
    )
    print("h,lambda,log2_bound")
    if args.continuation:
        if args.live_blocks_after_first is None or args.remaining_h is None or args.remaining_d is None:
            raise ValueError("--continuation requires --live-blocks-after-first, --remaining-h, and --remaining-d")
        q_values = parse_int_list(args.initial_q_values or "0")
        print("initial_q,lambda,log2_bound")
        for q in q_values:
            for lam in [float(x) for x in args.lambdas.split(",") if x.strip()]:
                bound = chernoff_continuation_fast(
                    profile=profile,
                    b=b,
                    live_blocks_after_first=args.live_blocks_after_first,
                    remaining_h=args.remaining_h,
                    initial_q=q,
                    remaining_d=args.remaining_d,
                    lam=lam,
                )
                print(f"{q},{lam:.8g},{bound:.6f}")
        return 0

    for h in parse_int_list(args.h_values):
        for lam in [float(x) for x in args.lambdas.split(",") if x.strip()]:
            if args.log_domain:
                bound = chernoff_for_h(
                    profile=profile,
                    b=b,
                    blocks=blocks,
                    h=h,
                    d=d,
                    lam=lam,
                    stop_after_blocks=args.stop_after_blocks,
                )
            else:
                bound = chernoff_for_h_fast(
                    profile=profile,
                    b=b,
                    blocks=blocks,
                    h=h,
                    d=d,
                    lam=lam,
                    stop_after_blocks=args.stop_after_blocks,
                )
            print(f"{h},{lam:.8g},{bound:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
