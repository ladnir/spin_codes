"""Burst-restricted fixed-tap banded outer generating functions.

For an allowed support set S consisting of one burst interval plus the final
late tail, this computes

    E[ sum_x z^{wt(outer(x))} 1{supp(outer(x)) subset S} ]

for the active fixed-tap systematic banded outer.  The calculation is exact for
a fixed interval and z, using a transfer matrix over the distance since the last
message one.

The current scan over burst starts is a diagnostic/max-sampled layer.  It is the
right object for the next proof step, but a formal theorem still needs a clean
argument identifying or upper-bounding the worst start class.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from dense_largek_eval import binom_cdf_half_entropy_log2, log2_binom, log2add


def format_log2(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.6f}"


def z_grid(z_min: float, z_max: float, count: int) -> list[float]:
    if count <= 1:
        return [z_max]
    lo = math.log(z_min)
    hi = math.log(z_max)
    return [math.exp(lo + (hi - lo) * i / (count - 1)) for i in range(count)]


def normalize_vec(v: np.ndarray) -> tuple[np.ndarray, float]:
    scale = float(np.max(v))
    if scale <= 0.0:
        return v, float("-inf")
    return v / scale, math.log2(scale)


def normalize_mat(mtx: np.ndarray) -> tuple[np.ndarray, float]:
    scale = float(np.max(mtx))
    if scale <= 0.0:
        return mtx, float("-inf")
    return mtx / scale, math.log2(scale)


def apply_power(v: np.ndarray, v_log: float, matrix: np.ndarray, exp: int) -> tuple[np.ndarray, float]:
    base = matrix
    base_log = 0.0
    while exp > 0:
        if exp & 1:
            v = base @ v
            v, inc = normalize_vec(v)
            v_log += base_log + inc
        exp >>= 1
        if exp:
            base = base @ base
            base, inc = normalize_mat(base)
            base_log = 2.0 * base_log + inc
    return v, v_log


def step_matrix(*, sigma: int, z: float, sys_allowed: bool, parity_allowed: bool) -> np.ndarray:
    states = sigma + 1
    inactive = sigma
    a_allowed = (1.0 + z) / 2.0
    a_forced_zero = 0.5
    mtx = np.zeros((states, states), dtype=np.float64)
    for d in range(states):
        prev_active = d < sigma

        # x_t = 1. The systematic output must be allowed.
        if sys_allowed:
            sys_factor = z
            if prev_active:
                par_factor = a_allowed if parity_allowed else a_forced_zero
            else:
                # First active parity coordinate in a fixed-tap cluster is a
                # deterministic one.
                par_factor = z if parity_allowed else 0.0
            if par_factor:
                mtx[0, d] += sys_factor * par_factor

        # x_t = 0.
        nd = min(sigma, d + 1)
        if nd < sigma:
            par_factor = a_allowed if parity_allowed else a_forced_zero
        else:
            par_factor = 1.0
        mtx[nd, d] += par_factor
    return mtx


def msg_only_matrix(*, sigma: int, z: float, sys_allowed: bool) -> np.ndarray:
    states = sigma + 1
    mtx = np.zeros((states, states), dtype=np.float64)
    for d in range(states):
        if sys_allowed:
            mtx[0, d] += z
        mtx[min(sigma, d + 1), d] += 1.0
    return mtx


def parity_only_matrix(*, sigma: int, z: float, parity_allowed: bool) -> np.ndarray:
    states = sigma + 1
    a_allowed = (1.0 + z) / 2.0
    a_forced_zero = 0.5
    mtx = np.zeros((states, states), dtype=np.float64)
    for d in range(states):
        nd = min(sigma, d + 1)
        if nd < sigma:
            par_factor = a_allowed if parity_allowed else a_forced_zero
        else:
            par_factor = 1.0
        mtx[nd, d] += par_factor
    return mtx


def interval_allowed(pos: int, intervals: list[tuple[int, int]]) -> bool:
    return any(lo <= pos < hi for lo, hi in intervals)


def restricted_outer_gf_log2(
    *,
    k: int,
    parity_n: int,
    sigma: int,
    z: float,
    burst_start: int,
    burst_len: int,
    tail_len: int,
) -> float:
    total_n = k + parity_n
    burst = (burst_start, min(total_n, burst_start + burst_len))
    tail = (max(0, total_n - tail_len), total_n)
    intervals = [burst, tail]
    states = sigma + 1
    inactive = sigma
    v = np.zeros(states, dtype=np.float64)
    v[inactive] = 1.0
    v_log = 0.0

    max_t = max(k, parity_n)
    breaks = {0, max_t}
    for lo, hi in intervals:
        breaks.add(min(max_t, max(0, lo)))
        breaks.add(min(max_t, max(0, hi)))
        breaks.add(min(max_t, max(0, lo - k)))
        breaks.add(min(max_t, max(0, hi - k)))
    breaks.add(k)
    breaks.add(parity_n)
    sorted_breaks = sorted(b for b in breaks if 0 <= b <= max_t)
    for t, end in zip(sorted_breaks, sorted_breaks[1:]):
        if end <= t:
            continue
        sys_allowed = t < k and interval_allowed(t, intervals)
        parity_allowed = t < parity_n and interval_allowed(k + t, intervals)

        if t < k and t < parity_n:
            matrix = step_matrix(sigma=sigma, z=z, sys_allowed=sys_allowed, parity_allowed=parity_allowed)
        elif t < k:
            matrix = msg_only_matrix(sigma=sigma, z=z, sys_allowed=sys_allowed)
        else:
            matrix = parity_only_matrix(sigma=sigma, z=z, parity_allowed=parity_allowed)
        v, v_log = apply_power(v, v_log, matrix, end - t)

    v_sum = float(np.sum(v))
    if v_sum <= 0.0:
        return float("-inf")
    log_all = v_log + math.log2(v_sum)
    # Zero message is always allowed and contributes one.
    if log_all <= 1e-10:
        val = max(0.0, (2.0**log_all) - 1.0)
        return math.log2(val) if val > 0.0 else float("-inf")
    return log_all + math.log2(1.0 - 2.0 ** (-log_all))


def burst_starts(total_n: int, burst_len: int, tail_len: int, samples: int) -> list[int]:
    hi = max(0, total_n - tail_len - burst_len)
    anchors = {0, hi, total_n // 2 - burst_len // 2}
    # Around the systematic/parity boundary.
    k = total_n // 2
    for s in (k - burst_len, k - burst_len // 2, k - 1, k, k + 1):
        anchors.add(min(hi, max(0, s)))
    if samples > 1:
        for i in range(samples):
            anchors.add(round(hi * i / (samples - 1)))
    return sorted(s for s in anchors if 0 <= s <= hi)


def separated_suffix_lower_log2(*, h: int, sigma: int, burst_sys_tail_len: int) -> tuple[float, int]:
    """Lower bound from separated message ones in the systematic burst suffix.

    Split the usable suffix into sigma-spaced slots and select w slots. Each
    selected message one contributes one systematic one and one deterministic
    fixed-tap parity one; the remaining w*(sigma-1) active parity coordinates
    are fair. This gives a concrete subfamily contained in the burst+late-tail
    support when the parity windows of the suffix lie inside the late tail.
    """
    slots = burst_sys_tail_len // sigma
    best = float("-inf")
    best_w = 0
    for w in range(1, min(h, slots) + 1):
        fair = w * (sigma - 1)
        need = h - 2 * w
        if 0 <= need <= fair:
            val = log2_binom(slots, w) + log2_binom(fair, need) - fair
            if val > best:
                best = val
                best_w = w
    return best, best_w


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--sigma", type=int, required=True)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--burst-len", type=int, default=None)
    parser.add_argument("--tail-len", type=int, default=None)
    parser.add_argument("--small-xi", type=float, default=0.1)
    parser.add_argument("--h-max", type=int, default=512)
    parser.add_argument("--z-min", type=float, default=0.02)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=128)
    parser.add_argument("--start-samples", type=int, default=65)
    parser.add_argument("--out-prefix", default=None)
    parser.add_argument("--no-png", action="store_true")
    args = parser.parse_args()

    parity_n = args.k
    total_n = args.k + parity_n
    d = math.floor(args.delta * total_n)
    burst_len = args.burst_len if args.burst_len is not None else 2 * d
    tail_len = args.tail_len if args.tail_len is not None else math.ceil((2.0 + args.small_xi) * args.delta * total_n)
    starts = burst_starts(total_n, burst_len, tail_len, args.start_samples)
    num_possible_starts = max(1, total_n - tail_len - burst_len + 1)
    inner_cost = -(args.sigma - 1) + binom_cdf_half_entropy_log2(max(0, burst_len - args.sigma), d)
    tail_parity_start = max(0, total_n - tail_len - args.k)
    zs = z_grid(args.z_min, args.z_max, args.z_count)

    print(
        f"k={args.k} sigma={args.sigma} total_n={total_n} burst_len={burst_len} "
        f"tail_len={tail_len} sampled_starts={len(starts)} possible_starts={num_possible_starts}"
    )

    rows: list[dict[str, str | int | float]] = []
    best = [float("-inf")] * (args.h_max + 1)
    best_start = [-1] * (args.h_max + 1)
    best_z = [float("nan")] * (args.h_max + 1)
    for s in starts:
        start_best = [float("inf")] * (args.h_max + 1)
        start_z = [float("nan")] * (args.h_max + 1)
        for z in zs:
            log_w = restricted_outer_gf_log2(
                k=args.k,
                parity_n=parity_n,
                sigma=args.sigma,
                z=z,
                burst_start=s,
                burst_len=burst_len,
                tail_len=tail_len,
            )
            log_z = math.log2(z)
            for h in range(1, args.h_max + 1):
                val = log_w - h * log_z
                if val < start_best[h]:
                    start_best[h] = val
                    start_z[h] = z
        for h in range(1, args.h_max + 1):
            val = start_best[h]
            if val != float("inf") and (best[h] == float("-inf") or val > best[h]):
                best[h] = val
                best_start[h] = s
                best_z[h] = start_z[h]

    union_log = math.log2(num_possible_starts)
    for h in range(1, args.h_max + 1):
        burst = (best_start[h], min(total_n, best_start[h] + burst_len))
        usable_sys_lo = max(0, burst[0], tail_parity_start)
        usable_sys_hi = min(args.k, burst[1])
        usable_sys_len = max(0, usable_sys_hi - usable_sys_lo)
        lower, lower_w = separated_suffix_lower_log2(
            h=h,
            sigma=args.sigma,
            burst_sys_tail_len=usable_sys_len,
        )
        rows.append(
            {
                "h": h,
                "sampled_max_coeff_log2": format_log2(best[h]),
                "union_start_bound_log2": format_log2(best[h] + union_log),
                "inner_one_burst_cost_log2": format_log2(inner_cost),
                "union_outer_inner_log2": format_log2(best[h] + union_log + inner_cost),
                "separated_suffix_lower_log2": format_log2(lower),
                "separated_suffix_lower_plus_inner_log2": format_log2(lower + inner_cost),
                "separated_suffix_best_w": lower_w,
                "usable_sys_tail_len": usable_sys_len,
                "best_start": best_start[h],
                "best_z": f"{best_z[h]:.12g}",
            }
        )
        if h == 1 or h % 32 == 0 or h == args.h_max:
            print(
                f"h={h:4d} sampled={format_log2(best[h])} "
                f"union+inner={format_log2(best[h] + union_log + inner_cost)} "
                f"lower+inner={format_log2(lower + inner_cost)} start={best_start[h]}"
            )

    if args.out_prefix is None:
        out_base = Path(__file__).resolve().parent / (
            f"burst_restricted_outer_k{args.k}_sig{args.sigma}_h{args.h_max}"
        )
    else:
        out_base = Path(args.out_prefix)
        if not out_base.is_absolute():
            out_base = Path.cwd() / out_base
    out_base.parent.mkdir(parents=True, exist_ok=True)

    csv_path = out_base.with_suffix(".csv")
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "h",
                "sampled_max_coeff_log2",
                "union_start_bound_log2",
                "inner_one_burst_cost_log2",
                "union_outer_inner_log2",
                "separated_suffix_lower_log2",
                "separated_suffix_lower_plus_inner_log2",
                "separated_suffix_best_w",
                "usable_sys_tail_len",
                "best_start",
                "best_z",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {csv_path}")

    if args.no_png:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not available; skipped PNG")
        return
    hs = [int(r["h"]) for r in rows]
    sampled = [float(r["sampled_max_coeff_log2"]) for r in rows]
    unioned = [float(r["union_start_bound_log2"]) for r in rows]
    product = [float(r["union_outer_inner_log2"]) for r in rows]
    lower = [
        float(r["separated_suffix_lower_plus_inner_log2"])
        if r["separated_suffix_lower_plus_inner_log2"] != "-inf"
        else float("nan")
        for r in rows
    ]
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    ax.plot(hs, sampled, label="sampled max interval")
    ax.plot(hs, unioned, label="sampled max + log2 possible starts")
    ax.plot(hs, product, label="union outer + one-burst inner cost")
    ax.plot(hs, lower, label="separated suffix lower + inner", linestyle="--")
    ax.set_xlabel("outer weight h")
    ax.set_ylabel("log2 burst-restricted outer coefficient bound")
    ax.set_title(f"Burst-restricted fixed-tap banded outer, sigma={args.sigma}")
    ax.legend()
    png_path = out_base.with_suffix(".png")
    fig.savefig(png_path, dpi=180)
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
