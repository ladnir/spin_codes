#!/usr/bin/env python3
"""Convert a dual-weight window into a BCH local-spectrum envelope.

For a binary [n,k] code C with redundancy r=n-k, MacWilliams gives

    A_h(C) = 2^{-r} sum_j A_j(C^perp) K_h(j).

If all nonzero dual words lie in a window W, then

    A_h(C) <= 2^{-r} binom(n,h) + max_{j in W} |K_h(j)|.

Equivalently, relative to the random-like main term,

    A_h(C) / (2^{-r} binom(n,h))
        <= 1 + 2^r max_{j in W} |K_h(j)| / binom(n,h).

This script evaluates that finite implication.  It is intended to audit
whether classical BCH dual-width bounds, such as Carlitz--Uchiyama or its
improvements, are strong enough for the block-outer first-moment certificate.
"""

from __future__ import annotations

import argparse
import math


def norm_krawtchouk(n: int, h: int, x: int) -> float:
    """Return K_h(x) / binom(n,h) by a stable normalized recurrence."""
    if h == 0:
        return 1.0
    l_prev = 1.0
    l_cur = (n - 2 * x) / n
    if h == 1:
        return l_cur
    for i in range(1, h):
        l_next = ((n - 2 * x) * l_cur - i * l_prev) / (n - i)
        l_prev, l_cur = l_cur, l_next
    return l_cur


def max_window_norm_krawtchouk(n: int, h: int, lo: int, hi: int) -> tuple[float, int]:
    best = 0.0
    best_j = lo
    for j in range(lo, hi + 1):
        val = abs(norm_krawtchouk(n, h, j))
        if val > best:
            best = val
            best_j = j
    return best, best_j


def log2_one_plus_pow(log_x: float) -> float:
    if log_x < -40.0:
        return math.log1p(2.0**log_x) / math.log(2.0)
    return math.log2(1.0 + 2.0**log_x)


def parse_weights(text: str) -> list[int]:
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=511)
    parser.add_argument("--k", type=int, default=250)
    parser.add_argument("--designed-distance", type=int, default=63)
    parser.add_argument("--weights", default="63,70,89,100,150,220")
    parser.add_argument(
        "--window-half-width",
        type=float,
        default=None,
        help="Dual half-window around (n+1)/2. Default: Carlitz-Uchiyama (t-1)*sqrt(n+1).",
    )
    parser.add_argument(
        "--target-slack",
        type=float,
        default=62.71,
        help="Reference uniform spectrum slack budget, in bits.",
    )
    args = parser.parse_args()

    t = (args.designed_distance - 1) // 2
    center = (args.n + 1) / 2.0
    half = args.window_half_width
    if half is None:
        half = (t - 1) * math.sqrt(args.n + 1)
    lo = max(1, math.ceil(center - half))
    hi = min(args.n, math.floor(center + half))
    r = args.n - args.k
    weights = parse_weights(args.weights)

    print("Dual-window Krawtchouk envelope")
    print(f"n={args.n}, k={args.k}, r={r}, designed_distance={args.designed_distance}, t={t}")
    print(f"center={(args.n + 1) / 2:.3f}, half_width={half:.6f}, window=[{lo},{hi}]")
    print(f"target_slack_bits={args.target_slack:.3f}")
    print("h,max_j,log2_max_abs_K_over_binom,log2_ratio_correction,log2_total_envelope_factor,slack_gap")
    for h in weights:
        max_norm, max_j = max_window_norm_krawtchouk(args.n, h, lo, hi)
        log_norm = math.log2(max_norm) if max_norm > 0.0 else float("-inf")
        ratio_corr = r + log_norm if log_norm != float("-inf") else float("-inf")
        total_factor = log2_one_plus_pow(ratio_corr) if ratio_corr != float("-inf") else 0.0
        print(
            f"{h},{max_j},{log_norm:.6f},{ratio_corr:.6f},{total_factor:.6f},"
            f"{args.target_slack - total_factor:.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
