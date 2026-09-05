#!/usr/bin/env python3
"""Check the conditional first-run late law from Lemma first-run-late-given-runs."""

from __future__ import annotations

import argparse
import math


def exact_tail(n: int, w: int, s: int, t: int) -> float:
    return math.comb(n - t - w + 1, s) / math.comb(n - w + 1, s)


def simple_bound(n: int, t: int, s: int) -> float:
    return (1.0 - t / n) ** s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--w", type=int, required=True)
    ap.add_argument("--s-max", type=int, default=6)
    ap.add_argument("--t-fracs", type=str, default="0.10,0.20,0.24")
    args = ap.parse_args()

    t_fracs = [float(x) for x in args.t_fracs.split(",") if x.strip()]
    for frac in t_fracs:
        t = int(math.floor(frac * args.n))
        print(f"t=floor({frac} * n) = {t}")
        for s in range(1, min(args.s_max, args.w) + 1):
            val = exact_tail(args.n, args.w, s, t)
            bd = simple_bound(args.n, t, s)
            print(
                f"  s={s}: exact={val:.6e}  bound={(bd):.6e}  gap={math.log2(bd/val):.3f} bits"
            )
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
