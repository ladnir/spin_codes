#!/usr/bin/env python3
"""Compare exact inner-only slice probabilities against the exact late-placement base."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

from check_inner_smallw import exact_slice_prob, parse_enum
from dense_largek_eval import late_start_base_prob, run_count_weight


def parse_n_from_filename(path: Path) -> int:
    m = re.search(r"WC\d+_WC(\d+)\.", path.name)
    if not m:
        raise ValueError(f"Could not parse N from {path.name}")
    return int(m.group(1))


def late_start_mixture(n: int, h: int, cut: int) -> float:
    l = math.ceil(2.0 * (cut / n) * n)
    p = 0.0
    for s in range(1, h + 1):
        p += run_count_weight(n, h, s) * late_start_base_prob(n, h, s, l)
    return p


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enum", type=Path, required=True)
    ap.add_argument("--cut", type=int, required=True)
    ap.add_argument("--h-max", type=int, default=10)
    args = ap.parse_args()

    rows = parse_enum(args.enum)
    n = parse_n_from_filename(args.enum)
    print(f"enum   : {args.enum}")
    print(f"n      : {n}")
    print(f"cut    : {args.cut}")
    print(f"delta  : {args.cut / n:.6f}")
    print(f"L      : {math.ceil(2.0 * (args.cut / n) * n)}")
    print()
    print("h, exact, late_base, late/exact, log2(late/exact)")
    for h in range(1, min(args.h_max, len(rows) - 1) + 1):
        exact = exact_slice_prob(rows, n, h, args.cut)
        late = late_start_mixture(n, h, args.cut)
        ratio = late / exact if exact > 0.0 else float("inf")
        log_ratio = math.log2(ratio) if ratio > 0.0 else float("-inf")
        print(f"{h}, {exact:.6e}, {late:.6e}, {ratio:.6e}, {log_ratio:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
