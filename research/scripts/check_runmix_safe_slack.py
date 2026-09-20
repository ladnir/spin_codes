#!/usr/bin/env python3
"""Measure the slack needed for the safe pairwise run-mixture law.

This script compares exact inner-only slice probabilities against the
theorem-shaped safe pairwise model

    p_h(delta) <= sum_s pi_{N,h}(s) q^s exp(-gamma_safe(q) C(s,2))

on the available exact wrap-convolution tables. It reports the worst observed
log2(exact / model) gap, which is the multiplicative slack (in bits) needed to
turn the model into an empirical upper envelope on the tested region.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Iterable

from dense_largek_eval import run_count_weight, runmix_gamma_safe_from_q


def parse_enum(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        vals: list[float] = []
        for x in line.split(","):
            x = x.strip()
            if not x:
                continue
            vals.append(float("-inf") if x == "-inf" else float(x))
        if vals:
            rows.append(vals)
    return rows


def exact_slice_prob(rows: list[list[float]], n: int, w: int, cut: int) -> float:
    row = rows[w]
    num = sum((0.0 if v == float("-inf") else 2.0**v) for v in row[: cut + 1])
    return num / math.comb(n, w)


def safe_pairwise_prob(n: int, w: int, delta: float) -> float:
    q_single = min(1.0, 2.0 * delta)
    gamma = runmix_gamma_safe_from_q(q_single)
    p = 0.0
    for s in range(1, w + 1):
        p += (
            run_count_weight(n, w, s)
            * (q_single**s)
            * math.exp(-gamma * s * (s - 1) / 2.0)
        )
    return p


def default_tables(root: Path) -> list[Path]:
    # Focus on the larger non-DP wrap-convolution tables; these are the most
    # relevant overlap points for the current theorem-shaped dense model.
    names = [
        "enum_WC12_WC120.120.12.txt",
        "enum_WC16_WC128.128.16.txt",
    ]
    return [root / name for name in names]


def parse_n_from_filename(path: Path) -> int:
    m = re.search(r"\.(\d+)\.\d+\.txt$", path.name)
    if not m:
        raise ValueError(f"Could not parse code length N from filename: {path.name}")
    return int(m.group(1))


def iter_deltas(text: str) -> Iterable[float]:
    for piece in text.split(","):
        piece = piece.strip()
        if piece:
            yield float(piece)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--root",
        type=Path,
        default=Path(r"C:\Users\peter\repo\libOTe"),
        help="Directory containing the exact inner-only enumerator dumps.",
    )
    ap.add_argument(
        "--tables",
        type=str,
        default="",
        help="Comma-separated table filenames to test. Defaults to the large non-DP overlap tables.",
    )
    ap.add_argument(
        "--deltas",
        type=str,
        default="0.10,0.11,0.12,0.125",
        help="Comma-separated delta values to scan.",
    )
    ap.add_argument(
        "--w-max",
        type=int,
        default=12,
        help="Maximum slice weight to check.",
    )
    args = ap.parse_args()

    tables = [args.root / name.strip() for name in args.tables.split(",") if name.strip()]
    if not tables:
        tables = default_tables(args.root)

    deltas = list(iter_deltas(args.deltas))
    global_worst = float("-inf")
    global_info: tuple[str, float, int] | None = None

    for path in tables:
        rows = parse_enum(path)
        n = parse_n_from_filename(path)
        local_worst = float("-inf")
        local_info: tuple[float, int] | None = None
        for delta in deltas:
            cut = math.floor(delta * n)
            for w in range(1, min(args.w_max, n) + 1):
                exact = exact_slice_prob(rows, n, w, cut)
                model = safe_pairwise_prob(n, w, delta)
                if exact <= 0.0 or model <= 0.0:
                    continue
                gap = math.log2(exact / model)
                if gap > local_worst:
                    local_worst = gap
                    local_info = (delta, w)
                if gap > global_worst:
                    global_worst = gap
                    global_info = (path.name, delta, w)

        if local_info is None:
            print(f"{path.name}: no positive points found")
            continue
        delta, w = local_info
        print(
            f"{path.name}: worst log2(exact/safe) = {local_worst:.3f} bits "
            f"at delta={delta:.3f}, w={w}"
        )

    if global_info is not None:
        name, delta, w = global_info
        print()
        print(
            f"GLOBAL worst-case slack = {max(0.0, global_worst):.3f} bits "
            f"on {name} at delta={delta:.3f}, w={w}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
