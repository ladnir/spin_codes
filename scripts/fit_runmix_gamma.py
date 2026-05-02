#!/usr/bin/env python3
"""Fit the pairwise run-mixture damping law gamma(q) from exact inner tables.

This script uses exact inner-only wrapConv tables, inverts the low-weight slice
probabilities into run-conditioned factors q_s, and extracts the effective
pairwise damping parameter from the two-run case

    q_2 ~= q_1^2 exp(-gamma).

It then fits gamma(q_1) with a simple linear law, which is the current default
working model used by the large-k evaluator.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from check_inner_smallw import exact_slice_prob, parse_enum


def infer_qs(rows: list[list[float]], n: int, cut: int, w_max: int) -> list[float]:
    p = [0.0] * (w_max + 1)
    for w in range(1, w_max + 1):
        p[w] = exact_slice_prob(rows, n, w, cut)

    q = [0.0] * (w_max + 1)
    for w in range(1, w_max + 1):
        acc = p[w]
        for s in range(1, w):
            acc -= (
                math.comb(w - 1, s - 1)
                * math.comb(n - w + 1, s)
                / math.comb(n, w)
                * q[s]
            )
        diag = math.comb(n - w + 1, w) / math.comb(n, w)
        q[w] = acc / diag
    return q


def fit_linear(xs: list[float], ys: list[float]) -> tuple[float, float]:
    n = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    det = n * sxx - sx * sx
    if det == 0.0:
        raise ValueError("Degenerate fit")
    a = (sy * sxx - sx * sxy) / det
    b = (n * sxy - sx * sy) / det
    return a, b


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    cases = [
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC12_WC120.120.12.txt"), 120, [10, 11, 12, 13, 14]),
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC16_WC128.128.16.txt"), 128, [12, 13, 14, 15, 16]),
    ]

    rows_out: list[dict[str, float | str]] = []
    xs: list[float] = []
    ys: list[float] = []

    for path, n, cuts in cases:
        rows = parse_enum(path)
        for cut in cuts:
            q = infer_qs(rows, n, cut, 10)
            q1 = q[1]
            q2 = q[2]
            gamma = math.log((q1 * q1) / q2)
            xs.append(q1)
            ys.append(gamma)
            rows_out.append(
                {
                    "file": path.name,
                    "n": n,
                    "cut": cut,
                    "delta": cut / n,
                    "q1": q1,
                    "q2": q2,
                    "gamma_from_q2": gamma,
                }
            )

    a, b = fit_linear(xs, ys)
    rmse = math.sqrt(sum((a + b * x - y) ** 2 for x, y in zip(xs, ys)) / len(xs))

    print(f"fit gamma(q) = a + b q")
    print(f"  a = {a:.8f}")
    print(f"  b = {b:.8f}")
    print(f"  rmse = {rmse:.8f}")
    print()
    print("samples:")
    for row in rows_out:
        pred = a + b * float(row["q1"])
        print(
            f"  {row['file']}: delta={float(row['delta']):.6f} q1={float(row['q1']):.6f} "
            f"gamma={float(row['gamma_from_q2']):.6f} pred={pred:.6f}"
        )

    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["file", "n", "cut", "delta", "q1", "q2", "gamma_from_q2", "pred_gamma"],
            )
            writer.writeheader()
            for row in rows_out:
                row = dict(row)
                row["pred_gamma"] = a + b * float(row["q1"])
                writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
