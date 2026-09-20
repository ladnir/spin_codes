#!/usr/bin/env python3
"""Fit the residual on the isolated-run slice after exact late placement."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from check_q1_law import parse_enum
from dense_largek_eval import late_start_base_prob


def exact_slice_prob(rows: list[list[float]], n: int, w: int, cut: int) -> float:
    row = rows[w]
    num = sum((0.0 if v == float("-inf") else 2.0**v) for v in row[: cut + 1])
    return num / math.comb(n, w)


def isolated_full_run_prob(n: int, h: int) -> float:
    return math.comb(n - h + 1, h) / math.comb(n, h)


def main() -> int:
    samples: list[tuple[str, int, float, int, float]] = []
    configs = [
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC12_WC120.120.12.txt"), 120, [10, 11, 12, 13, 14]),
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC16_WC128.128.16.txt"), 128, [12, 13, 14, 15, 16]),
    ]

    for path, n, cuts in configs:
        rows = parse_enum(path)
        for cut in cuts:
            delta = cut / n
            q1 = 2.0 * delta
            l = math.ceil(2.0 * delta * n)
            for h in range(2, 9):
                exact = exact_slice_prob(rows, n, h, cut)
                base = isolated_full_run_prob(n, h) * late_start_base_prob(n, h, h, l)
                bits = math.log2(exact / base)
                bits_per_pair = bits / (h * (h - 1) / 2.0)
                samples.append((path.name, n, q1, h, bits_per_pair))

    xs = [row[2] for row in samples]
    ys = [row[4] for row in samples]
    npts = len(samples)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    slope = (npts * sxy - sx * sy) / (npts * sxx - sx * sx)
    intercept = (sy - slope * sx) / npts
    residuals = [y - (intercept + slope * x) for x, y in zip(xs, ys)]

    # Second fit: allow a mild affine dependence on h itself.
    # bits_per_pair = a + b q1 + c (h-2)
    X = np.array([[1.0, row[2], row[3] - 2] for row in samples], dtype=float)
    y_vec = np.array([row[4] for row in samples], dtype=float)
    coef_h, *_ = np.linalg.lstsq(X, y_vec, rcond=None)
    a, b, c = map(float, coef_h)
    residuals_h = [y - (a + b * x + c * (h - 2)) for _, _, x, h, y in samples]

    out_csv = Path(__file__).with_name("isolated_residual_fit.csv")
    out_txt = Path(__file__).with_name("isolated_residual_fit_report.txt")

    with out_csv.open("w", newline="") as f:
        f.write("enum,n,q1,h,bits_per_pair,pred_bits_per_pair,pred_bits_per_pair_hdep\n")
        for enum_name, nn, q1, h, bits_per_pair in samples:
            pred = intercept + slope * q1
            pred_h = a + b * q1 + c * (h - 2)
            f.write(f"{enum_name},{nn},{q1},{h},{bits_per_pair},{pred},{pred_h}\n")

    with out_txt.open("w") as f:
        f.write("isolated residual fit\n")
        f.write(f"  bits_per_pair = {intercept:.8f} + ({slope:.8f}) * q1\n")
        f.write(f"  max residual  = {max(residuals):.6f} bits/pair\n")
        f.write(f"  min residual  = {min(residuals):.6f} bits/pair\n")
        f.write(f"  h-dependent fit: bits_per_pair = {a:.8f} + ({b:.8f}) * q1 + ({c:.8f}) * (h-2)\n")
        f.write(f"  max residual h = {max(residuals_h):.6f} bits/pair\n")
        f.write(f"  min residual h = {min(residuals_h):.6f} bits/pair\n")
        f.write("  samples:\n")
        for enum_name, nn, q1, h, bits_per_pair in samples:
            pred = intercept + slope * q1
            pred_h = a + b * q1 + c * (h - 2)
            f.write(
                f"    {enum_name}: q1={q1:.6f} h={h} bits/pair={bits_per_pair:.6f} pred={pred:.6f} pred_h={pred_h:.6f}\n"
            )

    print(out_txt.read_text(), end="")
    print(f"\nWrote {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
