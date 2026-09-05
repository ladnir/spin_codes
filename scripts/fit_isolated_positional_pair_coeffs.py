#!/usr/bin/env python3
"""Fit q1,m coefficient laws for the positional early-pair residual."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from check_q1_law import parse_enum
from check_isolated_positional import exact_slice_prob
from dense_largek_eval import isolated_early_pair_mean, isolated_position_single_model_prob


def main() -> int:
    configs = [
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC12_WC120.120.12.txt"), 120, 12, [10, 11, 12, 13, 14]),
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC16_WC128.128.16.txt"), 128, 16, [12, 13, 14, 15, 16]),
    ]

    rows_out: list[tuple[float, float, float, float]] = []  # q1, m, a, b

    for path, n, m, cuts in configs:
        rows = parse_enum(path)
        for cut in cuts:
            q1 = 2.0 * cut / n
            X: list[list[float]] = []
            y: list[float] = []
            l = 2 * cut
            for h in range(4, 11):
                exact = exact_slice_prob(rows, n, h, cut)
                model = isolated_position_single_model_prob(n, h, cut, m)
                gap_bits = math.log2(exact / model)
                epair = isolated_early_pair_mean(n, h, l)
                X.append([1.0, epair])
                y.append(gap_bits)
            X_np = np.array(X, dtype=float)
            y_np = np.array(y, dtype=float)
            coef, *_ = np.linalg.lstsq(X_np, y_np, rcond=None)
            a, b = map(float, coef)
            rows_out.append((q1, float(m), a, b))

    q = np.array([r[0] for r in rows_out], dtype=float)
    m = np.array([r[1] for r in rows_out], dtype=float)
    a = np.array([r[2] for r in rows_out], dtype=float)
    b = np.array([r[3] for r in rows_out], dtype=float)

    Xa = np.column_stack([np.ones_like(q), q, 1.0 / m])
    coef_a, *_ = np.linalg.lstsq(Xa, a, rcond=None)
    pred_a = Xa @ coef_a
    resid_a = a - pred_a

    Xb = np.column_stack([np.ones_like(q), q, 1.0 / m])
    coef_b, *_ = np.linalg.lstsq(Xb, b, rcond=None)
    pred_b = Xb @ coef_b
    resid_b = b - pred_b

    out = Path(__file__).with_name("isolated_positional_pair_coeffs_report.txt")
    with out.open("w") as f:
        f.write("isolated positional pair coefficient fit\n")
        f.write(
            "  a(q1,m) ~= "
            f"{coef_a[0]:.8f} + ({coef_a[1]:.8f}) q1 + ({coef_a[2]:.8f}) / m\n"
        )
        f.write(f"  a residual max = {max(resid_a):.6f}\n")
        f.write(f"  a residual min = {min(resid_a):.6f}\n")
        f.write(
            "  b(q1,m) ~= "
            f"{coef_b[0]:.8f} + ({coef_b[1]:.8f}) q1 + ({coef_b[2]:.8f}) / m\n"
        )
        f.write(f"  b residual max = {max(resid_b):.6f}\n")
        f.write(f"  b residual min = {min(resid_b):.6f}\n")
        f.write("  local fits:\n")
        for q1, mm, aa, bb in rows_out:
            f.write(f"    q1={q1:.6f} m={int(mm)} a={aa:.6f} b={bb:.6f}\n")

    print(out.read_text(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
