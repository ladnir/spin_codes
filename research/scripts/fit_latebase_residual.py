#!/usr/bin/env python3
"""Fit the residual correction after factoring out the exact late-placement base."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from check_inner_smallw import parse_enum
from dense_largek_eval import late_start_base_prob
from fit_runmix_gamma import infer_qs


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
    ys_nat: list[float] = []

    for path, n, cuts in cases:
        rows = parse_enum(path)
        for cut in cuts:
            q = infer_qs(rows, n, cut, 8)
            l = math.ceil(2.0 * (cut / n) * n)
            q1 = q[1]
            for s in range(2, 7):
                base = late_start_base_prob(n, s, s, l)
                corr = q[s] / base
                pairs = s * (s - 1) / 2.0
                beta_nat = math.log(corr) / pairs
                xs.append(q1)
                ys_nat.append(beta_nat)
                rows_out.append(
                    {
                        "file": path.name,
                        "n": n,
                        "cut": cut,
                        "delta": cut / n,
                        "q1": q1,
                        "s": s,
                        "pairs": pairs,
                        "q_s": q[s],
                        "late_base": base,
                        "corr": corr,
                        "beta_nat": beta_nat,
                        "beta_bits": beta_nat / math.log(2.0),
                    }
                )

    a, b = fit_linear(xs, ys_nat)
    rmse = math.sqrt(sum((a + b * x - y) ** 2 for x, y in zip(xs, ys_nat)) / len(xs))

    print("fit beta(q1) = a + b q1   [natural-log units per pair]")
    print(f"  a = {a:.8f}")
    print(f"  b = {b:.8f}")
    print(f"  rmse = {rmse:.8f}")
    print(f"  beta_bits(q1=0.24) ~= {(a + b * 0.24) / math.log(2.0):.6f}")
    print()
    print("samples:")
    for row in rows_out:
        pred = a + b * float(row["q1"])
        print(
            f"  {row['file']}: delta={float(row['delta']):.6f} s={int(row['s'])} "
            f"beta_bits={float(row['beta_bits']):.6f} pred_bits={pred / math.log(2.0):.6f}"
        )

    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "file",
                    "n",
                    "cut",
                    "delta",
                    "q1",
                    "s",
                    "pairs",
                    "q_s",
                    "late_base",
                    "corr",
                    "beta_nat",
                    "beta_bits",
                    "pred_beta_nat",
                    "pred_beta_bits",
                ],
            )
            writer.writeheader()
            for row in rows_out:
                row = dict(row)
                pred = a + b * float(row["q1"])
                row["pred_beta_nat"] = pred
                row["pred_beta_bits"] = pred / math.log(2.0)
                writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
