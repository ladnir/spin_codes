#!/usr/bin/env python3
"""Fit the residual after the position-sensitive isolated one-run model."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from check_q1_law import parse_enum
from check_isolated_positional import exact_slice_prob
from dense_largek_eval import isolated_early_pair_mean, isolated_position_single_model_prob


def main() -> int:
    samples: list[tuple[str, int, int, float]] = []
    configs = [
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC12_WC120.120.12.txt"), 120, 0.115217, 12),
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC16_WC128.128.16.txt"), 128, 0.125, 16),
    ]

    for path, n, delta, m in configs:
        rows = parse_enum(path)
        cut = math.floor(delta * n)
        l = 2 * cut
        for h in range(4, 11):
            exact = exact_slice_prob(rows, n, h, cut)
            model = isolated_position_single_model_prob(n, h, cut, m)
            gap_bits = math.log2(exact / model)
            epair = isolated_early_pair_mean(n, h, l)
            samples.append((path.name, h, epair, gap_bits))

    X = np.array([[1.0, row[2]] for row in samples], dtype=float)
    y = np.array([row[3] for row in samples], dtype=float)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b = map(float, coef)
    residuals = [gap - (a + b * ep) for _, _, ep, gap in samples]

    safe_b = 0.186
    safe_a = max(gap - safe_b * ep for _, _, ep, gap in samples)
    safe_residuals = [safe_a + safe_b * ep - gap for _, _, ep, gap in samples]

    out = Path(__file__).with_name("isolated_positional_pair_report.txt")
    with out.open("w") as f:
        f.write("isolated positional residual fit\n")
        f.write(f"  fit: gap_bits = {a:.8f} + ({b:.8f}) * E[pairs]\\n")
        f.write(f"  max residual = {max(residuals):.6f}\\n")
        f.write(f"  min residual = {min(residuals):.6f}\\n")
        f.write(f"  safe envelope: gap_bits <= {safe_a:.8f} + ({safe_b:.8f}) * E[pairs]\\n")
        f.write(f"  max safe slack = {max(safe_residuals):.6f}\\n")
        f.write("  samples:\\n")
        for name, h, ep, gap in samples:
            pred = a + b * ep
            safe = safe_a + safe_b * ep
            f.write(
                f"    {name}: h={h} Epair={ep:.6f} gap={gap:.6f} pred={pred:.6f} safe={safe:.6f}\\n"
            )

    print(out.read_text(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
