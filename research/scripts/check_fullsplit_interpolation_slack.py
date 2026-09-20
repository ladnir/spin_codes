#!/usr/bin/env python3
"""Audit interpolation slack for full-split episode inner envelopes.

The envelope summation uses a piecewise-linear interpolation of log2 inner
bounds over H.  For an upper-bound certificate, the interpolated value must be
at least as large as the true/audited value.  This helper compares a coarse knot
CSV against a denser audit CSV and reports the maximum positive lift needed.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_points(path: Path, *, h_column: str, value_column: str) -> list[tuple[int, float]]:
    points: list[tuple[int, float]] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            points.append((int(row[h_column]), float(row[value_column])))
    points.sort()
    return points


def interpolate(points: list[tuple[int, float]], x: int) -> float:
    if x < points[0][0] or x > points[-1][0]:
        raise ValueError(f"H={x} is outside [{points[0][0]}, {points[-1][0]}]")
    if x == points[0][0]:
        return points[0][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x <= x1:
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    return points[-1][1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse-csv", type=Path, required=True)
    parser.add_argument("--audit-csv", type=Path, required=True)
    parser.add_argument("--label", default="")
    parser.add_argument("--h-column", default="H")
    parser.add_argument("--value-column", default="total_log2")
    args = parser.parse_args()

    coarse = load_points(args.coarse_csv, h_column=args.h_column, value_column=args.value_column)
    audit = load_points(args.audit_csv, h_column=args.h_column, value_column=args.value_column)

    worst = None
    max_positive = 0.0
    for H, actual in audit:
        env = interpolate(coarse, H)
        needed = actual - env
        if needed > max_positive:
            max_positive = needed
        if worst is None or needed > worst[2]:
            worst = (H, actual, needed, env)

    print("Full-split interpolation slack audit")
    if args.label:
        print(f"label,{args.label}")
    print(f"coarse_csv,{args.coarse_csv}")
    print(f"audit_csv,{args.audit_csv}")
    print(f"audited_points,{len(audit)}")
    print(f"max_positive_lift_bits,{max_positive:.6f}")
    if worst is not None:
        H, actual, needed, env = worst
        print(f"worst_H,{H}")
        print(f"worst_actual_log2,{actual:.6f}")
        print(f"worst_linear_log2,{env:.6f}")
        print(f"worst_needed_lift_bits,{needed:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
