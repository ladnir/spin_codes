#!/usr/bin/env python3
"""Check the one-run law p_1(delta) against the 2*delta anchor on exact tables."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


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


def parse_n(path: Path) -> int:
    m = re.search(r"\.(\d+)\.\d+\.txt$", path.name)
    if not m:
        raise ValueError(f"Could not parse N from filename: {path.name}")
    return int(m.group(1))


def exact_p1(rows: list[list[float]], n: int, cut: int) -> float:
    row = rows[1]
    num = sum((0.0 if v == float("-inf") else 2.0**v) for v in row[: cut + 1])
    return num / n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--root",
        type=Path,
        default=Path(r"C:\Users\peter\repo\libOTe"),
        help="Directory containing exact inner-only enumerator dumps.",
    )
    ap.add_argument(
        "--tables",
        type=str,
        default="enum_WC12_WC120.120.12.txt,enum_WC16_WC128.128.16.txt",
        help="Comma-separated exact-table filenames to test.",
    )
    ap.add_argument(
        "--deltas",
        type=str,
        default="0.10,0.11,0.12,0.125",
        help="Comma-separated delta values to scan.",
    )
    args = ap.parse_args()

    deltas = [float(x) for x in args.deltas.split(",") if x.strip()]
    tables = [args.root / x.strip() for x in args.tables.split(",") if x.strip()]

    for path in tables:
        rows = parse_enum(path)
        n = parse_n(path)
        print(path.name)
        for delta in deltas:
            cut = math.floor(delta * n)
            p1 = exact_p1(rows, n, cut)
            anchor = (2.0 * cut) / n
            print(
                f"  delta={delta:.3f} cut={cut:2d} "
                f"p1={p1:.6f} 2cut/n={anchor:.6f} diff={p1-anchor:+.6e}"
            )
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
