#!/usr/bin/env python3
"""Audit selected-gap e-tail ratios in full-split episode CSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_log2(text: str) -> float:
    text = text.strip()
    if text == "-inf":
        return float("-inf")
    return float(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--start-e", type=int, default=1)
    parser.add_argument("--stop-e", type=int)
    parser.add_argument("--label", default="")
    args = parser.parse_args()

    worst = None
    count = 0
    with args.csv.open(newline="") as f:
        reader = csv.DictReader(f)
        fields = set(reader.fieldnames or [])
        stop_e = args.stop_e
        if stop_e is None:
            e_cols = [
                int(name[1:-5])
                for name in fields
                if name.startswith("e") and name.endswith("_log2") and name[1:-5].isdigit()
            ]
            stop_e = max(e_cols) - 1
        for row in reader:
            H = int(row["H"])
            for e in range(args.start_e, stop_e + 1):
                a_col = f"e{e}_log2"
                b_col = f"e{e + 1}_log2"
                if a_col not in fields or b_col not in fields:
                    continue
                a = parse_log2(row[a_col])
                b = parse_log2(row[b_col])
                if a == float("-inf") or b == float("-inf"):
                    continue
                diff = b - a
                count += 1
                if worst is None or diff > worst["diff"]:
                    worst = {"H": H, "e": e, "diff": diff, "prev": a, "next": b}

    print("Full-split e-tail ratio audit")
    if args.label:
        print(f"label,{args.label}")
    print(f"csv,{args.csv}")
    print(f"start_e,{args.start_e}")
    print(f"stop_e,{stop_e}")
    print(f"checked_pairs,{count}")
    if worst is not None:
        print(f"worst_H,{worst['H']}")
        print(f"worst_e,{worst['e']}")
        print(f"worst_next_minus_prev_log2,{worst['diff']:.6f}")
        print(f"worst_prev_log2,{worst['prev']:.6f}")
        print(f"worst_next_log2,{worst['next']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
