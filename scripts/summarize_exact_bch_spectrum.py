#!/usr/bin/env python3
"""Summarize exact small-BCH spectra against the random-like baseline."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from dense_largek_eval import log2_binom


def summarize(path: Path, n: int, k: int) -> dict[str, str]:
    total = 0
    dmin = None
    worst_loss = float("-inf")
    worst_weight = -1
    rows = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            weight = int(row["weight"])
            count = int(row["count"])
            rows.append((weight, count))
            total += count
            if weight > 0 and dmin is None:
                dmin = weight
            if 0 < weight < n:
                randomlike = log2_binom(n, weight) + k - n
                loss = math.log2(count) - randomlike
                if loss > worst_loss:
                    worst_loss = loss
                    worst_weight = weight
    return {
        "path": str(path),
        "n": str(n),
        "k": str(k),
        "total_log2": f"{math.log2(total):.6f}",
        "dmin": str(dmin),
        "worst_weight": str(worst_weight),
        "worst_loss_bits": f"{worst_loss:.6f}",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("items", nargs="+", help="Triples path:n:k")
    args = parser.parse_args()
    summaries = []
    for item in args.items:
        path_s, n_s, k_s = item.split(":")
        summaries.append(summarize(Path(path_s), int(n_s), int(k_s)))
    fields = ["path", "n", "k", "total_log2", "dmin", "worst_weight", "worst_loss_bits"]
    writer = csv.DictWriter(__import__("sys").stdout, fieldnames=fields)
    writer.writeheader()
    writer.writerows(summaries)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
