#!/usr/bin/env python3
"""Import a simple .wd weight-distribution file into weight,count CSV format.

The Okayama weight-distribution tables are plain text with an initial comment
followed by alternating weight/count tokens, for example:

    # EBCH128_64.wd desaki 0 1 22 243840 ...

This converts that format to the spectrum CSV consumed by
block_outer_upgrade_probe.py.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_wd(text: str) -> list[tuple[int, int]]:
    tokens = text.split()
    while tokens and tokens[0].startswith("#"):
        tokens.pop(0)
        while tokens and not tokens[0].lstrip("-").isdigit():
            tokens.pop(0)
    if len(tokens) % 2 != 0:
        raise ValueError("expected an even number of weight/count tokens")
    rows: list[tuple[int, int]] = []
    for i in range(0, len(tokens), 2):
        rows.append((int(tokens[i]), int(tokens[i + 1])))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--expect-dim", type=int, default=None)
    parser.add_argument("--expect-length", type=int, default=None)
    args = parser.parse_args()

    rows = parse_wd(args.input.read_text())
    total = sum(count for _, count in rows)
    if args.expect_dim is not None and total != 1 << args.expect_dim:
        raise ValueError(f"counts sum to {total}, expected 2^{args.expect_dim}")
    max_weight = max(weight for weight, _ in rows)
    if args.expect_length is not None and max_weight != args.expect_length:
        raise ValueError(f"max weight is {max_weight}, expected {args.expect_length}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["weight", "count"])
        writer.writerows(rows)

    dmin = min(weight for weight, count in rows if weight > 0 and count)
    print(f"wrote {args.output}")
    print(f"weights={len(rows)}")
    print(f"sum={total}")
    print(f"log2_sum={total.bit_length() - 1 if total and total & (total - 1) == 0 else 'not_power_of_two'}")
    print(f"max_weight={max_weight}")
    print(f"min_nonzero_weight={dmin}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
