#!/usr/bin/env python3
"""Freeze a grid atlas of generalized packet-profile outer witnesses."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from packet_group_outer_profile import atom_count
from packet_group_profile_bound import outer_probability


def integral_profile(group_bits: int, counts: tuple[int, ...], denominator: int) -> list[int]:
    atoms = atom_count(group_bits)
    raw = np.asarray(counts, dtype=np.float64) * atoms / denominator
    profile = np.floor(raw).astype(np.int64)
    remainder = atoms - int(np.sum(profile))
    order = np.argsort(-(raw - profile))
    profile[order[:remainder]] += 1
    return profile.tolist()


def compositions(total: int, parts: int, prefix: tuple[int, ...] = ()):
    if parts == 1:
        yield prefix + (total,)
        return
    for first in range(total + 1):
        yield from compositions(total - first, parts - 1, prefix + (first,))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--denominator", type=int, default=16)
    parser.add_argument("--full-outer", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.denominator <= 0:
        raise SystemExit("outer atlas: denominator must be positive")

    rows = []
    grid = list(compositions(args.denominator, args.group + 1))
    for index, counts in enumerate(grid, 1):
        profile = integral_profile(args.group, counts, args.denominator)
        value, details = outer_probability(
            args.group, profile, fast=not args.full_outer
        )
        log_variables = np.asarray(details["log_variables"], dtype=np.float64)
        charge = log_variables / math.log(2.0)
        constant = value + float(np.asarray(profile) @ charge)
        row = {
            "name": "outer_" + "_".join(map(str, counts)),
            "grid_counts": list(counts),
            "grid_denominator": args.denominator,
            "profile": profile,
            "constant_log2": constant,
            "charge": charge.tolist(),
            "anchor_value_log2": value,
            "outer_details": details,
        }
        rows.append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(
                {
                    "status": "DIAGNOSTIC_BINARY64_GENERALIZED_OUTER_FIXED_ATLAS",
                    "group_bits": args.group,
                    "denominator": args.denominator,
                    "rows": rows,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        if index % 20 == 0 or index == len(grid):
            print(
                f"outer_witnesses={index}/{len(grid)} last={row['name']} "
                f"anchor={value:.9f}",
                flush=True,
            )
    print("status=DIAGNOSTIC_BINARY64_GENERALIZED_OUTER_FIXED_ATLAS")


if __name__ == "__main__":
    main()
