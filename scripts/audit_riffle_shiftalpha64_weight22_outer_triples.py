#!/usr/bin/env python3
"""Audit the optimized outer ratio scanner on a reduced support set."""

from __future__ import annotations

import argparse
import itertools
import json
import struct
import subprocess
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as field  # noqa: E402


DEFAULT_RATIOS = (
    SCRIPT_DIRECTORY.parent
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
    / "goal11_weight22_ratio_counts.bin"
)


def read_ratios(path: Path) -> dict[int, int]:
    data = path.read_bytes()
    if len(data) % 16:
        raise RuntimeError("ratio file contains a partial record")
    result = {}
    for ratio, count in struct.iter_unpack("<QQ", data):
        result[ratio] = count
    if len(result) != len(data) // 16 or sum(result.values()) != 1_365_504:
        raise RuntimeError("ratio file failed mass validation")
    return result


def direct_scan(data_blocks: int, ratios: dict[int, int]) -> dict[str, int]:
    columns = [0] + [
        field.field_power(2, 64 + index) for index in range(data_blocks)
    ]
    result = {
        "data_supports_with_hits": 0,
        "data_outer_words": 0,
        "p0_supports_with_hits": 0,
        "p0_outer_words": 0,
    }
    for support in itertools.combinations(range(len(columns)), 3):
        first, second, third = (columns[index] for index in support)
        first_value = second ^ third
        second_value = first ^ third
        ratio = field.field_multiply(
            second_value, field.field_inverse(first_value)
        )
        count = ratios.get(ratio, 0)
        family = "p0" if support[0] == 0 else "data"
        if count:
            result[f"{family}_supports_with_hits"] += 1
            result[f"{family}_outer_words"] += count
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ratios", type=Path, default=DEFAULT_RATIOS)
    parser.add_argument("--data-blocks", type=int, default=12)
    parser.add_argument(
        "--scanner",
        type=Path,
        default=SCRIPT_DIRECTORY
        / "scan_riffle_shiftalpha64_weight22_outer_triples.exe",
    )
    args = parser.parse_args()
    ratios = read_ratios(args.ratios)
    direct = direct_scan(args.data_blocks, ratios)
    completed = subprocess.run(
        [str(args.scanner), str(args.ratios), str(args.data_blocks)],
        check=True,
        capture_output=True,
        text=True,
    )
    optimized = json.loads(completed.stdout)
    comparisons = {
        "data_supports_with_hits": optimized["data_only"]["supports_with_hits"],
        "data_outer_words": optimized["data_only"][
            "exact_weight22_triple_outer_words"
        ],
        "p0_supports_with_hits": optimized["p0_and_two_data"]["supports_with_hits"],
        "p0_outer_words": optimized["p0_and_two_data"][
            "exact_weight22_triple_outer_words"
        ],
    }
    if direct != comparisons:
        raise SystemExit(f"outer ratio audit mismatch: {direct} != {comparisons}")
    payload = {
        "schema": "riffle-shiftalpha64-weight22-outer-triples-audit-v1",
        "data_blocks": args.data_blocks,
        "direct": direct,
        "optimized": comparisons,
        "ratio_records": len(ratios),
        "ratio_ordered_pair_mass": sum(ratios.values()),
        "status": "PASS",
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
