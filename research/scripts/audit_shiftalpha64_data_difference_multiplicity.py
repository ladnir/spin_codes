#!/usr/bin/env python3
"""Audit the optimized ShiftAlpha64 difference scan on a reduced instance."""

from __future__ import annotations

import argparse
import collections
import json
import subprocess
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as field  # noqa: E402


def direct(data_blocks: int) -> dict[str, int]:
    coefficients = [
        field.field_power(2, 64 + index) for index in range(data_blocks)
    ]
    counts = collections.Counter(
        coefficients[first] ^ coefficients[second]
        for second in range(1, data_blocks)
        for first in range(second)
    )
    return {
        "unordered_data_pairs": data_blocks * (data_blocks - 1) // 2,
        "distinct_nonzero_differences": len(counts),
        "maximum_difference_multiplicity": max(counts.values()),
        "runs_with_repetitions": sum(count > 1 for count in counts.values()),
        "repeated_records": sum(count - 1 for count in counts.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=12)
    parser.add_argument(
        "--scanner",
        type=Path,
        default=SCRIPT_DIRECTORY
        / "certify_shiftalpha64_data_difference_multiplicity.exe",
    )
    args = parser.parse_args()
    expected = direct(args.data_blocks)
    completed = subprocess.run(
        [str(args.scanner), str(args.data_blocks)],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    observed = {key: int(payload[key]) for key in expected}
    if expected != observed:
        raise SystemExit(f"difference audit mismatch: {expected} != {observed}")
    print(
        json.dumps(
            {
                "schema": "shiftalpha64-data-difference-multiplicity-audit-v1",
                "data_blocks": args.data_blocks,
                "direct": expected,
                "optimized": observed,
                "status": "PASS",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
