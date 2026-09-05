#!/usr/bin/env python3
"""Audit the optimized p1 endpoint scanner on a reduced instance."""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as field  # noqa: E402
from analyze_riffle_bch_alpha_joint_spectrum import encode_weights  # noqa: E402


ALL_ONE_MESSAGE = 0x8E63CF44EFD4FA21


def weight(value: int) -> int:
    return int(encode_weights(np.asarray([np.uint64(value)]))[0])


def add(histogram: dict[int, int], value: int) -> None:
    encoded_weight = weight(value)
    histogram[encoded_weight] = histogram.get(encoded_weight, 0) + 1


def direct_histograms(data_blocks: int) -> tuple[dict[int, int], dict[int, int], int]:
    finite = [0] + [
        field.field_power(2, 64 + index) for index in range(data_blocks)
    ]
    finite_histogram: dict[int, int] = {}
    p1_histogram: dict[int, int] = {}
    overlaps = 0
    for left, right in itertools.combinations(finite, 2):
        difference = left ^ right
        finite_value = field.field_multiply(ALL_ONE_MESSAGE, difference)
        p1_value = field.field_multiply(
            ALL_ONE_MESSAGE, field.field_inverse(difference)
        )
        add(finite_histogram, finite_value)
        add(p1_histogram, p1_value)
        overlaps += int(finite_value == ALL_ONE_MESSAGE)
    return finite_histogram, p1_histogram, overlaps


def normalize_histogram(payload: dict[str, int]) -> dict[int, int]:
    return {int(weight): int(count) for weight, count in payload.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=12)
    parser.add_argument(
        "--scanner",
        type=Path,
        default=SCRIPT_DIRECTORY / "scan_riffle_shiftalpha64_p1_endpoints.exe",
    )
    args = parser.parse_args()
    completed = subprocess.run(
        [str(args.scanner), str(args.data_blocks)],
        check=True,
        capture_output=True,
        text=True,
    )
    optimized = json.loads(completed.stdout)
    finite_histogram, p1_histogram, overlaps = direct_histograms(args.data_blocks)
    if normalize_histogram(
        optimized["finite_coordinate_normalization_histogram"]
    ) != finite_histogram:
        raise SystemExit("finite-coordinate histogram mismatch")
    if normalize_histogram(
        optimized["p1_coordinate_normalization_histogram"]
    ) != p1_histogram:
        raise SystemExit("p1-coordinate histogram mismatch")
    if int(optimized["all_three_all_one_overlap_supports"]) != overlaps:
        raise SystemExit("all-one overlap mismatch")
    payload = {
        "schema": "riffle-shiftalpha64-p1-endpoint-audit-v1",
        "data_blocks": args.data_blocks,
        "finite_coordinate_normalization_histogram": finite_histogram,
        "p1_coordinate_normalization_histogram": p1_histogram,
        "all_three_all_one_overlap_supports": overlaps,
        "status": "PASS",
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
