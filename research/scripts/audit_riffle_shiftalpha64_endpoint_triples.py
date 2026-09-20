#!/usr/bin/env python3
"""Audit the optimized endpoint scanner against direct small enumeration."""

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


def direct_count(data_blocks: int) -> tuple[int, dict[int, int]]:
    finite = [field.field_power(2, 64 + index) for index in range(data_blocks)]
    columns: list[int | None] = finite + [0, None]
    result = 0
    finite_histogram: dict[int, int] = {}
    for support in itertools.combinations(range(len(columns)), 3):
        selected = [columns[index] for index in support]
        if selected.count(None):
            infinity = selected.index(None)
            others = [index for index in range(3) if index != infinity]
            bases = [0, 0, 0]
            bases[others[0]] = 1
            bases[others[1]] = 1
            bases[infinity] = int(selected[others[0]]) ^ int(selected[others[1]])
        else:
            first, second, third = (int(value) for value in selected)
            bases = [second ^ third, first ^ third, first ^ second]
        for endpoint in range(3):
            scale = field.field_multiply(
                ALL_ONE_MESSAGE, field.field_inverse(bases[endpoint])
            )
            other = next(index for index in range(3) if index != endpoint)
            candidate_weight = weight(field.field_multiply(scale, bases[other]))
            if None not in selected:
                canonical_weight = min(candidate_weight, 128 - candidate_weight)
                finite_histogram[canonical_weight] = (
                    finite_histogram.get(canonical_weight, 0) + 1
                )
            if candidate_weight in (22, 106):
                result += 1
    return result, finite_histogram


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=12)
    parser.add_argument(
        "--scanner",
        type=Path,
        default=SCRIPT_DIRECTORY / "scan_riffle_shiftalpha64_endpoint_triples.exe",
    )
    args = parser.parse_args()
    completed = subprocess.run(
        [str(args.scanner), str(args.data_blocks)],
        check=True,
        capture_output=True,
        text=True,
    )
    optimized = json.loads(completed.stdout)
    direct, direct_histogram = direct_count(args.data_blocks)
    if optimized["exact_22_106_128_outer_words"] != direct:
        raise SystemExit(
            "endpoint audit mismatch: "
            f"optimized={optimized['exact_22_106_128_outer_words']} direct={direct}"
        )
    optimized_histogram: dict[int, int] = {}
    for key, value in optimized["candidate_weight_histogram"].items():
        weight_value = int(key)
        canonical_weight = min(weight_value, 128 - weight_value)
        optimized_histogram[canonical_weight] = (
            optimized_histogram.get(canonical_weight, 0) + int(value)
        )
    if optimized_histogram != direct_histogram:
        raise SystemExit("endpoint audit candidate histogram mismatch")
    payload = {
        "schema": "riffle-shiftalpha64-endpoint-triples-audit-v1",
        "data_blocks": args.data_blocks,
        "direct_endpoint_profiles": direct,
        "optimized_endpoint_profiles": optimized[
            "exact_22_106_128_outer_words"
        ],
        "finite_support_candidate_histogram": direct_histogram,
        "status": "PASS",
        "scope": "Direct enumeration audit at the declared reduced data-block count.",
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
