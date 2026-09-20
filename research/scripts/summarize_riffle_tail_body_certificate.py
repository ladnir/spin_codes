#!/usr/bin/env python3
"""Combine the five disjoint tail/body bounds in a Riffle certificate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def log2_sum(log2_values: list[float]) -> float:
    pivot = max(log2_values)
    return pivot + math.log2(sum(2.0 ** (value - pivot) for value in log2_values))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--one-tail", type=Path, required=True)
    parser.add_argument("--tail-body", type=Path, required=True)
    parser.add_argument("--two-tail", type=Path, required=True)
    parser.add_argument("--mixed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    one_tail = load(args.one_tail)
    tail_body = load(args.tail_body)
    two_tail = load(args.two_tail)
    mixed = load(args.mixed)
    component_log2 = {
        "exactly_one_tail_no_body": float(one_tail["aggregate_log2_upper"]),
        "exactly_two_tail_no_body": float(two_tail["aggregate_log2_upper"]),
        "three_or_more_tail_no_body": float(
            tail_body["tail_only"]["three_plus_count_only_log2_upper"]
        ),
        "body_no_tail": float(tail_body["body_only_log2_upper"]),
        "body_and_tail": float(mixed["aggregate_log2_upper"]),
    }
    aggregate = log2_sum(list(component_log2.values()))
    components = {
        name: {
            "log2_upper": value,
            "margin_bits": -value,
            "fraction_of_aggregate": 2.0 ** (value - aggregate),
        }
        for name, value in component_log2.items()
    }
    payload = {
        "schema": "riffle-tail-body-certificate-summary-v1",
        "aggregate_log2_upper": aggregate,
        "aggregate_margin_bits": -aggregate,
        "components": components,
        "inputs": {
            "one_tail": str(args.one_tail),
            "tail_body": str(args.tail_body),
            "two_tail": str(args.two_tail),
            "mixed": str(args.mixed),
        },
        "scope": (
            "The five components are disjoint and exhaustive under the recorded "
            "endpoint-tail split. The component analyzers use nearest binary64; "
            "this summary is not an outward-rounded certificate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"margin_bits,{payload['aggregate_margin_bits']:.12f}")
    for name, row in components.items():
        print(
            f"component,{name},margin_bits,{row['margin_bits']:.12f},"
            f"fraction,{row['fraction_of_aggregate']:.12g}"
        )


if __name__ == "__main__":
    main()
