#!/usr/bin/env python3
"""Verify the closed boundary formula on packet weights zero, one, and three."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from probe_riffle_small_boundary_saddle import boundary_transition
from probe_riffle_small_typewise_bound import DEFAULT_RECEIPTS


def weak_compositions(items: int, boxes: int) -> int:
    if boxes == 0:
        return int(items == 0)
    return math.comb(items + boxes - 1, boxes - 1)


def closed_count(zero_packets: int, weight_one: int, weight_three: int) -> int:
    if (weight_one | weight_three) & 1:
        return 0
    one_excursions = weight_one // 2
    half_weight_three = weight_three // 2
    result = 0
    for three_excursions in range(half_weight_three + 1):
        internal_pairs = half_weight_three - three_excursions
        internal_allocations = weak_compositions(
            internal_pairs, one_excursions
        )
        if not internal_allocations:
            continue
        excursions = one_excursions + three_excursions
        excursion_orders = math.comb(excursions, three_excursions)
        zero_allocations = math.comb(zero_packets + excursions, excursions)
        result += (
            4**excursions
            * 6**internal_pairs
            * internal_allocations
            * excursion_orders
            * zero_allocations
        )
    return result


def audit(maximum_length: int) -> tuple[int, int]:
    current: dict[tuple[int, int, int, int], int] = {(0, 0, 0, 0): 1}
    comparisons = 0
    maximum_count_bits = 1
    for length in range(maximum_length + 1):
        for zero_packets in range(length + 1):
            for weight_one in range(length - zero_packets + 1):
                weight_three = length - zero_packets - weight_one
                exact = current.get(
                    (zero_packets, weight_one, weight_three, 0), 0
                )
                expected = closed_count(
                    zero_packets, weight_one, weight_three
                )
                if exact != expected:
                    raise RuntimeError(
                        "odd-face formula mismatch: "
                        f"N={length}, h0={zero_packets}, h1={weight_one}, "
                        f"h3={weight_three}, exact={exact}, expected={expected}"
                    )
                comparisons += 1
                maximum_count_bits = max(maximum_count_bits, exact.bit_length())

        if length == maximum_length:
            break
        following: defaultdict[tuple[int, int, int, int], int] = defaultdict(int)
        for (zero_packets, weight_one, weight_three, old), ways in current.items():
            for packet, coordinate in ((0, 0), (1, 1), (3, 2)):
                for new in range(5):
                    multiplicity = boundary_transition(old, new, packet)
                    if not multiplicity:
                        continue
                    counts = [zero_packets, weight_one, weight_three]
                    counts[coordinate] += 1
                    following[(*counts, new)] += ways * multiplicity
        current = dict(following)
    return comparisons, maximum_count_bits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-length", type=int, default=48)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.maximum_length < 0:
        raise ValueError("maximum length must be nonnegative")

    comparisons, maximum_count_bits = audit(args.maximum_length)
    payload = {
        "schema": "riffle-boundary-odd-face-formula-audit-v1",
        "evidence_label": "EXACT_INTEGER_DYNAMIC_PROGRAM_AUDIT",
        "packet_weights": [0, 1, 3],
        "maximum_length": args.maximum_length,
        "comparisons": comparisons,
        "maximum_exact_count_bits": maximum_count_bits,
        "formula_parameters": {
            "one_excursions": "b=h1/2",
            "three_excursions": "a=0,...,h3/2",
            "internal_weight_three_pairs": "t=h3/2-a",
        },
        "status": "PASS",
    }
    output = args.output or (
        DEFAULT_RECEIPTS / "goal24_boundary_odd_face_formula_audit.json"
    )
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output": str(output), **payload}, indent=2))


if __name__ == "__main__":
    main()
