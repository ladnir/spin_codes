#!/usr/bin/env python3
"""Verify the closed boundary formula on every one-dimensional packet face."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from probe_riffle_small_boundary_saddle import boundary_transition
from probe_riffle_small_typewise_bound import DEFAULT_RECEIPTS


def closed_count(zero_packets: int, packet_weight: int, active_packets: int) -> int:
    if active_packets & 1:
        return 0
    excursions = active_packets // 2
    return (
        math.comb(4, packet_weight) ** excursions
        * math.comb(zero_packets + excursions, excursions)
    )


def audit_face(packet_weight: int, maximum_length: int) -> tuple[int, int]:
    current: dict[tuple[int, int, int], int] = {(0, 0, 0): 1}
    comparisons = 0
    maximum_count_bits = 1
    for length in range(maximum_length + 1):
        for zero_packets in range(length + 1):
            active_packets = length - zero_packets
            exact = current.get((zero_packets, active_packets, 0), 0)
            expected = closed_count(zero_packets, packet_weight, active_packets)
            if exact != expected:
                raise RuntimeError(
                    "one-face formula mismatch: "
                    f"k={packet_weight}, N={length}, h0={zero_packets}, "
                    f"hk={active_packets}, exact={exact}, expected={expected}"
                )
            comparisons += 1
            maximum_count_bits = max(maximum_count_bits, exact.bit_length())

        if length == maximum_length:
            break
        following: defaultdict[tuple[int, int, int], int] = defaultdict(int)
        for (zero_packets, active_packets, old), ways in current.items():
            zero_multiplicity = boundary_transition(old, old, 0)
            if zero_multiplicity:
                following[(zero_packets + 1, active_packets, old)] += (
                    ways * zero_multiplicity
                )
            for new in range(5):
                multiplicity = boundary_transition(old, new, packet_weight)
                if multiplicity:
                    following[(zero_packets, active_packets + 1, new)] += (
                        ways * multiplicity
                    )
        current = dict(following)
    return comparisons, maximum_count_bits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-length", type=int, default=64)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.maximum_length < 0:
        raise ValueError("maximum length must be nonnegative")

    rows = []
    for packet_weight in range(1, 5):
        comparisons, maximum_count_bits = audit_face(
            packet_weight, args.maximum_length
        )
        rows.append(
            {
                "packet_weight": packet_weight,
                "outbound_multiplicity": math.comb(4, packet_weight),
                "comparisons": comparisons,
                "maximum_exact_count_bits": maximum_count_bits,
            }
        )

    payload = {
        "schema": "riffle-boundary-one-face-formula-audit-v1",
        "evidence_label": "EXACT_INTEGER_DYNAMIC_PROGRAM_AUDIT",
        "maximum_length": args.maximum_length,
        "total_comparisons": sum(row["comparisons"] for row in rows),
        "formula": (
            "A((h0,0,...,hk,...,0),(k*hk)/2) = "
            "0 for odd hk; otherwise binom(4,k)^(hk/2) "
            "* binom(h0+hk/2,hk/2)"
        ),
        "faces": rows,
        "status": "PASS",
    }
    output = args.output or (
        DEFAULT_RECEIPTS / "goal23_boundary_one_face_formula_audit.json"
    )
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output": str(output), **payload}, indent=2))


if __name__ == "__main__":
    main()
