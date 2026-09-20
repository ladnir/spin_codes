#!/usr/bin/env python3
"""Exact geometric envelope for Riffle DP g=4 reset patterns."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WINDOW_LEDGER = ROOT / "explorations" / "riffle_dp_g4_g2_autonomous_window_bound.json"
OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_reset_pattern_placements.json"
INNER_NODES = 32_772
PACKETS_PER_NODE = 16
PACKET_COUNT = INNER_NODES * PACKETS_PER_NODE
DISTANCE_THRESHOLD = 188_766


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def node_slot_subsets(packet_support: int, occupied_nodes: int) -> int:
    return sum(
        (-1) ** (occupied_nodes - used_nodes)
        * math.comb(occupied_nodes, used_nodes)
        * math.comb(PACKETS_PER_NODE * used_nodes, packet_support)
        for used_nodes in range(occupied_nodes + 1)
    )


def reset_pattern_count(
    node_count: int,
    occupied_nodes: int,
    resets: int,
    maximum_live_length: int,
) -> int:
    live = occupied_nodes - resets
    if live <= 0 or resets > live:
        return 0
    maximum_sum = min(maximum_live_length, node_count - resets)
    if maximum_sum < live:
        return 0
    term = math.comb(node_count - live, resets)
    subtotal = 0
    for live_sum in range(live, maximum_sum + 1):
        subtotal += term
        if live_sum == maximum_sum:
            break
        numerator = term * live_sum * (node_count - live_sum - resets)
        denominator = (live_sum - live + 1) * (node_count - live_sum)
        term, remainder = divmod(numerator, denominator)
        if remainder:
            raise RuntimeError("reset envelope: nonintegral recurrence")
    return math.comb(live, resets) * subtotal


def brute_reset_pattern_count(
    node_count: int,
    occupied_nodes: int,
    resets: int,
    maximum_live_length: int,
) -> int:
    result = 0
    for positions in itertools.combinations(range(node_count), occupied_nodes):
        lengths = tuple(
            positions[index + 1] - positions[index]
            for index in range(occupied_nodes - 1)
        ) + (node_count - positions[-1],)
        for mask in range(1 << occupied_nodes):
            if mask & 1:
                continue
            reset_indices = [index for index in range(occupied_nodes) if (mask >> index) & 1]
            if len(reset_indices) != resets:
                continue
            if any(right == left + 1 for left, right in zip(reset_indices, reset_indices[1:])):
                continue
            live_sum = sum(
                length for index, length in enumerate(lengths) if index not in reset_indices
            )
            result += live_sum <= maximum_live_length
    return result


def log2_fraction(value: Fraction) -> float:
    if not value:
        return float("-inf")
    return math.log2(value.numerator) - math.log2(value.denominator)


def main() -> None:
    window = json.loads(WINDOW_LEDGER.read_text())
    if window["distance_consequences"]["distance_threshold"] != DISTANCE_THRESHOLD:
        raise RuntimeError("reset envelope: distance threshold mismatch")
    for node_count in range(3, 9):
        for occupied_nodes in range(1, min(4, node_count) + 1):
            for resets in range(occupied_nodes // 2 + 1):
                for maximum_live_length in range(occupied_nodes - resets, node_count + 1):
                    exact = reset_pattern_count(
                        node_count,
                        occupied_nodes,
                        resets,
                        maximum_live_length,
                    )
                    brute = brute_reset_pattern_count(
                        node_count,
                        occupied_nodes,
                        resets,
                        maximum_live_length,
                    )
                    if exact != brute:
                        raise RuntimeError("reset envelope: small-instance check failed")

    pattern_counts = {}
    for occupied_nodes in range(1, 39):
        for resets in range(occupied_nodes // 2 + 1):
            live = occupied_nodes - resets
            maximum_live_length = (
                3 * DISTANCE_THRESHOLD + 16 * live
            ) // 19
            pattern_counts[(occupied_nodes, resets)] = (
                maximum_live_length,
                reset_pattern_count(
                    INNER_NODES,
                    occupied_nodes,
                    resets,
                    maximum_live_length,
                ),
            )

    support_rows = []
    for packet_support in range(33, 39):
        denominator = math.comb(PACKET_COUNT, packet_support)
        reset_totals: dict[int, int] = {}
        occupied_rows = []
        for occupied_nodes in range(
            (packet_support + PACKETS_PER_NODE - 1) // PACKETS_PER_NODE,
            packet_support + 1,
        ):
            slot_subsets = node_slot_subsets(packet_support, occupied_nodes)
            if slot_subsets <= 0:
                continue
            for resets in range(occupied_nodes // 2 + 1):
                maximum_live_length, node_patterns = pattern_counts[
                    (occupied_nodes, resets)
                ]
                numerator = slot_subsets * node_patterns
                reset_totals[resets] = reset_totals.get(resets, 0) + numerator
                if numerator:
                    occupied_rows.append(
                        {
                            "occupied_nodes": occupied_nodes,
                            "resets": resets,
                            "live_nodes": occupied_nodes - resets,
                            "maximum_total_live_gap_length": maximum_live_length,
                            "node_slot_subsets": str(slot_subsets),
                            "node_set_reset_patterns": str(node_patterns),
                        }
                    )
        reset_rows = []
        total = Fraction()
        for resets, numerator in sorted(reset_totals.items()):
            probability = Fraction(numerator, denominator)
            total += probability
            reset_rows.append(
                {
                    "resets": resets,
                    "probability_numerator": str(probability.numerator),
                    "probability_denominator": str(probability.denominator),
                    "log2_upper": log2_fraction(probability),
                }
            )
        support_rows.append(
            {
                "packet_support": packet_support,
                "reset_rows": reset_rows,
                "union_upper": {
                    "numerator": str(total.numerator),
                    "denominator": str(total.denominator),
                    "log2": log2_fraction(total),
                    "capped_log2": min(0.0, log2_fraction(total)),
                },
                "occupied_node_rows": occupied_rows,
            }
        )

    payload = {
        "schema": "riffle-dp-g4-g2-reset-pattern-placement-envelope-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "RIGOROUS_BOUND",
        "autonomous_window_ledger_sha256": digest(WINDOW_LEDGER),
        "necessary_condition": (
            "If an output has weight at most d, then 19 times the total length "
            "of live autonomous gaps is at most 3*d+16*live_gap_count."
        ),
        "reset_model": (
            "The first occupied node is live. A reset node cannot follow a reset "
            "node. Reset compatibility with packet values is deliberately ignored."
        ),
        "support_rows": support_rows,
        "scope_limitation": (
            "This is a geometric union envelope, not the complete G2 bound. "
            "It must be multiplied by certified algebraic reset costs without "
            "double counting before use in the first-moment ledger."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    for row in support_rows:
        zero_reset = next(item for item in row["reset_rows"] if item["resets"] == 0)
        print(
            f"support_{row['packet_support']}_zero_reset_log2_upper="
            f"{zero_reset['log2_upper']:.12f}"
        )
        print(
            f"support_{row['packet_support']}_all_reset_patterns_capped_log2="
            f"{row['union_upper']['capped_log2']:.12f}"
        )
    print(f"output={OUTPUT}")
    print("status=RIGOROUS_G2_RESET_PATTERN_PLACEMENT_ENVELOPE")


if __name__ == "__main__":
    main()
