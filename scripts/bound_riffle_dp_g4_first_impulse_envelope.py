#!/usr/bin/env python3
"""Use the exact first impulse to strengthen the zero-reset G2 envelope."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ORBIT_RECEIPT = ROOT / "explorations" / "riffle_dp_g4_g2_single_packet_orbits.json"
COLLISION_S2 = ROOT / "explorations" / "riffle_dp_g4_g2_initial_collision_orbits_s2.json"
COLLISION_S3 = ROOT / "explorations" / "riffle_dp_g4_g2_initial_collision_orbits_s3.json"
WINDOW_LEDGER = ROOT / "explorations" / "riffle_dp_g4_g2_autonomous_window_bound.json"
OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_first_impulse_envelope.json"
INNER_NODES = 32_772
PACKETS_PER_NODE = 16
PACKET_COUNT = INNER_NODES * PACKETS_PER_NODE
DISTANCE_THRESHOLD = 188_766


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def node_slot_subsets(packet_support: int, occupied_nodes: int) -> int:
    if packet_support < occupied_nodes or packet_support > 16 * occupied_nodes:
        return 0
    return sum(
        (-1) ** (occupied_nodes - used_nodes)
        * math.comb(occupied_nodes, used_nodes)
        * math.comb(16 * used_nodes, packet_support)
        for used_nodes in range(occupied_nodes + 1)
    )


def log2_fraction(value: Fraction) -> float:
    if not value:
        return float("-inf")
    return math.log2(value.numerator) - math.log2(value.denominator)


def main() -> None:
    orbit = json.loads(ORBIT_RECEIPT.read_text())
    collision_s2 = json.loads(COLLISION_S2.read_text())
    collision_s3 = json.loads(COLLISION_S3.read_text())
    window = json.loads(WINDOW_LEDGER.read_text())
    minimum_prefix_s1 = orbit["minimum_autonomous_prefix_weights"]
    if len(minimum_prefix_s1) != INNER_NODES:
        raise RuntimeError("first impulse: prefix table length mismatch")
    first_crossing = next(
        index + 1
        for index, weight in enumerate(minimum_prefix_s1)
        if weight > DISTANCE_THRESHOLD
    )
    if first_crossing != orbit["first_step_above_distance_max"]:
        raise RuntimeError("first impulse: worst crossing mismatch")
    if window["distance_consequences"]["distance_threshold"] != DISTANCE_THRESHOLD:
        raise RuntimeError("first impulse: window threshold mismatch")
    if collision_s2["initial_packet_support"] != 2 or collision_s3["initial_packet_support"] != 3:
        raise RuntimeError("first impulse: collision support mismatch")
    prefix_tables = {
        1: minimum_prefix_s1[: first_crossing - 1],
        2: collision_s2["minimum_prefix_weights_through_last_safe_length"],
        3: collision_s3["minimum_prefix_weights_through_last_safe_length"],
    }
    crossing_by_support = {
        1: first_crossing,
        2: collision_s2["first_step_above_distance_max"],
        3: collision_s3["first_step_above_distance_max"],
    }
    if any(len(prefix_tables[support]) != crossing_by_support[support] - 1 for support in (1, 2, 3)):
        raise RuntimeError("first impulse: safe prefix length mismatch")

    bad_node_sets = {}
    for first_packets, minimum_prefix in prefix_tables.items():
        for occupied_nodes in range(2, 39):
            remaining_segments = occupied_nodes - 1
            count = 0
            for first_length, first_weight in enumerate(minimum_prefix, start=1):
                remaining_budget = DISTANCE_THRESHOLD - first_weight
                maximum_remaining_length = (
                    3 * remaining_budget + 16 * remaining_segments
                ) // 19
                maximum_remaining_length = min(
                    maximum_remaining_length,
                    INNER_NODES - first_length,
                )
                if maximum_remaining_length >= remaining_segments:
                    count += math.comb(maximum_remaining_length, remaining_segments)
            bad_node_sets[(first_packets, occupied_nodes)] = count

    support_rows = []
    for packet_support in range(33, 39):
        denominator = math.comb(PACKET_COUNT, packet_support)
        prefix_numerators = {1: 0, 2: 0, 3: 0}
        first_ge4_numerator = 0
        occupied_rows = []
        for occupied_nodes in range(3, packet_support + 1):
            prefix_rows = []
            for first_packets in (1, 2, 3):
                slots = math.comb(16, first_packets) * node_slot_subsets(
                    packet_support - first_packets,
                    occupied_nodes - 1,
                )
                count = slots * bad_node_sets[(first_packets, occupied_nodes)]
                prefix_numerators[first_packets] += count
                prefix_rows.append(
                    {
                        "first_packets": first_packets,
                        "bad_node_sets": str(
                            bad_node_sets[(first_packets, occupied_nodes)]
                        ),
                        "slot_subsets": str(slots),
                    }
                )
            ge4_slots = sum(
                math.comb(16, first_packets)
                * node_slot_subsets(
                    packet_support - first_packets,
                    occupied_nodes - 1,
                )
                for first_packets in range(4, min(16, packet_support) + 1)
            )
            ge4_count = math.comb(INNER_NODES, occupied_nodes) * ge4_slots
            first_ge4_numerator += ge4_count
            if any(int(row["slot_subsets"]) for row in prefix_rows) or ge4_count:
                occupied_rows.append(
                    {
                        "occupied_nodes": occupied_nodes,
                        "prefix_rows": prefix_rows,
                        "first_ge4_slot_subsets": str(ge4_slots),
                    }
                )
        prefix_probabilities = {
            support: Fraction(numerator, denominator)
            for support, numerator in prefix_numerators.items()
        }
        ge4_probability = Fraction(first_ge4_numerator, denominator)
        combined = sum(prefix_probabilities.values(), ge4_probability)
        support_rows.append(
            {
                "packet_support": packet_support,
                "exact_prefix_rows": [
                    {
                        "first_packets": first_packets,
                        "probability_numerator": str(probability.numerator),
                        "probability_denominator": str(probability.denominator),
                        "log2": log2_fraction(probability),
                    }
                    for first_packets, probability in prefix_probabilities.items()
                ],
                "first_node_support_at_least_4_unconditional_upper": {
                    "numerator": str(ge4_probability.numerator),
                    "denominator": str(ge4_probability.denominator),
                    "log2": log2_fraction(ge4_probability),
                },
                "combined_upper": {
                    "numerator": str(combined.numerator),
                    "denominator": str(combined.denominator),
                    "log2": log2_fraction(combined),
                },
                "occupied_node_rows": occupied_rows,
            }
        )

    payload = {
        "schema": "riffle-dp-g4-g2-first-impulse-envelope-v2",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "RIGOROUS_BOUND",
        "orbit_receipt_sha256": digest(ORBIT_RECEIPT),
        "collision_support_2_receipt_sha256": digest(COLLISION_S2),
        "collision_support_3_receipt_sha256": digest(COLLISION_S3),
        "autonomous_window_ledger_sha256": digest(WINDOW_LEDGER),
        "scope": (
            "Zero-reset placements. If the first occupied node contains one, "
            "two, or three packets, use its exact worst prefix and the autonomous "
            "linear relaxation afterward. Count first-node support at least four "
            "without claiming an output bound."
        ),
        "first_crossing_by_packet_support": {
            str(support): crossing for support, crossing in crossing_by_support.items()
        },
        "support_rows": support_rows,
        "scope_limitation": (
            "The support-at-least-four row is unconditional and the later live "
            "segments use the weak autonomous relaxation. This artifact does not close G2."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(
        "first_crossings="
        + ",".join(f"s{support}:{crossing}" for support, crossing in crossing_by_support.items())
    )
    for row in support_rows:
        prefix_text = ",".join(
            f"s{item['first_packets']}:{item['log2']:.12f}"
            for item in row["exact_prefix_rows"]
        )
        print(f"support_{row['packet_support']}_prefix_log2_uppers={prefix_text}")
        print(
            f"support_{row['packet_support']}_first_ge4_log2_upper="
            f"{row['first_node_support_at_least_4_unconditional_upper']['log2']:.12f}"
        )
        print(
            f"support_{row['packet_support']}_combined_log2_upper="
            f"{row['combined_upper']['log2']:.12f}"
        )
    print(f"output={OUTPUT}")
    print("status=RIGOROUS_G2_FIRST_IMPULSE_ENVELOPE")


if __name__ == "__main__":
    main()
